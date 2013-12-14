// The implementation is intended to obey only the following limitions:
//  block_size <= 4GB
//  alphabet_size <= 254
//
// length could be anything up to 2^63 because we will encode partial SA using
// 4-byte integers (this restricts the block size to 4GB) and gap array using
// v-byte encoding. This requires slightly more complicated merging, but this
// we can afford.
//
// TODO: use uint40 for output on disk, long (8 bytes) for integers (in
//       general) in memory and unsigned (4 bytes) for inmemory block-related
//       integers.
// TODO: change the encoding of partial SA on disk, so they can *always* use
//       4-byte integers (even for huge texts), and during the merging we add
//       the missing bits.
// TODO: better memory management (less allocation).
// TODO: streaming left-to-right, not right-to-left, in general: test whether
//       it makes a differece.
// TODO: make the rightmost block the smallest.
// TODO: implement Juha's modified rank (and also test other ranks).
// TODO: check mem usage using 'memusage' and check leaks using valgrind
//       is it really using 5b bytes for a block of legth b?

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <algorithm>

#include "divsufsort.h"
#include "utils.h"
#include "rank.h"
#include "srank.h"
#include "gap_array.h"
#include "merge.h"
#include "bitvector.h"
#include "stream.h"

void FGM(std::string filename, long max_block_size) {
  long length = utils::file_size(filename);
  long end = length, block_id = 0, prev_end = length;
  text_reader *reader = new text_reader(filename);
  long double start = utils::wclock();
  while (end > 0) {
    long beg = std::max(end - max_block_size, 0L);
    long block_size = end - beg; // B = text[beg..end), current block
    fprintf(stderr, "Processing block %ld/%ld [%ld..%ld):\n", block_id + 1,
      (length + max_block_size - 1) / max_block_size, beg, end);

    // 1. Read current and previously processed block.
    fprintf(stderr, "  Reading blocks: ");
    long double read_start = utils::wclock();
    unsigned char *B = new unsigned char[block_size];
    reader->read_block(beg, block_size, B);
    unsigned char last = B[block_size - 1];
    unsigned char *extprevB = new unsigned char[max_block_size + 1];
    long ext_prev_block_size = prev_end - end + 1;
    reader->read_block(end - 1, ext_prev_block_size, extprevB);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - read_start);

    // 2. Compute gt_eof.
    fprintf(stderr, "  Compute gt_eof: ");
    long double gt_eof_start = utils::wclock();
    bitvector *gt_eof = new bitvector(block_size);
    bitvector *gt_head_bv = end < length ? new bitvector("gt_head") : NULL;
    compute_gt_eof(extprevB, ext_prev_block_size, B, block_size, gt_head_bv, gt_eof);
    delete gt_head_bv;
    delete[] extprevB;
    fprintf(stderr, "%.2Lf\n", utils::wclock() - gt_eof_start);

    // 2. Remap symbols of B to compute ordering of suffixes of BA
    // starting in B and then restore original B.
    fprintf(stderr, "  Computing corrected SA: ");
    long double sa_start = utils::wclock();
    for (long j = 0; j < block_size; ++j) B[j] += gt_eof->get(j);
    int *SA = new int[block_size];
    divsufsort(B, SA, (int)block_size);
    for (long j = 0; j < block_size; ++j) B[j] -= gt_eof->get(j);
    delete gt_eof;
    fprintf(stderr, "%.2Lf\n", utils::wclock() - sa_start);

    // 3. Store the head of new gt bitvector to disk.
    fprintf(stderr, "  Store new_gt_head to disk: ");
    long double new_gt_store_start = utils::wclock();
    long whole_suffix_pos = 0;
    for (long k = 0; k < block_size; ++k)
      if (!SA[k]) { whole_suffix_pos = k; break; }
    bitvector *new_gt_head = new bitvector(block_size);
    for (long k = whole_suffix_pos + 1; k < block_size; ++k)
      new_gt_head->set(block_size - 1 - SA[k]);
    new_gt_head->save("new_gt_head");
    delete new_gt_head;
    fprintf(stderr, "%.2Lf\n", utils::wclock() - new_gt_store_start);

    // 4. Store partial SA on disk.
    fprintf(stderr, "  Write partial SA to disk: ");
    long double write_sa_start = utils::wclock();
    for (long k = 0; k < block_size; ++k) SA[k] += beg;
    utils::write_ints_to_file(SA, block_size, "sparseSA." + utils::intToStr(block_id));
    for (long k = 0; k < block_size; ++k) SA[k] -= beg;
    fprintf(stderr, "%.2Lf\n", utils::wclock() - write_sa_start);

    // 5. Compute the BWT from SA and build rank on top of it.
    fprintf(stderr, "  Compute BWT and rank: ");
    long double compute_bwt_start = utils::wclock();
    long count[256] = {0}, dollar_pos = 0;
    if (block_id) {
      for (long j = 0; j < block_size; ++j) count[(int)B[j] + 1]++;
      for (int j = 1; j < 256; ++j) count[j] += count[j - 1];
      unsigned char *tmpBWT = (unsigned char *)SA;
      for (long j = 0, jj = 0; j < block_size; ++j)
        if (SA[j] == 0) dollar_pos = j;
        else tmpBWT[jj++] = B[SA[j] - 1];
      std::copy(tmpBWT, tmpBWT + block_size - 1, B);
    }
    delete[] SA;
    rank_4n *rank = block_id ? new rank_4n(B, block_size - 1) : NULL;
    delete[] B;
    fprintf(stderr, "%.2Lf\n", utils::wclock() - compute_bwt_start);

    // 6. Allocate the gap array, do the streaming and store gap to disk.
    fprintf(stderr, "  Stream:\r");
    long double stream_start = utils::wclock();
    buffered_gap_array *gap = new buffered_gap_array(block_size + 1);
    bit_stream_writer *new_gt_tail = new bit_stream_writer("new_gt_tail");
    bit_stream_reader *gt_tail = prev_end < length ? new bit_stream_reader("gt_tail") : NULL;
    long i = 0;
    unsigned char next_gt = 0;
    reader->init_backward_streaming();
    for (long j = length - 1, dbg = 0; j >= prev_end; --j, ++dbg) { // stream the tail
      if (dbg == (1 << 20)) {
        long double elapsed = utils::wclock() - stream_start;
        long double streamed_mib = (1.L * (length - j)) / (1 << 20);
        fprintf(stderr,"  Stream: %.1Lf%%. Time: %.2Lf. Speed: %.2LfMiB/s\r",
          (100.L * (length - j)) / (length - end), elapsed, streamed_mib / elapsed);
        dbg = 0;
      }
      unsigned char c = reader->read_next();        // c = text[j]
      i = count[c] + rank->rank(i - (i > dollar_pos), c);
      if (c == last && next_gt) ++i;               // next_gt = gt[j + 1]
      new_gt_tail->write(i > whole_suffix_pos);
      gap->increment(i);
      next_gt = gt_tail->read();
    }
    delete gt_tail;
    bit_stream_reader *gt_head = end < length ? new bit_stream_reader("gt_head") : NULL;
    for (long j = prev_end - 1, dbg = 0; j >= end; --j, ++dbg) { // stream the head
      if (dbg == (1 << 20)) {
        long double elapsed = utils::wclock() - stream_start;
        long double streamed_mib = (1.L * (length - j)) / (1 << 20);
        fprintf(stderr,"  Stream: %.1Lf%%. Time: %.2Lf. Speed: %.2LfMiB/s\r",
          (100.L * (length - j)) / (length - end), elapsed, streamed_mib / elapsed);
        dbg = 0;
      }
      unsigned char c = reader->read_next();    // c = text[j]
      i = count[c] + rank->rank(i - (i > dollar_pos), c);
      if (c == last && next_gt) ++i;            // next_gt = gt[j + 1]
      new_gt_tail->write(i > whole_suffix_pos);
      gap->increment(i);
      next_gt = gt_head->read();
    }
    long double stream_time = utils::wclock() - stream_start;
    long double streamed_mib = (1.L * (length - end)) / (1 << 20);
    fprintf(stderr,"  Stream: 100.0%%. Time: %.2Lf. Speed: %.2LfMiB/s\n",
      stream_time, streamed_mib / stream_time);
    delete gt_head;
    delete new_gt_tail;
    fprintf(stderr, "  gap->excess.size() = %d\n", (int)gap->excess.size());
    gap->save_to_file("gap." + utils::intToStr(block_id));

    // 7. Clean up.
    delete gap;
    delete rank;
    prev_end = end;
    end = beg;
    ++block_id;
    utils::execute("mv new_gt_head gt_head");
    utils::execute("mv new_gt_tail gt_tail");
  }

  delete reader;

  // Merge gap and sparseSA arrays into final SA.
  std::string out_filename = filename + ".sa";
  merge(length, max_block_size, out_filename);

  // Delete auxiliary files.
  for (int i = 0; i < block_id; ++i) {
    utils::file_delete("sparseSA." + utils::intToStr(i));
    utils::file_delete("gap." + utils::intToStr(i));
  }
  if (utils::file_exists("gt_head")) utils::file_delete("gt_head");
  if (utils::file_exists("gt_tail")) utils::file_delete("gt_tail");

  long double total_time = utils::wclock() - start;
  long double speed = total_time / ((1.L * length) / (1 << 20));
  fprintf(stderr, "Total time: %.2Lfs. Speed: %.2Lfs/MiB\n",
      total_time, speed);
}

