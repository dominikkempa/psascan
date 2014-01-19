#ifndef __PARTIAL_SUFSORT_H
#define __PARTIAL_SUFSORT_H

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

void partial_sufsort(std::string filename, long length, long max_block_size) {
  long n_block = (length + max_block_size - 1) / max_block_size;
  long block_id = n_block - 1, prev_end = length;

  while (block_id >= 0) {
    long beg = max_block_size * block_id;
    long end = std::min(length, beg + max_block_size);
    long block_size = end - beg; // B = text[beg..end), current block
    fprintf(stderr, "Processing block %ld/%ld [%ld..%ld):\n",
      n_block - block_id, n_block, beg, end);

    // 1. Read current and previously processed block.
    fprintf(stderr, "  Reading blocks: ");
    long double read_start = utils::wclock();
    unsigned char *B = new unsigned char[block_size];
    utils::read_block(filename, length - block_size - beg, block_size, B);
    std::reverse(B, B + block_size);
    unsigned char last = B[block_size - 1];
    unsigned char *extprevB = new unsigned char[max_block_size + 1];
    long ext_prev_block_size = prev_end - end + 1;
    utils::read_block(filename, length - ext_prev_block_size - (end - 1),
        ext_prev_block_size, extprevB);
    std::reverse(extprevB, extprevB + ext_prev_block_size);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - read_start);

    // 2. Compute symbols counts in the current block.
    fprintf(stderr, "  Compute counts: ");
    long double compute_counts_start = utils::wclock();
    long count[256] = {0};
    if (block_id + 1 != n_block) {
      for (long j = 0; j < block_size; ++j) count[(int)B[j] + 1]++;
      for (int j = 1; j < 256; ++j) count[j] += count[j - 1];
    }
    fprintf(stderr, "%.2Lf\n", utils::wclock() - compute_counts_start);

    // 3. Compute gt_eof.
    fprintf(stderr, "  Compute gt_eof_bv: ");
    long double gt_eof_start = utils::wclock();
    bitvector *gt_eof_bv = new bitvector(block_size);
    bitvector *gt_head_bv = end < length ? new bitvector(filename + ".gt_head") : NULL;
    compute_gt_eof_bv(extprevB, ext_prev_block_size, B, block_size, gt_head_bv, gt_eof_bv);
    delete gt_head_bv;
    delete[] extprevB;
    fprintf(stderr, "%.2Lf\n", utils::wclock() - gt_eof_start);

    // 4. Remap symbols of B.
    fprintf(stderr, "  Remapping B: ");
    long double remap_start = utils::wclock();
    for (long j = 0; j < block_size; ++j) B[j] += gt_eof_bv->get(j);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - remap_start);

    // 5. Compute the head of the new gt bitvector.
    fprintf(stderr, "  Compute new_gt_head_bv: ");
    long double new_gt_head_bv_start = utils::wclock();
    bitvector *new_gt_head_bv = new bitvector(block_size);
    long whole_suffix_rank = compute_new_gt_head_bv(B, block_size, new_gt_head_bv);
    new_gt_head_bv->save(filename + ".new_gt_head");
    delete new_gt_head_bv;
    fprintf(stderr, "%.2Lf\n", utils::wclock() - new_gt_head_bv_start);

    // 6. Compute ordering of suffixes of BA starting in B.
    fprintf(stderr, "  Computing partial SA: ");
    long double sa_start = utils::wclock();
    int *SA = new int[block_size];
    divsufsort(B, SA, (int)block_size);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - sa_start);

    // 7. Remapping B back.
    fprintf(stderr, "  Re-remapping B: ");
    long double reremap_start = utils::wclock();
    for (long j = 0; j < block_size; ++j) B[j] -= gt_eof_bv->get(j);
    delete gt_eof_bv;
    fprintf(stderr, "%.2Lf\n", utils::wclock() - reremap_start);

    // 8. Store partial SA on disk and compute BWT.
    fprintf(stderr, "  Write partial SA to disk and compute BWT: ");
    long double write_sa_and_compute_bwt_start = utils::wclock();
    std::string sa_fname = filename + ".partial_sa." + utils::intToStr(block_id);
    utils::write_ints_to_file(SA, block_size, sa_fname);
    if (block_id + 1 != n_block) {
      unsigned char *tmp = (unsigned char *)SA;
      for (long j = 0, jj = 0; j < block_size; ++j)
        if (SA[j]) tmp[jj++] = B[SA[j] - 1];
      std::copy(tmp, tmp + block_size - 1, B);
    }
    delete[] SA;
    fprintf(stderr, "%.2Lf\n", utils::wclock() - write_sa_and_compute_bwt_start);
    
    // 9. Build the rank query support BWT.
    fprintf(stderr, "  Building the rank data structure: ");
    long double building_rank_start = utils::wclock();
    context_rank_4n *rank = (block_id + 1 != n_block) ? new context_rank_4n(B, block_size - 1) : NULL;
    delete[] B;
    fprintf(stderr, "%.2Lf\n", utils::wclock() - building_rank_start);

    // 10. Allocate the gap array, do the streaming and store gap to disk.
    fprintf(stderr, "  Stream:\r");
    long double stream_start = utils::wclock();
    buffered_gap_array *gap = new buffered_gap_array(block_size + 1);
    bit_stream_writer *new_gt_tail = new bit_stream_writer(filename + ".new_gt_tail");
    bit_stream_reader *gt_tail = prev_end < length ? new bit_stream_reader(filename + ".gt_tail") : NULL;
    long i = 0;
    unsigned char next_gt = 0;
    stream_reader<unsigned char> *streamer = new stream_reader<unsigned char>(filename, 1 << 20);
    for (long j = length - 1, dbg = 0; j >= prev_end; --j, ++dbg) { // stream the tail
      if (dbg == (1 << 23)) {
        long double elapsed = utils::wclock() - stream_start;
        long double streamed_mib = (1.L * (length - j)) / (1 << 20);
        fprintf(stderr,"  Stream: %.1Lf%%. Time: %.2Lf. Speed: %.2LfMiB/s\r",
          (100.L * (length - j)) / (length - end), elapsed, streamed_mib / elapsed);
        dbg = 0;
      }
      unsigned char c = streamer->read();          // c = text^R[j]
      i = count[c] + rank->rank(i - (i > whole_suffix_rank), c);
      if (c == last && next_gt) ++i;               // next_gt = gt[j + 1]
      new_gt_tail->write(i > whole_suffix_rank);
      gap->increment(i);
      next_gt = gt_tail->read();
    }
    delete gt_tail;
    bit_stream_reader *gt_head = end < length ? new bit_stream_reader(filename + ".gt_head") : NULL;
    for (long j = prev_end - 1, dbg = 0; j >= end; --j, ++dbg) { // stream the head
      if (dbg == (1 << 23)) {
        long double elapsed = utils::wclock() - stream_start;
        long double streamed_mib = (1.L * (length - j)) / (1 << 20);
        fprintf(stderr,"  Stream: %.1Lf%%. Time: %.2Lf. Speed: %.2LfMiB/s\r",
          (100.L * (length - j)) / (length - end), elapsed, streamed_mib / elapsed);
        dbg = 0;
      }
      unsigned char c = streamer->read();       // c = text^R[j]
      i = count[c] + rank->rank(i - (i > whole_suffix_rank), c);
      if (c == last && next_gt) ++i;            // next_gt = gt[j + 1]
      new_gt_tail->write(i > whole_suffix_rank);
      gap->increment(i);
      next_gt = gt_head->read();
    }
    delete streamer;
    long double stream_time = utils::wclock() - stream_start;
    long double streamed_mib = (1.L * (length - end)) / (1 << 20);
    fprintf(stderr,"  Stream: 100.0%%. Time: %.2Lf. Speed: %.2LfMiB/s\n",
      stream_time, streamed_mib / stream_time);
    delete gt_head;
    delete new_gt_tail;
    

    // 11. Save gap.
    fprintf(stderr, "  gap->excess.size() = %lu\n", (size_t)gap->excess.size());
    gap->save_to_file(filename + ".gap." + utils::intToStr(block_id));
    delete gap;

    // 12. Clean up.
    delete rank;
    prev_end = end;
    end = beg;
    --block_id;
    utils::execute("mv " + filename + ".new_gt_head " + filename + ".gt_head");
    utils::execute("mv " + filename + ".new_gt_tail " + filename + ".gt_tail");
  }

  if (utils::file_exists(filename + ".gt_head")) utils::file_delete(filename + ".gt_head");
  if (utils::file_exists(filename + ".gt_tail")) utils::file_delete(filename + ".gt_tail");
}


#endif // __PARTIAL_SUFSORT_H
