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

void compute_partial_sa_and_bwt(unsigned char*, long, long, std::string, long,
    long, bitvector*, unsigned char**);

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
    fprintf(stderr, "  Reading block: ");
    long double read_start = utils::wclock();
    unsigned char *B = new unsigned char[block_size];
    utils::read_block(filename, length - block_size - beg, block_size, B);
    std::reverse(B, B + block_size);
    unsigned char last = B[block_size - 1];
    fprintf(stderr, "%.2Lf\n", utils::wclock() - read_start);

    long count[256] = {0};
    bitvector *gt_eof_bv = NULL;
    if (block_id + 1 != n_block) {
      // 2a. Compute symbols counts of B.
      fprintf(stderr, "  Compute counts: ");
      long double compute_counts_start = utils::wclock();
      for (long j = 0; j < block_size; ++j) count[(int)B[j] + 1]++;
      for (int j = 1; j < 256; ++j) count[j] += count[j - 1];
      fprintf(stderr, "%.2Lf\n", utils::wclock() - compute_counts_start);

      // 2b. Read previous block.
      fprintf(stderr, "  Read previous block: ");
      long double prev_block_reading_start = utils::wclock();
      unsigned char *extprevB = new unsigned char[max_block_size + 1];
      long ext_prev_block_size = prev_end - end + 1;
      utils::read_block(filename, length - ext_prev_block_size - (end - 1),
          ext_prev_block_size, extprevB);
      std::reverse(extprevB, extprevB + ext_prev_block_size);
      fprintf(stderr, "%.2Lf\n", utils::wclock() - prev_block_reading_start);

      // 2c. Compute gt_eof.
      fprintf(stderr, "  Compute gt_eof_bv: ");
      long double gt_eof_start = utils::wclock();
      gt_eof_bv = new bitvector(block_size);
      bitvector *gt_head_bv = new bitvector(filename + ".gt_head");
      compute_gt_eof_bv(extprevB, ext_prev_block_size, B, block_size, gt_head_bv, gt_eof_bv);
      delete gt_head_bv;
      delete[] extprevB;
      fprintf(stderr, "%.2Lf\n", utils::wclock() - gt_eof_start);

      // 2d. Remap symbols of B.
      fprintf(stderr, "  Remapping B: ");
      long double remap_start = utils::wclock();
      for (long j = 0; j < block_size; ++j) B[j] += gt_eof_bv->get(j);
      fprintf(stderr, "%.2Lf\n", utils::wclock() - remap_start);
    }

    // 3. Compute the head of the new gt bitvector.
    long whole_suffix_rank = 0;
    if (n_block > 1) {
      fprintf(stderr, "  Compute new_gt_head_bv: ");
      long double new_gt_head_bv_start = utils::wclock();
      bitvector *new_gt_head_bv = new bitvector(block_size);
      whole_suffix_rank = compute_new_gt_head_bv(B, block_size, new_gt_head_bv);
      new_gt_head_bv->save(filename + ".new_gt_head");
      delete new_gt_head_bv;
      fprintf(stderr, "%.2Lf\n", utils::wclock() - new_gt_head_bv_start);
    }

    // 4. Compute/save partial SA. Compute BWT if it's not the last block.
    unsigned char *BWT = NULL;
    compute_partial_sa_and_bwt(B, block_size, max_block_size, filename, block_id, n_block, gt_eof_bv, &BWT);

    if (block_id + 1 != n_block) {
      // 5a. Build the rank support for BWT.
      fprintf(stderr, "  Building the rank data structure:\n");
      long double building_rank_start = utils::wclock();
      context_rank_4n *rank = new context_rank_4n(BWT, block_size - 1);
      delete[] BWT;
      fprintf(stderr, "  Building the rank data structure: %.2Lf\n", utils::wclock() - building_rank_start);

      // 5b. Allocate the gap array, do the streaming and store gap to disk.
      fprintf(stderr, "  Stream:\r");
      long double stream_start = utils::wclock();
      buffered_gap_array *gap = new buffered_gap_array(block_size + 1);
      bit_stream_writer *new_gt_tail = new bit_stream_writer(filename + ".new_gt_tail");
      bit_stream_reader *gt_tail = prev_end < length ? new bit_stream_reader(filename + ".gt_tail") : NULL;
      long i = 0;
      unsigned char next_gt = 0;
      stream_reader<unsigned char> *text_streamer = new stream_reader<unsigned char>(filename, 1 << 20);
      for (long j = length - 1, dbg = 0; j >= prev_end; --j, ++dbg) { // stream the tail
        if (dbg == (1 << 23)) {
          long double elapsed = utils::wclock() - stream_start;
          long double streamed_mib = (1.L * (length - j)) / (1 << 20);
          fprintf(stderr,"  Stream: %.1Lf%%. Time: %.2Lf. Speed: %.2LfMiB/s\r",
            (100.L * (length - j)) / (length - end), elapsed, streamed_mib / elapsed);
          dbg = 0;
        }
        unsigned char c = text_streamer->read();     // c = text^R[j]
        i = count[c] + rank->rank(i - (i > whole_suffix_rank), c);
        if (c == last && next_gt) ++i;               // next_gt = gt[j + 1]
        new_gt_tail->write(i > whole_suffix_rank);
        gap->increment(i);
        next_gt = gt_tail->read();
      }
      delete gt_tail;
      bit_stream_reader *gt_head = new bit_stream_reader(filename + ".gt_head");
      for (long j = prev_end - 1, dbg = 0; j >= end; --j, ++dbg) { // stream the head
        if (dbg == (1 << 23)) {
          long double elapsed = utils::wclock() - stream_start;
          long double streamed_mib = (1.L * (length - j)) / (1 << 20);
          fprintf(stderr,"  Stream: %.1Lf%%. Time: %.2Lf. Speed: %.2LfMiB/s\r",
            (100.L * (length - j)) / (length - end), elapsed, streamed_mib / elapsed);
          dbg = 0;
        }
        unsigned char c = text_streamer->read();  // c = text^R[j]
        i = count[c] + rank->rank(i - (i > whole_suffix_rank), c);
        if (c == last && next_gt) ++i;            // next_gt = gt[j + 1]
        new_gt_tail->write(i > whole_suffix_rank);
        gap->increment(i);
        next_gt = gt_head->read();
      }
      long double stream_time = utils::wclock() - stream_start;
      long double streamed_mib = (1.L * (length - end)) / (1 << 20);
      fprintf(stderr,"  Stream: 100.0%%. Time: %.2Lf. Speed: %.2LfMiB/s\n",
        stream_time, streamed_mib / stream_time);
      delete text_streamer;
      delete gt_head;
      delete new_gt_tail;
      delete rank;
      utils::execute("mv " + filename + ".new_gt_tail " + filename + ".gt_tail");

      // 5c. Save gap to file.
      gap->save_to_file(filename + ".gap." + utils::intToStr(block_id));
      delete gap;
    }

    prev_end = end;
    end = beg;
    --block_id;
    if (n_block > 1)
      utils::execute("mv " + filename + ".new_gt_head " + filename + ".gt_head");
  }

  if (utils::file_exists(filename + ".gt_head")) utils::file_delete(filename + ".gt_head");
  if (utils::file_exists(filename + ".gt_tail")) utils::file_delete(filename + ".gt_tail");
}

// Compute partial SA of B[0..block_size) and store on disk. Also
// compute the BWT iff block_id != n_block.
void compute_partial_sa_and_bwt(unsigned char *B, long block_size, long max_block_size,
    std::string text_fname, long block_id, long n_block,
    bitvector *gt_eof_bv, unsigned char **BWT) {
  if (block_size <= 2147483647L /*2^31-1*/) {
    // Easy case, just use use 32-bit divsufsort to sufsort.
    fprintf(stderr, "  Computing partial SA: ");
    long double sa_start = utils::wclock();
    int *SA = new int[block_size];
    divsufsort(B, SA, (int)block_size);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - sa_start);

    fprintf(stderr, "  Writing partial SA to disk: ");
    long double writing_sa_start = utils::wclock();
    std::string sa_fname = text_fname + ".partial_sa." + utils::intToStr(block_id);
    if (n_block > 1 && max_block_size <= 2147483647L) utils::write_objects_to_file<int>(SA, block_size, sa_fname);
    else {
      stream_writer<uint40> *sa_writer = new stream_writer<uint40>(sa_fname, 1 << 20);
      for (long j = 0; j < block_size; ++j) sa_writer->write(uint40((unsigned long)SA[j]));
      delete sa_writer;
    }
    fprintf(stderr, "%.2Lf\n", utils::wclock() - writing_sa_start);

    if (block_id + 1 != n_block) {
      // Remap symbols of B back to original.
      fprintf(stderr, "  Re-remapping B: ");
      long double reremap_start = utils::wclock();
      for (long j = 0; j < block_size; ++j) B[j] -= gt_eof_bv->get(j);
      delete gt_eof_bv;
      fprintf(stderr, "%.2Lf\n", utils::wclock() - reremap_start);

      // Compute BWT
      fprintf(stderr, "  Compute BWT: ");
      long double bwt_start = utils::wclock();
      unsigned char *tmp = (unsigned char *)SA;
      for (long j = 0, jj = 0; j < block_size; ++j)
        if (SA[j]) tmp[jj++] = B[SA[j] - 1];
      std::copy(tmp, tmp + block_size - 1, B);
      *BWT = B;
      delete[] SA;
      fprintf(stderr, "%.2Lf\n", utils::wclock() - bwt_start);
    } else {
      delete[] SA;
      delete[] B;
    }
  } else {
    /* NOT IMPLEMENTED YET */
    fprintf(stderr, "\nNOT IMPLEMENTED YET, exiting...\n");
    std::exit(EXIT_FAILURE);
  }
}

#endif // __PARTIAL_SUFSORT_H
