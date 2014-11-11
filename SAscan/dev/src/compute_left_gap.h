#ifndef __COMPUTE_LEFT_GAP_H_INCLUDED
#define __COMPUTE_LEFT_GAP_H_INCLUDED

#include <cstdio>
#include <vector>

#include "bitvector.h"
#include "gap_array.h"

//==============================================================================
// Given the gap array of the block (representation using 2 bytes per elements)
// and the gap array of the left half-block wrt right half-block (bitvector
// representation), compute the gap array (wrt tail) of the left half-block
// and write to a given file using v-byte encoding.
//
// The whole computation is performed under given ram budget. It is fully
// parallelized and uses asynchronous I/O as much as possible.
//==============================================================================
void compute_left_gap(long left_block_size, long right_block_size,
    gap_array_2n *block_gap, bitvector *left_block_gap_bv,
    std::string out_filename, long max_threads, long ram_budget) {
  long block_size = left_block_size + right_block_size;

  fprintf(stderr, "    Compute gap for left half-block:\n");
  fprintf(stderr, "      Computation: ");
  long double left_block_gap_compute_start = utils::wclock();
  gap_array_2n *left_gap = new gap_array_2n(left_block_size + 1);

  long left_gap_filled = 0L;
  long excess_ptr = 0L;
  long gap_val = block_gap->m_count[0];

  while (excess_ptr < (long)block_gap->m_excess.size() && block_gap->m_excess[excess_ptr] == 0L) {
    ++excess_ptr;
    gap_val += (1L << 16);
  }

  long left_half_block_current_gap = 0L;
  for (long i = 0; i < block_size; ++i) {
    left_half_block_current_gap += gap_val;

    if (left_block_gap_bv->get(i)) {
      ++left_half_block_current_gap;
    } else {
      left_gap->set_count(left_gap_filled++, left_half_block_current_gap);
      left_half_block_current_gap = 0L;
    }

    gap_val = block_gap->m_count[i + 1];
    while (excess_ptr < (long)block_gap->m_excess.size() && block_gap->m_excess[excess_ptr] == i + 1) {
      ++excess_ptr;
      gap_val += (1L << 16);
    }
  }

  left_half_block_current_gap += gap_val;
  left_gap->set_count(left_gap_filled++, left_half_block_current_gap);

  long double left_block_gap_compute_time = utils::wclock() - left_block_gap_compute_start;
  long double left_block_gap_compute_speed = (left_block_size / (1024.L * 1024)) / left_block_gap_compute_time;
  fprintf(stderr, "%.2Lf (%.2LfMiB/s)\n", left_block_gap_compute_time, left_block_gap_compute_speed);

  fprintf(stderr, "      Write to disk: ");
  long double left_gap_write_start = utils::wclock();
  left_gap->write_to_disk(out_filename);
  long double left_gap_write_time = utils::wclock() - left_gap_write_start;
  long double left_gap_write_io = (left_block_size / (1024.L * 1024)) / left_gap_write_time;
  fprintf(stderr, "%.2Lf (I/O: %.2LfMiB/s)\n", left_gap_write_time, left_gap_write_io);

  delete left_gap;
}

#endif  // __COMPUTE_LEFT_GAP_H_INCLUDED
