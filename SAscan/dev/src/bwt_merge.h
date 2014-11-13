// XXX The merge_bwt function should probably be tested more thoroughly.

#ifndef __BWT_MERGE_H
#define __BWT_MERGE_H

#include <thread>
#include <algorithm>

#include "bitvector.h"
#include "ranksel_support.h"


//==============================================================================
// Compute bwt[beg..end).
//==============================================================================
void merge_bwt_aux(long beg, long end, long left_ptr, long right_ptr,
    unsigned char *left_bwt, unsigned char *right_bwt, unsigned char *bwt,
    bitvector *bv) {
  for (long i = beg; i < end; ++i) {
    if (bv->get(i)) bwt[i] = right_bwt[right_ptr++];
    else bwt[i] = left_bwt[left_ptr++];
  }
}

void compute_initial_rank(long i, ranksel_support *ranksel, long &result) {
  result = ranksel->rank(i);
}

//==============================================================================
// Merge partial bwt of half-blocks (of size left_size and right_size) into
// partial bwt of the whole block.
//==============================================================================
long merge_bwt(unsigned char *left_bwt, unsigned char *right_bwt,
    long left_size, long right_size, long left_block_i0, long right_block_i0,
    unsigned char left_block_last, unsigned char *bwt, bitvector *bv,
    long max_threads) {
  long block_size = left_size + right_size;

  ranksel_support *bv_ranksel = new ranksel_support(bv, block_size, max_threads);

  long max_range_size = (block_size + max_threads - 1) / max_threads;
  long n_ranges = (block_size + max_range_size - 1) / max_range_size;

  // 3
  //
  // Compute starting parameters for each thread.
  long *left_ptr = new long[n_ranges];
  long *right_ptr = new long[n_ranges];
  long *rank_at_range_beg = new long[n_ranges];

  std::thread **threads = new std::thread*[n_ranges];
  for (long t = 0; t < n_ranges; ++t) {
    long range_beg = t * max_range_size;
    threads[t] = new std::thread(compute_initial_rank,
        range_beg, bv_ranksel, std::ref(rank_at_range_beg[t]));
  }

  for (long t = 0; t < n_ranges; ++t) threads[t]->join();
  for (long t = 0; t < n_ranges; ++t) delete threads[t];

  for (long t = 0; t < n_ranges; ++t) {
    long range_beg = t * max_range_size;
    left_ptr[t] = range_beg - rank_at_range_beg[t];
    right_ptr[t] = rank_at_range_beg[t];
  }
  delete[] rank_at_range_beg;

  // 4
  //
  // Proper computation of bwt.
  for (long t = 0; t < n_ranges; ++t) {
    long range_beg = max_range_size * t;
    long range_end = std::min(range_beg + max_range_size, block_size);

    threads[t] = new std::thread(merge_bwt_aux, range_beg, range_end,
        left_ptr[t], right_ptr[t], left_bwt, right_bwt, bwt, bv);
  }

  for (long t = 0; t < n_ranges; ++t) threads[t]->join();
  for (long t = 0; t < n_ranges; ++t) delete threads[t];
  delete[] threads;
  delete[] left_ptr;
  delete[] right_ptr;

  // 5
  //
  // Find position j = select_1(bv, right_block_i0) and replace bwt[j] with
  // left_block_last. To speed up the search for j, we use sparse_rank.
  bwt[bv_ranksel->select1(right_block_i0)] = left_block_last;

  // 6
  //
  // Compute the returned value.
  long block_i0 = bv_ranksel->select0(left_block_i0);

  // 7
  //
  // Clean up and exit
  delete bv_ranksel;
  return block_i0;
}


#endif  // __BWT_MERGE_H
