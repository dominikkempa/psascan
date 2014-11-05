// XXX The merge_bwt function should probably be tested more thoroughly.

#ifndef __BWT_MERGE_H
#define __BWT_MERGE_H

#include <thread>
#include <algorithm>

#include "bitvector.h"


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


//==============================================================================
// Compute the number of 1-bits in bv[0..i) with the help of sparse_rank.
// Note:
// - i is an integer in the range from 0 to length of bv (inclusive),
// - sparse_rank[k] = number of 1-bits in bv[0..k * chunk_size),
//==============================================================================
void compute_rank_using_sparse_rank(long chunk_size, long i,
    long *sparse_rank, bitvector *bv, long &result) {
  long j = i / chunk_size;
  result = sparse_rank[j];
  j *= chunk_size;

  while (j < i)
    result += bv->get(j++);
}


//==============================================================================
// Find the largest position j such that the number of 1s in bv[0..j) is <= i.
// In other words, find the position of i-th 1-bit in bv (i = 0, 1, ..).
//==============================================================================
long compute_select1_using_sparse_rank(long i, long chunk_size, long n_chunks,
    long *sparse_rank, bitvector *bv) {
  // Fast-forward through chunks preceding the chunk with the answer.
  long j = 0L;
  while (j < n_chunks && sparse_rank[j + 1] <= i)
    ++j;

  long rank_j = sparse_rank[j];
  j *= chunk_size;

  // Slowly find the final position in a single chunk.
  while (rank_j + bv->get(j) <= i)
    rank_j += bv->get(j++);

  return j;
}

//==============================================================================
// Find the largest position j such that the number of 0s in bv[0..j) is <= i.
// In other words, find the position of i-th 0-bit in bv (i = 0, 1, ..).
//==============================================================================
long compute_select0_using_sparse_rank(long i, long chunk_size, long n_chunks,
    long *sparse_rank, bitvector *bv) {
  // Fast forward through chunks preceding the chunk with the answer.
  long j = 0L;
  while (j < n_chunks && ((j + 1) * chunk_size) - sparse_rank[j + 1] <= i)
    ++j;

  long zero_cnt_j = (j * chunk_size) - sparse_rank[j];
  j *= chunk_size;

  // Slowly find the final position in a single chunk.
  while (zero_cnt_j + (1 - bv->get(j)) <= i)
    zero_cnt_j += (1 - bv->get(j++));

  return j;
}


//==============================================================================
// Compute sparse_rank[group_beg..group_end).
//==============================================================================
void process_group_of_chunks(long group_beg, long group_end, long chunk_size,
    long *sparse_rank, bitvector *bv) {
  for (long chunk_id = group_beg; chunk_id < group_end; ++chunk_id) {
    long chunk_beg = chunk_id * chunk_size;
    long chunk_end = chunk_beg + chunk_size;

    sparse_rank[chunk_id] = bv->range_sum(chunk_beg, chunk_end);
  }
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

  long chunk_size = std::min((1L << 20), (block_size + max_threads - 1) / max_threads);
  long n_chunks = block_size / chunk_size;  // we exclude the last partial chunk
  long *sparse_rank = (long *)malloc((n_chunks + 1) * sizeof(long));

  // 1
  //
  // Compute the values of sparse_rank, that is, the sum of 1-bits inside
  // each chunk. Since there can be more chunks than threads, we split chunks
  // into groups and let each thread compute the sum of 1-bits inside the
  // group of chunks.
  long chunk_max_group_size = (n_chunks + max_threads - 1) / max_threads;
  long n_chunk_groups = (n_chunks + chunk_max_group_size - 1) / chunk_max_group_size;

  std::thread **threads = new std::thread*[n_chunk_groups];
  for (long t = 0; t < n_chunk_groups; ++t) {
    long chunk_group_beg = t * chunk_max_group_size;
    long chunk_group_end = std::min(chunk_group_beg + chunk_max_group_size, n_chunks);
    threads[t] = new std::thread(process_group_of_chunks, chunk_group_beg,
        chunk_group_end, chunk_size, sparse_rank, bv);
  }

  for (long t = 0; t < n_chunk_groups; ++t) threads[t]->join();
  for (long t = 0; t < n_chunk_groups; ++t) delete threads[t];
  delete[] threads;

  // 2
  //
  // Compute cumulative sum of sparse_rank.
  for (long i = 0, sum = 0L; i <= n_chunks; ++i) {
    long temp = sparse_rank[i];
    sparse_rank[i] = sum;
    sum += temp;
  }

  long max_range_size = (block_size + max_threads - 1) / max_threads;
  long n_ranges = (block_size + max_range_size - 1) / max_range_size;

  // 3
  //
  // Compute starting parameters for each thread.
  long *left_ptr = new long[n_ranges];
  long *right_ptr = new long[n_ranges];
  long *rank_at_range_beg = new long[n_ranges];

  threads = new std::thread*[n_ranges];
  for (long t = 0; t < n_ranges; ++t) {
    long range_beg = t * max_range_size;
    threads[t] = new std::thread(compute_rank_using_sparse_rank,
        chunk_size, range_beg, sparse_rank, bv, std::ref(rank_at_range_beg[t]));
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
  long j = compute_select1_using_sparse_rank(right_block_i0, chunk_size, n_chunks, sparse_rank, bv);
  bwt[j] = left_block_last;

  // 6
  //
  // Compute the returned value.
  long block_i0 = compute_select0_using_sparse_rank(left_block_i0, chunk_size, n_chunks, sparse_rank, bv);

  // 7
  //
  // Clean up and exit
  free(sparse_rank);
  return block_i0;
}


#endif  // __BWT_MERGE_H
