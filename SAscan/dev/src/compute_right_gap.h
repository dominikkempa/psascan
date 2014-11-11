#ifndef __COMPUTE_RIGHT_GAP_H_INCLUDED
#define __COMPUTE_RIGHT_GAP_H_INCLUDED

#include <cstdio>
#include <vector>
#include <thread>
#include <mutex>
#include <future>
#include <condition_variable>
#include <algorithm>

#include "bitvector.h"
#include "gap_array.h"
#include "parallel_utils.h"


//==============================================================================
// Compute the number of 1-bits in bv[0..i) with the help of sparse_rank.
// Note:
// - i is an integer in the range from 0 to length of bv (inclusive),
// - sparse_rank[k] = number of 1-bits in bv[0..k * chunk_size),
//==============================================================================
// XXX remove the x from the function name, right not it is like that only
//     because of name conflict
// XXX maybe we should have class: rank support that can be built over bitvector
//     in parallel to support fast rank queries.
//==============================================================================
void xcompute_rank_using_sparse_rank(long chunk_size, long i,
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
long xcompute_select1_using_sparse_rank(long i, long chunk_size, long n_chunks,
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
long xcompute_select0_using_sparse_rank(long i, long chunk_size, long n_chunks,
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
void xprocess_group_of_chunks(long group_beg, long group_end, long chunk_size,
    long *sparse_rank, bitvector *bv) {
  for (long chunk_id = group_beg; chunk_id < group_end; ++chunk_id) {
    long chunk_beg = chunk_id * chunk_size;
    long chunk_end = chunk_beg + chunk_size;

    sparse_rank[chunk_id] = bv->range_sum(chunk_beg, chunk_end);
  }
}


//==============================================================================
// Compute the range_gap values corresponging to bv[part_beg..part_end).
//==============================================================================
void handle_bv_part(long part_beg, long part_end, long range_beg, long chunk_size,
    long *range_gap, long *sparse_rank, gap_array_2n *block_gap, bitvector *bv,
    long &res_sum, long &res_rank) {
  size_t excess_ptr = std::lower_bound(block_gap->m_excess.begin(),
      block_gap->m_excess.end(), part_beg) - block_gap->m_excess.begin();

  // Initialize j.
  long j = part_beg;

  // Compute gap[j].
  long gap_j = block_gap->m_count[j];
  while (excess_ptr < block_gap->m_excess.size() && block_gap->m_excess[excess_ptr] == j) {
    ++excess_ptr;
    gap_j += (1L << 16);
  }

  // Initialize sum.
  long sum = gap_j;

  while (j != part_end - 1 && bv->get(j) == 0) {
    // Update j.
    ++j;

    // Compute gap[j].
    gap_j = block_gap->m_count[j];
    while (excess_ptr < block_gap->m_excess.size() && block_gap->m_excess[excess_ptr] == j) {
      ++excess_ptr;
      gap_j += (1L << 16);
    }

    // Update sum.
    sum += gap_j;
  }

  // Store gap[part_beg] + .. + gap[j] and bv.rank(part_beg) (== bv.rank(j)).
  res_sum = sum;
  xcompute_rank_using_sparse_rank(chunk_size, part_beg, sparse_rank, bv, res_rank);

  if (j == part_end - 1)
    return;

  sum = 0L;
  long range_gap_ptr = res_rank + 1;
  while (j != part_end - 1) {
    // Update j.
    ++j;

    // Compute gap[j].
    gap_j = block_gap->m_count[j];
    while (excess_ptr < block_gap->m_excess.size() && block_gap->m_excess[excess_ptr] == j) {
      ++excess_ptr;
      gap_j += (1L << 16);
    }

    // Update sum.
    sum += gap_j;

    // Update range_gap.
    if (bv->get(j) == 1) {
      range_gap[range_gap_ptr - range_beg] = sum;
      ++range_gap_ptr;
      sum = 0L;
    }
  }

  if (bv->get(j) == 0)
    range_gap[range_gap_ptr - range_beg] = sum;
}


void async_write_code(unsigned char* &slab, long &length, std::mutex &mtx,
    std::condition_variable &cv, bool &avail, bool &finished, std::string filename) {
  while (true) {
    // Wait until the passive buffer is available.
    std::unique_lock<std::mutex> lk(mtx);
    while (!avail && !finished)
      cv.wait(lk);

    if (!avail && finished) {
      // We're done, terminate the thread.
      lk.unlock();
      return;
    }
    lk.unlock();

    // Safely write the data to disk.
    utils::add_objects_to_file(slab, length, filename);

    // Let the caller know what the I/O thread finished writing.
    lk.lock();
    avail = false;
    lk.unlock();
    cv.notify_one();
  }
}


//==============================================================================
// Given the gap array of the block (representation using 2 bytes per elements)
// and the gap array of the left half-block wrt right half-block (bitvector
// representation), compute the gap array (wrt tail) of the right half-block
// and write to a given file using v-byte encoding.
//
// The whole computation is performed under given ram budget. It is fully
// parallelized and uses asynchronous I/O as much as possible.
//==============================================================================
void compute_right_gap(long left_block_size, long right_block_size,
    gap_array_2n *block_gap, bitvector *bv, std::string out_filename,
    long max_threads, long ram_budget) {
  long block_size = left_block_size + right_block_size;
  long bv_size = block_size + 1;
  long right_gap_size = right_block_size + 1;

  fprintf(stderr, "  Compute gap for right half-block: ");
  long compute_gap_start = utils::wclock();

  //----------------------------------------------------------------------------
  // STEP 1: Preprocess left_block_gap_bv for rank and select queries,
  //         i.e., compute sparse_gap.
  //----------------------------------------------------------------------------
  
  // 1.a
  //
  // Compute chunk size and allocate spase rank.
  long chunk_size = std::min((1L << 20), (bv_size + max_threads - 1) / max_threads);
  long n_chunks = bv_size / chunk_size;  // we exclude the last partial chunk
  long *sparse_rank = (long *)malloc((n_chunks + 1) * sizeof(long));

  // 1.b
  //
  // Compute the values of sparse_rank, that is, the sum of 1-bits inside
  // each chunk. Since there can be more chunks than threads, we split chunks
  // into groups and let each thread handle the group of chunks.
  long chunk_max_group_size = (n_chunks + max_threads - 1) / max_threads;
  long n_chunk_groups = (n_chunks + chunk_max_group_size - 1) / chunk_max_group_size;

  std::thread **threads = new std::thread*[n_chunk_groups];
  for (long t = 0; t < n_chunk_groups; ++t) {
    long chunk_group_beg = t * chunk_max_group_size;
    long chunk_group_end = std::min(chunk_group_beg + chunk_max_group_size, n_chunks);
    threads[t] = new std::thread(xprocess_group_of_chunks, chunk_group_beg,
        chunk_group_end, chunk_size, sparse_rank, bv);
  }

  for (long t = 0; t < n_chunk_groups; ++t) threads[t]->join();
  for (long t = 0; t < n_chunk_groups; ++t) delete threads[t];
  delete[] threads;

  // 1.c
  //
  // Compute cumulative sum of sparse_rank. From now on we can quickly
  // answer rank, select0 and select1 queries on the bitvector.
  for (long i = 0, sum = 0L; i <= n_chunks; ++i) {
    long temp = sparse_rank[i];
    sparse_rank[i] = sum;
    sum += temp;
  }


  //============================================================================
  // STEP 2: compute the values of the right gap array, one range at a time.
  //============================================================================
  long max_range_size = std::max(1L, ram_budget / (3L * (long)sizeof(long)));
  long n_ranges = (right_gap_size + max_range_size - 1) / max_range_size;

  // To ensure that asynchronous I/O is really taking
  // place, we try to make 8 parts.
  if (n_ranges < 8L) {
    max_range_size = (right_gap_size + 7L) / 8L;
    n_ranges = (right_gap_size + max_range_size - 1) / max_range_size;
  }

  long *range_gap = (long *)malloc(max_range_size * sizeof(long));
  unsigned char *active_vbyte_slab = (unsigned char *)malloc(max_range_size * sizeof(long));
  unsigned char *passive_vbyte_slab = (unsigned char *)malloc(max_range_size * sizeof(long));
  long active_vbyte_slab_length;
  long passive_vbyte_slab_length;

  // Used for communication with thread doing asynchronous writes.  
  std::mutex mtx;
  std::condition_variable cv;
  bool avail = false;
  bool finished = false;
  
  // Start the thread doing asynchronius writes.
  std::thread *async_writer = new std::thread(async_write_code, std::ref(passive_vbyte_slab),
    std::ref(passive_vbyte_slab_length), std::ref(mtx), std::ref(cv), std::ref(avail),
    std::ref(finished), out_filename);

  for (long range_id = 0L; range_id < n_ranges; ++range_id) {
    // Compute the range [range_beg..range_end) of values in the right gap
    // (which if indexed [0..right_gap_size)).
    long range_beg = range_id * max_range_size;
    long range_end = std::min(range_beg + max_range_size, right_gap_size);
    long range_size = range_end - range_beg;

    // 2.a
    //
    // Find the section in the bitvector that contains
    // the bits necessary to compute the answer.
    long bv_section_beg = 0L;
    long bv_section_end = 0L;
    if (range_beg > 0)
      bv_section_beg = xcompute_select1_using_sparse_rank(range_beg - 1, chunk_size, n_chunks, sparse_rank, bv) + 1;
    bv_section_end = xcompute_select1_using_sparse_rank(range_end - 1, chunk_size, n_chunks, sparse_rank, bv) + 1;
    long bv_section_size = bv_section_end - bv_section_beg;

    // We split the current bitvector section into
    // equal parts. Each thread handles one part.
    long max_part_size = (bv_section_size + max_threads - 1) / max_threads;
    long n_parts = (bv_section_size + max_part_size - 1) / max_part_size;

    parallel_utils::parallel_fill<long>(range_gap, range_size, 0L, max_threads);

    // Allocate arrays used to store the answers for part boundaries.
    long *res_sum = new long[n_parts];
    long *res_rank = new long[n_parts];

    threads = new std::thread*[n_parts];
    for (long t = 0; t < n_parts; ++t) {
      long part_beg = bv_section_beg + t * max_part_size;
      long part_end = std::min(part_beg + max_part_size, bv_section_end);

      threads[t] = new std::thread(handle_bv_part, part_beg, part_end, range_beg,
          chunk_size, range_gap, sparse_rank, block_gap, bv, std::ref(res_sum[t]),
          std::ref(res_rank[t]));
    }

    for (long t = 0; t < n_parts; ++t) threads[t]->join();
    for (long t = 0; t < n_parts; ++t) delete threads[t];
    delete[] threads;

    // Update range_gap with values computed at part boundaries.
    for (long t = 0; t < n_parts; ++t)
      range_gap[res_rank[t] - range_beg] += res_sum[t];
    delete[] res_sum;
    delete[] res_rank;

    // 2.c
    //
    // Convert the range_gap to the slab of vbyte encoding.
    active_vbyte_slab_length = parallel_utils::convert_array_to_vbyte_slab(range_gap, range_size, active_vbyte_slab, max_threads);


    // 2.d
    //
    // Asynchronously schedule the write of the slab.

    // Wait for the async I/O thread to finish writing.
    std::unique_lock<std::mutex> lk(mtx);
    while (avail == true)
      cv.wait(lk);

    // Set the new passive slab.
    std::swap(active_vbyte_slab, passive_vbyte_slab);
    passive_vbyte_slab_length = active_vbyte_slab_length;

    // Let the I/O thread know that the slab is waiting.
    avail = true;
    lk.unlock();
    cv.notify_one();
  }

  // Let the I/O thread know that we're done.
  std::unique_lock<std::mutex> lk(mtx);
  finished = true;
  lk.unlock();
  cv.notify_one();
  
  // Wait for the thread to actually finish and delete it.
  async_writer->join();
  delete async_writer;


  free(sparse_rank);
  free(range_gap);
  free(active_vbyte_slab);
  free(passive_vbyte_slab);

  long double compute_gap_time = utils::wclock() - compute_gap_start;
  long double compute_gap_speed = (block_size / (1024.L * 1024)) / compute_gap_time;
  fprintf(stderr, "%.2Lf (%.2LfMiB/s)\n", compute_gap_time, compute_gap_speed);
}

#endif  // __COMPUTE_RIGHT_GAP_H_INCLUDED
