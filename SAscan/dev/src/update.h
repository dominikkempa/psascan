#ifndef __UPDATE_H_INCLUDED
#define __UPDATE_H_INCLUDED

#include <thread>
#include <mutex>
#include <algorithm>

#include "utils.h"
#include "buffer.h"
#include "gap_array.h"
#include "stream_info.h"

const long n_updaters = 24;
const long n_counters = 24;

long double housekeeping_time;
long double counting_time;
long double permuting_time;
long double increasing_time;

static const long max_buckets = 1048; // not necessarily a power of 2
long bucket_size_bits;

template<typename block_offset_type>
void parallel_count(buffer<block_offset_type> *b,
    long block_beg, long block_end, int *count) {
  for (long i = block_beg; i < block_end; ++i) {
    long x = (long)b->m_content[i];
    ++count[x >> bucket_size_bits];
  }
}

template<typename block_offset_type>
void parallel_permute(buffer<block_offset_type> *b, block_offset_type *temp,
    long block_beg, long block_end, int *placement_ptrs, long *sbucket_beg, long *oracle) {
  for (long i = block_beg; i < block_end; ++i) {
    long x = (long)b->m_content[i];
    long sblock_id = n_updaters - 1;
    while (sbucket_beg[sblock_id] > x) --sblock_id;
    oracle[i - block_beg] = placement_ptrs[sblock_id]++;
  }

  for (long i = block_beg; i < block_end; ++i) {
    long addr = oracle[i - block_beg];
    temp[addr] = b->m_content[i];
  }
}

template<typename block_offset_type>
void parallel_increase(block_offset_type *temp, buffered_gap_array *gap,
    long sblock_beg, long sblock_end) {
  for (long i = sblock_beg; i < sblock_end; ++i) {
    block_offset_type x = temp[i];
    gap->m_count[x]++;

    // check is gap values wrapped-around
    if (gap->m_count[x] == 0) {
      gap->m_excess_mutex.lock();
      gap->add_excess(x);
      gap->m_excess_mutex.unlock();
    }
  }
}

template<typename block_offset_type>
void update_gap(buffer<block_offset_type> *b, buffered_gap_array *gap, block_offset_type *temp, long *oracle) {
  if (b->m_filled == 0) {
    fprintf(stderr, "Error: trying to do update on empty buffer\n");
    std::exit(EXIT_FAILURE);
  }

  //----------------------------------------------------------------------------
  // STEP 1: compute superbuckets. All elements belonging to a superbucket
  //         make up a superblock. The goal is to choose superbucket sizes
  //         so that superblock are as even as possible.
  //----------------------------------------------------------------------------
  //   a) First, we compute the bucket size. Its size is a power of two
  //      (to allows fast computation of bucket id based in the value) and
  //      there can be at most max_bucket buckets.
  //   b) We divide the buffer into n_counters block of equal size
  //      (block_size).
  //   c) For each block we run a thread counting the number of elements
  //      in that block falling into each bucket.
  //   d) Based on the counts computed in the previous step we compute the
  //      superbucket sizes. Also, we compute the contribution of each
  //      thread from step c into each superbucket. The latter is needed in
  //      STEP 2.
  //----------------------------------------------------------------------------

  // 1.a
  //
  // Compute the smallest bucket_size that results in
  // distributing elements from range [0..gap->m_length)
  // into at most max_buckets buckets.
  long double tt = utils::wclock();
  long bucket_size = 1;
  bucket_size_bits = 0;
  while ((gap->m_length + bucket_size - 1) / bucket_size > max_buckets)
    bucket_size <<= 1, ++bucket_size_bits;
  long n_buckets = (gap->m_length + bucket_size - 1) / bucket_size;


  // 1.b
  //
  // Compute block size.
  long block_size = (b->m_filled + n_counters - 1) / n_counters;


  // 1.c
  //
  // For each block compute the number of elements in each bucket.
  int **block_counts = new int*[n_counters];
  for (long t = 0; t < n_counters; ++t) {
    block_counts[t] = new int[n_buckets];
    std::fill(block_counts[t], block_counts[t] + n_buckets, 0);
  }
  std::thread **counters = new std::thread*[n_counters];
  housekeeping_time += utils::wclock() - tt;
  tt = utils::wclock();
  for (long t = 0, beg = 0, end; t < n_counters; ++t, beg = end) {
    end = std::min(b->m_filled, beg + block_size);
    counters[t] = new std::thread(parallel_count<block_offset_type>,
        b, beg, end, block_counts[t]);
  }
  for (long i = 0; i < n_counters; ++i) counters[i]->join();
  counting_time += utils::wclock() - tt;
  tt = utils::wclock();
  for (long i = 0; i < n_counters; ++i) delete counters[i];
  delete[] counters;


  // 1.d
  //
  // Compute superblock sizes.
  // In addition, we compute the contribution (i.e., the number of elements)
  // of each block into each superblock: block_to_sblock_contrib[i][t]
  // = contribution of block i into superblock t.
  long *sbucket_beg = new long[n_updaters];
  long *sblock_size = new long[n_updaters];

  int *total_count = new int[n_buckets]; // accumulate the counts
  std::fill(total_count, total_count + n_buckets, 0);
  for (long i = 0L; i < n_counters; ++i)
    for (long j = 0L; j < n_buckets; ++j)
      total_count[j] += block_counts[i][j];

  int **block_to_sblock_contrib = new int*[n_counters];
  for (int i = 0; i < n_counters; ++i)
    block_to_sblock_contrib[i] = new int[n_updaters];

  long ideal_sblock_size = (b->m_filled + n_updaters - 1) / n_updaters;
  long bucket_id_beg = 0;
  for (long t = 0; t < n_updaters; ++t) {
    long bucket_id_end = bucket_id_beg;

    long size = 0L;
    while (bucket_id_end < n_buckets && size < ideal_sblock_size)
      size += total_count[bucket_id_end++];

    sblock_size[t] = size;
    sbucket_beg[t] = std::min(gap->m_length, bucket_size * bucket_id_beg);

    // Compute the contribution of each block to this superblock.
    for (long i = 0; i < n_counters; ++i) {
      block_to_sblock_contrib[i][t] = 0;
      for (long id = bucket_id_beg; id < bucket_id_end; ++id)
        block_to_sblock_contrib[i][t] += block_counts[i][id];
    }

    bucket_id_beg = bucket_id_end;
  }
  delete[] total_count;
  for (long i = 0; i < n_counters; ++i)
    delete[] block_counts[i];
  delete[] block_counts;

  //----------------------------------------------------------------------------
  // STEP 2: permute elements from the buffer into superblocks. Superblock
  //         is a contugous range of elements in the buffer containing all
  //         elements from one superbucket (in any order).
  //----------------------------------------------------------------------------
  //    a) compute superblock starting positions. This is easy compute if we
  //       know what was the contribution of each block into each superbucket.
  //       we computed that in step 1.d
  //    b) Permute elements in the buffer so that elements from the same
  //       superbucket are in the contiguous range. We do this in parallel,
  //       one thread per block.
  //----------------------------------------------------------------------------

  // 2.a
  //
  // For each block and superblock compute the pointer where
  // the elements from that block should be placed.
  long *sblock_beg = new long[n_updaters];
  for (long t = 0, curbeg = 0; t < n_updaters; curbeg += sblock_size[t++])
    sblock_beg[t] = curbeg;

  int **placement_pointers = new int*[n_counters];
  for (long i = 0; i < n_counters; ++i) {
    placement_pointers[i] = new int[n_updaters];
    for (long t = 0; t < n_updaters; ++t) {
      placement_pointers[i][t] = sblock_beg[t];
      sblock_beg[t] += block_to_sblock_contrib[i][t];
    }
  }

  for (long i = 0; i < n_counters; ++i)
    delete[] block_to_sblock_contrib[i];
  delete[] block_to_sblock_contrib;

  // 2.b
  //
  // Create and run threads performing the placement.
  std::thread **permuters = new std::thread*[n_counters];
  housekeeping_time += utils::wclock() - tt;
  tt = utils::wclock();
  for (long i = 0, beg = 0, end; i < n_counters; ++i, beg = end) {
    end = std::min(b->m_filled, beg + block_size);
    permuters[i] = new std::thread(parallel_permute<block_offset_type>,
        b, temp, beg, end, placement_pointers[i], sbucket_beg, oracle + beg);
  }
  for (long i = 0; i < n_counters; ++i) permuters[i]->join();
  permuting_time += utils::wclock() - tt;
  tt = utils::wclock();
  for (long i = 0; i < n_counters; ++i) {
    delete permuters[i];
    delete[] placement_pointers[i];
  }
  delete[] permuters;
  delete[] placement_pointers;
  delete[] sbucket_beg;

  //----------------------------------------------------------------------------
  // STEP 3: increate the gap array values. Create one thread for every
  // superblock and perform the updates.
  //----------------------------------------------------------------------------

  for (long t = 0, curbeg = 0; t < n_updaters; curbeg += sblock_size[t++])
    sblock_beg[t] = curbeg;

  std::thread **increasers = new std::thread*[n_updaters];
  housekeeping_time += utils::wclock() - tt;
  tt = utils::wclock();
  for (long t = 0; t < n_updaters; ++t) {
    increasers[t] = new std::thread(parallel_increase<block_offset_type>,
        temp, gap, sblock_beg[t], sblock_beg[t] + sblock_size[t]);
  }
  for (long t = 0; t < n_updaters; ++t) increasers[t]->join();
  increasing_time += utils::wclock() - tt;
  tt = utils::wclock();
  for (long t = 0; t < n_updaters; ++t) delete increasers[t];
  delete[] increasers;

  delete[] sblock_size;
  delete[] sblock_beg;
  housekeeping_time += utils::wclock() - tt;
}

template<typename block_offset_type>
void gap_updater(buffer_poll<block_offset_type> *full_buffers,
    buffer_poll<block_offset_type> *empty_buffers,
    buffered_gap_array *gap, long stream_buf_size) {

  housekeeping_time = 0.L;
  counting_time = 0.L;
  permuting_time = 0.L;
  increasing_time = 0.L;

  long max_buffer_elems = stream_buf_size / sizeof(block_offset_type);
  block_offset_type *temp = new block_offset_type[max_buffer_elems];
  long *oracle = new long[max_buffer_elems];

  while (true) {
    // Get a buffer from the poll of full buffers.
    std::unique_lock<std::mutex> lk(full_buffers->m_mutex);
    while (!full_buffers->available() && !full_buffers->finished())
      full_buffers->m_cv.wait(lk);

    if (!full_buffers->available() && full_buffers->finished()) {
      // All workers finished. We're exiting, but before, we're letting
      // other updating threads know that they also should check for exit.
      lk.unlock();
      full_buffers->m_cv.notify_one();
      break;
    }

    buffer<block_offset_type> *b = full_buffers->get();
    lk.unlock();
    full_buffers->m_cv.notify_one(); // let others know they should try

    // Process buffer.
    update_gap<block_offset_type>(b, gap, temp, oracle);

    // Add the buffer to the poll of empty buffers and notify waiting thread.
    std::unique_lock<std::mutex> lk2(empty_buffers->m_mutex);
    empty_buffers->add(b);
    lk2.unlock();
    empty_buffers->m_cv.notify_one();
  }
  
  delete[] temp;
  delete[] oracle;

  fprintf(stderr, "\nhousekeeping_time: %.2Lfs\n", housekeeping_time);
  fprintf(stderr, "counting_time: %.2Lfs\n", counting_time);
  fprintf(stderr, "permuting_time: %.2Lfs\n", permuting_time);
  fprintf(stderr, "increasing_time: %.2Lfs\n", increasing_time);
}

#endif // __UPDATE_H_INCLUDED
