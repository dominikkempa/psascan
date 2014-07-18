
//
// TODO:
//   I don't understand why this version is slower than the one in update.h
//   should be exactly as fast, the houskeeping cannot make such a big difference.
//

#ifndef __UPDATE_H_INCLUDED
#define __UPDATE_H_INCLUDED

#include <thread>
#include <mutex>
#include <algorithm>

#include "utils.h"
#include "buffer.h"
#include "gap_array.h"
#include "stream_info.h"

const long n_buckets = 24; // also number of threads
const long buffer_sample_size = 1000;

long threads_created;
long double thread_time_creation;

long double housekeeping_time;
long double counting_time;
long double permuting_time;
long double increasing_time;

template<typename block_offset_type> // XXX: try using block_offset_type for bucket_beg
void parallel_count(buffer<block_offset_type> *b, long block_beg,
    long block_end, long *bucket_lbound, int *block2bucket_contrib,
    int *bucket_id) {
  for (long i = block_beg; i < block_end; ++i) {
    long x = (long)b->m_content[i];

    long id = n_buckets;
    while (bucket_lbound[id] > x) --id;

    bucket_id[i] = id;
    ++block2bucket_contrib[id];
  }
}

// XXX: check if oracle helps
// XXX: oracle[i] vs oracle[i - block_beg]???
template<typename block_offset_type>
void parallel_permute(buffer<block_offset_type> *b, block_offset_type *dest,
    long block_beg, long block_end, int *placement_ptrs, int *bucket_id,
    long *oracle) {
  for (long i = block_beg; i < block_end; ++i) {
    long id = bucket_id[i];
    oracle[i] = placement_ptrs[id]++;
  }
  for (long i = block_beg; i < block_end; ++i) {
    long addr = oracle[i];
    dest[addr] = b->m_content[i];
  }
}

template<typename block_offset_type>
void parallel_increase(block_offset_type *temp, buffered_gap_array *gap,
    long bucket_beg, long bucket_end) {
  for (long i = bucket_beg; i < bucket_end; ++i) {
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

  //============================================================================
  // OVERVIEW
  //============================================================================
  // Each bucket is characterized by two ranges:
  //   * lower and upper bound determine the range of elements that fall into
  //     the bucket. All ranges are disjoint and altogether sup up to
  //     [0..gap->m_length).
  //   * bucket beg and bucket end are the pointers to the buffer, in which
  //     the elements of that buffer will be placed.
  //
  // The computation goes as follows:
  //   1) we guess bucket boundaries (based on the random sample of buffer
  //      elements) so that all buckets will contain roungly equal number of
  //      elements from the buffer (do not confuse with the number of elements
  //      in the range, which is usually much larger)
  //   2) then we compute (in parallel) actual bucket sizes
  //   3) then we are able to compute bucket beg and end
  //   4) we place the buffer elements in their buckets. We permute the
  //      elements out of place (into 'temp' array).
  //   5) finally, we update the gap array, each thread handles one bucket
  //============================================================================


  //----------------------------------------------------------------------------
  // STEP1 1: compute random sample of elements in the buffer
  //----------------------------------------------------------------------------
  long double tt = utils::wclock();
  std::vector<block_offset_type> samples(buffer_sample_size);
  for (long i = 0; i < buffer_sample_size; ++i)
    samples[i] = b->m_content[utils::random_long(0L, b->m_filled - 1)];
  std::sort(samples.begin(), samples.end());
  samples.erase(std::unique(samples.begin(), samples.end()), samples.end());

  //----------------------------------------------------------------------------
  // STEP 2: compute bucket boundaries (only lower bound is sufficient).
  //
  // To simplify parallel code we allow empty buckets.
  // After computation the boundaries of i-th (i = 0,..,n_bucket-1)
  // will be [bucket_lbound[i], bucket_lbound[i + 1]).
  //----------------------------------------------------------------------------
  long *bucket_lbound = new long[n_buckets + 1]; // with sentinel
  std::fill(bucket_lbound, bucket_lbound + n_buckets + 1, gap->m_length);

  long step = (samples.size() + n_buckets - 1) / n_buckets;
  for (size_t i = 1, p = step; p < samples.size(); ++i, p += step)
    bucket_lbound[i] = (samples[p - 1] + samples[p] + 1) / 2;
  bucket_lbound[0] = 0;

  //----------------------------------------------------------------------------
  // STEP 3: divide the buffer into blocks and compute bucket counts in each.
  //
  // Again, we allow empty blocks just to simplify the parallel code.
  // In addition, we compute for each element to which bucket it belongs.
  //----------------------------------------------------------------------------

  // 3.a
  //
  // Allocate block id and bucket counts.
  int *bucket_id = new int[b->m_filled];
  int **block2bucket_contrib = new int*[n_buckets];
  for (long i = 0; i < n_buckets; ++i) {
    block2bucket_contrib[i] = new int[n_buckets];
    std::fill(block2bucket_contrib[i], block2bucket_contrib[i] + n_buckets, 0);
  }

  // 3.b
  //
  // Compute bucket counts and block ids in parallel.
  long block_size = (b->m_filled + n_buckets - 1) / n_buckets;
  std::thread **counters = new std::thread*[n_buckets];
  housekeeping_time += utils::wclock() - tt;
  tt = utils::wclock();
  for (long i = 0, beg = 0, end; i < n_buckets; ++i, beg = end) {
    end = std::min(b->m_filled, beg + block_size);
    counters[i] = new std::thread(parallel_count<block_offset_type>,
        b, beg, end, bucket_lbound, block2bucket_contrib[i], bucket_id);
  }
  for (long i = 0; i < n_buckets; ++i) counters[i]->join();
  counting_time += utils::wclock() - tt;
  tt = utils::wclock();
  for (long i = 0; i < n_buckets; ++i) delete counters[i];
  delete[] counters;
  delete[] bucket_lbound;

  //----------------------------------------------------------------------------
  // STEP 4: place the elements in the buckets.
  //
  // For that we first compute placement pointers and then perform the
  // permutation (with the help of oracle) in parallel.
  //----------------------------------------------------------------------------

  // 4.a
  //
  // Compute bucket size and beg.
  int *bucket_size = new int[n_buckets];
  int *bucket_beg = new int[n_buckets];
  std::fill(bucket_size, bucket_size + n_buckets, 0L);
  for (long i = 0; i < n_buckets; ++i)
    for (long j = 0; j < n_buckets; ++j)
      bucket_size[j] += block2bucket_contrib[i][j];
  for (long t = 0, curbeg = 0; t < n_buckets; curbeg += bucket_size[t++])
    bucket_beg[t] = curbeg;

  // 4.b
  //
  // Allocate placement pointers.
  int **placement_ptrs = new int*[n_buckets];
  for (long i = 0; i < n_buckets; ++i) {
    placement_ptrs[i] = new int[n_buckets];
    for (long j = 0; j < n_buckets; ++j) {
      placement_ptrs[i][j] = bucket_beg[j];
      bucket_beg[j] += block2bucket_contrib[i][j];
    }
    delete[] block2bucket_contrib[i];
  }
  delete[] block2bucket_contrib;

  // 4.c
  //
  // Place the elements in buckets (out of place).
  std::thread **permuters = new std::thread*[n_buckets];
  housekeeping_time += utils::wclock() - tt;
  tt = utils::wclock();
  for (long i = 0, beg = 0, end; i < n_buckets; ++i, beg = end) {
    end = std::min(b->m_filled, beg + block_size);
    permuters[i] = new std::thread(parallel_permute<block_offset_type>,
        b, temp, beg, end, placement_ptrs[i], bucket_id, oracle);
  }
  for (long i = 0; i < n_buckets; ++i) permuters[i]->join();
  permuting_time += utils::wclock() - tt;
  tt = utils::wclock();
  for (long i = 0; i < n_buckets; ++i) {
    delete permuters[i];
    delete[] placement_ptrs[i];
  }
  delete[] placement_ptrs;
  delete[] permuters;
  delete[] bucket_id;

  //----------------------------------------------------------------------------
  // STEP 5: update the gap array in parallel.
  //
  // One thread per bucket (note again, buckets are allowed to be empty).
  //----------------------------------------------------------------------------

  // 5.a
  //
  // Recompute bucket beg.
  for (long t = 0, curbeg = 0; t < n_buckets; curbeg += bucket_size[t++])
    bucket_beg[t] = curbeg;

  // 5.b
  //
  // Increase gap array entries.
  std::thread **increasers = new std::thread*[n_buckets];
  housekeeping_time += utils::wclock() - tt;
  tt = utils::wclock();
  for (long i = 0; i < n_buckets; ++i) {
    increasers[i] = new std::thread(parallel_increase<block_offset_type>,
        temp, gap, bucket_beg[i], bucket_beg[i] + bucket_size[i]);
  }
  for (long i = 0; i < n_buckets; ++i) increasers[i]->join();
  increasing_time += utils::wclock() - tt;
  tt = utils::wclock();
  for (long i = 0; i < n_buckets; ++i) delete increasers[i];
  delete[] increasers;
  delete[] bucket_size;
  delete[] bucket_beg;
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

  fprintf(stderr, "\nhousekeeping time: %.2Lfs\n", housekeeping_time);
  fprintf(stderr, "counting_time: %.2Lfs\n", counting_time);
  fprintf(stderr, "permuting_time: %.2Lfs\n", permuting_time);
  fprintf(stderr, "increasing_time: %.2Lfs\n", increasing_time);
}

#endif // __UPDATE_H_INCLUDED
