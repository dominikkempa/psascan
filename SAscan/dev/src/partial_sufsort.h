#ifndef __PARTIAL_SUFSORT_H_INCLUDED
#define __PARTIAL_SUFSORT_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <iostream>
#include <queue>
#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <algorithm>
#include <vector>

#include "divsufsort.h"
#include "divsufsort64.h"
#include "utils.h"
#include "rank.h"
#include "srank.h"
#include "gap_array.h"
#include "merge.h"
#include "bitvector.h"
#include "stream.h"
#include "sascan.h"
#include "settings.h"
#include "smaller_suffixes.h"
#include "radixsort.h"

// Stores information about a contigous subsequence of bits
// of gt bitvector produced by a single thread. The bits are
// stored in file m_fname and their range is [m_start..m_end).
struct gt_substring_info {
  gt_substring_info() {}
  gt_substring_info(long start, long end, std::string fname)
    : m_start(start),
      m_end(end),
      m_fname(fname) {}

  long m_start;
  long m_end;
  std::string m_fname;
};

template<typename T>
struct buffer {  
  buffer(long size_bytes)
      : m_filled(0L),
        m_size(size_bytes / sizeof(T)) {
    m_content = new T[m_size];
  }
  
  ~buffer() {
    delete[] m_content;
  }

  long m_filled, m_size;
  T *m_content;
};

// Same class for the poll of empty and full buffers.  
template<typename T>
struct buffer_poll {
  buffer_poll(long worker_threads = 0L) {
    m_worker_threads = worker_threads; // unused for the poll of empty buffers.
    m_worker_threads_finished = 0L;
  }
  
  void add(buffer<T> *b) {
    m_queue.push(b);
  }
  
  bool available() const {
    return m_queue.size() > 0;
  }

  buffer<T> *get() {
    if (m_queue.empty()) {
      fprintf(stderr, "Error: requesting a buffer from empty poll!\n");
      std::exit(EXIT_FAILURE);
    }

    buffer<T> *ret = m_queue.front();
    m_queue.pop();

    return ret;
  }

  bool finished() const {
    return m_worker_threads_finished == m_worker_threads;
  }

  void increment_finished_workers() {
    ++m_worker_threads_finished;
  }

  std::condition_variable m_cv;
  std::mutex m_mutex;

private:
  long m_worker_threads; // used to detect that worker threads are done
  long m_worker_threads_finished;

  std::queue<buffer<T>* > m_queue;
};

long n_streamers;
long n_updaters;
long stream_buffer_size;
long n_stream_buffers;
long max_gap_sections;

std::mutex stdout_mutex;

//-----------------------------------------------------------------------------
// On Platform-L 6 is a good value (better than 4 and 8).
// Bigger values make the program run slower.
//-----------------------------------------------------------------------------

template<typename output_type> void SAscan(std::string input_filename, long ram_use);
template<typename output_type> distributed_file<output_type> *partial_SAscan(std::string input_filename,
    bool compute_bwt, long ram_use, unsigned char **BWT, std::string text_filename, long text_offset);

//=============================================================================
// Compute partial SA of B[0..block_size) and store on disk.
// If block_id != n_block also compute the BWT of B.
//
// INVARIANT: on entry to the function it holds: 5.2 * block_size <= ram_use
//=============================================================================
template<typename block_offset_type>
distributed_file<block_offset_type> *compute_partial_sa_and_bwt(
    unsigned char *B,
    long block_size,
    long block_id,
    long ram_use, // XXX is this always the same throughout the algorithm?
    std::string text_fname,
    std::string sa_fname,
    bool compute_bwt,
    bitvector *gt_eof_bv,
    unsigned char **BWT,
    long block_offset) {

  distributed_file<block_offset_type> *result = new distributed_file<block_offset_type>(sa_fname.c_str(),
      std::max((long)sizeof(block_offset_type), ram_use / 10L));
  fprintf(stderr, "  compute-bwt = %s\n", compute_bwt ? "TRUE" : "FALSE");

  if (block_size <= MAX_32BIT_DIVSUFSORT_LENGTH) {
    // Easy case, just use use 32-bit divsufsort.
    fprintf(stderr, "  Computing partial SA (divsufsort): ");
    long double sa_start = utils::wclock();
    int *SA = new int[block_size];
    divsufsort(B, SA, (int)block_size);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - sa_start);

    fprintf(stderr, "  Writing partial SA to disk (using %lu-byte ints): ", sizeof(block_offset_type));
    long double writing_sa_start = utils::wclock();
    result->initialize_writing(4 << 20);
    for (long i = 0; i < block_size; ++i)
      result->write((block_offset_type)SA[i]);
    result->finish_writing();
    fprintf(stderr, "%.2Lf\n", utils::wclock() - writing_sa_start);

    if (compute_bwt) {
      // Remap symbols of B back to original.
      fprintf(stderr, "  Re-remapping B: ");
      long double reremap_start = utils::wclock();
      for (long j = 0; j < block_size; ++j) B[j] -= gt_eof_bv->get(j);
      fprintf(stderr, "%.2Lf\n", utils::wclock() - reremap_start);

      // Compute BWT
      fprintf(stderr, "  Compute BWT: ");
      long double bwt_start = utils::wclock();
      unsigned char *tmp = (unsigned char *)SA;
      for (long j = 0, jj = 0; j < block_size; ++j)
        if (SA[j]) tmp[jj++] = B[SA[j] - 1];
      std::copy(tmp, tmp + block_size - 1, B);
      *BWT = B;
      fprintf(stderr, "%.2Lf\n", utils::wclock() - bwt_start);
    } else delete[] B;
    
    delete[] SA;
    delete gt_eof_bv;
  } else if (9L * block_size <= ram_use) {
    // Easy case: block_size >= 2GiB but enough RAM to use divsufsort64.
    fprintf(stderr, "  Computing partial SA (divsufsort64): ");
    long double sa_start = utils::wclock();
    long *SA = new long[block_size];
    divsufsort64(B, SA, block_size);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - sa_start);
    
    fprintf(stderr, "  Writing partial SA to disk (using %lu-byte ints): ", sizeof(block_offset_type));
    long double writing_sa_start = utils::wclock();
    result->initialize_writing(4 << 20);
    for (long i = 0; i < block_size; ++i)
      result->write((block_offset_type)SA[i]);
    result->finish_writing();
    fprintf(stderr, "%.2Lf\n", utils::wclock() - writing_sa_start);

    if (compute_bwt) {
      // Remap symbols of B back to original.
      fprintf(stderr, "  Re-remapping B: ");
      long double reremap_start = utils::wclock();
      for (long j = 0; j < block_size; ++j) B[j] -= gt_eof_bv->get(j);
      fprintf(stderr, "%.2Lf\n", utils::wclock() - reremap_start);

      // Compute BWT
      fprintf(stderr, "  Compute BWT: ");
      long double bwt_start = utils::wclock();
      unsigned char *tmp = (unsigned char *)SA;
      for (long j = 0, jj = 0; j < block_size; ++j)
        if (SA[j]) tmp[jj++] = B[SA[j] - 1];
      std::copy(tmp, tmp + block_size - 1, B);
      *BWT = B;
      fprintf(stderr, "%.2Lf\n", utils::wclock() - bwt_start);
    } else  delete[] B;
    
    delete gt_eof_bv;
    delete[] SA;
  } else {
    //-------------------------------------------------------------------------
    // (block_size >= 2GiB and 9*block_size > ram_use) => use recursion.
    // To save I/O from recursion we also get the BWT, which is obtained
    // during the storage of the resulting partial SA to disk.
    //-------------------------------------------------------------------------
    fprintf(stderr, "  Recursively computing partial SA:\n");
    long double rec_partial_sa_start = utils::wclock();

    // Save the remapped block to temp file.
    std::string B_fname = text_fname + ".block" + utils::intToStr(block_id);
    utils::write_objects_to_file(B, block_size, B_fname);

    // Free all memory.
    delete gt_eof_bv;
    delete[] B;

    result = partial_SAscan<block_offset_type>(B_fname, compute_bwt, ram_use, BWT, text_fname, block_offset);

    utils::file_delete(B_fname);
    fprintf(stderr, "  Recursively computing partial SA: %.2Lf\n",
        utils::wclock() - rec_partial_sa_start);
  }

  return result;
}

//=============================================================================
// Used to store progress information for different threads during streaming.
//=============================================================================
struct stream_info {
  stream_info(long thread_count, long tostream)
    : m_update_count(0L),
      m_thread_count(thread_count),
      m_tostream(tostream) {
    m_streamed = new long[thread_count];
    std::fill(m_streamed, m_streamed + thread_count, 0L);

    m_idle_update = new long double[thread_count];
    m_idle_work  = new long double[thread_count];
    std::fill(m_idle_update, m_idle_update + thread_count, 0.L);
    std::fill(m_idle_work, m_idle_work + thread_count, 0.L);

    m_timestamp = utils::wclock();
  }

  ~stream_info() {
    delete[] m_streamed;
    delete[] m_idle_work;
    delete[] m_idle_update;
  }

  long m_update_count;     // number of updates
  long m_thread_count;     // number of threads
  long m_tostream;         // total text length to stream
  long double m_timestamp; // when the streaming started
  long *m_streamed;        // how many bytes streamed by each thread
  long double *m_idle_update;
  long double *m_idle_work;
};

std::mutex stream_info_mutex;

template<typename block_offset_type>
void parallel_stream(
    buffer_poll<block_offset_type> *full_buffers,
    buffer_poll<block_offset_type> *empty_buffers,
    long j_start,
    long j_end,
    block_offset_type i,
    long *count,
    block_offset_type whole_suffix_rank,
    context_rank_4n *rank,
    unsigned char last,
    std::string text_filename,
    long length,
    std::string &tail_gt_filename,
    stream_info *info,
    int thread_id) {
  gt_accessor *gt_in = new gt_accessor(text_filename + std::string(".gt")); // 1MiB buffer
  bool next_gt = (j_start + 1 == length) ? 0 : (*gt_in)[length - j_start - 2];
  backward_skip_stream_reader<unsigned char> *text_streamer
    = new backward_skip_stream_reader<unsigned char>(text_filename, length - 1 - j_start, 1 << 20); // 1MiB buffer

  tail_gt_filename = text_filename + std::string(".gt_tail.") + utils::random_string_hash();
  bit_stream_writer *gt_out = new bit_stream_writer(tail_gt_filename); // 1MiB buffer

  long j = j_start, dbg = 0L;
  while (j > j_end) {
    if (dbg > (1 << 26)) {
      stream_info_mutex.lock();
      info->m_streamed[thread_id] = j_start - j;
      info->m_update_count += 1;
      if (info->m_update_count == info->m_thread_count) {
        info->m_update_count = 0L;
        long double elapsed = utils::wclock() - info->m_timestamp;
        long total_streamed = 0L;

        for (long t = 0; t < info->m_thread_count; ++t)
          total_streamed += info->m_streamed[t];
        long double speed = (total_streamed / (1024.L * 1024)) / elapsed;

        stdout_mutex.lock();
        fprintf(stderr, "\r  [PARALLEL]Stream: %.2Lf%%. Time: %.2Lf. Threads: %ld. "
            "Speed: %.2LfMiB/s (avg), %.2LfMiB/s (total)",
            (total_streamed * 100.L) / info->m_tostream, elapsed,
            info->m_thread_count, speed / info->m_thread_count, speed);
        stdout_mutex.unlock();
      }
      stream_info_mutex.unlock();
      dbg = 0L;
    }

    // Get a buffer from the poll of empty buffers.
    std::unique_lock<std::mutex> lk(empty_buffers->m_mutex);
    while (!empty_buffers->available())
      empty_buffers->m_cv.wait(lk);
    buffer<block_offset_type> *b = empty_buffers->get();
    lk.unlock();
    empty_buffers->m_cv.notify_one(); // let others know they should re-check

    // Process buffer -- fill with gap values.
    long left = j - j_end;
    b->m_filled = std::min(left, b->m_size);
    dbg += b->m_filled;
    for (long t = 0L; t < b->m_filled; ++t, --j) {
      unsigned char c = text_streamer->read();
      i = (block_offset_type)(count[c] + rank->rank((long)(i - (i > whole_suffix_rank)), c));
      if (c == last && next_gt) ++i;
      gt_out->write(i > whole_suffix_rank);
      next_gt = (*gt_in)[length - j - 1];
      b->m_content[t] = i;
    }

    // Add the buffer to the poll of full buffers and notify waiting thread.
    std::unique_lock<std::mutex> lk2(full_buffers->m_mutex);
    full_buffers->add(b);
    lk2.unlock();
    full_buffers->m_cv.notify_one();    
  }

  delete text_streamer;
  delete gt_in;
  delete gt_out;
  
  // Report that another worker thread has finished.
  std::unique_lock<std::mutex> lk(full_buffers->m_mutex);
  full_buffers->increment_finished_workers();
  lk.unlock();
  
  // Notify waiting update threads in case no more buffers
  // are going to be produces by worker threads.
  full_buffers->m_cv.notify_one();
}

//=============================================================================
// Parallel streaming
//-----------------------------------------------------------------------------
// Currently each thread is using about 11MiB of RAM.
//=============================================================================
template<typename block_offset_type>
void gap_update(buffer_poll<block_offset_type> *full_buffers,
    buffer_poll<block_offset_type> *empty_buffers,
    buffered_gap_array *gap) {
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
    gap->increment(b->m_content, b->m_filled);

    // Add the buffer to the poll of empty buffers and notify waiting thread.
    std::unique_lock<std::mutex> lk2(empty_buffers->m_mutex);
    empty_buffers->add(b);
    lk2.unlock();
    empty_buffers->m_cv.notify_one();
  }
}

//=============================================================================
// Compute partial SAs and gap arrays and write to disk.
// Return the array of handlers to distributed files as a result.
//=============================================================================
template<typename block_offset_type>
distributed_file<block_offset_type> **partial_sufsort(std::string filename, long length, long max_block_size, long ram_use) {
  long n_block = (length + max_block_size - 1) / max_block_size;
  long block_id = n_block - 1, prev_end = length;

  distributed_file<block_offset_type> **distrib_files = new distributed_file<block_offset_type>*[n_block];

  while (block_id >= 0) {
    long beg = max_block_size * block_id;
    long end = std::min(length, beg + max_block_size);
    long block_size = end - beg; // B = text[beg..end), current block
    bool need_streaming = (block_id + 1 != n_block);

    fprintf(stderr, "Processing block %ld/%ld [%ld..%ld):\n",
      n_block - block_id, n_block, beg, end);
    fprintf(stderr, "  need_streaming = %s\n", need_streaming ? "TRUE" : "FALSE");
    fprintf(stderr, "  block_size = %ld (%.2LfMiB)\n", block_size,
        (long double)block_size / (1 << 20));

    //--------------------------------------------------------------------------
    // Invariant: if the current block [beg..end) is not the last block of
    // text, there is a file called filename.gt on disk containing bitvector of
    // size length - end defined as follows:
    //
    // gt[i] == 1
    //
    //   iff
    //
    // the suffix of text of length i + 1 is lexicographically greater
    // than the current tail of the text, text[end..length).
    //--------------------------------------------------------------------------
    // The bitvector is required during the streaming and other computations,
    // e.g. during the computation of gt_eof bitvector required for renaming
    // the block.
    //
    // Typically the bitvector is accessed via the 'gt_accessor' class (see
    // smaller_suffixes.h). This class implements random access via
    // operator []. Access is efficient (essentially streaming) if its
    // sequential.
    //
    // Example:
    //
    // gt_accessor *gt = new gt_accessor(filename + ".gt");
    // printf("%d", gt[5]);
    //--------------------------------------------------------------------------

    // 1. Read current block.
    fprintf(stderr, "  Reading block: ");
    long double read_start = utils::wclock();
    unsigned char *B = new unsigned char[block_size];
    utils::read_block(filename, beg, block_size, B);
    unsigned char last = B[block_size - 1];
    fprintf(stderr, "%.2Lf\n", utils::wclock() - read_start);

    int starting_positions = std::min(n_streamers, length - end);
    std::vector<gt_substring_info> gt_info(starting_positions);
    std::vector<long> initial_rank(starting_positions);

    long count[256] = {0};
    bitvector *gt_eof_bv = NULL;
    if (need_streaming) {
      // Compute initial ranks for streaming.
      //-----------------------------------------------------------------------
      fprintf(stderr, "  [PARALLEL]Computing initial ranks: ");
      long double initial_ranks_start = utils::wclock();

      // Evenly distribute the starting positions across the range [end..length).
      //
      // Thread t streams symbols of text from gt_info[t].m_start to
      // gt_info[t].m_end - 1 (note that m_end < m_start since streaming goes
      // backwards).
      long parallel_block = (length - end) / starting_positions;
      long prev_start = end - 1;
      for (int t = 0; t < starting_positions; ++t) {
        gt_info[t].m_end = prev_start;
        gt_info[t].m_start = std::min(prev_start + parallel_block, length - 1);
        prev_start = gt_info[t].m_start;
      }
      gt_info[starting_positions - 1].m_start = length - 1;

      // Compute initial ranks for all threads in parallel.
      std::thread **threads = new std::thread*[starting_positions];
      for (int t = 0; t < starting_positions; ++t) {
        threads[t] = new std::thread(parallel_smaller_suffixes,
            B, block_size, filename, gt_info[t].m_start + 1, std::ref(initial_rank[t]));
      }
      for (int t = 0; t < starting_positions; ++t) threads[t]->join();
      for (int t = 0; t < starting_positions; ++t) delete threads[t];
      delete[] threads;
      fprintf(stderr, "%.2Lf\n", utils::wclock() - initial_ranks_start);
      //-----------------------------------------------------------------------
 
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
      utils::read_block(filename, end - 1, ext_prev_block_size, extprevB);
      fprintf(stderr, "%.2Lf\n", utils::wclock() - prev_block_reading_start);

      // 2c. Compute gt_eof.
      fprintf(stderr, "  Compute gt_eof_bv: ");
      long double gt_eof_start = utils::wclock();
      gt_eof_bv = new bitvector(block_size);
      gt_accessor *gt = new gt_accessor(filename + ".gt");
      long gt_length = length - end;
      compute_gt_eof_bv(extprevB, ext_prev_block_size, B, block_size, *gt, gt_length, gt_eof_bv);
      delete gt;
      delete[] extprevB;
      fprintf(stderr, "%.2Lf\n", utils::wclock() - gt_eof_start);

      // FIXME
      bool used = false;
      for (long j = 0; j < block_size; ++j)
        if (B[j] == 255) used = true;
      if (used) {
        printf("\nError: Input text cannot contain symbol 255.\n");
        std::exit(EXIT_FAILURE);
      }

      // 2d. Remap symbols of B.
      fprintf(stderr, "  Remapping B: ");
      long double remap_start = utils::wclock();
      for (long j = 0; j < block_size; ++j) B[j] += gt_eof_bv->get(j);
      fprintf(stderr, "%.2Lf\n", utils::wclock() - remap_start);
    }

    // 3. Compute the head of the new gt bitvector.
    long whole_suffix_rank = 0;
    fprintf(stderr, "  Compute new_gt_head_bv: ");
    long double new_gt_head_bv_start = utils::wclock();
    bitvector *new_gt_head_bv = new bitvector(block_size);
    whole_suffix_rank = compute_new_gt_head_bv(B, block_size, new_gt_head_bv);
    new_gt_head_bv->save(filename + ".gt_head");
    delete new_gt_head_bv;
    fprintf(stderr, "%.2Lf\n", utils::wclock() - new_gt_head_bv_start);

    // 4. Compute/save partial SA. Compute BWT if it's not the last block.
    unsigned char *BWT = NULL;
    std::string sa_fname = filename + ".partial_sa." + utils::intToStr(block_id);

    distrib_files[block_id] =
      compute_partial_sa_and_bwt<block_offset_type>(B, block_size, block_id,
          ram_use, filename, sa_fname, need_streaming, gt_eof_bv, &BWT, beg);

    if (need_streaming) {
      // 5a. Build the rank support for BWT.
      fprintf(stderr, "  Building the rank data structure: ");
      long double building_rank_start = utils::wclock();
      context_rank_4n *rank = new context_rank_4n(BWT, block_size - 1);
      delete[] BWT;
      fprintf(stderr, "%.2Lf\n", utils::wclock() - building_rank_start);

      // 5b. Allocate the gap array, do the streaming and store gap to disk.
      fprintf(stderr, "  [PARALLEL]Stream:");
      long double stream_start = utils::wclock();

      long gap_sections = std::min(block_size + 1, max_gap_sections);
      buffered_gap_array *gap = new buffered_gap_array(block_size + 1, gap_sections);

      // Allocate buffers.
      buffer<block_offset_type> **buffers = new buffer<block_offset_type>*[n_stream_buffers];
      for (long i = 0L; i < n_stream_buffers; ++i)
        buffers[i] = new buffer<block_offset_type>(stream_buffer_size);

      // Create poll of empty and full buffers.
      buffer_poll<block_offset_type> *empty_buffers = new buffer_poll<block_offset_type>();
      buffer_poll<block_offset_type> *full_buffers = new buffer_poll<block_offset_type>(starting_positions);

      // Add empty buffers to empty poll.
      for (long i = 0L; i < n_stream_buffers; ++i)
        empty_buffers->add(buffers[i]);
      
      // Start workers.
      stream_info info(starting_positions, length - end);
      std::thread **streamers = new std::thread*[starting_positions];
      for (long t = 0L; t < starting_positions; ++t) {
        long j_start = gt_info[t].m_start;
        long j_end = gt_info[t].m_end;

        streamers[t] = new std::thread(parallel_stream<block_offset_type>,
            full_buffers, empty_buffers, j_start, j_end, initial_rank[t],
            count, whole_suffix_rank, rank, last, filename, length,
            std::ref(gt_info[t].m_fname), &info, t);
      }

      // Start updaters.
      std::thread **updaters = new std::thread*[n_updaters];
      for (long i = 0L; i < n_updaters; ++i) {
        updaters[i] = new std::thread(gap_update<block_offset_type>,
            full_buffers, empty_buffers, gap);
      }

      // Wait for all threads to finish.        
      for (long i = 0L; i < starting_positions; ++i) streamers[i]->join();
      for (long i = 0L; i < n_updaters; ++i) updaters[i]->join();

      // Clean up.
      for (long i = 0L; i < starting_positions; ++i) delete streamers[i];
      for (long i = 0L; i < n_updaters; ++i) delete updaters[i];
      for (long i = 0L; i < n_stream_buffers; ++i) delete buffers[i];
      delete[] updaters;
      delete[] streamers;
      delete[] buffers;
      delete empty_buffers;
      delete full_buffers;
      delete rank;
      utils::file_delete(filename + ".gt");

      long double stream_time = utils::wclock() - stream_start;
      long double speed = ((length - end) / (1024.L * 1024)) / stream_time;
      fprintf(stderr,"\r  [PARALLEL]Stream: 100.0%%. Time: %.2Lf. Threads: %ld. "
          "Speed: %.2LfMiB/s (avg), %.2LfMiB/s (total)\n",
          stream_time, info.m_thread_count, speed / starting_positions, speed);

      //-----------------------------------------------------------------------
      // for (long t = 0L; t < starting_points; ++t)
      //   fprintf(stderr, "  IDLE(%ld) = %.3Lfsec (work), %.3Lfsec (update)\n",
      //       t, info.m_idle_work[t], info.m_idle_update[t]);
      //-----------------------------------------------------------------------

      // 5c. Save gap to file.
      gap->save_to_file(filename + ".gap." + utils::intToStr(block_id));
      delete gap;
    }

    // Concatenate all bitvectors into one.
    //-------------------------------------------------------------------------
    fprintf(stderr, "  Concatenating gt bitvectors: ");
    long double gt_concat_start = utils::wclock();
    bit_stream_writer *gt = new bit_stream_writer(filename + ".gt");
    if (need_streaming) {
      for (long jj = starting_positions - 1; jj >= 0; --jj) {
        long j_start = gt_info[jj].m_start;
        long j_end = gt_info[jj].m_end;
        bit_stream_reader *gt_tail = new bit_stream_reader(gt_info[jj].m_fname);
        for (long tt = 0; tt < j_start - j_end; ++tt) gt->write(gt_tail->read());
        delete gt_tail;
        utils::file_delete(gt_info[jj].m_fname);
      }
    }
    bit_stream_reader *gt_head = new bit_stream_reader(filename + ".gt_head");
    for (long tt = end - 1; tt >= beg; --tt) gt->write(gt_head->read());
    delete gt_head;
    utils::file_delete(filename + ".gt_head");
    delete gt;
    fprintf(stderr, "%.2Lf\n", utils::wclock() - gt_concat_start);
    //--------------------------------------------------------------------------

    prev_end = end;
    end = beg;
    --block_id;
  }

  if (utils::file_exists(filename + ".gt"))
    utils::file_delete(filename + ".gt");

  return distrib_files;
}

#endif // __PARTIAL_SUFSORT_H_INCLUDED
