#ifndef __COMPUTE_GAP_H
#define __COMPUTE_GAP_H

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <string>
#include <thread>
#include <algorithm>
#include <vector>

#include "utils.h"
#include "rank.h"
#include "gap_array.h"
#include "io_streamer.h"
#include "buffer.h"
#include "stream.h"
#include "update.h"
#include "stream_info.h"
#include "multifile.h"


//==============================================================================
// Compute the gap for an arbitrary range of suffixes of tail. This version is
// more general, and can be used also when processing half-blocks.
//==============================================================================
template<typename block_offset_type>
void compute_gap(rank4n<> *rank, buffered_gap_array *gap,
    long tail_begin, long tail_end, long text_length, long max_threads,
    long block_isa0, long stream_bufsize, unsigned char block_last_symbol,
    std::vector<long> initial_ranks, std::string text_filename,
    multifile *tail_gt_begin_rev, multifile *newtail_gt_begin_rev) {

  long tail_length = tail_end - tail_begin;
  long stream_max_block_size = (tail_length + max_threads - 1) / max_threads;
  long n_threads = (tail_length + stream_max_block_size - 1) / stream_max_block_size;

  // Get symbol counts of a block and turn into exclusive partial sum.
  long *count = new long[256];
  std::copy(rank->m_count, rank->m_count + 256, count);
  ++count[block_last_symbol];
  --count[0];
  for (long j = 0, s = 0, t; j < 256; ++j) {
    t = count[j];
    count[j] = s;
    s += t;
  }

  // Allocate the gap array, do the streaming and store gap to disk.
  fprintf(stderr, "    Stream:");
  long double stream_start = utils::wclock();

  // Allocate buffers.
  long n_stream_buffers = 2 * max_threads;
  buffer<block_offset_type> **buffers = new buffer<block_offset_type>*[n_stream_buffers];
  for (long i = 0L; i < n_stream_buffers; ++i)
    buffers[i] = new buffer<block_offset_type>(stream_bufsize, max_threads);

  // Create poll of empty and full buffers.
  buffer_poll<block_offset_type> *empty_buffers = new buffer_poll<block_offset_type>();
  buffer_poll<block_offset_type> *full_buffers = new buffer_poll<block_offset_type>(n_threads);

  // Add empty buffers to empty poll.
  for (long i = 0L; i < n_stream_buffers; ++i)
    empty_buffers->add(buffers[i]);

  // Start workers.
  stream_info info(n_threads, tail_length);
  std::thread **streamers = new std::thread*[n_threads];
  std::vector<std::string> gt_filenames(n_threads);

  for (long t = 0L; t < n_threads; ++t) {
    long stream_block_beg = tail_begin + t * stream_max_block_size;
    long stream_block_end = std::min(stream_block_beg + stream_max_block_size, tail_end);

    gt_filenames[t] = text_filename + ".gt_tail." + utils::random_string_hash();
    newtail_gt_begin_rev->add_file(text_length - stream_block_end, text_length - stream_block_beg, gt_filenames[t]);

    streamers[t] = new std::thread(parallel_stream<block_offset_type>, full_buffers, empty_buffers, stream_block_beg,
        stream_block_end, initial_ranks[t], count, block_isa0, rank, block_last_symbol, text_filename, text_length,
        std::ref(gt_filenames[t]), &info, t, gap->m_length, stream_bufsize, tail_gt_begin_rev, max_threads);
  }

  // Start updaters.
  std::thread *updater = new std::thread(gap_updater<block_offset_type>,
        full_buffers, empty_buffers, gap, max_threads);

  // Wait for all threads to finish.
  for (long i = 0L; i < n_threads; ++i) streamers[i]->join();
  updater->join();

  // Clean up.
  for (long i = 0L; i < n_threads; ++i) delete streamers[i];
  for (long i = 0L; i < n_stream_buffers; ++i) delete buffers[i];
  delete updater;
  delete[] streamers;
  delete[] buffers;
  delete empty_buffers;
  delete full_buffers;
  delete[] count;

  long double stream_time = utils::wclock() - stream_start;
  long double speed = (tail_length / (1024.L * 1024)) / stream_time;
  fprintf(stderr,"\r    Stream: 100.0%%. Time: %.2Lf. Threads: %ld. "
      "Speed: %.2LfMiB/s (avg), %.2LfMiB/s (total)\n",
      stream_time, info.m_thread_count, speed / n_threads, speed);
}


#endif  // __COMPUTE_GAP_H
