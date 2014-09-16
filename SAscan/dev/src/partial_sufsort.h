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
#include "io_streamer.h"
#include "sascan.h"
#include "settings.h"
#include "smaller_suffixes.h"
#include "radixsort.h"
#include "buffer.h"
#include "update.h"
#include "stream_info.h"
#include "aux_parallel.h"
#include "multifile_bitvector.h"
#include "inmem_sascan/inmem_sascan.h"



template<typename saidx_t>
void compute_initial_ranks(unsigned char *block, long block_beg,
    long block_end, long text_length, saidx_t *block_partial_sa,
    std::string text_filename, std::vector<long> &result, multifile *tail_gt_begin_reversed,
    long max_threads, long tail_begin, long tail_end) {
  fprintf(stderr, "  Computing initial ranks: ");
  long double start = utils::wclock();

  long tail_length = tail_end - tail_begin;
  long stream_max_block_size = (tail_length + max_threads - 1) / max_threads;
  long n_threads = (tail_length + stream_max_block_size - 1) / stream_max_block_size;

  result.resize(n_threads);
  std::thread **threads = new std::thread*[n_threads];
  for (int t = 0; t < n_threads; ++t) {
    long stream_block_beg = tail_begin + t * stream_max_block_size;
    long stream_block_end = std::min(stream_block_beg + stream_max_block_size, tail_end);

    threads[t] = new std::thread(parallel_smaller_suffixes2<saidx_t>, block, block_beg, block_end, 
        text_length, block_partial_sa,
        text_filename, stream_block_end, std::ref(result[t]), tail_gt_begin_reversed);
  }

  for (int t = 0; t < n_threads; ++t) threads[t]->join();
  for (int t = 0; t < n_threads; ++t) delete threads[t];
  delete[] threads;
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
}


// compute the gap for an arbitrary range of suffixes of tail.
// this version is more general, and can be used also when processing half-blocks.

template<typename block_offset_type>
void compute_gap(
    rank4n<> *rank,
    buffered_gap_array *gap,
    long tail_begin,
    long tail_end,
    long text_length,
    std::string text_filename,
    std::vector<long> initial_ranks,
    long max_threads,
    unsigned char block_last_symbol,
    long block_isa0,
    multifile *tail_gt_begin_reversed,
    multifile *newtail_gt_begin_reversed,
    long stream_buffer_size) {

  long tail_length = tail_end - tail_begin;
  long stream_max_block_size = (tail_length + max_threads - 1) / max_threads;
  long n_threads = (tail_length + stream_max_block_size - 1) / stream_max_block_size;


  // Get symbol counts of a block and turn into exclusive partial sum.
  long *count = new long[256];
  std::copy(rank->c_rank, rank->c_rank + 256, count);
  ++count[block_last_symbol];
  --count[0];
  for (long j = 0, s = 0, t; j < 256; ++j) {
    t = count[j];
    count[j] = s;
    s += t;
  }


  // 5b. Allocate the gap array, do the streaming and store gap to disk.
  fprintf(stderr, "  Stream:");
  long double stream_start = utils::wclock();

  // Allocate buffers.
  long n_stream_buffers = 2 * max_threads;
  buffer<block_offset_type> **buffers = new buffer<block_offset_type>*[n_stream_buffers];
  for (long i = 0L; i < n_stream_buffers; ++i)
    buffers[i] = new buffer<block_offset_type>(stream_buffer_size, max_threads);

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
    newtail_gt_begin_reversed->add_file(text_length - stream_block_end, text_length - stream_block_beg, gt_filenames[t]);

    streamers[t] = new std::thread(parallel_stream<block_offset_type>, full_buffers, empty_buffers, stream_block_beg,
        stream_block_end, initial_ranks[t], count, block_isa0, rank, block_last_symbol, text_filename, text_length,
        std::ref(gt_filenames[t]), &info, t, gap->m_length, stream_buffer_size, tail_gt_begin_reversed, max_threads);
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
  fprintf(stderr,"\r  Stream: 100.0%%. Tail length = %ld. Time: %.2Lf. Threads: %ld. "
      "Speed: %.2LfMiB/s (avg), %.2LfMiB/s (total)\n", tail_length,
      stream_time, info.m_thread_count, speed / n_threads, speed);
}




//=============================================================================
// Note: on entry to the function it holds: 5.2 * block_size <= ram_use
//==============================================================================
template<typename block_offset_type>
distributed_file<block_offset_type> *process_block(
    long block_beg,
    long block_end,
    long max_block_size,
    long text_length,
    long ram_use,
    std::string text_filename,
    long max_threads,
    multifile *newtail_gt_begin_reversed,
    multifile *tail_gt_begin_reversed,
    long stream_buffer_size) {

  long block_size = block_end - block_beg;
  long block_id = block_beg / max_block_size;

  fprintf(stderr, "  Reading block: ");
  long double read_start = utils::wclock();
  unsigned char *block = (unsigned char *)malloc(block_size);
  utils::read_block(text_filename, block_beg, block_size, block);
  fprintf(stderr, "%.2Lf\n", utils::wclock() - read_start);


  unsigned char block_last_symbol = block[block_size - 1];

  bool last_block = (block_end == text_length);
  bool first_block = (block_beg == 0);

  fprintf(stderr, "  block_size = %ld (%.2LfMiB)\n", block_size, 1.L * block_size / (1 << 20));

  fprintf(stderr, "\n******************** Running inmem SAscan *********************\n");
  unsigned char *sabwt = (unsigned char *)malloc(block_size * (sizeof(block_offset_type) + 1) + 1);
  block_offset_type *partial_sa = (block_offset_type *)sabwt;
  unsigned char *bwt = (unsigned char *)(partial_sa + block_size);
  bitvector *block_gt_begin = NULL;
  if (!first_block) block_gt_begin = new bitvector(block_size, max_threads);
  long block_isa0;
  inmem_sascan<block_offset_type>(block, block_size, sabwt, max_threads, !last_block,
      !first_block, block_gt_begin, -1, block_beg, block_end, text_length, text_filename,
      tail_gt_begin_reversed, &block_isa0);
  fprintf(stderr, "***************************************************************\n\n");

  if (!first_block) {
    fprintf(stderr, "  Saving block gt_begin to file: ");
    long double block_gt_begin_save_start = utils::wclock();
    std::string block_gt_begin_reversed_filename = text_filename + ".block_gt_begin_reversed." + utils::random_string_hash();
    block_gt_begin->save_reversed(block_gt_begin_reversed_filename, block_size);
    newtail_gt_begin_reversed->add_file(text_length - block_end, text_length - block_beg, block_gt_begin_reversed_filename);
    delete block_gt_begin;
    fprintf(stderr, "%.2Lf\n", utils::wclock() - block_gt_begin_save_start);
  }

  fprintf(stderr, "  Writing partial SA to disk (using %lu-byte ints): ", sizeof(block_offset_type));
  long double writing_sa_start = utils::wclock();
  std::string sa_fname = text_filename + ".partial_sa." + utils::intToStr(block_id); // base should be SA file.
  distributed_file<block_offset_type> *result = new distributed_file<block_offset_type>(
      sa_fname.c_str(), std::max((long)sizeof(block_offset_type), ram_use / 10L));
  result->initialize_writing(4 << 20);
  for (long i = 0; i < block_size; ++i)
    result->write(partial_sa[i]);
  result->finish_writing();
  fprintf(stderr, "%.2Lf\n", utils::wclock() - writing_sa_start);


  if (!last_block) {
#if 0
    std::vector<long> initial_ranks;
    compute_initial_ranks<block_offset_type>(block, block_beg, block_end, text_length, partial_sa,
        text_filename, initial_ranks, tail_gt_begin_reversed, max_threads, block_end, text_length);
    free(block);

    fprintf(stderr, "  Building the rank data structure: ");
    long double building_rank_start = utils::wclock();
    rank4n<> *rank = new rank4n<>(bwt, block_size, max_threads);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - building_rank_start);

    buffered_gap_array *gap = new buffered_gap_array(block_size + 1, bwt);
    compute_gap<block_offset_type>(rank, gap, block_end, text_length, text_length, text_filename, initial_ranks,
        max_threads, block_last_symbol, block_isa0, tail_gt_begin_reversed, newtail_gt_begin_reversed, stream_buffer_size);
    gap->save_to_file(text_filename + ".gap." + utils::intToStr(block_id));

    delete gap;
    delete rank;
#else
    rank4n<> *rank = new rank4n<>(bwt, block_size, max_threads);
    buffered_gap_array *gap = new buffered_gap_array(block_size + 1, bwt);

    long tail_beg = block_end;
    while (tail_beg != text_length) {
      long tmp_tail_size = utils::random_long(1L, text_length - tail_beg);
      long tail_end = tail_beg + tmp_tail_size;

      std::vector<long> initial_ranks;
      compute_initial_ranks<block_offset_type>(block, block_beg, block_end, text_length, partial_sa,
          text_filename, initial_ranks, tail_gt_begin_reversed, max_threads, tail_beg, tail_end);

      unsigned char *temp = (unsigned char *)malloc(block_size + 1);
      buffered_gap_array *temp_gap = new buffered_gap_array(block_size + 1, temp);

      compute_gap<block_offset_type>(rank, temp_gap, tail_beg, tail_end, text_length, text_filename, initial_ranks,
         max_threads, block_last_symbol, block_isa0, tail_gt_begin_reversed, newtail_gt_begin_reversed, stream_buffer_size);

      for (long i = 0; i <= block_size; ++i) {
        long x = (long)gap->m_count[i] + (long)temp_gap->m_count[i];
        while (x >= 256) {
          gap->add_excess(i);
          x -= 256;
        }
        gap->m_count[i] = x;
      }
      for (long jj = 0; jj < temp_gap->m_excess_filled; ++jj)
        gap->add_excess(temp_gap->m_excess[jj]);

      free(temp);
      delete temp_gap;

      tail_beg = tail_end;
    }

    free(block);

    gap->save_to_file(text_filename + ".gap." + utils::intToStr(block_id));
    delete gap;
    delete rank;
#endif
  } else free(block);
  
  free(sabwt);

  return result;
}




//=============================================================================
// Compute partial SAs and gap arrays and write to disk.
// Return the array of handlers to distributed files as a result.
//=============================================================================
template<typename block_offset_type>
distributed_file<block_offset_type> **partial_sufsort(std::string text_filename,
    long text_length, long max_block_size, long ram_use, long max_threads, long stream_buffer_size) {
  fprintf(stderr, "sizeof(block_offset_type) = %lu\n\n", sizeof(block_offset_type));

  long n_blocks = (text_length + max_block_size - 1) / max_block_size;
  multifile *tail_gt_begin_reversed = NULL;
  distributed_file<block_offset_type> **distrib_files = new distributed_file<block_offset_type>*[n_blocks];

  for (long block_id = n_blocks - 1; block_id >= 0; --block_id) {
    long block_beg = max_block_size * block_id;
    long block_end = std::min(block_beg + max_block_size, text_length);
    fprintf(stderr, "Processing block %ld/%ld [%ld..%ld):\n", n_blocks - block_id, n_blocks, block_beg, block_end);

    multifile *newtail_gt_begin_reversed = new multifile();
    distrib_files[block_id] = process_block<block_offset_type>(block_beg, block_end, max_block_size, text_length,
        ram_use, text_filename, max_threads, newtail_gt_begin_reversed, tail_gt_begin_reversed, stream_buffer_size);

    delete tail_gt_begin_reversed;
    tail_gt_begin_reversed = newtail_gt_begin_reversed;
  }

  delete tail_gt_begin_reversed;
  return distrib_files;
}

#endif // __PARTIAL_SUFSORT_H_INCLUDED
