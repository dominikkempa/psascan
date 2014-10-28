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


//------------------------------------------------------------------------------
// For this to work you need to run sascan e.g. as follows:
//
//   $ nohup sudo -b numactl -i 0-1 ./sascan /data01/sac-corpus/wiki.6G -m 20000
//
//------------------------------------------------------------------------------
void drop_cache() {
/*  long double start = utils::wclock();
  fprintf(stderr, "  Clearing cache: ");
  fprintf(stderr, "Before:\n");
  utils::execute("free -m");
  utils::execute("echo 3 | tee /proc/sys/vm/drop_caches");
  fprintf(stderr, "After:\n");
  utils::execute("free -m");
  fprintf(stderr, "Clearing time: %.2Lf\n", utils::wclock() - start);*/
}

template<typename saidx_t>
void compute_initial_ranks(unsigned char *block, long block_beg,
    long block_end, long text_length, saidx_t *block_partial_sa,
    std::string text_filename, std::vector<long> &result,
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
        text_filename, stream_block_end, std::ref(result[t]));
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

  long block_tail_beg = block_end;
  long block_tail_end = text_length;

  long block_id = block_beg / max_block_size;
  long block_size = block_end - block_beg;

  bool last_block = (block_end == text_length);
  bool first_block = (block_beg == 0);

  long left_block_size;
  if (last_block) left_block_size = std::min(block_size, std::max(1L, ram_use / 10L));
  else left_block_size = std::max(1L, block_size / 2L);
  long right_block_size = block_size - left_block_size;

  long left_block_beg = block_beg;
  long left_block_end = block_beg + left_block_size;
  long right_block_beg = left_block_end;
  long right_block_end = block_end;

  fprintf(stderr, "  Block size = %ld (%.2LfMiB)\n", block_size, (1.L * block_size / (1 << 20)));
  fprintf(stderr, "  Left block size = %ld (%.2LfMiB)\n", left_block_size, 1.L * left_block_size / (1 << 20));
  fprintf(stderr, "  Right block size = %ld (%.2LfMiB)\n", right_block_size, 1.L * right_block_size / (1 << 20));

  std::vector<long> block_initial_ranks;



  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Process right block -------------------------------------------------------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------
  long right_block_i0;
  std::string right_block_psa_filename = text_filename + ".right_block.partial_sa"; // XXX different base name!
  std::string right_block_bwt_filename = text_filename + ".right_block.bwt";
  std::string right_block_gt_begin_filename = text_filename + ".right_block.gt_begin";

  distributed_file<block_offset_type> *right_block_psa = NULL;
  multifile *right_block_gt_begin_reversed = NULL;
  unsigned char right_block_last_symbol = 0;
  unsigned char block_last_symbol = 0;

  if (right_block_size > 0) {
    fprintf(stderr, "  Processing right half-block [%ld..%ld):\n", right_block_beg, right_block_end);

    // At this point, we don't have any memory allocated.
    // 1. Allocate space and read right block into memory.
    unsigned char *right_block = (unsigned char *)malloc(right_block_size);
    fprintf(stderr, "  Reading right block: ");
    long double right_block_read_start = utils::wclock();
    utils::read_block(text_filename, right_block_beg, right_block_size, right_block);
    block_last_symbol = right_block_last_symbol = right_block[right_block_size - 1];
    fprintf(stderr, "%.2Lf\n", utils::wclock() - right_block_read_start);
 
    // 2. Allocate the space to hold the suffix array and BWT (if needed) of right block.
    unsigned char *right_block_sabwt = (unsigned char *)malloc(right_block_size * (sizeof(block_offset_type) + 1));
    block_offset_type *right_block_partial_sa_ptr = (block_offset_type *)right_block_sabwt;
    unsigned char *right_block_bwt = (unsigned char *)(right_block_partial_sa_ptr + right_block_size);

    // 3. Allocate the space for right_block_gt_begin (always necessary).
    bitvector *right_block_gt_begin_bv = new bitvector(right_block_size, max_threads);

    // 4. Compute partial suffix array for right block.
    fprintf(stderr, "\n******************** Running inmem SAscan *********************\n");
    inmem_sascan<block_offset_type>(right_block, right_block_size, right_block_sabwt, max_threads,
        !last_block, true, right_block_gt_begin_bv, -1, right_block_beg, right_block_end, text_length,
        text_filename, tail_gt_begin_reversed, &right_block_i0);
    fprintf(stderr, "****************************************************************\n\n");

    // compute the first term of initial_ranks for the block.
    if (!last_block) {
      compute_initial_ranks<block_offset_type>(right_block, right_block_beg, right_block_end, text_length, right_block_partial_sa_ptr,
          text_filename, block_initial_ranks, max_threads, block_tail_beg, block_tail_end);
    }

    // 5. Save the partial suffix array of the right block to disk using distributed file.
    fprintf(stderr, "  Saving partial sa of right block to disk (using %lu-byte ints): ", sizeof(block_offset_type));
    long double right_block_psa_save_start = utils::wclock();
    right_block_psa = new distributed_file<block_offset_type>(right_block_psa_filename.c_str(), 100L << 20);
    right_block_psa->initialize_writing(4L << 20);
    for (long i = 0; i < right_block_size; ++i) {
      //fprintf(stderr, "Writing %ld\n", (long)right_block_partial_sa_ptr[i]);
      right_block_psa->write(right_block_partial_sa_ptr[i]);
    }
    right_block_psa->finish_writing();
    fprintf(stderr, "%.2Lf\n", utils::wclock() - right_block_psa_save_start);

    // 6. If necessary, save the bwt of the right block on disk.
    if (!last_block) {
      fprintf(stderr, "  Saving bwt of right block to disk: ");
      long double right_block_bwt_save_start = utils::wclock();
      utils::write_objects_to_file(right_block_bwt, right_block_size, right_block_bwt_filename);
      fprintf(stderr, "%.2Lf\n", utils::wclock() - right_block_bwt_save_start);
    }

    // 7. Write reversed gt_begin for the right block to disk.
    fprintf(stderr, "  Writing gt_begin of right block to disk: ");
    long double right_block_gt_begin_save_start = utils::wclock();
    right_block_gt_begin_bv->save_reversed(right_block_gt_begin_filename, right_block_size);
    right_block_gt_begin_reversed = new multifile();
    right_block_gt_begin_reversed->add_file(text_length - right_block_end, text_length - right_block_beg,
      right_block_gt_begin_filename);
    delete right_block_gt_begin_bv;
    fprintf(stderr, "%.2Lf\n", utils::wclock() - right_block_gt_begin_save_start);

    free(right_block_sabwt);
    free(right_block);

    drop_cache();
  }







  //---------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Process left block -------------------------------------------------------------------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // At this point no memory is allocated.
  fprintf(stderr, "  Processing left half-block [%ld..%ld):\n", left_block_beg, left_block_end);

  // 1. allocate the space and read the left block.
  fprintf(stderr, "  Reading left block: ");
  long double left_block_reading_start = utils::wclock();
  unsigned char *left_block = (unsigned char *)malloc(left_block_size);
  utils::read_block(text_filename, left_block_beg, left_block_size, left_block);
  fprintf(stderr, "%.2Lf\n", utils::wclock() - left_block_reading_start);

  // 2. allocate space to hold sa and bwt of left block.
  unsigned char *left_block_sabwt = (unsigned char *)malloc(left_block_size * (sizeof(block_offset_type) + 1) + 1);
  block_offset_type *left_block_psa_ptr = (block_offset_type *)left_block_sabwt;
  unsigned char *left_block_bwt_ptr = (unsigned char *)(left_block_psa_ptr + left_block_size);

  long left_block_i0;
  bitvector *left_block_gt_beg_bv = NULL;
  if (!first_block) left_block_gt_beg_bv = new bitvector(left_block_size, max_threads);

  // Run inmem SAscan.
  fprintf(stderr, "\n******************** Running inmem SAscan *********************\n");
  inmem_sascan<block_offset_type>(left_block, left_block_size, left_block_sabwt, max_threads,
      (right_block_size > 0), !first_block, left_block_gt_beg_bv, -1,
      left_block_beg, left_block_end, text_length, text_filename,
      right_block_gt_begin_reversed, &left_block_i0);
  fprintf(stderr, "****************************************************************\n\n");

  // Save gt_begin for the left block as a part of newtail_gt_begin_reversed.
  if (!first_block) {
    fprintf(stderr, "  Writing gt_begin of left block to disk: ");
    long double left_block_gt_begin_save_start = utils::wclock();
    std::string left_block_gt_begin_reversed_filename = text_filename + ".left_block_gt_beg_rev" + utils::random_string_hash();
    left_block_gt_beg_bv->save_reversed(left_block_gt_begin_reversed_filename, left_block_size);

    newtail_gt_begin_reversed->add_file(text_length - left_block_end, text_length - left_block_beg,
        left_block_gt_begin_reversed_filename);
    delete left_block_gt_beg_bv;
    fprintf(stderr, "%.2Lf\n", utils::wclock() - left_block_gt_begin_save_start);
  }

  if (!last_block) {
    // Compute the second term for block_initial_ranks
    std::vector<long> block_initial_ranks_second_term;
    compute_initial_ranks<block_offset_type>(left_block, left_block_beg, left_block_end, text_length,
        left_block_psa_ptr, text_filename, block_initial_ranks_second_term, max_threads, block_tail_beg, block_tail_end);

    for (size_t j = 0; j < block_initial_ranks_second_term.size(); ++j)
      block_initial_ranks[j] += block_initial_ranks_second_term[j];
  }

  if (right_block_size > 0) {
    // Compute initial ranks for streaming right block.
    std::vector<long> initial_ranks2;
    compute_initial_ranks<block_offset_type>(left_block, left_block_beg, left_block_end, text_length, left_block_psa_ptr,
        text_filename, initial_ranks2, max_threads, right_block_beg, right_block_end);

    unsigned char left_block_last_symbol = left_block[left_block_size - 1];
    free(left_block);

    // Save bwt of the left block to disk.
    std::string left_block_bwt_filename = text_filename + ".left_block.bwt";
    if (!last_block) {
      fprintf(stderr, "  Saving bwt of left block to disk: ");
      long double left_block_bwt_save_start = utils::wclock();
      utils::write_objects_to_file(left_block_bwt_ptr, left_block_size, left_block_bwt_filename);
      fprintf(stderr, "%.2Lf\n", utils::wclock() - left_block_bwt_save_start);
    }

    drop_cache();

    // 5. Build the rank over bwt of left block.
    fprintf(stderr, "  Building the rank data structure for left block: ");
    long double left_block_rank_build_start = utils::wclock();
    rank4n<> *left_block_rank = new rank4n<>(left_block_bwt_ptr, left_block_size, max_threads);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - left_block_rank_build_start);

    // At this point we are using 10 * left_block_size of RAM. The maximal amount.
    // 6. Overwrite bwt of left block with gap array.
    buffered_gap_array *left_block_gap = new buffered_gap_array(left_block_size + 1, left_block_bwt_ptr);

    // Compute gap array for the left block.
    compute_gap<block_offset_type>(left_block_rank, left_block_gap, right_block_beg, right_block_end,
        text_length, text_filename, initial_ranks2, max_threads, left_block_last_symbol, left_block_i0, 
        right_block_gt_begin_reversed, newtail_gt_begin_reversed, stream_buffer_size);

    delete left_block_rank;
    delete right_block_gt_begin_reversed;

    // merge partial suffix array of the left block (in memory) with the partial suffix array of the right block (on disk)
    // to obtain partial suffix array of the block.
    fprintf(stderr, "  Merging partial SA of left and right block and writing to disk (using %lu-byte ints): ", sizeof(block_offset_type));
    long double writing_sa_start = utils::wclock();
    std::string sa_fname = text_filename + ".partial_sa." + utils::intToStr(block_id); // base should be SA file.
    distributed_file<block_offset_type> *block_psa = new distributed_file<block_offset_type>(sa_fname.c_str(), std::max((long)sizeof(block_offset_type), ram_use / 10L));

    block_psa->initialize_writing(4 << 20);
    right_block_psa->initialize_reading(4 << 20);

    left_block_gap->start_sequential_access();
    for (long i = 0; i <= left_block_size; ++i) {
      long gap_i = left_block_gap->get_next();

      for (long j = 0; j < gap_i; ++j) {
        long next_value = left_block_size + right_block_psa->read();
        block_psa->write(next_value);
      }
      if (i < left_block_size)
        block_psa->write(left_block_psa_ptr[i]);
    }
    right_block_psa->finish_reading();
    block_psa->finish_writing();
    delete right_block_psa;
    fprintf(stderr, "%.2Lf\n", utils::wclock() - writing_sa_start);

    if (!last_block) {
      // Merge BWT of the right and left block (both are on disk).
      fprintf(stderr, "  Merging partial BWT of left and right block: ");
      long double left_right_bwt_merge_start = utils::wclock();
      unsigned char *block_bwt = (unsigned char *)malloc(block_size);
      long block_bwt_filled = 0;

      stream_reader<unsigned char> *left_block_bwt_reader  = new stream_reader<unsigned char> (left_block_bwt_filename);
      stream_reader<unsigned char> *right_block_bwt_reader = new stream_reader<unsigned char>(right_block_bwt_filename);

      left_block_gap->start_sequential_access();
      long block_i0 = 0;
      long right_block_extracted = 0;
      for (long i = 0; i <= left_block_size; ++i) {
        long gap_i = left_block_gap->get_next();

        for (long j = 0; j < gap_i; ++j)  {
          long ccc = right_block_bwt_reader->read();
          if (right_block_extracted == right_block_i0) ccc = left_block_last_symbol;
          block_bwt[block_bwt_filled++] = ccc;
          ++right_block_extracted;
        }
        if (i == left_block_i0) block_i0 = block_bwt_filled;
        if (i < left_block_size) block_bwt[block_bwt_filled++] = left_block_bwt_reader->read();
      }
      left_block_gap->stop_sequential_access();
      delete left_block_bwt_reader;
      delete right_block_bwt_reader;
      utils::file_delete(left_block_bwt_filename);
      utils::file_delete(right_block_bwt_filename);

      free(left_block_sabwt);

      fprintf(stderr, "%.2Lf\n", utils::wclock() - left_right_bwt_merge_start);

      drop_cache();

      fprintf(stderr, "  Building rank for the block: ");
      long double whole_block_rank_build_start = utils::wclock();
      rank4n<> *block_rank = new rank4n<>(block_bwt, block_size, max_threads);
      fprintf(stderr, "%.2Lf\n", utils::wclock() - whole_block_rank_build_start);

      free(block_bwt);
      unsigned char *block_gap_count = (unsigned char *)malloc(block_size + 1);
      buffered_gap_array *block_gap = new buffered_gap_array(block_size + 1, block_gap_count);

      // compute gap for the block.
      compute_gap<block_offset_type>(block_rank, block_gap, block_tail_beg, block_tail_end, text_length, text_filename,
          block_initial_ranks, max_threads, block_last_symbol, block_i0, tail_gt_begin_reversed, newtail_gt_begin_reversed, 
          stream_buffer_size);

      delete block_rank;

      block_gap->save_to_file(text_filename + ".gap." + utils::intToStr(block_id));
      delete block_gap;
      free(block_gap_count);
    } else {
      free(left_block_sabwt);
      left_block_gap->stop_sequential_access(); // XXX we delete it right away, can we skip that?
    }

    delete left_block_gap;

    drop_cache();

    return block_psa;
  } else {
    free(left_block);

    fprintf(stderr, "  Writing partial SA of block to disk: ");
    long double block_partial_sa_writing_start = utils::wclock();
    std::string sa_fname = text_filename + ".partial_sa." + utils::intToStr(block_id); // base should be SA file.
    distributed_file<block_offset_type> *block_psa = new distributed_file<block_offset_type>(sa_fname.c_str(), std::max((long)sizeof(block_offset_type), ram_use / 10L));

    block_psa->initialize_writing(4 << 20);
    for (long j = 0; j < left_block_size; ++j) {
      block_psa->write(left_block_psa_ptr[j]);
    }
    block_psa->finish_writing();
    fprintf(stderr, "%.2Lf\n", utils::wclock() - block_partial_sa_writing_start);

    free(left_block_sabwt);

    drop_cache();
    return block_psa;
  }
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
