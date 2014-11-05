// XXX fix the interface issue of inmem_sascan
// XXX check if inmem_sascan can compute gt_begin reversed

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

#include <sys/stat.h>
#include <fcntl.h>

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
#include "smaller_suffixes.h"
#include "radixsort.h"
#include "buffer.h"
#include "update.h"
#include "stream_info.h"
#include "aux_parallel.h"
#include "multifile_bitvector.h"
#include "inmem_sascan/inmem_sascan.h"
#include "half_block_info.h"
#include "bwt_merge.h"


template<typename saidx_t>
void compute_initial_ranks(unsigned char *block, long block_beg,
    long block_end, long text_length, saidx_t *block_partial_sa,
    std::string text_filename, std::vector<long> &result,
    long max_threads, long tail_begin, long tail_end) {
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
}


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
  std::copy(rank->c_rank, rank->c_rank + 256, count);
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


//=============================================================================
// The main function processing the block.
//=============================================================================
// INPUT
// In addition to the block boundaries, ram budget and other integer
// parameters, the function requires that gt_begin was computed for the tail
// of the text and the multifile representation of the bitvector on disk is
// given as an input.
//
// OUTPUT
// The function produces the following output:
// - partial suffix array of the block, stored on disk as a distributed
//   file. The handle to this file is returned from the function. Note that
//   any distributed file is unnamed (the actual files created on disk have
//   names decided by the implementation of the distributed_file class), so
//   the handle is the only way of accessing distributed files (well, this
//   is actually not yet implemented).
// - gap array of the block, stored on disk as a regular file using v-byte
//   encoding. The name of the file is output_filename + ".gap." + block_id.
// - gt_begin of the new tail (which consits of the block and the old tail)
//   stored on disk as a multifile. The handle to this multifile is returned
//   via the reference.
//
// NOTE
// * The multifile representation is different from distributed file!
// * On entry to the function it holds: 5.2 * block_size <= ram_use
// * Next version of this function will return two distributed files tather
//   than one -- each holding the partial suffix of the half-block. For this
//   to work, we will have to change few this here and there, but overall it
//   should save some I/O.
//==============================================================================
template<typename block_offset_type>
void process_block(long block_beg, long block_end,
    long text_length, long ram_use, long max_threads, long stream_buffer_size,
    std::string text_filename, std::string output_filename,
    multifile *newtail_gt_begin_rev, multifile *tail_gt_begin_rev,
    std::vector<half_block_info<block_offset_type> > &hblock_info) {

  long block_size = block_end - block_beg;

  if (block_end != text_length && block_size <= 1) {
    fprintf(stderr, "Error: any block other than the last one has to be of length at least two.\n");
    std::exit(EXIT_FAILURE);
  }

  long block_tail_beg = block_end;
  long block_tail_end = text_length;

  bool last_block = (block_end == text_length);
  bool first_block = (block_beg == 0);

  // Note: left_block_size > 0.
  long left_block_size;
  if (!last_block) left_block_size = std::max(1L, block_size / 2L);
  else left_block_size = std::min(block_size, std::max(1L, ram_use / 10L));
  long right_block_size = block_size - left_block_size;
  long left_block_beg = block_beg;
  long left_block_end = block_beg + left_block_size;
  long right_block_beg = left_block_end;
  long right_block_end = block_end;

  fprintf(stderr, "  Block size = %ld (%.2LfMiB)\n", block_size, 1.L * block_size / (1 << 20));
  fprintf(stderr, "  Left half-block size = %ld (%.2LfMiB)\n", left_block_size, 1.L * left_block_size / (1 << 20));
  fprintf(stderr, "  Right half-block size = %ld (%.2LfMiB)\n", right_block_size, 1.L * right_block_size / (1 << 20));

  std::vector<long> block_initial_ranks;
  unsigned char block_last_symbol = 0;

  long right_block_i0 = 0;
  long left_block_i0 = 0;

  std::string right_block_pbwt_fname = output_filename + "." + utils::random_string_hash();
  std::string right_block_gt_begin_rev_fname = output_filename + "." + utils::random_string_hash();

  half_block_info<block_offset_type> info_left;
  half_block_info<block_offset_type> info_right;

  info_left.beg = left_block_beg;
  info_left.end = left_block_end;
  if (right_block_size > 0) {
    info_right.beg = right_block_beg;
    info_right.end = right_block_end;
  }


  //----------------------------------------------------------------------------
  // STEP 1: Process right half-block.
  //
  // The output of this step if the following:
  // - if right_block_size == 0, then nothing is produced and pointers
  //   right_block_psa and right_block_gt_begin_rev remain NULL,
  // - otherwise right_block_psa and right_block_gt_begin_rev are not NULL.
  //   right_block_psa is a handler to the distributed file containing partial
  //   suffix array of the right half-block. right_block_gt_begin_rev is a
  //   handler to multifile containing gt_begin or the right block. Both
  //   structures are of negligible size. The actual data is store on disk.
  //
  // NOTE: The right half-block can be empty ONLY if the block under
  // consideration is the last block of text.
  //----------------------------------------------------------------------------
  multifile *right_block_gt_begin_rev = NULL;

  if (right_block_size > 0) {
    fprintf(stderr, "  Process right half-block:\n");

    // 1.a
    //
    // Read the right half-block from disk.
    fprintf(stderr, "    Read: ");
    unsigned char *right_block = (unsigned char *)malloc(right_block_size);
    long double right_block_read_start = utils::wclock();
    utils::read_block(text_filename, right_block_beg, right_block_size, right_block);
    block_last_symbol = right_block[right_block_size - 1];
    long double right_block_read_time = utils::wclock() - right_block_read_start;
    long double right_block_read_io = (right_block_size / (1024.L * 1024)) / right_block_read_time;
    fprintf(stderr, "%.2Lf (I/O: %.2LfMiB/s)\n", right_block_read_time, right_block_read_io);
 
    // 1.b
    //
    // Compute partial sa, bwt and gt_begin of the right half-block.

    // Allocate suffix array, bwt and gt_begin.
    unsigned char *right_block_sabwt = (unsigned char *)malloc(right_block_size * (sizeof(block_offset_type) + 1));
    block_offset_type *right_block_psa_ptr = (block_offset_type *)right_block_sabwt;
    unsigned char *right_block_bwt = (unsigned char *)(right_block_psa_ptr + right_block_size);
    bitvector *right_block_gt_begin_bv = new bitvector(right_block_size);
    fprintf(stderr, "    Internal memory sufsort: ");
    long double right_block_sascan_start = utils::wclock();

    // Close stderr.
    int stderr_backup;
    int stderr_temp;
    std::fflush(stderr);
    stderr_backup = dup(2);
    stderr_temp = open("/dev/null", O_WRONLY);
    dup2(stderr_temp, 2);
    close(stderr_temp);

    // Run in-memory SAscan.
    inmem_sascan<block_offset_type>(right_block, right_block_size, right_block_sabwt, max_threads,
        !last_block, true, right_block_gt_begin_bv, -1, right_block_beg, right_block_end, text_length,
        text_filename, tail_gt_begin_rev, &right_block_i0);

    // Restore stderr.
    std::fflush(stderr);
    dup2(stderr_backup, 2);
    close(stderr_backup);

    // Print summary.
    long double right_block_sascan_time = utils::wclock() - right_block_sascan_start;
    long double right_block_sascan_speed = (right_block_size / (1024.L * 1024)) / right_block_sascan_time;
    fprintf(stderr, "%.2Lf (%.2LfMiB/s)\n", right_block_sascan_time, right_block_sascan_speed);
 
    // 1.c
    //
    // Compute the first term of initial_ranks for the block.
    if (!last_block) {
      fprintf(stderr, "    Compute initial tail ranks (part 1): ");
      long double initial_ranks_first_term_start = utils::wclock();
      compute_initial_ranks<block_offset_type>(right_block, right_block_beg, right_block_end, text_length,
          right_block_psa_ptr, text_filename, block_initial_ranks, max_threads, block_tail_beg, block_tail_end);
      fprintf(stderr, "%.2Lf\n", utils::wclock() - initial_ranks_first_term_start);
    }
    free(right_block);

    // 1.d
    //
    // Write the partial SA of the right half-block to disk.
    fprintf(stderr, "    Write partial SA to disk: ");
    long double right_psa_save_start = utils::wclock();
    long right_psa_max_part_length = std::max((long)sizeof(block_offset_type), ram_use / 10L);
    info_right.psa = new distributed_file<block_offset_type>(output_filename,
        right_psa_max_part_length, right_block_psa_ptr, right_block_psa_ptr + right_block_size);
    long double right_psa_save_time = utils::wclock() - right_psa_save_start;
    long double right_psa_save_io = ((right_block_size * sizeof(block_offset_type)) / (1024.L * 1024)) / right_psa_save_time;
    fprintf(stderr, "%.2Lf (I/O: %.2LfMiB/s)\n", right_psa_save_time, right_psa_save_io);

    // 1.e
    //
    // Write the BWT of the right half-block on disk.
    if (!last_block) {
      fprintf(stderr, "    Write BWT to disk: ");
      long double right_bwt_save_start = utils::wclock();
      utils::write_objects_to_file(right_block_bwt, right_block_size, right_block_pbwt_fname);
      long double right_bwt_save_time = utils::wclock() - right_bwt_save_start;
      long double right_bwt_save_io = (right_block_size / (1024.L * 1024)) / right_bwt_save_time;
      fprintf(stderr, "%.2Lf (I/O: %.2LfMiB/s)\n", right_bwt_save_time, right_bwt_save_io);
    }
    free(right_block_sabwt);

    // 1.f
    //
    // Write gt_begin of the right half-block to disk.
    fprintf(stderr, "    Write gt_begin to disk: ");
    long double right_gt_begin_save_start = utils::wclock();
    right_block_gt_begin_bv->save_reversed(right_block_gt_begin_rev_fname, right_block_size);
    right_block_gt_begin_rev = new multifile();
    right_block_gt_begin_rev->add_file(text_length - right_block_end, text_length - right_block_beg, right_block_gt_begin_rev_fname);
    delete right_block_gt_begin_bv;
    long double right_gt_begin_save_time = utils::wclock() - right_gt_begin_save_start;
    long double right_gt_begin_save_io = (right_block_size / (8.L * (1 << 20))) / right_gt_begin_save_time;
    fprintf(stderr, "%.2Lf (I/O: %.2LfMiB/s)\n", right_gt_begin_save_time, right_gt_begin_save_io);
  }


  //----------------------------------------------------------------------------
  // STEP 2: Process left half-block.
  //
  // At this point in RAM reside only handlers to gt_begin and partial SA of
  // the right half-block (if it was empty they are both NULL). Both handles
  // take negligible space, the actual data structures are stored on disk.
  //----------------------------------------------------------------------------
  fprintf(stderr, "  Process left half-block:\n");

  // 2.a
  //
  // Read the left half-block from disk.
  fprintf(stderr, "    Read: ");
  long double left_block_read_start = utils::wclock();
  unsigned char *left_block = (unsigned char *)malloc(left_block_size);
  utils::read_block(text_filename, left_block_beg, left_block_size, left_block);
  unsigned char left_block_last = left_block[left_block_size - 1];
  long double left_block_read_time = utils::wclock() - left_block_read_start;
  long double left_block_read_io = (left_block_size / (1024.L * 1024)) / left_block_read_time;
  fprintf(stderr, "%.2Lf (I/O: %.2LfMiB/s)\n", left_block_read_time, left_block_read_io);

  // 2.b
  //
  // Compute partial SA, BWT and gt_begin for left half-block.

  // Allocate suffix array, bwt and gt_begin.
  unsigned char *left_block_sabwt = (unsigned char *)malloc(left_block_size * (sizeof(block_offset_type) + 1) + 1);
  block_offset_type *left_block_psa_ptr = (block_offset_type *)left_block_sabwt;
  unsigned char *left_block_bwt_ptr = (unsigned char *)(left_block_psa_ptr + left_block_size);
  bitvector *left_block_gt_beg_bv = NULL;
  if (!first_block) left_block_gt_beg_bv = new bitvector(left_block_size);  
  fprintf(stderr, "    Internal memory sufsort: ");
  long double left_block_sascan_start = utils::wclock();

  // Close stderr.
  int stderr_backup;
  int stderr_temp;
  std::fflush(stderr);
  stderr_backup = dup(2);
  stderr_temp = open("/dev/null", O_WRONLY);
  dup2(stderr_temp, 2);
  close(stderr_temp);

  // Run in-memory SAscan.
  inmem_sascan<block_offset_type>(left_block, left_block_size, left_block_sabwt, max_threads,
      (right_block_size > 0), !first_block, left_block_gt_beg_bv, -1, left_block_beg,
      left_block_end, text_length, text_filename, right_block_gt_begin_rev, &left_block_i0);

  // Restore stderr.
  std::fflush(stderr);
  dup2(stderr_backup, 2);
  close(stderr_backup);

  // Print summary.
  long double left_block_sascan_time = utils::wclock() - left_block_sascan_start;
  long double left_block_sascan_speed = (left_block_size / (1024.L * 1024)) / left_block_sascan_time;
  fprintf(stderr, "%.2Lf (%.2LfMiB/s)\n", left_block_sascan_time, left_block_sascan_speed);

  // 2.c
  //
  // Compute the second term for block_initial_ranks.
  if (!last_block) {
    fprintf(stderr, "    Compute initial tail ranks (part 2): ");
    long double initial_ranks_second_term_start = utils::wclock();
    std::vector<long> block_initial_ranks_second_term;
    compute_initial_ranks<block_offset_type>(left_block, left_block_beg, left_block_end, text_length, left_block_psa_ptr,
        text_filename, block_initial_ranks_second_term, max_threads, block_tail_beg, block_tail_end);
    for (size_t j = 0; j < block_initial_ranks_second_term.size(); ++j)
      block_initial_ranks[j] += block_initial_ranks_second_term[j];
    fprintf(stderr, "%.2Lf\n", utils::wclock() - initial_ranks_second_term_start);
  }

  // 2.d
  //
  // Write the partial SA of the left half-block to disk.
  fprintf(stderr, "    Write partial SA to disk: ");
  long double left_psa_save_start = utils::wclock();
  long left_psa_max_part_length = std::max((long)sizeof(block_offset_type), ram_use / 10L);
  info_left.psa = new distributed_file<block_offset_type>(output_filename,
      left_psa_max_part_length, left_block_psa_ptr, left_block_psa_ptr + left_block_size);
  long double left_psa_save_time = utils::wclock() - left_psa_save_start;
  long double left_psa_save_io = ((left_block_size * sizeof(block_offset_type)) / (1024.L * 1024)) / left_psa_save_time;
  fprintf(stderr, "%.2Lf (I/O: %.2LfMiB/s)\n", left_psa_save_time, left_psa_save_io);

  // 2.e
  //
  // Write the BWT of the left half-block to disk.
  unsigned char *left_block_bwt = NULL;
  if (right_block_size > 0) {
    // XXX when I fix the in-memory SAscan, this step will be obsolete.
    fprintf(stderr, "    Copy left BWT to separate array: ");
    long double left_bwt_copy_start = utils::wclock();
    left_block_bwt = (unsigned char *)malloc(left_block_size);
    std::copy(left_block_bwt_ptr, left_block_bwt_ptr + left_block_size, left_block_bwt);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - left_bwt_copy_start);
  }

  // 2.f
  //
  // Write gt_begin of the left half-block to disk.
  // Note that gt_begin of the block is the part of the output of this function.
  // It is only needed if the block is not the first block of the text.
  if (!first_block) {
    fprintf(stderr, "    Write gt_begin to disk: ");
    long double left_gt_begin_save_start = utils::wclock();
    std::string left_block_gt_begin_rev_fname = output_filename + "." + utils::random_string_hash();
    left_block_gt_beg_bv->save_reversed(left_block_gt_begin_rev_fname, left_block_size);
    newtail_gt_begin_rev->add_file(text_length - left_block_end, text_length - left_block_beg, left_block_gt_begin_rev_fname);
    delete left_block_gt_beg_bv;
    long double left_gt_begin_save_time = utils::wclock() - left_gt_begin_save_start;
    long double left_gt_begin_save_io = (left_block_size / (8.L * (1 << 20))) / left_gt_begin_save_time;
    fprintf(stderr, "%.2Lf (I/O: %.2LfMiB/s)\n", left_gt_begin_save_time, left_gt_begin_save_io);
  }


  //----------------------------------------------------------------------------
  // STEP 3: Compute the partial SA of the block.
  //
  // At this point in RAM we still have the left half-block, its partial SA and
  // BWT occupying 7 * left_block_size bytes in total (assuming 5-byte integers).
  //
  // There are 2 cases.
  // I  if the right half-block was empty, we just write the partial suffix
  //    array of the left half-block to disk (as a distributed file).
  // II otherwise, we first we compute the gap array of the left half-block wrt
  //    to the right half-block and then merge the partial suffix arrays of the
  //    half-blocks. Note that the partial SA of the left half-block is already
  //    in memory.
  //----------------------------------------------------------------------------
  if (right_block_size == 0) {  // STEP 3, case I
    hblock_info.push_back(info_left);
    free(left_block);
    free(left_block_sabwt);
    return;
  }


  //----------------------------------------------------------------------------
  // From now on we assume right_block_size > 0.
  //----------------------------------------------------------------------------

  fprintf(stderr, "  Compute partial gap array for left half-block:\n");
  buffered_gap_array *left_block_gap = NULL;

  // STEP 3, case II
    // 3.a
    //
    // Compute initial ranks for streaming of the right half-block.
    fprintf(stderr, "    Compute initial ranks: ");
    long double initial_ranks_right_half_block_start = utils::wclock();
    std::vector<long> initial_ranks2;
    compute_initial_ranks<block_offset_type>(left_block, left_block_beg, left_block_end, text_length,
        left_block_psa_ptr, text_filename, initial_ranks2, max_threads, right_block_beg, right_block_end);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - initial_ranks_right_half_block_start);
    free(left_block);
    free(left_block_sabwt);

    // 3.b
    //
    // Build the rank over BWT of left half-block.
    // RAM: left_block_sabwt, handles to right block psa and gt_begin.
    fprintf(stderr, "    Construct rank: ");
    long double left_block_rank_build_start = utils::wclock();
    rank4n<> *left_block_rank = new rank4n<>(left_block_bwt, left_block_size, max_threads);
    long double left_block_rank_build_time = utils::wclock() - left_block_rank_build_start;
    long double left_block_rank_build_speed = (left_block_size / (1024.L * 1024)) / left_block_rank_build_time;
    fprintf(stderr, "%.2Lf (%.2LfMiB/s)\n", left_block_rank_build_time, left_block_rank_build_speed);

    // 3.c
    //
    // Compute gap array of the left half-block wrt to the right half-block.
    // RAM: left_block_rank, left_block_sabwt, handles to right block psa and gt_begin.
    left_block_gap = new buffered_gap_array(left_block_size + 1, output_filename);
    compute_gap<block_offset_type>(left_block_rank, left_block_gap, right_block_beg, right_block_end,
        text_length, max_threads, left_block_i0, stream_buffer_size, left_block_last,
        initial_ranks2, text_filename, right_block_gt_begin_rev, newtail_gt_begin_rev);
    delete left_block_rank;
    delete right_block_gt_begin_rev;


  //----------------------------------------------------------------------------
  // The computation continues only if the block under consideration is not the
  // last block of the text.
  //----------------------------------------------------------------------------
  if (last_block) {
    free(left_block_bwt);

    // INVARIANT:
    //   The gap of the left half-block wrt to the right half-block is
    //   the gap of the left half-block wrt to the whole tail and right_block_size > 0.
    // What we should do in this situation, is to write the gap to disk
    // and update the information about the gap array filename in info_left.
    info_left.gap_filename = output_filename + ".gap." + utils::random_string_hash();
    left_block_gap->save_to_file(info_left.gap_filename);
    delete left_block_gap;

    hblock_info.push_back(info_left);
    hblock_info.push_back(info_right);
    return;
  }


  //----------------------------------------------------------------------------
  // STEP 4: Compute the BWT for the block. This step is only performed
  //         if the block under consideration is not the last block of text
  //
  // RAM:
  // - handle to block_psa,
  // - left_block_sabwt containing the count array of the left_block_gap (which
  //   is the gap array of the left half-block wrt to the right half-block).
  //
  // DISK:
  // - if the block under consideration is not the last block of the text,
  //   at this point the partial BWT of the left and right half-block is
  //   stored in disk as a regular file and is accesible via filenames
  //   left_block_pbwt_fname and right_block_pbwt_fname. Otherwise, the
  //   files were not created.
  //
  // ADDITIONAL VALUES:
  // - when computing the BWT of the block, we also need to know i0 values
  //   for the BWTs of the left and right half-blocks. This is necessary to
  //     1) replace the occurrence of 0 in the BWT of the right half-block
  //        with left_block_last,
  //     2) identify position i0 for the BWT of the block.
  //
  // INVARIANT:
  // - at this point we have right_block_size > 0. The other case would not
  //   reach this point,
  // - sequential access of the left_block_gap is still initialized (i.e.,
  //   excess values are in RAM).
  //----------------------------------------------------------------------------
  fprintf(stderr, "  Compute block gap:\n");

  // 4.a
  //
  // Convert the partial gap of the left half-block into bitvector.
  fprintf(stderr, "    Convert partial gap of left half-block to bitvector: ");
  long double convert_to_bitvector_start = utils::wclock();
  bitvector *left_block_gap_bv = left_block_gap->convert_to_bitvector(max_threads);
  long double convert_to_bitvector_time = utils::wclock() - convert_to_bitvector_start;
  long double convert_to_bitvector_speed = (block_size / (1024.L * 1024)) / convert_to_bitvector_time;
  fprintf(stderr, "%.2Lf (%.2LfMiB/s)\n", convert_to_bitvector_time, convert_to_bitvector_speed);


  delete left_block_gap;


  // 4.b
  //
  // Read the BWT of the right half-block into RAM.
  fprintf(stderr, "    Read BWT of right half-block: ");
  long double right_block_bwt_read_start = utils::wclock();
  unsigned char *right_block_bwt = (unsigned char *)malloc(right_block_size);
  utils::read_objects_from_file(right_block_bwt, right_block_size, right_block_pbwt_fname);
  long double right_block_bwt_read_time = utils::wclock() - right_block_bwt_read_start;
  long double right_block_bwt_read_io = (right_block_size / (1024.L * 1024)) / right_block_bwt_read_time;
  fprintf(stderr, "%.2Lf (%.2LfMiB/s)\n", right_block_bwt_read_time, right_block_bwt_read_io);


  utils::file_delete(right_block_pbwt_fname);


  unsigned char *block_pbwt = (unsigned char *)malloc(block_size);
  long block_i0 = 0;

  // 4.c
  //
  // Actual merging of BWT.
  fprintf(stderr, "    Merge BWT of half-blocks: ");
  long double bwt_merge_start = utils::wclock();
  block_i0 = merge_bwt(left_block_bwt, right_block_bwt, left_block_size, right_block_size,
      left_block_i0, right_block_i0, left_block_last, block_pbwt, left_block_gap_bv, max_threads);
  long double bwt_merge_time = utils::wclock() - bwt_merge_start;
  long double bwt_merge_speed = (block_size / (1024.L * 1024)) / bwt_merge_time;
  fprintf(stderr, "%.2Lf (%.2LfMiB/s)\n", bwt_merge_time, bwt_merge_speed);

  free(left_block_bwt);
  free(right_block_bwt);


  //----------------------------------------------------------------------------
  // STEP 5: Compute the gap array of the block.
  //
  // RAM:
  // - BWT of the block (block_size bytes)
  // - additional values necessary for streaming:
  //   * last symbol of the block (block_last),
  //   * block i0,
  //   * initial_ranks.
  //----------------------------------------------------------------------------

  // 5.a
  //
  // Construct the rank data structure over BWT of the block.
  fprintf(stderr, "    Construct rank: ");
  long double whole_block_rank_build_start = utils::wclock();
  rank4n<> *block_rank = new rank4n<>(block_pbwt, block_size, max_threads);
  free(block_pbwt);
  long double whole_block_rank_build_time = utils::wclock() - whole_block_rank_build_start;
  long double whole_block_rank_build_io = (block_size / (1024.L * 1024)) / whole_block_rank_build_time;
  fprintf(stderr, "%.2Lf (%.2LfMiB/s)\n", whole_block_rank_build_time, whole_block_rank_build_io);

  buffered_gap_array *block_gap = new buffered_gap_array(block_size + 1, output_filename);

  // 5.b
  //
  // Compute gap for the block. During this step we also compute gt_begin
  // for the new tail (i.e., the block and the old tail).
  // RAM: block_rank, block_gap_array, block_gap.
  compute_gap<block_offset_type>(block_rank, block_gap, block_tail_beg, block_tail_end, text_length,
      max_threads, block_i0, stream_buffer_size, block_last_symbol, block_initial_ranks, text_filename,
      tail_gt_begin_rev, newtail_gt_begin_rev);
  delete block_rank;


  //----------------------------------------------------------------------------
  // STEP 6: compute the gap array of the right half-block
  //
  // At this point we know that right_block_size > 0 and the current block
  // is not the last one in the text. The task now is to compute the gap
  // arrays for the left and right half-blocks from the pseudo gap array for
  // the left half-block and the gap array of the whole block. They are streamed
  // directly to disk and after we're done we update the information about the
  // location of gap arrays into info_left and info_right structures.
  //----------------------------------------------------------------------------
  fprintf(stderr, "  Compute gap array for half-blocks: ");
  long double compute_hb_gaps_start = utils::wclock();
  info_left.gap_filename = output_filename + ".gap." + utils::random_string_hash();
  info_right.gap_filename = output_filename + ".gap." + utils::random_string_hash();
  typedef vbyte_stream_writer<long> vbyte_writer_t;
  vbyte_writer_t *left_gap_writer = new vbyte_writer_t(info_left.gap_filename);
  vbyte_writer_t *right_gap_writer = new vbyte_writer_t(info_right.gap_filename);

  block_gap->start_sequential_access();
  long gap_val = block_gap->get_next();
  long right_half_block_current_gap = 0L;
  long left_half_block_current_gap = 0L;

  for (long i = 0; i < block_size; ++i) {
    // Invariant: val = gap[i].
    right_half_block_current_gap += gap_val;
    left_half_block_current_gap += gap_val;

    if (left_block_gap_bv->get(i)) {
      right_gap_writer->write(right_half_block_current_gap);
      right_half_block_current_gap = 0L;
      ++left_half_block_current_gap;
    } else {
      left_gap_writer->write(left_half_block_current_gap);
      left_half_block_current_gap = 0L;
    }

    gap_val = block_gap->get_next();
  }

  right_half_block_current_gap += gap_val;
  left_half_block_current_gap += gap_val;
  right_gap_writer->write(right_half_block_current_gap);
  left_gap_writer->write(left_half_block_current_gap);

  block_gap->stop_sequential_access();
  delete block_gap;
  delete left_gap_writer;
  delete right_gap_writer;
  delete left_block_gap_bv;

  long double compute_hb_gaps_time = utils::wclock() - compute_hb_gaps_start;
  long double compute_hb_gaps_io = (block_size / (1024.L * 1024)) / compute_hb_gaps_time;
  fprintf(stderr, "%.2Lf (I/O: %.2LfMiB/s)\n", compute_hb_gaps_time, compute_hb_gaps_io);

  hblock_info.push_back(info_left);
  hblock_info.push_back(info_right);
}


//=============================================================================
// Compute partial SAs and gap arrays and write to disk.
// Return the array of handlers to distributed files as a result.
//=============================================================================
template<typename block_offset_type>
std::vector<half_block_info<block_offset_type> > partial_sufsort(std::string text_filename, std::string output_filename,
    long text_length, long max_block_size, long ram_use, long max_threads, long stream_buffer_size) {
  fprintf(stderr, "sizeof(block_offset_type) = %lu\n\n", sizeof(block_offset_type));

  long n_blocks = (text_length + max_block_size - 1) / max_block_size;
  multifile *tail_gt_begin_reversed = NULL;

  std::vector<half_block_info<block_offset_type> > hblock_info;
  for (long block_id = n_blocks - 1; block_id >= 0; --block_id) {
    long block_beg = max_block_size * block_id;
    long block_end = std::min(block_beg + max_block_size, text_length);
    fprintf(stderr, "Processing block %ld/%ld [%ld..%ld):\n", n_blocks - block_id, n_blocks, block_beg, block_end);

    multifile *newtail_gt_begin_reversed = new multifile();
    process_block<block_offset_type>(block_beg, block_end, text_length, ram_use, max_threads, stream_buffer_size,
        text_filename, output_filename, newtail_gt_begin_reversed, tail_gt_begin_reversed, hblock_info);

    delete tail_gt_begin_reversed;
    tail_gt_begin_reversed = newtail_gt_begin_reversed;
  }

  delete tail_gt_begin_reversed;
  return hblock_info;
}

#endif // __PARTIAL_SUFSORT_H_INCLUDED
