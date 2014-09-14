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

long stream_buffer_size;

template<typename output_type> void SAscan(std::string input_filename, long ram_use, long max_threads);
template<typename output_type> distributed_file<output_type> *partial_SAscan(std::string input_filename,
    bool compute_bwt, long ram_use, unsigned char **BWT, std::string text_filename, long text_offset,
    long max_threads);


void compute_initial_ranks(unsigned char *block, long block_beg,
    long block_end, long text_length, std::string text_filename, long stream_max_block_size,
    std::vector<long> &result, multifile *tail_gt_begin_reversed) {
  fprintf(stderr, "  [PARALLEL]Computing initial ranks: ");
  long double start = utils::wclock();

  long tail_length = text_length - block_end;
  long n_threads = (tail_length + stream_max_block_size - 1) / stream_max_block_size;

  result.resize(n_threads);
  std::thread **threads = new std::thread*[n_threads];
  for (int t = 0; t < n_threads; ++t) {
    long stream_block_beg = block_end + t * stream_max_block_size;
    long stream_block_end = std::min(stream_block_beg + stream_max_block_size, text_length);

    threads[t] = new std::thread(parallel_smaller_suffixes, block, block_beg, block_end, 
        text_length, text_filename, stream_block_end, std::ref(result[t]), tail_gt_begin_reversed);
  }

  for (int t = 0; t < n_threads; ++t) threads[t]->join();
  for (int t = 0; t < n_threads; ++t) delete threads[t];
  delete[] threads;
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
}


void rename_block(unsigned char *block, long block_beg, long block_end,
    long text_length, std::string text_filename, multifile *tail_gt_begin_reversed) {
  long block_size = block_end - block_beg;


  //----------------------------------------------------------------------------
  // STEP 1: compute block_gt_end from tail_gt_begin.
  //----------------------------------------------------------------------------
  fprintf(stderr, "  Compute gt_end for block: ");
  long double start = utils::wclock();
  bitvector *block_gt_end = new bitvector(block_size);
  compute_block_gt_end(block, block_beg, block_end, text_length, text_filename, tail_gt_begin_reversed, block_gt_end);
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);


  //----------------------------------------------------------------------------
  // STEP 3: remap the block.
  //----------------------------------------------------------------------------
  fprintf(stderr, "  Rename block: ");
  start = utils::wclock();
  unsigned char block_last = block[block_size - 1];
  for (long j = 0; j + 1 < block_size; ++j)
    if (block[j] > block_last || (block[j] == block_last && block_gt_end->get(j + 1)))
      block[j] += 1;
  block[block_size - 1] += 1;
  delete block_gt_end;
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
}


//=============================================================================
// Compute partial SA of block[0..block_size) and store on disk.
// If block_id != n_block also compute the BWT of block.
//
// INVARIANT: on entry to the function it holds: 5.2 * block_size <= ram_use
//=============================================================================
//
//
// Notes for the future:
//   - this function should be called process_block and there should be
//     exactly specified what it does: computes the partial suffix array
//     of the given block, if requested returns the BWT (in memory)
//     and gt_begin (on disk), again, only if requested.
//     essentially it calles parallel inmem_sascan twice.
//   - gt as a result is not neded when block_beg == 0, i.e., when
//     processing the first block
//   - BWT as a result if not needed for the last block
//     (we will never build a rank from it to stream anything).
//
//==============================================================================
template<typename block_offset_type>
distributed_file<block_offset_type> *compute_partial_sa_and_bwt(
    unsigned char *block,
    long block_beg,
    long block_end,
    long max_block_size,
    long text_length,
    long ram_use,
    std::string text_filename,
    unsigned char **BWT,
    long max_threads,
    long &isa0,
    multifile *newtail_gt_begin_reversed,
    multifile *tail_gt_begin_reversed
    ) {

  long block_id = block_beg / max_block_size;
  long block_size = block_end - block_beg;
  bool compute_bwt = (block_end != text_length);



  fprintf(stderr, "  Compute gt_begin for block: ");
  long double start = utils::wclock();
  bitvector *block_gt_begin_reversed = new bitvector(block_size);

  isa0 = compute_block_gt_begin_reversed(block, block_beg, block_end, text_length, text_filename,
      tail_gt_begin_reversed, block_gt_begin_reversed);

  std::string block_gt_begin_reversed_filename = text_filename + ".block_gt_begin_reversed." + utils::random_string_hash();
  block_gt_begin_reversed->save(block_gt_begin_reversed_filename);
  newtail_gt_begin_reversed->add_file(text_length - block_end, text_length - block_beg, block_gt_begin_reversed_filename);
  delete block_gt_begin_reversed;
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);


 

  // with parallel SAscan we will not do that.
  if (block_end != text_length)
    rename_block(block, block_beg, block_end, text_length, text_filename, tail_gt_begin_reversed);


  std::string sa_fname = text_filename + ".partial_sa." + utils::intToStr(block_id); // XXX base of filename should be SA file.
  distributed_file<block_offset_type> *result = new distributed_file<block_offset_type>(sa_fname.c_str(),
      std::max((long)sizeof(block_offset_type), ram_use / 10L));
  fprintf(stderr, "  compute-bwt = %s\n", compute_bwt ? "TRUE" : "FALSE");

  if (block_size <= MAX_32BIT_DIVSUFSORT_LENGTH) {
    // Easy case, just use use 32-bit divsufsort.
    fprintf(stderr, "  Computing partial SA (divsufsort): ");
    long double sa_start = utils::wclock();
    int *SA = new int[block_size];
    divsufsort(block, SA, (int)block_size);
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
      unsigned char block_last = block[block_size - 1] - 1;
      for (long j = 0; j < block_size; ++j)
        if (block[j] > block_last) block[j] -= 1;
      fprintf(stderr, "%.2Lf\n", utils::wclock() - reremap_start);

      // Compute BWT
      fprintf(stderr, "  Compute BWT: ");
      long double bwt_start = utils::wclock();
      bwt_from_sa_replace_text(SA, block, block_size, max_threads);
      *BWT = block;
      fprintf(stderr, "%.2Lf\n", utils::wclock() - bwt_start);
    } else free(block);
    
    delete[] SA;
  } else if (9L * block_size <= ram_use) {
    // Easy case: block_size >= 2GiB but enough RAM to use divsufsort64.
    fprintf(stderr, "  Computing partial SA (divsufsort64): ");
    long double sa_start = utils::wclock();
    long *SA = new long[block_size];
    divsufsort64(block, SA, block_size);
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
      unsigned char block_last = block[block_size - 1] - 1;
      for (long j = 0; j < block_size; ++j)
        if (block[j] > block_last) block[j] -= 1;
      fprintf(stderr, "%.2Lf\n", utils::wclock() - reremap_start);

      // Compute BWT
      fprintf(stderr, "  Compute BWT: ");
      long double bwt_start = utils::wclock();
      bwt_from_sa_replace_text(SA, block, block_size, max_threads); /// XXX this we also change.
      *BWT = block;
      fprintf(stderr, "%.2Lf\n", utils::wclock() - bwt_start);
    } else delete[] block;
    
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
    std::string B_fname = text_filename + ".block" + utils::intToStr(block_id);
    utils::write_objects_to_file(block, block_size, B_fname);

    // Free all memory.
    delete[] block;

    result = partial_SAscan<block_offset_type>(B_fname, compute_bwt, ram_use,
        BWT, text_filename, block_beg, max_threads);

    utils::file_delete(B_fname);
    fprintf(stderr, "  Recursively computing partial SA: %.2Lf\n",
        utils::wclock() - rec_partial_sa_start);
  }

  return result;
}

template<typename block_offset_type>
buffered_gap_array* compute_gap(
    unsigned char *BWT,
    long block_beg,
    long block_end,
    long text_length,
    std::string text_filename,
    std::vector<long> initial_ranks,
    long max_threads,
    unsigned char block_last_symbol,
    long isa0,
    multifile *tail_gt_begin_reversed,
    multifile *newtail_gt_begin_reversed) {

  long block_size = block_end - block_beg;

  long tail_length = text_length - block_end;
  long stream_max_block_size = (tail_length + max_threads - 1) / max_threads;
  long n_threads = (tail_length + stream_max_block_size - 1) / stream_max_block_size;

  std::vector<std::string> gt_filenames(n_threads);

  // 5a. Build the rank support for BWT.
  fprintf(stderr, "  Building the rank data structure: ");
  long double building_rank_start = utils::wclock();
  rank4n<> *rank = new rank4n<>(BWT, block_size, max_threads);
  free(BWT);
  fprintf(stderr, "%.2Lf\n", utils::wclock() - building_rank_start);


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
  fprintf(stderr, "  [PARALLEL]Stream:");
  long double stream_start = utils::wclock();

  buffered_gap_array *gap = new buffered_gap_array(block_size + 1);

  // Allocate buffers.
  long n_stream_buffers = 2 * max_threads;
  buffer<block_offset_type> **buffers = new buffer<block_offset_type>*[n_stream_buffers];
  for (long i = 0L; i < n_stream_buffers; ++i)
    buffers[i] = new buffer<block_offset_type>(stream_buffer_size);

  // Create poll of empty and full buffers.
  buffer_poll<block_offset_type> *empty_buffers = new buffer_poll<block_offset_type>();
  buffer_poll<block_offset_type> *full_buffers = new buffer_poll<block_offset_type>(n_threads);

  // Add empty buffers to empty poll.
  for (long i = 0L; i < n_stream_buffers; ++i)
    empty_buffers->add(buffers[i]);

  // Start workers.
  stream_info info(n_threads, tail_length);
  std::thread **streamers = new std::thread*[n_threads];
  for (long t = 0L; t < n_threads; ++t) {
    long stream_block_beg = block_end + t * stream_max_block_size;
    long stream_block_end = std::min(stream_block_beg + stream_max_block_size, text_length);
    gt_filenames[t] = text_filename + ".gt_tail." + utils::random_string_hash();

    newtail_gt_begin_reversed->add_file(text_length - stream_block_end,
        text_length - stream_block_beg, gt_filenames[t]);

    streamers[t] = new std::thread(parallel_stream<block_offset_type>,
        full_buffers, empty_buffers, stream_block_beg, stream_block_end,
        initial_ranks[t], count, isa0, rank, block_last_symbol, text_filename, text_length,
        std::ref(gt_filenames[t]), &info, t, gap->m_length, stream_buffer_size, tail_gt_begin_reversed);
  }

  // Start updaters.
  std::thread *updater = new std::thread(gap_updater<block_offset_type>,
        full_buffers, empty_buffers, gap);

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
  delete rank;
  delete[] count;

  long double stream_time = utils::wclock() - stream_start;
  long double speed = ((text_length - block_end) / (1024.L * 1024)) / stream_time;
  fprintf(stderr,"\r  [PARALLEL]Stream: 100.0%%. Time: %.2Lf. Threads: %ld. "
      "Speed: %.2LfMiB/s (avg), %.2LfMiB/s (total)\n",
      stream_time, info.m_thread_count, speed / n_threads, speed);

  // 5c. Save gap to file.
  return gap;
}


//=============================================================================
// Compute partial SAs and gap arrays and write to disk.
// Return the array of handlers to distributed files as a result.
//=============================================================================
template<typename block_offset_type>
distributed_file<block_offset_type> **partial_sufsort(std::string text_filename,
    long text_length, long max_block_size, long ram_use, long max_threads) {
  long n_blocks = (text_length + max_block_size - 1) / max_block_size;
  distributed_file<block_offset_type> **distrib_files = new distributed_file<block_offset_type>*[n_blocks];

  multifile *tail_gt_begin_reversed = NULL;

  for (long block_id = n_blocks - 1; block_id >= 0; --block_id) {
    long block_beg = max_block_size * block_id;
    long block_end = std::min(block_beg + max_block_size, text_length);
    long block_size = block_end - block_beg;

    long tail_length = text_length - block_end;
    long stream_max_block_size = (tail_length + max_threads - 1) / max_threads;
    bool need_streaming = (block_end != text_length);
    multifile *newtail_gt_begin_reversed = new multifile();

    fprintf(stderr, "Processing block %ld/%ld [%ld..%ld):\n", n_blocks - block_id, n_blocks, block_beg, block_end);
    fprintf(stderr, "  need_streaming = %s\n", need_streaming ? "TRUE" : "FALSE");
    fprintf(stderr, "  block_size = %ld (%.2LfMiB)\n", block_size, (long double)block_size / (1 << 20));



    // 1. Read current block.
    fprintf(stderr, "  Reading block: ");
    long double read_start = utils::wclock();
    unsigned char *block = (unsigned char *)malloc(block_size);
    utils::read_block(text_filename, block_beg, block_size, block);
    unsigned char block_last_symbol = block[block_size - 1];
    fprintf(stderr, "%.2Lf\n", utils::wclock() - read_start);




    std::vector<long> initial_ranks;
    if (block_end != text_length)
      compute_initial_ranks(block, block_beg, block_end, text_length, text_filename, stream_max_block_size,
        initial_ranks, tail_gt_begin_reversed);





    unsigned char *BWT = NULL;
    long isa0;
    distrib_files[block_id] = compute_partial_sa_and_bwt<block_offset_type>(block, block_beg,
        block_end, max_block_size, text_length, ram_use, text_filename, &BWT, max_threads, isa0,
        newtail_gt_begin_reversed, tail_gt_begin_reversed);



    if (block_end != text_length) {
      buffered_gap_array *gap = compute_gap<block_offset_type>(BWT, block_beg,
          block_end, text_length, text_filename, initial_ranks,
          max_threads, block_last_symbol, isa0, tail_gt_begin_reversed, newtail_gt_begin_reversed);


      gap->save_to_file(text_filename + ".gap." + utils::intToStr(block_id));
      delete gap;
    }

    delete tail_gt_begin_reversed;
    tail_gt_begin_reversed = newtail_gt_begin_reversed;
  }

  delete tail_gt_begin_reversed;

  return distrib_files;
}

#endif // __PARTIAL_SUFSORT_H_INCLUDED
