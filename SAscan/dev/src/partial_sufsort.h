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

long stream_buffer_size;

template<typename output_type> void SAscan(std::string input_filename, long ram_use, long max_threads);
template<typename output_type> distributed_file<output_type> *partial_SAscan(std::string input_filename,
    bool compute_bwt, long ram_use, unsigned char **BWT, std::string text_filename, long text_offset,
    long max_threads);

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
    long ram_use,
    std::string text_fname,
    std::string sa_fname,
    bool compute_bwt,
    bitvector *gt_eof_bv,
    unsigned char **BWT,
    long block_offset,
    long max_threads) {

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
      bwt_from_sa_replace_text(SA, B, block_size, max_threads);
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
      bwt_from_sa_replace_text(SA, B, block_size, max_threads);
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

    result = partial_SAscan<block_offset_type>(B_fname, compute_bwt, ram_use,
        BWT, text_fname, block_offset, max_threads);

    utils::file_delete(B_fname);
    fprintf(stderr, "  Recursively computing partial SA: %.2Lf\n",
        utils::wclock() - rec_partial_sa_start);
  }

  return result;
}

//=============================================================================
// Compute partial SAs and gap arrays and write to disk.
// Return the array of handlers to distributed files as a result.
//=============================================================================
template<typename block_offset_type>
distributed_file<block_offset_type> **partial_sufsort(std::string filename,
    long length, long max_block_size, long ram_use, long max_threads) {
  long n_stream_buffers = 2 * max_threads;
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

    long right_part_length = length - end;
    long stream_block_size = (right_part_length + max_threads - 1) / max_threads;
    long n_threads = 0L;
    if (end != length)
      n_threads = (right_part_length + stream_block_size - 1) / stream_block_size;
    std::vector<std::string> gt_filenames(n_threads);
    std::vector<long> initial_rank(n_threads);

    long count[256] = {0};
    bitvector *gt_eof_bv = NULL;
    if (need_streaming) {
      //------------------------------------------------------------------------
      // Compute initial ranks for streaming.
      //------------------------------------------------------------------------
      fprintf(stderr, "  [PARALLEL]Computing initial ranks: ");
      long double initial_ranks_start = utils::wclock();
      std::thread **threads = new std::thread*[n_threads];
      for (int t = 0; t < n_threads; ++t) {
        long stream_block_beg = end + t * stream_block_size;
        long stream_block_end = std::min(stream_block_beg + stream_block_size, length);
        threads[t] = new std::thread(parallel_smaller_suffixes, B, block_size, filename,
            stream_block_end, std::ref(initial_rank[t]));
      }
      for (int t = 0; t < n_threads; ++t) threads[t]->join();
      for (int t = 0; t < n_threads; ++t) delete threads[t];
      delete[] threads;
      fprintf(stderr, "%.2Lf\n", utils::wclock() - initial_ranks_start);
 
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

    distrib_files[block_id] = compute_partial_sa_and_bwt<block_offset_type>(
        B, block_size, block_id, ram_use, filename, sa_fname, need_streaming,
        gt_eof_bv, &BWT, beg, max_threads);

    if (need_streaming) {
      // 5a. Build the rank support for BWT.
      fprintf(stderr, "  Building the rank data structure: ");
      long double building_rank_start = utils::wclock();
      rank4n<> *rank = new rank4n<>(BWT, block_size - 1, max_threads);
      delete[] BWT;
      fprintf(stderr, "%.2Lf\n", utils::wclock() - building_rank_start);

      // 5b. Allocate the gap array, do the streaming and store gap to disk.
      fprintf(stderr, "  [PARALLEL]Stream:");
      long double stream_start = utils::wclock();

      buffered_gap_array *gap = new buffered_gap_array(block_size + 1);

      // Allocate buffers.
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
      stream_info info(n_threads, right_part_length);
      std::thread **streamers = new std::thread*[n_threads];
      for (long t = 0L; t < n_threads; ++t) {
        long stream_block_beg = end + t * stream_block_size;
        long stream_block_end = std::min(stream_block_beg + stream_block_size, length);
        streamers[t] = new std::thread(parallel_stream<block_offset_type>,
            full_buffers, empty_buffers, stream_block_beg, stream_block_end,
            initial_rank[t], count, whole_suffix_rank, rank, last, filename, length,
            std::ref(gt_filenames[t]), &info, t, gap->m_length, stream_buffer_size);
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
      utils::file_delete(filename + ".gt");

      long double stream_time = utils::wclock() - stream_start;
      long double speed = ((length - end) / (1024.L * 1024)) / stream_time;
      fprintf(stderr,"\r  [PARALLEL]Stream: 100.0%%. Time: %.2Lf. Threads: %ld. "
          "Speed: %.2LfMiB/s (avg), %.2LfMiB/s (total)\n",
          stream_time, info.m_thread_count, speed / n_threads, speed);

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
      for (long jj = n_threads - 1; jj >= 0; --jj) {
        long stream_block_beg = end + jj * stream_block_size;
        long stream_block_end = std::min(stream_block_beg + stream_block_size, length);
        long this_stream_block_size = stream_block_end - stream_block_beg;

        bit_stream_reader *gt_tail = new bit_stream_reader(gt_filenames[jj]);
        for (long tt = 0; tt < this_stream_block_size; ++tt) gt->write(gt_tail->read());
        delete gt_tail;
        utils::file_delete(gt_filenames[jj]);
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
