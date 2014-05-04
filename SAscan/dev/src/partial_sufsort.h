#ifndef __PARTIAL_SUFSORT_H_INCLUDED
#define __PARTIAL_SUFSORT_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <iostream>
#include <string>
#include <thread>
#include <mutex>
#include <algorithm>

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

std::mutex stdout_mutex;
std::mutex gap_mutex;

template<typename output_type> void SAscan(std::string input_filename, long ram_use);
template<typename output_type> distributed_file<output_type> *partial_SAscan(std::string input_filename,
  long ram_use, unsigned char **BWT, std::string text_filename, long text_offset);

// Compute partial SA of B[0..block_size) and store on disk.
// If block_id != n_block also compute the BWT of B.
//
// INVARIANT: on entry to the function it holds: 5 * block_size <= ram_use
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
    // (block_size >= 2GiB and 9*block_size > ram_use) => use recursion.
    // To save I/O from recursion we also get the BWT, which is obtained
    // during the storage of the resulting partial SA to disk.
    fprintf(stderr, "  Recursively computing partial SA:\n");
    long double rec_partial_sa_start = utils::wclock();

    // Save the remapped block to temp file.
    std::string B_fname = text_fname + ".block" + utils::intToStr(block_id);
    utils::write_objects_to_file(B, block_size, B_fname);

    // Free all memory.
    delete gt_eof_bv;
    delete[] B;

    result = partial_SAscan<block_offset_type>(B_fname, ram_use, BWT, text_fname, block_offset);

    utils::file_delete(B_fname);
    fprintf(stderr, "  Recursively computing partial SA: %.2Lf\n",
        utils::wclock() - rec_partial_sa_start);
  }

  return result;
}

/*********************** CODE EXECUTED IN PARALLEL *****************************/
void parallel_stream(long j_beg, long j_end, long i, buffered_gap_array *gap,
    long *count, long whole_suffix_rank, context_rank_4n *rank,
    unsigned char last, std::string text_filename, long length,
    std::string *tail_gt_filename) {
  long double start = utils::wclock();
  stdout_mutex.lock();
  std::cerr << "Thread " << std::this_thread::get_id() << " begins.\n";
  std::cerr << "Thread " << std::this_thread::get_id() << " j-range: [" << j_beg << ", " << j_end << "]\n";
  stdout_mutex.unlock();

  static const int gap_buf_size = (1 << 19);
  long *gap_buf = new long[gap_buf_size];
  int gap_buf_filled = 0;

  gt_accessor *gt_in = new gt_accessor(text_filename + std::string(".gt"));
  bool next_gt = (j_end + 1 == length) ? 0 : (*gt_in)[length - j_end - 2];
  backward_skip_stream_reader<unsigned char> *text_streamer
    = new backward_skip_stream_reader<unsigned char>(text_filename, length - 1 - j_end);

  *tail_gt_filename = text_filename + std::string(".gt_tail.") + utils::random_string_hash();
  bit_stream_writer *gt_out = new bit_stream_writer(*tail_gt_filename);

  for (long j = j_end, dbg = 0; j > j_beg; --j, ++dbg) {
    if (dbg == (1 << 26)) {
      long double elapsed = utils::wclock() - start;
      long double streamed_mib = (1.L * (j_end - j)) / (1 << 20);
      stdout_mutex.lock();
      std::cerr << "Thread " << std::this_thread::get_id() << ": ";
      fprintf(stderr,"  stream: %.1Lf%%. Time: %.2Lf. Speed: %.2LfMiB/s\n",
        (100.L * (j_end - j)) / (j_end - j_beg), elapsed, streamed_mib / elapsed);
      stdout_mutex.unlock();
      dbg = 0;
    }

    unsigned char c = text_streamer->read();
    i = count[c] + rank->rank(i - (i > whole_suffix_rank), c);
    if (c == last && next_gt) ++i;
    gt_out->write(i > whole_suffix_rank);
    next_gt = (*gt_in)[length - j - 1];

    gap_buf[gap_buf_filled++] = i;
    if (gap_buf_filled == gap_buf_size) {
      gap_mutex.lock();
      gap->increment(gap_buf, gap_buf_filled);
      gap_mutex.unlock();
      gap_buf_filled = 0;    
    }
  }

  if (gap_buf_filled) {
    gap_mutex.lock();
    gap->increment(gap_buf, gap_buf_filled);
    gap_mutex.unlock();
  }
  
  delete[] gap_buf;
  delete text_streamer;

  delete gt_in;
  delete gt_out;

  long double elapsed = utils::wclock() - start;
  stdout_mutex.lock();
  std::cerr << "Thread " << std::this_thread::get_id() << " finished in "
    << elapsed << "secs.\n";
  stdout_mutex.unlock();
}
/********************* END OF CODE EXECUTED IN PARALLEL ************************/


// Compute partial SAs and gap arrays and write to disk.
// Return the array of handlers to distributed files as a result.
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

    // 1. Read current block.
    fprintf(stderr, "  Reading block: ");
    long double read_start = utils::wclock();
    unsigned char *B = new unsigned char[block_size];
    utils::read_block(filename, beg, block_size, B);
    unsigned char last = B[block_size - 1];
    fprintf(stderr, "%.2Lf\n", utils::wclock() - read_start);

#define MAX_THREADS 16L
    int starting_points = std::min(MAX_THREADS, length - end);
    long *start_j = NULL, *start_i = NULL;
    std::string *tail_gt_filenames = NULL;

    long count[256] = {0};
    bitvector *gt_eof_bv = NULL;
    if (need_streaming) {
      //-----------------------------------------------------------------------
      fprintf(stderr, "  Computing starting points: ");
      long double starting_points_start = utils::wclock();
      start_j = new long[starting_points];
      start_i = new long[starting_points];
      tail_gt_filenames = new std::string[starting_points];

      // Evenly distribute the starting points (start_j).
      long parallel_block = (length - end) / starting_points;
      long prev_start_j = end - 1;
      for (int jj = 0; jj < starting_points; ++jj) {
        start_j[jj] = std::min(prev_start_j + parallel_block, length - 1);
        prev_start_j = start_j[jj];
      }
      start_j[starting_points - 1] = length - 1;

      // Run string range matching in parallel
      std::thread **threads = new std::thread*[starting_points];
      for (int jj = 0; jj < starting_points; ++jj)
        threads[jj] = new std::thread(parallel_smaller_suffixes,
            B, block_size, filename, start_j[jj] + 1, start_i + jj);
      for (int jj = 0; jj < starting_points; ++jj) threads[jj]->join();
      for (int jj = 0; jj < starting_points; ++jj) delete threads[jj];
      delete[] threads;
      fprintf(stderr, "%.2Lf\n", utils::wclock() - starting_points_start);
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
      bitvector *gt_head_bv = new bitvector(filename + ".gt_head");
      compute_gt_eof_bv(extprevB, ext_prev_block_size, B, block_size, gt_head_bv, gt_eof_bv);
      delete gt_head_bv;
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

    distrib_files[block_id] = compute_partial_sa_and_bwt<block_offset_type>
      (B, block_size, block_id, ram_use, filename, sa_fname, need_streaming, gt_eof_bv, &BWT, beg);

    if (need_streaming) {
      // 5a. Build the rank support for BWT.
      fprintf(stderr, "  Building the rank data structure: ");
      long double building_rank_start = utils::wclock();
      context_rank_4n *rank = new context_rank_4n(BWT, block_size - 1);
      delete[] BWT;
      fprintf(stderr, "%.2Lf\n", utils::wclock() - building_rank_start);

      // 5b. Allocate the gap array, do the streaming and store gap to disk.
      //------------------------------------------------------------------------
      fprintf(stderr, "Stream: \r");
      long double stream_start = utils::wclock();
      buffered_gap_array *gap = new buffered_gap_array(block_size + 1);
      std::thread **threads = new std::thread*[starting_points];
      for (int jj = 0; jj < starting_points; ++jj) {
        long j_beg = (jj > 0) ? start_j[jj - 1] : end - 1;
        long j_end = start_j[jj];
        threads[jj] = new std::thread(parallel_stream,
            j_beg, j_end, start_i[jj], gap, count, whole_suffix_rank,
            rank, last, filename, length, tail_gt_filenames + jj);
      }
      for (int jj = 0; jj < starting_points; ++jj) threads[jj]->join();
      for (int jj = 0; jj < starting_points; ++jj) delete threads[jj];
      delete[] threads;
      delete rank;
      long double stream_time = utils::wclock() - stream_start;
      long double streamed_mib = (1.L * (length - end)) / (1 << 20);
      fprintf(stderr,"  Stream: 100.0%%. Time: %.2Lf. Speed: %.2LfMiB/s\n",
        stream_time, streamed_mib / stream_time);

      // 5c. Save gap to file.
      gap->save_to_file(filename + ".gap." + utils::intToStr(block_id));
      delete gap;
    }

    // Concatenate all bitvectors into one (but leave filename.gt_head, for now)
    //-------------------------------------------------------------------------
    fprintf(stderr, "  Concatenating gt bitvectors: ");
    long double gt_concat_start = utils::wclock();
    bit_stream_writer *gt = new bit_stream_writer(filename + ".gt");
    if (need_streaming) {
      for (long jj = starting_points - 1; jj >= 0; --jj) {
        long j_beg = (jj > 0) ? start_j[jj - 1] : end - 1, j_end = start_j[jj];
        bit_stream_reader *gt_tail = new bit_stream_reader(tail_gt_filenames[jj]);
        for (long tt = 0; tt < j_end - j_beg; ++tt) gt->write(gt_tail->read());
        delete gt_tail;
        utils::file_delete(tail_gt_filenames[jj]);
      }
      delete[] start_i;
      delete[] start_j;
      delete[] tail_gt_filenames;
    }
    bit_stream_reader *gt_head = new bit_stream_reader(filename + ".gt_head");
    for (long tt = end - 1; tt >= beg; --tt) gt->write(gt_head->read());
    delete gt_head;
    delete gt;
    fprintf(stderr, "%.2Lf\n", utils::wclock() - gt_concat_start);
    //--------------------------------------------------------------------------

    prev_end = end;
    end = beg;
    --block_id;
  }

  if (utils::file_exists(filename + ".gt_head")) utils::file_delete(filename + ".gt_head");
  if (utils::file_exists(filename + ".gt")) utils::file_delete(filename + ".gt");

  return distrib_files;
}

#endif // __PARTIAL_SUFSORT_H_INCLUDED
