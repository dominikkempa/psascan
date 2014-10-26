#ifndef __PARTIAL_SUFSORT_H_INCLUDED
#define __PARTIAL_SUFSORT_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>
#include <algorithm>

#include "divsufsort.h"
#include "divsufsort64.h"
#include "utils.h"
#include "rank.h"
#include "srank.h"
#include "gap_array.h"
#include "merge.h"
#include "multifile_bitvector.h"
#include "bitvector.h"
#include "stream.h"
#include "sascan.h"


//==============================================================================
// Compute partial SA of B[0..block_size) and store on disk.
//==============================================================================
template<typename block_offset_type>
distributed_file<block_offset_type> *compute_partial_sa_and_bwt(
    unsigned char *B,
    long block_size,
    long ram_use,
    std::string sa_fname,
    bool compute_bwt,
    unsigned char **BWT,
    long &whole_suffix_rank) {

  distributed_file<block_offset_type> *result = new distributed_file<block_offset_type>(sa_fname.c_str(),
      std::max((long)sizeof(block_offset_type), ram_use / 10L));

  if (block_size < (1L << 31)) {
    // Easy case, just use use 32-bit divsufsort.
    fprintf(stderr, "  Compute partial sa (divsufsort32): ");
    long double sa_start = utils::wclock();
    int *SA = new int[block_size];
    divsufsort(B, SA, (int)block_size);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - sa_start);

    fprintf(stderr, "  Write partial sa to disk: ");
    long double writing_sa_start = utils::wclock();
    result->initialize_writing(4 << 20);
    for (long i = 0; i < block_size; ++i)
      result->write((block_offset_type)SA[i]);
    result->finish_writing();
    fprintf(stderr, "%.2Lf\n", utils::wclock() - writing_sa_start);

    if (compute_bwt) {
      // Remap symbols of B back to original.
      fprintf(stderr, "  Reremap block: ");
      long double reremap_start = utils::wclock();
      unsigned char block_last = B[block_size - 1] - 1;
      for (long j = 0; j < block_size; ++j)
        if (B[j] > block_last) B[j] -= 1;
      fprintf(stderr, "%.2Lf\n", utils::wclock() - reremap_start);

      // Compute BWT
      fprintf(stderr, "  Compute bwt: ");
      long double bwt_start = utils::wclock();
      unsigned char *tmp = (unsigned char *)SA;
      for (long j = 0; j < block_size; ++j) {
        if (SA[j]) tmp[j] = B[SA[j] - 1];
        else {
          tmp[j] = 0;
          whole_suffix_rank = j;
        }
      }
      std::copy(tmp, tmp + block_size, B);
      *BWT = B;
      fprintf(stderr, "%.2Lf\n", utils::wclock() - bwt_start);
    } else delete[] B;
    
    delete[] SA;
  } else {
    // Easy case: block_size >= 2GiB but enough RAM to use divsufsort64.
    fprintf(stderr, "  Compute partial sa (divsufsort64): ");
    long double sa_start = utils::wclock();
    long *SA = new long[block_size];
    divsufsort64(B, SA, block_size);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - sa_start);
    
    fprintf(stderr, "  Write partial sa to disk: ");
    long double writing_sa_start = utils::wclock();
    result->initialize_writing(4 << 20);
    for (long i = 0; i < block_size; ++i)
      result->write((block_offset_type)SA[i]);
    result->finish_writing();
    fprintf(stderr, "%.2Lf\n", utils::wclock() - writing_sa_start);

    if (compute_bwt) {
      // Remap symbols of B back to original.
      fprintf(stderr, "  Reremap B: ");
      long double reremap_start = utils::wclock();
      unsigned char block_last = B[block_size - 1] - 1;
      for (long j = 0; j < block_size; ++j)
        if (B[j] > block_last) B[j] -= 1;
      fprintf(stderr, "%.2Lf\n", utils::wclock() - reremap_start);

      // Compute BWT
      fprintf(stderr, "  Compute BWT: ");
      long double bwt_start = utils::wclock();
      unsigned char *tmp = (unsigned char *)SA;
      for (long j = 0; j < block_size; ++j) {
        if (SA[j]) tmp[j] = B[SA[j] - 1];
        else{
          tmp[j] = 0;
          whole_suffix_rank = j;
        }
      }
      std::copy(tmp, tmp + block_size, B);
      *BWT = B;
      fprintf(stderr, "%.2Lf\n", utils::wclock() - bwt_start);
    } else  delete[] B;
    delete[] SA;
  }

  return result;
}

template<typename block_offset_type>
void compute_gap(
    rank4n<> *rank,
    buffered_gap_array *gap,
    long tail_begin,
    long tail_end,
    long text_length,
    std::string text_filename,
    std::string output_filename,
    long i,
    unsigned char block_last_symbol,
    long block_isa0,
    multifile *tail_gt_begin_reversed,
    multifile *newtail_gt_begin_reversed) {
  long tail_length = tail_end - tail_begin;

  // Obtain symbol counts from the rank data structure.
  long count[256] = {0};
  std::copy(rank->m_count, rank->m_count + 256, count);
  count[block_last_symbol]++;
  count[0]--;
  if (count[255]) {
    fprintf(stderr, "Error: input cannot contain byte 255\n");
    std::exit(EXIT_FAILURE);
  }
  for (long j = 0, s = 0, t; j < 256; ++j) {
    t = count[j];
    count[j] = s;
    s += t;
  }

  fprintf(stderr, "  Stream:\r");
  long double stream_start = utils::wclock();

  // Update the information about the bitvector on disk.
  std::string newtail_gt_begin_tail_filename = output_filename + ".gt." + utils::random_string_hash();
  newtail_gt_begin_reversed->add_file(text_length - tail_end, text_length - tail_begin, newtail_gt_begin_tail_filename);
  bit_stream_writer *gt_out = new bit_stream_writer(newtail_gt_begin_tail_filename);

  // Initialize reading of old bitvectors.
  multifile_bitvector_reader gt_in(tail_gt_begin_reversed);
  gt_in.initialize_sequential_reading(0L);

  backward_stream_reader<unsigned char> *text_streamer = new backward_stream_reader<unsigned char>(text_filename);

  for (long j = tail_end, dbg = 0; j > tail_begin; --j, ++dbg) {
    if (dbg == (1 << 24)) {
      long double elapsed = utils::wclock() - stream_start;
      long double streamed_mib = (1.L * (tail_end - j)) / (1 << 20);
      fprintf(stderr,"  Stream: %.1Lf%%. Time: %.2Lf. Speed: %.2LfMiB/s\r",
        (100.L * (tail_end - j)) / tail_length,
        elapsed,
        streamed_mib / elapsed);
      dbg = 0;
    }

    unsigned char c = text_streamer->read();
    unsigned char next_gt = gt_in.read();
    gt_out->write(i > block_isa0);

    int delta = (i > block_isa0 && c == 0);
    i = count[c] + rank->rank(i, c) - delta;
    if (c == block_last_symbol && next_gt) ++i;
    gap->increment(i);
  }

  long double stream_time = utils::wclock() - stream_start;
  long double streamed_mib = (1.L * tail_length) / (1 << 20);
  fprintf(stderr,"  Stream: 100.0%%. Time: %.2Lf. Speed: %.2LfMiB/s\n",
      stream_time, streamed_mib / stream_time);

  delete text_streamer;
  delete gt_out;
}

template<typename block_offset_type>
distributed_file<block_offset_type> *process_block(long block_beg,
    long block_end, long max_block_size, long text_length,
    std::string text_filename, std::string output_filename, long ram_use,
    multifile *newtail_gt_begin_rev, multifile *tail_gt_begin_rev) {
  long block_size = block_end - block_beg;
  long block_id = block_beg / max_block_size;
  bool last_block = (block_end == text_length);
  bool first_block = (block_beg == 0);

  // Read current block.
  fprintf(stderr, "  Read block: ");
  long double start = utils::wclock();
  unsigned char *block = new unsigned char[block_size];
  utils::read_block(text_filename, block_beg, block_size, block);
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
  unsigned char block_last = block[block_size - 1];

  bitvector *bv = new bitvector(block_size);
  if (!last_block) {
    // Compute block_sm_end.
    fprintf(stderr, "  Compute block sm_end: ");
    long double gt_eof_start = utils::wclock();
    bitvector *block_sm_end = bv;
    compute_sm_end(block, block_beg, block_end, text_length, text_filename, tail_gt_begin_rev, block_sm_end);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - gt_eof_start);

    // Remap symbols of the block.
    fprintf(stderr, "  Remap block: ");
    long double remap_start = utils::wclock();
    for (long j = 0; j + 1 < block_size; ++j)
      if (block[j] > block_last || (block[j] == block_last && !block_sm_end->get(j + 1)))
        ++block[j];;
    ++block[block_size - 1];
    fprintf(stderr, "%.2Lf\n", utils::wclock() - remap_start);
  }

  // Compute the head of the new gt bitvector.
  if (!first_block) {
    fprintf(stderr, "  Compute block gt_begin: ");
    start = utils::wclock();
    bitvector *block_gt_begin_rev = bv;
    transform_sm_end_into_gt_begin_reversed(block, block_size, block_gt_begin_rev);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

    fprintf(stderr, "  Write block gt_begin to disk: ");
    start = utils::wclock();
    std::string newtail_gt_begin_head_filename = output_filename + ".gt." + utils::random_string_hash();
    block_gt_begin_rev->save(newtail_gt_begin_head_filename);
    newtail_gt_begin_rev->add_file(text_length - block_end, text_length - block_beg, newtail_gt_begin_head_filename);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
  }

  delete bv;

  // Compute and save partial SA.
  unsigned char *bwt = NULL;
  std::string sa_fname = output_filename + ".partial_sa." + utils::intToStr(block_id);
  long whole_suffix_rank = 0;

  distributed_file<block_offset_type> *result = compute_partial_sa_and_bwt<block_offset_type>(
      block, block_size, ram_use,  sa_fname, !last_block, &bwt, whole_suffix_rank);

  if (!last_block) {
    // Build the rank support for BWT.
    fprintf(stderr, "  Build rank over bwt: ");
    start = utils::wclock();
    rank4n<> *rank = new rank4n<>(bwt, block_size);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
    delete[] bwt;

    // Allocate and compute the gap array.
    std::string gap_filename = output_filename + ".excess";
    buffered_gap_array *gap = new buffered_gap_array(block_size + 1, gap_filename);
    compute_gap<block_offset_type>(rank, gap, block_end, text_length, text_length, text_filename,
        output_filename, 0, block_last, whole_suffix_rank, tail_gt_begin_rev, newtail_gt_begin_rev);

    delete rank;

    // Save gap to file.
    gap->save_to_file(output_filename + ".gap." + utils::intToStr(block_id));
    delete gap;
  }

  return result;
}


//==============================================================================
// Compute partial SAs and gap arrays and write to disk.
// Return the array of handlers to distributed files as a result.
//==============================================================================
template<typename block_offset_type>
distributed_file<block_offset_type> **partial_sufsort(
    std::string text_filename, long text_length, std::string output_filename,
    long max_block_size, long ram_use) {
  fprintf(stderr, "sizeof(block_offset_type) = %lu\n\n", sizeof(block_offset_type));

  long n_blocks = (text_length + max_block_size - 1) / max_block_size;
  distributed_file<block_offset_type> **distrib_files = new distributed_file<block_offset_type>*[n_blocks];
  multifile *tail_gt_begin_rev = NULL;

  for (long block_id = n_blocks - 1; block_id >= 0; --block_id) {
    long block_beg = max_block_size * block_id;
    long block_end = std::min(text_length, block_beg + max_block_size);
    fprintf(stderr, "Process block %ld/%ld [%ld..%ld):\n", n_blocks - block_id, n_blocks, block_beg, block_end);

    multifile *newtail_gt_begin_rev = new multifile();
    distrib_files[block_id] = process_block<block_offset_type>(block_beg, block_end, max_block_size,
        text_length, text_filename, output_filename, ram_use, newtail_gt_begin_rev, tail_gt_begin_rev);

    delete tail_gt_begin_rev;
    tail_gt_begin_rev = newtail_gt_begin_rev;
  }

  delete tail_gt_begin_rev;
  return distrib_files;
}

#endif // __PARTIAL_SUFSORT_H_INCLUDED
