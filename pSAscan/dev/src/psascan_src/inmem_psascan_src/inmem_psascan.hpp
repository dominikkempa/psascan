/**
 * @file    src/psascan_src/inmem_psascan_src/inmem_psascan.hpp
 * @section LICENCE
 *
 * This file is part of pSAscan v0.2.0
 * See: http://www.cs.helsinki.fi/group/pads/
 *
 * Copyright (C) 2014-2017
 *   Juha Karkkainen <juha.karkkainen (at) cs.helsinki.fi>
 *   Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/

#ifndef __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_PSASCAN_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_PSASCAN_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>

#include "../bitvector.hpp"
#include "../io/multifile.hpp"
#include "../io/background_block_reader.hpp"
#include "inmem_gap_array.hpp"
#include "compute_initial_gt_bitvectors.hpp"
#include "initial_partial_sufsort.hpp"
#include "change_gt_reference_point.hpp"
#include "inmem_bwt_from_sa.hpp"
#include "inmem_compute_initial_ranks.hpp"
#include "parallel_merge.hpp"
#include "inmem_bwtsa_merge.hpp"
#include "pagearray.hpp"
#include "bwtsa.hpp"
#include "parallel_shrink.hpp"
#include "merge_schedule.hpp"


namespace psascan_private {
namespace inmem_psascan_private {

//==============================================================================
// We assume sa is already allocated.
//
// What should be done is some kind of optimization to the initial gt
// bitvectors computation, so that the last bitvector is not computed, and
// also that the last block is not renamed and so on. The end result should
// be, that when used with one thread, the procedure simply runs divsufsort
// and there is no overhead of running this function over running divsufsort.
//==============================================================================
template<typename text_offset_type, unsigned pagesize_log = 12>
void inmem_psascan(
    std::uint8_t *text,
    std::uint64_t text_length,
    std::uint8_t *sa_bwt,
    std::uint64_t max_threads = 1,
    bool compute_bwt = false,
    bool compute_gt_begin = false,
    bitvector *gt_begin = NULL,
    std::uint64_t max_blocks = 0,
    std::uint64_t text_beg = 0,
    std::uint64_t text_end = 0,
    std::uint64_t supertext_length = 0,
    std::string supertext_filename = "",
    const multifile *tail_gt_begin_reversed = NULL,
    std::uint64_t *i0 = NULL,
    std::uint8_t *tail_prefix_preread = NULL) {
  static const std::uint32_t pagesize = (1U << pagesize_log);
  long double absolute_start = utils::wclock();
  long double start;

  if ((std::uint64_t)std::numeric_limits<text_offset_type>::max() < text_length) {
    fprintf(stderr, "\n\nError: text is too long (%lu bytes),\n", text_length);
    fprintf(stderr, "       std::numeric_limits<text_offset_type>::max() = %lu\n",
        (std::uint64_t)std::numeric_limits<text_offset_type>::max());
    std::exit(EXIT_FAILURE);
  }

  if (max_blocks == 0)
    max_blocks = max_threads;

  if (text_end == 0) {
    supertext_length = text_length;
    text_end = text_length;
    text_beg = 0;
    supertext_filename = "";
    tail_gt_begin_reversed = NULL;
  }

  bool has_tail = (text_end != supertext_length);

  if (!has_tail && tail_prefix_preread != NULL) {
    fprintf(stderr, "\n\nError: has_tail == false but tail_prefix_preread != NULL\n");
    std::exit(EXIT_FAILURE);
  }

  // long max_block_size = (text_length + max_blocks - 1) / max_blocks;
  // while ((max_block_size & 7) || (max_block_size & pagesize_mask)) ++max_block_size;
  // long n_blocks = (text_length + max_block_size - 1) / max_block_size;

  //----------------------------------------------------------------------------
  // min_block_size must be a multiple ot alignment unit. Alignement unit is to
  // simplify the computation involving bitvectors and page arrays. Note: it
  // may happen then min_block_size > text_length. This is perfectly fine due
  // to the way we compute block boundaries (always separate if for the last
  // block).
  //----------------------------------------------------------------------------

  std::uint64_t alignment_unit = std::max(pagesize, 8U);
  std::uint64_t max_block_size = (text_length + max_blocks - 1) / max_blocks;
  while ((max_block_size & (alignment_unit - 1)) && max_block_size < text_length)
    ++max_block_size;

  std::uint64_t n_blocks = (text_length + max_block_size - 1) / max_block_size;

  if (!compute_gt_begin) {
    if (gt_begin) {
      fprintf(stderr, "\n\nError: check gt_begin == NULL failed\n");
      std::exit(EXIT_FAILURE);
    }
    if (n_blocks > 1 || has_tail)
      gt_begin = new bitvector(text_length);
  } else {
    if (!gt_begin) {
      fprintf(stderr, "inmem_sascan: gt_begin was requested but is not allocated!\n");
      std::exit(EXIT_FAILURE);
    }
  }

  fprintf(stderr, "Text length = %lu (%.2LfMiB)\n", text_length, text_length / (1024.L * 1024));
  fprintf(stderr, "Max block size = %lu (%.2LfMiB)\n", max_block_size, max_block_size / (1024.L * 1024));
  fprintf(stderr, "Max blocks = %lu\n", max_blocks);
  fprintf(stderr, "Number of blocks = %ld\n", n_blocks);
  fprintf(stderr, "Max threads = %lu\n", max_threads);
  fprintf(stderr, "sizeof(text_offset_type) = %lu\n", sizeof(text_offset_type));
  fprintf(stderr, "Pagesize = %u\n", (1U << pagesize_log));
  fprintf(stderr, "Compute bwt = %s\n", compute_bwt ? "true" : "false");
  fprintf(stderr, "Compute gt begin = %s\n", compute_gt_begin ? "true" : "false");
  fprintf(stderr, "Text beg = %lu\n", text_beg);
  fprintf(stderr, "Text end = %lu\n", text_end);
  fprintf(stderr, "Supertext length = %lu (%.2LfMiB)\n", supertext_length, supertext_length / (1024.L * 1024));
  fprintf(stderr, "Supertext filename = %s\n", supertext_filename.c_str());
  fprintf(stderr, "Has tail = %s\n", has_tail ? "true" : "false");
  fprintf(stderr, "\n");

  bwtsa_t<text_offset_type> *bwtsa = (bwtsa_t<text_offset_type> *)sa_bwt;

  // Initialize reading of the tail prefix in the background.
  std::uint64_t tail_length = supertext_length - text_end;
  std::uint64_t tail_prefix_length = std::min(text_length, tail_length);

  background_block_reader *tail_prefix_background_reader = NULL;
  if (has_tail && tail_prefix_preread == NULL)
    tail_prefix_background_reader =
      new background_block_reader(supertext_filename, text_end, tail_prefix_length);

  //----------------------------------------------------------------------------
  // STEP 1: compute initial bitvectors, and partial suffix arrays.
  //----------------------------------------------------------------------------
  if (n_blocks > 1 || compute_gt_begin || has_tail) {
    fprintf(stderr, "Compute initial bitvectors:\n");
    start = utils::wclock();
    compute_initial_gt_bitvectors(text, text_length, gt_begin, max_block_size,
        max_threads, text_end, supertext_length, tail_gt_begin_reversed,
        tail_prefix_background_reader, tail_prefix_preread);
    fprintf(stderr, "Total time: %.2Lfs\n\n", utils::wclock() - start);
  }

  fprintf(stderr, "Initial sufsort:\n");
  start = utils::wclock();
  initial_partial_sufsort(text, text_length, gt_begin, bwtsa, max_block_size, max_threads, has_tail);
  fprintf(stderr, "Total time: %.2Lfs\n\n", utils::wclock() - start);

  //----------------------------------------------------------------------------
  // STEP 2: compute matrix of block ranks.
  //----------------------------------------------------------------------------
  fprintf(stderr, "Compute matrix of initial ranks: ");
  start = utils::wclock();
  std::uint64_t **block_rank_matrix = new std::uint64_t*[n_blocks];
  for (std::uint64_t j = 0; j < n_blocks; ++j)
    block_rank_matrix[j] = new std::uint64_t[n_blocks];
  compute_block_rank_matrix<text_offset_type>(text, text_length, bwtsa,
      max_block_size, text_beg, supertext_length, supertext_filename,
      tail_gt_begin_reversed, tail_prefix_background_reader,
      tail_prefix_preread, block_rank_matrix);

  // Stop reading next block in the background or free memory taken by next block.
  if (has_tail) {
    if (tail_prefix_background_reader != NULL) {
      tail_prefix_background_reader->stop();
      delete tail_prefix_background_reader;
    } else free(tail_prefix_preread);
  }

  fprintf(stderr, "%.2Lfs\n\n", utils::wclock() - start);

  //----------------------------------------------------------------------------
  // STEP 3: compute the gt bitvectors for blocks that will be on the right
  //         side during the merging.
  //----------------------------------------------------------------------------
  if (n_blocks > 1 || compute_gt_begin) {
    fprintf(stderr, "Overwrite gt_end with gt_begin: ");
    start = utils::wclock();
    gt_end_to_gt_begin(text, text_length, gt_begin, max_block_size);
    fprintf(stderr, "%.2Lfs\n\n", utils::wclock() - start);
  }

  float rl_ratio = 10.L;  // estimated empirically
  // Note that 9n for the 32-bit version and 10n for 40-bit version are the most reasonable
  // space usages we can get. In the worst case there are two blocks, thus during the
  // merging the rank + gap array for the left block will take 2.5n. This added to the 7.125n
  // (for 40-bit) and 6.125n (for 32-bit) gives 9.6125n and 8.625n space usages.
  std::uint64_t max_ram_usage_per_input_byte = 10L;  // peak ram usage = 10n
  std::uint32_t max_left_size = (std::uint32_t)std::max(1,
      (std::int32_t)floor(n_blocks * (((long double)max_ram_usage_per_input_byte - (2.125L + sizeof(text_offset_type))) / 5.L)));
  fprintf(stderr, "Assumed rl_ratio: %.2f\n", rl_ratio);
  fprintf(stderr, "Max left size = %u\n", max_left_size);
  fprintf(stderr, "Peak memory usage during last merge = %.3Lfn\n",
      (2.125L + sizeof(text_offset_type)) + (5.L * max_left_size) / n_blocks);
  MergeSchedule schedule(n_blocks, rl_ratio, max_left_size);

  fprintf(stderr, "Skewed merge schedule:\n");
  print_schedule(schedule, n_blocks);
  fprintf(stderr, "\n");

  std::int64_t *i0_array = new std::int64_t[n_blocks];
  if (n_blocks > 1 || compute_bwt) {
    for (std::uint64_t block_id = 0; block_id < n_blocks; ++block_id) {
      std::uint64_t block_end = text_length - (n_blocks - 1 - block_id) * max_block_size;
      std::uint64_t block_beg = (std::uint64_t)std::max(0L, (std::int64_t)block_end - (std::int64_t)max_block_size);
      std::uint64_t block_size = block_end - block_beg;

      if (block_id + 1 != n_blocks || compute_bwt) {
        fprintf(stderr, "Compute BWT for block %ld: ", block_id + 1);
        long double bwt_start = utils::wclock();
        compute_bwt_in_bwtsa<text_offset_type>(text + block_beg, block_size,
            bwtsa + block_beg, max_threads, i0_array[block_id]);
        fprintf(stderr, "%.2Lfs\n", utils::wclock() - bwt_start);
      }
    }
    fprintf(stderr, "\n");
  }

  if (n_blocks > 1) {
    std::uint64_t i0_result;
    pagearray<bwtsa_t<text_offset_type>, pagesize_log> *result =
      inmem_bwtsa_merge<text_offset_type, pagesize_log>(text, text_length, bwtsa,
          gt_begin, max_block_size, 0, n_blocks, max_threads, compute_gt_begin,
          compute_bwt, i0_result, schedule, text_beg, text_end,
          supertext_length, supertext_filename, tail_gt_begin_reversed,
          i0_array, block_rank_matrix);
    if (i0) *i0 = i0_result;

    // Permute SA to plain array.
    fprintf(stderr, "\nPermute the resulting SA to plain array: ");
    start = utils::wclock();
    result->permute_to_plain_array(max_threads);
    fprintf(stderr, "%.2Lfs\n", utils::wclock() - start);

    delete result;
  } else if (compute_bwt) {
    if (i0) *i0 = i0_array[0];
  }
  delete[] i0_array;
  for (std::uint64_t j = 0; j < n_blocks; ++j)
    delete[] block_rank_matrix[j];
  delete[] block_rank_matrix;

  if (!compute_gt_begin && (n_blocks > 1 || has_tail)) {
    delete gt_begin;
    gt_begin = NULL;
  }

  std::uint8_t *bwt = NULL;
  if (compute_bwt) {

    // Allocate aux, copy bwt into aux.
    fprintf(stderr, "Copy bwtsa.bwt into aux memory: ");
    start = utils::wclock();
    bwt = (std::uint8_t *)malloc(text_length);

#ifdef _OPENMP
    #pragma omp parallel for
    for (std::uint64_t j = 0; j < text_length; ++j) {
      bwt[j] = bwtsa[j].m_bwt;
    }
#else
    for (std::uint64_t j = 0; j < text_length; ++j)
      bwt[j] = bwtsa[j].m_bwt;
#endif

    fprintf(stderr, "%.2Lfs\n", utils::wclock() - start);
  }

  fprintf(stderr, "Shrink bwtsa.sa into sa: ");
  start = utils::wclock();

  parallel_shrink<bwtsa_t<text_offset_type>, text_offset_type>(bwtsa, text_length, max_threads);

  fprintf(stderr, "%.2Lfs\n", utils::wclock() - start);

  if (compute_bwt) {

    // Copy from aux into the end of bwtsa.
    fprintf(stderr, "Copy bwt from aux memory to the end of bwtsa: ");
    start = utils::wclock();
    std::uint8_t *dest = (std::uint8_t *)(((text_offset_type *)bwtsa) + text_length);
#ifdef _OPENMP
    #pragma omp parallel for
    for (std::uint64_t j = 0; j < text_length; ++j) {
      dest[j] = bwt[j];
    }
#else
    for (std::uint64_t j = 0; j < text_length; ++j)
      dest[j] = bwt[j];
#endif

    free(bwt);
    fprintf(stderr, "%.2Lfs\n", utils::wclock() - start);
  }

  long double total_sascan_time = utils::wclock() - absolute_start;
  fprintf(stderr, "\n\nComputation finished. Summary:\n");
  fprintf(stderr, "  elapsed time: %.2Lfs (%.4Lfs/MiB)\n", total_sascan_time,
      total_sascan_time / ((long double)text_length / (1 << 20)));
  fprintf(stderr, "  speed: %.2LfMiB/s\n", ((long double)text_length / (1 << 20)) / total_sascan_time);
}

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_PSASCAN_HPP_INCLUDED
