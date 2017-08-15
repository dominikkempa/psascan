/**
 * @file    src/psascan_src/inmem_psascan_src/inmem_bwtsa_merge.hpp
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

#ifndef __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_BWTSA_MERGE_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_BWTSA_MERGE_HPP_INCLUDED

#include <cstdio>
#include <cstdint>
#include <vector>
#include <string>
#include <algorithm>

#include "../bitvector.hpp"
#include "../io/multifile.hpp"
#include "inmem_gap_array.hpp"
#include "inmem_compute_gap.hpp"
#include "parallel_merge.hpp"
#include "pagearray.hpp"
#include "bwtsa.hpp"
#include "merge_schedule.hpp"


namespace psascan_private {
namespace inmem_psascan_private {

template<typename block_offset_type, unsigned pagesize_log>
pagearray<bwtsa_t<block_offset_type>, pagesize_log> *inmem_bwtsa_merge(
    const std::uint8_t *text,
    std::uint64_t text_length,
    bwtsa_t<block_offset_type> *bwtsa,
    bitvector *gt,
    std::uint64_t max_block_size,
    std::uint64_t range_beg,
    std::uint64_t range_end,
    std::uint64_t max_threads,
    bool need_gt,
    bool need_bwt,
    std::uint64_t &result_i0,
    MergeSchedule &schedule,
    std::uint64_t text_beg,
    std::uint64_t text_end,
    std::uint64_t supertext_length,
    std::string supertext_filename,
    const multifile *tail_gt_begin_reversed,
    std::int64_t *i0_array,
    std::uint64_t **block_rank_matrix) {
  typedef pagearray<bwtsa_t<block_offset_type>, pagesize_log> pagearray_type;

  std::uint64_t shift = (max_block_size - text_length % max_block_size) % max_block_size;
  std::uint64_t range_size = range_end - range_beg;

  if (range_size == 1) {
    std::uint64_t block_beg = range_beg * max_block_size;
    std::uint64_t block_end = block_beg + max_block_size;
    block_beg = (std::uint64_t)std::max(0L, (std::int64_t)block_beg - (std::int64_t)shift);
    block_end -= shift;

    result_i0 = i0_array[range_beg];
    pagearray_type *bwtsa_pagearray =
      new pagearray_type(bwtsa + block_beg, bwtsa + block_end);
    return bwtsa_pagearray;
  }

  //----------------------------------------------------------------------------
  // STEP 1: Split the blocks in the left and right group.
  //----------------------------------------------------------------------------
  std::uint64_t lrange_size = schedule.left_size(range_size);
  std::uint64_t rrange_size = range_size - lrange_size;

  std::uint64_t lrange_beg = range_beg;
  std::uint64_t lrange_end = range_beg + lrange_size;
  std::uint64_t rrange_beg = lrange_end;
  std::uint64_t rrange_end = rrange_beg + rrange_size;

  std::uint64_t lbeg = lrange_beg * max_block_size;
  std::uint64_t rbeg = rrange_beg * max_block_size;
  std::uint64_t lend = rbeg;
  std::uint64_t rend = rbeg + rrange_size * max_block_size;
  lbeg = (std::uint64_t)std::max(0L, (std::int64_t)lbeg - (std::int64_t)shift);
  rbeg -= shift;
  lend -= shift;
  rend -= shift;

  std::uint64_t lsize = lend - lbeg;
  std::uint64_t rsize = rend - rbeg;

  //----------------------------------------------------------------------------
  // STEP 2: Compute partial SAs and BWTs for left and right block.
  //----------------------------------------------------------------------------

  // 2.a
  //
  // Left block
  std::uint64_t left_i0;
  pagearray_type *l_bwtsa = inmem_bwtsa_merge<block_offset_type, pagesize_log>(text,
      text_length, bwtsa, gt, max_block_size, lrange_beg, lrange_end,
      max_threads, need_gt, true, left_i0, schedule, text_beg, text_end,
      supertext_length, supertext_filename, tail_gt_begin_reversed, i0_array,
      block_rank_matrix);

  // 2.b
  // 
  // Right block
  std::uint64_t right_i0;
  pagearray_type *r_bwtsa = inmem_bwtsa_merge<block_offset_type, pagesize_log>(text,
      text_length, bwtsa, gt, max_block_size, rrange_beg, rrange_end,
      max_threads, true, need_bwt, right_i0, schedule, text_beg, text_end,
      supertext_length, supertext_filename, tail_gt_begin_reversed, i0_array,
      block_rank_matrix);

  //----------------------------------------------------------------------------
  // STEP 3: Merge partial SAs and BWTs.
  //----------------------------------------------------------------------------
  fprintf(stderr, "Merge blocks %ld-%ld with %ld-%ld\n",
      lrange_beg + 1, lrange_end, rrange_beg + 1, rrange_end);
  long double start = utils::wclock();

  // 3.a
  //
  // Compute gap
  fprintf(stderr, "  Compute gap:\n");
  inmem_gap_array *gap;
  long double rank_init_time;
  long double streaming_time;
  long double start1 = utils::wclock();
  inmem_compute_gap<block_offset_type, pagesize_log>(text, text_length, lbeg, lsize,
      rsize, *l_bwtsa, gt, gap, max_threads, need_gt, left_i0, (1L << 21),
      rank_init_time, streaming_time, block_rank_matrix, lrange_beg,
      lrange_size, rrange_size);
  fprintf(stderr, "  Total time: %.2Lfs\n", utils::wclock() - start1);

  // 3.b
  //
  // Merge partial SAs and BWTs
  fprintf(stderr, "  Merge SA/BWT:  ");
  start1 = utils::wclock();
  std::uint64_t delta_i0;
  if (need_bwt)
    (*r_bwtsa)[right_i0].m_bwt = text[rbeg - 1];
  pagearray_type *result = parallel_merge(l_bwtsa, r_bwtsa, gap,
      max_threads, left_i0, delta_i0, lsize);
  result_i0 = left_i0 + delta_i0;
  long double merging_time = utils::wclock() - start1;
  fprintf(stderr, "total: %.2Lfs\n", merging_time);

  // 3.c
  //
  // Clean up.
  start1 = utils::wclock();
  delete l_bwtsa;
  delete r_bwtsa;
  delete gap;
  long double cleaning_time = utils::wclock() - start1;
  if (cleaning_time > 0.2L)
    fprintf(stderr, "Clean: %.2Lfs\n", cleaning_time);

  long double time_per_elem_left = merging_time / (lsize + rsize) + rank_init_time / lsize;
  long double time_per_elem_right = merging_time / (lsize + rsize) + streaming_time / rsize;
  long double ratio = time_per_elem_right / time_per_elem_left;
  fprintf(stderr, "Total time: %.2Lfs (rl_ratio = %.3Lf)\n",
      utils::wclock() - start, ratio);

  return result;
}

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_BWTSA_MERGE_HPP_INCLUDED
