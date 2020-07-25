/**
 * @file    src/psascan_src/merge_bwt.hpp
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

#ifndef __SRC_PSASCAN_SRC_MERGE_BWT_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_MERGE_BWT_HPP_INCLUDED

#include <cstdint>
#include <algorithm>
#include <omp.h>

#include "bitvector.hpp"
#include "ranksel_support.hpp"


namespace psascan_private {

//==============================================================================
// Let block_size = left_block_size + right_block_size. Given two arrays
// left_bwt[0..left_block_size) and right_bwt[0..right_block_size) and
// a bitvector bv[0..block_size), compute the array bwt[0..block_size)
// defined as, for 0 <= i < block_size, bwt[i] = left_bwt[j] if bv[i] = 0
// and j is the number of 0-bits in bv[0..i), and bwt[i] = right_bwt[k] if
// bv[i] = 1 and k is the number of 1-bits in bv[0..i).
//
// The function additionally is given two positions: left_block_i0 in
// [0..left_block_size) and right_block_i0 in [0..right_block_size).
// After computing bwt, the function changes bwt[j] to left_block_last,
// where j is the position of (right_block_i0 + 1)-th 1-bit in bv and
// returns the position of (left_block_i0 + 1)-th 0-bit in bv.
//==============================================================================
std::uint64_t merge_bwt(
    const std::uint8_t * const left_bwt,
    const std::uint8_t * const right_bwt,
    const bitvector * const bv,
    const std::uint64_t left_block_size,
    const std::uint64_t right_block_size,
    const std::uint64_t left_block_i0,
    const std::uint64_t right_block_i0,
    const std::uint8_t left_block_last,
    std::uint8_t * const bwt) {

  // Initialize basic parameters.
  const std::uint64_t block_size = left_block_size + right_block_size;

  // Initialize rank/select queries support for bv.
  typedef ranksel_support<> ranksel_support_type;
  ranksel_support_type * const bv_ranksel =
    new ranksel_support_type(bv, block_size);

  // Merge left and right bwt.
  {

#ifdef _OPENMP

    // Compute the number of ranges.
    const std::uint64_t max_threads = omp_get_max_threads();
    const std::uint64_t max_range_size =
      (block_size + max_threads - 1) / max_threads;
    const std::uint64_t n_ranges =
      (block_size + max_range_size - 1) / max_range_size;

    // Each thread handles on range.
    #pragma omp parallel num_threads(n_ranges)
    {

      // Compute the range boundaries.
      const std::uint64_t range_id = omp_get_thread_num();
      const std::uint64_t range_beg = range_id * max_range_size;
      const std::uint64_t range_end = std::min(block_size,
          range_beg + max_range_size);

      // Compute the pointers to
      // left and right BWT.
      std::uint64_t left_bwt_ptr =
        bv_ranksel->rank0(range_beg);
      std::uint64_t right_bwt_ptr =
        range_beg - left_bwt_ptr;

      // Actual merging follows.
      for (std::uint64_t i = range_beg; i < range_end; ++i) {
        if (bv->get(i))
          bwt[i] = right_bwt[right_bwt_ptr++];
        else bwt[i] = left_bwt[left_bwt_ptr++];
      }
    }
#else

    // Sequential version.
    std::uint64_t left_bwt_ptr = 0;
    std::uint64_t right_bwt_ptr = 0;
    for (std::uint64_t i = 0; i < block_size; ++i) {
      if (bv->get(i))
        bwt[i] = right_bwt[right_bwt_ptr++];
      else bwt[i] = left_bwt[left_bwt_ptr++];
    }

#endif  // _OPENMP
  }

  // Find position j = select_1(bv, right_block_i0)
  // and replace bwt[j] with left_block_last.
  const std::uint64_t j = bv_ranksel->select1(right_block_i0);
  bwt[j] = left_block_last;

  // Compute the returned value.
  const std::uint64_t block_i0 =
    bv_ranksel->select0(left_block_i0);

  // Clean up.
  delete bv_ranksel;

  // Return the result.
  return block_i0;
}

}  // namespace psascan_private

#endif  // __SRC_PSASCAN_SRC_MERGE_BWT_HPP_INCLUDED
