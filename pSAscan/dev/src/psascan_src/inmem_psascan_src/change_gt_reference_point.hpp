/**
 * @file    psascan_src/inmem_psascan_src/change_gt_reference_point.hpp
 * @section LICENCE
 *
 * This file is part of pSAscan v0.2.0
 * See: http://www.cs.helsinki.fi/group/pads/
 *
 * Copyright (C) 2014-2016
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

#ifndef __PSASCAN_SRC_INMEM_PSASCAN_SRC_CHANGE_GT_REFERENCE_POINT_HPP_INCLUDED
#define __PSASCAN_SRC_INMEM_PSASCAN_SRC_CHANGE_GT_REFERENCE_POINT_HPP_INCLUDED

#include <cstdint>
#include <cstring>
#include <algorithm>
#include <thread>

#include "../bitvector.hpp"
#include "srank_aux.hpp"


namespace psascan_private {
namespace inmem_psascan_private {

/**
 * In-place computation of gt_begin bitvector from gt_end bitvector
 * (reversed). The procedure uses the string range matching algorithm
 * described in
 *
 *   Juha Karkkainen, Dominik Kempa, Simon J. Puglisi:
 *   String Range Matching.
 *   In Proc. CPM 2014, p. 232-241.
 **/

//==============================================================================
// Compute range [microblock_beg..microblock_end) of bits in the output
// bitvector gt_out.
//==============================================================================
void gt_end_to_gt_begin_aux(const std::uint8_t *text, std::uint64_t text_length,
    std::uint64_t block_beg, std::uint64_t block_end, bitvector *gt) {
  std::uint64_t block_size = block_end - block_beg;
  const std::uint8_t *pat = text + block_beg, *txt = pat;

  std::uint64_t i = 1, el = 0, s = 0, p = 0;
  std::uint64_t i_max = i, el_max = 0, s_max = 0, p_max = 0;
  std::uint64_t rev_end = text_length - block_beg;

  while (i < block_size) {
    // Compute lcp(text[left_block_beg..), text[left_block_beg+i..),
    // but compare not more than left_block_size symbols (we have gt
    // to resolve the long comparisons).
    while (block_beg + i + el < block_end && txt[i + el] == pat[el])
      update_ms(pat, ++el, s, p);

    if (((block_beg + i + el != block_end && txt[i + el] > pat[el]) ||
         (block_beg + i + el == block_end && !gt->get(rev_end - i))))
      gt->set(rev_end - i);
    else gt->reset(rev_end - i);

    std::uint64_t j = i_max;
    if (el > el_max) {
      std::swap(el, el_max);
      std::swap(s, s_max);
      std::swap(p, p_max);
      i_max = i;
    }

    if (el < 100) {
      ++i;
      el = 0;
    } else if (p > 0 && (p << 2) <= el && !memcmp(pat, pat + p, s)) {
      std::uint64_t maxk = std::min(block_size - i, p);
      for (std::uint64_t k = 1; k < maxk; ++k) {
        if (gt->get(rev_end - (j + k))) gt->set(rev_end - (i + k));
        else gt->reset(rev_end - (i + k));
      }

      i += p;
      el -= p;
    } else {
      std::uint64_t h = (el >> 2) + 1;
      std::uint64_t maxk = std::min(h, block_size - i);
      for (std::uint64_t k = 1; k < maxk; ++k) {
        if (gt->get(rev_end - (j + k))) gt->set(rev_end - (i + k));
        else gt->reset(rev_end - (i + k));
      }

      i += h;
      el = 0;
      s = 0;
      p = 0;
    }
  }
}


//==============================================================================
// Change gt_end bitvector into gt_begin using string range matching.
//==============================================================================
void gt_end_to_gt_begin(const std::uint8_t *text, std::uint64_t text_length,
    bitvector *gt, std::uint64_t max_block_size) {
  std::uint64_t n_blocks = (text_length + max_block_size - 1) / max_block_size;

  //----------------------------------------------------------------------------
  // STEP 1: Compute the last bit in every block.
  //----------------------------------------------------------------------------
  for (std::uint64_t i = 0; i < n_blocks; ++i) {
    std::uint64_t block_end = text_length - (n_blocks - 1 - i) * max_block_size;
    std::uint64_t rev_beg = text_length - block_end;
    gt->flip(rev_beg);
  }

  //----------------------------------------------------------------------------
  // STEP 2: compute remaining bits in every block.
  //----------------------------------------------------------------------------
  std::thread **threads = new std::thread*[n_blocks];
  for (std::uint64_t i = 0; i < n_blocks; ++i) {
    std::uint64_t block_end = text_length - (n_blocks - 1 - i) * max_block_size;
    std::uint64_t block_beg = (std::uint64_t)std::max(0L, (std::int64_t)block_end - (std::int64_t)max_block_size);

    threads[i] = new std::thread(gt_end_to_gt_begin_aux,
        text, text_length, block_beg, block_end, gt);
  }

  for (std::uint64_t i = 0; i < n_blocks; ++i) threads[i]->join();
  for (std::uint64_t i = 0; i < n_blocks; ++i) delete threads[i];
  delete[] threads;
}

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __PSASCAN_SRC_INMEM_PSASCAN_SRC_CHANGE_GT_REFERENCE_POINT_HPP_INCLUDED
