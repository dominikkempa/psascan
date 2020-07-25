/**
 * @file    src/psascan_src/inmem_psascan_src/inmem_bwt_from_sa.hpp
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

#ifndef __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_BWT_FROM_SA_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_BWT_FROM_SA_HPP_INCLUDED

#include <cstdint>
#include <algorithm>
#include <omp.h>

#include "../utils.hpp"
#include "bwtsa.hpp"


namespace psascan_private {
namespace inmem_psascan_private {


template<typename block_offset_type>
void compute_bwt_in_bwtsa(
    const std::uint8_t *text,
    std::uint64_t text_length,
    bwtsa_t<block_offset_type> *dest,
    std::uint64_t &longest_suffix_rank) {

  longest_suffix_rank = 0;
  std::uint64_t max_threads = omp_get_max_threads();
  std::uint64_t max_block_size =
    (text_length + max_threads - 1) / max_threads;
  std::uint64_t n_blocks =
    (text_length + max_block_size - 1) / max_block_size;

  #pragma omp parallel num_threads(n_blocks)
  {
    std::uint64_t thread_id = omp_get_thread_num();
    std::uint64_t block_beg = thread_id * max_block_size;
    std::uint64_t block_end = std::min(text_length,
        block_beg + max_block_size);

    for (std::uint64_t i = block_beg; i < block_end; ++i) {
      std::uint64_t sa_value = dest[i].m_sa;

      if (sa_value > 0) {
        dest[i].m_bwt = text[sa_value - 1];
      } else {
        dest[i].m_bwt = 0;
        longest_suffix_rank = i;
      }
    }
  }
}

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_BWT_FROM_SA_HPP_INCLUDED
