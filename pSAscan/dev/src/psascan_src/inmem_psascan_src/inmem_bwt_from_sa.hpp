/**
 * @file    psascan_src/inmem_psascan_src/inmem_bwt_from_sa.hpp
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

#ifndef __PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_BWT_FROM_SA_HPP_INCLUDED
#define __PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_BWT_FROM_SA_HPP_INCLUDED

#include <cstdint>
#include <algorithm>
#include <thread>

#include "../utils.hpp"
#include "bwtsa.hpp"


namespace psascan_private {
namespace inmem_psascan_private {

template<typename saidx_t>
void compute_bwt_in_bwtsa_aux(const std::uint8_t *text, std::uint64_t beg,
    std::uint64_t end, bwtsa_t<saidx_t> *dest, std::int64_t *i0) {
  *i0 = -1;
  for (std::uint64_t j = beg; j < end; ++j) {
    if (dest[j].m_sa > saidx_t(0)) dest[j].m_bwt = text[dest[j].m_sa - 1];
    else { dest[j].m_bwt = 0; *i0 = j; }
  }
}

template<typename saidx_t>
void compute_bwt_in_bwtsa(const std::uint8_t *text, std::uint64_t length,
  bwtsa_t<saidx_t> *dest, std::uint64_t max_threads, std::int64_t &result) {
  std::uint64_t max_block_size = (length + max_threads - 1) / max_threads;
  std::uint64_t n_blocks = (length + max_block_size - 1) / max_block_size;
  std::int64_t *index_0 = new std::int64_t[n_blocks];

  // Compute bwt and find i0, where sa[i0] == 0.
  std::thread **threads = new std::thread*[n_blocks];
  for (std::uint64_t i = 0; i < n_blocks; ++i) {
    std::uint64_t block_beg = i * max_block_size;
    std::uint64_t block_end = std::min(block_beg + max_block_size, length);

    threads[i] = new std::thread(compute_bwt_in_bwtsa_aux<saidx_t>,
        text, block_beg, block_end, dest, index_0 + i);
  }

  for (std::uint64_t i = 0; i < n_blocks; ++i) threads[i]->join();
  for (std::uint64_t i = 0; i < n_blocks; ++i) delete threads[i];
  delete[] threads;

  // Find and return i0.
  result = -1;
  for (std::uint64_t i = 0; i < n_blocks; ++i)
    if (index_0[i] != -1) result = index_0[i];
  delete[] index_0;
}

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_BWT_FROM_SA_HPP_INCLUDED
