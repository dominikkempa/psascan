/**
 * @file    src/psascan_src/parallel_utils.hpp
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

#ifndef __SRC_PSASCAN_SRC_PARALLEL_UTILS_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_PARALLEL_UTILS_HPP_INCLUDED

#include <vector>
#include <algorithm>
#include <omp.h>


namespace psascan_private {
namespace parallel_utils {

std::uint64_t convert_array_to_vbyte_slab(const std::uint64_t *tab,
    std::uint64_t length,
    std::uint8_t *dest) {

#ifdef _OPENMP
  std::uint64_t max_threads = omp_get_num_threads();
  std::uint64_t max_block_size = (length + max_threads - 1) / max_threads;
  std::uint64_t n_blocks = (length + max_block_size - 1) / max_block_size;
  std::vector<std::uint64_t> block_slab_length(n_blocks + 1, (std::uint64_t)0);

  // Compute the size of slab for each block.
  #pragma omp parallel num_threads(n_blocks)
  {
    std::uint64_t thread_id = omp_get_thread_num();
    std::uint64_t block_beg = thread_id * max_block_size;
    std::uint64_t block_end = std::min(block_beg + max_block_size, length);
    std::uint64_t vbyte_encoding_len = 0;

    for (std::uint64_t i = block_beg; i < block_end; ++i) {
      std::uint64_t value = tab[i];
      while (value > 127) {
        ++vbyte_encoding_len;
        value >>= 7;
      }
      ++vbyte_encoding_len;
    }
    block_slab_length[thread_id + 1] = vbyte_encoding_len;
  }

  // Partial sum over block slab length.
  std::partial_sum(block_slab_length.begin(),
      block_slab_length.end(), block_slab_length.begin());

  // Encode the slab of each block.
  #pragma omp parallel num_threads(n_blocks)
  {
    std::uint64_t thread_id = omp_get_thread_num();
    std::uint64_t block_beg = thread_id * max_block_size;
    std::uint64_t block_end = std::min(block_beg + max_block_size, length);
    std::uint64_t ptr = block_slab_length[thread_id];

    for (std::uint64_t i = block_beg; i < block_end; ++i) {
      std::uint64_t value = tab[i];
      while (value > 127) {
        dest[ptr++] = ((value & 0x7F) | 0x80);
        value >>= 7;
      }
      dest[ptr++] = value;
    }
  }

  return block_slab_length[n_blocks];

#else
  std::uint64_t slab_length = 0;
  for (std::uint64_t i = 0; i < length; ++i) {
    std::uint64_t value = tab[i];
    while (value > 127) {
      dest[slab_length++] = ((value & 0x7F) | 0x80);
      value >>= 7;
    }
    dest[slab_length++] = value;
  }

  return slab_length;
#endif

}

}  // namespace parallel_utils
}  // namespace psascan_private

#endif  // __SRC_PSASCAN_SRC_PARALLEL_UTILS_HPP_INCLUDED
