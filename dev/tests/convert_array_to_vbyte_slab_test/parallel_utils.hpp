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

#include <cstdint>
#include <algorithm>
#include <omp.h>

#include "utils.hpp"


namespace psascan_private {
namespace parallel_utils {

//=============================================================================
// Compute the vbyte encoding of tab[0..length) and store
// into dest array. The dest array has been preallocated
// and is sufficiently large to accomodate the result.
//=============================================================================
std::uint64_t convert_array_to_vbyte_slab(
    const std::uint64_t * const tab,
    const std::uint64_t length,
    std::uint8_t * const dest) {

#ifdef _OPENMP

  // Compute the number of blocks.
  const std::uint64_t max_threads = omp_get_max_threads();
  const std::uint64_t max_block_size =
    (length + max_threads - 1) / max_threads;
  const std::uint64_t n_blocks =
    (length + max_block_size - 1) / max_block_size;

  // Allocate the array holding
  // the size of slab for each block.
  std::uint64_t * const block_slab_length =
    utils::allocate_array<std::uint64_t>(n_blocks);
  std::fill(block_slab_length, block_slab_length + n_blocks, 0);

  // Compute the size of slab for each block.
  #pragma omp parallel num_threads(n_blocks)
  {

    // Compute block boundaries.
    const std::uint64_t block_id = omp_get_thread_num();
    const std::uint64_t block_beg = block_id * max_block_size;
    const std::uint64_t block_end =
      std::min(block_beg + max_block_size, length);

    // Compute the size of slab.
    std::uint64_t vbyte_encoding_len = 0;
    for (std::uint64_t i = block_beg; i < block_end; ++i) {
      std::uint64_t value = tab[i];
      while (value > 127) {
        ++vbyte_encoding_len;
        value >>= 7;
      }
      ++vbyte_encoding_len;
    }

    // Store the result.
    block_slab_length[block_id] =
      vbyte_encoding_len;
  }

  // Exclusive partial sum
  // over block slab lengths.
  std::uint64_t total_slab_length = 0;
  for (std::uint64_t block_id = 0;
      block_id < n_blocks; ++block_id) {
    const std::uint64_t temp = block_slab_length[block_id];
    block_slab_length[block_id] = total_slab_length;
    total_slab_length += temp;
  }

  // Encode the slab of each block.
  #pragma omp parallel num_threads(n_blocks)
  {

    // Compute block boundaries.
    const std::uint64_t block_id = omp_get_thread_num();
    const std::uint64_t block_beg = block_id * max_block_size;
    const std::uint64_t block_end =
      std::min(block_beg + max_block_size, length);

    // Compute the slab.
    std::uint64_t pos = block_slab_length[block_id];
    for (std::uint64_t i = block_beg; i < block_end; ++i) {
      std::uint64_t value = tab[i];
      while (value > 127) {
        dest[pos++] = ((value & 0x7F) | 0x80);
        value >>= 7;
      }
      dest[pos++] = value;
    }
  }

  // Clean up.
  utils::deallocate(block_slab_length);

  // Return the result.
  return total_slab_length;

#else

  // Sequential version.
  std::uint64_t slab_length = 0;
  for (std::uint64_t i = 0; i < length; ++i) {
    std::uint64_t value = tab[i];
    while (value > 127) {
      dest[slab_length++] = ((value & 0x7F) | 0x80);
      value >>= 7;
    }
    dest[slab_length++] = value;
  }

  // Return the result.
  return slab_length;
#endif

}

}  // namespace parallel_utils
}  // namespace psascan_private

#endif  // __SRC_PSASCAN_SRC_PARALLEL_UTILS_HPP_INCLUDED
