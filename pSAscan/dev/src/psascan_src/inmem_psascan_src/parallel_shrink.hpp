/**
 * @file    psascan_src/inmem_psascan_src/parallel_shrink.hpp
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

#ifndef __PSASCAN_SRC_INMEM_PSASCAN_SRC_PARALLEL_SHRINK_HPP_INCLUDED
#define __PSASCAN_SRC_INMEM_PSASCAN_SRC_PARALLEL_SHRINK_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <thread>


namespace psascan_private {
namespace inmem_psascan_private {

template<typename src_type, typename dest_type>
void parallel_shrink_aux(src_type *src, dest_type *dest, std::uint64_t length) {
  for (std::uint64_t i = 0; i < length; ++i)
    dest[i] = (dest_type)src[i];
}

// Specialization
template<>
void parallel_shrink_aux(bwtsa_t<uint40> *src, uint40 *dest, std::uint64_t length) {
  for (std::uint64_t i = 0; i < length; ++i)
    dest[i] = src[i].m_sa;
}

// Specialization
template<>
void parallel_shrink_aux(bwtsa_t<std::int64_t> *src, std::int64_t *dest, std::uint64_t length) {
  for (std::uint64_t i = 0; i < length; ++i)
    dest[i] = src[i].m_sa;
}


// Specialization
template<>
void parallel_shrink_aux(bwtsa_t<std::int32_t> *src, std::int32_t *dest, std::uint64_t length) {
  for (std::uint64_t i = 0; i < length; ++i)
    dest[i] = src[i].m_sa;
}

template<typename src_type, typename dest_type>
dest_type *parallel_shrink(src_type *tab, std::uint64_t length, std::uint64_t max_threads) {
  dest_type *result = (dest_type *)tab;
  if (length < (1UL << 20)) {
    parallel_shrink_aux<src_type, dest_type>(tab, result, length);
    return result;
  }

  // Compute the index of the smallest element (of type src_type)
  // that lies past the end of the last element of tab
  // (after converting all elemeents to type dest_type).
  std::uint64_t bytes_after_shrinking = length * sizeof(dest_type);
  std::uint64_t split = (bytes_after_shrinking + sizeof(src_type) - 1) / sizeof(src_type);

  // Recursively shrink the part up to (but excluding) split.
  parallel_shrink<src_type, dest_type>(tab, split, max_threads);

  // Move the elements in the range [split, length) in parallel.
  // This is safe (no element overwriting) because of how we
  // computed the split.
  std::uint64_t elems = length - split;
  std::uint64_t max_block_size = (elems + max_threads - 1) / max_threads;
  std::uint64_t n_blocks = (elems + max_block_size - 1) / max_block_size;

  std::thread **threads = new std::thread*[n_blocks];
  for (std::uint64_t i = 0; i < n_blocks; ++i) {
    std::uint64_t block_beg = split + i * max_block_size;
    std::uint64_t block_end = std::min(block_beg + max_block_size, length);
    std::uint64_t block_size = block_end - block_beg;

    threads[i] = new std::thread(parallel_shrink_aux<src_type, dest_type>,
        tab + block_beg, result + block_beg, block_size);
  }

  for (std::uint64_t i = 0; i < n_blocks; ++i) threads[i]->join();
  for (std::uint64_t i = 0; i < n_blocks; ++i) delete threads[i];
  delete[] threads;

  return result;
}

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __PSASCAN_SRC_INMEM_PSASCAN_SRC_PARALLEL_SHRINK_HPP_INCLUDED
