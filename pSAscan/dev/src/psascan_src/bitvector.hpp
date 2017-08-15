/**
 * @file    src/psascan_src/bitvector.hpp
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

#ifndef __SRC_PSASCAN_SRC_BITVECTOR_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_BITVECTOR_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>

#include "utils.hpp"


namespace psascan_private {

struct bitvector {
  private:
    std::uint64_t m_alloc_bytes;
    std::uint8_t *m_data;

  public:
    bitvector(std::string filename) {
      m_alloc_bytes = utils::file_size(filename);
      m_data = (std::uint8_t *)malloc(m_alloc_bytes);
      utils::read_from_file(m_data, m_alloc_bytes, filename);
    }

    bitvector(std::uint64_t length) {
      m_alloc_bytes = (length + 7) / 8;
      m_data = (std::uint8_t *)calloc(m_alloc_bytes, 1);
    }

    inline bool get(std::uint64_t i) const {
      return m_data[i >> 3] & (1 << (i & 7));
    }

    inline void set(std::uint64_t i) {
      m_data[i >> 3] |= (1 << (i & 7));
    }

    inline void reset(std::uint64_t i) {
      m_data[i >> 3] &= (~(1 << (i & 7)));
    }

    inline void flip(std::uint64_t i) {
      if (get(i)) reset(i);
      else set(i);
    }

    inline void save(std::string filename) const {
      utils::write_to_file(m_data, m_alloc_bytes, filename);
    }

    // Number of 1 bits in the range [beg..end).
    std::uint64_t range_sum(std::uint64_t beg, std::uint64_t end) const {
      std::uint64_t result = 0;

      std::uint64_t j = beg;
      while (j < end && (j & 63))
        result += get(j++);

      uint64_t *ptr64 = (uint64_t *)(m_data + (j >> 3));
      while (j + 64 <= end) {
        result += __builtin_popcountll(*ptr64++);
        j += 64;
      }

      while (j < end)
        result += get(j++);

      return result;
    }

    ~bitvector() {
      free(m_data);
    }
};

}  // namespace psascan_private

#endif  // __SRC_PSASCAN_SRC_BITVECTOR_HPP_INCLUDED
