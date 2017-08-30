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

class bitvector {
  private:
    std::uint64_t m_alloc_bytes;
    std::uint8_t *m_data;

  public:

    //=========================================================================
    // Constructor.
    //=========================================================================
    bitvector(std::string filename) {
      m_alloc_bytes = utils::file_size(filename);
      m_data = utils::allocate_array<std::uint8_t>(m_alloc_bytes);
      utils::read_from_file(m_data, m_alloc_bytes, filename);
    }

    //=========================================================================
    // Constructor.
    //=========================================================================
    bitvector(const std::uint64_t length) {
      m_alloc_bytes = (length + 7) / 8;
      m_data = utils::allocate_array<std::uint8_t>(m_alloc_bytes);
      std::fill(m_data, m_data + m_alloc_bytes, 0);
    }

    //=========================================================================
    // Return the i-th bit,
    // 0 <= i < length.
    //=========================================================================
    inline bool get(const std::uint64_t i) const {
      return m_data[i >> 3] & (1 << (i & 7));
    }

    //=========================================================================
    // Set the i-th bit to 1,
    // 0 <= i < length.
    //=========================================================================
    inline void set(const std::uint64_t i) {
      m_data[i >> 3] |= (1 << (i & 7));
    }

    //=========================================================================
    // Set the i-th bit to 0,
    // 0 <= i < length.
    //=========================================================================
    inline void reset(const std::uint64_t i) {
      m_data[i >> 3] &= (~(1 << (i & 7)));
    }

    //=========================================================================
    // Flip the i-th bit,
    // 0 <= i < length.
    //=========================================================================
    inline void flip(const std::uint64_t i) {
      if (get(i))
        reset(i);
      else set(i);
    }

    //=========================================================================
    // Store the bitvector to a given file.
    //=========================================================================
    inline void save(const std::string filename) const {
      utils::write_to_file(m_data, m_alloc_bytes, filename);
    }

    //=========================================================================
    // Return the number of 1-bits
    // in the range [beg..end).
    //=========================================================================
    inline std::uint64_t range_sum(
        const std::uint64_t beg,
        const std::uint64_t end) const {

      // Initialize the result.
      std::uint64_t result = 0;

      // Accumulate the bits from the first word.
      std::uint64_t j = beg;
      while (j < end && (j & 63))
        result += get(j++);

      // Quickly counts the bits inside
      // words that are entirely in range.
      const std::uint64_t *ptr64 =
        (const std::uint64_t *)(m_data + (j >> 3));
      while (j + 64 <= end) {
        result += __builtin_popcountll(*ptr64++);
        j += 64;
      }

      // Add the bits from the last word.
      while (j < end)
        result += get(j++);

      // Return the result.
      return result;
    }

    //=========================================================================
    // Return the (i+1)-th 0-bit (i = 0, ..) in the range
    // [beg..length), where length is the bitvector size (in bits).
    //=========================================================================
    inline std::uint64_t select0(
        const std::uint64_t beg,
        const std::uint64_t i) const {

      // Slow select in the
      // first word.
      std::uint64_t j = beg;
      std::uint64_t left = i + (std::uint64_t)get(j);
      while (left > 0 && ((j + 1) & 63))
        left -= (1 - get(++j));

      // Fast forward through
      // the rest of the words.
      const std::uint64_t *ptr64 =
        (const std::uint64_t *)(m_data + ((j + 1) >> 3));
      while (left >= 64) {
        j += 64;
        left -= (64 - __builtin_popcountll(*ptr64++));
      }

      // Slow select in the
      // last word.
      while (left > 0)
        left -= (1 - get(++j));

      // Return the result.
      return j;
    }

    //=========================================================================
    // Return the (i+1)-th 1-bit (i = 0, ..) in the range
    // [beg..length), where length is the bitvector size (in bits).
    //=========================================================================
    inline std::uint64_t select1(
        const std::uint64_t beg,
        const std::uint64_t i) const {

      // Slow select in the
      // first word.
      std::uint64_t j = beg;
      std::uint64_t left = (i + 1) - (std::uint64_t)get(j);
      while (left > 0 && ((j + 1) & 63))
        left -= get(++j);

      // Fast forward through
      // the rest of the words.
      const std::uint64_t *ptr64 =
        (const std::uint64_t *)(m_data + ((j + 1) >> 3));
      while (left >= 64) {
        j += 64;
        left -= __builtin_popcountll(*ptr64++);
      }

      // Slow select in the
      // last word.
      while (left > 0)
        left -= get(++j);

      // Return the result.
      return j;
    }

    //=========================================================================
    // Destuctor.
    //=========================================================================
    ~bitvector() {
      utils::deallocate(m_data);
    }
};

}  // namespace psascan_private

#endif  // __SRC_PSASCAN_SRC_BITVECTOR_HPP_INCLUDED
