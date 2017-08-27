/**
 * @file    src/psascan_src/space_efficient_isa.hpp
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

#ifndef __SRC_PSASCAN_SRC_SPACE_EFFICIENT_ISA_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_SPACE_EFFICIENT_ISA_HPP_INCLUDED

#include <cstdint>
#include <algorithm>
#include <omp.h>

#include "utils.hpp"


namespace psascan_private {

/**
 * Space efficient ISA representation based on the
 * LZ-ISA algorithm computing Lempel-Ziv (LZ77)
 * factorization described in
 *
 *   Dominik Kempa, Simon J. Puglisi:
 *   Lempel-Ziv Factorization: Simple, Fast, Practical.
 *   In Proc. ALENEX 2013, p. 103-112.
 **/

template<
  typename approx_rank_type,
  typename psa_offset_type,
  std::uint64_t k_sampling_rate_log>
class space_efficient_isa {

  private:
    std::uint64_t m_text_length;
    std::uint64_t m_last_isa;
    std::uint64_t m_i0;
    std::uint8_t *m_mem;
    std::uint64_t *m_count;
    std::uint64_t *m_sparse_isa;

    const psa_offset_type *m_psa;
    const std::uint8_t *m_text;
    const approx_rank_type *m_rank;

    static const std::uint64_t k_sampling_rate;
    static const std::uint64_t k_sampling_rate_mask;
    static const std::uint64_t k_sigma;

  public:

    //=========================================================================
    // Constructor.
    //=========================================================================
    space_efficient_isa(
        const psa_offset_type *psa,
        const std::uint8_t *text,
        const approx_rank_type *rank,
        std::uint64_t text_length,
        std::uint64_t i0) {

      // Initialize members.
      m_psa = psa;
      m_text_length = text_length;
      m_rank = rank;
      m_text = text;
      m_i0 = i0;

      // Allocate all memory.
      std::uint64_t sparse_isa_size =
        (m_text_length + k_sampling_rate - 1) >> k_sampling_rate_log;
      std::uint64_t toalloc =
        k_sigma * sizeof(std::uint64_t) +
        sparse_isa_size * sizeof(std::uint64_t);
      m_mem = utils::allocate_array<std::uint8_t>(toalloc);

      // Assign memory for symbol counts.
      std::uint8_t *mem_ptr = m_mem;
      m_count = (std::uint64_t *)mem_ptr;
      mem_ptr += k_sigma * sizeof(std::uint64_t);

      // Assign memory for sparse ISA.
      m_sparse_isa = (std::uint64_t *)mem_ptr;
      mem_ptr += sparse_isa_size * sizeof(std::uint64_t);

      // Compute sparse ISA.
#ifdef _OPENMP

      // Parallel version.
      #pragma omp parallel for
      for (std::uint64_t i = 0; i < m_text_length; ++i) {
        std::uint64_t sa_value = psa[i];
        if (!(sa_value & k_sampling_rate_mask))
          m_sparse_isa[sa_value >> k_sampling_rate_log] = i;
        if (sa_value + 1 == m_text_length)
          m_last_isa = i;
      }
#else

      // Sequential version.
      for (std::uint64_t i = 0; i < m_text_length; ++i) {
        std::uint64_t sa_value = psa[i];
        if (!(sa_value & k_sampling_rate_mask))
          m_sparse_isa[sa_value >> k_sampling_rate_log] = i;
        if (sa_value + 1 == m_text_length)
          m_last_isa = i;
      }
#endif

      // Zero-initialize symbol counts.
      std::fill(m_count, m_count + k_sigma, 0);

      // Compute symbol counts.
      for (std::uint64_t c = 0; c < k_sigma; ++c)
        m_count[c] = m_rank->count((std::uint8_t)c);

      // Modify symbol counts.
      ++m_count[text[m_text_length - 1]];
      --m_count[0];

      // Exclusive partial sum over symbol counts.
      for (std::uint64_t i = 0, sum = 0; i < k_sigma; ++i) {
        std::uint64_t temp = m_count[i];
        m_count[i] = sum;
        sum += temp;
      }
    }

    //=========================================================================
    // Return the position p such that m_psa[p] = j.
    //=========================================================================
    inline std::uint64_t query(std::uint64_t j) const {

      // Obtain the initial approximation of
      // the final position from ISA samples.
      std::uint64_t isa_i = 0;
      std::uint64_t i =
        ((j + k_sampling_rate - 1) >> k_sampling_rate_log);
      if ((i << k_sampling_rate_log) < m_text_length) {
        isa_i = m_sparse_isa[i];
        i <<= k_sampling_rate_log;
      } else {
        isa_i = m_last_isa;
        i = m_text_length - 1;
      }

      // Apply backward search
      // until the answer is found.
      while (i != j) {

        // Compute ISA[i - 1] from ISA[i].
        // Invariants: i > j, isa_i = ISA[i].
        std::uint8_t c = m_text[i - 1];
        std::uint64_t delta =
          (isa_i > m_i0 && c == 0);
        isa_i = m_count[c] + m_rank->rank(isa_i, c);
        if (isa_i > 0)
          isa_i -= delta;
        while ((std::uint64_t)m_psa[isa_i] + 1 != i)
          ++isa_i;
        --i;
      }

      // Return the answer.
      return isa_i;
    }

    //=========================================================================
    // Destructor.
    //=========================================================================
    ~space_efficient_isa() {

      // Clean up.
      utils::deallocate(m_mem);
    }
};

template<
  typename approx_rank_type,
  typename psa_offset_type,
  std::uint64_t k_sampling_rate_log>
const std::uint64_t space_efficient_isa<approx_rank_type,
      psa_offset_type, k_sampling_rate_log>
  ::k_sampling_rate = ((std::uint64_t)1 << k_sampling_rate_log);

template<
  typename approx_rank_type,
  typename psa_offset_type,
  std::uint64_t k_sampling_rate_log>
const std::uint64_t space_efficient_isa<approx_rank_type,
      psa_offset_type, k_sampling_rate_log>
  ::k_sampling_rate_mask = ((std::uint64_t)1 << k_sampling_rate_log) - 1;

template<
  typename approx_rank_type,
  typename psa_offset_type,
  std::uint64_t k_sampling_rate_log>
const std::uint64_t space_efficient_isa<approx_rank_type,
      psa_offset_type, k_sampling_rate_log>
  ::k_sigma = (std::uint64_t)256;

}  // namespace psascan_private

#endif  // __SRC_PSASCAN_SRC_SPACE_EFFICIENT_ISA_HPP_INCLUDED
