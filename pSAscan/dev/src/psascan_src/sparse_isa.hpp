/**
 * @file    src/psascan_src/sparse_isa.hpp
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

#ifndef __SRC_PSASCAN_SRC_SPARSE_ISA_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_SPARSE_ISA_HPP_INCLUDED

#include <cstdint>
#include <algorithm>
#include <omp.h>

#include "utils.hpp"


namespace psascan_private {

//-----------------------------------------------------------------------------
// Sparse ISA encoding based on the ISAs algorithm computing
// Lempel-Ziv (LZ77) factorization described in
//
//   Dominik Kempa, Simon J. Puglisi:
//   Lempel-Ziv Factorization: Simple, Fast, Practical.
//   In Proc. ALENEX 2013, p. 103-112.
//-----------------------------------------------------------------------------
template<
  typename approx_rank_type,
  typename sa_int_type,
  std::uint64_t k_sampling_rate_log>
class sparse_isa {
  private:
    std::uint64_t m_length;
    std::uint64_t m_last_isa;
    std::uint64_t m_i0;

    std::uint64_t *m_count;
    std::uint64_t *m_sparse_isa;

    const sa_int_type *m_psa;
    const std::uint8_t *m_text;
    const approx_rank_type *m_rank;
  
    static const std::uint64_t k_sampling_rate;
    static const std::uint64_t k_sampling_rate_mask;
    static const std::uint64_t k_sigma = 256;

  public:
    sparse_isa(
        const sa_int_type *psa,
        const std::uint8_t *text,
        std::uint64_t length,
        std::uint64_t i0,
        const approx_rank_type *rank) {

      // Initialize members.
      m_psa = psa;
      m_length = length;
      m_rank = rank;
      m_text = text;
      m_i0 = i0;

      // Allocate sparse ISA.
      std::uint64_t sparse_isa_size =
        (m_length + k_sampling_rate - 1) / k_sampling_rate + 1;
      m_sparse_isa =
        utils::allocate_array<std::uint64_t>(sparse_isa_size);

#ifdef _OPENMP

      // Compute sparse ISA.
      #pragma omp parallel for
      for (std::uint64_t i = 0; i < m_length; ++i) {
        std::uint64_t psa_value = psa[i];
        if (!(psa_value & k_sampling_rate_mask))
          m_sparse_isa[psa_value >> k_sampling_rate_log] = i;
        if (psa_value + 1 == m_length)
          m_last_isa = i;
      }
#else

      // Compute sparse ISA.
      for (std::uint64_t i = 0; i < m_length; ++i) {
        std::uint64_t psa_value = psa[i];
        if (!(psa_value & k_sampling_rate_mask))
          m_sparse_isa[psa_value >> k_sampling_rate_log] = i;
        if (psa_value + 1 == m_length)
          m_last_isa = i;
      }
#endif

      // Compute symbol counts.
      m_count = utils::allocate_array<std::uint64_t>(k_sigma);
      for (std::uint64_t i = 0; i < k_sigma; ++i)
        m_count[i] = rank->query(m_length, (std::uint8_t)i);
      ++m_count[text[m_length - 1]];
      --m_count[0];

      // Exclusive partial sum over symbol counts.
      for (std::uint64_t i = 0, sum = 0; i < k_sigma; ++i) {
        std::uint64_t temp = m_count[i];
        m_count[i] = sum;
        sum += temp;
      }
    }

    inline std::uint64_t query(std::uint64_t j) const {
      std::uint64_t isa_i;
      std::uint64_t i = ((j + k_sampling_rate - 1) >> k_sampling_rate_log);
      if ((i << k_sampling_rate_log) < m_length) {
        isa_i = m_sparse_isa[i];
        i <<= k_sampling_rate_log;
      } else {
        isa_i = m_last_isa;
        i = m_length - 1;
      }

      while (i != j) {

        // Compute ISA[i - 1] from ISA[i].
        // Invariant: isa_i = ISA[i], j <= i.
        std::uint8_t c = m_text[i - 1];
        std::uint64_t delta = (isa_i > m_i0 && c == 0);
        std::int64_t temp = (m_count[c] + m_rank->query(isa_i, c)) - delta;
        while (temp < 0 || (std::uint64_t)m_psa[temp] + 1 != i)
          ++temp;

        isa_i = temp;
        --i;
      }

      return isa_i;
    }

    ~sparse_isa() {
      utils::deallocate(m_count);
      utils::deallocate(m_sparse_isa);
    }
};

template<
  typename approx_rank_type,
  typename sa_int_type,
  std::uint64_t k_sampling_rate_log>
const std::uint64_t sparse_isa<approx_rank_type,
      sa_int_type, k_sampling_rate_log>
  ::k_sampling_rate = ((std::uint64_t)1 << k_sampling_rate_log);

template<
  typename approx_rank_type,
  typename sa_int_type,
  std::uint64_t k_sampling_rate_log>
const std::uint64_t sparse_isa<approx_rank_type,
      sa_int_type, k_sampling_rate_log>
  ::k_sampling_rate_mask = ((std::uint64_t)1 << k_sampling_rate_log) - 1;

}  // namespace psascan_private

#endif  // __SRC_PSASCAN_SRC_SPARSE_ISA_HPP_INCLUDED
