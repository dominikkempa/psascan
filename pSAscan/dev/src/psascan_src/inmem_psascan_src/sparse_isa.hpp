/**
 * @file    psascan_src/inmem_psascan_src/sparse_isa.hpp
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

#ifndef __PSASCAN_SRC_INMEM_PSASCAN_SRC_SPARSE_ISA_HPP_INCLUDED
#define __PSASCAN_SRC_INMEM_PSASCAN_SRC_SPARSE_ISA_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <omp.h>


namespace psascan_private {
namespace inmem_psascan_private {

/**
 * Sparse ISA encoding based on the ISAs algorithm computing
 * Lempel-Ziv (LZ77) factorization described in
 *
 *   Dominik Kempa, Simon J. Puglisi:
 *   Lempel-Ziv factorization: Simple, fast, practical.
 *   In Proc. ALENEX 2013, p. 103-112.
 **/

template<typename pagearray_type, typename rank_type, std::uint64_t k_isa_sampling_rate_log>
class sparse_isa {
  private:
    static const std::uint64_t k_isa_sampling_rate;
    static const std::uint64_t k_isa_sampling_rate_mask;
    static const std::uint64_t k_sigma;

  public:
    sparse_isa(const pagearray_type *bwtsa, const std::uint8_t *text,
        const rank_type *rank, std::uint64_t length, std::uint64_t i0) {
      m_bwtsa = bwtsa;
      m_length = length;
      m_rank = rank;
      m_text = text;
      m_i0 = i0;

      if (m_length <= 0) {
        fprintf(stderr, "\nError: m_length in the constructor of sparse_isa\n");
        std::exit(EXIT_FAILURE);
      }

      std::uint64_t items = (m_length + k_isa_sampling_rate - 1) / k_isa_sampling_rate + 1;
      m_sparse_isa = (std::uint64_t *)malloc(items * sizeof(std::uint64_t));

#ifdef _OPENMP
      #pragma omp parallel for
      for (std::uint64_t j = 0; j < m_length; ++j) {
        std::uint64_t sa_j = (*m_bwtsa)[j].m_sa;
        if (!(sa_j & k_isa_sampling_rate_mask))
          m_sparse_isa[sa_j >> k_isa_sampling_rate_log] = j;
        if (sa_j + 1 == m_length) m_last_isa = j;
      }
#else
      for (std::uint64_t j = 0; j < m_length; ++j) {
        std::uint64_t sa_j = (*m_bwtsa)[j].m_sa;
        if (!(sa_j & k_isa_sampling_rate_mask))
          m_sparse_isa[sa_j >> k_isa_sampling_rate_log] = j;
        if (sa_j + 1 == m_length) m_last_isa = j;
      }
#endif

      m_count = (std::uint64_t *)malloc(k_sigma * sizeof(std::uint64_t));
      for (std::uint64_t j = 0; j < k_sigma; ++j)
        m_count[j] = rank->rank(length, (std::uint8_t)j);

      ++m_count[text[length - 1]];
      --m_count[0];

      for (std::uint64_t i = 0, s = 0; i < k_sigma; ++i) {
        std::uint64_t t = m_count[i];
        m_count[i] = s;
        s += t;
      }
    }

    inline std::uint64_t query(std::uint64_t j) const {
      std::int64_t isa_i;
      std::uint64_t i = ((j + k_isa_sampling_rate - 1) >> k_isa_sampling_rate_log);
      if ((i << k_isa_sampling_rate_log) < m_length) {
        isa_i = m_sparse_isa[i];
        i <<= k_isa_sampling_rate_log;
      } else {
        isa_i = m_last_isa;
        i = m_length - 1;
      }

      while (i != j) {
        // Compute ISA[i - 1] from ISA[i].
        // Invariant:
        //   i > 0
        //   isa_i >= 0
        //   isa_i = ISA[i]
        //   j <= i
        std::uint8_t c = m_text[i - 1];
        std::int64_t delta = ((std::uint64_t)isa_i > m_i0 && c == 0);

        isa_i = (std::int64_t)m_count[c] + (std::int64_t)m_rank->rank(isa_i, c) - delta;
        if (isa_i < 0 || ((std::uint64_t)((*m_bwtsa)[isa_i].m_sa)) + 1 != i)
          ++isa_i;

        --i;
      }

      return (std::uint64_t)isa_i;
    }

    ~sparse_isa() {
      free(m_sparse_isa);
      free(m_count);
    }

  private:
    std::uint64_t m_length;
    std::uint64_t m_last_isa;
    std::uint64_t m_i0;

    std::uint64_t *m_count;
    std::uint64_t *m_sparse_isa;

    const std::uint8_t *m_text;
    const pagearray_type *m_bwtsa;
    const rank_type *m_rank;
};

template<typename pagearray_type, typename rank_type, std::uint64_t k_isa_sampling_rate_log>
  const std::uint64_t sparse_isa<pagearray_type, rank_type, k_isa_sampling_rate_log>
  ::k_isa_sampling_rate = (1UL << k_isa_sampling_rate_log);

template<typename pagearray_type, typename rank_type, std::uint64_t k_isa_sampling_rate_log>
  const std::uint64_t sparse_isa<pagearray_type, rank_type, k_isa_sampling_rate_log>
  ::k_isa_sampling_rate_mask = (1UL << k_isa_sampling_rate_log) - 1;

template<typename pagearray_type, typename rank_type, std::uint64_t k_isa_sampling_rate_log>
  const std::uint64_t sparse_isa<pagearray_type, rank_type, k_isa_sampling_rate_log>
  ::k_sigma = 256UL;


}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __PSASCAN_SRC_INMEM_PSASCAN_SRC_SPARSE_ISA_HPP_INCLUDED
