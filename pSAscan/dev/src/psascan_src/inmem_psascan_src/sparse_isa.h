/**
 * @file    src/psascan_src/inmem_psascan_src/sparse_isa.h
 * @author  Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * @section DESCRIPTION
 *
 * Sparse ISA encoding based on the ISAs algorithm computing
 * Lempel-Ziv (LZ77) factorization described in
 *
 *   Dominik Kempa, Simon J. Puglisi:
 *   Lempel-Ziv factorization: Simple, fast, practical.
 *   In Proc. ALENEX 2013, p. 103-112.
 *
 * @section LICENCE
 *
 * This file is part of pSAscan v0.1.0
 * See: http://www.cs.helsinki.fi/group/pads/
 *
 * Copyright (C) 2014-2015
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

#ifndef __PSASCAN_SRC_INMEM_PSASCAN_SRC_SPARSE_ISA_H_INCLUDED
#define __PSASCAN_SRC_INMEM_PSASCAN_SRC_SPARSE_ISA_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <thread>


namespace psascan_private {
namespace inmem_psascan_private {

template<typename pagearray_type, typename rank_type, std::uint64_t k_isa_sampling_rate_log>
class sparse_isa {
  private:
    static const std::uint64_t k_isa_sampling_rate;
    static const std::uint64_t k_isa_sampling_rate_mask;
    static const std::uint64_t k_sigma;

  private:
    static void compute_sparse_isa_aux(const pagearray_type &bwtsa,
        std::uint64_t block_beg, std::uint64_t block_end,
        std::uint64_t psa_size, std::uint64_t *sparse_isa,
        std::uint64_t &last) {
      for (std::uint64_t j = block_beg; j < block_end; ++j) {
        std::uint64_t sa_j = bwtsa[j].m_sa;
        if (!(sa_j & k_isa_sampling_rate_mask))
          sparse_isa[sa_j >> k_isa_sampling_rate_log] = j;
        if (sa_j + 1 == psa_size) last = j;
      }
    }

  public:
    sparse_isa(const pagearray_type *bwtsa, const unsigned char *text,
        const rank_type *rank, std::uint64_t length, std::uint64_t i0,
        std::uint64_t max_threads) {
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

      std::uint64_t max_block_size = (m_length + max_threads - 1) / max_threads;
      std::uint64_t n_blocks = (m_length + max_block_size - 1) / max_block_size;

      std::thread **threads = new std::thread*[n_blocks];
      for (std::uint64_t t = 0; t < n_blocks; ++t) {
        std::uint64_t block_beg = t * max_block_size;
        std::uint64_t block_end = std::min(block_beg + max_block_size, m_length);

        threads[t] = new std::thread(compute_sparse_isa_aux, std::ref(*m_bwtsa),
            block_beg, block_end, m_length, m_sparse_isa, std::ref(m_last_isa));
      }

      for (std::uint64_t t = 0; t < n_blocks; ++t) threads[t]->join();
      for (std::uint64_t t = 0; t < n_blocks; ++t) delete threads[t];
      delete[] threads;

      m_count = (std::uint64_t *)malloc(k_sigma * sizeof(std::uint64_t));
      std::copy(rank->m_count, rank->m_count + k_sigma, m_count);
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
        unsigned char c = m_text[i - 1];
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

    const unsigned char *m_text;
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

#endif  // __PSASCAN_SRC_INMEM_PSASCAN_SRC_SPARSE_ISA_H_INCLUDED
