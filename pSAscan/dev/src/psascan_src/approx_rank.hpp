/**
 * @file    src/psascan_src/approx_rank.hpp
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

#ifndef __SRC_PSASCAN_SRC_APPROX_RANK_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_APPROX_RANK_HPP_INCLUDED

#include <cstdint>
#include <algorithm>
#include <omp.h>

#include "utils.hpp"


namespace psascan_private {

/**
 * Data structure answering approximate rank queries. based on the
 * LZ-ISA algorithm computing Lempel-Ziv (LZ77) factorization described in
 *
 *   Dominik Kempa, Simon J. Puglisi:
 *   Lempel-Ziv Factorization: Simple, Fast, Practical.
 *   In Proc. ALENEX 2013, p. 103-112.
 **/

template<std::uint64_t k_sampling_rate_log>
class approx_rank {
  private:
    static const std::uint64_t k_sampling_rate;
    static const std::uint64_t k_sampling_rate_mask;

    std::uint8_t *m_mem;
    std::uint64_t *m_lists;
    std::uint64_t *m_cum_symbol_count;

  public:

    //=========================================================================
    // Constructor
    //=========================================================================
    approx_rank(
        const std::uint8_t * const text,
        const std::uint64_t text_length) {

      // Set the alphabet size.
      static const std::uint64_t k_sigma = 256;

      // Allocate all memory.
      const std::uint64_t max_samples =
        (text_length >> k_sampling_rate_log);
      const std::uint64_t toalloc =
        (k_sigma + 1) * sizeof(std::uint64_t) +
        max_samples * sizeof(std::uint64_t);
      m_mem = utils::allocate_array<std::uint8_t>(toalloc);

      // Assign memory for symbol counts.
      const std::uint8_t *mem_ptr = m_mem;
      m_cum_symbol_count = (std::uint64_t *)mem_ptr;
      mem_ptr += (k_sigma + 1) * sizeof(std::uint64_t);

      // Assign memory for lists.
      m_lists = (std::uint64_t *)mem_ptr;
      mem_ptr += max_samples * sizeof(std::uint64_t);

      // Compute lists.
      {
#ifdef _OPENMP

        // Compute basic parameters.
        const std::uint64_t max_threads = omp_get_max_threads();
        const std::uint64_t max_block_size =
          (text_length + max_threads - 1) / max_threads;
        const std::uint64_t n_blocks =
          (text_length + max_block_size - 1) / max_block_size;

        // Allocate a set of counts for each block.
        std::uint64_t * const blocks_counts =
          utils::allocate_array<std::uint64_t>(k_sigma * n_blocks);

        // Compute counts in each block.
        #pragma omp parallel num_threads(n_blocks)
        {

          // Compute basic parameters
          // (block boundaries, etc).
          const std::uint64_t block_id = omp_get_thread_num();
          const std::uint64_t block_beg = block_id * max_block_size;
          const std::uint64_t block_end = std::min(text_length,
              block_beg + max_block_size);
          std::uint64_t * const block_counts =
            blocks_counts + block_id * k_sigma;

          // Zero-initialize block counts.
          std::fill(block_counts, block_counts + k_sigma, 0);

          // Fill in block counts.
          for (std::uint64_t i = block_beg; i < block_end; ++i)
            ++block_counts[text[i]];
        }

        // Compute global cumulative
        // symbol counts.
        {
          std::uint64_t total_symbol_count = 0;
          for (std::uint64_t c = 0; c < k_sigma; ++c) {

            // Prefix sum for symbol c.
            std::uint64_t total_c_count = 0;
            for (std::uint64_t block_id = 0;
                block_id < n_blocks; ++block_id) {
              const std::uint64_t temp =
                blocks_counts[block_id * k_sigma + c];
              blocks_counts[block_id * k_sigma + c] = total_c_count;
              total_c_count += temp;
            }

            // Update m_symbol_count and
            // total_symbol_count.
            m_cum_symbol_count[c] = total_symbol_count;
            total_symbol_count += total_c_count;
          }

          // Set the sentinel count.
          m_cum_symbol_count[k_sigma] =
            total_symbol_count;
        }

        // Fill in lists.
        #pragma omp parallel num_threads(n_blocks)
        {

          // Compute basic parameters.
          const std::uint64_t block_id = omp_get_thread_num();
          const std::uint64_t block_beg = block_id * max_block_size;
          const std::uint64_t block_end = std::min(text_length,
              block_beg + max_block_size);
          std::uint64_t * const block_counts =
            blocks_counts + block_id * k_sigma;

          // Fill in lists.
          for (std::uint64_t i = block_beg; i < block_end; ++i) {
            const std::uint8_t c = text[i];
            ++block_counts[c];

            // Store the sample.
            if (!(block_counts[c] & k_sampling_rate_mask)) {
              const std::uint64_t list_beg =
                (m_cum_symbol_count[c] >> k_sampling_rate_log);
              const std::uint64_t offset =
                (block_counts[c] >> k_sampling_rate_log) - 1;
              m_lists[list_beg + offset] = i;
            }
          }
        }

        // Clean up.
        utils::deallocate(blocks_counts);
#else

        // Zero-initialize symbol counts.
        std::fill(m_cum_symbol_count,
            m_cum_symbol_count + k_sigma, 0);

        // Compute symbol counts.
        for (std::uint64_t i = 0; i < text_length; ++i)
          ++m_cum_symbol_count[text[i]];

        // Compute exclusive prefix
        // sum over symbol counts.
        {
          std::uint64_t total_symbol_count = 0;
          for (std::uint64_t c = 0; c < k_sigma; ++c) {
            const std::uint64_t total_c_count =
              m_cum_symbol_count[c];

            // Compute list size for c.
            m_cum_symbol_count[c] = total_symbol_count;
            total_symbol_count += total_c_count;
          }

          // Set the sentinel count.
          m_cum_symbol_count[k_sigma] =
            total_symbol_count;
        }

        // Fill in lists.
        {

          // Allocate and  zero-initialize
          // temporary symbol counts.
          std::uint64_t * const symbol_count =
            utils::allocate_array<std::uint64_t>(k_sigma);
          std::fill(symbol_count, symbol_count + k_sigma, 0);

          // Compute lists.
          for (std::uint64_t i = 0; i < text_length; ++i) {
            const std::uint8_t c = text[i];
            ++symbol_count[c];

            // Store the sample.
            if (!(symbol_count[c] & k_sampling_rate_mask)) {
              const std::uint64_t list_beg =
                (m_cum_symbol_count[c] >> k_sampling_rate_log);
              const std::uint64_t offset =
                (symbol_count[c] >> k_sampling_rate_log) - 1;
              m_lists[list_beg + offset] = i;
            }
          }

          // Clean up.
          utils::deallocate(symbol_count);
        }

#endif  // _OPENMP
      }
    }

    //=========================================================================
    // Return a lower bound for the number of occurrences of c in
    // text[0..i), i.e., a value x such that x <= rank(i, c). It is
    // guaranteed to be a "good" lower bound, i.e., to satisfy
    // rank(i, c) - x < k_sampling_rate.
    //=========================================================================
    inline std::uint64_t rank(
        const std::uint64_t i,
        const std::uint8_t c) const {

      // Compute the list size and pointer.
      const std::uint64_t list_beg =
        (m_cum_symbol_count[c] >> k_sampling_rate_log);
      const std::uint64_t count_c =
        m_cum_symbol_count[(std::uint64_t)c + 1] - m_cum_symbol_count[c];
      const std::uint64_t list_size = (count_c >> k_sampling_rate_log);
      const std::uint64_t * const m_list = m_lists + list_beg;

      // Handle special case.
      if (list_size == 0 ||
          m_list[0] >= i)
        return 0;

      // Find the number of items smaller than i.
      std::uint64_t left = 0;
      std::uint64_t right = list_size;
      while (left + 1 != right) {

        // Invariant: the answer is in range (left, right].
        const std::uint64_t mid = (left + right) / 2;
        if (m_list[mid] < i) left = mid;
        else right = mid;
      }

      // Return the answer.
      return (right << k_sampling_rate_log);
    }

    //=========================================================================
    // Return the number of occurrences of c in text.
    //=========================================================================
    inline std::uint64_t count(const std::uint8_t c) const {
      return
        m_cum_symbol_count[(std::uint64_t)c + 1] -
        m_cum_symbol_count[c];
    }

    //=========================================================================
    // Destructor
    //=========================================================================
    ~approx_rank() {
      utils::deallocate(m_mem);
    }
};

template<std::uint64_t k_sampling_rate_log>
  const std::uint64_t approx_rank<k_sampling_rate_log>::
  k_sampling_rate = ((std::uint64_t)1 << k_sampling_rate_log);

template<std::uint64_t k_sampling_rate_log>
  const std::uint64_t approx_rank<k_sampling_rate_log>::
  k_sampling_rate_mask = ((std::uint64_t)1 << k_sampling_rate_log) - 1;

}  // namespace psascan_private

#endif // __SRC_PSASCAN_SRC_APPROX_RANK_HPP_INCLUDED
