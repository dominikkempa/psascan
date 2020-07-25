/**
 * @file    src/psascan_src/ranksel_support.hpp
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

#ifndef __SRC_PSASCAN_SRC_RANKSEL_SUPPORT_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_RANKSEL_SUPPORT_HPP_INCLUDED

#include <cstdint>
#include <algorithm>
#include <omp.h>

#include "bitvector.hpp"
#include "utils.hpp"


namespace psascan_private {

template<
  typename offset_type = std::uint64_t,
  std::uint64_t k_sampling_rate_log = 20>
class ranksel_support {

  private:
    static const std::uint64_t k_sampling_rate;
    static const std::uint64_t k_sampling_rate_mask;

    const bitvector * const m_bv;
    const std::uint64_t m_length;
    std::uint64_t n_blocks;
    std::uint64_t m_total_ones;
    offset_type *m_sparse_rank;

  public:

    //=========================================================================
    // Constructor, length is the size of bv in bits.
    //=========================================================================
    ranksel_support(
        const bitvector * const bv,
        const std::uint64_t length) :
          m_bv(bv),
          m_length(length) {

      // Allocate sparse rank.
      n_blocks = (m_length + k_sampling_rate - 1) >> k_sampling_rate_log;
      m_sparse_rank = utils::allocate_array<offset_type>(n_blocks);

      // Compute sparse rank.
      {

#ifdef _OPENMP

        // Compute the number of groups.
        const std::uint64_t max_threads = omp_get_max_threads();
        const std::uint64_t max_group_size =
          (n_blocks + max_threads - 1) / max_threads;
        const std::uint64_t n_groups =
          (n_blocks + max_group_size - 1) / max_group_size;

        // Each thread handles a group of samples.
        #pragma omp parallel num_threads(n_groups)
        {

          // Compute the range of blocks to handle.
          const std::uint64_t group_id = omp_get_thread_num();
          const std::uint64_t group_beg = group_id * max_group_size;
          const std::uint64_t group_end = std::min(n_blocks,
              group_beg + max_group_size);

          // Compute the sum of 1-bits inside
          // each of the handled blocks.
          for (std::uint64_t block_id = group_beg;
              block_id < group_end; ++block_id) {

            // Compute block boundaries.
            const std::uint64_t block_beg =
              (block_id << k_sampling_rate_log);
            const std::uint64_t block_end =
              std::min(m_length, block_beg + k_sampling_rate);

            // Store a number of 1-bits in the block.
            m_sparse_rank[block_id] =
              m_bv->range_sum(block_beg, block_end);
          }
        }
#else

        // Sequential version.
        for (std::uint64_t block_id = 0;
            block_id < n_blocks; ++block_id) {

          // Compute block boundaries.
          const std::uint64_t block_beg =
            (block_id << k_sampling_rate_log);
          const std::uint64_t block_end =
            std::min(m_length, block_beg + k_sampling_rate);

          // Store a number of 1-bits in the block.
          m_sparse_rank[block_id] =
            m_bv->range_sum(block_beg, block_end);
        }

#endif  // _OPENMP
      }

      // Exclusive partial sum over m_sparse_rank.
      std::uint64_t bit_count = 0;
      for (std::uint64_t block_id = 0;
          block_id < n_blocks; ++block_id) {
        const std::uint64_t temp = m_sparse_rank[block_id];
        m_sparse_rank[block_id] = bit_count;
        bit_count += temp;
      }

      // Store the total number of ones.
      m_total_ones = bit_count;
    }

    //=========================================================================
    // Find the largest position j such that the number of 0s in bv[0..j)
    // is <= i. In other words, find the position of (i+1)-th 0-bit in bv
    // (i = 0, 1, ..). 0 <= i < number of 0-bits in bv.
    //==========================================================================
    inline std::uint64_t select0(const std::uint64_t i) const {

      // Fast-forward through blocks
      // preceding the block with answer.
      std::uint64_t block_id = 0;
      std::uint64_t block_beg = 0;
      while (block_id + 2 < n_blocks &&
          (block_beg + k_sampling_rate) -
          (std::uint64_t)m_sparse_rank[block_id + 1] <= i) {
        ++block_id;
        block_beg += k_sampling_rate;
      }

      // Compute the final position
      // in a single block.
      const std::uint64_t zero_count =
        block_beg - (std::uint64_t)m_sparse_rank[block_id];
      return m_bv->select0(block_beg, i - zero_count);
    }

    //=========================================================================
    // Find the largest position j such that the number of 1s in bv[0..j)
    // is <= i. In other words, find the position of (i+1)-th 1-bit in bv
    // (i = 0, 1, ..). 0 <= i < number of 1-bits in bv.
    //==========================================================================
    inline std::uint64_t select1(const std::uint64_t i) const {

      // Fast-forward through blocks
      // preceding the block with answer.
      std::uint64_t block_id = 0;
      std::uint64_t block_beg = 0;
      while (block_id + 2 < n_blocks &&
          (std::uint64_t)m_sparse_rank[block_id + 1] <= i) {
        ++block_id;
        block_beg += k_sampling_rate;
      }

      // Compute the final position
      // in a single block.
      const std::uint64_t one_count = m_sparse_rank[block_id];
      return m_bv->select1(block_beg, i - one_count);
    }

    //==========================================================================
    // Compute the number of 1-bits in bv[0..i),
    // i is in range [0..m_length].
    //=========================================================================
    inline std::uint64_t rank(const std::uint64_t i) const {
      if (i == m_length)
        return m_total_ones;

      // Compute the rank from the sample array.
      const std::uint64_t block_id = (i >> k_sampling_rate_log);
      const std::uint64_t ones_count = m_sparse_rank[block_id];

      // Add the number of 1s in the last block
      // and return the answer.
      const std::uint64_t block_beg = (block_id << k_sampling_rate_log);
      return ones_count + m_bv->range_sum(block_beg, i);
    }

    //=========================================================================
    // Compute the number of 0-bits in bv[0..i),
    // i is in range [0..m_length].
    //=========================================================================
    inline std::uint64_t rank0(const std::uint64_t i) const {
      return i - rank(i);
    }

    //=========================================================================
    // Destructor.
    //=========================================================================
    ~ranksel_support() {
      utils::deallocate(m_sparse_rank);
    }
};

template<typename offset_type, std::uint64_t k_sampling_rate_log>
  const std::uint64_t ranksel_support<offset_type, k_sampling_rate_log>
  ::k_sampling_rate = ((std::uint64_t)1 << k_sampling_rate_log);

}  // psascan_private

#endif  // __SRC_PSASCAN_SRC_RANKSEL_SUPPORT_HPP_INCLUDED
