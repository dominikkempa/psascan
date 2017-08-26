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

#include <thread>
#include <algorithm>

#include "bitvector.hpp"


namespace psascan_private {

class ranksel_support {
  private:
    std::uint64_t m_length;      // length of bitvector
    std::uint64_t m_chunk_size;  // chunk size
    std::uint64_t n_chunks;      // number of chunks
    std::uint64_t *m_sparse_rank;

    const bitvector *m_bv;

  private:

    //=========================================================================
    // Compute sparse_rank[group_beg..group_end).
    //=========================================================================
    static void process_group_of_chunks(
        std::uint64_t group_beg,
        std::uint64_t group_end,
        std::uint64_t chunk_size,
        std::uint64_t *sparse_rank,
        const bitvector *bv) {

      for (std::uint64_t chunk_id = group_beg;
          chunk_id < group_end; ++chunk_id) {
        std::uint64_t chunk_beg = chunk_id * chunk_size;
        std::uint64_t chunk_end = chunk_beg + chunk_size;

        sparse_rank[chunk_id] = bv->range_sum(chunk_beg, chunk_end);
      }
    }


  public:

    //=========================================================================
    // Constructor.
    //=========================================================================
    ranksel_support(
        const bitvector *bv,
        std::uint64_t length,
        std::uint64_t max_threads) {

      m_bv = bv;
      m_length = length;

      // Compute chunk size and allocate m_sparse_rank.
      // We excluse the last partial chunk.
      m_chunk_size = std::min(((std::uint64_t)1 << 20),
          (m_length + max_threads - 1) / max_threads);
      n_chunks = m_length / m_chunk_size;
      m_sparse_rank = utils::allocate_array<std::uint64_t>(n_chunks + 1);

      // Compute the sum of 1-bits inside each chunk and write to
      // m_sparse_rank. Since there can be more chunks than threads,
      // we split chunks into groups and let each thread handle the
      // group of chunks.
      std::uint64_t chunk_max_group_size =
        (n_chunks + max_threads - 1) / max_threads;
      std::uint64_t n_chunk_groups =
        (n_chunks + chunk_max_group_size - 1) / chunk_max_group_size;

      std::thread **threads = new std::thread*[n_chunk_groups];
      for (std::uint64_t t = 0; t < n_chunk_groups; ++t) {
        std::uint64_t chunk_group_beg = t * chunk_max_group_size;
        std::uint64_t chunk_group_end =
          std::min(chunk_group_beg + chunk_max_group_size, n_chunks);

        threads[t] =
          new std::thread(process_group_of_chunks, chunk_group_beg,
              chunk_group_end, m_chunk_size, m_sparse_rank, m_bv);
      }

      // Stop all threads and clean up.
      for (std::uint64_t t = 0; t < n_chunk_groups; ++t) threads[t]->join();
      for (std::uint64_t t = 0; t < n_chunk_groups; ++t) delete threads[t];
      delete[] threads;
    
      // Partial exclusive sum on m_sparse_rank.
      std::uint64_t ones = 0;
      for (std::uint64_t i = 0; i < n_chunks; ++i) {
        std::uint64_t temp = m_sparse_rank[i];
        m_sparse_rank[i] = ones;
        ones += temp;
      }
      m_sparse_rank[n_chunks] = ones;
    }

    //=========================================================================
    // Find the largest position j such that the number of 0s in bv[0..j)
    // is <= i. In other words, find the position of i-th 0-bit in bv
    // (i = 0, 1, ..). 0 <= i < number of 0-bits in bv.
    //==========================================================================
    inline std::uint64_t select0(std::uint64_t i) const {

      // Fast-forward through chunks preceding
      // the chunk with the answer.
      std::uint64_t j = 0;
      while (j < n_chunks &&
        ((j + 1) * m_chunk_size) - m_sparse_rank[j + 1] <= i) ++j;
      std::uint64_t zero_cnt_j =
        (j * m_chunk_size) - m_sparse_rank[j];
      j *= m_chunk_size;

      // Find the final position in a single chunk.
      while (zero_cnt_j + (1 - m_bv->get(j)) <= i)
        zero_cnt_j += (1 - m_bv->get(j++));

      // Return result.
      return j;
    }

    //=========================================================================
    // Find the largest position j such that the number of 1s in bv[0..j)
    // is <= i In other words, find the position of i-th 1-bit in bv
    // (i = 0, 1, ..). 0 <= i < number of 1-bits in bv.
    //=========================================================================
    inline std::uint64_t select1(std::uint64_t i) const {

      // Fast-forward through chunks preceding
      // the chunk with the answer.
      std::uint64_t j = 0;
      while (j < n_chunks &&
          m_sparse_rank[j + 1] <= i) ++j;
      std::uint64_t rank_j = m_sparse_rank[j];
      j *= m_chunk_size;

      // Find the final position in a single chunk.
      while (rank_j + m_bv->get(j) <= i)
        rank_j += m_bv->get(j++);

      // Return result.
      return j;
    }
  
    //==========================================================================
    // Compute the number of 1-bits in bv[0..i) with the help of sparse_rank.
    // Note:
    // - i is an integer in the range from 0 to length of bv (inclusive),
    // - sparse_rank[k] = number of 1-bits in bv[0..k * chunk_size).
    //=========================================================================
    inline std::uint64_t rank(std::uint64_t i) const {

      // Compute the rank from the sample array.
      std::uint64_t j = i / m_chunk_size;
      std::uint64_t result = m_sparse_rank[j];

      // Scan the chunk.
      j *= m_chunk_size;
      while (j < i)
        result += m_bv->get(j++);

      // Return result.
      return result;
    }

    //=========================================================================
    // Compute the number of 0-bits in bv[0..i).
    // 0 <= i <= m_length.
    //=========================================================================
    inline std::uint64_t rank0(std::uint64_t i) const {
      return i - rank(i);
    }

    ~ranksel_support() {
      utils::deallocate(m_sparse_rank);
    }
};

}  // psascan_private

#endif  // __SRC_PSASCAN_SRC_RANKSEL_SUPPORT_HPP_INCLUDED
