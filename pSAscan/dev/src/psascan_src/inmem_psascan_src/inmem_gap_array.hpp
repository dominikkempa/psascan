/**
 * @file    src/psascan_src/inmem_psascan_src/inmem_gap_array.hpp
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

#ifndef __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_SASCAN_INMEM_GAP_ARRAY_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_SASCAN_INMEM_GAP_ARRAY_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <mutex>
#include <stack>
#include <thread>


namespace psascan_private {
namespace inmem_psascan_private {

class inmem_gap_array {
  public:
    std::uint8_t *m_count;
    std::uint64_t m_length;

    std::vector<std::uint64_t> m_excess;
    std::mutex m_excess_mutex;

  public:
    inmem_gap_array(std::uint64_t length)
      : m_length(length) {
      m_count = (std::uint8_t *)calloc(m_length, 1);
    }

    ~inmem_gap_array() {
      free(m_count);
    }

    //=========================================================================
    // Find and smallest j such that j + gap[0] + .. + gap[j] >= a. Store
    // the value of j into b and gap[0] + .. + gap[j] into c. To speed up the
    // algorithm, we have array gapsum defined as
    //
    //    gapsum[i] = gap[0] + .. + gap[i * block_size - 1].
    //
    //=========================================================================
    static void answer_single_gap_query(
        const inmem_gap_array *gap,
        std::uint64_t block_size,
        const std::uint64_t *gapsum,
        std::uint64_t a,
        std::uint64_t &b,
        std::uint64_t &c) {

      std::uint64_t n_blocks =
        (gap->m_length + block_size - 1) / block_size;

      // Find the block containing the correct index.
      // To do that find the largest j such that
      // gapsum[j] + block_size * j - 1 < a and start
      // searching from j * block_size.
      std::uint64_t j = 0;
      while (j + 1 < n_blocks &&
          gapsum[j + 1] + block_size * (j + 1) < a + 1) ++j;

      // Invariant: the j we are searching
      // for is > j * block_size - 1.

      std::uint64_t sum = gapsum[j];
      j = block_size * j;
      std::uint64_t excess_ptr =
        std::lower_bound(gap->m_excess.begin(),
            gap->m_excess.end(), j) - gap->m_excess.begin();
      while (true) {

        // Invariant: sum = gap[0] + .. + gap[j - 1].
        // Compute gap[j] using small gap array
        // representation.
        std::uint64_t gap_j = gap->m_count[j];
        while (excess_ptr < gap->m_excess.size() &&
            gap->m_excess[excess_ptr] == j) {
          gap_j += 256;
          ++excess_ptr;
        }

        if (j + sum + gap_j >= a) {
          b = j;
          c = sum + gap_j;
          return;
        } else {
          sum += gap_j;
          ++j;
        }
      }
    }

    //=========================================================================
    // Compute gap[0] + gap[1] + .. + gap[j - 1] with the help of gapsum array.
    //=========================================================================
    static std::uint64_t compute_sum3(
        const inmem_gap_array *gap,
        std::uint64_t j,
        std::uint64_t max_block_size,
        std::uint64_t *gapsum) {

      // Invariant: j > 0.
      if (j == 0) {
        fprintf(stderr, "\nError: j <= 0 in compute_sum3\n");
        std::exit(EXIT_FAILURE);
      }

      std::uint64_t block_id = j / max_block_size;
      std::uint64_t result = gapsum[block_id];
      std::uint64_t scan_beg = block_id * max_block_size;
      std::uint64_t scan_end = j;
      std::int64_t occ =  // XXX really needed signed int?
        std::upper_bound(
            gap->m_excess.begin(),
            gap->m_excess.end(),
            scan_end - 1) -
        std::lower_bound(
            gap->m_excess.begin(),
            gap->m_excess.end(),
            scan_beg);

      result +=
        (std::uint64_t)256 *
        (std::uint64_t)std::max((std::int64_t)0, occ);

      for (std::uint64_t i = block_id * max_block_size; i < j; ++i)
        result += gap->m_count[i];

      return result;
    }

    //=========================================================================
    // Compute sum of gap values for blocks in range [range_beg..range_end).
    // The sum for each block is stored in gapsum array.
    //=========================================================================
    static void compute_sum2(
        const inmem_gap_array *gap,
        std::uint64_t range_beg,
        std::uint64_t range_end,
        std::uint64_t max_block_size,
        std::uint64_t *gapsum) {

      for (std::uint64_t block_id = range_beg;
          block_id < range_end; ++block_id) {
        std::uint64_t block_beg = block_id * max_block_size;
        std::uint64_t block_end =
          std::min(block_beg + max_block_size, gap->m_length);

        // Invariant: block_end > 0.
        if (block_end == 0) {
          fprintf(stderr, "\nblock_end <= 0 in compute_sum2\n");
          std::exit(EXIT_FAILURE);
        }

        // Process block.
        std::int64_t occ =
          std::upper_bound(
              gap->m_excess.begin(),
              gap->m_excess.end(),
              block_end - 1) -
          std::lower_bound(
              gap->m_excess.begin(),
              gap->m_excess.end(),
              block_beg);
        std::uint64_t block_gap_sum =
          (std::uint64_t)256 *
          (std::uint64_t)std::max((std::int64_t)0, occ);
        for (std::uint64_t j = block_beg; j < block_end; ++j)
          block_gap_sum += gap->m_count[j];
        gapsum[block_id] = block_gap_sum;
      }
    }

    //=========================================================================
    // Parallel computaton of answers to n_queries queries of the form:
    // What is the smallest j such that j + gap[0] + .. + gap[j] >= a[i]"
    //   - the answer to i-th query is stored in b[i]
    //   - in addition we also return gap[0] + .. + gap[j] in c[i]
    //
    // To do that we first split the gap array into blocks of size of
    // about length / max_threads and (in parallel) compute sums of gap
    // values inside these blocks. We the accumulate these sums into
    // array of prefix sums.
    //
    // To answer each of the queries we start a separate thread. Each
    // thread uses the partial sums of gap array at block boundaries to
    // find a good starting point for search and then scans the gap array
    // from there.
    //=========================================================================
    std::uint64_t answer_queries(
        std::uint64_t n_queries,
        const std::uint64_t *a,
        std::uint64_t *b,
        std::uint64_t *c,
        std::uint64_t max_threads,
        std::uint64_t i0) const {

      //-----------------------------------------------------------------------
      // STEP 1: split gap array into at most max_threads blocks
      // and in parallel compute sum of values inside each block.
      //-----------------------------------------------------------------------
      std::uint64_t max_block_size = std::min(
          ((std::uint64_t)4) << 20,
          (m_length + max_threads - 1) / max_threads);
      std::uint64_t n_blocks =
        (m_length + max_block_size - 1) / max_block_size;
      std::uint64_t *gapsum = new std::uint64_t[n_blocks];
  
      // Each thread handles range of blocks.
      std::uint64_t range_size = (n_blocks + max_threads - 1) / max_threads;
      std::uint64_t n_ranges = (n_blocks + range_size - 1) / range_size;
      std::thread **threads = new std::thread*[max_threads];
      for (std::uint64_t range_id = 0; range_id < n_ranges; ++range_id) {
        std::uint64_t range_beg = range_id * range_size;
        std::uint64_t range_end = std::min(range_beg + range_size, n_blocks);

        threads[range_id] = new std::thread(compute_sum2, this,
            range_beg, range_end, max_block_size, gapsum);
      }
      for (std::uint64_t i = 0; i < n_ranges; ++i) threads[i]->join();
      for (std::uint64_t i = 0; i < n_ranges; ++i) delete threads[i];
      delete[] threads;

      // Compute partial sum from block counts.
      for (std::uint64_t i = 0, s = 0, t; i < n_blocks; ++i) {
        t = gapsum[i];
        gapsum[i] = s;
        s += t;
      }

      // Answer the queries in parallel.
      threads = new std::thread*[n_queries];
      for (std::uint64_t i = 0; i < n_queries; ++i)
        threads[i] = new std::thread(
            answer_single_gap_query,
            this, max_block_size, gapsum,
            a[i], std::ref(b[i]), std::ref(c[i]));

      for (std::uint64_t i = 0; i < n_queries; ++i) threads[i]->join();
      for (std::uint64_t i = 0; i < n_queries; ++i) delete threads[i];
      delete[] threads;

      std::uint64_t result =
        compute_sum3(this, i0 + 1, max_block_size, gapsum);

      delete[] gapsum;
      return result;
    }
};

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_GAP_ARRAY_HPP_INCLUDED
