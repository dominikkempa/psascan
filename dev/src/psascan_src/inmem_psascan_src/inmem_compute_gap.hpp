/**
 * @file    src/psascan_src/inmem_psascan_src/inmem_compute_gap.hpp
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

#ifndef __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_COMPUTE_GAP_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_COMPUTE_GAP_HPP_INCLUDED

#include <cstdio>
#include <cstdint>
#include <map>
#include <vector>
#include <thread>
#include <algorithm>

#include "../bitvector.hpp"
#include "../gap_buffer.hpp"
#include "../io/multifile.hpp"
#include "rank.hpp"
#include "inmem_gap_array.hpp"
#include "inmem_compute_initial_ranks.hpp"
#include "inmem_stream.hpp"
#include "inmem_update.hpp"
#include "inmem_bwt_from_sa.hpp"
#include "pagearray.hpp"
#include "bwtsa.hpp"
#include "space_efficient_isa.hpp"


namespace psascan_private {
namespace inmem_psascan_private {

template<typename block_offset_type, unsigned pagesize_log>
void inmem_compute_gap(
    const std::uint8_t *text,
    std::uint64_t text_length,
    std::uint64_t left_block_beg,
    std::uint64_t left_block_size,
    std::uint64_t right_block_size,
    const pagearray<bwtsa_t<block_offset_type>, pagesize_log> &bwtsa,
    bitvector *gt,
    inmem_gap_array* &gap,
    std::uint64_t max_threads,
    bool need_gt,
    std::uint64_t i0,
    std::uint64_t gap_buf_size,
    long double &rank_init_time,
    long double &streaming_time,
    std::uint64_t **block_rank_matrix,
    std::uint64_t lrange_beg,
    std::uint64_t lrange_size,
    std::uint64_t rrange_size) {

  std::uint64_t lrange_end = lrange_beg + lrange_size;
  std::uint64_t rrange_end = lrange_end + rrange_size;

  //---------------------------------------------------------------------------
  // STEP 1: build rank data structure over BWT.
  //---------------------------------------------------------------------------
  fprintf(stderr, "    Build rank: ");
  long double start = utils::wclock();
  typedef rank4n<block_offset_type, pagesize_log> rank_type;
  rank_type *rank = new rank_type(&bwtsa, left_block_size, max_threads);
  rank_init_time = utils::wclock() - start;
  fprintf(stderr, "total: %.2Lfs\n", rank_init_time);

  //---------------------------------------------------------------------------
  // STEP 2: compute symbol counts and the last symbol of the left block.
  //---------------------------------------------------------------------------
  static const std::uint64_t k_sigma = 256;
  std::uint64_t *count = new std::uint64_t[k_sigma];
  for (std::uint64_t c = 0; c < k_sigma; ++c)
    count[c] = rank->query(left_block_size, (std::uint8_t)c);

  const std::uint8_t *left_block = text + left_block_beg;
  std::uint8_t last = left_block[left_block_size - 1];
  ++count[last];
  --count[0];

  for (std::uint64_t i = 0, sum = 0, temp = 0; i < k_sigma; ++i) {
    temp = count[i];
    count[i] = sum;
    sum += temp;
  }

  //---------------------------------------------------------------------------
  // STEP 3: compute starting positions for all streaming threads.
  //---------------------------------------------------------------------------
  std::uint64_t left_block_end = left_block_beg + left_block_size;
  std::uint64_t right_block_beg = left_block_end;
  std::uint64_t right_block_end = left_block_end + right_block_size;

  std::uint64_t max_stream_block_size =
    (right_block_size + max_threads - 1) / max_threads;
  while (max_stream_block_size & 7) ++max_stream_block_size;
  std::uint64_t n_threads =
    (right_block_size + max_stream_block_size - 1) / max_stream_block_size;

  fprintf(stderr, "    Compute initial ranks: ");
  start = utils::wclock();
  std::vector<std::uint64_t> initial_ranks(n_threads);
  std::vector<std::pair<std::uint64_t, std::uint64_t> >
    initial_ranges(n_threads);
  std::thread **threads = new std::thread*[n_threads];

  // 3.a
  //
  // Compute the last starting position using the matrix of initial ranks.
  typedef bwtsa_t<block_offset_type> bwtsa_type;
  typedef pagearray<bwtsa_type, pagesize_log> pagearray_bwtsa_type;
  std::uint64_t last_stream_block_beg =
    right_block_beg + (n_threads - 1) * max_stream_block_size;
  std::uint64_t last_stream_block_end = right_block_end;

  initial_ranks[n_threads - 1] = 0;
  for (std::uint64_t j = lrange_beg; j < lrange_end; ++j)
    initial_ranks[n_threads - 1] += block_rank_matrix[j][rrange_end - 1];

  // 3.b
  //
  // Compute the starting position for all
  // starting positions other than the last one.
  std::uint64_t prev_stream_block_size =
    last_stream_block_end - last_stream_block_beg;
  for (std::uint64_t i_plusplus = n_threads; i_plusplus > 1; --i_plusplus) {
    std::uint64_t i = i_plusplus - 2;
    std::uint64_t stream_block_beg =
      right_block_beg + i * max_stream_block_size;
    std::uint64_t stream_block_end =
      std::min(stream_block_beg + max_stream_block_size, right_block_end);
    std::uint64_t stream_block_size = stream_block_end - stream_block_beg;
    const std::uint8_t *pat = text + stream_block_end;

    threads[i] = new std::thread(compute_range<pagearray_bwtsa_type>,
        text, left_block_beg, left_block_size, pat, prev_stream_block_size,
        std::ref(bwtsa), std::ref(initial_ranges[i]));

    prev_stream_block_size = stream_block_size;
  }

  for (std::uint64_t i = 0; i + 1 < n_threads; ++i) threads[i]->join();
  for (std::uint64_t i = 0; i + 1 < n_threads; ++i) delete threads[i];
  delete[] threads;
  fprintf(stderr, "%.2Lfs ", utils::wclock() - start);

  bool nontrivial_range = false;
  for (std::uint64_t j = 0; j + 1 < n_threads; ++j)
    if (initial_ranges[j].first != initial_ranges[j].second)
      nontrivial_range = true;

  if (nontrivial_range) {

    // 3.c
    //
    // Build the data structure allowing answering ISA queries.
    start = utils::wclock();
    typedef pagearray<bwtsa_type, pagesize_log> pagearray_type;

#ifdef INMEM_PSASCAN_DEBUG
    typedef space_efficient_isa<pagearray_type, rank_type, 2> isa_type;
#else
    typedef space_efficient_isa<pagearray_type, rank_type, 12> isa_type;
#endif

    isa_type *sp_isa = new isa_type(&bwtsa,
        text + left_block_beg, rank, left_block_size, i0);
    fprintf(stderr, "%.3Lfs ", utils::wclock() - start);

    // 3.d
    //
    // Narrow nontrivial ranges to single elements.
    start = utils::wclock();
    prev_stream_block_size = last_stream_block_end - last_stream_block_beg;
    std::uint64_t prev_rank = initial_ranks[n_threads - 1];
    for (std::uint64_t i_plus = n_threads - 1; i_plus > 0; --i_plus) {
      std::uint64_t i = i_plus - 1;
      std::uint64_t stream_block_beg =
        right_block_beg + i * max_stream_block_size;
      std::uint64_t stream_block_end =
        std::min(stream_block_beg + max_stream_block_size, right_block_end);
      std::uint64_t stream_block_size = stream_block_end - stream_block_beg;
      std::uint64_t suf_start = stream_block_end;

      std::uint64_t left = initial_ranges[i].first;
      std::uint64_t right = initial_ranges[i].second;

      // Keep refining the range [left..right) until it's empty.
      while (left != right) {

        // Valid values for mid are in [left..right).
        std::uint64_t mid = (left + right) / 2;

        // Check if suffix starting at position suf_start
        // is larger than the one starting at block_beg +
        // bwtsa[mid].sa in the text. We know they have
        // a common prefix of length prev_stream_block_size.
        if ((std::uint64_t)bwtsa[mid].m_sa +
            prev_stream_block_size >= left_block_size) {
          if (gt->get(text_length - 1 -
                (suf_start + left_block_size -
                 (std::uint64_t)bwtsa[mid].m_sa - 1)))
            left = mid + 1;
          else right = mid;
        } else {
          std::uint64_t j = bwtsa[mid].m_sa + prev_stream_block_size;
          if (sp_isa->query(j) < prev_rank)
            left = mid + 1;
          else right = mid;
        }
      }

      initial_ranks[i] = left;
      prev_rank = left;
      prev_stream_block_size = stream_block_size;
    }

    delete sp_isa;
    fprintf(stderr, "%.3Lfs ", utils::wclock() - start);
  } else {
    for (std::uint64_t j = 0; j + 1 < n_threads; ++j)
      initial_ranks[j] = initial_ranges[j].first;
  }
  fprintf(stderr, "\n");

  //---------------------------------------------------------------------------
  // STEP 4: allocate gap array. The gap array is indexed from 0 to
  //         left_block_size so the number of elements is left_block_size + 1.
  //---------------------------------------------------------------------------
  start = utils::wclock();
  gap = new inmem_gap_array(left_block_size + 1);

  //---------------------------------------------------------------------------
  // STEP 5: allocate buffers, buffer polls and auxiliary arrays.
  //---------------------------------------------------------------------------

  // Allocate gap buffers.
  std::uint64_t n_gap_buffers = 2 * n_threads;
  gap_buffer<block_offset_type> **gap_buffers =
    new gap_buffer<block_offset_type>*[n_gap_buffers];
  for (std::uint64_t i = 0; i < n_gap_buffers; ++i)
    gap_buffers[i] =
      new gap_buffer<block_offset_type>(gap_buf_size, max_threads);

  // Create poll of empty and full buffers.
  gap_buffer_poll<block_offset_type> *empty_gap_buffers =
    new gap_buffer_poll<block_offset_type>();
  gap_buffer_poll<block_offset_type> *full_gap_buffers =
    new gap_buffer_poll<block_offset_type>(n_threads);

  // Add empty buffers to empty poll.
  for (std::uint64_t i = 0; i < n_gap_buffers; ++i)
    empty_gap_buffers->add(gap_buffers[i]);

  // Allocate temp arrays and oracles.
  std::uint64_t max_buffer_elems = gap_buf_size / sizeof(block_offset_type);
  block_offset_type *temp =
    utils::allocate_array<block_offset_type>(max_buffer_elems * n_threads);
  std::uint32_t *oracle =
    utils::allocate_array<std::uint32_t>(max_buffer_elems * n_threads);
  long double allocations_time = utils::wclock() - start;
  if (allocations_time > 0.05L)
    fprintf(stderr, "    Allocations: %.2Lfs\n", allocations_time);

  //---------------------------------------------------------------------------
  // STEP 6: run the parallel streaming.
  //---------------------------------------------------------------------------

  // Start streaming threads.
  fprintf(stderr, "    Stream: ");
  start = utils::wclock();
  threads = new std::thread*[n_threads];
  for (std::uint64_t t = 0; t < n_threads; ++t) {
    std::uint64_t beg = right_block_beg + t * max_stream_block_size;
    std::uint64_t end =
      std::min(beg + max_stream_block_size, right_block_end);

    threads[t] = new std::thread(
        inmem_parallel_stream<rank_type, block_offset_type>,
        text, text_length, beg, end, last, count, full_gap_buffers,
        empty_gap_buffers, initial_ranks[t], i0, rank, gap->m_length,
        max_threads, gt, temp + t * max_buffer_elems,
        oracle + t * max_buffer_elems, need_gt);
  }

  // Start updating thread.
  std::thread *updater =
    new std::thread(inmem_gap_updater<block_offset_type>,
        full_gap_buffers, empty_gap_buffers, gap, max_threads);

  // Wait to all threads to finish.
  for (std::uint64_t t = 0; t < n_threads; ++t) threads[t]->join();
  updater->join();
  streaming_time = utils::wclock() - start;
  long double streaming_speed =
    (right_block_size / (1024.L * 1024)) / streaming_time;
  fprintf(stderr, "%.2Lfs (%.2LfMiB/s)\n", streaming_time,
      streaming_speed);

  //---------------------------------------------------------------------------
  // STEP 7: clean up and sort gap->m_excess.
  // XXX Consider using gnu parallel sort.
  //---------------------------------------------------------------------------
  start = utils::wclock();
  utils::deallocate(oracle);
  utils::deallocate(temp);
  for (std::uint64_t i = 0; i < n_threads; ++i) delete threads[i];
  for (std::uint64_t i = 0; i < n_gap_buffers; ++i) delete gap_buffers[i];
  delete updater;
  delete[] threads;
  delete[] gap_buffers;
  delete empty_gap_buffers;
  delete full_gap_buffers;
  delete rank;
  delete[] count;

  std::sort(gap->m_excess.begin(), gap->m_excess.end());

  long double cleaning_time = utils::wclock() - start;
  if (cleaning_time > 0.1L)
    fprintf(stderr, "    Clean: %.2Lfs\n", cleaning_time);
}

}  // namespace inmem_psascan_private
}  // namespace psascan_private
                 
#endif  // __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_COMPUTE_GAP_HPP_INCLUDED
