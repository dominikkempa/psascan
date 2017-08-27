/**
 * @file    src/psascan_src/compute_gap.hpp
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

#ifndef __SRC_PSASCAN_SRC_COMPUTE_GAP_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_COMPUTE_GAP_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <thread>
#include <algorithm>
#include <vector>

#include "utils.hpp"
#include "rank.hpp"
#include "gap_array.hpp"
#include "gap_buffer.hpp"
#include "stream.hpp"
#include "update.hpp"
#include "stream_info.hpp"
#include "io/multifile.hpp"
#include "io/async_multifile_bit_writer.hpp"


namespace psascan_private {

//==============================================================================
// Compute the gap for an arbitrary range of suffixes of tail. This version is
// more general, and can be used also when processing half-blocks.
//==============================================================================
template<typename block_offset_type,
         typename rank_type>
void compute_gap(
    const rank_type *rank,
    std::uint64_t block_size,
    buffered_gap_array *gap,
    std::uint64_t tail_begin,
    std::uint64_t tail_end,
    std::uint64_t text_length,
    std::uint64_t max_threads,
    std::uint64_t block_isa0,
    std::uint64_t gap_buf_size,
    std::uint8_t block_last_symbol,
    std::vector<std::uint64_t> initial_ranks,
    std::string text_filename,
    std::string output_filename,
    const multifile *tail_gt_begin_rev,
    multifile *newtail_gt_begin_rev) {

  std::uint64_t tail_length = tail_end - tail_begin;
  std::uint64_t stream_max_block_size =
    (tail_length + max_threads - 1) / max_threads;
  std::uint64_t n_threads =
    (tail_length + stream_max_block_size - 1) / stream_max_block_size;

  fprintf(stderr, "    Stream:");
  long double stream_start = utils::wclock();

  // 1
  //
  // Get symbol counts of a block and turn into exclusive partial sum.
  static const std::uint64_t k_sigma = 256;
  std::uint64_t *count = new std::uint64_t[k_sigma];
  for (std::uint64_t j = 0; j < k_sigma; ++j)
    count[j] = rank->rank(block_size, (std::uint8_t)j);
  ++count[block_last_symbol];
  --count[0];

  // Exclusive partial sum over the count array.
  for (std::uint64_t j = 0, sum = 0, temp = 0; j < k_sigma; ++j) {
    temp = count[j];
    count[j] = sum;
    sum += temp;
  }

  // 2
  //
  // Allocate gap buffers.
  std::uint64_t n_gap_buffers = 2 * max_threads;
  gap_buffer<block_offset_type> **gap_buffers =
    new gap_buffer<block_offset_type>*[n_gap_buffers];
  for (std::uint64_t i = 0; i < n_gap_buffers; ++i)
    gap_buffers[i] =
      new gap_buffer<block_offset_type>(gap_buf_size, max_threads);

  // 3
  //
  // Create poll of empty and full buffers.
  gap_buffer_poll<block_offset_type> *empty_gap_buffers =
    new gap_buffer_poll<block_offset_type>();
  gap_buffer_poll<block_offset_type> *full_gap_buffers =
    new gap_buffer_poll<block_offset_type>(n_threads);

  // 4
  //
  // Add all buffers to the poll of empty buffers.
  for (std::uint64_t i = 0; i < n_gap_buffers; ++i)
    empty_gap_buffers->add(gap_buffers[i]);

  // Create the async multifile bit writer.
  typedef async_multifile_bit_writer gt_writer_type;
  gt_writer_type *gt_bit_writer = new gt_writer_type();
  for (std::uint64_t t = 0; t < n_threads; ++t) {
    std::uint64_t stream_block_beg = tail_begin + t * stream_max_block_size;
    std::uint64_t stream_block_end = std::min(tail_end,
        stream_block_beg + stream_max_block_size);

    std::string filename = output_filename +
      ".gt_tail." + utils::random_string_hash();
    newtail_gt_begin_rev->add_file(
        text_length - stream_block_end,
        text_length - stream_block_beg, filename);
    gt_bit_writer->add_file(filename);
  }

  // 5
  //
  // Start threads doing the backward search.
  stream_info info(n_threads, tail_length);
  std::thread **streamers = new std::thread*[n_threads];

  for (std::uint64_t t = 0; t < n_threads; ++t) {
    std::uint64_t stream_block_beg = tail_begin + t * stream_max_block_size;
    std::uint64_t stream_block_end = std::min(tail_end,
      stream_block_beg + stream_max_block_size);

    streamers[t] = new std::thread(
        parallel_stream<block_offset_type, rank_type>,
        full_gap_buffers, empty_gap_buffers, stream_block_beg,
        stream_block_end, initial_ranks[t], count, block_isa0,
        rank, block_last_symbol, text_filename, text_length,
        &info, t, gap->m_length, gap_buf_size, tail_gt_begin_rev,
        max_threads, gt_bit_writer);
  }

  // 6
  //
  // Start threads doing the gap array updates.
  std::thread *updater = new std::thread(gap_updater<block_offset_type>,
        full_gap_buffers, empty_gap_buffers, gap, max_threads);

  // 7
  //
  // Wait for all threads to finish.
  for (std::uint64_t i = 0; i < n_threads; ++i)
    streamers[i]->join();
  updater->join();

  // 8
  //
  // Clean up.
  for (std::uint64_t i = 0; i < n_threads; ++i) delete streamers[i];
  for (std::uint64_t i = 0; i < n_gap_buffers; ++i) delete gap_buffers[i];
  delete updater;
  delete[] streamers;
  delete[] gap_buffers;
  delete empty_gap_buffers;
  delete full_gap_buffers;
  delete[] count;
  delete gt_bit_writer;

  // 9
  //
  // Print summary and exit.
  long double stream_time = utils::wclock() - stream_start;
  long double speed = (tail_length / (1024.L * 1024)) / stream_time;
  fprintf(stderr,"\r    Stream: 100.0%%. Time: %.2Lfs. Speed: %.2LfMiB/s\n",
      stream_time, speed);
}

}  // namespace psascan_private

#endif  // __SRC_PSASCAN_SRC_COMPUTE_GAP_HPP_INCLUDED
