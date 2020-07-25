/**
 * @file    src/psascan_src/stream.hpp
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

#ifndef __SRC_PSASCAN_SRC_STREAM_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_STREAM_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <mutex>
#include <algorithm>

#include "utils.hpp"
#include "rank.hpp"
#include "gap_buffer.hpp"
#include "update.hpp"
#include "stream_info.hpp"
#include "io/multifile.hpp"
#include "io/multifile_bit_stream_reader.hpp"
#include "io/async_scatterfile_bit_reader.hpp"
#include "io/async_backward_skip_stream_reader.hpp"
#include "io/async_multifile_bit_writer.hpp"


namespace psascan_private {

std::mutex stdout_mutex;

template<
  typename block_offset_type,
  typename rank_type>
void parallel_stream(
    const std::uint64_t stream_block_beg,
    const std::uint64_t stream_block_end,
    const std::uint64_t initial_rank,
    const std::uint64_t thread_id,
    const std::uint64_t gap_range_size,
    const std::uint64_t gap_buf_size,
    const std::uint64_t n_increasers,
    const std::uint64_t text_length,
    const std::uint64_t longest_block_suffix_rank,
    const std::uint8_t block_last_symbol,
    const std::string text_filename,
    const std::uint64_t * const block_symbol_count,
    const rank_type * const block_bwt_rank,
    const multifile * const tail_gt_begin,
    stream_info * const info,
    gap_buffer_poll<block_offset_type> * const full_gap_buffers,
    gap_buffer_poll<block_offset_type> * const empty_gap_buffers,
    async_multifile_bit_writer * const gt_bit_writer) {

  static const std::uint64_t max_buckets = 4096;
  std::uint32_t * const block_id_to_sblock_id = new std::uint32_t[max_buckets];

  std::uint64_t bucket_size = 1;
  std::uint64_t bucket_size_bits = 0;
  while ((gap_range_size + bucket_size - 1) / bucket_size > max_buckets) {
    bucket_size <<= 1;
    ++bucket_size_bits;
  }

  // Add description XXX.
  const std::uint64_t n_buckets =
    (gap_range_size + bucket_size - 1) / bucket_size;
  std::uint32_t * const block_count = new std::uint32_t[n_buckets];

  const std::uint64_t max_buffer_elems =
    gap_buf_size / sizeof(block_offset_type);
  block_offset_type * const temp =
    new block_offset_type[max_buffer_elems];
  std::uint32_t * const oracle =
    new std::uint32_t[max_buffer_elems];

  static const std::uint64_t buffer_sample_size = 512;
  std::vector<block_offset_type> samples(buffer_sample_size);
  std::uint64_t * const ptr =
    new std::uint64_t[n_increasers];
  block_offset_type * const bucket_lbound =
    new block_offset_type[n_increasers + 1];

  typedef async_scatterfile_bit_reader bit_stream_reader_type;
  typedef async_backward_skip_stream_reader<std::uint8_t> text_reader_type;

  text_reader_type *text_streamer = new text_reader_type(
      text_filename, text_length - stream_block_end, 4 << 20);
  bit_stream_reader_type gt_in(tail_gt_begin,
      text_length - stream_block_end, 1 << 20);

  std::uint64_t current_rank = initial_rank;
  std::uint64_t j = stream_block_end;
  std::uint64_t dbg = 0;
  while (j > stream_block_beg) {

    // Update streaming speed.
    if (dbg > ((std::uint64_t)1 << 26)) {
      info->m_mutex.lock();
      info->m_streamed[thread_id] = stream_block_end - j;
      info->m_update_count += 1;

      if (info->m_update_count == info->m_thread_count) {
        info->m_update_count = 0L;
        long double elapsed = utils::wclock() - info->m_timestamp;
        std::uint64_t total_streamed = 0L;

        for (long t = 0; t < info->m_thread_count; ++t)
          total_streamed += info->m_streamed[t];
        long double speed = (total_streamed / (1024.L * 1024)) / elapsed;

        stdout_mutex.lock();
        fprintf(stderr, "\r    Stream: %.2Lf%%. Time: %.2Lfs. "
            "Speed: %.2LfMiB/s",
            (total_streamed * 100.L) / info->m_tostream,
            elapsed, speed);
        stdout_mutex.unlock();
      }

      info->m_mutex.unlock();
      dbg = 0L;
    }

    // Get a gap buffer from the poll of empty buffers.
    std::unique_lock<std::mutex> lk(empty_gap_buffers->m_mutex);
    while (!empty_gap_buffers->available())
      empty_gap_buffers->m_cv.wait(lk);

    gap_buffer<block_offset_type> * const b = empty_gap_buffers->get();
    lk.unlock();
    empty_gap_buffers->m_cv.notify_one(); // let others know they should re-check

    // Process buffer -- fill with gap values.
    const std::uint64_t left = j - stream_block_beg;
    b->m_filled = std::min(left, b->m_size);
    dbg += b->m_filled;
    std::fill(block_count, block_count + n_buckets, (std::uint32_t)0);

    for (std::uint64_t t = 0; t < b->m_filled; ++t, --j) {
      const std::uint8_t c = text_streamer->read();
      const std::uint8_t gt_bit =
        (current_rank > longest_block_suffix_rank);

      gt_bit_writer->write_to_ith_file(thread_id, gt_bit);
      const bool next_gt = gt_in.read();

      const std::uint8_t delta = (c == 0 &&
          current_rank > longest_block_suffix_rank);

      current_rank =
        block_symbol_count[c] +
        block_bwt_rank->rank(current_rank, c);

      if (c == block_last_symbol && next_gt)
        ++current_rank;
      current_rank -= delta;

      temp[t] = current_rank;
      block_count[current_rank >> bucket_size_bits]++;
    }

    // Compute super-buckets.
    std::uint64_t ideal_sblock_size =
      (b->m_filled + n_increasers - 1) / n_increasers;
    std::uint64_t max_sbucket_size = 0;
    std::uint64_t bucket_id_beg = 0;
    for (std::uint64_t t = 0; t < n_increasers; ++t) {

      std::uint64_t bucket_id_end = bucket_id_beg;
      std::uint64_t size = 0;
      while (bucket_id_end < n_buckets && size < ideal_sblock_size)
        size += block_count[bucket_id_end++];

      b->sblock_size[t] = size;
      max_sbucket_size = std::max(max_sbucket_size, size);

      for (std::uint64_t id = bucket_id_beg; id < bucket_id_end; ++id)
        block_id_to_sblock_id[id] = t;
      bucket_id_beg = bucket_id_end;
    }

    if (max_sbucket_size < 4 * ideal_sblock_size) {
      for (std::uint64_t t = 0, curbeg = 0; t < n_increasers;
          curbeg += b->sblock_size[t++])
        b->sblock_beg[t] = ptr[t] = curbeg;

      // Permute the elements of the buffer.
      for (std::uint64_t t = 0; t < b->m_filled; ++t) {
        const std::uint64_t id = ((std::uint64_t)temp[t] >> bucket_size_bits);
        const std::uint64_t sblock_id = block_id_to_sblock_id[id];
        oracle[t] = ptr[sblock_id]++;
      }

      for (std::uint64_t t = 0; t < b->m_filled; ++t) {
        const std::uint64_t addr = oracle[t];
        b->m_content[addr] = temp[t];
      }
    } else {

      // Repeat the partition into sbuckets, this time using random sample.
      // This is a fallback mechanism in case the quick partition failed.
      // It is not suppose to happen to often.

      // Compute random sample of elements in the buffer.
      for (std::uint64_t t = 0; t < buffer_sample_size; ++t)
        samples[t] = temp[utils::random_int64(0, b->m_filled - 1)];
      std::sort(samples.begin(), samples.end());
      samples.erase(std::unique(samples.begin(),
            samples.end()), samples.end());

      // Compute bucket boundaries (lower bound is enough).
      std::fill(bucket_lbound,
          bucket_lbound + n_increasers + 1, gap_range_size);

      std::uint64_t step =
        (samples.size() + n_increasers - 1) / n_increasers;
      for (std::uint64_t t = 1, p = step; p < samples.size(); ++t, p += step)
        bucket_lbound[t] = (samples[p - 1] + samples[p] + 1) / 2;
      bucket_lbound[0] = 0;

      // Compute bucket sizes and sblock id into oracle array.
      std::fill(b->sblock_size, b->sblock_size + n_increasers, 0L);
      for (std::uint64_t t = 0; t < b->m_filled; ++t) {
        const std::uint64_t x = temp[t];

        std::uint64_t id = n_increasers;
        while ((std::uint64_t)bucket_lbound[id] > x)
          --id;

        oracle[t] = id;
        b->sblock_size[id]++;
      }

      // Permute elements into their own buckets using oracle.
      for (std::uint64_t t = 0, curbeg = 0; t < n_increasers;
          curbeg += b->sblock_size[t++])
        b->sblock_beg[t] = ptr[t] = curbeg;

      for (std::uint64_t t = 0; t < b->m_filled; ++t) {
        const std::uint64_t sblock_id = oracle[t];
        oracle[t] = ptr[sblock_id]++;
      }

      for (std::uint64_t t = 0; t < b->m_filled; ++t) {
        const std::uint64_t addr = oracle[t];
        b->m_content[addr] = temp[t];
      }
    }

    // Add the buffer to the poll of full
    // buffers and notify waiting thread.
    std::unique_lock<std::mutex>
      lk2(full_gap_buffers->m_mutex);
    full_gap_buffers->add(b);
    lk2.unlock();
    full_gap_buffers->m_cv.notify_one();
  }

  delete text_streamer;
  
  // Report that another worker thread has finished.
  std::unique_lock<std::mutex> lk(full_gap_buffers->m_mutex);
  full_gap_buffers->increment_finished_workers();
  lk.unlock();

  // Notify waiting update threads in case no
  // more buffers are going to be produces by
  // worker threads.
  full_gap_buffers->m_cv.notify_one();

  delete[] block_count;
  delete[] block_id_to_sblock_id;
  delete[] temp;
  delete[] oracle;
  delete[] ptr;
  delete[] bucket_lbound;
}

}  // namespace psascan_private

#endif  // __SRC_PSASCAN_SRC_STREAM_HPP_INCLUDED
