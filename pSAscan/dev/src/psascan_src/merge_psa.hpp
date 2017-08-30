/*
 * @file    src/psascan_src/merge_psa.hpp
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

#ifndef __SRC_PSASCAN_SRC_MERGE_PSA_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_MERGE_PSA_HPP_INCLUDED

#include <cstdio>
#include <cstdint>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

#include "utils.hpp"
#include "half_block_info.hpp"
#include "io/scatterfile.hpp"
#include "io/async_stream_writer.hpp"
#include "io/async_stream_vbyte_reader.hpp"
#include "io/async_scatterfile_reader.hpp"


namespace psascan_private {

// Merge partial suffix arrays into final suffix array.
template<typename text_offset_type, typename block_offset_type>
void merge(std::string output_filename, std::uint64_t ram_use,
    std::vector<half_block_info<block_offset_type> > &hblock_info) {
  std::uint64_t n_block = hblock_info.size();
  std::uint64_t text_length = 0;

  std::sort(hblock_info.begin(), hblock_info.end());
  for (std::uint64_t j = 0; j < hblock_info.size(); ++j)
    text_length += hblock_info[j].end - hblock_info[j].beg;

  std::uint64_t pieces = (1 + sizeof(block_offset_type)) * n_block -
    1 + sizeof(text_offset_type);
  std::uint64_t buffer_size = (ram_use + pieces - 1) / pieces;

  fprintf(stderr, "\nMerge partial suffix arrays:\n");
  fprintf(stderr, "  buffer size per block = %lu (%.2LfMiB)\n",
      sizeof(block_offset_type) * buffer_size,
      (1.L * sizeof(block_offset_type) * buffer_size) / (1 << 20));

  typedef async_stream_vbyte_reader vbyte_reader_type;
  typedef async_stream_writer<text_offset_type> output_writer_type;
  typedef async_scatterfile_reader<block_offset_type> psa_reader_type;

  psa_reader_type **psa_readers = new psa_reader_type*[n_block];
  output_writer_type *output = new output_writer_type(output_filename,
      sizeof(text_offset_type) * buffer_size, 4UL, "w");
  vbyte_reader_type **gap = new vbyte_reader_type*[n_block - 1];

  for (std::uint64_t i = 0; i < n_block; ++i) {
    psa_readers[i] = new psa_reader_type(&hblock_info[i].psa,
        buffer_size * sizeof(block_offset_type));

    if (i + 1 != n_block)
      gap[i] = new vbyte_reader_type(
          hblock_info[i].gap_filename, buffer_size,
          std::max(4UL, buffer_size / (2UL << 20)));
  }

  std::uint64_t *gap_head = new std::uint64_t[n_block];
  for (std::uint64_t i = 0; i + 1 < n_block; ++i)
    gap_head[i] = gap[i]->read();
  gap_head[n_block - 1] = 0;

  std::uint64_t tmp = sqrtl((long double)n_block);
  std::uint64_t sblock_size = 1L;
  std::uint64_t sblock_size_log = 0;
  while (sblock_size * 2L <= tmp) {
    sblock_size *= 2L;
    ++sblock_size_log;
  }

  std::uint64_t n_sblocks = (n_block + sblock_size - 1) / sblock_size;
  std::pair<std::uint64_t, std::uint64_t> *sblock_info =
    new std::pair<std::uint64_t, std::uint64_t>[n_sblocks];

  for (std::uint64_t i = 0; i < n_sblocks; ++i) {
    std::uint64_t sblock_beg = i * sblock_size;
    std::uint64_t sblock_end = std::min(n_block, sblock_beg + sblock_size);

    sblock_info[i].second = 0;
    sblock_info[i].first = gap_head[sblock_beg];
    for (std::uint64_t j = sblock_beg + 1; j < sblock_end; ++j)
      sblock_info[i].first = std::min(sblock_info[i].first, gap_head[j]);
  }

  long double merge_start = utils::wclock();
  for (std::uint64_t i = 0, dbg = 0; i < text_length; ++i, ++dbg) {
    if (dbg == (1 << 23)) {
      long double elapsed = utils::wclock() - merge_start;
      std::uint64_t inp_vol = (1L + sizeof(block_offset_type)) * i;
      std::uint64_t out_vol = sizeof(text_offset_type) * i;
      std::uint64_t tot_vol = inp_vol + out_vol;
      long double tot_vol_m = tot_vol / (1024.L * 1024);
      long double io_speed = tot_vol_m / elapsed;
      fprintf(stderr, "\r  %.1Lf%%. Time = %.2Lfs. I/O: %2.LfMiB/s",
          (100.L * i) / text_length, elapsed, io_speed);
      dbg = 0;
    }

    // Find the superblock containing gap head equal to zero.
    std::uint64_t k = 0;
    while (sblock_info[k].first != 0) {
      sblock_info[k].first--;
      sblock_info[k].second++;
      ++k;
    }

    // Find the block with the gap head equal to zero.
    std::uint64_t sblock_beg = (k << sblock_size_log);
    std::uint64_t sblock_end = std::min(n_block, sblock_beg + sblock_size);

    std::uint64_t new_min = text_length;
    std::uint64_t j = sblock_beg;
    while (gap_head[j] != sblock_info[k].second) {
      gap_head[j] -= (sblock_info[k].second + 1);
      new_min = std::min(new_min, gap_head[j]);
      ++j;
    }

    std::uint64_t SA_i =
      (std::uint64_t)psa_readers[j]->read() +
      (std::uint64_t)hblock_info[j].beg;

    if (j != n_block - 1) gap_head[j] = gap[j]->read();
    new_min = std::min(new_min, gap_head[j]);
    ++j;

    while (j < sblock_end) {
      gap_head[j] -= sblock_info[k].second;
      new_min = std::min(new_min, gap_head[j]);
      ++j;
    }

    sblock_info[k].first = new_min;
    sblock_info[k].second = 0;

    output->write(SA_i);
  }

  // Stop I/O threads.
  for (std::uint64_t i = 0; i + 1 < n_block; ++i)
    gap[i]->stop_reading();

  // Update I/O volume.
  std::uint64_t io_volume =
    (1 + sizeof(block_offset_type) + sizeof(text_offset_type)) * text_length;

  // Print summary.
  long double merge_time = utils::wclock() - merge_start;
  fprintf(stderr, "\r  100.0%%. Time: %.2Lfs. I/O: %.2LfMiB/s\n",
      merge_time, (io_volume / (1024.L * 1024)) / merge_time);

  // Clean up.
  delete output;
  for (std::uint64_t i = 0; i < n_block; ++i) {
    delete psa_readers[i];
    if (i + 1 != n_block)
      delete gap[i];
  }

  delete[] psa_readers;
  delete[] gap;
  delete[] gap_head;
  delete[] sblock_info;
  
  for (std::uint64_t i = 0; i + 1 < n_block; ++i)
    utils::file_delete(hblock_info[i].gap_filename);
}

}  // namespace psascan_private

#endif  // __SRC_PSASCAN_SRC_MERGE_PSA_HPP_INCLUDED
