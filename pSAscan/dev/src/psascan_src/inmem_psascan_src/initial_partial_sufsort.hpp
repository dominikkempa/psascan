/**
 * @file    src/psascan_src/inmem_psascan_src/initial_partial_sufsort.hpp
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

#ifndef __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_INITIAL_PARTIAL_SUFSORT_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_INITIAL_PARTIAL_SUFSORT_HPP_INCLUDED

#include <algorithm>
#include <thread>

#include "../bitvector.hpp"
#include "divsufsort_template.hpp"
#include "bwtsa.hpp"
#include "parallel_shrink.hpp"
#include "parallel_expand.hpp"


namespace psascan_private {
namespace inmem_psascan_private {

//=============================================================================
// Rename the given block using its gt bitvector.
//=============================================================================
void rename_block(
    std::uint8_t *text,
    std::uint64_t text_length,
    std::uint64_t block_beg,
    std::uint64_t block_length,
    bitvector *gt,
    bool &renaming_error) {

  std::uint64_t block_end = block_beg + block_length;
  std::uint64_t beg_rev = text_length - block_end;
  std::uint8_t *block = text + block_beg;
  std::uint8_t last = block[block_length - 1];
  bool err = false;
  for (std::uint64_t i = 0; i + 1 < block_length; ++i)
    if (block[i] > last ||
        (block[i] == last && gt->get(beg_rev + i + 1))) {
      if (block[i] == 255)
        err = true;
      ++block[i];
    }
  if (block[block_length - 1] == 255)
    err = true;
  ++block[block_length - 1];

  if (err)
    renaming_error = true;
}

//=============================================================================
// Re-rename block back to original.
//=============================================================================
void rerename_block(std::uint8_t *block, std::uint64_t block_length) {
  std::uint8_t last = block[block_length - 1] - 1;
  for (std::uint64_t i = 0; i < block_length; ++i)
    if (block[i] > last) --block[i];
}

//==============================================================================
// Given gt bitvectors, compute partial suffix arrays of blocks.
//==============================================================================
template<typename block_offset_type>
void initial_partial_sufsort(
    std::uint8_t *,
    std::uint64_t,
    bitvector *,
    bwtsa_t<block_offset_type> *,
    std::uint64_t,
    std::uint64_t, bool) {
  fprintf(stderr, "\n\nError: initial_partial_sufsort: "
      "given block_offset_type is not supported, "
      "sizeof(block_offset_type) = %lu\n", sizeof(block_offset_type));
  std::exit(EXIT_FAILURE);
}

template<>
void initial_partial_sufsort(
    std::uint8_t *text,
    std::uint64_t text_length,
    bitvector* gt,
    bwtsa_t<std::int64_t> *bwtsa,
    std::uint64_t max_block_size,
    std::uint64_t max_threads,
    bool has_tail) {

  long double start = utils::wclock();
  std::uint64_t n_blocks =
    (text_length + max_block_size - 1) / max_block_size;

  //----------------------------------------------------------------------------
  // STEP 1: Rename the blocks in parallel.
  //----------------------------------------------------------------------------
  if (n_blocks > 1 || has_tail) {
    fprintf(stderr, "  Rename blocks: ");
    start = utils::wclock();
    bool *renaming_error = new bool[n_blocks];
    std::fill(renaming_error, renaming_error + n_blocks, false);
    std::thread **threads = new std::thread*[n_blocks];
    for (std::uint64_t i = 0; i < n_blocks; ++i) {
      std::uint64_t block_end =
        text_length - (n_blocks - 1 - i) * max_block_size;
      std::uint64_t block_beg =
        std::max(0L, (std::int64_t)block_end - (std::int64_t)max_block_size);
      std::uint64_t block_size = block_end - block_beg;

      threads[i] = new std::thread(rename_block, text, text_length,
          block_beg, block_size, gt, std::ref(renaming_error[i]));
    }

    for (std::uint64_t i = 0; i < n_blocks; ++i) threads[i]->join();
    for (std::uint64_t i = 0; i < n_blocks; ++i) delete threads[i];
    delete[] threads;

    fprintf(stderr, "%.2Lfs\n", utils::wclock() - start);

    bool err = false;
    for (std::uint64_t i = 0; i < n_blocks; ++i)
      if (renaming_error[i]) err = true;
    delete[] renaming_error;

    if (err) {
      fprintf(stdout, "\n\nError: byte with value 255 was "
          "detected in the input text!\n"
          "See the section on limitations in the README "
          "for more information.\n");
      std::fflush(stdout);
      std::exit(EXIT_FAILURE);
    }
  }

  {

    // Use 64-bit divsufsort.
    std::int64_t *temp_sa = (std::int64_t *)bwtsa;

    //--------------------------------------------------------------------------
    // STEP 2: Compute suffix arrays in parallel.
    //--------------------------------------------------------------------------
    fprintf(stderr, "  Run divsufsort32: ");
    start = utils::wclock();
    std::thread **threads = new std::thread*[n_blocks];
    for (std::uint64_t i = 0; i < n_blocks; ++i) {
      std::uint64_t block_end =
        text_length - (n_blocks - 1 - i) * max_block_size;
      std::uint64_t block_beg = std::max((std::int64_t)0,
          (std::int64_t)block_end - (std::int64_t)max_block_size);
      std::uint64_t block_size = block_end - block_beg;

      threads[i] = new std::thread(run_divsufsort<std::int64_t>,
          text + block_beg, temp_sa + block_beg, block_size);
    }

    for (std::uint64_t i = 0; i < n_blocks; ++i) threads[i]->join();
    for (std::uint64_t i = 0; i < n_blocks; ++i) delete threads[i];
    delete[] threads;

    fprintf(stderr, "%.2Lfs\n", utils::wclock() - start);

    fprintf(stderr, "  Expand 64-bit integers to bwtsa objects: ");
    start = utils::wclock();
    parallel_expand<std::int64_t, bwtsa_t<int64_t> >(temp_sa,
        text_length, max_threads);
    fprintf(stderr, "%.2Lfs\n", utils::wclock() - start);
  }

  //----------------------------------------------------------------------------
  // STEP 3: Restore the original text.
  //----------------------------------------------------------------------------
  if (n_blocks > 1 || has_tail) {
    fprintf(stderr, "  Rerename blocks: ");
    start = utils::wclock();
    std::thread **threads = new std::thread*[n_blocks];
    for (std::uint64_t i = 0; i < n_blocks; ++i) {
      std::uint64_t block_end =
        text_length - (n_blocks - 1 - i) * max_block_size;
      std::uint64_t block_beg = std::max((std::int64_t)0,
          (std::int64_t)block_end - (std::int64_t)max_block_size);
      std::uint64_t block_size = block_end - block_beg;

      threads[i] = new std::thread(rerename_block,
          text + block_beg, block_size);
    }

    for (std::uint64_t i = 0; i < n_blocks; ++i) threads[i]->join();
    for (std::uint64_t i = 0; i < n_blocks; ++i) delete threads[i];
    delete[] threads;

    fprintf(stderr, "%.2Lfs\n", utils::wclock() - start);
  }
}


template<>
void initial_partial_sufsort(
    std::uint8_t *text,
    std::uint64_t text_length,
    bitvector* gt,
    bwtsa_t<uint40> *bwtsa,
    std::uint64_t max_block_size,
    std::uint64_t max_threads,
    bool has_tail) {

  long double start = utils::wclock();
  std::uint64_t n_blocks =
    (text_length + max_block_size - 1) / max_block_size;

  //----------------------------------------------------------------------------
  // STEP 1: Rename the blocks in parallel.
  //----------------------------------------------------------------------------

  // XXX change this parallelism to vertical!
  if (n_blocks > 1 || has_tail) {
    fprintf(stderr, "  Rename blocks: ");
    start = utils::wclock();
    bool *renaming_error = new bool[n_blocks];
    std::fill(renaming_error, renaming_error + n_blocks, false);
    std::thread **threads = new std::thread*[n_blocks];
    for (std::uint64_t i = 0; i < n_blocks; ++i) {
      std::uint64_t block_end =
        text_length - (n_blocks - 1 - i) * max_block_size;
      std::uint64_t block_beg = std::max((std::int64_t)0,
          (std::int64_t)block_end - (std::int64_t)max_block_size);
      std::uint64_t block_size = block_end - block_beg;

      threads[i] = new std::thread(rename_block, text, text_length,
          block_beg, block_size, gt, std::ref(renaming_error[i]));
    }

    for (std::uint64_t i = 0; i < n_blocks; ++i) threads[i]->join();
    for (std::uint64_t i = 0; i < n_blocks; ++i) delete threads[i];
    delete[] threads;

    fprintf(stderr, "%.2Lfs\n", utils::wclock() - start);

    bool err = false;
    for (std::uint64_t i = 0; i < n_blocks; ++i)
      if (renaming_error[i]) err = true;
    delete[] renaming_error;

    if (err) {
      fprintf(stdout, "\n\nError: byte with value 255 "
          "was detected in the input text!\n"
          "See the section on limitations in the README "
          "for more information.\n");
      std::fflush(stdout);
      std::exit(EXIT_FAILURE);
    }
  }

  if (max_block_size >= (2UL << 30)) {  // Use 64-bit divsufsort.
    fprintf(stdout, "\n\nError: 2GiB+ partial suffix arrays are not "
        "yet supported by the internal-memory pSAscan.\n");
    std::fflush(stdout);
    std::exit(EXIT_FAILURE);
  } else {  // Use 32-bit divsufsort.
    std::int32_t *temp_sa = (std::int32_t *)bwtsa;

    //--------------------------------------------------------------------------
    // STEP 2: Compute suffix arrays in parallel.
    //--------------------------------------------------------------------------
    fprintf(stderr, "  Run divsufsort32: ");
    start = utils::wclock();
    std::thread **threads = new std::thread*[n_blocks];
    for (std::uint64_t i = 0; i < n_blocks; ++i) {
      std::uint64_t block_end =
        text_length - (n_blocks - 1 - i) * max_block_size;
      std::uint64_t block_beg = std::max((std::int64_t)0,
          (std::int64_t)block_end - (std::int64_t)max_block_size);
      std::uint64_t block_size = block_end - block_beg;

      threads[i] = new std::thread(run_divsufsort<std::int32_t>,
          text + block_beg, temp_sa + block_beg, block_size);
    }

    for (std::uint64_t i = 0; i < n_blocks; ++i) threads[i]->join();
    for (std::uint64_t i = 0; i < n_blocks; ++i) delete threads[i];
    delete[] threads;

    fprintf(stderr, "%.2Lfs\n", utils::wclock() - start);

    fprintf(stderr, "  Expand 32-bit integers to bwtsa objects: ");
    start = utils::wclock();
    parallel_expand<std::int32_t, bwtsa_t<uint40> >(temp_sa,
        text_length, max_threads);
    fprintf(stderr, "%.2Lfs\n", utils::wclock() - start);
  }

  //----------------------------------------------------------------------------
  // STEP 3: Restore the original text.
  //----------------------------------------------------------------------------
  if (n_blocks > 1 || has_tail) {
    fprintf(stderr, "  Rerename blocks: ");
    start = utils::wclock();
    std::thread **threads = new std::thread*[n_blocks];
    for (std::uint64_t i = 0; i < n_blocks; ++i) {
      std::uint64_t block_end =
        text_length - (n_blocks - 1 - i) * max_block_size;
      std::uint64_t block_beg = std::max((std::int64_t)0,
          (std::int64_t)block_end - (std::int64_t)max_block_size);
      std::uint64_t block_size = block_end - block_beg;

      threads[i] = new std::thread(rerename_block,
          text + block_beg, block_size);
    }

    for (std::uint64_t i = 0; i < n_blocks; ++i) threads[i]->join();
    for (std::uint64_t i = 0; i < n_blocks; ++i) delete threads[i];
    delete[] threads;

    fprintf(stderr, "%.2Lfs\n", utils::wclock() - start);
  }
}

template<>
void initial_partial_sufsort(
    std::uint8_t *text,
    std::uint64_t text_length,
    bitvector* gt,
    bwtsa_t<int> *bwtsa,
    std::uint64_t max_block_size,
    std::uint64_t max_threads,
    bool has_tail) {

  long double start = utils::wclock();
  std::uint64_t n_blocks =
    (text_length + max_block_size - 1) / max_block_size;

  //----------------------------------------------------------------------------
  // STEP 1: Rename the blocks in parallel.
  //----------------------------------------------------------------------------
  // XXX change this parallelism to vertical!
  if (n_blocks > 1 || has_tail) {
    fprintf(stderr, "  Rename blocks: ");
    start = utils::wclock();
    bool *renaming_error = new bool[n_blocks];
    std::fill(renaming_error, renaming_error + n_blocks, false);
    std::thread **threads = new std::thread*[n_blocks];
    for (std::uint64_t i = 0; i < n_blocks; ++i) {
      std::uint64_t block_end =
        text_length - (n_blocks - 1 - i) * max_block_size;
      std::uint64_t block_beg = std::max((std::int64_t)0,
          (std::int64_t)block_end - (std::int64_t)max_block_size);
      std::uint64_t block_size = block_end - block_beg;

      threads[i] = new std::thread(rename_block, text,
          text_length, block_beg,
          block_size, gt, std::ref(renaming_error[i]));
    }

    for (std::uint64_t i = 0; i < n_blocks; ++i) threads[i]->join();
    for (std::uint64_t i = 0; i < n_blocks; ++i) delete threads[i];
    delete[] threads;

    fprintf(stderr, "%.2Lfs\n", utils::wclock() - start);

    bool err = false;
    for (std::uint64_t i = 0; i < n_blocks; ++i)
      if (renaming_error[i]) err = true;
    delete[] renaming_error;

    if (err) {
      fprintf(stdout, "\n\nError: byte with value 255 "
          "was detected in the input text!\n"
          "See the section on limitations in the README "
          "for more information.\n");
      std::fflush(stdout);
      std::exit(EXIT_FAILURE);
    }
  }
  
  int *temp_sa = (int *)bwtsa;

  //----------------------------------------------------------------------------
  // STEP 2: Compute suffix arrays in parallel.
  //----------------------------------------------------------------------------
  fprintf(stderr, "  Run divsufsort32: ");
  start = utils::wclock();
  std::thread **threads = new std::thread*[n_blocks];
  for (std::uint64_t i = 0; i < n_blocks; ++i) {
    std::uint64_t block_end =
      text_length - (n_blocks - 1 - i) * max_block_size;
    std::uint64_t block_beg = std::max((std::int64_t)0,
        (std::int64_t)block_end - (std::int64_t)max_block_size);
    std::uint64_t block_size = block_end - block_beg;

    threads[i] = new std::thread(run_divsufsort<int>,
        text + block_beg, temp_sa + block_beg, block_size);
  }

  for (std::uint64_t i = 0; i < n_blocks; ++i) threads[i]->join();
  for (std::uint64_t i = 0; i < n_blocks; ++i) delete threads[i];
  delete[] threads;
  fprintf(stderr, "%.2Lfs\n", utils::wclock() - start);

  fprintf(stderr, "  Expand 32-bit integers to bwtsa objects: ");
  start = utils::wclock();
  parallel_expand<int, bwtsa_t<int> >(temp_sa, text_length, max_threads);
  fprintf(stderr, "%.2Lfs\n", utils::wclock() - start);

  //----------------------------------------------------------------------------
  // STEP 3: Restore the original text.
  //----------------------------------------------------------------------------
  // XXX: change parallelism to vertical.
  if (n_blocks > 1 || has_tail) {
    fprintf(stderr, "  Rerename blocks: ");
    start = utils::wclock();
    threads = new std::thread*[n_blocks];
    for (std::uint64_t i = 0; i < n_blocks; ++i) {
      std::uint64_t block_end =
        text_length - (n_blocks - 1 - i) * max_block_size;
      std::uint64_t block_beg = std::max((std::int64_t)0,
          (std::int64_t)block_end - (std::int64_t)max_block_size);
      std::uint64_t block_size = block_end - block_beg;

      threads[i] = new std::thread(rerename_block,
          text + block_beg, block_size);
    }

    for (std::uint64_t i = 0; i < n_blocks; ++i) threads[i]->join();
    for (std::uint64_t i = 0; i < n_blocks; ++i) delete threads[i];
    delete[] threads;

    fprintf(stderr, "%.2Lfs\n", utils::wclock() - start);
  }
}

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_INITIAL_PARTIAL_SUFSORT_HPP_INCLUDED
