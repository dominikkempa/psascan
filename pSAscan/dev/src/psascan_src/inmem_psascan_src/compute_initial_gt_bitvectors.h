/**
 * @file    src/psascan_src/inmem_psascan_src/compute_initial_gt_bitvectors.h
 * @author  Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * @section DESCRIPTION
 *
 * Parallel computation of gt_end bitvectors. The procedure uses the
 * string range matching algorithm described in
 *
 *   Juha Karkkainen, Dominik Kempa, Simon J. Puglisi:
 *   String Range Matching.
 *   In Proc. CPM 2014, p. 232-241.
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

#ifndef __PSASCAN_SRC_INMEM_PSASCAN_SRC_COMPUTE_INITIAL_GT_BITVECTORS_H_INCLUDED
#define __PSASCAN_SRC_INMEM_PSASCAN_SRC_COMPUTE_INITIAL_GT_BITVECTORS_H_INCLUDED

#include <cstdio>
#include <cstring>
#include <algorithm>
#include <thread>

#include "../bitvector.h"
#include "../multifile.h"
#include "../multifile_bit_stream_reader.h"
#include "../background_block_reader.h"
#include "srank_aux.h"


namespace psascan_private {
namespace inmem_psascan_private {

//==============================================================================
// Compute bitvectors bv[0..ref_pos) and undecided[0..ref_pos), where:
//   undecided[i] == 0 iff lcp(text[0..txtlen), text[ref_pos..txtlen)) < max_lcp
//   gt[i] == 1 iff undecided[i] == 0 and text[0..txtlen) > text[ref_pos..txtlen)
//
// NOTE: we store the gt bitvector reversed, so that later we can overwrite
// it with gt_begin in place.
//==============================================================================
void compute_partial_gt_end(const unsigned char *text, std::uint64_t text_length,
    std::uint64_t begin, std::uint64_t end, std::uint64_t max_lcp,
    bitvector *gt, bitvector *undecided,
    bool &all_decided, std::uint64_t text_end, std::uint64_t supertext_length,
    const multifile *tail_gt_begin_rev,
    background_block_reader *tail_prefix_background_reader,
    const std::uint8_t *tail_prefix_preread) {
  bool res = true;
  all_decided = true;
  std::uint64_t revbeg = text_length - end;

  if (end == text_length) {
    // It's ok if tail_gt_begin_rev is NULL
    multifile_bit_stream_reader tail_gt_beg_rev(tail_gt_begin_rev);
    std::uint64_t tail_length = supertext_length - text_end;
    std::uint64_t range_size = end - begin;
    std::uint64_t tail_prefix_length = std::min(text_length, tail_length);
    std::uint64_t tail_prefix_fetched = 0;

    const std::uint8_t *txt = text + begin;
    const std::uint8_t *tail_prefix = NULL;

    if (tail_prefix_length > 0) {
      if (tail_prefix_preread != NULL) {
        // Whole tail prefix is already in memory.
        tail_prefix = tail_prefix_preread;
        tail_prefix_fetched = tail_prefix_length;
      } else {
        // Tail prefix will be fetched asynchronously in the background.
        tail_prefix = tail_prefix_background_reader->m_data;
        tail_prefix_fetched = 0;
      }
    }

    std::uint64_t i = 0, el = 0, s = 0, p = 0;
    std::uint64_t i_max = 0, el_max = 0, s_max = 0, p_max = 0;

    static const std::uint64_t chunk_size = (1UL << 20);

    while (i < range_size) {
      while (i + el < range_size && el < tail_length) {
        if (el == tail_prefix_fetched) {
          std::uint64_t next_chunk = std::min(chunk_size,
              tail_prefix_length - tail_prefix_fetched);
          tail_prefix_fetched += next_chunk;
          tail_prefix_background_reader->wait(tail_prefix_fetched);
        }
        while (i + el < range_size && el < tail_length &&
            el < tail_prefix_fetched && txt[i + el] == tail_prefix[el])
          update_ms(tail_prefix, ++el, s, p);
        if (el < tail_prefix_fetched) break;
      }

      if ((el == tail_length) ||
          (i + el == range_size && !tail_gt_beg_rev.access(tail_length - el)) ||
          (i + el < range_size && txt[i + el] > tail_prefix[el]))
        gt->set(revbeg + i);

      std::uint64_t j = i_max;
      if (el > el_max) {
        std::swap(el, el_max);
        std::swap(s, s_max);
        std::swap(p, p_max);
        i_max = i;
      }

      if (el < 100) {
        ++i;
        el = 0;
      } else if (p > 0 && (p << 2) <= el &&
          !memcmp(tail_prefix, tail_prefix + p, s)) {
        std::uint64_t maxk = std::min(p, range_size - i);
        for (std::uint64_t k = 1; k < maxk; ++k)
          if (gt->get(revbeg + j + k)) gt->set(revbeg + i + k);
        i += p;
        el -= p;
      } else {
        std::uint64_t h = (el >> 2) + 1;
        std::uint64_t maxk = std::min(h, range_size - i);
        for (std::uint64_t k = 1; k < maxk; ++k)
          if (gt->get(revbeg + j + k)) gt->set(revbeg + i + k);
        i += h;
        el = 0;
        p = 0;
        s = 0;
      }
    }
  } else {
    std::uint64_t i = 0, el = 0, s = 0, p = 0;
    std::uint64_t i_max = 0, el_max = 0, s_max = 0, p_max = 0;

    const std::uint8_t *txt = text + begin;
    const std::uint8_t *pat = text + end;
    std::uint64_t range_size = end - begin;

    while (i < range_size) {
      while (el < max_lcp && txt[i + el] == pat[el])
        update_ms(pat, ++el, s, p);

      if (el < max_lcp) {
        if (txt[i + el] > pat[el]) gt->set(revbeg + i);
      } else {
        undecided->set(revbeg + i);
        res = false;
      }

      std::uint64_t j = i_max;
      if (el > el_max) {
        std::swap(el, el_max);
        std::swap(s, s_max);
        std::swap(p, p_max);
        i_max = i;
      }

      if (el < 100) {
        ++i;
        el = 0;
      } else if (p > 0 && (p << 2) <= el && !memcmp(pat, pat + p, s)) {
        std::uint64_t maxk = std::min(p, range_size - i);
        for (std::uint64_t k = 1; k < maxk; ++k) {
          if (undecided->get(revbeg + (j + k))) undecided->set(revbeg + (i + k));
          if (gt->get(revbeg + (j + k))) gt->set(revbeg + (i + k));
        }

        i += p;
        el -= p;
      } else {
        std::uint64_t h = (el >> 2) + 1;
        std::uint64_t maxk = std::min(h, range_size - i);
        for (std::uint64_t k = 1; k < maxk; ++k) {
          if (undecided->get(revbeg + (j + k))) undecided->set(revbeg + (i + k));
          if (gt->get(revbeg + (j + k))) gt->set(revbeg + (i + k));
        }

        i += h;
        el = 0;
        s = 0;
        p = 0;
      }
    }
  }

  all_decided = res;
}

//==============================================================================
// Set all undecided bits inside the given microblock (that is, the range
// [mb_beg..mb_end)) of all gt bitvectors to their correct values.
//==============================================================================
void compute_final_gt(std::uint64_t text_length, std::uint64_t max_block_size,
    std::uint64_t mb_beg, std::uint64_t mb_end, bitvector *gt,
    const bitvector *undecided, const bool *all_decided) {
  std::uint64_t n_blocks = (text_length + max_block_size - 1) / max_block_size;

  // Go through blocks right to left.
  for (std::uint64_t t = n_blocks - 1; t > 0; --t) {
    std::uint64_t block_end = text_length - (n_blocks - t) * max_block_size;
    std::uint64_t block_beg = (std::uint64_t)std::max(0L, (std::int64_t)block_end - (std::int64_t)max_block_size);
    std::uint64_t this_block_size = block_end - block_beg;
    std::uint64_t this_mb_beg = mb_beg;
    std::uint64_t this_mb_end = std::min(mb_end, this_block_size);

    std::uint64_t rev_beg = text_length - block_end;
    std::uint64_t rev_end = text_length - block_beg;

    if (!all_decided[t - 1]) {
      // This eliminates the problem with accessing bits located in the same
      // byte in the bitvector. Skipped bits are later updated sequentially.
      while (((rev_end - 1 - this_mb_beg) & 7) != 7) ++this_mb_beg;
      for (std::uint64_t j = this_mb_beg; j < this_mb_end; ++j)
        if (undecided->get(rev_end - 1 - j) && gt->get(rev_beg - 1 - j))
          gt->set(rev_end - 1 - j);
    }
  }
}

//==============================================================================
// Update the bits omitted in compute_final_gt.
//==============================================================================
void compute_final_gt_last_bits(std::uint64_t text_length, std::uint64_t max_block_size,
    std::uint64_t mb_beg, std::uint64_t mb_end, bitvector *gt, const bitvector *undecided,
    bool *all_decided) {
  std::uint64_t n_blocks = (text_length + max_block_size - 1) / max_block_size;
  if (!all_decided[0]) {
    std::uint64_t block_end = text_length - (n_blocks - 1) * max_block_size;
    std::uint64_t block_beg = (std::uint64_t)std::max(0L, (std::int64_t)block_end - (std::int64_t)max_block_size);
    std::uint64_t this_block_size = block_end - block_beg;
    std::uint64_t this_mb_beg = mb_beg;
    std::uint64_t this_mb_end = std::min(mb_end, this_block_size);

    std::uint64_t rev_beg = text_length - block_end;
    std::uint64_t rev_end = text_length - block_beg;

    std::uint64_t temp_this_mb_beg = this_mb_beg;
    while (((rev_end - 1 - temp_this_mb_beg) & 7) != 7) ++temp_this_mb_beg;
    this_mb_end = temp_this_mb_beg;

    // [this_mb_beg..this_mb_end) were omitted.
    for (std::uint64_t j = this_mb_beg; j < this_mb_end; ++j)
      if (undecided->get(rev_end - 1 - j) && gt->get(rev_beg - 1 - j))
        gt->set(rev_end - 1 - j);
  }
}

//==============================================================================
// Fully parallel computation of gt bitvectors.
//==============================================================================
void compute_initial_gt_bitvectors(const std::uint8_t *text, std::uint64_t text_length,
    bitvector *gt, std::uint64_t max_block_size, std::uint64_t max_threads, std::uint64_t text_end,
    std::uint64_t supertext_length, const multifile *tail_gt_begin_reversed,
    background_block_reader *tail_prefix_background_reader,
    const std::uint8_t *tail_prefix_preread) {
  long double start;
  std::uint64_t n_blocks = (text_length + max_block_size - 1) / max_block_size;

  //----------------------------------------------------------------------------
  // STEP 1: compute gt bitvectors, some bits may still be undecided after this.
  //----------------------------------------------------------------------------

  // Allocate ane zero-initialize (in parallel) bitvectors.
  fprintf(stderr, "  Allocate: ");
  start = utils::wclock();
  bitvector *undecided = new bitvector(text_length);
  fprintf(stderr, "%.2Lfs\n", utils::wclock() - start);

  // all_decided[i] == true, if all bits inside block i were
  // decided in the first stage. This can be used by threads in the
  // second stage to completely skip inspecting some blocks.
  bool *all_decided = new bool[n_blocks];

  // Process blocks right-to-left.
  fprintf(stderr, "  Compute decided bits: ");
  start = utils::wclock();
  std::thread **threads = new std::thread*[n_blocks];
  for (std::uint64_t i = 0; i < n_blocks; ++i) {
    std::uint64_t block_end = text_length - (n_blocks - 1 - i) * max_block_size;
    std::uint64_t block_beg = (std::uint64_t)std::max(0L, (std::int64_t)block_end - (std::int64_t)max_block_size);

    // Compute bitvectors 'gt' and 'undecided' for block i.
    threads[i] = new std::thread(compute_partial_gt_end,
        text, text_length, block_beg, block_end, max_block_size, gt,
        undecided, std::ref(all_decided[i]), text_end, supertext_length,
        tail_gt_begin_reversed, tail_prefix_background_reader,
        tail_prefix_preread);
  }

  // Wait for the threads to finish and clean up.
  for (std::uint64_t i = 0; i < n_blocks; ++i) threads[i]->join();
  for (std::uint64_t i = 0; i < n_blocks; ++i) delete threads[i];
  delete[] threads;
  fprintf(stderr, "%.2Lfs\n", utils::wclock() - start);

  //----------------------------------------------------------------------------
  // STEP 2: compute the undecided bits in the gt bitvectors.
  //----------------------------------------------------------------------------
  
  // The size of micro block has to be a multiple of 8, otherwise two
  // threads might try to update the same char inside bitvector.
  std::uint64_t max_microblock_size = (max_block_size + max_threads - 1) / max_threads;
  while ((max_microblock_size & 7) && max_microblock_size < max_block_size)
    ++max_microblock_size;
  std::uint64_t n_microblocks = (max_block_size + max_microblock_size - 1) / max_microblock_size;

  fprintf(stderr, "  Compute undecided bits: ");
  start = utils::wclock();
  threads = new std::thread*[n_microblocks];
  for (std::uint64_t i = 0; i < n_microblocks; ++i) {
    std::uint64_t mb_beg = i * max_microblock_size;
    std::uint64_t mb_end = std::min(mb_beg + max_microblock_size, max_block_size);

    threads[i] = new std::thread(compute_final_gt, text_length, max_block_size,
        mb_beg, mb_end, std::ref(gt), std::ref(undecided), all_decided);
  }

  // Wait for the threads to finish and clean up.
  for (std::uint64_t i = 0; i < n_microblocks; ++i) threads[i]->join();
  for (std::uint64_t i = 0; i < n_microblocks; ++i) delete threads[i];
  fprintf(stderr, "%.2Lfs\n", utils::wclock() - start);
  
  // Fill in the skipped (due to parallel byte access issue) undecided bits.
  for (std::uint64_t i = 0; i < n_microblocks; ++i) {
    std::uint64_t mb_beg = i * max_microblock_size;
    std::uint64_t mb_end = std::min(mb_beg + max_microblock_size, max_block_size);

    compute_final_gt_last_bits(text_length, max_block_size, mb_beg, mb_end,
        gt, undecided, all_decided);
  }

  fprintf(stderr, "  Deallocate: ");
  start = utils::wclock();
  delete[] threads;
  delete undecided;
  delete[] all_decided;
  fprintf(stderr, "%.2Lfs\n", utils::wclock() - start);
}

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __PSASCAN_SRC_INMEM_PSASCAN_SRC_COMPUTE_INITIAL_GT_BITVECTORS_H_INCLUDED
