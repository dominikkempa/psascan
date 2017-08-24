/**
 * @file    src/psascan_src/inmem_psascan_src/inmem_compute_initial_ranks.hpp
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

#ifndef __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_COMPUTE_INITIAL_RANKS_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_COMPUTE_INITIAL_RANKS_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>
#include <thread>
#include <mutex>

#include "../io/background_block_reader.hpp"
#include "../io/multifile_bit_stream_reader.hpp"
#include "../io/multifile.hpp"
#include "../utils.hpp"
#include "bwtsa.hpp"
#include "pagearray.hpp"


namespace psascan_private {
namespace inmem_psascan_private {

inline int lcp_compare(
    const std::uint8_t *text,
    std::uint64_t text_length,
    const std::uint8_t *pat,
    std::uint64_t pat_length,
    std::uint64_t gt_begin_length,
    std::uint64_t j,
    multifile_bit_stream_reader &rev_gt_begin_reader,
    std::uint64_t &lcp) {

  while (lcp < pat_length &&
      j + lcp < text_length &&
      pat[lcp] == text[j + lcp])
    ++lcp;

  if (j + lcp >= text_length) {
    if (rev_gt_begin_reader.access(gt_begin_length - (text_length - j)))
      return 1;
    else return -1;
  } else if (lcp == pat_length) return 0;
  else {
    if (pat[lcp] < text[j + lcp]) return -1;
    else return 1;
  }
}

inline int lcp_compare(
    const std::uint8_t *text,
    const std::uint8_t *pat,
    std::uint64_t pat_length,
    std::uint64_t j,
    std::uint64_t &lcp) {

  while (lcp < pat_length && pat[lcp] == text[j + lcp])
    ++lcp;

  if (lcp == pat_length) return 0;
  else if (pat[lcp] < text[j + lcp]) return -1;
  else return 1;
}

//-----------------------------------------------------------------------------
// Find the range [left..right) of suffixes starting inside the block that are
// prefixed with pat[0..pat_length). In case there is no such suffix, left ==
// right and they both point to the first suffix larger than the pattern.
//-----------------------------------------------------------------------------
template<typename pagearray_type>
void compute_range(
    const std::uint8_t *text,
    std::uint64_t block_beg,
    std::uint64_t block_size,
    const std::uint8_t *pat,
    std::uint64_t pat_length,
    const pagearray_type &bwtsa,
    std::pair<std::uint64_t, std::uint64_t> &ret) {

#ifdef BLOCK_MATRIX_MODULE_DEBUG_MODE
  std::uint64_t min_discrepancy =
    utils::random_int64((std::int64_t)0, (std::int64_t)10);
  std::uint64_t balancing_factor =
    utils::random_int64((std::int64_t)1, (std::int64_t)10);
#else
  static const std::uint64_t min_discrepancy = (1 << 16);
  static const std::uint64_t balancing_factor = 64;
#endif

  // Find left.
  std::int64_t low = -1;
  std::int64_t high = block_size;
  std::uint64_t llcp = 0;
  std::uint64_t rlcp = 0;
  while (low + 1 != high) {

    // Invariant: left is in the range (low..high].
    std::uint64_t lcp = std::min(llcp, rlcp);

    // Compute mid.
    // Valid values for mid are: low + 1, .., high - 1.
    std::uint64_t mid = 0;
    if (llcp + min_discrepancy < rlcp) {

      // Choose the pivot that split the range into two
      // parts of sizes with ratio equal to logd / d.
      std::uint64_t d = rlcp - llcp;
      std::uint64_t logd = utils::log2ceil(d);
      mid = low + 1 +
        ((high - low - 1) * balancing_factor * logd) /
        (d + balancing_factor * logd);
    } else if (rlcp + min_discrepancy < llcp) {
      std::uint64_t d = llcp - rlcp;
      std::uint64_t logd = utils::log2ceil(d);
      mid = high - 1 -
        ((high - low - 1) * balancing_factor * logd) /
        (d + balancing_factor * logd);
    } else {

      // Discrepancy is too small, use standard binary search.
      mid = (low + high) / 2;
    }

    if (lcp_compare(text, pat, pat_length,
          block_beg + (std::uint64_t)bwtsa[mid].m_sa, lcp) <= 0) {
      high = mid;
      rlcp = lcp;
    } else {
      low = mid;
      llcp = lcp;
    }
  }
  std::uint64_t left = high;

  // Find right.
  if (rlcp == pat_length) {
    high = block_size;
    rlcp = 0;

    while (low + 1 != high) {

      // Invariant: right is in the range (low..high].
      std::uint64_t lcp = std::min(llcp, rlcp);
      std::uint64_t mid = 0;
      if (llcp + min_discrepancy < rlcp) {
        std::uint64_t d = rlcp - llcp;
        std::uint64_t logd = utils::log2ceil(d);
        mid = low + 1 +
          ((high - low - 1) * balancing_factor * logd) /
          (d + balancing_factor * logd);
      } else if (rlcp + min_discrepancy < llcp) {
        std::uint64_t d = llcp - rlcp;
        std::uint64_t logd = utils::log2ceil(d);
        mid = high - 1 -
          ((high - low - 1) * balancing_factor * logd) /
          (d + balancing_factor * logd);
      } else mid = (low + high) / 2;

      if (lcp_compare(text, pat, pat_length,
            block_beg + (std::uint64_t)bwtsa[mid].m_sa, lcp) < 0) {
        high = mid;
        rlcp = lcp;
      } else {
        low = mid;
        llcp = lcp;
      }
    }
  }
  std::uint64_t right = high;

  ret = std::make_pair(left, right);
}

//------------------------------------------------------------------------------
// On the entry to the function:
// - all suffixes in the range [0..left) are smaller than pat[0..old_pat_length),
// - all suffixes in the range [right..text_length) are larger than the pattern,
// - suffixes in the range [left..right) are unknown -- they can either be
//   larger or smaller than the pattern, or equal -- in any case, they have a
//   common prefix of length `old_pat_length' with the pattern.
//------------------------------------------------------------------------------
template<typename block_offset_type>
void refine_range(
    const std::uint8_t *text,
    std::uint64_t block_beg,
    const bwtsa_t<block_offset_type> *block_psa,
    std::uint64_t left,
    std::uint64_t right,
    std::uint64_t old_pat_length,
    std::uint64_t pat_length,
    const std::uint8_t *pat,
    std::uint64_t &newleft,
    std::uint64_t &newright) {

  std::int64_t low = left - 1;
  std::int64_t high = right;
  std::uint64_t llcp = old_pat_length;
  std::uint64_t rlcp = old_pat_length;

#ifdef BLOCK_MATRIX_MODULE_DEBUG_MODE
  std::uint64_t min_discrepancy =
    utils::random_int64((std::int64_t)0, (std::int64_t)10);
  std::uint64_t balancing_factor =
    utils::random_int64((std::int64_t)1, (std::int64_t)10);
#else
  static const std::uint64_t min_discrepancy = (1 << 16);
  static const std::uint64_t balancing_factor = 64;
#endif

  while (low + 1 != high) {

    // Invariant: newleft is in the range (low, high].
    std::uint64_t lcp = std::min(llcp, rlcp);

    // Compute mid.
    // Valid values for mid are: low + 1, .., high - 1.
    std::uint64_t mid = 0;
    if (llcp + min_discrepancy < rlcp) {

      // Choose the pivot that split the range into two
      // parts of sizes with ratio equal to logd / d.
      std::uint64_t d = rlcp - llcp;
      std::uint64_t logd = utils::log2ceil(d);
      mid = low + 1 +
        ((high - low - 1) * balancing_factor * logd) /
        (d + balancing_factor * logd);
    } else if (rlcp + min_discrepancy < llcp) {
      std::uint64_t d = llcp - rlcp;
      std::uint64_t logd = utils::log2ceil(d);
      mid = high - 1 -
        ((high - low - 1) * balancing_factor * logd) /
        (d + balancing_factor * logd);
    } else {

      // Discrepancy is too small, use standard binary search.
      mid = (low + high) / 2;
    }

    if (lcp_compare(text, pat, pat_length,
          block_beg + (std::uint64_t)block_psa[mid].m_sa, lcp) <= 0) {
      high = mid;
      rlcp = lcp;
    } else {
      low = mid;
      llcp = lcp;
    }
  }

  newleft = high;

  if (rlcp >= pat_length) {
    high = right;
    rlcp = old_pat_length;

    while (low + 1 != high) {

      // Invariant: newright is in the range (low, high].
      std::uint64_t lcp = std::min(llcp, rlcp);
      std::uint64_t mid = 0;
      if (llcp + min_discrepancy < rlcp) {
        std::uint64_t d = rlcp - llcp;
        std::uint64_t logd = utils::log2ceil(d);
        mid = low + 1 +
          ((high - low - 1) * balancing_factor * logd) /
          (d + balancing_factor * logd);
      } else if (rlcp + min_discrepancy < llcp) {
        std::uint64_t d = llcp - rlcp;
        std::uint64_t logd = utils::log2ceil(d);
        mid = high - 1 -
          ((high - low - 1) * balancing_factor * logd) /
          (d + balancing_factor * logd);
      } else mid = (low + high) / 2;

      if (lcp_compare(text, pat, pat_length,
            block_beg + (std::uint64_t)block_psa[mid].m_sa, lcp) < 0) {
        high = mid;
        rlcp = lcp;
      } else {
        low = mid;
        llcp = lcp;
      }
    }
  }

  newright = high;
}


//=============================================================================
// NOTE: range returned by this function is NOT the range of suffixes
// starting in the block that has pat[0..pat_length) as a prefix. This would
// only be the case, if we didn't use gt_begin for the tail.
// Instead what we get is a subrange of the above: the range of suffixes
// starting inside the block that have pat[0..pat_length) as a prefix and
// we could not determine whether they are larger than the whole pattern
// or not. This is because we stop the symbol comparisons (when computing
// lcp value) at the end of the text and then use gt_begin for the tail.
// This may cause the we classify the suffix with pat[0..pat_length) as a
// prefix already at this point as larger than the whole pattern because
// gt_begin encodes that information.
//=============================================================================
template<typename block_offset_type>
void refine_range(
    const std::uint8_t *text,
    std::uint64_t text_length,
    std::uint64_t tail_gt_begin_reversed_length,
    std::uint64_t block_beg,
    const bwtsa_t<block_offset_type> *block_psa,
    std::uint64_t left,
    std::uint64_t right,
    const multifile *tail_gt_begin_reversed,
    std::uint64_t old_pat_length,
    std::uint64_t pat_length,
    const std::uint8_t *pat,
    std::uint64_t &newleft,
    std::uint64_t &newright) {

  multifile_bit_stream_reader reader(tail_gt_begin_reversed);

  std::int64_t low = (std::int64_t)left - 1;
  std::int64_t high = right;
  std::uint64_t llcp = old_pat_length;
  std::uint64_t rlcp = old_pat_length;

#ifdef BLOCK_MATRIX_MODULE_DEBUG_MODE
  std::uint64_t min_discrepancy =
    utils::random_int64((std::int64_t)0, (std::int64_t)10);
  std::uint64_t balancing_factor =
    utils::random_int64((std::int64_t)1, (std::int64_t)10);
#else
  static const std::uint64_t min_discrepancy = (1 << 16);
  static const std::uint64_t balancing_factor = 64;
#endif

  while (low + 1 != high) {

    // Invariant: newleft is in the range (low, high].
    std::uint64_t lcp = std::min(llcp, rlcp);
    std::uint64_t mid = 0;
    if (llcp + min_discrepancy < rlcp) {
      std::uint64_t d = rlcp - llcp;
      std::uint64_t logd = utils::log2ceil(d);
      mid = low + 1 +
        ((high - low - 1) * balancing_factor * logd) /
        (d + balancing_factor * logd);
    } else if (rlcp + min_discrepancy < llcp) {
      std::uint64_t d = llcp - rlcp;
      std::uint64_t logd = utils::log2ceil(d);
      mid = high - 1 -
        ((high - low - 1) * balancing_factor * logd) /
        (d + balancing_factor * logd);
    } else mid = (low + high) / 2;

    if (lcp_compare(text, text_length, pat, pat_length,
          tail_gt_begin_reversed_length, block_beg +
          (std::uint64_t)block_psa[mid].m_sa, reader, lcp) <= 0) {
      high = mid;
      rlcp = lcp;
    } else {
      low = mid;
      llcp = lcp;
    }
  }
  newleft = high;

  if (rlcp >= pat_length) {
    high = right;
    rlcp = old_pat_length;

    while (low + 1 != high) {

      // Invariant: newright is in the range (low, high].
      std::uint64_t lcp = std::min(llcp, rlcp);
      std::uint64_t mid = 0;
      if (llcp + min_discrepancy < rlcp) {
        std::uint64_t d = rlcp - llcp;
        std::uint64_t logd = utils::log2ceil(d);
        mid = low + 1 +
          ((high - low - 1) * balancing_factor * logd) /
          (d + balancing_factor * logd);
      } else if (rlcp + min_discrepancy < llcp) {
        std::uint64_t d = llcp - rlcp;
        std::uint64_t logd = utils::log2ceil(d);
        mid = high - 1 -
          ((high - low - 1) * balancing_factor * logd) /
          (d + balancing_factor * logd);
      } else mid = (low + high) / 2;

      if (lcp_compare(text, text_length, pat, pat_length,
            tail_gt_begin_reversed_length, block_beg +
            (std::uint64_t)block_psa[mid].m_sa, reader, lcp) < 0) {
        high = mid;
        rlcp = lcp;
      } else {
        low = mid;
        llcp = lcp;
      }
    }
  }
  newright = high;
}

//=============================================================================
// Variant 1: compute ranges for columns other than the last two.
//=============================================================================
template<typename block_offset_type>
void compute_ranges_1(
    const std::uint8_t *text,
    std::uint64_t text_length,
    const bwtsa_t<block_offset_type> *bwtsa,
    std::uint64_t max_block_size,
    std::pair<std::uint64_t, std::uint64_t> **primary_range,
    std::pair<std::uint64_t, std::uint64_t> **secondary_range,
    std::uint64_t row,
    std::uint64_t column) {

  std::uint64_t n_blocks = (text_length + max_block_size - 1) / max_block_size;
  std::uint64_t block_end = text_length - (n_blocks - 1 - row) * max_block_size;
  std::uint64_t block_begin = std::max((std::int64_t)0,
      (std::int64_t)block_end - (std::int64_t)max_block_size);
  std::uint64_t block_size = block_end - block_begin;
  std::uint64_t pat_start = text_length -
    (n_blocks - 1 - column) * max_block_size;

  const std::uint8_t *pat = text + pat_start;
  const bwtsa_t<block_offset_type> *block_psa = bwtsa + block_begin;

  // Check that 0 <= row < column < n_blocks - 2 and
  // pat_start + 2 * max_block_size <= text_length.
  if (0 > row || row >= column || column >= n_blocks - 2 ||
      pat_start + 2L * max_block_size > text_length) {
    fprintf(stdout, "\n\nError: invariant in compute_ranges_1 failed.\n");
    std::fflush(stdout);
    std::exit(EXIT_FAILURE);
  }

  std::uint64_t left = 0;
  std::uint64_t right = block_size;
  std::uint64_t cur_pat_length = 0;

  // Compute the primary range.
  {
    std::uint64_t new_pat_length = max_block_size;
    if (left != right && cur_pat_length < new_pat_length) {
      std::uint64_t newleft = 0;
      std::uint64_t newright = 0;
      refine_range(text, block_begin, block_psa, left, right,
          cur_pat_length, new_pat_length, pat, newleft, newright);
      left = newleft;
      right = newright;
    }
    cur_pat_length = new_pat_length;
  }
  primary_range[row][column] = std::make_pair(left, right);

#ifdef BLOCK_MATRIX_MODULE_DEBUG_MODE

  // Verify the primary range.
  {
    std::uint64_t smaller = 0;
    std::uint64_t equal = 0;
    for (std::uint64_t j = block_begin; j < block_end; ++j) {
      std::uint64_t lcp = 0L;
      while (lcp < max_block_size && text[j + lcp] == pat[lcp]) ++lcp;
      if (lcp == max_block_size) ++equal;
      else if (text[j + lcp] < pat[lcp]) ++smaller;
    }
    std::uint64_t check_left = smaller;
    std::uint64_t check_right = smaller + equal;
    if (primary_range[row][column] !=
        std::make_pair(check_left, check_right)) {
      fprintf(stdout, "\n\nError: incorrect primary range!\n");
      std::fflush(stdout);
      std::exit(EXIT_FAILURE);
    }
  }
#endif

  // Compute secondary range.
  {
    std::uint64_t new_pat_length = cur_pat_length + max_block_size;
    if (left != right && cur_pat_length < new_pat_length) {
      std::uint64_t newleft = 0;
      std::uint64_t newright = 0;
      refine_range(text, block_begin, block_psa, left, right,
          cur_pat_length, new_pat_length, pat, newleft, newright);
      left = newleft;
      right = newright;
    }
    cur_pat_length = new_pat_length;
  }
  secondary_range[row][column] = std::make_pair(left, right);

#ifdef BLOCK_MATRIX_MODULE_DEBUG_MODE

  // Verify the secondary range.
  {
    std::uint64_t smaller = 0;
    std::uint64_t equal = 0;
    for (std::uint64_t j = block_begin; j < block_end; ++j) {
      std::uint64_t lcp = 0;
      while (lcp < cur_pat_length && text[j + lcp] == pat[lcp]) ++lcp;
      if (lcp == cur_pat_length) ++equal;
      else if (text[j + lcp] < pat[lcp]) ++smaller;
    }
    std::uint64_t check_left = smaller;
    std::uint64_t check_right = smaller + equal;
    if (secondary_range[row][column] !=
        std::make_pair(check_left, check_right)) {
      fprintf(stdout, "\n\nError: incorrect secondary range!\n");
      std::fflush(stdout);
      std::exit(EXIT_FAILURE);
    }
  }
#endif
}

//=============================================================================
// Variant 2: compute primary and secondary range for second to last column.
//=============================================================================
template<typename block_offset_type>
void compute_ranges_2(
    const std::uint8_t *text,
    std::uint64_t text_length,
    std::uint64_t text_beg,
    std::uint64_t supertext_length,
    const bwtsa_t<block_offset_type> *bwtsa,
    std::uint64_t max_block_size,
    background_block_reader *reader,
    const std::uint8_t *next_block,
    std::pair<std::uint64_t, std::uint64_t> **primary_range,
    std::pair<std::uint64_t, std::uint64_t> **secondary_range,
    std::uint64_t row,
    std::uint64_t column) {

  std::uint64_t text_end = text_beg + text_length;
  std::uint64_t tail_length = supertext_length - text_end;
  std::uint64_t n_blocks = (text_length + max_block_size - 1) / max_block_size;
  std::uint64_t block_end = text_length - (n_blocks - 1 - row) * max_block_size;
  std::uint64_t block_begin = std::max((std::int64_t)0,
      (std::int64_t)block_end - (std::int64_t)max_block_size);
  std::uint64_t block_size = block_end - block_begin;
  std::uint64_t pat_start =
    text_length - (n_blocks - 1 - column) * max_block_size;

  const std::uint8_t *pat = text + pat_start;
  const bwtsa_t<block_offset_type> *block_psa = bwtsa + block_begin;

  // Check that 0 <= row < column and column == n_blocks - 2
  // and pat_start + max_block_size == text_length.
  if (0 > row || row >= column || column != n_blocks - 2 ||
        pat_start + max_block_size != text_length) {
    fprintf(stdout, "\n\nError: invariant in compute_ranges_2 failed.\n");
    std::fflush(stdout);
    std::exit(EXIT_FAILURE);
  }

  std::uint64_t left = 0;
  std::uint64_t right = block_size;
  std::uint64_t cur_pat_length = 0;

  // Compute primary range.
  {
    std::uint64_t new_pat_length = max_block_size;
    if (left != right && cur_pat_length < new_pat_length) {
      std::uint64_t newleft = 0;
      std::uint64_t newright = 0;
      refine_range(text, block_begin, block_psa, left, right,
          cur_pat_length, new_pat_length, pat, newleft, newright);
      left = newleft;
      right = newright;
    }
    cur_pat_length = new_pat_length;
  }
  primary_range[row][column] = std::make_pair(left, right);

#ifdef BLOCK_MATRIX_MODULE_DEBUG_MODE

  // Verify the primary range.
  {
    std::uint64_t smaller = 0;
    std::uint64_t equal = 0;
    for (std::uint64_t j = block_begin; j < block_end; ++j) {
      std::uint64_t lcp = 0;
      while (lcp < cur_pat_length && text[j + lcp] == pat[lcp]) ++lcp;
      if (lcp == cur_pat_length) ++equal;
      else if (text[j + lcp] < pat[lcp]) ++smaller;
    }
    std::uint64_t check_left = smaller;
    std::uint64_t check_right = smaller + equal;
    if (primary_range[row][column] !=
        std::make_pair(check_left, check_right)) {
      fprintf(stdout, "\n\nError: incorrect primary range!\n");
      std::fflush(stdout);
      std::exit(EXIT_FAILURE);
    }
  }
#endif

  static const std::uint64_t chunk_size = (1 << 20);

  // Compute secondary range.
  std::uint64_t pat_length = cur_pat_length +
    std::min(tail_length, max_block_size);

  if (reader) {

    // The reader != NULL, meaning that
    // we have to gradually refine the range.
    while (left != right && cur_pat_length < pat_length) {
      std::uint64_t next_chunk =
        std::min(chunk_size,
            pat_length - cur_pat_length);
      std::uint64_t new_pat_length = cur_pat_length + next_chunk;
      reader->wait(new_pat_length - max_block_size);

      std::uint64_t newleft = 0;
      std::uint64_t newright = 0;
      refine_range(text, block_begin, block_psa,
          left, right, cur_pat_length, new_pat_length,
          reader->m_data - max_block_size, newleft, newright);
      left = newleft;
      right = newright;
      cur_pat_length = new_pat_length;
    }
  } else {

#ifdef BLOCK_MATRIX_MODULE_DEBUG_MODE

    // Debug version extends the range chunk by
    // chunk (using random chunk lengths) even if
    // the whole next block is available.
    while (left != right && cur_pat_length < pat_length) {
      std::uint64_t next_chunk =
        utils::random_int64((std::int64_t)1,
            (std::int64_t)(pat_length - cur_pat_length));
      std::uint64_t new_pat_length = cur_pat_length + next_chunk;

      std::uint64_t newleft = 0;
      std::uint64_t newright = 0;
      refine_range(text, block_begin, block_psa,
          left, right, cur_pat_length, new_pat_length,
          next_block - max_block_size, newleft, newright);
      left = newleft;
      right = newright;
      cur_pat_length = new_pat_length;
    }
#else

    // The whole next block is available,
    // we can just do one binary search.
    std::uint64_t new_pat_length = pat_length;
    if (left != right && cur_pat_length < new_pat_length) {
      std::uint64_t newleft = 0;
      std::uint64_t newright = 0;
      refine_range(text, block_begin, block_psa,
          left, right, cur_pat_length, new_pat_length,
          next_block - max_block_size, newleft, newright);
      left = newleft;
      right = newright;
    }
    cur_pat_length = new_pat_length;
#endif

  }

  secondary_range[row][column] =
    std::make_pair(left, right);
}

//=============================================================================
// Variant 3: compute primary and secondary range for the last column.
//=============================================================================
template<typename block_offset_type>
void compute_ranges_3(
    const std::uint8_t *text,
    std::uint64_t text_length,
    std::uint64_t text_beg,
    std::uint64_t supertext_length,
    const bwtsa_t<block_offset_type> *bwtsa,
    std::uint64_t max_block_size,
    const multifile *tail_gt_begin_reversed,
    background_block_reader *reader,
    const std::uint8_t *next_block,
    std::pair<std::uint64_t, std::uint64_t> **primary_range,
    std::pair<std::uint64_t, std::uint64_t> **secondary_range,
    std::uint64_t row,
    std::uint64_t column) {

  std::uint64_t text_end = text_beg + text_length;
  std::uint64_t tail_length = supertext_length - text_end;
  std::uint64_t n_blocks = (text_length + max_block_size - 1) / max_block_size;
  std::uint64_t block_end = text_length - (n_blocks - 1 - row) * max_block_size;
  std::uint64_t block_beg = std::max((std::int64_t)0,
      (std::int64_t)block_end - (std::int64_t)max_block_size);
  std::uint64_t block_size = block_end - block_beg;
  const bwtsa_t<block_offset_type> *block_psa = bwtsa + block_beg;
  std::uint64_t first_range_pat_length =
    std::min(max_block_size, tail_length);

  // Length of text stored in next_block (if not NULL).
  std::uint64_t pat_length = std::min(text_length, tail_length);

  // Note: max_block_size <= text_length thus
  // first_range_pat_length <= pat_length
  // Invariant: one of the following cases holds:
  // (1) next_block != NULL and reader == NULL and next_block stores
  //     std::min(text_length, tail_length) symbols after text
  // (2) next_block == NULL and reader != NULL and reader will read
  //     std::min(text_length, tail_length) symbols after text

  // Check that 0 <= row < colum and column == n_blocks - 1.
  if (0 > row || row >= column || column != n_blocks - 1) {
    fprintf(stdout, "\n\nError: invariant 1 in compute_ranges_3 failed.\n");
    std::fflush(stdout);
    std::exit(EXIT_FAILURE);
  }

  std::uint64_t left = 0;
  std::uint64_t right = block_size;
  std::uint64_t cur_pat_length = 0;

  static const std::uint64_t chunk_size = (1 << 20);

  // Compute the primary range.
  if (reader) {

    // The reader != NULL, meaning that we
    // have to gradually refine the range.
    while (left != right && cur_pat_length < first_range_pat_length) {
      std::uint64_t next_chunk =
        std::min(chunk_size,
            first_range_pat_length - cur_pat_length);
      std::uint64_t new_pat_length = cur_pat_length + next_chunk;
      reader->wait(new_pat_length);

      std::uint64_t newleft = 0;
      std::uint64_t newright = 0;
      refine_range(text, text_length, tail_length, block_beg,
          block_psa, left, right, tail_gt_begin_reversed,
          cur_pat_length, new_pat_length, reader->m_data,
          newleft, newright);
      left = newleft;
      right = newright;
      cur_pat_length = new_pat_length;
    }
  } else {

#ifdef BLOCK_MATRIX_MODULE_DEBUG_MODE

    // Debug version extends the range chunk by
    // chunk (using random chunk lengths) even if
    // the whole next block is available.
    while (left != right && cur_pat_length < first_range_pat_length) {
      std::uint64_t next_chunk = utils::random_int64((std::int64_t)1,
          (std::int64_t)(first_range_pat_length - cur_pat_length));
      std::uint64_t new_pat_length = cur_pat_length + next_chunk;

      std::uint64_t newleft = 0;
      std::uint64_t newright = 0;
      refine_range(text, text_length, tail_length, block_beg,
          block_psa, left, right, tail_gt_begin_reversed,
          cur_pat_length, new_pat_length, next_block, newleft,
          newright);
      left = newleft;
      right = newright;
      cur_pat_length = new_pat_length;
    }
#else

    // The whole next block is available,
    // we can just do one binary search.
    std::uint64_t new_pat_length = first_range_pat_length;
    if (left != right && cur_pat_length < new_pat_length) {
      std::uint64_t newleft = 0;
      std::uint64_t newright = 0;
      refine_range(text, text_length, tail_length, block_beg,
          block_psa, left, right, tail_gt_begin_reversed,
          cur_pat_length, new_pat_length, next_block, newleft,
          newright);
      left = newleft;
      right = newright;
    }
    cur_pat_length = new_pat_length;
#endif

  }

  primary_range[row][column] =
    std::make_pair(left, right);

  // Compute the secondary range.
  if (reader) {

    // The reader != NULL, meaning that
    // we have to gradually refine the range.
    while (left != right && cur_pat_length < pat_length) {
      std::uint64_t next_chunk =
        std::min(chunk_size,
            pat_length - cur_pat_length);
      std::uint64_t new_pat_length = cur_pat_length + next_chunk;
      reader->wait(new_pat_length);

      std::uint64_t newleft = 0;
      std::uint64_t newright = 0;
      refine_range(text, text_length, tail_length, block_beg,
          block_psa, left, right, tail_gt_begin_reversed,
          cur_pat_length, new_pat_length, reader->m_data,
          newleft, newright);
      left = newleft;
      right = newright;
      cur_pat_length = new_pat_length;
    }
  } else {

#ifdef BLOCK_MATRIX_MODULE_DEBUG_MODE

    // Debug version extends the range chunk by
    // chunk (using random chunk lengths) even
    // if the whole next block is available.
    while (left != right && cur_pat_length < pat_length) {
      std::uint64_t next_chunk = utils::random_int64((std::int64_t)1,
         (std::int64_t)(pat_length - cur_pat_length));
      std::uint64_t new_pat_length = cur_pat_length + next_chunk;

      std::uint64_t newleft = 0;
      std::uint64_t newright = 0;
      refine_range(text, text_length, tail_length, block_beg,
          block_psa, left, right, tail_gt_begin_reversed,
          cur_pat_length, new_pat_length, next_block, newleft,
          newright);
      left = newleft;
      right = newright;
      cur_pat_length = new_pat_length;
    }
#else

    // The whole next block is available,
    // we can just do one binary search.
    std::uint64_t new_pat_length = pat_length;
    if (left != right && cur_pat_length < new_pat_length) {
      std::uint64_t newleft = 0;
      std::uint64_t newright = 0;
      refine_range(text, text_length, tail_length, block_beg,
          block_psa, left, right, tail_gt_begin_reversed,
          cur_pat_length, new_pat_length, next_block, newleft,
          newright);
      left = newleft;
      right = newright;
    }
    cur_pat_length = new_pat_length;
#endif

  }

  secondary_range[row][column] =
    std::make_pair(left, right);

  if (left != right && text_length <= tail_length) {
    fprintf(stdout, "\n\nError: left != right && "
        "text_length <= tail_length.\n");
    std::fflush(stdout);
    std::exit(EXIT_FAILURE);
  }
}

template<typename block_offset_type>
void task_solver_code(
    const std::uint8_t *text,
    std::uint64_t text_length,
    const bwtsa_t<block_offset_type> *bwtsa,
    std::uint64_t max_block_size,
    std::pair<std::uint64_t, std::uint64_t> **primary_range,
    std::pair<std::uint64_t, std::uint64_t> **secondary_range,
    std::vector<std::pair<std::uint64_t, std::uint64_t> > &tasks,
    std::mutex &tasks_mutex) {

  while (true) {

    // Get a task from the task collection.
    std::pair<std::uint64_t, std::uint64_t> task;
    bool task_avail = true;
    std::unique_lock<std::mutex> lk(tasks_mutex);
    if (tasks.empty()) task_avail = false;
    else {
      task = tasks.back();
      tasks.pop_back();
    }
    lk.unlock();

    if (!task_avail) break;

    // Solve the task and save the answer.
    compute_ranges_1(text, text_length, bwtsa, max_block_size,
        primary_range, secondary_range, task.first, task.second);
  }
}

template<typename block_offset_type>
void compute_block_rank_matrix(
    const std::uint8_t *text,
    std::uint64_t text_length,
    const bwtsa_t<block_offset_type> *bwtsa,
    std::uint64_t max_block_size,
    std::uint64_t text_beg,
    std::uint64_t supertext_length,
    const multifile *tail_gt_begin_reversed,
    background_block_reader *reader,
    const std::uint8_t *next_block,
    std::uint64_t **block_rank_matrix) {

  std::uint64_t n_blocks = (text_length + max_block_size - 1) / max_block_size;
  std::uint64_t text_end = text_beg + text_length;
  std::uint64_t tail_length = supertext_length - text_end;

  // Allocate primary and secondary ranges.
  std::pair<std::uint64_t, std::uint64_t> **primary_range =
    new std::pair<std::uint64_t, std::uint64_t>*[n_blocks];
  std::pair<std::uint64_t, std::uint64_t> **secondary_range =
    new std::pair<std::uint64_t, std::uint64_t>*[n_blocks];
  for (std::uint64_t row = 0; row < n_blocks; ++row) {
    primary_range[row] = new std::pair<std::uint64_t, std::uint64_t>[n_blocks];
    secondary_range[row] = new std::pair<std::uint64_t, std::uint64_t>[n_blocks];
  }

  // Start the threads computing
  // ranges for the last column.
  std::thread **threads_last_col = NULL;
  if (n_blocks > 1) {
    threads_last_col = new std::thread*[n_blocks - 1];
    for (std::uint64_t row = 0; row + 1 < n_blocks; ++row) {
      std::uint64_t column = n_blocks - 1;
      threads_last_col[row] = new std::thread(
          compute_ranges_3<block_offset_type>, text, text_length,
          text_beg, supertext_length, bwtsa, max_block_size,
          tail_gt_begin_reversed, reader, next_block, primary_range,
          secondary_range, row, column);
    }
  }

  // Start the threads computing ranges
  // for the second-to-last column.
  std::thread **threads_second_last_col = NULL;
  if (n_blocks > 2) {
    threads_second_last_col = new std::thread*[n_blocks - 2];
    for (std::uint64_t row = 0; row + 2 < n_blocks; ++row) {
      std::uint64_t column = n_blocks - 2;
      threads_second_last_col[row] = new std::thread(
          compute_ranges_2<block_offset_type>, text, text_length,
          text_beg, supertext_length, bwtsa, max_block_size,
          reader, next_block, primary_range, secondary_range,
          row, column);
    }
  }

  // Start threads computing columns
  // other than the last two.
  std::vector<std::pair<std::uint64_t, std::uint64_t> > tasks;
  std::mutex tasks_mutex;
  for (std::uint64_t row = 0; row < n_blocks; ++row)
    for (std::uint64_t col = row + 1; col + 2 < n_blocks; ++col)
      tasks.push_back(std::make_pair(row, col));
  std::random_shuffle(tasks.begin(), tasks.end());  // solve in any order
  std::thread **threads_other = new std::thread*[n_blocks];
  for (std::uint64_t t = 0; t < n_blocks; ++t)
    threads_other[t] = new std::thread(
        task_solver_code<block_offset_type>, text, text_length,
        bwtsa, max_block_size, primary_range, secondary_range,
        std::ref(tasks), std::ref(tasks_mutex));

  // Wait for the threads computing
  // columns other than last two.
  for (std::uint64_t t = 0; t < n_blocks; ++t) threads_other[t]->join();
  for (std::uint64_t t = 0; t < n_blocks; ++t) delete threads_other[t];
  delete[] threads_other;

  // Wait for the threads computing
  // second-to-last column to finish.
  if (n_blocks > 2) {
    for (std::uint64_t row = 0; row + 2 < n_blocks; ++row)
      threads_second_last_col[row]->join();
    for (std::uint64_t row = 0; row + 2 < n_blocks; ++row)
      delete threads_second_last_col[row];
    delete[] threads_second_last_col;
  }

  // Wait for the threads computing
  // the last column to finish.
  if (n_blocks > 1) {
    for (std::uint64_t row = 0; row + 1 < n_blocks; ++row)
      threads_last_col[row]->join();
    for (std::uint64_t row = 0; row + 1 < n_blocks; ++row)
      delete threads_last_col[row];
    delete[] threads_last_col;
  }

  // Compute the rank values from
  // primary and secondary ranges.
  for (std::uint64_t row_plus = n_blocks; row_plus > 0; --row_plus) {
    std::uint64_t row = row_plus - 1;
    for (std::uint64_t col = n_blocks - 1; col > row; --col) {
      std::uint64_t left = secondary_range[row][col].first;
      std::uint64_t right = secondary_range[row][col].second;

      if (col != n_blocks - 1 &&
          (col != n_blocks - 2 || tail_length >= max_block_size)) {

        std::uint64_t cur_block_end =
          text_length - (n_blocks - 1 - row) * max_block_size;
        std::uint64_t cur_block_beg = std::max((std::int64_t)0,
            (std::int64_t)cur_block_end - (std::int64_t)max_block_size);
        std::uint64_t cur_block_size = cur_block_end - cur_block_beg;
        std::uint64_t shift = max_block_size - cur_block_size;
        std::uint64_t next_block_end =
          text_length - (n_blocks - 1 - (row + 1)) * max_block_size;
        std::uint64_t next_block_beg = std::max((std::int64_t)0,
            (std::int64_t)next_block_end - (std::int64_t)max_block_size);

        const bwtsa_t<block_offset_type> *cur_block_psa =
          bwtsa + cur_block_beg;
        const bwtsa_t<block_offset_type> *next_block_psa =
          bwtsa + next_block_beg;

        // Compute the ranges.
        std::uint64_t next_primary_range_beg =
          primary_range[row + 1][col + 1].first;
        std::uint64_t next_primary_range_end =
          primary_range[row + 1][col + 1].second;
        std::uint64_t next_primary_range_size =
          next_primary_range_end - next_primary_range_beg;

        // Compute the difference of the arithmetic progression.
        std::int64_t delta = 0;
        std::uint64_t next_psa_first = 0;
        std::uint64_t next_psa_second = 0;
        if (next_primary_range_size > 1) {
          next_psa_first = next_block_psa[next_primary_range_beg].m_sa;
          next_psa_second = next_block_psa[next_primary_range_beg + 1].m_sa;
          delta = (std::int64_t)next_psa_second - (std::int64_t)next_psa_first;
        }

        // Invariant:
        // 1. the primary range of next block contains (possibly
        //    zero) values forming an arithmetic progression,
        // 2. elements in the range [left..right) of the psa of the
        //    current block incremented by `shift' appear in the primary
        //    range of the next block.

#ifdef BLOCK_MATRIX_MODULE_DEBUG_MODE

        // Check that both invariants hold.
        for (std::uint64_t j = next_primary_range_beg;
            j + 1 < next_primary_range_end; ++j) {
          if ((std::int64_t)next_block_psa[j + 1].m_sa -
              (std::int64_t)next_block_psa[j].m_sa != delta) {
            fprintf(stdout, "\n\nError: invariant 1 in "
                "compute_block_rank_matrix failed.\n");
            std::fflush(stdout);
            std::exit(EXIT_FAILURE);
          }
        }

        for (std::uint64_t j = left; j < right; ++j) {
          std::uint64_t suf = (std::uint64_t)cur_block_psa[j].m_sa + shift;
          bool found = false;
          for (std::uint64_t jj = next_primary_range_beg;
              jj < next_primary_range_end; ++jj) {
            if ((std::uint64_t)next_block_psa[jj].m_sa == suf) {
              found = true;
              break;
            }
          }

          if (!found) {
            fprintf(stdout, "\n\nError: invariant 2 in "
                "compute_block_rank_matrix failed.\n");
            std::fflush(stdout);
            std::exit(EXIT_FAILURE);
          }
        }
#endif

        // Keep refining the range [left..right) until it's empty.
        while (left != right) {

          // Valid values for mid are in [left..right).
          std::uint64_t mid = (left + right) / 2;
          std::uint64_t suf = (std::uint64_t)cur_block_psa[mid].m_sa + shift;

          // Locate suf in next_block_psa using invariants 1. and 2.
          std::int64_t pos = next_primary_range_beg;
          if (next_primary_range_size > 1)
            pos += ((std::int64_t)suf - (std::int64_t)next_psa_first) / delta;

          // Refine the range.
          if ((std::uint64_t)pos <
              block_rank_matrix[row + 1][col + 1]) left = mid + 1;
          else right = mid;
        }
      }

      block_rank_matrix[row][col] = left;
    }
  }

  // Clean up.
  for (std::uint64_t row = 0; row < n_blocks; ++row) {
    delete[] primary_range[row];
    delete[] secondary_range[row];
  }
  delete[] primary_range;
  delete[] secondary_range;
}

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_INMEM_COMPUTE_INITIAL_RANKS_HPP_INCLUDED
