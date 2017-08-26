/**
 * @file    src/psascan_src/em_compute_initial_ranks.hpp
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

#ifndef __SRC_PSASCAN_SRC_EM_COMPUTE_INITIAL_RANKS_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_EM_COMPUTE_INITIAL_RANKS_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>
#include <thread>

#include "approx_rank.hpp"
#include "sparse_isa.hpp"
#include "utils.hpp"
#include "io/background_block_reader.hpp"
#include "io/background_chunk_reader.hpp"
#include "io/multifile_bit_stream_reader.hpp"


namespace psascan_private {

// #define EM_STARTING_POS_MODULE_DEBUG_MODE

inline int lcp_compare(
    const std::uint8_t *text,    // only text[block_suf_beg..block_end)
    std::uint64_t text_length,   //   can be accessed
    std::uint64_t block_end,     // wrt to text beg
    std::uint64_t block_suf_beg, // wrt to text beg
    const std::uint8_t *pat,     // only pat[lcp..pat_length) can be accessed
    std::uint64_t pat_beg,       // wrt to text beg
    std::uint64_t pat_length,
    multifile_bit_stream_reader &gt_reader,
    std::uint64_t &lcp) {

  // Continue matching until we reach the end of
  // block of the end of the pattern.
  while (block_suf_beg + lcp < block_end &&
      lcp < pat_length &&
      text[block_suf_beg + lcp] == pat[lcp])
    ++lcp;

  if (block_suf_beg + lcp >= block_end) {
    std::uint64_t block_suf_len = block_end - block_suf_beg;

    // Invariant: lcp >= block_suf_len.
    std::uint64_t pat_ptr = pat_beg + block_suf_len;
    if (gt_reader.access(text_length - pat_ptr)) return 1;
    else return -1;
  } else if (lcp == pat_length) {
    if (pat_beg + pat_length >= text_length) return -1;
    else return 0;
  } else {
    if (pat[lcp] > text[block_suf_beg + lcp]) return 1;
    else return -1;
  } 
}

template<typename block_offset_type>
void refine_range(
    const std::uint8_t *block,
    const block_offset_type *block_psa,
    std::uint64_t block_beg,  // wrt to text beg
    std::uint64_t block_end,  // same here
    std::uint64_t pat_beg,    // same here
    std::uint64_t text_length,
    std::uint64_t left,
    std::uint64_t right,
    std::uint64_t old_lcp,
    std::uint64_t new_lcp,
    const std::uint8_t *pat,  // only pat[old_lcp..new_lcp) can be accessed
    multifile_bit_stream_reader &gt_reader,
    std::uint64_t &newleft,
    std::uint64_t &newright) {

  std::int64_t low = (std::int64_t)left - 1;
  std::int64_t high = right;
  std::uint64_t llcp = old_lcp;
  std::uint64_t rlcp = old_lcp;

#ifdef EM_STARTING_POS_MODULE_DEBUG_MODE
  std::uint64_t min_discrepancy =
    utils::random_int64((std::int64_t)0, (std::int64_t)10);
  std::uint64_t balancing_factor =
    utils::random_int64((std::int64_t)1, (std::int64_t)10);
#else
  static const std::uint64_t min_discrepancy = (1 << 16);
  static const std::uint64_t balancing_factor = 64;
#endif

  const std::uint8_t *text = block - block_beg;
  while (low + 1 != high) {

    // Invariant: newleft is in the range (low, high].
    std::uint64_t lcp = std::min(llcp, rlcp);
    std::uint64_t mid = 0L;
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

    if (lcp_compare(text, text_length, block_end,
          block_beg + (std::uint64_t)block_psa[mid],
          pat, pat_beg, new_lcp, gt_reader, lcp) <= 0) {
      high = mid;
      rlcp = lcp;
    } else {
      low = mid;
      llcp = lcp;
    }
  }
  newleft = high;

  if (rlcp >= new_lcp) {
    high = right;
    rlcp = old_lcp;

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

      if (lcp_compare(text, text_length, block_end,
            block_beg + (std::uint64_t)block_psa[mid],
            pat, pat_beg, new_lcp, gt_reader, lcp) < 0) {
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

template<typename block_offset_type>
void em_compute_single_initial_rank(
    const std::uint8_t *block,
    const block_offset_type *block_psa,
    std::uint64_t block_beg,  // wrt to text beg
    std::uint64_t block_end,  // same here
    std::uint64_t pat_beg,    // same here
    std::uint64_t text_length,
    std::uint64_t max_lcp,
    std::string text_filename,
    const multifile *tail_gt_begin_reversed,
    std::pair<std::uint64_t, std::uint64_t> &result) {

  if (pat_beg == text_length) {
    result = std::make_pair((std::uint64_t)0, (std::uint64_t)0);
    return;
  }

  std::uint64_t block_size = block_end - block_beg;
  std::uint64_t pat_end = pat_beg + max_lcp;

  multifile_bit_stream_reader gt_reader(tail_gt_begin_reversed);

  // Reads text[pat_beg..pat_end) in chunks.
#ifdef EM_STARTING_POS_MODULE_DEBUG_MODE
  std::uint64_t chunk_length =
    utils::random_int64((std::int64_t)1, (std::int64_t)10);
  background_chunk_reader *chunk_reader =
    new background_chunk_reader(
        text_filename, pat_beg, pat_end, chunk_length);
#else
  background_chunk_reader *chunk_reader =
    new background_chunk_reader(text_filename, pat_beg, pat_end);
#endif

  // The current range is [left, right).
  std::uint64_t left = 0;
  std::uint64_t right = block_size;
  std::uint64_t lcp = 0;

  while (left != right && lcp < max_lcp) {
    std::uint64_t this_chunk_length =
      std::min(max_lcp - lcp, chunk_reader->get_chunk_size());
    std::uint64_t new_lcp = lcp + this_chunk_length;
    chunk_reader->wait(pat_beg + new_lcp);

    // Invariant:
    // reader->chunk[0..chunk_length) = pattern[lcp..new_lcp).
    std::uint64_t newleft = 0;
    std::uint64_t newright = 0;
    refine_range(block, block_psa, block_beg, block_end,
        pat_beg, text_length, left, right, lcp, new_lcp,
        chunk_reader->m_chunk - lcp, gt_reader, newleft,
        newright);

    left = newleft;
    right = newright;
    lcp = new_lcp;
  }

  // Clean up.
  delete chunk_reader;

  // Store the answer.
  result = std::make_pair(left, right);
}

template<typename block_offset_type>
void em_compute_initial_ranks(
    const std::uint8_t *block,
    const block_offset_type *block_psa,
    const std::uint8_t *block_pbwt,
    std::uint64_t i0,
    std::uint64_t block_beg,  // wrt to text beg
    std::uint64_t block_end,  // same here
    std::uint64_t text_length,
    std::string text_filename,
    const multifile *tail_gt_begin_reversed,
    std::vector<std::uint64_t> &result,
    std::uint64_t max_threads,
    std::uint64_t tail_end,
    std::uint64_t initial_rank_after_tail) {

  // Note, that bits of tail_gt_begin_reversed are indexed in the
  // range [text_length - tail_end.. text_length - block_end). This
  // is because the same multifile is then used in the streaming and
  // for streaming is much more natural to use this indexing.
  std::uint64_t block_length = block_end - block_beg;
  std::uint64_t tail_length = tail_end - block_end;
  std::uint64_t stream_max_block_size =
    (tail_length + max_threads - 1) / max_threads;
  std::uint64_t n_threads =
    (tail_length + stream_max_block_size - 1) / stream_max_block_size;

  std::vector<std::pair<std::uint64_t, std::uint64_t> > ranges(n_threads);
  std::thread **threads = new std::thread*[n_threads];

  for (std::uint64_t t_plus = n_threads; t_plus > 0; --t_plus) {
    std::uint64_t t = t_plus - 1;
    std::uint64_t stream_block_beg = block_end + t * stream_max_block_size;
    std::uint64_t stream_block_end =
      std::min(stream_block_beg + stream_max_block_size, tail_end);
    std::uint64_t stream_block_size = stream_block_end - stream_block_beg;

    threads[t] = new std::thread(
        em_compute_single_initial_rank<block_offset_type>,
        block, block_psa, block_beg, block_end,
        stream_block_beg, text_length, stream_block_size,
        text_filename, tail_gt_begin_reversed, std::ref(ranges[t]));
  }

  for (std::uint64_t t = 0; t < n_threads; ++t) threads[t]->join();
  for (std::uint64_t t = 0; t < n_threads; ++t) delete threads[t];
  delete[] threads;

  // Refine ranges until all are single elements.
  result.resize(n_threads);

  bool nontrivial_range = false;
  for (std::uint64_t t = 0; t < n_threads; ++t)
    if (ranges[t].first != ranges[t].second)
      nontrivial_range = true;

  if (nontrivial_range) {
    multifile_bit_stream_reader *gt_reader =
      new multifile_bit_stream_reader(tail_gt_begin_reversed);

#ifdef EM_STARTING_POS_MODULE_DEBUG_MODE
    typedef approx_rank<1L> rank_type;
    typedef sparse_isa<rank_type, block_offset_type, 1L> isa_type;
#else
    typedef approx_rank<8L> rank_type;
    typedef sparse_isa<rank_type, block_offset_type, 8L> isa_type;
#endif
    rank_type *pbwt_rank =
      new rank_type(block_pbwt, block_length);
    isa_type *block_sparse_isa =
      new isa_type(block_psa, block, block_length, i0, pbwt_rank);

    std::uint64_t prev_rank = initial_rank_after_tail;
    for (std::uint64_t t_plus = n_threads; t_plus > 0; --t_plus) {
      std::uint64_t t = t_plus - 1;

      // Compute block boundaries.
      std::uint64_t stream_block_beg = block_end + t * stream_max_block_size;
      std::uint64_t stream_block_end =
        std::min(stream_block_beg + stream_max_block_size, tail_end);
      std::uint64_t stream_block_size = stream_block_end - stream_block_beg;

      std::uint64_t left = ranges[t].first;
      std::uint64_t right = ranges[t].second;

      while (left != right) {

        // Valid values for mid are in [left..right).
        std::uint64_t mid = (left + right) / 2;

        if ((std::uint64_t)block_psa[mid] +
            stream_block_size >= block_length) {
          std::uint64_t suf_len = block_length - (std::uint64_t)block_psa[mid];
          if (gt_reader->access(text_length - (stream_block_beg + suf_len)))
            left = mid + 1;
          else right = mid; 
        } else {
          std::uint64_t j = (std::uint64_t)block_psa[mid] + stream_block_size;
          if (block_sparse_isa->query(j) < prev_rank) left = mid + 1;
          else right = mid;
        }
      }

      result[t] = left;
      prev_rank = result[t];
    }

    delete pbwt_rank;
    delete block_sparse_isa;
    delete gt_reader;
  } else {
    for (std::uint64_t t = 0; t < n_threads; ++t)
      result[t] = ranges[t].first;
  }
}

inline int lcp_compare_2(
    const std::uint8_t *text,     // only text[block_suf_beg..block_end)
    std::uint64_t text_length,    //   can be accessed
    std::uint64_t block_end,      // wrt to text beg
    std::uint64_t block_suf_beg,  // wrt to text beg
    const std::uint8_t *pat,      // only pat[lcp..pat_length) can be accessed
    std::uint64_t pat_beg,        // wrt to text beg
    std::uint64_t pat_length,
    std::uint64_t tail_begin,     // wrt to text beg
    background_block_reader *mid_block_reader,
    multifile_bit_stream_reader &gt_reader,
    std::uint64_t &lcp) {

  while (block_suf_beg + lcp < block_end &&
      lcp < pat_length &&
      text[block_suf_beg + lcp] == pat[lcp])
    ++lcp;

  if (block_suf_beg + lcp < block_end && lcp < pat_length) {
    if (pat[lcp] > text[block_suf_beg + lcp]) return 1;
    else return -1;
  }

  if (block_suf_beg + lcp >= block_end &&
      block_end < tail_begin && lcp < pat_length) {

    // To finish the comparison, we need to access symbols from
    // the mid block. First, wait until enough symbols are available.
    mid_block_reader->wait(
        std::min(tail_begin, block_suf_beg + pat_length) - block_end);

    // Now continue the comparison.
    const std::uint8_t *text2 = mid_block_reader->m_data - block_end;
    while (block_suf_beg + lcp < tail_begin && lcp < pat_length &&
        text2[block_suf_beg + lcp] == pat[lcp])
      ++lcp;

    if (block_suf_beg + lcp < tail_begin && lcp < pat_length) {
      if (pat[lcp] > text2[block_suf_beg + lcp]) return 1;
      else return -1;
    }
  }

  if (block_suf_beg + lcp >= tail_begin) {

    // Use gt to resolve comparison.
    if (gt_reader.access(text_length -
          (pat_beg + (tail_begin - block_suf_beg))))
      return 1;
    else return -1;

  } else {

    // lcp == pat_length
    if (pat_beg + pat_length >= text_length) return -1;
    else return 0;
  }
}

template<typename block_offset_type>
void refine_range_2(
    const std::uint8_t *block,
    const block_offset_type *block_psa,
    std::uint64_t block_beg,  // wrt to text beg
    std::uint64_t block_end,  // same here
    std::uint64_t pat_beg,    // same here
    std::uint64_t tail_begin,
    background_block_reader *mid_block_reader,
    std::uint64_t text_length,
    std::uint64_t left,
    std::uint64_t right,
    std::uint64_t old_lcp,
    std::uint64_t new_lcp,
    const std::uint8_t *pat,  // only pat[old_lcp..new_lcp) can be accessed
    multifile_bit_stream_reader &gt_reader,
    std::uint64_t &newleft,
    std::uint64_t &newright) {

  std::int64_t low = (std::int64_t)left - 1;
  std::int64_t high = right;
  std::uint64_t llcp = old_lcp;
  std::uint64_t rlcp = old_lcp;

#ifdef EM_STARTING_POS_MODULE_DEBUG_MODE
  std::uint64_t min_discrepancy = utils::random_int64(0L, 10L);
  std::uint64_t balancing_factor = utils::random_int64(1L, 10L);
#else
  static const std::uint64_t min_discrepancy = (1 << 16);
  static const std::uint64_t balancing_factor = 64;
#endif

  const std::uint8_t *text = block - block_beg;
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

    if (lcp_compare_2(text, text_length, block_end,
          block_beg + (std::uint64_t)block_psa[mid],
        pat, pat_beg, new_lcp, tail_begin,
        mid_block_reader, gt_reader, lcp) <= 0) {
      high = mid;
      rlcp = lcp;
    } else {
      low = mid;
      llcp = lcp;
    }
  }
  newleft = high;

  if (rlcp >= new_lcp) {
    high = right;
    rlcp = old_lcp;

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

      if (lcp_compare_2(text, text_length, block_end,
            block_beg + (std::uint64_t)block_psa[mid],
            pat, pat_beg, new_lcp, tail_begin,
            mid_block_reader, gt_reader, lcp) < 0) {
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

template<typename block_offset_type>
void em_compute_single_initial_rank_2(
    const std::uint8_t *block,
    const block_offset_type *block_psa,
    std::uint64_t block_beg,  // wrt to text beg
    std::uint64_t block_end,  // same here
    std::uint64_t pat_beg,    // same here
    std::uint64_t text_length,
    std::uint64_t max_lcp,
    std::uint64_t tail_begin,
    background_block_reader *mid_block_reader,
    std::string text_filename,
    const multifile *tail_gt_begin_reversed,
    std::uint64_t &result) {

  if (pat_beg == text_length) {
    result = 0;
    return;
  }

  std::uint64_t block_size = block_end - block_beg;
  std::uint64_t pat_end = std::min(text_length, pat_beg + max_lcp);

  multifile_bit_stream_reader gt_reader(tail_gt_begin_reversed);

  // Reads text[pat_beg..pat_end) in chunks.
#ifdef EM_STARTING_POS_MODULE_DEBUG_MODE
  std::uint64_t chunk_length = utils::random_int64(1L, 10L);
  background_chunk_reader *chunk_reader =
    new background_chunk_reader(text_filename, pat_beg, pat_end, chunk_length);
#else
  background_chunk_reader *chunk_reader =
    new background_chunk_reader(text_filename, pat_beg, pat_end);
#endif

  // The current range is [left, right).
  std::uint64_t left = 0;
  std::uint64_t right = block_size;
  std::uint64_t lcp = 0;

  while (left != right && lcp < max_lcp) {
    std::uint64_t this_chunk_length =
      std::min(max_lcp - lcp, chunk_reader->get_chunk_size());
    std::uint64_t new_lcp = lcp + this_chunk_length;
    chunk_reader->wait(pat_beg + new_lcp);

    // Invariant:
    // reader->chunk[0..chunk_length) = pattern[lcp..new_lcp).
    std::uint64_t newleft = 0;
    std::uint64_t newright = 0;
    refine_range_2(block, block_psa, block_beg, block_end,
        pat_beg, tail_begin, mid_block_reader, text_length,
        left, right, lcp, new_lcp, chunk_reader->m_chunk - lcp,
        gt_reader, newleft, newright);
    left = newleft;
    right = newright;
    lcp = new_lcp;
  }
  result = left;

  delete chunk_reader;
}

template<typename block_offset_type>
void em_compute_initial_ranks(
    const std::uint8_t *block,
    const block_offset_type *block_psa,
    std::uint64_t block_beg,  // wrt to text beg
    std::uint64_t block_end,  // same here
    std::uint64_t text_length,
    std::string text_filename,
    const multifile *tail_gt_begin_reversed,
    std::vector<std::uint64_t> &result,
    std::uint64_t max_threads,
    std::uint64_t tail_begin) {

  // Compute some initial parameters.
  std::uint64_t block_length = block_end - block_beg;
  std::uint64_t tail_length = text_length - tail_begin;
  std::uint64_t mid_block_beg = block_end;
  std::uint64_t mid_block_end = tail_begin;
  std::uint64_t mid_block_size = mid_block_end - mid_block_beg;
  std::uint64_t stream_max_block_size =
    (tail_length + max_threads - 1) / max_threads;
  std::uint64_t n_threads =
    (tail_length + stream_max_block_size - 1) / stream_max_block_size;

  // Start reading the text between the block and the tail in the backgrond.
  background_block_reader *mid_block_reader =
    new background_block_reader(text_filename, mid_block_beg, mid_block_size);

  // Compute the initial ranks.
  std::vector<std::uint64_t> res(n_threads);
  std::thread **threads = new std::thread*[n_threads];

  for (std::uint64_t t = 0; t < n_threads; ++t) {
    std::uint64_t stream_block_beg = tail_begin + t * stream_max_block_size;
    std::uint64_t max_lcp = std::min(block_length + mid_block_size,
        text_length - stream_block_beg);

    threads[t] = new std::thread(
        em_compute_single_initial_rank_2<block_offset_type>, block,
        block_psa, block_beg, block_end, stream_block_beg,
        text_length, max_lcp, tail_begin, mid_block_reader,
        text_filename, tail_gt_begin_reversed, std::ref(res[t]));
  }

  for (std::uint64_t t = 0; t < n_threads; ++t) threads[t]->join();
  for (std::uint64_t t = 0; t < n_threads; ++t) delete threads[t];
  delete[] threads;

  mid_block_reader->stop();
  delete mid_block_reader;

  result = res;
}

}  // namespace psascan_private

#endif  // __SRC_PSASCAN_SRC_EM_COMPUTE_INITIAL_RANKS_HPP_INCLUDED
