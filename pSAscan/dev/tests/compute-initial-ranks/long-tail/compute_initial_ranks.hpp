/**
 * @file    src/psascan_src/compute_initial_ranks.hpp
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

#ifndef __SRC_PSASCAN_SRC_COMPUTE_INITIAL_RANKS_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_COMPUTE_INITIAL_RANKS_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>
#include <thread>
#include <omp.h>

#include "approx_rank.hpp"
#include "space_efficient_isa.hpp"
#include "utils.hpp"
#include "io/background_block_reader.hpp"
#include "io/background_chunk_reader.hpp"
#include "io/multifile_bit_stream_reader.hpp"


namespace psascan_private {


//=============================================================================
// Determine the lexicographical order between pat[0..pat_length) and
// text[block_suf_beg..text_length). We already know that they share a
// common prefix of length lcp. The task is to finish the comparison.
// We can (and will) only read symbols from text[block_suf_beg..block_end)
// because: (1) other symbols are not in fact available in text array, and
// (2) if we indeed reached text[block_end] we can resolve the comparison
// using the gt bitvector, to which we have access via the provided reader.
//
// XXX flaky description, because we completely ignore that fact, that
//     pat is somehow the suffix of text. Thus function, should not thus
//     receive a pointer to pat separatelly, but should acknowledge that
//     the pattern is a suffix of text. Then, the usage of gt bitvecotr
//     will be more clear.
//=============================================================================
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

  // Continue matching until we reach the
  // end of block or the end of the pattern.
  while (block_suf_beg + lcp < block_end &&
      lcp < pat_length &&
      text[block_suf_beg + lcp] == pat[lcp])
    ++lcp;

  // Return the answer depending on the lcp.
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
    typedef space_efficient_isa<rank_type, block_offset_type, 1L> isa_type;
#else
    typedef approx_rank<8L> rank_type;
    typedef space_efficient_isa<rank_type, block_offset_type, 8L> isa_type;
#endif

    rank_type *pbwt_rank =
      new rank_type(block_pbwt, block_length);
    isa_type *block_isa =
      new isa_type(block_psa, block, pbwt_rank, block_length, i0);

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
          std::uint64_t suf_len =
            block_length - (std::uint64_t)block_psa[mid];
          if (gt_reader->access(text_length -
                (stream_block_beg + suf_len)))
            left = mid + 1;
          else right = mid; 
        } else {
          std::uint64_t j =
            (std::uint64_t)block_psa[mid] + stream_block_size;
          if (block_isa->query(j) < prev_rank)
            left = mid + 1;
          else right = mid;
        }
      }

      result[t] = left;
      prev_rank = result[t];
    }

    delete pbwt_rank;
    delete block_isa;
    delete gt_reader;
  } else {
    for (std::uint64_t t = 0; t < n_threads; ++t)
      result[t] = ranges[t].first;
  }
}

//=============================================================================
// Return the result of the lexicographical comparison between
// two suffixes of text: text[block_suf_beg..text_length) and
// text[pat_beg..text_length) assuming that:
// * 0 <= lcp < pat_length <= text_length - pat_beg and the
//   two suffixes share a common prefix of length lcp.
// * 0 <= block_beg <= block_suf_beg < block_end <=
//   tail_beg <= pat_beg.
// * only substrings text[pat_beg+lcp..pat_beg+pat_length)
//   and text[block_suf_beg..block_end) are stored in RAM
//   and can be accessed, respectively, via variables
//   text[block_suf_beg..block_end) and pat[lcp..pat_length).
//   Any accesses to these variables outside the specified
//   ranges are undefined.
//
// Other information about the order between suffixes of the
// text is provided in the form of the "gt" bitvector (see the
// description of compute_tail_ranks below).
//
// The function is not required to resolve the comparison: the
// return value of -1 indicates that text[pat_length..text_length)
// is smaller than the other suffix, +1 indicates the opposite.
// 0 indicates that based on all information available the
// determination could not be made. Note that the return value
// of 0 implies that text[pat_beg..pat_beg+pat_length) is a prefix
// of suffix text[block_suf_beg..text_length), but the opposite is
// not true, i.e., the function could return != 0 even if the
// substring is a prefix of text[block_suf_beg..text_length),
// because the determination could be made on the information
// provided by the gt bitvector or the boundary conditions
// (e.g., if pat_beg + pat_length = text_length).
//
// NOTE: there are no constraints on the provided "lcp" value on
// the entry to the function, in particular it is even possible
// that block_suf_beg + lcp > block_end.
//=============================================================================
inline int lcp_compare(
    const std::uint8_t *text,
    const std::uint8_t *pat,
    background_block_reader *mid_block_reader,
    multifile_bit_stream_reader &gt_reader,
    std::uint64_t text_length,
    std::uint64_t block_end,
    std::uint64_t block_suf_beg,
    std::uint64_t pat_beg,
    std::uint64_t pat_length,
    std::uint64_t tail_begin,
    std::uint64_t &lcp) {

  // Continue the comparison using the
  // symbols explicitly available in RAM.
  while (block_suf_beg + lcp < block_end &&
      lcp < pat_length &&
      text[block_suf_beg + lcp] == pat[lcp])
    ++lcp;

  // At this point, one of the three conditions above
  // failed. Consider the case if it was the third,
  // i.e., we found a mismatch.
  if (block_suf_beg + lcp < block_end &&
      lcp < pat_length) {
    if (pat[lcp] > text[block_suf_beg + lcp])
      return +1;
    else return -1;
  }

  // Either we reached the end of block (or possibly, due to large
  // initial lcp value, we are already past the end of block), or
  // we reached the end of pattern. Consider the first case but only
  // if block_end < tail_begin, otherwise we can use gt to decide.
  if (block_suf_beg + lcp >= block_end &&
      block_end < tail_begin && lcp < pat_length) {

    // To continue the comparison, we need to access symbols in the
    // range text[block_end..tail_begin). Wait until enough symbols
    // are available.
    mid_block_reader->wait(
        std::min(tail_begin, block_suf_beg + pat_length) - block_end);

    // Now continue the comparison using
    // the symbols from the mid block.
    const std::uint8_t *text2 = mid_block_reader->m_data - block_end;
    while (block_suf_beg + lcp < tail_begin &&
        lcp < pat_length &&
        text2[block_suf_beg + lcp] == pat[lcp])
      ++lcp;

    // If the comparison ended with the
    // mismatch, we know the final answer
    // (either -1 or +1).
    if (block_suf_beg + lcp < tail_begin &&
        lcp < pat_length) {

      // Invariant: pat[lcp] != text2[block_suf_beg + lcp].
      if (pat[lcp] < text2[block_suf_beg + lcp])
        return -1;
      else return +1;
    }
  }

  // No more information can be gained
  // using symbol comparisons.
  if (block_suf_beg + lcp >= tail_begin) {

    // Use gt to resolve comparison.
    // Note that gt is stored reversed.
    std::uint64_t gt_pos =
      pat_beg + (tail_begin - block_suf_beg);
    if (gt_reader.access(text_length - gt_pos))
      return +1;
    else return -1;

  } else {

    // The only remaining option is that
    // lcp == pat_length. In this case we
    // can resolve the comparison only due
    // to boundary case, i.e., if we have
    // pat_beg + pat_length == text_length.
    if (pat_beg + pat_length == text_length)
      return -1;
    else return 0;
  }
}

//=============================================================================
// Auxiliary function. See the
// description of compute_tail_ranks.
//=============================================================================
template<typename block_offset_type>
void refine_range(
    const std::uint8_t *block,
    const std::uint8_t *pat,
    const block_offset_type *block_psa,
    background_block_reader *mid_block_reader,
    multifile_bit_stream_reader &gt_reader,
    std::uint64_t block_beg,
    std::uint64_t block_end,
    std::uint64_t pat_beg,
    std::uint64_t tail_begin,
    std::uint64_t text_length,
    std::uint64_t left,
    std::uint64_t right,
    std::uint64_t range_lcp,
    std::uint64_t pat_length,
    std::uint64_t &newleft,
    std::uint64_t &newright) {

  // Initialize the pointer to text. Due to the
  // construction of this function we will never
  // access symbols in text which are not in RAM.
  const std::uint8_t *text = block - block_beg;

  // Compute the tuning parameters
  // for string binary search.
#ifdef EM_STARTING_POS_MODULE_DEBUG_MODE
  std::uint64_t min_discrepancy = utils::random_int64(0L, 10L);
  std::uint64_t balancing_factor = utils::random_int64(1L, 10L);
#else
  static const std::uint64_t min_discrepancy = (1 << 16);
  static const std::uint64_t balancing_factor = 64;
#endif

  // Compute the lower bound of the new range.
  // Invariant: the answer is always in (low, high].
  std::int64_t low = (std::int64_t)left - 1;
  std::int64_t high = right;
  std::uint64_t llcp = range_lcp;
  std::uint64_t rlcp = range_lcp;
  while (low + 1 != high) {

    // Compute the pivot. We use the skewed
    // string binary search. The pivot
    // satisfies: low < mid < high.
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

    // Compare the pivot to the pattern. See the
    // lcp_compare for the exact explanation.
    std::uint64_t lcp = std::min(llcp, rlcp);
    std::uint64_t block_suf_beg =
      block_beg + (std::uint64_t)block_psa[mid];
    if (lcp_compare(text, pat, mid_block_reader, gt_reader,
          text_length, block_end, block_suf_beg, pat_beg,
          pat_length, tail_begin, lcp) <= 0) {
      high = mid;
      rlcp = lcp;
    } else {
      low = mid;
      llcp = lcp;
    }
  }

  // Store the lower bound
  // of the new range.
  newleft = high;

  // rlcp < pat_length means no suffix in block_psa
  // has pat[0..pat_length) as a prefix or the ones
  // that do were already classified as larger. We
  // only need to do second binary search if rlcp
  // >= pat_length.
  if (rlcp >= pat_length) {
    high = right;
    rlcp = range_lcp;

    // Compute the upper bound of the new range.
    // Invariant: the answer is always in (low, high].
    while (low + 1 != high) {

      // Compute the pivot. We use the skewed
      // string binary search. The pivot
      // satisfies: low < mid < high.
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

      // Compare the pivot to the pattern. See the
      // lcp_compare for the exact explanation.
      std::uint64_t lcp = std::min(llcp, rlcp);
      std::uint64_t block_suf_beg =
        block_beg + (std::uint64_t)block_psa[mid];
      if (lcp_compare(text, pat, mid_block_reader, gt_reader,
            text_length, block_end, block_suf_beg, pat_beg,
            pat_length, tail_begin, lcp) < 0) {
        high = mid;
        rlcp = lcp;
      } else {
        low = mid;
        llcp = lcp;
      }
    }
  }

  // Store the upper bound
  // of the new range.
  newright = high;
}

//=============================================================================
// Auxiliary function. See the
// description of compute_tail_ranks.
//=============================================================================
template<typename block_offset_type>
void compute_pattern_rank(
    const std::uint8_t *block,
    const block_offset_type *block_psa,
    const multifile *tail_gt_begin_reversed,
    background_block_reader *mid_block_reader,
    std::string text_filename,
    std::uint64_t block_begin,
    std::uint64_t block_end,
    std::uint64_t pat_begin,
    std::uint64_t text_length,
    std::uint64_t tail_begin,
    std::uint64_t &result) {

  // Handle special case.
  if (pat_begin == text_length) {
    result = 0;
    return;
  }

  // Compute basic parameters.
  std::uint64_t block_size = block_end - block_begin;
  std::uint64_t mid_block_size = tail_begin - block_end;
  std::uint64_t max_pat_symbols_needed = std::min(
      text_length - pat_begin, block_size + mid_block_size);

  // Create the reader of
  // tail_gt_begin_reversed.
  multifile_bit_stream_reader
    gt_reader(tail_gt_begin_reversed);

  // The computation is done with the binary search over the array
  // 'psa'. However, a comparison of the pattern with a single suffix
  // from psa could require accessing nearly whole pattern. We cannot
  // afford to store it, so in the worst case we would have to stream
  // the whole pattern many times. To prevent this, we perform the
  // computation in chunks. Each time we read a small chunk of the
  // pattern into RAM, we refine the range of suffixes in the suffix
  // array using symbols in the current chunk. This way we use little
  // RAM and are guaranteed to only read the pattern once.
#ifdef EM_STARTING_POS_MODULE_DEBUG_MODE
  std::uint64_t chunk_length =
    utils::random_int64(1, 10);
  background_chunk_reader *chunk_reader =
    new background_chunk_reader(text_filename, pat_begin,
        pat_begin + max_pat_symbols_needed, chunk_length);
#else
  background_chunk_reader *chunk_reader =
    new background_chunk_reader(text_filename, pat_begin,
        pat_begin + max_pat_symbols_needed);
#endif

  // Find the answer using binary search. Invariant:
  // the answer (i.e., final rank) is always in the
  // range [left, right]. Each step of the loop below
  // processes a single chunk and refines the range
  // untils it contains one item.
  std::uint64_t range_lcp = 0;
  std::uint64_t left = 0;
  std::uint64_t right = block_size;
  while (left != right) {

    // Read the next chunk of the pattern.
    std::uint64_t this_chunk_length = std::min(
        max_pat_symbols_needed - range_lcp,
        chunk_reader->get_chunk_size());
    std::uint64_t pat_length = range_lcp + this_chunk_length;
    chunk_reader->wait(pat_begin + pat_length);

    // Compute new range boundaries.
    // Invariant: all suffixes in psa[left..right]
    // share a common prefix of length range_lcp
    // with the pattern.
    // Invariant:
    // chunk_reader->m_chunk[0..this_chunk_length)
    // == text[pat_begin+range_lcp..pat_begin+pat_length).
    std::uint64_t newleft = 0;
    std::uint64_t newright = 0;
    refine_range(block, chunk_reader->m_chunk - range_lcp, block_psa,
        mid_block_reader, gt_reader, block_begin, block_end,
        pat_begin, tail_begin, text_length, left, right, range_lcp,
        pat_length, newleft, newright);

    // Update range boundaries.
    left = newleft;
    right = newright;
    range_lcp = pat_length;
  }

  // Store the result.
  result = left;

  // Clean up.
  delete chunk_reader;
}

//=============================================================================
// The aim of this function is to compute the rank (i.e., the number
// of smaller suffixes) among the suffixes of text starting inside
// text[block_beg..block_end) for each of the strings in the sparse set
// of suffixes of text called here "patterns". All patterns start inside
// text[tail_begin..text_length). It is guaranteed that tail_begin >=
// block_end, but not that tail_begin == block_end. The symbols inside
// block are provided explicitly in the array block[0..block_size). The
// text in between the block and the tail (text[block_end..tail_begin)
// is referred to as "mid block" and is not provided explicitly. The
// number and the exact starting positions of patterns are computed inside
// the function. The distances between the starting positions of patterns
// are equal, but the function does not take advantage of this fact. The
// rank of each pattern is computed independently by a separate thread.
//
// To determine ties in case of very long common prefixes between suffixes
// starting inside the block and the pattern, the function is provided
// with a bitvector (called here "gt", for technical reason it is reversed)
// that tells, for every position i in the range [tail_begin..text_length)
// whether the suffix text[i..text_length) is lexicographically greater
// than the suffix text[tail_begin..text_length). This way we never have
// to compare more than tail_begin - block_beg symbols to resolve the
// comparison between any pattern and any suffix starting in text[block_beg
// ..block_end). This prevents the worst case occurring on artificial
// inputs, in practice the bitvector is almost never used.
//=============================================================================
template<typename block_offset_type>
void compute_tail_ranks(
    const std::uint8_t *block,
    const block_offset_type *block_psa,
    const multifile *tail_gt_begin_reversed,
    std::string text_filename,
    std::uint64_t block_beg,
    std::uint64_t block_end,
    std::uint64_t text_length,
    std::uint64_t tail_begin,
    std::uint64_t max_threads,
    std::vector<std::uint64_t> &result) {

  // Compute some initial parameters.
  std::uint64_t tail_length = text_length - tail_begin;
  std::uint64_t mid_block_beg = block_end;
  std::uint64_t mid_block_end = tail_begin;
  std::uint64_t mid_block_size = mid_block_end - mid_block_beg;
  std::uint64_t stream_max_block_size =
    (tail_length + max_threads - 1) / max_threads;
  std::uint64_t n_patterns =
    (tail_length + stream_max_block_size - 1) / stream_max_block_size;

  // Start reading the mid block in the
  // background. We only need it in rare cases.
  background_block_reader *mid_block_reader =
    new background_block_reader(text_filename,
        mid_block_beg, mid_block_size);

  // Compute the ranks.
  std::vector<std::uint64_t> res(n_patterns);
  {

#ifdef _OPENMP

    // Parallel version.
    #pragma omp parallel num_threads(n_patterns)
    {

      // Compute pattern id and starting position.
      std::uint64_t pattern_id = omp_get_thread_num();
      std::uint64_t pattern_beg = tail_begin +
        pattern_id * stream_max_block_size;

      // Compute the rank.
      compute_pattern_rank<block_offset_type>(
          block, block_psa, tail_gt_begin_reversed,
          mid_block_reader, text_filename, block_beg,
          block_end, pattern_beg, text_length,
          tail_begin, res[pattern_id]);
    }
#else

    // Sequential version.
    for (std::uint64_t pattern_id = 0;
        pattern_id < n_patterns; ++pattern_id) {

      // Compute pattern starting position.
      std::uint64_t pattern_beg = tail_begin +
        pattern_id * stream_max_block_size;

      // Compute the rank.
      compute_pattern_rank<block_offset_type>(
          block, block_psa, tail_gt_begin_reversed,
          mid_block_reader, text_filename, block_beg,
          block_end, pattern_beg, text_length,
          tail_begin, res[pattern_id]);
    }

#endif  // _OPENMP
  }

  // Stop and delete the
  // mid block reader.
  mid_block_reader->stop();
  delete mid_block_reader;

  // Store the result.
  result = res;
}

}  // namespace psascan_private

#endif  // __SRC_PSASCAN_SRC_COMPUTE_INITIAL_RANKS_HPP_INCLUDED
