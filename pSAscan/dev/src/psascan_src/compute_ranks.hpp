/**
 * @file    src/psascan_src/compute_ranks.hpp
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

#ifndef __SRC_PSASCAN_SRC_COMPUTE_RANKS_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_COMPUTE_RANKS_HPP_INCLUDED

#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>
#include <omp.h>

#include "approx_rank.hpp"
#include "space_efficient_isa.hpp"
#include "utils.hpp"
#include "io/background_block_reader.hpp"
#include "io/background_chunk_reader.hpp"
#include "io/multifile_bit_stream_reader.hpp"


namespace psascan_private {


//=============================================================================
// Let 0 <= init_lcp <= pattern_length, 0 <= block_beg <= block_suf_begin
// < block_end <= pattern_begin < tail_end <= text_length, pattern_begin
// + pattern_length <= tail_end. Let x = text[pattern_begin..text_length)
// and y = text[block_suf_begin..text_length). Assume that x and y share
// a common prefix of length init_lcp. To explain the aim of this function,
// observe that at least one of the following three cases holds: (1) x > y,
// (2) x < y, (3) x and y share a common prefix of length >= pattern_length.
//
// The aim of this function is to determine which of the cases hold and
// return +1 for case (1), -1 for case (2) and 0 for case (3). If more
// than one case holds, any integer indicating one of the holding cases
// is a correct return value. If the function returns 0, ret_lcp has to be
// set to pattern_length. In other cases, the ret_lcp must be a true
// lower bound for lcp(x, y) not smaller than init_lcp. No other requirement
// is placed on ret_lcp in those cases.
//
// The function is not provided with the full text, nor the filename with
// the text. It is given the following information about the text and the
// ordering of its suffixes:
// * text[block_begin..block_end) and text[pattern_begin+init_lcp..
//   pattern_begin+pattern_length) are the accessible substrings of text.
//   They are available, respectively, via text[block_beg..block_end) and
//   pat[init_lcp..pattern_length) variables. Any accesses to these
//   variables outside the specified ranges are undefined.
// * gt bitvector telling for every i in the range (block_end..tail_end]
//   whether the suffix text[i..text_length) is greater (the bit is 1)
//   than text[block_end..text_length).
// The above information is sufficient to compute the answer.
//
// Note that gt bitvector is reversed (for technical reasons) and shifted
// (to provide convenient indexing), i.e., the bit of gt corresponding
// position i in (block_end..tail_end] (see definition above) is accessed
// as gt[text_length - i].
//=============================================================================
template<typename char_type>
inline std::int32_t lcp_compare(
    const std::uint64_t text_length,
    const std::uint64_t block_end,
    const std::uint64_t block_suf_begin,
    const std::uint64_t pattern_begin,
    const std::uint64_t pattern_length,
    const std::uint64_t init_lcp,
    const char_type * const text,
    const char_type * const pat,
    multifile_bit_stream_reader * const gt_reader,
    std::uint64_t * const ret_lcp) {

  // Continue matching until we reach the
  // end of block or the end of the pattern.
  std::uint64_t new_lcp = init_lcp;
  while (block_suf_begin + new_lcp < block_end &&
      new_lcp < pattern_length &&
      text[block_suf_begin + new_lcp] == pat[new_lcp])
    ++new_lcp;

  // Store the new lcp.
  *ret_lcp = new_lcp;

  // Return the answer depending on new_lcp.
  if (block_suf_begin + new_lcp >= block_end) {
    const std::uint64_t block_suf_len = block_end - block_suf_begin;
    const std::uint64_t pat_ptr = pattern_begin + block_suf_len;

    // We reached the end of the block. The answer
    // can be computed from the gt bitvector.
    if (gt_reader->access(text_length - pat_ptr))
      return 1;
    else return -1;
  } else if (new_lcp == pattern_length) {

    // We have reached the end of pattern.
    // Return -1 but only if we also reached
    // the end of text. Otherwise, return 0.
    if (pattern_begin + pattern_length == text_length)
      return -1;
    else return 0;
  } else {

    // We found a mismatch. Return the value
    // depending on the symbol comparison.
    if (pat[new_lcp] > text[block_suf_begin + new_lcp])
      return 1;
    else return -1;
  } 
}

//=============================================================================
// Let 0 <= block_begin < block_end <= text_length and define block_size
// = block_end - block_begin. Further, let block_psa[0..block_size) be
// an array that contains a permutation of the set {0, .., block_size - 1}
// that describes the ascending lexicographical order of suffixes of text
// starting inside block (and extending beyond the end of the block).
//
// Assume that block_end <= pattern_begin < tail_end <= text_length and
// that 0 <= range_lcp <= pattern_length <= tail_end - pattern_begin.
// Further, assume that integers 0 <= left <= right <= block_size satisfy
// the following conditions:
// * any suffix of text starting at position block_begin + block_psa[i],
//   for 0 <= i < left is smaller than text[pattern_begin..text_length).
// * any suffix of text starting at position block_begin + block_psa[i],
//   for right <= i < block_size is larger than text[pattern_begin..
//   text_length).
// * any suffix of text starting at position block_begin + block_psa[i]
//   for left <= i < right has text[pattern_begin..pattern_begin +
//   range_lcp) as a prefix.
// Note that we do not assume the range [left, right) to be the maximal
// range satisfying the above conditions.
//
// The aim of this function is to compute integers newleft, newright
// that satisfy left <= newleft <= newright <= right and that satisfy
// all three conditions listed above (with left and right replaced with,
// respectively, newleft and newright), except with range_lcp in the
// third condition replaced by pattern_length. Note again that we do
// not require the range [newleft, newright) to the be the maximal range
// satisfying the specified conditions, hence the result of this function
// is not uniquely defined.
//
// Note that the function is not provided with the full text, nor the
// filename with the text. It is given the following information about
// the text and the ordering of its suffixes:
// * text[block_begin..block_end) and text[pattern_begin+range_lcp..
//   pattern_begin+pattern_length) are the accessible substrings of
//   text. They are available, respectively, via block[0..block_size)
//   and pat[range_lcp..pattern_length) variables. Any accesses to
//   these variables outside the specified ranges are undefined.
// * block_psa[0..block_size) as define above.
// * gt bitvector telling for every i in range (block_end..tail_end]
//   whether the suffix text[i..text_length) is greater (then bit is
//   1) than text[block_end..text_length).
// The above information is sufficient to compute the answer.
//
// Note that gt bitvector is reversed (for technical reasons) and shifted
// (to provide convenient indexing), i.e., the bit of gt corresponding
// position i in (block_end..tail_end] (see definition above) is accessed
// as gt[text_length - i].
//=============================================================================
template<
  typename char_type,
  typename block_offset_type>
void refine_range(
    const std::uint64_t block_begin,
    const std::uint64_t block_end,
    const std::uint64_t pattern_begin,
    const std::uint64_t pattern_length,
    const std::uint64_t text_length,
    const std::uint64_t left,
    const std::uint64_t right,
    const std::uint64_t range_lcp,
    const char_type * const block,
    const char_type * const pat,
    const block_offset_type * const block_psa,
    multifile_bit_stream_reader * const gt_reader,
    std::uint64_t * const newleft,
    std::uint64_t * const newright) {

  // Initialize the pointer to text. Due to the
  // construction of this function we will never
  // access symbols in text which are not in RAM.
  const char_type * const text = block - block_begin;

  // Compute the tuning parameters
  // for string binary search.
#ifdef EM_STARTING_POS_MODULE_DEBUG_MODE
  std::uint64_t min_discrepancy = utils::random_int64(0, 10);
  std::uint64_t balancing_factor = utils::random_int64(1, 10);
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
      const std::uint64_t d = rlcp - llcp;
      const std::uint64_t logd = utils::log2ceil(d);
      mid = low + 1 +
        ((high - low - 1) * balancing_factor * logd) /
        (d + balancing_factor * logd);
    } else if (rlcp + min_discrepancy < llcp) {
      const std::uint64_t d = llcp - rlcp;
      const std::uint64_t logd = utils::log2ceil(d);
      mid = high - 1 -
        ((high - low - 1) * balancing_factor * logd) /
        (d + balancing_factor * logd);
    } else mid = (low + high) / 2;

    // Compare the pivot to the pattern.
    std::uint64_t new_lcp = 0;
    const std::uint64_t init_lcp = std::min(llcp, rlcp);
    const std::uint64_t block_suf_begin =
      block_begin + (std::uint64_t)block_psa[mid];
    if (lcp_compare(text_length, block_end,
          block_suf_begin, pattern_begin, pattern_length,
          init_lcp, text, pat, gt_reader, &new_lcp) <= 0) {
      high = mid;
      rlcp = new_lcp;
    } else {
      low = mid;
      llcp = new_lcp;
    }
  }

  // Store the lower bound
  // of the new range.
  *newleft = high;

  // rlcp < pattern_length means no suffix in block_psa
  // has pat[0..pattern_length) as a prefix or the ones
  // that do were already classified as larger. We
  // only need to do second binary search if rlcp
  // >= pattern_length.
  if (rlcp >= pattern_length) {
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
        const std::uint64_t d = rlcp - llcp;
        const std::uint64_t logd = utils::log2ceil(d);
        mid = low + 1 +
          ((high - low - 1) * balancing_factor * logd) /
          (d + balancing_factor * logd);
      } else if (rlcp + min_discrepancy < llcp) {
        const std::uint64_t d = llcp - rlcp;
        const std::uint64_t logd = utils::log2ceil(d);
        mid = high - 1 -
          ((high - low - 1) * balancing_factor * logd) /
          (d + balancing_factor * logd);
      } else mid = (low + high) / 2;

      // Compare the pivot to the pattern.
      std::uint64_t new_lcp = 0;
      const std::uint64_t init_lcp = std::min(llcp, rlcp);
      const std::uint64_t block_suf_begin =
        block_begin + (std::uint64_t)block_psa[mid];
      if (lcp_compare(text_length, block_end,
            block_suf_begin, pattern_begin, pattern_length,
            init_lcp, text, pat, gt_reader, &new_lcp) < 0) {
        high = mid;
        rlcp = new_lcp;
      } else {
        low = mid;
        llcp = new_lcp;
      }
    }
  }

  // Store the upper bound
  // of the new range.
  *newright = high;
}

//=============================================================================
// Let 0 <= block_begin < block_end <= text_length and define block_size
// = block_end - block_begin. Further, let block_psa[0..block_size) be
// an array that contains a permutation of the set {0, .., block_size - 1}
// that describes the ascending lexicographical order of suffixes of text
// starting inside block (and extending beyond the end of the block).
//
// Assume block_end <= pattern_begin < tail_end, max_pat_length > 0, and
// 0 <= pattern_begin <= tail_end - max_pat_length. This function aims
// to compute the integers left and right such that 0 <= left <= right
// <= block_size satisfying the following conditions:
// * any suffix of text starting at position block_begin + block_psa[i],
//   for 0 <= i < left is smaller than text[pattern_begin..text_length).
// * any suffix of text starting at position block_begin + block_psa[i],
//   for right <= i < block_size is larger than text[pattern_begin..
//   text_length).
// * any suffix of text starting at position block_begin + block_psa[i]
//   for left <= i < right has text[pattern_begin..pattern_begin +
//   max_pat_length) as a prefix.
// Note that we do not require the range [left, right) to be the maximal
// range satisfying the conditions above, thus the result of this function
// is not uniquely defined.
//
// The function is given the name of the file containing the text.
// The text itself is in principle sufficient to compute the answer.
// However, to perform the computation efficiently, the function receives
// the following information about the text and ordering of its suffixes:
// * block[0..block_size) = text[block_begin..block_end).
// * block_psa[0..block_size), as described above.
// * gt bitvector telling for every i in range (block_end..tail_end]
//   whether the suffix text[i..text_length) is greater (then bit is
//   1) than text[block_end..text_length).
//
// Note that gt bitvector is reversed (for technical reasons) and shifted
// (to provide convenient indexing), i.e., the bit of gt corresponding
// position i in (block_end..tail_end] (see definition above) is accessed
// as gt[text_length - i].
//=============================================================================
template<
  typename char_type,
  typename block_offset_type>
std::pair<std::uint64_t, std::uint64_t>
compute_range(
    const std::uint64_t block_begin,
    const std::uint64_t block_end,
    const std::uint64_t pattern_begin,
    const std::uint64_t text_length,
    const std::uint64_t max_pat_length,
    const char_type * const block,
    const block_offset_type * const block_psa,
    const multifile * const tail_gt_begin_reversed,
    const std::string text_filename) {

  // Handle special case.
  if (pattern_begin == text_length)
    return std::make_pair(0, 0);

  // Compute initial parameters.
  const std::uint64_t block_size =
    block_end - block_begin;

  // Create the reader of the gt bitvector.
  multifile_bit_stream_reader *gt_reader =
    new multifile_bit_stream_reader(tail_gt_begin_reversed);

  // Read text[pattern_begin..pattern_begin+max_pat_length) in chunks.
  // The computation is done with the binary search over the block_psa
  // array. However, a comparison of the pattern with a single suffix
  // from block_psa could require accessing nearly whole pattern. We
  // cannot afford to store it, so in the worst case we would have to
  // stream the whole pattern many times. To prevent this, we perform
  // the computation in chunks. Each time we read a small chunk of the
  // pattern into RAM, we refine the range of suffixes in the suffix
  // array using symbols in the current chunk. This way we use little
  // extra RAM and are we guaranteed to only read the pattern once.
  // Further, in most cases, refining a range based on some initial
  // chunk already reduced the range to a single element and the
  // computation can be finished.
#ifdef EM_STARTING_POS_MODULE_DEBUG_MODE
  std::uint64_t chunk_length =
    utils::random_int64(1, 10);
  background_chunk_reader * const chunk_reader =
    new background_chunk_reader(text_filename,
        pattern_begin, pattern_begin + max_pat_length,
        chunk_length);
#else

  // Use default chunk length.
  background_chunk_reader * const chunk_reader =
    new background_chunk_reader(text_filename,
        pattern_begin, pattern_begin + max_pat_length);
#endif

  // Initialize the output range to [0, block_size).
  // At each step we refine the range until it
  // becomes empty or the common prefix between
  // text[pattern_begin..text_length) and all
  // suffixes in the range [left, right) is equal
  // to max_pat_length.
  std::uint64_t range_lcp = 0;
  std::uint64_t left = 0;
  std::uint64_t right = block_size;
  while (left != right && range_lcp < max_pat_length) {

    // Read the next chunk of the pattern.
    const std::uint64_t this_chunk_length = std::min(
        max_pat_length - range_lcp,
        chunk_reader->get_chunk_size());
    const std::uint64_t pattern_length =
      range_lcp + this_chunk_length;
    chunk_reader->wait(pattern_begin + pattern_length);

    // Refine the range using the next chunk. Invariant:
    // all suffixes in block_psa[left..right) share a
    // common prefix of length range_lcp with the pattern.
    // Invariant:
    // chunk_reader->m_chunk[0..this_chunk_length) == text[
    // pattern_begin+range_lcp..pattern_begin+pattern_length).
    std::uint64_t newleft = 0;
    std::uint64_t newright = 0;
    refine_range(block_begin, block_end, pattern_begin,
        pattern_length, text_length, left, right, range_lcp,
        block, chunk_reader->m_chunk - range_lcp, block_psa,
        gt_reader, &newleft, &newright);

    // Update range boundaries.
    left = newleft;
    right = newright;
    range_lcp = pattern_length;
  }

  // Clean up.
  delete chunk_reader;
  delete gt_reader;

  // Return the result.
  return std::make_pair(left, right);
}

//=============================================================================
// The aim of this function is to compute the rank (i.e., the number
// of smaller suffixes) among the suffixes of text starting inside
// text[block_begin..block_end) (0 <= block_begin < block_end <= text_length)
// for each of the strings in the sparse set of suffixes of text called
// patterns. The set of patterns is defined as the set of all suffixes
// of text that start at positions i such that block_end <= i < tail_end
// and (i - block_end) is a multiple of patterns_dist > 0.
//
// The function is given the name of the file containing the text.
// The text itself is in principle sufficient to compute all the ranks.
// However, to perform the computation efficiently, the function receives
// the following information about the text and ordering of its suffixes.
// Let block_size = block_end - block_begin. The additional info is:
// * block[0..block_size) = text[block_begin..block_end).
// * block_psa[0..block_size) contains permutation of integer set
//   {0, .., block_size-1} that describes the ascending lexicographical
//   order of suffixes of text starting inside block (and extending
//   beyond the end of the block).
// * gt bitvector telling for every i in range (block_end..tail_end]
//   whether the suffix text[i..text_length) is greater (then bit is
//   1) than text[block_end..text_length).
// * rank_tail_end contains the number of suffixes of text starting
//   inside text[block_begin..block_end) that are smaller than the suffix
//   text[tail_end..text_length).
//
// Furthermore, the function is given auxiliary information (which
// can be quickly computed from above, but is provided for efficiency):
// * block_pbwt[0..block_size) which contains the Burrows-Wheeler
//   transform of block. Namely, for any i such that block_psa[i] != 0,
//   block_pbwt[i] = block[block_psa[i] - 1). If block_psa[i] == 0,
//   then block_pbwt[i] = 0.
// * i0 is the position such that block_psa[i0] = 0.
//
// Note that gt bitvector is reversed (for technical reasons) and shifted
// (to provide convenient indexing), i.e., the bit of gt corresponding
// position i in (block_end..tail_end] (see definition above) is accessed
// as gt[text_length - i].
//=============================================================================
template<
  typename char_type,
  typename block_offset_type>
void compute_ranks(
    const std::uint64_t block_begin,
    const std::uint64_t block_end,
    const std::uint64_t text_length,
    const std::uint64_t patterns_dist,
    const std::uint64_t tail_end,
    const std::uint64_t rank_tail_end,
    const std::uint64_t i0,
    const char_type * const block,
    const char_type * const block_pbwt,
    const block_offset_type * const block_psa,
    const multifile * const tail_gt_begin_reversed,
    const std::string text_filename,
    std::vector<std::uint64_t> &result) {

  // Compute the number of patterns.
  const std::uint64_t block_length = block_end - block_begin;
  const std::uint64_t tail_length = tail_end - block_end;
  const std::uint64_t n_patterns =
    (tail_length + patterns_dist - 1) / patterns_dist;

  // Allocate initial ranges.
  std::vector<std::pair<std::uint64_t, std::uint64_t> >
    ranges(n_patterns);

  // Compute initial ranges. See the
  // description of compute_range.
  {

#ifdef _OPENMP

    // Parallel version.
    #pragma omp parallel num_threads(n_patterns)
    {

      // Compute pattern id and boundaries.
      const std::uint64_t pattern_id = omp_get_thread_num();
      const std::uint64_t pattern_begin =
        block_end + pattern_id * patterns_dist;
      const std::uint64_t pattern_prefix_end =
        std::min(tail_end, pattern_begin + patterns_dist);
      const std::uint64_t pattern_prefix_length =
        pattern_prefix_end - pattern_begin;

      // Compute the range.
      ranges[pattern_id] = compute_range(
          block_begin, block_end, pattern_begin, text_length,
          pattern_prefix_length, block, block_psa,
          tail_gt_begin_reversed, text_filename);
    }

#else

    // Sequential version.
    for (std::uint64_t pattern_id = 0;
        pattern_id < n_patterns; ++pattern_id) {

      // Compute pattern boundaries.
      const std::uint64_t pattern_begin =
        block_end + pattern_id * patterns_dist;
      const std::uint64_t pattern_prefix_end =
        std::min(tail_end, pattern_begin + patterns_dist);
      const std::uint64_t pattern_prefix_length =
        pattern_prefix_end - pattern_begin;

      // Compute the range.
      ranges[pattern_id] = compute_range(
          block_begin, block_end, pattern_begin, text_length,
          pattern_prefix_length, block, block_psa,
          tail_gt_begin_reversed, text_filename);
    }

#endif  // _OPENMP
  }

  // Allocate space for output.
  result.resize(n_patterns);

  // Determine if there is at
  // least one nontrivial range.
  bool nontrivial_range = false;
  for (std::uint64_t pattern_id = 0;
      pattern_id < n_patterns; ++pattern_id)
    if (ranges[pattern_id].first !=
        ranges[pattern_id].second)
      nontrivial_range = true;

  // If the nontrivial range was found
  // refine all ranges using the fact
  // that patterns are evenly spaced.
  if (nontrivial_range == true) {

    // Create the reader of the gt bitvector.
    multifile_bit_stream_reader * const gt_reader =
      new multifile_bit_stream_reader(tail_gt_begin_reversed);

    // Decide on the sampling rates for
    // auxiliary data structures.
#ifdef EM_STARTING_POS_MODULE_DEBUG_MODE
    typedef approx_rank<1> rank_type;
    typedef space_efficient_isa<rank_type,
            block_offset_type, 1> isa_type;
#else
    typedef approx_rank<8> rank_type;
    typedef space_efficient_isa<rank_type,
            block_offset_type, 8> isa_type;
#endif

    // Create data structures able to quickly
    // locate, for any p, a position i such
    // that block_psa[i] = p.
    rank_type * const pbwt_rank =
      new rank_type(block_pbwt, block_length);
    isa_type * const block_isa =
      new isa_type(block_psa, block, pbwt_rank, block_length, i0);

    // Process all patterns in right-to-left order.
    std::uint64_t prev_rank = rank_tail_end;
    for (std::uint64_t pattern_id_plus = n_patterns;
        pattern_id_plus > 0; --pattern_id_plus) {
      const std::uint64_t pattern_id = pattern_id_plus - 1;

      // Compute block boundaries.
      const std::uint64_t pattern_begin =
        block_end + pattern_id * patterns_dist;
      const std::uint64_t pattern_prefix_end =
        std::min(tail_end, pattern_begin + patterns_dist);
      const std::uint64_t pattern_prefix_length =
        pattern_prefix_end - pattern_begin;

      // Binary search for the final rank of the
      // current pattern.
      std::uint64_t left = ranges[pattern_id].first;
      std::uint64_t right = ranges[pattern_id].second;
      while (left != right) {

        // Invariant: all suffixes of text starting at
        // positions block_begin + block_psa[i] for
        // i in [left, right) have a common prefix
        // of length pattern_prefix_length with the
        // suffix text[pattern_begin..text_length).
        // Invariant: the answer (the final rank
        // of suffix text[pattern_begin..text_length)
        // is in range [left, right].
        // Choose pivot mid in [left, right).
        const std::uint64_t mid = (left + right) / 2;

        // Determine if text[block_begin+block_psa[mid]..text_length)
        // is smaller or larger than text[pattern_begin..text_length).
        if ((std::uint64_t)block_psa[mid] +
            pattern_prefix_length >= block_length) {

          // The comparison can be resolved
          // by looking up the gt bitvector.
          const std::uint64_t suf_len =
            block_length - (std::uint64_t)block_psa[mid];
          if (gt_reader->access(text_length -
                (pattern_begin + suf_len)))
            left = mid + 1;
          else right = mid; 
        } else {

          // The final comparison can be resolved by
          // looking up the answer for the pattern to
          // the right and taking advantage of the
          // fact that all patterns are evenly spaced.
          // This case requires looking up the inverse
          // suffix array of block_psa.
          const std::uint64_t j =
            (std::uint64_t)block_psa[mid] +
            pattern_prefix_length;
          if (block_isa->query(j) < prev_rank)
            left = mid + 1;
          else right = mid;
        }
      }

      // Store the final rank.
      result[pattern_id] = left;
      prev_rank = result[pattern_id];
    }

    // Clean up.
    delete pbwt_rank;
    delete block_isa;
    delete gt_reader;
  } else {

    // If all ranges are trivial,
    // we can simply copy the answers.
    for (std::uint64_t pattern_id = 0;
        pattern_id < n_patterns; ++pattern_id)
      result[pattern_id] = ranges[pattern_id].first;
  }
}

//=============================================================================
// Let 0 <= init_lcp <= pattern_length, 0 <= block_beg <= block_suf_begin
// < block_end <= tail_begin <= pattern_begin < text_length, pattern_begin
// + pattern_length <= text_length. Let x = text[pattern_begin..text_length)
// and y = text[block_suf_begin..text_length). Assume that x and y share
// a common prefix of length init_lcp. To explain the aim of this function,
// observe that at least one of the following three cases holds: (1) x > y,
// (2) x < y, (3) x and y share a common prefix of length >= pattern_length.
//
// The aim of this function is to determine which of the cases hold and
// return +1 for case (1), -1 for case (2) and 0 for case (3). If more
// than one case holds, any integer indicating one of the holding cases
// is a correct return value. If the function returns 0, ret_lcp has to be
// set to pattern_length. In other cases, the ret_lcp must be a true
// lower bound for lcp(x, y) not smaller than init_lcp. No other requirement
// is placed on ret_lcp in those cases.
//
// The function is not provided with the full text, nor the filename with
// the text. It is given the following information about the text and the
// ordering of its suffixes:
// * text[block_begin..block_end) and text[pattern_begin+init_lcp..
//   pattern_begin+pattern_length) are the accessible substrings of text.
//   They are available, respectively, via text[block_beg..block_end) and
//   pat[init_lcp..pattern_length) variables. Any accesses to these
//   variables outside the specified ranges are undefined.
// * gt bitvector telling for every i in range (tail_begin..text_length]
//   whether the suffix text[i..text_length) is greater (then bit is 1)
//   than text[tail_begin..text_length).
// * mid_block_reader allows accessing (after first calling the wait
//   method to check that enough text has been read into RAM) the
//   substring text[block_end..tail_begin).
// The above information is sufficient to compute the answer.
//
// Note that gt bitvector is reversed (for technical reasons) and shifted
// (to provide convenient indexing), i.e., the bit of gt corresponding
// position i in (tail_begin..text_length] (see definition above) is
// accessed as gt[text_length - i].
//=============================================================================
template<typename char_type>
inline std::int32_t lcp_compare(
    const std::uint64_t text_length,
    const std::uint64_t block_end,
    const std::uint64_t block_suf_begin,
    const std::uint64_t pattern_begin,
    const std::uint64_t pattern_length,
    const std::uint64_t tail_begin,
    const std::uint64_t init_lcp,
    const char_type * const text,
    const char_type * const pat,
    background_block_reader * const mid_block_reader,
    multifile_bit_stream_reader * const gt_reader,
    std::uint64_t * const ret_lcp) {

  // Continue the comparison using the
  // symbols explicitly available in RAM.
  std::uint64_t new_lcp = init_lcp;
  while (block_suf_begin + new_lcp < block_end &&
      new_lcp < pattern_length &&
      text[block_suf_begin + new_lcp] == pat[new_lcp])
    ++new_lcp;

  // Store the result.
  *ret_lcp = new_lcp;

  // At this point, one of the three conditions above
  // failed. Consider the case if it was the third,
  // i.e., we found a mismatch.
  if (block_suf_begin + new_lcp < block_end &&
      new_lcp < pattern_length) {
    if (pat[new_lcp] > text[block_suf_begin + new_lcp])
      return +1;
    else return -1;
  }

  // Either we reached the end of block (or possibly, due to large
  // initial lcp value, we are already past the end of block), or
  // we reached the end of pattern. Consider the first case but only
  // if block_end < tail_begin, otherwise we can use gt to decide.
  if (block_suf_begin + new_lcp >= block_end &&
      block_end < tail_begin &&
      new_lcp < pattern_length) {

    // To continue the comparison, we need to access symbols in the
    // range text[block_end..tail_begin). Wait until enough symbols
    // are available.
    const std::uint64_t mid_block_needed = std::min(tail_begin,
        block_suf_begin + pattern_length) - block_end;
    mid_block_reader->wait(mid_block_needed);

    // Now continue the comparison using the symbols from the mid block.
    // Invariant: mid_block_reader->m_data[0..mid_block_needed) ==
    // text[block_end..block_end+mid_block_needed).
    const char_type * const text2 =
      mid_block_reader->m_data - block_end;
    while (block_suf_begin + new_lcp < tail_begin &&
        new_lcp < pattern_length &&
        text2[block_suf_begin + new_lcp] == pat[new_lcp])
      ++new_lcp;

    // Update ret_lcp.
    *ret_lcp = new_lcp;

    // If the comparison ended with the
    // mismatch, we know the final answer
    // (either -1 or +1).
    if (block_suf_begin + new_lcp < tail_begin &&
        new_lcp < pattern_length) {

      // Invariant: pat[new_lcp] !=
      // text2[block_suf_begin + new_lcp].
      if (pat[new_lcp] < text2[block_suf_begin + new_lcp])
        return -1;
      else return +1;
    }
  }

  // No more information can be gained
  // using symbol comparisons.
  if (block_suf_begin + new_lcp >= tail_begin) {

    // Use gt to resolve comparison.
    const std::uint64_t gt_pos =
      pattern_begin + (tail_begin - block_suf_begin);
    if (gt_reader->access(text_length - gt_pos))
      return +1;
    else return -1;

  } else {

    // The only remaining option is that
    // new_lcp == pattern_length. In this case
    // we can resolve the comparison only due
    // to boundary case, i.e., if we have
    // pattern_begin + pattern_length == text_length.
    if (pattern_begin + pattern_length == text_length)
      return -1;
    else return 0;
  }
}

//=============================================================================
// Let 0 <= block_begin < block_end <= text_length and define block_size
// = block_end - block_begin. Further, let block_psa[0..block_size) be
// an array that contains a permutation of the set {0, .., block_size - 1}
// that describes the ascending lexicographical order of suffixes of text
// starting inside block (and extending beyond the end of the block).
//
// Assume that block_end <= tail_begin <= pattern_begin < text_length and
// that 0 <= range_lcp <= pattern_length <= text_length - pattern_begin.
// Further, assume that integers 0 <= left <= right <= block_size satisfy
// the following conditions:
// * any suffix of text starting at position block_begin + block_psa[i],
//   for 0 <= i < left is smaller than text[pattern_begin..text_length).
// * any suffix of text starting at position block_begin + block_psa[i],
//   for right <= i < block_size is larger than text[pattern_begin..
//   text_length).
// * any suffix of text starting at position block_begin + block_psa[i]
//   for left <= i < right has text[pattern_begin..pattern_begin +
//   range_lcp) as a prefix.
// Note that we do not assume the range [left, right) to be the maximal
// range satisfying the above conditions.
//
// The aim of this function is to compute integers newleft, newright
// that satisfy left <= newleft <= newright <= right and that satisfy
// all three conditions listed above (with left and right replaced with,
// respectively, newleft and newright), except with range_lcp in the
// third condition replaced by pattern_length. Note again that we do
// not require the range [newleft, newright) to the be the maximal range
// satisfying the specified conditions, hence the result of this function
// is not uniquely defined.
//
// Note that the function is not provided with the full text, nor the
// filename with the text. It is given the following information about
// the text and the ordering of its suffixes:
// * text[block_begin..block_end) and text[pattern_begin+range_lcp..
//   pattern_begin+pattern_length) are the accessible substrings of
//   text. They are available, respectively, via block[0..block_size)
//   and pat[range_lcp..pattern_length) variables. Any accesses to
//   these variables outside the specified ranges are undefined.
// * block_psa[0..block_size) as define above.
// * gt bitvector telling for every i in range (tail_begin..text_length]
//   whether the suffix text[i..text_length) is greater (then bit is 1)
//   than text[tail_begin..text_length).
// * mid_block_reader allows accessing (after first calling the wait
//   method to check that enough text has been read into RAM) the
//   substring text[block_end..tail_begin).
// The above information is sufficient to compute the answer.
//
// Note that gt bitvector is reversed (for technical reasons) and shifted
// (to provide convenient indexing), i.e., the bit of gt corresponding
// position i in (block_end..tail_end] (see definition above) is accessed
// as gt[text_length - i].
//=============================================================================
template<
  typename char_type,
  typename block_offset_type>
void refine_range(
    const std::uint64_t block_begin,
    const std::uint64_t block_end,
    const std::uint64_t pattern_begin,
    const std::uint64_t pattern_length,
    const std::uint64_t tail_begin,
    const std::uint64_t text_length,
    const std::uint64_t left,
    const std::uint64_t right,
    const std::uint64_t range_lcp,
    const char_type * const block,
    const char_type * const pat,
    const block_offset_type * const block_psa,
    background_block_reader * const mid_block_reader,
    multifile_bit_stream_reader * const gt_reader,
    std::uint64_t * const newleft,
    std::uint64_t * const newright) {

  // Initialize the pointer to text. Due to the
  // construction of this function we will never
  // access symbols in text which are not in RAM.
  const char_type * const text = block - block_begin;

  // Compute the tuning parameters
  // for string binary search.
#ifdef EM_STARTING_POS_MODULE_DEBUG_MODE
  std::uint64_t min_discrepancy = utils::random_int64(0, 10);
  std::uint64_t balancing_factor = utils::random_int64(1, 10);
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
      const std::uint64_t d = rlcp - llcp;
      const std::uint64_t logd = utils::log2ceil(d);
      mid = low + 1 +
        ((high - low - 1) * balancing_factor * logd) /
        (d + balancing_factor * logd);
    } else if (rlcp + min_discrepancy < llcp) {
      const std::uint64_t d = llcp - rlcp;
      const std::uint64_t logd = utils::log2ceil(d);
      mid = high - 1 -
        ((high - low - 1) * balancing_factor * logd) /
        (d + balancing_factor * logd);
    } else mid = (low + high) / 2;

    // Compare the pivot to the pattern.
    std::uint64_t new_lcp = 0;
    const std::uint64_t init_lcp = std::min(llcp, rlcp);
    const std::uint64_t block_suf_begin =
      block_begin + (std::uint64_t)block_psa[mid];
    if (lcp_compare(text_length, block_end, block_suf_begin,
          pattern_begin, pattern_length, tail_begin, init_lcp,
          text, pat, mid_block_reader, gt_reader, &new_lcp) <= 0) {
      high = mid;
      rlcp = new_lcp;
    } else {
      low = mid;
      llcp = new_lcp;
    }
  }

  // Store the lower bound
  // of the new range.
  *newleft = high;

  // rlcp < pattern_length means no suffix in block_psa
  // has pat[0..pattern_length) as a prefix or the ones
  // that do were already classified as larger. We
  // only need to do second binary search if rlcp
  // >= pattern_length.
  if (rlcp >= pattern_length) {
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
        const std::uint64_t d = rlcp - llcp;
        const std::uint64_t logd = utils::log2ceil(d);
        mid = low + 1 +
          ((high - low - 1) * balancing_factor * logd) /
          (d + balancing_factor * logd);
      } else if (rlcp + min_discrepancy < llcp) {
        const std::uint64_t d = llcp - rlcp;
        const std::uint64_t logd = utils::log2ceil(d);
        mid = high - 1 -
          ((high - low - 1) * balancing_factor * logd) /
          (d + balancing_factor * logd);
      } else mid = (low + high) / 2;

      // Compare the pivot to the pattern.
      std::uint64_t new_lcp = 0;
      const std::uint64_t init_lcp = std::min(llcp, rlcp);
      const std::uint64_t block_suf_begin =
        block_begin + (std::uint64_t)block_psa[mid];
      if (lcp_compare(text_length, block_end, block_suf_begin,
            pattern_begin, pattern_length, tail_begin, init_lcp,
            text, pat, mid_block_reader, gt_reader, &new_lcp) < 0) {
        high = mid;
        rlcp = new_lcp;
      } else {
        low = mid;
        llcp = new_lcp;
      }
    }
  }

  // Store the upper bound
  // of the new range.
  *newright = high;
}

//=============================================================================
// Let 0 <= block_begin < block_end <= tail_begin <= pattern_begin <
// text_length. The aim of this function is to compute the rank (i.e.,
// the number of smaller suffixes) of the suffix text[pattern_begin..
// text_length) among the suffixes of text starting inside text[
// block_begin..block_end) (and extending beyong the end of the block).
//
// The function is given the name of the file containing the text.
// The text itself is in principle sufficient to compute the answer.
// However, to perform the computation efficiently, the function receives
// the following information about the text and ordering of its suffixes.
// Let block_size = block_end - block_begin. The additional info is:
// * block[0..block_size) = text[block_begin..block_end).
// * block_psa[0..block_size) contains permutation of integer set
//   {0, .., block_size-1} that describes the ascending lexicographical
//   order of suffixes of text starting inside block (and extending
//   beyond the end of the block).
// * gt bitvector telling for every i in range (tail_begin..text_length]
//   whether the suffix text[i..text_length) is greater (then bit is 1)
//   than text[tail_begin..text_length).
// * mid_block_reader allows accessing (after first calling the wait
//   method to check that enough text has been read into RAM) the
//   substring text[block_end..tail_begin).
//
// Note that gt bitvector is reversed (for technical reasons) and shifted
// (to provide convenient indexing), i.e., the bit of gt corresponding
// position i in (tail_begin..text_length] (see definition above) is
// accessed as gt[text_length - i].
//=============================================================================
template<
  typename char_type,
  typename block_offset_type>
std::uint64_t compute_rank(
    const std::uint64_t block_begin,
    const std::uint64_t block_end,
    const std::uint64_t pattern_begin,
    const std::uint64_t text_length,
    const std::uint64_t tail_begin,
    const char_type * const block,
    const block_offset_type * const block_psa,
    const multifile * const tail_gt_begin_reversed,
    background_block_reader * const mid_block_reader,
    const std::string text_filename) {

  // Handle special case.
  if (pattern_begin == text_length)
    return 0;

  // Compute basic parameters.
  const std::uint64_t block_size = block_end - block_begin;
  const std::uint64_t max_pat_symbols_needed = std::min(
      text_length - pattern_begin, tail_begin - block_begin);

  // Create the reader of the gt bitvector.
  multifile_bit_stream_reader * const gt_reader
    = new multifile_bit_stream_reader(tail_gt_begin_reversed);

  // The computation is done with the binary search over the block_psa
  // array. However, a comparison of the pattern with a single suffix
  // from block_psa could require accessing nearly whole pattern. We
  // cannot afford to store it, so in the worst case we would have to
  // stream the whole pattern many times. To prevent this, we perform
  // the computation in chunks. Each time we read a small chunk of the
  // pattern into RAM, we refine the range of suffixes in the suffix
  // array using symbols in the current chunk. This way we use little
  // extra RAM and are we guaranteed to only read the pattern once.
  // Further, in most cases, refining a range based on some initial
  // chunk already reduced the range to a single element and the
  // computation can be finished.
#ifdef EM_STARTING_POS_MODULE_DEBUG_MODE
  std::uint64_t chunk_length =
    utils::random_int64(1, 10);
  background_chunk_reader *chunk_reader =
    new background_chunk_reader(text_filename, pattern_begin,
        pattern_begin + max_pat_symbols_needed, chunk_length);
#else
  background_chunk_reader *chunk_reader =
    new background_chunk_reader(text_filename, pattern_begin,
        pattern_begin + max_pat_symbols_needed);
#endif

  // Find the answer using binary search. Each step
  // of the loop below processes a single chunk and
  // refines the range until left == right.
  std::uint64_t range_lcp = 0;
  std::uint64_t left = 0;
  std::uint64_t right = block_size;
  while (left != right) {

    // Read the next chunk of the pattern.
    const std::uint64_t this_chunk_length = std::min(
        max_pat_symbols_needed - range_lcp,
        chunk_reader->get_chunk_size());
    const std::uint64_t pattern_length = range_lcp + this_chunk_length;
    chunk_reader->wait(pattern_begin + pattern_length);

    // Compute new range boundaries. Invariant: the answer
    // (i.e., final rank) is in in the range [left, right].
    // Invariant: all suffixes in psa[left..right) share a
    // common prefix of length range_lcp with the pattern.
    // Invariant:
    // chunk_reader->m_chunk[0..this_chunk_length) == text[
    // pattern_begin+range_lcp..pattern_begin+pattern_length).
    std::uint64_t newleft = 0;
    std::uint64_t newright = 0;
    refine_range(block_begin, block_end, pattern_begin,
        pattern_length, tail_begin, text_length, left, right,
        range_lcp, block, chunk_reader->m_chunk - range_lcp,
        block_psa, mid_block_reader, gt_reader, &newleft, &newright);

    // Update range boundaries.
    left = newleft;
    right = newright;
    range_lcp = pattern_length;
  }

  // Clean up.
  delete chunk_reader;
  delete gt_reader;

  // Return the result.
  return left;
}

//=============================================================================
// The aim of this function is to compute the rank (i.e., the number
// of smaller suffixes) among the suffixes of text starting inside
// text[block_begin..block_end) (0 <= block_begin < block_end <= text_length)
// for each of the strings in the sparse set of suffixes of text called
// patterns. The set of patterns is defined as the set of all suffixes
// of text that start at positions i such that tail_begin <= i and
// (i - tail_begin) is a multiple of patterns_dist > 0. It is guaranteed
// that block_end <= tail_begin.
//
// The function is given the name of the file containing the text.
// The text itself is in principle sufficient to compute all the ranks.
// However, to perform the computation efficiently, the function receives
// the following information about the text and ordering of its suffixes.
// Let block_size = block_end - block_begin. The additional info is:
// * block[0..block_size) = text[block_begin..block_end).
// * block_psa[0..block_size) contains permutation of integer set
//   {0, .., block_size-1} that describes the ascending lexicographical
//   order of suffixes of text starting inside block (and extending
//   beyond the end of the block).
// * gt bitvector telling for every i in range (tail_begin..text_length]
//   whether the suffix text[i..text_length) is greater (then bit is 1)
//   than text[tail_begin..text_length).
//
// Note that gt bitvector is reversed (for technical reasons) and shifted
// (to provide convenient indexing), i.e., the bit of gt corresponding
// position i in (tail_begin..text_length] (see definition above) is
// accessed as gt[text_length - i].
//=============================================================================
template<
  typename char_type,
  typename block_offset_type>
void compute_ranks(
    const std::uint64_t block_begin,
    const std::uint64_t block_end,
    const std::uint64_t text_length,
    const std::uint64_t tail_begin,
    const std::uint64_t patterns_dist,
    const char_type * const block,
    const block_offset_type * const block_psa,
    const multifile * const tail_gt_begin_reversed,
    const std::string text_filename,
    std::vector<std::uint64_t> &result) {

  // Compute the number of patterns.
  const std::uint64_t tail_length = text_length - tail_begin;
  const std::uint64_t n_patterns = (tail_length +
      patterns_dist - 1) / patterns_dist;

  // During the computation we might require access to the
  // substring of text: text[block_end..tail_begin). It is
  // allocated right away and then, since it is rarely used,
  // we start reading its symbols in the background while the
  // computation continues. Each thread that wants to access
  // those symbols first checks whether sufficient prefix has
  // been already read and waits if necessary.
  background_block_reader * const mid_block_reader =
    new background_block_reader(text_filename,
        block_end, tail_begin - block_end);

  // Compute the ranks.
  std::vector<std::uint64_t> res(n_patterns);
  {

#ifdef _OPENMP

    // Parallel version.
    #pragma omp parallel num_threads(n_patterns)
    {

      // Compute pattern id and starting position.
      const std::uint64_t pattern_id = omp_get_thread_num();
      const std::uint64_t pattern_begin = tail_begin +
        pattern_id * patterns_dist;

      // Compute the rank.
      res[pattern_id] = compute_rank(
          block_begin, block_end, pattern_begin, text_length,
          tail_begin, block, block_psa, tail_gt_begin_reversed,
          mid_block_reader, text_filename);
    }
#else

    // Sequential version.
    for (std::uint64_t pattern_id = 0;
        pattern_id < n_patterns; ++pattern_id) {

      // Compute pattern starting position.
      const std::uint64_t pattern_begin = tail_begin +
        pattern_id * patterns_dist;

      // Compute the rank.
      res[pattern_id] = compute_rank(
          block_begin, block_end, pattern_begin, text_length,
          tail_begin, block, block_psa, tail_gt_begin_reversed,
          mid_block_reader, text_filename);
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

#endif  // __SRC_PSASCAN_SRC_COMPUTE_RANKS_HPP_INCLUDED
