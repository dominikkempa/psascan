#ifndef __EM_COMPUTE_INITIAL_RANKS_INCLUDED
#define __EM_COMPUTE_INITIAL_RANKS_INCLUDED

#include <algorithm>

#include "background_chunk_reader.h"
#include "multifile_bit_stream_reader.h"
#include "utils.h"


#define EM_STARTING_POS_MODULE_DEBUG_MODE

inline int lcp_compare_special(
    const unsigned char *text,  // only text[block_suf_beg..block_end) will be accessed
    long text_length,
    long block_end,             // wrt to text beg
    long block_suf_beg,         // wrt to text beg
    const unsigned char *pat,   // only pat[lcp..pat_length) will be accessed
    long pat_beg,               // wrt to text beg
    long pat_length,
    multifile_bit_stream_reader &gt_reader,
    long &lcp) {
  while (block_suf_beg + lcp < block_end && lcp < pat_length &&
      text[block_suf_beg + lcp] == pat[lcp]) ++lcp;
  if (block_suf_beg + lcp >= block_end) {
    if (gt_reader.access(text_length - (pat_beg + (block_end - block_suf_beg)))) return 1;
    else return -1;
  } else if (lcp == pat_length) return 0;
  else {
    if (pat[lcp] > text[block_suf_beg + lcp]) return 1;
    else return -1;
  } 
}

template<typename saidx_t>
void refine_range(
    const unsigned char *block,
    const saidx_t *block_psa,
    long block_beg,  // wrt to text beg
    long block_end,  // same here
    long pat_beg,    // same here
    long text_length,
    long left,
    long right,
    long old_lcp,
    long new_lcp,
    const unsigned char *pat,  // only pat[old_lcp..new_lcp) can and will be accessed
    multifile_bit_stream_reader &gt_reader,
    long &newleft,
    long &newright) {
  long low = left - 1;
  long high = right;
  long llcp = old_lcp;
  long rlcp = old_lcp;

#ifdef EM_STARTING_POS_MODULE_DEBUG_MODE
  long min_discrepancy = utils::random_long(0L, 10L);
  long balancing_factor = utils::random_long(1L, 10L);
#else
  static const long min_discrepancy = (1L << 16);
  static const long balancing_factor = 64L;
#endif

  const unsigned char *text = block - block_beg;
  while (low + 1 != high) {
    // Invariant: newleft is in the range (low, high].
    long lcp = std::min(llcp, rlcp);
    long mid = 0L;
    if (llcp + min_discrepancy < rlcp) {
      long d = rlcp - llcp;
      long logd = utils::log2ceil(d);
      mid = low + 1 + ((high - low - 1) * balancing_factor * logd) / (d + balancing_factor * logd);
    } else if (rlcp + min_discrepancy < llcp) {
      long d = llcp - rlcp;
      long logd = utils::log2ceil(d);
      mid = high - 1 - ((high - low - 1) * balancing_factor * logd) / (d + balancing_factor * logd);
    } else mid = (low + high) / 2;

    if (lcp_compare_special(text, text_length, block_end, block_beg + (long)block_psa[mid],
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
      long lcp = std::min(llcp, rlcp);
      long mid = 0L;
      if (llcp + min_discrepancy < rlcp) {
        long d = rlcp - llcp;
        long logd = utils::log2ceil(d);
        mid = low + 1 + ((high - low - 1) * balancing_factor * logd) / (d + balancing_factor * logd);
      } else if (rlcp + min_discrepancy < llcp) {
        long d = llcp - rlcp;
        long logd = utils::log2ceil(d);
        mid = high - 1 - ((high - low - 1) * balancing_factor * logd) / (d + balancing_factor * logd);
      } else mid = (low + high) / 2;

      if (lcp_compare_special(text, text_length, block_end, block_beg + (long)block_psa[mid],
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

template<typename saidx_t>
void compute_single_starting_position(
    const unsigned char *block,
    const saidx_t *block_psa,
    long block_beg,  // wrt to text beg
    long block_end,  // same here
    long pat_beg,    // same here
    long text_length,
    long max_lcp,
    std::string text_filename,
    const multifile *tail_gt_begin_reversed,
    std::pair<long, long> &result) {
  long block_size = block_end - block_beg;
  long pat_end = std::min(text_length, pat_beg + max_lcp);

  multifile_bit_stream_reader gt_reader(tail_gt_begin_reversed);

  // Reads text[pat_beg..pat_end) in chunks.
#ifdef EM_STARTING_POS_MODULE_DEBUG_MODE
  long chunk_length = utils::random_long(1L, 10L); 
  background_chunk_reader *chunk_reader =
    new background_chunk_reader(text_filename, pat_beg, pat_end, chunk_length);
#else
  background_chunk_reader *chunk_reader =
    new background_chunk_reader(text_filename, pat_beg, pat_end);
#endif

  // The current range is [left, right).
  long left = 0;
  long right = block_size;
  long lcp = 0;

  while (left != right && lcp < max_lcp) {
    long this_chunk_length = std::min(max_lcp - lcp, chunk_reader->get_chunk_size());
    long new_lcp = lcp + this_chunk_length;
    chunk_reader->wait(pat_beg + new_lcp);

    // Invariant:
    //   reader->chunk[0..chunk_length) = pattern[lcp..new_lcp).
    long newleft = 0;
    long newright = 0;
    refine_range(block, block_psa, block_beg, block_end, pat_beg, text_length, left,
        right, lcp, new_lcp, chunk_reader->m_chunk - lcp, gt_reader, newleft, newright);
    left = newleft;
    right = newright;
    lcp = new_lcp;
  }

  delete chunk_reader;

  result = std::make_pair(left, right);
}

#endif  // __EM_COMPUTE_INITIAL_RANKS_INCLUDED

