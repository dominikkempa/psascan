#ifndef __INMEM_SASCAN_INMEM_SMALLER_SUFFIXES_H_INCLUDED
#define __INMEM_SASCAN_INMEM_SMALLER_SUFFIXES_H_INCLUDED

#include "../../bitvector.h"
#include "pagearray.h"
#include "../../multifile.h"
#include "../../multifile_bit_stream_reader.h"


namespace inmem_sascan_private {


//==============================================================================
// Return true iff text[i..length) < text[j..length). To speed up the
// computation we are given a lower bound for the lcp values, that is
// lcp <= lcp(i, j).
//==============================================================================
// We assume that j < i.
//==============================================================================
int lcp_compare(unsigned char *text, long i, long j, long maxlcp, long &lcp) {
  while (lcp < maxlcp && text[i + lcp] == text[j + lcp])
    ++lcp;

  if (lcp == maxlcp) return 0;
  else if (text[i + lcp] < text[j + lcp]) return -1;
  else return 1;
}

template<typename pagearray_type>
void compute_other_starting_position(unsigned char *text, long block_beg, long block_end,
    long suf_start, long maxlcp, const pagearray_type &partial_bwtsa, std::pair<long, long> &ret) {

  long block_size = block_end - block_beg;

  // 1
  //
  // Find the smallest j (at the and will be stored in `right'),
  // such that
  // - text[suf_start..text_length) < text[psa[j]..text_length) or,
  // - lcp(suf_start, psa[j]) = maxlcp, i.e., the we cannot determine the order.
  // If there is no such j (i.e., all suffixes were smaller), the answer will
  // be equal to block_size.

  static const long min_discrepancy = (1L << 16);
  static const long balancing_constant = 64L;

  long left = -1L, right = block_size;
  long llcp = 0L, rlcp = 0L;
  while (left + 1 != right) {
    // Invariant: llcp = lcp(left, suf_start), rlcp = lcp(right, suf_start).
    // Invariant: the index we are searching for is in the range (beg..end].

    // Compute mid.
    // Valid values for mid are: left + 1, .., right - 1.
    long mid = 0L;
    if (llcp + min_discrepancy < rlcp) {
      // Choose the pivot that split the range into two parts of sizes
      // with ratio equal to logd / d.
      long d = rlcp - llcp;
      long logd = utils::log2ceil(d);
      mid = left + 1 + ((right - left - 1) * balancing_constant * logd) / (d + balancing_constant * logd);
    } else if (rlcp + min_discrepancy < llcp) {
      // Choose the pivot that split the range into two parts of sizes
      // with ratio equal to logd / d.
      long d = llcp - rlcp;
      long logd = utils::log2ceil(d);
      mid = right - 1 - ((right - left - 1) * balancing_constant * logd) / (d + balancing_constant * logd);
    } else {
      // Discrepancy between lcp values is small, use standard binary search.
      mid = (left + right) / 2;
    }

    long lcp = std::min(llcp, rlcp);

    if (lcp_compare(text, suf_start, block_beg + partial_bwtsa[mid].sa, maxlcp, lcp) <= 0) {
      right = mid;
      rlcp = lcp;
    } else {
      left = mid;
      llcp = lcp;
    }
  }

  long lower_bound = right - 1;
  long upper_bound = right;

  if (rlcp == maxlcp) {
    right = block_size;
    rlcp = 0L;

    while (left + 1 != right) {
      // Compute mid.
      // Valid values for mid are: left + 1, .., right - 1.
      long mid = 0L;
      if (llcp + min_discrepancy < rlcp) {
        // Choose the pivot that split the range into two parts of sizes
        // with ratio equal to logd / d.
        long d = rlcp - llcp;
        long logd = utils::log2ceil(d);
        mid = left + 1 + ((right - left - 1) * balancing_constant * logd) / (d + balancing_constant * logd);
      } else if (rlcp + min_discrepancy < llcp) {
        // Choose the pivot that split the range into two parts of sizes
        // with ratio equal to logd / d.
        long d = llcp - rlcp;
        long logd = utils::log2ceil(d);
        mid = right - 1 - ((right - left - 1) * balancing_constant * logd) / (d + balancing_constant * logd);
      } else {
        // Discrepancy between lcp values is small, use standard binary search.
        mid = (left + right) / 2;
      }

      long lcp = std::min(llcp, rlcp);

      if (lcp_compare(text, suf_start, block_beg + partial_bwtsa[mid].sa, maxlcp, lcp) == -1) {
        right = mid;
        rlcp = lcp;
      } else {
        left = mid;
        llcp = lcp;
      }
    }

    upper_bound = right;
  }

  // Invariant: the position we are looking for is in the range (lower_bound..upper_bound].
  // Moreover:  either the range is of size one, and then we just know the rank of the suffix,
  // or if it's of size > 1 and then the range (lower_bound..upper_bound) contains all
  // suffixes starting inside text[block_beg..block_end) that have text[suf_start..suf_start..maxlcp)
  // as a prefix.

  ret = std::make_pair(lower_bound, upper_bound);
}

}  // namespace inmem_sascan

#endif // __INMEM_SMALLER_SUFFIXES_H_INCLUDED
