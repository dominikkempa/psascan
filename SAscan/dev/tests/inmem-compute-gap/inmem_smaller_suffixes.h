#ifndef __INMEM_SMALLER_SUFFIXES_H_INCLUDED
#define __INMEM_SMALLER_SUFFIXES_H_INCLUDED

#include "bitvector.h"

//==============================================================================
// Return true iff text[i..length) < text[j..length). To speed up the
// computation we are given a lower bound for the lcp values, that is
// lcp <= lcp(i, j).
//==============================================================================
bool lcp_compare(unsigned char *text, long length, long i, long j, long &lcp) {
  while (i + lcp < length && j + lcp < length && text[i + lcp] == text[j + lcp])
    ++lcp;

  return (j + lcp != length &&
      (i + lcp == length || text[i + lcp] < text[j + lcp]));
}

template<typename T>
void inmem_smaller_suffixes(unsigned char *text, long text_length,
    long block_beg, long block_end, long suf_start,  T *partial_sa, long &ret) {

  //----------------------------------------------------------------------------
  // Find the smallest j such that text[partial_sa[j]..text_length) >
  // text[suf_start..text_length) using binary search or block_size if there
  // is not such j.
  //----------------------------------------------------------------------------
  long left = -1L, right = block_end - block_beg;
  long llcp = 0L, rlcp = 0L;
  while (left + 1 != right) {
    // Invariant: llcp = lcp(left, suf_start), rlcp = lcp(right, suf_start).
    long mid = (left + right);
    long lcp = std::min(llcp, rlcp);

    if (lcp_compare(text, text_length, suf_start, block_beg + partial_sa[mid], lcp)) {
      right = mid;
      rlcp = lcp;
    } else {
      left = mid;
      llcp = lcp;
    }
  }

  ret = right;
}

#endif // __INMEM_SMALLER_SUFFIXES_H_INCLUDED
