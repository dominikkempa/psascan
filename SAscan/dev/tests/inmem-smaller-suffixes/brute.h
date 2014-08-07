#ifndef __BRUTE_H
#define __BRUTE_H

#include "bitvector.h"

void inmem_smaller_suffixes_brute(unsigned char *text, long text_length,
    long block_beg, long block_end, long suf_start, long &res) {
  res = 0;
  for (long i = block_beg; i < block_end; ++i) {
    long lcp = 0;

    // Compute lcp(i, suf_start).
    while (i + lcp < text_length && suf_start + lcp < text_length &&
        text[i + lcp] == text[suf_start + lcp]) ++lcp;

    // Update the answer.
    if (suf_start + lcp != text_length && ((i + lcp == text_length) || 
        (i + lcp < text_length && text[i + lcp] < text[suf_start + lcp])))
      ++res;
  }
}

#endif  // __BRUTE_H
