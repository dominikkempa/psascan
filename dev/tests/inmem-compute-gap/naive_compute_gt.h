#ifndef __NAIVE_COMPUTE_GT_H
#define __NAIVE_COMPUTE_GT_H

#include "bitvector.h"

//==============================================================================
// Compute gt bitvector of length gt_length such that
// gt[i] == 1 iff text[suf_start + i .. text_length)
// is greater than the suffix text[suf_start .. text_lengtH).
//==============================================================================
void naive_compute_gt(unsigned char *text, long text_length,
    long suf_start, long gt_length, bitvector* &gt) {
  gt = new bitvector(gt_length);
  for (long i = suf_start; i < suf_start + gt_length; ++i) {
    long lcp = 0;

    // Compute lcp(suf_start, i).
    while (i + lcp < text_length &&
        text[suf_start + lcp] == text[i + lcp]) ++lcp;

    // Compute gt[i].
    if (i + lcp != text_length &&
        text[i + lcp] > text[suf_start + lcp]) {
      gt->set(i - suf_start);
    }
  }
}

#endif  // __NAIVE_COMPUTE_GT_H
