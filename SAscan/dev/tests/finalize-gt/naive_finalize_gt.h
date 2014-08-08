#ifndef __NAIVE_FINALIZE_GT_H
#define __NAIVE_FINALIZE_GT_H


#include "bitvector.h"


void naive_finalize_gt(unsigned char *text, long text_length, long left_block_beg,
    long left_block_size, bitvector *gt_out) {
  long left_block_end = left_block_beg + left_block_size;

  for (long i = left_block_beg; i < left_block_end; ++i) {
    // compute lcp(text[left_block_beg..), text[i..)).
    long lcp = 0L;
    while (i + lcp < text_length && text[left_block_beg + lcp] == text[i + lcp]) ++lcp;
    if (i + lcp != text_length && text[i + lcp] > text[left_block_beg + lcp])
      gt_out->set(i - left_block_beg);
  }
}

#endif  // __NAIVE_FINALIZE_GT_H
