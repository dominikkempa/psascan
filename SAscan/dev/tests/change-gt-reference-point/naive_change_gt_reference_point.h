#ifndef __NAIVE_CHANGE_GT_REFERENCE_POINT_H
#define __NAIVE_CHANGE_GT_REFERENCE_POINT_H

#include "bitvector.h"


void naive_change_gt_reference_point(unsigned char *text, long text_length,
    long block_beg, long block_end, bitvector* &gt_out) {
  gt_out = new bitvector(block_end - block_beg + 1);
  for (long i = block_beg; i <= block_end; ++i) {
    long lcp = 0L;
    while (i + lcp < text_length && text[i + lcp] == text[block_beg + lcp]) ++lcp;
    if (i + lcp < text_length && text[i + lcp] > text[block_beg + lcp])
      gt_out->set(i - block_beg);
  }
}


#endif  // __NAIVE_CHANGE_GT_REFERENCE_POINT_H
