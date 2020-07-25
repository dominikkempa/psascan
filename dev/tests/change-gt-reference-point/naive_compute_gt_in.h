#ifndef __NAIVE_COMPUTE_GT_IN_H
#define __NAIVE_COMPUTE_GT_IN_H

#include "bitvector.h"


void naive_compute_gt_in(unsigned char *text, long text_length,
    long block_beg, long block_end, bitvector* &gt_in) {
  long block_size = block_end - block_beg;
  gt_in = new bitvector(block_size);
  for (long i = block_beg; i < block_end; ++i) {
    long lcp = 0L;
    while (block_end + lcp < text_length && text[block_end + lcp] == text[i + lcp]) ++lcp;
    if (block_end + lcp == text_length ||
        (block_end + lcp != text_length && text[i + lcp] > text[block_end + lcp]))
      gt_in->set(i - block_beg);
  }
}


#endif  // __NAIVE_COMPUTE_GT_IN
