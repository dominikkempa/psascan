#ifndef __NAIVE_COMPUTE_PARTIAL_SA_H
#define __NAIVE_COMPUTE_PARTIAL_SA_H

#include "divsufsort.h"

void naive_compute_partial_sa(unsigned char *text, long length,
    long block_beg, long block_end, int* &partial_sa) {

  long block_size = block_end - block_beg;
  partial_sa = new int[block_size];

  int *sa = new int[length];
  divsufsort(text, sa, (int)length);

  for (long i = 0, ptr = 0; i < length; ++i)
    if (block_beg <= sa[i] && sa[i] < block_end)
      partial_sa[ptr++] = sa[i] - block_beg;

  delete[] sa;
}

#endif  // __NAIVE_COMPUTE_PARTIAL_SA_H
