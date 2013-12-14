#ifndef __SRANK_H
#define __SRANK_H

#include <algorithm>
#include <cstring>

#include "bitvector.h"
#include "stream.h"

// Update ms-decomposition of T[0..n) from T[0..n-1).
void next(unsigned char *T, long n, long &s, long &p) {
  if (n == 1) { s = 0; p = 1; return; }
  long i = n - 1, r = (i - s) % p;
  while (i < n) {
    unsigned char a = T[s + r], b = T[i];
    if (a > b) { p = i - s + 1; r = 0; }
    else if (a < b) { i -= r; s = i; p = 1; r = 0; }
    else { ++r; if (r == p) r = 0; } ++i;
  }
}

void compute_gt_eof(unsigned char *A, long A_length,
                    unsigned char *B, long B_length,
                    bitvector *gt_head_bv, bitvector *gt_eof) {
  long i = 0, el = 0, s = 0, p = 0, i_max = 0, el_max = 0, s_max = 0, p_max = 0;
  while (i < B_length) {
    while (i + el < B_length && el < A_length && B[i + el] == A[el]) next(A, ++el, s, p);
    if (el == A_length || (i + el == B_length && gt_head_bv->get(A_length - 1 - el) == false) ||
        (i + el < B_length && B[i + el] > A[el])) gt_eof->set(i);
    long j = i_max;
    if (el > el_max) {
      std::swap(el, el_max);
      std::swap(s, s_max);
      std::swap(p, p_max);
      i_max = i;
    }
    if (p && 3 * p <= el && !memcmp(A, A + p, s)) {
      for (long k = 1; k < p; ++k)
        if (gt_eof->get(j + k)) gt_eof->set(i + k);
      i += p; el -= p;
    } else {
      long h = (el / 3) + 1;
      for (long k = 1; k < h; ++k)
        if (gt_eof->get(j + k)) gt_eof->set(i + k);
      i += h; el = 0; s = 0; p = 0;
    }
  }
}

#endif // __SRANK_H
