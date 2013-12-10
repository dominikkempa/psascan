#ifndef __SUFFIX_RANKING_H
#define __SUFFIX_RANKING_H

#include <algorithm>
#include <cstring>

// Update ms-decomposition of T[0..n) from T[0..n-1).
void next(unsigned char *T, int n, int &s, int &p) {
  if (n == 1) { s = 0; p = 1; return; }
  int i = n - 1, r = (i - s) % p;
  while (i < n) {
    unsigned char a = T[s + r], b = T[i];
    if (a > b) { p = i - s + 1; r = 0; }
    else if (a < b) { i -= r; s = i; p = 1; r = 0; }
    else { ++r; if (r == p) r = 0; } ++i;
  }
}

void compute_bitmap(unsigned char *A, int A_length,
                    unsigned char *B, int B_length,
                    unsigned char *gt_A, unsigned char *gt_eof) {
  std::fill(gt_eof, gt_eof + B_length, 0);
  int i = 0, el = 0, s = 0, p = 0, i_max = 0, el_max = 0, s_max = 0, p_max = 0;
  while (i < B_length) {
    while (i + el < B_length && el < A_length && B[i + el] == A[el]) next(A, ++el, s, p);
    if (el == A_length || (i + el == B_length && !gt_A[el]) ||
        (i + el < B_length && B[i + el] > A[el])) gt_eof[i] = 1;
    int j = i_max;
    if (el > el_max) {
      std::swap(el, el_max);
      std::swap(s, s_max);
      std::swap(p, p_max);
      i_max = i;
    }
    if (p && 3 * p <= el && !memcmp(A, A + p, s)) {
      std::copy(gt_eof + j + 1, gt_eof + j + p, gt_eof + i + 1);
      i += p; el -= p;
    } else {
      int h = (el / 3) + 1;
      std::copy(gt_eof + j + 1, gt_eof + j + h, gt_eof + i + 1);
      i += h; el = 0; s = 0; p = 0;
    }
  }
}

#endif // __SUFFIX_RANKING_H
