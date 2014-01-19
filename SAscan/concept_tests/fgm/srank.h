#ifndef __SRANK_H
#define __SRANK_H

#include <algorithm>
#include <cstring>

#include "bitvector.h"
#include "stream.h"

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

void compute_gt_eof_bv(unsigned char *A, int A_length,
                    unsigned char *B, int B_length,
                    bitvector *gt_head_bv, bitvector *gt_eof_bv) {
  int i = 0, el = 0, s = 0, p = 0, i_max = 0, el_max = 0, s_max = 0, p_max = 0;
  while (i < B_length) {
    while (i + el < B_length && el < A_length && B[i + el] == A[el]) next(A, ++el, s, p);
    if (el == A_length || (i + el == B_length && gt_head_bv->get(A_length - 1 - el) == false) ||
        (i + el < B_length && B[i + el] > A[el])) gt_eof_bv->set(i);
    int j = i_max;
    if (el > el_max) {
      std::swap(el, el_max);
      std::swap(s, s_max);
      std::swap(p, p_max);
      i_max = i;
    }
    if (p && 3 * p <= el && !memcmp(A, A + p, s)) {
      for (int k = 1; k < p; ++k)
        if (gt_eof_bv->get(j + k)) gt_eof_bv->set(i + k);
      i += p; el -= p;
    } else {
      int h = (el / 3) + 1;
      for (int k = 1; k < h; ++k)
        if (gt_eof_bv->get(j + k)) gt_eof_bv->set(i + k);
      i += h; el = 0; s = 0; p = 0;
    }
  }
}

// Given T[0..n), compute bitmap of length n such that new_gt_head_bv[i] == 1
// iff T[i..n) > T. Returns the number of suffixes T[i..n) < T.
int compute_new_gt_head_bv(unsigned char *T, int n, bitvector *new_gt_head_bv) {
  int whole_suffix_rank = n - 1;
  int i = 0, el = 0, s = 0, p = 0;
  int i_max = 0, el_max = 0, s_max = 0, p_max = 0;
  while (i < n) {
    while (i + el < n && el < n && T[i + el] == T[el]) next(T, ++el, s, p);
    if (i + el < n && (el == n || T[i + el] > T[el])) {
      new_gt_head_bv->set(i);
      --whole_suffix_rank;
    }
    int j = i_max;
    if (el > el_max) {
      std::swap(el, el_max);
      std::swap(s, s_max);
      std::swap(p, p_max);
      i_max = i;
    }
    if (p && 3 * p <= el && !memcmp(T, T + p, s)) {
      for (int k = 1; k < p; ++k)
        if (new_gt_head_bv->get(j + k)) {
          new_gt_head_bv->set(i + k);
          --whole_suffix_rank;
        }
      i += p; el -= p;
    } else {
      int h = (el / 3) + 1;
      for (int k = 1; k < h; ++k)
        if (new_gt_head_bv->get(j + k)) {
          new_gt_head_bv->set(i + k);
          --whole_suffix_rank;
        }
      i += h; el = 0; s = 0; p = 0;
    }
  }

  return whole_suffix_rank;
}

#endif // __SRANK_H
