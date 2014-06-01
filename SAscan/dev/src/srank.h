//==============================================================================
// String range matching algorithms described in
//
// Juha Kärkkäinen, Dominik Kempa, Simon J. Puglisi:
// String Range Matching. In Proc. CPM 2014.
//==============================================================================

#ifndef __SRANK_H_INCLUDED
#define __SRANK_H_INCLUDED

#include <cstring>
#include <algorithm>

#include "bitvector.h"
#include "smaller_suffixes.h"

// Update ms-decomposition of T[0..n) from T[0..n-1).
void next(unsigned char *T, long n, long &s, long &p, long &r) {
  if (n == 1) { s = 0; p = 1; r = 0; return; }
  long i = n - 1;
  while (i < n) {
    unsigned char a = T[s + r], b = T[i];
    if (a > b) { p = i - s + 1; r = 0; }
    else if (a < b) { i -= r; s = i; p = 1; r = 0; }
    else { ++r; if (r == p) r = 0; } ++i;
  }
}


void compute_gt_eof_bv(unsigned char *A, long A_length,
                       unsigned char *B, long B_length,
                       gt_accessor &gt, long gt_length, bitvector *gt_eof_bv) {
  long i = 0, el = 0, s = 0, p = 0, r = 0;
  long i_max = 0, el_max = 0, s_max = 0, p_max = 0, r_max = 0;
  while (i < B_length) {
    while (i + el < B_length && el < A_length && B[i + el] == A[el]) next(A, ++el, s, p, r);
    if (el == A_length || (i + el == B_length && !gt[gt_length - el]) ||
        (i + el < B_length && B[i + el] > A[el])) gt_eof_bv->set(i);
    long j = i_max;
    if (el > el_max) {
      std::swap(el, el_max);
      std::swap(s, s_max);
      std::swap(p, p_max);
      std::swap(r, r_max);
      i_max = i;
    }
    if (p && 3 * p <= el && !memcmp(A, A + p, s)) {
      for (long k = 1; k < p; ++k)
        if (gt_eof_bv->get(j + k)) gt_eof_bv->set(i + k);
      i += p; el -= p;
    } else {
      long h = (el / 3) + 1;
      for (long k = 1; k < h; ++k)
        if (gt_eof_bv->get(j + k)) gt_eof_bv->set(i + k);
      i += h; el = 0; s = 0; p = 0;
    }
  }
}

// Given T[0..n), compute bitmap of length n such that new_gt_head_bv[i] == 1
// iff T[i..n) > T. Returns the number of suffixes T[i..n) < T.
long compute_new_gt_head_bv(unsigned char *T, long n, bitvector *new_gt_head_bv) {
  long whole_suffix_rank = n - 1;
  long i = 1, el = 0, s = 0, p = 0, r = 0;
  long i_max = 0, el_max = 0, s_max = 0, p_max = 0, r_max = 0;
  while (i < n) {
    while (i + el < n && el < n && T[i + el] == T[el]) next(T, ++el, s, p, r);
    if (i + el < n && (el == n || T[i + el] > T[el])) {
      // To avoid reversing, fill new_gt_head_bv backwards.
      // new_gt_head_bv->set(i);
      new_gt_head_bv->set(n - 1 - i);
      --whole_suffix_rank;
    }
    long j = i_max;
    if (el > el_max) {
      std::swap(el, el_max);
      std::swap(s, s_max);
      std::swap(p, p_max);
      std::swap(r, r_max);
      i_max = i;
    }
    if (p && 3 * p <= el && !memcmp(T, T + p, s)) {
      for (long k = 1; k < p; ++k) // Optimized as above.
        if (new_gt_head_bv->get(n - 1 - (j + k))) {
          new_gt_head_bv->set(n - 1 - (i + k));
          --whole_suffix_rank;
        }
        //if (new_gt_head_bv->get(j + k)) {
        //  new_gt_head_bv->set(i + k);
        //  --whole_suffix_rank;
        //}
      i += p; el -= p;
    } else {
      long h = (el / 3) + 1;
      for (long k = 1; k < h; ++k) // Optimized as above.
        if (new_gt_head_bv->get(n - 1 - (j + k))) {
          new_gt_head_bv->set(n - 1 - (i + k));
          --whole_suffix_rank;
        }
        //if (new_gt_head_bv->get(j + k)) {
        //  new_gt_head_bv->set(i + k);
        //  --whole_suffix_rank;
        //}
      i += h; el = 0; s = 0; p = 0;
    }
  }

  return whole_suffix_rank;
}

#endif // __SRANK_H_INCLUDED
