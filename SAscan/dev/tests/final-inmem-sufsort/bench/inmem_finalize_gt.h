#ifndef __INMEM_FINALIZE_GT_H
#define __INMEM_FINALIZE_GT_H


#include <algorithm>
#include <cstring>
#include <thread>

#include "bitvector.h"


//==============================================================================
// Compute ms-decomposition of T[0..n) from ms-decomposition of T[0..n-1).
// The result is returned via updated values s, p, r.
//==============================================================================
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


//==============================================================================
// Compute gt_out[range_beg..range_end).
//
// gt_in is indexed from 0, with respect to the end of left block and is
// guaranteed to be long enough
//-----------------------------------------------------------------------------
// Rule of thumb in this kind of computation based on string range matching:
//
//   pat should be always a pointer to the beginning of
//   pattern, and should be 0-indexed.
//
//==============================================================================
void finalize_gt_aux(unsigned char *text, long text_length, long left_block_beg,
    long left_block_size, long range_beg, long range_end, bitvector *gt_in,
    bitvector *gt_out) {
  // We avoid the redundant long comparisong, gt_out[0] (correctly) stays 0.
  if (range_beg == left_block_beg) ++range_beg;

  // Run string range matching, i is relative to left_block_beg,
  // its valid values are inside range [range_beg - left_block_beg .. 
  // range_end - left_block_beg).
  unsigned char *pat = text + left_block_beg, *txt = pat;
  long i = range_beg - left_block_beg, el = 0L, s = 0L, p = 0L, r = 0L;
  long i_max = i, el_max = 0L, s_max = 0L, p_max = 0L, r_max = 0L;
  while (i < range_end - left_block_beg) {
    // Compute lcp(text[left_block_beg..), text[left_block_beg+i..),
    // but compare not more than left_block_size symbols (we have gt_in
    // to resolve the long comparisons).
    while (el < left_block_size && left_block_beg + i + el < text_length && txt[i + el] == pat[el])
      next(pat, ++el, s, p, r);

    if (left_block_beg + i + el != text_length &&
         ((el < left_block_size && txt[i + el] > pat[el]) ||
         (el == left_block_size && gt_in->get(i)))) gt_out->set(i);

    long j = i_max;
    if (el > el_max) {
      std::swap(el, el_max);
      std::swap(s, s_max);
      std::swap(p, p_max);
      std::swap(r, r_max);
      i_max = i;
    }

    if (p > 0L && 3L * p <= el && !memcmp(pat, pat + p, s)) {
      for (long k = 1L; k < std::min(range_end - left_block_beg - i, p); ++k)
        if (gt_out->get(j + k)) gt_out->set(i + k);

      i += p;
      el -= p;
    } else {
      long h = (el / 3L) + 1L;
      for (long k = 1L; k < std::min(range_end - left_block_beg - i, h); ++k)
        if (gt_out->get(j + k)) gt_out->set(i + k);

      i += h;
      el = 0;
      s = 0;
      p = 0;
    }
  }
}


//==============================================================================
// Compute the first left_block_size bits of gt_out bitvector. It is defined
// as follows. gt_out[i] == 1 iff the suffix of text starting at position
// left_block_beg + i if greater than suffix starting at position
// left_block_beg.
//
// We assume left_block_size and right_block_size are > 0 and that gt_in[i] == 1
// iff the suffix of text starting at right_block_beg + i is greater than suffix
// of text starting at right_block_beg, where right_block_beg = left_block_beg
// + left_block_size. gt_in is guaranteed to contain at least left_block_size
// bits, but also not more than right_block_size + 1. Only then we can always
// resolve the ties.
//
// max_threads is the intended number of threads used for computation.
//==============================================================================
void finalize_gt(unsigned char *text, long text_length, long left_block_beg,
    long left_block_size, bitvector *gt_in, bitvector *gt_out, long max_threads) {
  // Max range size has to be a multiple of 8. Otherwise two
  // threads might try to update bits in the same byte.
  long left_block_end = left_block_beg + left_block_size;
  long max_range_size = (left_block_size + max_threads - 1) / max_threads;
  while (max_range_size & 7) ++max_range_size;
  long n_blocks = (left_block_size + max_range_size - 1) / max_range_size;

  std::thread **threads = new std::thread*[n_blocks];
  for (long i = 0; i < n_blocks; ++i) {
    long range_beg = left_block_beg + i * max_range_size;
    long range_end = std::min(range_beg + max_range_size, left_block_end);
    
    threads[i] = new std::thread(finalize_gt_aux, text, text_length,
        left_block_beg, left_block_size, range_beg, range_end, gt_in, gt_out);
  }
  for (long i = 0; i < n_blocks; ++i) threads[i]->join();
  for (long i = 0; i < n_blocks; ++i) delete threads[i];
  delete[] threads;
}

#endif  // __INMEM_FINALIZE_GT_H
