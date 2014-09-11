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
#include "multifile_bitvector.h"
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


void compute_gt_end(unsigned char *block, long block_length,
                    unsigned char *prev_block, long prev_block_length,
                    long gt_length, multifile *gt_begin_multifile, bitvector *gt_end) {
  multifile_bitvector_reader gt_begin_reader(*gt_begin_multifile);

  long i = 0, el = 0, s = 0, p = 0, r = 0;
  long i_max = 0, el_max = 0, s_max = 0, p_max = 0, r_max = 0;
  while (i < block_length) {
    while (i + el < block_length && el < prev_block_length && block[i + el] == prev_block[el])
      next(prev_block, ++el, s, p, r);

    if ((el == prev_block_length && prev_block_length == gt_length) ||
        (i + el == block_length && (el == gt_length || (!gt_begin_reader.access(gt_length - el - 1)))) ||
        (i + el < block_length && block[i + el] > prev_block[el]))
      gt_end->set(i);

    long j = i_max;
    if (el > el_max) {
      std::swap(el, el_max);
      std::swap(s, s_max);
      std::swap(p, p_max);
      std::swap(r, r_max);
      i_max = i;
    }

    if (p && 3 * p <= el && !memcmp(prev_block, prev_block + p, s)) {
      for (long k = 1; k < std::min(p, block_length - i); ++k)
        if (gt_end->get(j + k)) gt_end->set(i + k);

      i += p;
      el -= p;
    } else {
      long h = (el / 3) + 1;
      for (long k = 1; k < std::min(h, block_length - i); ++k)
        if (gt_end->get(j + k)) gt_end->set(i + k);

      i += h;
      el = 0;
      s = 0;
      p = 0;
    }
  }
}


//==============================================================================
// Compute gt_begin for text.
//==============================================================================
long compute_gt_begin(unsigned char *text, long length, bitvector *gt_begin) {
  long whole_suffix_rank = length - 1;
  long i = 1, el = 0, s = 0, p = 0, r = 0;
  long i_max = 0, el_max = 0, s_max = 0, p_max = 0, r_max = 0;
  while (i < length) {
    while (i + el < length && el < length && text[i + el] == text[el])
      next(text, ++el, s, p, r);
 
    if (i + el < length && (el == length || text[i + el] > text[el])) {
      // To avoid reversing, fill gt_begin backwards.
      gt_begin->set(length - 1 - i);
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

    if (p && 3 * p <= el && !memcmp(text, text + p, s)) {
      for (long k = 1; k < std::min(length - i, p); ++k) { // Optimized as above.
        if (gt_begin->get(length - 1 - (j + k))) {
          gt_begin->set(length - 1 - (i + k));
          --whole_suffix_rank;
        }
      }

      i += p;
      el -= p;
    } else {
      long h = (el / 3) + 1;
      for (long k = 1; k < std::min(length - i, h); ++k) { // Optimized as above.
        if (gt_begin->get(length - 1 - (j + k))) {
          gt_begin->set(length - 1 - (i + k));
          --whole_suffix_rank;
        }
      }

      i += h;
      el = 0;
      s = 0;
      p = 0;
    }
  }

  return whole_suffix_rank;
}

#endif // __SRANK_H_INCLUDED
