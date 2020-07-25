//==============================================================================
// Implementation of string range matching, as descrined in
//
//   Juha Kärkkäinen, Dominik Kempa, Simon J. Puglisi:
//   String Range Matching. In Proc. CPM 2014
//=============================================================================

#ifndef __SRANK_H_INCLUDED
#define __SRANK_H_INCLUDED

#include <cstring>
#include <algorithm>

#include "bitvector.h"
#include "multifile_bitvector.h"
#include "disk_pattern.h"


//==============================================================================
// Update ms-decomposition of text[0..length) from text[0..length - 1).
//==============================================================================
void next(unsigned char *text, long length, long &s, long &p, long &r) {
  if (length == 1) { s = 0; p = 1; r = 0; return; }
  long i = length - 1;
  while (i < length) {
    unsigned char a = text[s + r], b = text[i];
    if (a > b) { p = i - s + 1; r = 0; }
    else if (a < b) { i -= r; s = i; p = 1; r = 0; }
    else { ++r; if (r == p) r = 0; } ++i;
  }
}



void compute_sm_end(unsigned char *block, long block_beg, long block_end,
    long text_length, std::string text_filename,
    multifile *tail_gt_begin_reversed, bitvector *block_sm_end) {
  long block_size = block_end - block_beg;

  multifile_bitvector_reader tail_gt_begin_reversed_reader(tail_gt_begin_reversed);
  pattern pat(text_filename, block_end);
  long pat_length = text_length - block_end;

  long i = 0, el = 0;
  while (i < block_size) {
    while (i + el < block_size && el < pat_length && block[i + el] == pat[el])
      ++el;

    if ((i + el < block_size && el < pat_length && block[i + el] < pat[el]) ||
        (i + el == block_size && tail_gt_begin_reversed_reader.access(pat_length - el)))
      block_sm_end->set(i);

    el = 0;
    ++i;
  }
}


//==============================================================================
// Given we are gt_end for the block -> very well defined and all.
// We want to compute gt_begin from it, and to do it in-place we want to
// have gt_begin reversed in fact.
//==============================================================================
void transform_sm_end_into_gt_begin_reversed(
    unsigned char *text, long length, bitvector *bv) {
  long i = 1, el = 0;
  while (i < length) {
    while (i + el < length && text[i + el] == text[el]) ++el;
    if (i + el < length) {
      if (text[i + el] > text[el]) bv->set(length - i);
      else bv->reset(length - i);
    }
    ++i;
    el = 0;
  }
}

#endif // __SRANK_H_INCLUDED
