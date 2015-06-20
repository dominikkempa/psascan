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
#include "multifile.h"
#include "multifile_bit_stream_reader.h"
#include "disk_pattern.h"


// Update ms-decomposition of T[0..n) from T[0..n-1).
void next(const unsigned char *T, long n, long &s, long &p, long &r) {
  if (n == 1) { s = 0; p = 1; r = 0; return; }
  long i = n - 1;
  while (i < n) {
    unsigned char a = T[s + r], b = T[i];
    if (a > b) { p = i - s + 1; r = 0; }
    else if (a < b) { i -= r; s = i; p = 1; r = 0; }
    else { ++r; if (r == p) r = 0; } ++i;
  }
}


// Update ms-decomposition of T[0..n) from T[0..n-1).
void next(pattern *T, long n, long &s, long &p, long &r) {
  if (n == 1) { s = 0; p = 1; r = 0; return; }
  long i = n - 1;
  while (i < n) {
    unsigned char a = (*T)[s + r], b = (*T)[i];
    if (a > b) { p = i - s + 1; r = 0; }
    else if (a < b) { i -= r; s = i; p = 1; r = 0; }
    else { ++r; if (r == p) r = 0; } ++i;
  }
}


void compute_block_gt_end(const unsigned char *block, long block_beg, long block_end,
    long text_length, std::string text_filename,
    const multifile *tail_gt_begin_reversed, bitvector *block_gt_end) {
  long block_size = block_end - block_beg;

  multifile_bit_stream_reader tail_gt_begin_reversed_reader(tail_gt_begin_reversed);
  pattern pat(text_filename, block_end);
  long pat_length = text_length - block_end;

  long i = 0, el = 0;
  while (i < block_size) {
    while (i + el < block_size && el < pat_length && block[i + el] == pat[el])
      ++el;

    if (el == pat_length ||
        (i + el == block_size && (!tail_gt_begin_reversed_reader.access(pat_length - el))) ||
        (i + el < block_size && block[i + el] > pat[el]))
      block_gt_end->set(i);

    el = 0;
    ++i;
  }
}


// Inplace transformation.
void compute_block_gt_begin_reversed_from_block_gt_end(const unsigned char *block,
    long block_beg, long block_end, bitvector *gt) {
  long block_size = block_end - block_beg;
    
  gt->flip(0);
  long i = 1, el = 0;
  
  while (i < block_size) {
    while (i + el < block_size && block[i + el] == block[el]) ++el;

    if ((i + el == block_size && !(gt->get(block_size - i))) ||
        (i + el < block_size && block[i + el] > block[el])) gt->set(block_size - i);
    else gt->reset(block_size - i);

    ++i;
    el = 0;
  }
}


//==============================================================================
// Compute (reversed) gt_begin for the given block. Returns the number of suffixes
// starting inside the block that are smaller than the suffix starting at
// the beginning of the block (which means that the maximal value is
// block_size - 1.
//==============================================================================
long compute_block_gt_begin_reversed(const unsigned char *block, long block_beg,
    long block_end, long text_length, std::string text_filename,
    const multifile *tail_gt_begin_reversed, bitvector *block_gt_begin_reversed) {
  long block_size = block_end - block_beg;

  compute_block_gt_end(block, block_beg, block_end, text_length, text_filename,
                       tail_gt_begin_reversed, block_gt_begin_reversed);
  compute_block_gt_begin_reversed_from_block_gt_end(block, block_beg,
                                                    block_end,
                                                    block_gt_begin_reversed);

  long result = block_size - 1;
  for (long i = 1; i < block_size; ++i)
    if (block_gt_begin_reversed->get(i)) --result;

  return result;
}


#endif  // __SRANK_H_INCLUDED
