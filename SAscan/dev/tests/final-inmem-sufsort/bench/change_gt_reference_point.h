#ifndef __CHANGE_GT_REFERENCE_POINT_H
#define __CHANGE_GT_REFERENCE_POINT_H

#include <cstring>
#include <algorithm>
#include <thread>

#include "bitvector.h"
#include "srank_aux.h"


//==============================================================================
// Compute range [microblock_beg..microblock_end) of bits in the output
// bitvector gt_out.
//==============================================================================
void change_gt_reference_point_aux(unsigned char *text, long text_length,
    long block_beg, long block_end, long microblock_beg, long microblock_end,
    bitvector *gt_in, bitvector *gt_out) {
  // Void comparing the suffix with itself, the bit in gt_out
  // should be 0 anyway, so we don't have ot do anything.
  if (block_beg == microblock_beg) ++microblock_beg;

  unsigned char *pat = text + block_beg, *txt = pat;
  long i = microblock_beg - block_beg, el = 0L, s = 0L, p = 0L, r = 0L;
  long i_max = i, el_max = 0L, s_max = 0L, p_max = 0L, r_max = 0L;
  while (i < microblock_end - block_beg) {
    // Compute lcp(text[left_block_beg..), text[left_block_beg+i..),
    // but compare not more than left_block_size symbols (we have gt_in
    // to resolve the long comparisons).
    while (block_beg + i + el < block_end && txt[i + el] == pat[el])
      next(pat, ++el, s, p, r);

    if (block_beg + i + el != text_length &&
        ((block_beg + i + el != block_end && txt[i + el] > pat[el]) ||
         (block_beg + i + el == block_end && !gt_in->get(el)))) gt_out->set(i);

    long j = i_max;
    if (el > el_max) {
      std::swap(el, el_max);
      std::swap(s, s_max);
      std::swap(p, p_max);
      std::swap(r, r_max);
      i_max = i;
    }

    if (p > 0L && 3L * p <= el && !memcmp(pat, pat + p, s)) {
      for (long k = 1L; k < std::min(microblock_end - block_beg - i, p); ++k)
        if (gt_out->get(j + k)) gt_out->set(i + k);

      i += p;
      el -= p;
    } else {
      long h = (el / 3L) + 1L;
      for (long k = 1L; k < std::min(microblock_end - block_beg - i, h); ++k)
        if (gt_out->get(j + k)) gt_out->set(i + k);

      i += h;
      el = 0;
      s = 0;
      p = 0;
    }
  }
}


//==============================================================================
// Given bitvector gt_in of length block_size such that for i = 0, ..,
// block_size - 1, gt_in[i] == 1 iff text[block_beg + i..) > text[block_end..)
// (where block_size = block_end - block_beg), compute bivctor gt_out of length
// block_size + 1such that for i = 0, .., block_size, gt_out[i] == 1 iff
// text[block_beg + i..) > text[block_beg..). Of course gt_out[0] is always 0.
//
// In other words, function takes the bitvector gt computed with respect to
// block end and from it computes the gt bitvector with respect to block start.
//==============================================================================
void change_gt_reference_point(unsigned char *text, long text_length,
    long block_beg, long block_end, bitvector *gt_in, bitvector* &gt_out,
    long max_threads) {
  // As usual, we make sure max_block_size is a multiple of 8 so that
  // no two threads will attempt to modify the bits inside the same byte.
  long block_size = block_end - block_beg;
  long max_microblock_size = (block_size + max_threads - 1) / max_threads;
  while (max_microblock_size & 7) ++max_microblock_size;
  long n_microblocks = (block_size + max_microblock_size - 1) / max_microblock_size;

  //----------------------------------------------------------------------------
  // STEP 1: allocate and zero-initialize (in parallel) the bitvector.
  //         Also, copy its last bit from the gt_in.
  //----------------------------------------------------------------------------
  gt_out = new bitvector(block_size + 1, max_threads);
  if (!gt_in->get(0)) gt_out->set(block_size);


  //----------------------------------------------------------------------------
  // STEP 2: compute the remaining block_size bits.
  //----------------------------------------------------------------------------
  std::thread **threads = new std::thread*[n_microblocks];
  for (long i = 0; i < n_microblocks; ++i) {
    long microblock_beg = block_beg + i * max_microblock_size;
    long microblock_end = std::min(microblock_beg + max_microblock_size, block_end);

    threads[i] = new std::thread(change_gt_reference_point_aux,
        text, text_length, block_beg, block_end, microblock_beg,
        microblock_end, gt_in, gt_out);
  }

  for (long i = 0; i < n_microblocks; ++i) threads[i]->join();
  for (long i = 0; i < n_microblocks; ++i) delete threads[i];
  delete[] threads;
}

#endif  // __CHANGE_GT_REFERENCE_POINT_H
