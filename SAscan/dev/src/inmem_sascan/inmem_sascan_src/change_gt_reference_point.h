#ifndef __INMEM_SASCAN_CHANGE_GT_REFERENCE_POINT_H
#define __INMEM_SASCAN_CHANGE_GT_REFERENCE_POINT_H

#include <cstring>
#include <algorithm>
#include <thread>

#include "../../bitvector.h"
#include "srank_aux.h"

namespace inmem_sascan_private {

//==============================================================================
// Compute range [microblock_beg..microblock_end) of bits in the output
// bitvector gt_out.
//==============================================================================
void gt_end_to_gt_begin_aux(unsigned char *text,
    long block_beg, long block_end, bitvector *gt) {
  long block_size = block_end - block_beg;
  unsigned char *pat = text + block_beg, *txt = pat;

  long i = 1, el = 0L, s = 0L, p = 0L, r = 0L;
  long i_max = i, el_max = 0L, s_max = 0L, p_max = 0L, r_max = 0L;

  while (i < block_size) {
    // Compute lcp(text[left_block_beg..), text[left_block_beg+i..),
    // but compare not more than left_block_size symbols (we have gt
    // to resolve the long comparisons).
    while (block_beg + i + el < block_end && txt[i + el] == pat[el])
      next(pat, ++el, s, p, r);

    if (((block_beg + i + el != block_end && txt[i + el] > pat[el]) ||
         (block_beg + i + el == block_end && !gt->get(block_beg + i - 1))))
      gt->set(block_beg + i - 1);
    else gt->reset(block_beg + i - 1);

    long j = i_max;
    if (el > el_max) {
      std::swap(el, el_max);
      std::swap(s, s_max);
      std::swap(p, p_max);
      std::swap(r, r_max);
      i_max = i;
    }

    if (el < 100) {
      ++i;
      el = 0;
    } else if (p > 0L && (p << 2) <= el && !memcmp(pat, pat + p, s)) {
      for (long k = 1L; k < std::min(block_size - i, p); ++k) {
        if (gt->get(block_beg + j + k - 1)) gt->set(block_beg + i + k - 1);
        else gt->reset(block_beg + i + k - 1);
      }

      i += p;
      el -= p;
    } else {
      long h = (el >> 2) + 1L;
      for (long k = 1L; k < std::min(block_size - i, h); ++k) {
        if (gt->get(block_beg + j + k - 1)) gt->set(block_beg + i + k - 1);
        else gt->reset(block_beg + i + k - 1);
      }

      i += h;
      el = 0;
      s = 0;
      p = 0;
    }
  }
}


//==============================================================================
// Change gt_end bitvector into gt_begin using string range matching.
//==============================================================================
void gt_end_to_gt_begin(unsigned char *text, long text_length,
    bitvector *gt, long min_block_size, long /*max_threads*/) {
  // NOTE: The called of this function *must* use the same max block size.
  long n_blocks = text_length / min_block_size;


  //----------------------------------------------------------------------------
  // STEP 1: Compute the last bit in every block.
  //----------------------------------------------------------------------------
  for (long i = 0; i < n_blocks; ++i) {
    long block_beg = i * min_block_size;
    long block_end = block_beg + min_block_size;
    if (block_end + min_block_size > text_length) block_end = text_length;

    gt->flip(block_end - 1);
  }


  //----------------------------------------------------------------------------
  // STEP 2: compute remaining bits in every block.
  //----------------------------------------------------------------------------
  std::thread **threads = new std::thread*[n_blocks];
  for (long i = 0; i < n_blocks; ++i) {
    long block_beg = i * min_block_size;
    long block_end = block_beg + min_block_size;
    if (block_end + min_block_size > text_length) block_end = text_length;

    threads[i] = new std::thread(gt_end_to_gt_begin_aux,
        text, block_beg, block_end, gt);
  }

  for (long i = 0; i < n_blocks; ++i) threads[i]->join();
  for (long i = 0; i < n_blocks; ++i) delete threads[i];
  delete[] threads;
}

}  // namespace inmem_sascan

#endif  // __CHANGE_GT_REFERENCE_POINT_H
