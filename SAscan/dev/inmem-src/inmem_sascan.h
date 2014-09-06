#ifndef __FINAL_INMEM_SUFSORT_H
#define __FINAL_INMEM_SUFSORT_H

#include <vector>

#include "bitvector.h"
#include "inmem_gap_array.h"
#include "compute_initial_gt_bitvectors.h"
#include "initial_partial_sufsort.h"
#include "change_gt_reference_point.h"
#include "inmem_compute_gap.h"
#include "parallel_merge.h"
#include "balanced_merge.h"
#include "pagearray.h"


//==============================================================================
// We assume sa is already allocated.
//
// What should be done is some kind of optimization to the initial gt
// bitvectors computation, so that the last bitvector is not computed, and
// also that the last block is not renamed and so on. The end result should
// be, that when used with one thread, the procedure simply runs divsufsort
// and there is no overhead of running this function over running divsufsort.
//==============================================================================
template<typename T, unsigned pagesize_log = 12>
void inmem_sascan(unsigned char *text, long text_length, T* sa,
    long max_threads = 1, long max_blocks = -1) {
  static const unsigned pagesize_mask = (1U << pagesize_log) - 1;
  long double start;
  if (max_blocks == -1)
    max_blocks = max_threads;

  long max_block_size = (text_length + max_blocks - 1) / max_blocks;
  while ((max_block_size & 7) || (max_block_size & pagesize_mask)) ++max_block_size;
  long n_blocks = (text_length + max_block_size - 1) / max_block_size;

  unsigned char *bwt = (unsigned char *)malloc(text_length);

  //----------------------------------------------------------------------------
  // STEP 1: compute initial bitvectors, and partial suffix arrays.
  //----------------------------------------------------------------------------
  bitvector *gt = NULL;
  if (n_blocks > 1) {
    fprintf(stderr, "Compute initial bitvectors:\n");
    start = utils::wclock();
    compute_initial_gt_bitvectors(text, text_length, gt, max_block_size, max_threads);
    fprintf(stderr, "Time: %.2Lf\n\n", utils::wclock() - start);
  }

  fprintf(stderr, "Initial sufsort:\n");
  start = utils::wclock();
  initial_partial_sufsort(text, text_length, gt, sa, max_block_size, max_threads);
  fprintf(stderr, "Time: %.2Lf\n\n", utils::wclock() - start);


  //----------------------------------------------------------------------------
  // STEP 2: compute the gt bitvectors for blocks that will be on the right
  //         side during the merging. Also, create block description array.
  //----------------------------------------------------------------------------
  if (n_blocks > 1) {
    fprintf(stderr, "Overwriting gt_end with gt_begin: ");
    start = utils::wclock();
    gt_end_to_gt_begin(text, text_length, gt, max_block_size, max_threads);
    fprintf(stderr, "%.2Lf\n\n", utils::wclock() - start);

    typedef pagearray<T, pagesize_log> pagearray_type;
    long i0;
    pagearray_type *result = 
      balanced_merge<T, pagesize_log>(text, text_length, sa, gt, max_block_size, 0, n_blocks, max_threads, bwt, false, i0);

    // We permute it to plain array.
    fprintf(stderr, "\nPermuting the resulting SA to plain array: ");
    start = utils::wclock();
    result->permute_to_plain_array(max_threads);
    delete result;
    fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
  }

  delete gt;
  free(bwt);
}


#endif  // __FINAL_INMEM_SUFSORT_H
