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


struct block_description {
  block_description(long beg = 0, long end = 0) {
    m_beg = beg;
    m_end = end;
  }

  long m_beg, m_end;
};


//==============================================================================
// We assume sa is already allocated.
//
// What should be done is some kind of optimization to the initial gt
// bitvectors computation, so that the last bitvector is not computed, and
// also that the last block is not renamed and so on. The end result should
// be, that when used with one thread, the procedure simply runs divsufsort
// and there is no overhead of running this function over running divsufsort.
//==============================================================================
template<typename T>
void inmem_sascan(unsigned char *text, long text_length, T* sa,
    long max_threads = 1, long max_blocks = -1) {
  long double start;
  if (max_blocks == -1)
    max_blocks = max_threads;

  //----------------------------------------------------------------------------
  // STEP 1: compute initial bitvectors, and partial suffix arrays.
  //----------------------------------------------------------------------------
  fprintf(stderr, "Compute initial bitvectors:\n");
  bitvector *gt;
  start = utils::wclock();
  compute_initial_gt_bitvectors(text, text_length, gt, max_blocks, max_threads);
  fprintf(stderr, "Time: %.2Lf\n", utils::wclock() - start);

  fprintf(stderr, "Initial sufsort:\n");
  start = utils::wclock();
  initial_partial_sufsort(text, text_length, gt, sa, max_blocks);
  fprintf(stderr, "Time: %.2Lf\n", utils::wclock() - start);


  //----------------------------------------------------------------------------
  // STEP 2: compute the gt bitvectors for blocks that will be on the right
  //         side during the merging. Also, create block description array.
  //----------------------------------------------------------------------------
  long max_block_size = (text_length + max_blocks - 1) / max_blocks;
  while (max_block_size & 7) ++max_block_size;
  long n_blocks = (text_length + max_block_size - 1) / max_block_size;

  std::vector<block_description> block_desc;
  fprintf(stderr, "Overwriting gt_end with gt_begin: ");
  start = utils::wclock();
  gt_end_to_gt_begin(text, text_length, gt, max_blocks);
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

  // Create initial blocks.
  for (long i = 0; i < n_blocks; ++i) {
    long block_beg = i * max_block_size;
    long block_end = std::min(block_beg + max_block_size, text_length);
    block_desc.push_back(block_description(block_beg, block_end));
  }

  //----------------------------------------------------------------------------
  // STEP 3: keep merging blocks until we have only one left.
  //----------------------------------------------------------------------------
  while (n_blocks != 1) {
    std::vector<block_description> new_block_desc;

    for (long i = 0; i < n_blocks; i += 2) {
      // Check if block has no right neighbour to merge
      // with. If yes, promote to the next merging roung.
      if (i + 1 == n_blocks) {
        new_block_desc.push_back(block_desc[i]);
        break;
      }

      //------------------------------------------------------------------------
      // Merge blocks with description in block_desc[i] and block_desc[i + 1].
      //------------------------------------------------------------------------
      fprintf(stderr, "Merging blocks [%ld..%ld) and [%ld..%ld):\n",
          block_desc[i].m_beg, block_desc[i].m_end,
          block_desc[i + 1].m_beg, block_desc[i + 1].m_end);
      start = utils::wclock();

      // 1
      //
      // Compute some initial parameters.
      long lbeg = block_desc[i].m_beg;
      long rend = block_desc[i + 1].m_end;
      long lsize = block_desc[i].m_end - lbeg;
      long rsize = rend - block_desc[i + 1].m_beg;

      // 2
      //
      // Compute gap and new gt bitvector.
      inmem_gap_array *gap;
      fprintf(stderr, "  Computing gap:\n");
      long double start1 = utils::wclock();
      inmem_compute_gap(text, text_length, lbeg, lsize, rsize,
          sa + lbeg, gt, gap, max_threads, (1L << 21));
      fprintf(stderr, "    Time: %.2Lf\n", utils::wclock() - start1);

      // 3
      //
      // Merge (in place) partial suffix arrays.
      fprintf(stderr, "  Merging: ");
      start1 = utils::wclock();
      merge<T, 12>(sa + lbeg, lsize, rsize, gap, max_threads);
      fprintf(stderr, "total: %.2Lf\n", utils::wclock() - start1);

      fprintf(stderr, "  Deleting gap: ");
      start1 = utils::wclock();
      delete gap;
      fprintf(stderr, "%.2Lf\n", utils::wclock() - start1);

      fprintf(stderr, "Time: %.2Lf\n", utils::wclock() - start);

      // 4
      //
      // Promote the merged block to the next round.
      new_block_desc.push_back(block_description(lbeg, rend));
    }

    block_desc = new_block_desc;
    n_blocks = (long)block_desc.size();
  }

  delete gt;
}


#endif  // __FINAL_INMEM_SUFSORT_H
