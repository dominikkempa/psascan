#ifndef __BALANCED_MERGE_H
#define __BALANCED_MERGE_H

#include <vector>

#include "bitvector.h"
#include "inmem_gap_array.h"
#include "inmem_compute_gap.h"
#include "inmem_bwt_from_sa.h"
#include "parallel_merge.h"


// Returns i0, or -1 if unknown.
template<typename T>
long balanced_merge(unsigned char *text, long text_length, T *sa, bitvector *gt,
    long max_block_size, long range_beg, long range_end, long max_threads,
    unsigned char *bwt, bool need_gt) {
  long range_size = range_end - range_beg;

  if (range_size == 1) {
    long block_beg = max_block_size * range_beg;
    long block_end = std::min(block_beg + max_block_size, text_length);
    long block_size = block_end - block_beg;
    return bwt_from_sa_into_dest<T>(sa + block_beg, text + block_beg, block_size, bwt + block_beg, max_threads);
  }

  // Split blocks almost equally into two groups.
  long lrange_size = (range_size + 1) / 2;
  long rrange_size = range_size - lrange_size;

  long lrange_beg = range_beg;
  long lrange_end = range_beg + lrange_size;
  long rrange_beg = lrange_end;
  long rrange_end = rrange_beg + rrange_size;

  long lbeg = lrange_beg * max_block_size;
  long rbeg = rrange_beg * max_block_size;
  long lend = rbeg;
  long rend = std::min(rbeg + rrange_size * max_block_size, text_length);
  long lsize = lend - lbeg;
  long rsize = rend - rbeg;

  // Compute partial sa for the blocks recursively.
  long left_i0  = balanced_merge(text, text_length, sa, gt, max_block_size,
      lrange_beg, lrange_end, max_threads, bwt, need_gt);
  // right
  long right_i0 = balanced_merge(text, text_length, sa, gt, max_block_size, rrange_beg,
      rrange_end, max_threads, bwt, true);

  // Merge partial suffix arrays.
  fprintf(stderr, "Merging blocks [%ld..%ld) and [%ld..%ld):\n", lbeg, lend, rbeg, rend);
  long double start = utils::wclock();

  inmem_gap_array *gap;
  fprintf(stderr, "  Computing gap:\n");
  long double start1 = utils::wclock();
  inmem_compute_gap(text, text_length, lbeg, lsize, rsize, sa + lbeg, bwt + lbeg, gt,
      gap, max_threads, need_gt, left_i0, (1L << 21));
  fprintf(stderr, "  Time: %.2Lf\n", utils::wclock() - start1);

  fprintf(stderr, "  Merging partial SAs: ");
  start1 = utils::wclock();
  long delta_i0 = merge<T, 12>(sa + lbeg, lsize, rsize, gap, max_threads, left_i0);
  fprintf(stderr, "total: %.2Lf\n", utils::wclock() - start1);

  fprintf(stderr, "  Merging BWTs: ");
  start1 = utils::wclock();
  bwt[rbeg + right_i0] = text[rbeg - 1];
  merge<unsigned char, 12>(bwt + lbeg, lsize, rsize, gap, max_threads, -1);
  fprintf(stderr, "total: %.2Lf\n", utils::wclock() - start1);

  fprintf(stderr, "  Deleting gap: ");
  start1 = utils::wclock();
  delete gap;
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start1);

  fprintf(stderr, "Time: %.2Lf\n", utils::wclock() - start);
  return left_i0 + delta_i0;
}

#endif  // __BALANCED_MERGE_H
