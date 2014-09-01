#ifndef __BALANCED_MERGE_H
#define __BALANCED_MERGE_H

#include <vector>

#include "bitvector.h"
#include "inmem_gap_array.h"
#include "inmem_compute_gap.h"
#include "parallel_merge.h"


template<typename T>
void balanced_merge(unsigned char *text, long text_length, T *sa, bitvector *gt,
    long max_block_size, long range_beg, long range_end, long max_threads, bool need_gt) {
  long range_size = range_end - range_beg;
  if (range_size == 1) return;

  // Split blocks almost equally into two groups.
  long lrange_size = (range_size + 1) / 2;
  long rrange_size = range_size - lrange_size;

  long lrange_beg = range_beg;
  long lrange_end = range_beg + lrange_size;
  long rrange_beg = lrange_end;
  long rrange_end = rrange_beg + rrange_size;


  // Compute partial sa for the blocks recursively.
  balanced_merge(text, text_length, sa, gt, max_block_size, lrange_beg, lrange_end, max_threads, need_gt);
  balanced_merge(text, text_length, sa, gt, max_block_size, rrange_beg, rrange_end, max_threads, true);


  // Merge partial suffix arrays.
  long lbeg = lrange_beg * max_block_size;
  long rbeg = rrange_beg * max_block_size;
  long lend = rbeg;
  long rend = std::min(rbeg + rrange_size * max_block_size, text_length);
  long lsize = lend - lbeg;
  long rsize = rend - rbeg;

  fprintf(stderr, "Merging blocks [%ld..%ld) and [%ld..%ld):\n", lbeg, lend, rbeg, rend);
  long double start = utils::wclock();

  inmem_gap_array *gap;
  fprintf(stderr, "  Computing gap:\n");
  long double start1 = utils::wclock();
  inmem_compute_gap(text, text_length, lbeg, lsize, rsize, sa + lbeg, gt, gap, max_threads, need_gt, (1L << 21));
  fprintf(stderr, "  Time: %.2Lf\n", utils::wclock() - start1);

  fprintf(stderr, "  Merging: ");
  start1 = utils::wclock();
  merge<T, 12>(sa + lbeg, lsize, rsize, gap, max_threads);
  fprintf(stderr, "total: %.2Lf\n", utils::wclock() - start1);

  fprintf(stderr, "  Deleting gap: ");
  start1 = utils::wclock();
  delete gap;
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start1);

  fprintf(stderr, "Time: %.2Lf\n", utils::wclock() - start);
}

#endif  // __BALANCED_MERGE_H
