#ifndef __BALANCED_MERGE_H
#define __BALANCED_MERGE_H

#include <vector>

#include "bitvector.h"
#include "inmem_gap_array.h"
#include "inmem_compute_gap.h"
#include "inmem_bwt_from_sa.h"
#include "parallel_merge.h"
#include "pagearray.h"


// Returns i0, or -1 if unknown.
template<typename T, unsigned pagesize_log>
pagearray<T, pagesize_log> *balanced_merge(unsigned char *text, long text_length, T *sa,
    bitvector *gt, long max_block_size, long range_beg, long range_end, long max_threads,
    unsigned char *bwt, bool need_gt, long &result_i0) {
  typedef pagearray<T, pagesize_log> pagearray_type;
  long range_size = range_end - range_beg;

  if (range_size == 1) {
    long block_beg = max_block_size * range_beg;
    long block_end = std::min(block_beg + max_block_size, text_length);
    long block_size = block_end - block_beg;

    fprintf(stderr, "Computing BWT for block [%ld..%ld): ", block_beg, block_end);
    long double start = utils::wclock();
    bwt_from_sa_into_dest<T>(sa + block_beg, text + block_beg, block_size, bwt + block_beg, max_threads, result_i0);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
    return new pagearray_type(sa + block_beg, sa + block_beg + block_size);
  }

  // Split blocks almost equally into two groups.
  long lrange_size = /*(range_size + 1) / 2*/range_size - 1;
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
  typedef pagearray<T, pagesize_log> pagearray_type;

  long left_i0;
  pagearray_type *l_pagearray = balanced_merge<T, pagesize_log>(text, text_length, sa, gt,
      max_block_size, lrange_beg, lrange_end, max_threads, bwt, need_gt, left_i0);

  // right
  long right_i0;
  pagearray_type *r_pagearray = balanced_merge<T, pagesize_log>(text, text_length, sa, gt,
      max_block_size, rrange_beg, rrange_end, max_threads, bwt, true, right_i0);

  // Merge partial suffix arrays.
  fprintf(stderr, "Merging blocks [%ld..%ld) and [%ld..%ld):\n", lbeg, lend, rbeg, rend);
  long double start = utils::wclock();

  inmem_gap_array *gap;
  fprintf(stderr, "  Computing gap:\n");
  long double start1 = utils::wclock();
  inmem_compute_gap(text, text_length, lbeg, lsize, rsize, l_pagearray,
      bwt + lbeg, gt, gap, max_threads, need_gt, left_i0, (1L << 21));
  fprintf(stderr, "  Time: %.2Lf\n", utils::wclock() - start1);

  fprintf(stderr, "  Merging partial SAs:  ");
  start1 = utils::wclock();
  long delta_i0;
  pagearray_type *merged_pagearray = parallel_merge(l_pagearray, r_pagearray, gap, max_threads, left_i0, lsize, delta_i0);
  delete l_pagearray;
  delete r_pagearray;

  fprintf(stderr, "total: %.2Lf\n", utils::wclock() - start1);

  fprintf(stderr, "  Merging partial BWTs: ");
  start1 = utils::wclock();
  bwt[rbeg + right_i0] = text[rbeg - 1];
  merge<unsigned char, pagesize_log>(bwt + lbeg, lsize, rsize, gap, max_threads, -1, 0);

  fprintf(stderr, "total: %.2Lf\n", utils::wclock() - start1);

  fprintf(stderr, "  Deleting gap: ");
  start1 = utils::wclock();
  delete gap;
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start1);

  fprintf(stderr, "Time: %.2Lf\n", utils::wclock() - start);
  result_i0 = left_i0 + delta_i0;

  return merged_pagearray;
}

#endif  // __BALANCED_MERGE_H
