#ifndef __INMEM_SASCAN_BALANCED_MERGE_H_INCLUDED
#define __INMEM_SASCAN_BALANCED_MERGE_H_INCLUDED

#include <vector>

#include "../../bitvector.h"
#include "inmem_gap_array.h"
#include "inmem_compute_gap.h"
#include "inmem_bwt_from_sa.h"
#include "parallel_merge.h"
#include "pagearray.h"
#include "bwtsa.h"
#include "skewed-merge.h"
#include "../../multifile.h"


namespace inmem_sascan_private {

template<typename saidx_t, unsigned pagesize_log>
pagearray<bwtsa_t<saidx_t>, pagesize_log> *balanced_merge(const unsigned char *text,
    long text_length, bwtsa_t<saidx_t> *bwtsa, bitvector *gt,
    long max_block_size, long range_beg, long range_end, long max_threads,
    bool need_gt, bool need_bwt, long &result_i0, MergeSchedule &schedule,
    long text_beg, long text_end, long supertext_length,
    std::string supertext_filename, const multifile *tail_gt_begin_reversed,
    long *i0_array, long **block_rank_matrix) {
  typedef pagearray<bwtsa_t<saidx_t>, pagesize_log> pagearray_type;

  long shift = (max_block_size - text_length % max_block_size) % max_block_size;
  long range_size = range_end - range_beg;

  if (range_size == 1) {
    long block_beg = range_beg * max_block_size;
    long block_end = block_beg + max_block_size;
    block_beg = std::max(0L, block_beg - shift);
    block_end -= shift;

    result_i0 = i0_array[range_beg];
    pagearray_type *bwtsa_pagearray = new pagearray_type(bwtsa + block_beg, bwtsa + block_end);
    return bwtsa_pagearray;
  }


  //----------------------------------------------------------------------------
  // STEP 1: Split the blocks in the left and right group.
  //----------------------------------------------------------------------------
  long lrange_size = schedule.left_size(range_size);
  long rrange_size = range_size - lrange_size;

  long lrange_beg = range_beg;
  long lrange_end = range_beg + lrange_size;
  long rrange_beg = lrange_end;
  long rrange_end = rrange_beg + rrange_size;

  long lbeg = lrange_beg * max_block_size;
  long rbeg = rrange_beg * max_block_size;
  long lend = rbeg;
  long rend = rbeg + rrange_size * max_block_size;
  lbeg = std::max(0L, lbeg - shift);
  rbeg -= shift;
  lend -= shift;
  rend -= shift;

  long lsize = lend - lbeg;
  long rsize = rend - rbeg;


  //----------------------------------------------------------------------------
  // STEP 2: Compute partial SAs and BWTs for left and right block.
  //----------------------------------------------------------------------------

  // 1
  //
  // left block
  long left_i0;
  pagearray_type *l_bwtsa =
    balanced_merge<saidx_t, pagesize_log>(text, text_length, bwtsa, gt,
    max_block_size, lrange_beg, lrange_end, max_threads, need_gt, true, left_i0, schedule,
    text_beg, text_end, supertext_length, supertext_filename, tail_gt_begin_reversed, i0_array,
    block_rank_matrix);

  // 2
  // 
  // right block
  long right_i0;
  pagearray_type *r_bwtsa =
      balanced_merge<saidx_t, pagesize_log>(text, text_length, bwtsa, gt,
      max_block_size, rrange_beg, rrange_end, max_threads, true, need_bwt, right_i0, schedule,
      text_beg, text_end, supertext_length, supertext_filename, tail_gt_begin_reversed, i0_array,
      block_rank_matrix);


  //----------------------------------------------------------------------------
  // STEP 3: Merge partial suffix arrays.
  //----------------------------------------------------------------------------
  fprintf(stderr, "Merging blocks %ld-%ld with %ld-%ld\n", lrange_beg + 1, lrange_end, rrange_beg + 1, rrange_end);
  long double start = utils::wclock();


  // 1
  //
  // Compute gap
  fprintf(stderr, "  Computing gap:\n");
  inmem_gap_array *gap;
  long double rank_init_time;
  long double streaming_time;
  long double start1 = utils::wclock();
  inmem_compute_gap<saidx_t, pagesize_log>(text, text_length, lbeg, lsize, rsize,
      *l_bwtsa, gt, gap, max_threads, need_gt, left_i0, (1L << 21), rank_init_time, streaming_time,
      block_rank_matrix, lrange_beg, lrange_size, rrange_size);
  fprintf(stderr, "  Time: %.2Lf\n", utils::wclock() - start1);


  // 2
  //
  // Merge partial SAs and BWTs
  fprintf(stderr, "  Merging SA/BWT:  ");
  start1 = utils::wclock();
  long delta_i0;
  if (need_bwt)
    (*r_bwtsa)[right_i0].bwt = text[rbeg - 1];
  pagearray_type *result = parallel_merge(l_bwtsa, r_bwtsa, gap, max_threads, left_i0, delta_i0, lsize);
  result_i0 = left_i0 + delta_i0;
  long double merging_time = utils::wclock() - start1;
  fprintf(stderr, "total: %.2Lf\n", merging_time);


  // 3
  //
  // Clean up.
  start1 = utils::wclock();
  delete l_bwtsa;
  delete r_bwtsa;
  delete gap;
  long double cleaning_time = utils::wclock() - start1;
  if (cleaning_time > 0.2L)
    fprintf(stderr, "Cleaning: %.2Lf\n", cleaning_time);

  long double time_per_elem_left = merging_time / (lsize + rsize) + rank_init_time / lsize;
  long double time_per_elem_right = merging_time / (lsize + rsize) + streaming_time / rsize;
  long double ratio = time_per_elem_right / time_per_elem_left;
  fprintf(stderr, "Time: %.2Lf (rl_ratio = %.3Lf)\n",
      utils::wclock() - start, ratio);

  return result;
}

}  // namespace inmem_sascan

#endif  // __INMEM_SASCAN_BALANCED_MERGE_H_INCLUDED
