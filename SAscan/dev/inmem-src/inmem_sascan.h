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
#include "bwtsa.h"
#include "parallel_shrink.h"
#include "skewed-merge.h"


//==============================================================================
// We assume sa is already allocated.
//
// What should be done is some kind of optimization to the initial gt
// bitvectors computation, so that the last bitvector is not computed, and
// also that the last block is not renamed and so on. The end result should
// be, that when used with one thread, the procedure simply runs divsufsort
// and there is no overhead of running this function over running divsufsort.
//==============================================================================
template<typename saidx_t, unsigned pagesize_log = 12>
void inmem_sascan(unsigned char *text, long text_length, unsigned char *sa_bwt,
    long max_threads = 1, bool compute_bwt = false, bool compute_gt_begin = false,
    bitvector *gt_begin = NULL, long max_blocks = -1) {
  static const unsigned pagesize = (1U << pagesize_log);

  long double start;
  if (max_blocks == -1)
    max_blocks = max_threads;

  if (!compute_gt_begin) {
    if (gt_begin) {
      fprintf(stderr, "Error: check gt_begin == NULL failed\n");
      std::exit(EXIT_FAILURE);
    }
    gt_begin = new bitvector(text_length, max_threads);
  } else {
    if (!gt_begin) {
      fprintf(stderr, "inmem_sascan: gt_begin was requested but is not allocated!\n");
      std::exit(EXIT_FAILURE);
    }
  }

  // long max_block_size = (text_length + max_blocks - 1) / max_blocks;
  // while ((max_block_size & 7) || (max_block_size & pagesize_mask)) ++max_block_size;
  // long n_blocks = (text_length + max_block_size - 1) / max_block_size;

  //----------------------------------------------------------------------------
  // min_block_size must be a multiple ot alignment unit. Alignement unit is to
  // simplify the computation involving bitvectors and page arrays. Note: it
  // may happen then min_block_size > text_length. This is perfectly fine due
  // to the way we compute block boundaries (always separate if for the last
  // block).
  //----------------------------------------------------------------------------

  long alignment_unit = (long)std::max(pagesize, 8U);
  long min_block_size = text_length / max_blocks;
  while (min_block_size & (alignment_unit - 1)) --min_block_size;
  if (!min_block_size) min_block_size = text_length;
  long n_blocks = text_length / min_block_size;

  fprintf(stderr, "Text length = %ld (%.2LfMiB)\n", text_length, text_length / (1024.L * 1024));
  fprintf(stderr, "Min block size = %ld (%.2LfMiB)\n", min_block_size, min_block_size / (1024.L * 1024));
  fprintf(stderr, "Max blocks = %ld\n", max_blocks);
  fprintf(stderr, "Number of blocks = %ld\n", n_blocks);
  fprintf(stderr, "Max threads = %ld\n", max_threads);
  fprintf(stderr, "sizeof(saidx_t) = %lu\n", sizeof(saidx_t));
  fprintf(stderr, "pagesize = %u\n", (1U << pagesize_log));
  fprintf(stderr, "compute bwt = %s\n", compute_bwt ? "true" : "false");
  fprintf(stderr, "compute gt_begin = %s\n", compute_gt_begin ? "true" : "false");
  fprintf(stderr, "\n");

  bwtsa_t<saidx_t> *bwtsa = (bwtsa_t<saidx_t> *)sa_bwt;

  //----------------------------------------------------------------------------
  // STEP 1: compute initial bitvectors, and partial suffix arrays.
  //----------------------------------------------------------------------------
  if (compute_gt_begin || n_blocks > 1) {
    fprintf(stderr, "Compute initial bitvectors:\n");
    start = utils::wclock();
    compute_initial_gt_bitvectors(text, text_length, gt_begin, min_block_size, max_threads);
    fprintf(stderr, "Time: %.2Lf\n\n", utils::wclock() - start);
  }

  fprintf(stderr, "Initial sufsort:\n");
  start = utils::wclock();
  initial_partial_sufsort(text, text_length, gt_begin, bwtsa, min_block_size, max_threads);
  fprintf(stderr, "Time: %.2Lf\n\n", utils::wclock() - start);


  //----------------------------------------------------------------------------
  // STEP 2: compute the gt bitvectors for blocks that will be on the right
  //         side during the merging. Also, create block description array.
  //----------------------------------------------------------------------------
  if (compute_gt_begin || n_blocks > 1) {
    fprintf(stderr, "Overwriting gt_end with gt_begin: ");
    start = utils::wclock();
    gt_end_to_gt_begin(text, text_length, gt_begin, min_block_size, max_threads);
    fprintf(stderr, "%.2Lf\n\n", utils::wclock() - start);
  }

  float rl_ratio = 10.L; // estimated empirically
  int max_left_size = std::max(1, (int)floor(n_blocks * (2.L - (2.125L + sizeof(saidx_t)) / 5.L)));
  fprintf(stderr, "Assumed rl_ratio: %.2f\n", rl_ratio);
  fprintf(stderr, "Max left size = %d\n", max_left_size);
  fprintf(stderr, "Peak memory usage during last merging = %.3Lfn\n",
      (2.125L + sizeof(saidx_t)) + (5.L * max_left_size) / n_blocks);
  MergeSchedule schedule(n_blocks, rl_ratio, max_left_size);

  fprintf(stderr, "Skewed merge schedule:\n");
  print_schedule(schedule, n_blocks);
  fprintf(stderr, "\n");

  long i0;
  pagearray<bwtsa_t<saidx_t>, pagesize_log> *result = NULL;
  if (n_blocks > 1 || compute_bwt) {
    result = balanced_merge<saidx_t, pagesize_log>(text, text_length, bwtsa,
        gt_begin, min_block_size, 0, n_blocks, max_threads, compute_gt_begin, i0, schedule);
  }

  if (n_blocks > 1) {
    // We permute it to plain array.
    fprintf(stderr, "\nPermuting the resulting SA to plain array: ");
    start = utils::wclock();
    result->permute_to_plain_array(max_threads);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
  }

  delete result;
  if (!compute_gt_begin && n_blocks > 1) {
    delete gt_begin;
    gt_begin = NULL;
  }

  unsigned char *bwt = NULL;
  if (compute_bwt) {
    // Allocate aux, copy bwt into aux.
    fprintf(stderr, "Copying bwtsa.bwt into aux memory: ");
    start = utils::wclock();
    bwt = (unsigned char *)malloc(text_length);
    parallel_copy<bwtsa_t<saidx_t>, unsigned char>(bwtsa, bwt, text_length, max_threads);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
  }

  fprintf(stderr, "Shrinking bwtsa.sa into sa: ");
  start = utils::wclock();

  parallel_shrink<bwtsa_t<saidx_t>, saidx_t>(bwtsa, text_length, max_threads);

  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

  if (compute_bwt) {
    // Copy from aux into the end of bwtsa.
    fprintf(stderr, "Copying bwt from aux memory to the end of bwtsa: ");
    start = utils::wclock();
    unsigned char *dest = (unsigned char *)(((saidx_t *)bwtsa) + text_length);
    parallel_copy<unsigned char, unsigned char>(bwt, dest, text_length, max_threads);
    free(bwt);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
  }
}


#endif  // __FINAL_INMEM_SUFSORT_H
