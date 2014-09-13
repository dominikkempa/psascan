#ifndef __INITIAL_PARTIAL_SUFSORT_H
#define __INITIAL_PARTIAL_SUFSORT_H

#include <algorithm>
#include <thread>

#include "divsufsort_template.h"
#include "bitvector.h"
#include "bwtsa.h"
#include "parallel_shrink.h"
#include "parallel_expand.h"
#include "parallel_copy.h"


//==============================================================================
// Rename the given block using its gt bitvector.
//==============================================================================
void rename_block(unsigned char *text, long block_beg, long block_length,
    bitvector *gt) {
  long block_end = block_beg + block_length;

  unsigned char *block = text + block_beg;
  unsigned char last = block[block_length - 1];
  for (long i = 0; i + 1 < block_length; ++i)
    if (block[i] > last || (block[i] == last && gt->get(block_end - i - 2)))
      ++block[i];
  ++block[block_length - 1];
}


//==============================================================================
// Re-rename block back to original.
//==============================================================================
void rerename_block(unsigned char *block, long block_length) {
  unsigned char last = block[block_length - 1] - 1;
  for (long i = 0; i < block_length; ++i)
    if (block[i] > last) --block[i];
}


//==============================================================================
// Given gt bitvectors, compute partial suffix arrays of blocks.
// To do this, in parallel:
//   1) rename the blocks
//   2) run divsufsort on each block
//==============================================================================

template<typename saidx_t>
void initial_partial_sufsort(unsigned char *, long, bitvector* &, bwtsa_t<saidx_t> *, long, long) {
  fprintf(stderr, "Error: initial_partial_sufsort: given saidx_t is "
      "not supported, sizeof(saidx_t) = %ld\n", (long)sizeof(saidx_t));
  std::exit(EXIT_FAILURE);
}

template<>
void initial_partial_sufsort(unsigned char *text, long text_length,
    bitvector* &gt, bwtsa_t<uint40> *bwtsa, long min_block_size, long max_threads) {
  long double start = utils::wclock();
  long n_blocks = text_length / min_block_size;

  //----------------------------------------------------------------------------
  // STEP 1: Rename the blocks in parallel.
  //----------------------------------------------------------------------------
  // XXX change this parallelism to vertical!
  if (n_blocks > 1) {
    fprintf(stderr, "  Renaming blocks: ");
    start = utils::wclock();
    std::thread **threads = new std::thread*[n_blocks - 1];
    for (long i = 0; i + 1 < n_blocks; ++i) {
      long block_beg = i * min_block_size;
      long block_end = block_beg + min_block_size;
      long block_size = block_end - block_beg;

      threads[i] = new std::thread(rename_block,
          text, block_beg, block_size, gt);
    }

    for (long i = 0; i + 1 < n_blocks; ++i) threads[i]->join();
    for (long i = 0; i + 1 < n_blocks; ++i) delete threads[i];
    delete[] threads;

    fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
  }

  long max_block_size = text_length - (n_blocks - 1) * min_block_size;
  if (max_block_size >= (2L << 30)) {  // Use 64-bit divsufsort.


    fprintf(stderr, "Not implemented yet (time saving, this will never come up in experiments).\n");
    std::exit(EXIT_FAILURE);

    /*long *temp_space = (long *)malloc(text_length * sizeof(long));

    fprintf(stderr, "  Running divsufsort64 in parallel: ");
    start = utils::wclock();
    std::thread **threads = new std::thread*[n_blocks];
    for (long i = 0; i < n_blocks; ++i) {
      long block_beg = i * min_block_size;
      long block_end = block_beg + min_block_size;

      if (i == n_blocks - 1) block_end = text_length;
      long block_size = block_end - block_beg;

      threads[i] = new std::thread(run_divsufsort<long>, text + block_beg, temp_space + block_beg, block_size);
    }

    for (long i = 0; i < text_length; ++i)
      bwtsa[i].sa = (uint40)(temp_space[i]);

    free(temp_space);*/

#if 0
    //----------------------------------------------------------------------------
    // STEP 2: Compute and allocate extra space.
    //----------------------------------------------------------------------------
    long *temp_sa = (long *)bwtsa;
    long available_mem = (sizeof(uint40) + 1) * text_length;

    long accomodated_blocks = available_mem / (max_block_size * sizeof(long));
    long extra_blocks = n_blocks - accomodated_blocks;
    fprintf(stderr, "  Allocating memory:\n");
    fprintf(stderr, "    Accomodated blocks = %ld\n", accomodated_blocks);
    fprintf(stderr, "    Need extra space for %ld blocks\n", extra_blocks);
    long extra_integers_needed = text_length - std::min(text_length, accomodated_blocks * max_block_size);

    long *extra_mem = NULL;
    if (extra_integers_needed)
      extra_mem = (long *)malloc(sizeof(long) * extra_integers_needed);

    long double total_space_usage = available_mem +
      sizeof(long) * extra_integers_needed + text_length * 1.125L;
    fprintf(stderr, "    Total space usage: %.2Lfn\n",
        (long double)total_space_usage / text_length);


    //----------------------------------------------------------------------------
    // STEP 3: Run the threads. This parallelism has to be horizontal.
    //----------------------------------------------------------------------------
    fprintf(stderr, "  Running divsufsort64 in parallel: ");
    start = utils::wclock();
    std::thread **threads = new std::thread*[n_blocks];
    for (long i = 0; i < n_blocks; ++i) {
      long block_beg = i * max_block_size;
      long block_end = std::min(block_beg + max_block_size, text_length);
      long block_size = block_end - block_beg;

      if (i < accomodated_blocks)
        threads[i] = new std::thread(run_divsufsort<long>, text + block_beg, temp_sa + block_beg, block_size);
      else {
        long block_beg2 = block_beg - accomodated_blocks * max_block_size;
        threads[i] = new std::thread(run_divsufsort<long>, text + block_beg, extra_mem + block_beg2, block_size);
      }
    }

    for (long i = 0; i < n_blocks; ++i) threads[i]->join();
    for (long i = 0; i < n_blocks; ++i) delete threads[i];
    delete[] threads;
    fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

    fprintf(stderr, "  Shrinking 64-bit integers to bwtsa objects: ");
    start = utils::wclock();
    parallel_shrink<long, bwtsa_t<uint40> >(temp_sa, accomodated_blocks * max_block_size, max_threads);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

    fprintf(stderr, "  Moving data from extra space to bwtsa: ");
    start = utils::wclock();
    parallel_copy<long, bwtsa_t<uint40> >(extra_mem,
        bwtsa + accomodated_blocks * max_block_size, extra_integers_needed, max_threads);
    if (extra_mem) free(extra_mem);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
#endif
  } else {  // Use 32-bit divsufsort.
    int *temp_sa = (int *)bwtsa;

    //----------------------------------------------------------------------------
    // STEP 3: Run the threads. This parallelism has to be horizontal.
    //----------------------------------------------------------------------------
    fprintf(stderr, "  Running divsufsort32 in parallel: ");
    start = utils::wclock();
    std::thread **threads = new std::thread*[n_blocks];
    for (long i = 0; i < n_blocks; ++i) {
      long block_beg = i * min_block_size;
      long block_end = block_beg + min_block_size;

      if (i == n_blocks - 1) block_end = text_length;
      long block_size = block_end - block_beg;

      threads[i] = new std::thread(run_divsufsort<int>, text + block_beg, temp_sa + block_beg, block_size);
    }

    for (long i = 0; i < n_blocks; ++i) threads[i]->join();
    for (long i = 0; i < n_blocks; ++i) delete threads[i];
    delete[] threads;

    fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

    fprintf(stderr, "  Expanding 32-bit integers to bwtsa objects: ");
    start = utils::wclock();
    parallel_expand<int, bwtsa_t<uint40> >(temp_sa, text_length, max_threads);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
  }


  //----------------------------------------------------------------------------
  // STEP 4: Finally, we restore the text.
  //----------------------------------------------------------------------------
  // XXX: change parallelism to vertical.
  if (n_blocks > 1) {
    fprintf(stderr, "  Rerenaming blocks: ");
    start = utils::wclock();
    std::thread **threads = new std::thread*[n_blocks - 1];
    for (long i = 0; i + 1 < n_blocks; ++i) {
      long block_beg = i * min_block_size;
      long block_end = block_beg + min_block_size;
      long block_size = block_end - block_beg;

      threads[i] = new std::thread(rerename_block,
          text + block_beg, block_size);
    }

    for (long i = 0; i + 1 < n_blocks; ++i) threads[i]->join();
    for (long i = 0; i + 1 < n_blocks; ++i) delete threads[i];
    delete[] threads;

    fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
  }
}





template<>
void initial_partial_sufsort(unsigned char *text, long text_length,
    bitvector* &gt, bwtsa_t<int> *bwtsa, long min_block_size, long max_threads) {
  long double start = utils::wclock();
  long n_blocks = text_length / min_block_size;

  //----------------------------------------------------------------------------
  // STEP 1: Rename the blocks in parallel.
  //----------------------------------------------------------------------------
  // XXX change this parallelism to vertical!
  if (n_blocks > 1) {
    fprintf(stderr, "  Renaming blocks: ");
    start = utils::wclock();
    std::thread **threads = new std::thread*[n_blocks - 1];
    for (long i = 0; i + 1 < n_blocks; ++i) {
      long block_beg = i * min_block_size;
      long block_end = block_beg + min_block_size;
      long block_size = block_end - block_beg;

      threads[i] = new std::thread(rename_block,
          text, block_beg, block_size, gt);
    }

    for (long i = 0; i + 1 < n_blocks; ++i) threads[i]->join();
    for (long i = 0; i + 1 < n_blocks; ++i) delete threads[i];
    delete[] threads;

    fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
  }
  
  int *temp_sa = (int *)bwtsa;

  //----------------------------------------------------------------------------
  // STEP 2: Run the threads. This parallelism has to be horizontal.
  //----------------------------------------------------------------------------
  fprintf(stderr, "  Running divsufsort32 in parallel: ");
  start = utils::wclock();
  std::thread **threads = new std::thread*[n_blocks];
  for (long i = 0; i < n_blocks; ++i) {
    long block_beg = i * min_block_size;
    long block_end = block_beg + min_block_size;

    if (i == n_blocks - 1) block_end = text_length;
    long block_size = block_end - block_beg;

    threads[i] = new std::thread(run_divsufsort<int>, text + block_beg, temp_sa + block_beg, block_size);
  }

  for (long i = 0; i < n_blocks; ++i) threads[i]->join();
  for (long i = 0; i < n_blocks; ++i) delete threads[i];
  delete[] threads;

  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

  fprintf(stderr, "  Expanding 32-bit integers to bwtsa objects: ");
  start = utils::wclock();
  parallel_expand<int, bwtsa_t<int> >(temp_sa, text_length, max_threads);
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);


  //----------------------------------------------------------------------------
  // STEP 3: Finally, we restore the text.
  //----------------------------------------------------------------------------
  // XXX: change parallelism to vertical.
  if (n_blocks > 1) {
    fprintf(stderr, "  Rerenaming blocks: ");
    start = utils::wclock();
    threads = new std::thread*[n_blocks - 1];
    for (long i = 0; i + 1 < n_blocks; ++i) {
      long block_beg = i * min_block_size;
      long block_end = block_beg + min_block_size;
      long block_size = block_end - block_beg;

      threads[i] = new std::thread(rerename_block,
          text + block_beg, block_size);
    }

    for (long i = 0; i + 1 < n_blocks; ++i) threads[i]->join();
    for (long i = 0; i + 1 < n_blocks; ++i) delete threads[i];
    delete[] threads;

    fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
  }
}





#endif  // __INITIAL_PARTIAL_SUFSORT_H
