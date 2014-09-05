#ifndef __INITIAL_PARTIAL_SUFSORT_H
#define __INITIAL_PARTIAL_SUFSORT_H

#include <algorithm>
#include <thread>

#include "divsufsort_template.h"
#include "bitvector.h"


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
template<typename T>
void initial_partial_sufsort(unsigned char *text, long text_length,
    bitvector* &gt, T* &partial_sa, long max_block_size, long /*max_threads*/) {
  long double start = utils::wclock();
  long n_blocks = (text_length + max_block_size - 1) / max_block_size;

  // Rename the blocks in parallel.
  // XXX change this parallelism to vertical!
  if (n_blocks > 1) {
    fprintf(stderr, "    Renaming blocks: ");
    start = utils::wclock();
    std::thread **threads = new std::thread*[n_blocks - 1];
    for (long i = 0; i + 1 < n_blocks; ++i) {
      long block_beg = i * max_block_size;
      long block_end = std::min(block_beg + max_block_size, text_length);
      long block_size = block_end - block_beg;
      threads[i] = new std::thread(rename_block,
          text, block_beg, block_size, gt);
    }
    for (long i = 0; i + 1 < n_blocks; ++i) threads[i]->join();
    for (long i = 0; i + 1 < n_blocks; ++i) delete threads[i];
    delete[] threads;
    fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
  }

  // Now run the threads. This parallelism has to be horizontal.
  fprintf(stderr, "    Running divsufsort in parallel: ");
  start = utils::wclock();
  std::thread **threads = new std::thread*[n_blocks];
  for (long i = 0; i < n_blocks; ++i) {
    long block_beg = i * max_block_size;
    long block_end = std::min(block_beg + max_block_size, text_length);
    long block_size = block_end - block_beg;
    threads[i] = new std::thread(run_divsufsort<T>, text + block_beg,
        partial_sa + block_beg, (T)block_size);
  }
  for (long i = 0; i < n_blocks; ++i) threads[i]->join();
  for (long i = 0; i < n_blocks; ++i) delete threads[i];
  delete[] threads;
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

  // Finally, we restore the text.
  // XXX: change parallelism to vertical.
  if (n_blocks > 1) {
    fprintf(stderr, "    Rerenaming blocks: ");
    start = utils::wclock();
    threads = new std::thread*[n_blocks - 1];
    for (long i = 0; i + 1 < n_blocks; ++i) {
      long block_beg = i * max_block_size;
      long block_end = std::min(block_beg + max_block_size, text_length);
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
