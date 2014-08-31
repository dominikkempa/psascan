#ifndef __PARALLEL_UTILS_H
#define __PARALLEL_UTILS_H

#include <algorithm>
#include <thread>

#include "utils.h"


namespace parallel_utils {

template<typename T>
void fill_aux(T *tab, long length, T val) {
  for (long i = 0; i < length; ++i)
    tab[i] = val;
}

template<typename T>
void fill(T *tab, long length, T val, long max_threads) {
  long max_block_size = (length + max_threads - 1) / max_threads;
  long n_blocks = (length + max_block_size - 1) / max_block_size;

  std::thread **threads = new std::thread*[n_blocks];
  for (long i = 0; i < n_blocks; ++i) {
    long block_beg = i * max_block_size;
    long block_end = std::min(block_beg + max_block_size, length);
    long block_size = block_end - block_beg;

    threads[i] = new std::thread(fill_aux<T>, tab + block_beg,
        block_size, val);
  }

  for (long i = 0; i < n_blocks; ++i) threads[i]->join();
  for (long i = 0; i < n_blocks; ++i) delete threads[i];
  delete[] threads;
}

} // namespace parallel_utils

#endif  // __PARALLEL_UTILS_H
