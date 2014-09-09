#ifndef __PARALLEL_COPY_H_INCLUDED
#define __PARALLEL_COPY_H_INCLUDED

#include <algorithm>
#include <thread>

#include "bwtsa.h"
#include "uint40.h"

template<typename T, typename S>
void parallel_copy_aux(T *src, S *dest, long length) {
  for (long i = 0; i < length; ++i)
    dest[i] = (S)src[i];
}

// Specilization
template<>
void parallel_copy_aux(bwtsa_t<uint40> *src, unsigned char *dest, long length) {
  for (long i = 0; i < length; ++i)
    dest[i] = src[i].bwt;
}

// Specilization
template<>
void parallel_copy_aux(bwtsa_t<int> *src, unsigned char *dest, long length) {
  for (long i = 0; i < length; ++i)
    dest[i] = src[i].bwt;
}


// Conversion from T to S has to make sense.
template<typename T, typename S>
void parallel_copy(T *src, S *dest, long length, long max_threads) {
  long max_block_size = (length + max_threads - 1) / max_threads;
  long n_blocks = (length + max_block_size - 1) / max_block_size;

  std::thread **threads = new std::thread*[n_blocks];
  for (long i = 0; i < n_blocks; ++i) {
    long block_beg = i * max_block_size;
    long block_end = std::min(block_beg + max_block_size, length);
    long block_size = block_end - block_beg;

    threads[i] = new std::thread(parallel_copy_aux<T, S>,
        src + block_beg, dest + block_beg, block_size);
  }

  for (long i = 0; i < n_blocks; ++i) threads[i]->join();
  for (long i = 0; i < n_blocks; ++i) delete threads[i];
  delete[] threads;
}

// Specialization
template<>
void parallel_copy(bwtsa_t<uint40> *src, unsigned char *dest, long length, long max_threads) {
  long max_block_size = (length + max_threads - 1) / max_threads;
  long n_blocks = (length + max_block_size - 1) / max_block_size;

  std::thread **threads = new std::thread*[n_blocks];
  for (long i = 0; i < n_blocks; ++i) {
    long block_beg = i * max_block_size;
    long block_end = std::min(block_beg + max_block_size, length);
    long block_size = block_end - block_beg;

    threads[i] = new std::thread(parallel_copy_aux<bwtsa_t<uint40>, unsigned char>,
        src + block_beg, dest + block_beg, block_size);
  }

  for (long i = 0; i < n_blocks; ++i) threads[i]->join();
  for (long i = 0; i < n_blocks; ++i) delete threads[i];
  delete[] threads;
}

// Specialization
template<>
void parallel_copy(bwtsa_t<int> *src, unsigned char *dest, long length, long max_threads) {
  long max_block_size = (length + max_threads - 1) / max_threads;
  long n_blocks = (length + max_block_size - 1) / max_block_size;

  std::thread **threads = new std::thread*[n_blocks];
  for (long i = 0; i < n_blocks; ++i) {
    long block_beg = i * max_block_size;
    long block_end = std::min(block_beg + max_block_size, length);
    long block_size = block_end - block_beg;

    threads[i] = new std::thread(parallel_copy_aux<bwtsa_t<int>, unsigned char>,
        src + block_beg, dest + block_beg, block_size);
  }

  for (long i = 0; i < n_blocks; ++i) threads[i]->join();
  for (long i = 0; i < n_blocks; ++i) delete threads[i];
  delete[] threads;
}



#endif  // __PARALLEL_SHRINK_H_INCLUDED
