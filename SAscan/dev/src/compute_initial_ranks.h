#ifndef __COMPUTE_INITIAL_RANKS_H_INCLUDED
#define __COMPUTE_INITIAL_RANKS_H_INCLUDED

#include <string>
#include <thread>
#include <algorithm>
#include <vector>

#include "smaller_suffixes.h"


template<typename saidx_t>
void compute_initial_ranks(unsigned char *block, long block_beg,
    long block_end, long text_length, saidx_t *block_partial_sa,
    std::string text_filename, std::vector<long> &result,
    long max_threads, long tail_begin, long tail_end) {
  long tail_length = tail_end - tail_begin;
  long stream_max_block_size = (tail_length + max_threads - 1) / max_threads;
  long n_threads = (tail_length + stream_max_block_size - 1) / stream_max_block_size;

  result.resize(n_threads);
  std::thread **threads = new std::thread*[n_threads];
  for (int t = 0; t < n_threads; ++t) {
    long stream_block_beg = tail_begin + t * stream_max_block_size;
    long stream_block_end = std::min(stream_block_beg + stream_max_block_size, tail_end);

    threads[t] = new std::thread(parallel_smaller_suffixes2<saidx_t>, block,
        block_beg, block_end, text_length, block_partial_sa, text_filename,
        stream_block_end, std::ref(result[t]));
  }

  for (int t = 0; t < n_threads; ++t) threads[t]->join();
  for (int t = 0; t < n_threads; ++t) delete threads[t];
  delete[] threads;
}


#endif  // __COMPUTE_INITIAL_RANKS_H_INCLUDED
