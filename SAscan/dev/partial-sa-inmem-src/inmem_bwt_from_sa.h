////////////////////////////////////////////////////////////////////////////////
//  Implementation of two methods:
//
//  * bwt_from_sa_replace_text
//    Given partial sa and text, compute the bwt and replace the text. The
//    algorithm does not use any extra space and is fully parallelized.
//    Destroys the suffix array!
//
//  * bwt_from_sa_into_dest
//    Given partial sa and text, compute bwt and write into dest, where dest is
//    any place in memory different from sa and text. The algorithm is fully
//    parallelized.
//==============================================================================

#ifndef __INMEM_BWT_FROM_SA_H_INCLUDED
#define __INMEM_BWT_FROM_SA_H_INCLUDED

#include <algorithm>
#include <thread>

#include "utils.h"
#include "bwtsa.h"


template<typename saidx_t>
void bwt_from_sa_into_dest_aux(unsigned char *text, long beg, long end,
    bwtsa_t<saidx_t> *dest, long *i0) {
  *i0 = -1;
  for (long j = beg; j < end; ++j)
    if (dest[j].sa) dest[j].bwt = text[dest[j].sa - 1];
    else { dest[j].bwt = 0; *i0 = j; }
}


template<typename saidx_t>
void bwt_from_sa_into_dest(unsigned char *text, long length,
  bwtsa_t<saidx_t> *dest, long max_threads, long &result) {
  long max_block_size = (length + max_threads - 1) / max_threads;
  long n_blocks = (length + max_block_size - 1) / max_block_size;
  long *index_0 = new long[n_blocks];

  // Compute bwt and find i0, where sa[i0] == 0.
  std::thread **threads = new std::thread*[n_blocks];
  for (long i = 0; i < n_blocks; ++i) {
    long block_beg = i * max_block_size;
    long block_end = std::min(block_beg + max_block_size, length);

    threads[i] = new std::thread(bwt_from_sa_into_dest_aux<saidx_t>,
        text, block_beg, block_end, dest, index_0 + i);
  }

  for (long i = 0; i < n_blocks; ++i) threads[i]->join();
  for (long i = 0; i < n_blocks; ++i) delete threads[i];
  delete[] threads;

  // Find and return i0.
  result = -1;
  for (long i = 0; i < n_blocks; ++i)
    if (index_0[i] != -1) result = index_0[i];
  delete[] index_0;
}


#endif  // __INMEM_BWT_FROM_SA_INCLUDED
