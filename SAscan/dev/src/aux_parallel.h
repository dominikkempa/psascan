/*******************************************************************************
  Implementation of two methods:

  * bwt_from_sa_replace_text
    Given partial sa and text, compute the bwt and replace the text. The
    algorithm does not use any extra space and is fully parallelized.
    Destroys the suffix array!

  * bwt_from_sa_into_dest
    Given partial sa and text, compute bwt and write into dest, where dest is
    any place in memory different from sa and text. The algorithm is fully
    parallelized.
*******************************************************************************/

#ifndef __AUX_PARALLEL_H_INCLUDED
#define __AUX_PARALLEL_H_INCLUDED

#include <algorithm>
#include <thread>

#include "utils.h"


template<typename T>
void find_index_0(const T *tab, long beg, long end, long *answer) {
  *answer = -1;
  for (long j = beg; j < end; ++j)
    if (tab[j] == 0) *answer = j;
}

template<typename T>
void bwt_of_range0(const T *tab, const unsigned char *text, long beg, long end,
    unsigned char *dest) {
  for (long j = beg; j < end; ++j)
    dest[j] = text[tab[j] - 1];
}

//==============================================================================
// Compute in parallel for i = range_beg, .., range_end - 1:
//    dest[i - beg] := text[sa[i] - 1]
//==============================================================================
template<typename T>
void bwt_of_range(const T *sa, const unsigned char *text, long range_beg, long range_end,
    unsigned char *dest, long max_threads) {
  long range_size = range_end - range_beg;
  if (!range_size) return;
  long block_size = (range_size + max_threads - 1) / max_threads;
  long n_blocks = (range_size + block_size - 1) / block_size;

  std::thread **threads = new std::thread*[n_blocks];
  for (long j = 0; j < n_blocks; ++j) {
    long block_beg = range_beg + j * block_size;
    long block_end = std::min(block_beg + block_size, range_end);
    threads[j] = new std::thread(bwt_of_range0<T>,
        sa, text, block_beg, block_end, dest + block_beg);
  }
  for (long j = 0; j < n_blocks; ++j) threads[j]->join();
  for (long j = 0; j < n_blocks; ++j) delete threads[j];
  delete[] threads;
}

//==============================================================================
// Computation equivalent to:
//
//   for (long i = 0, ptr = 0; i < length; ++i)
//     if (SA[i] != 0) dest[ptr++] = text[SA[i] - 1];
//
// We assume length > 0.
//==============================================================================
template<typename T>
void bwt_from_sa_into_dest(const T *sa, const unsigned char *text, long length,
    unsigned char *dest, long max_threads) {
  // First, we find j such that SA[j] = 0 in parallel.
  long block_size = (length + max_threads - 1) / max_threads;
  long n_blocks = (length + block_size - 1) / block_size;

  long *index_0 = new long[n_blocks];
  std::thread **threads = new std::thread*[n_blocks];
  for (long j = 0; j < n_blocks; ++j) {
    long block_beg = j * block_size;
    long block_end = std::min(block_beg + block_size, length);
    threads[j] = new std::thread(find_index_0<T>, sa,
        block_beg, block_end, index_0 + j);
  }
  for (long j = 0; j < n_blocks; ++j) threads[j]->join();
  for (long j = 0; j < n_blocks; ++j) delete threads[j];
  delete[] threads;

  // Retreive the answer.
  long SA0 = 0;
  for (long j = 0; j < n_blocks; ++j)
    if (index_0[j] != -1) SA0 = index_0[j];
  delete[] index_0;

  // Compute BWT in parallel.
  bwt_of_range(sa, text, 0, SA0, dest, max_threads);
  bwt_of_range(sa, text, SA0 + 1, length, dest + SA0, max_threads);
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

template<typename T>
void replace_sa_with_bwt(T *sa, const unsigned char *text, long block_beg,
    long block_end, long *index_0) {
  *index_0 = -1;
  for (long i = block_beg; i < block_end; ++i) {
    long sai = sa[i];
    if (sai == 0) {
      *index_0 = i;
      sa[i] = 0;
    }
    else sa[i] = text[sai - 1];
  }
}

template<typename T>
void move_bwt_from_sa_to_dest_aux(const T *sa, unsigned char *dest, long length) {
  for (long i = 0; i < length; ++i)
    dest[i] = (unsigned char)sa[i];
}

//==============================================================================
// In parallel do:
//    for i = 0, .., length - 1: dest[i] = (unsigned char)sa[i];
//==============================================================================
template<typename T>
void move_bwt_from_sa_to_dest(const T *sa, long length, unsigned char *dest,
    long max_threads) {
  if (!length) return;
  long block_size = (length + max_threads - 1) / max_threads;
  long n_blocks = (length + block_size - 1) / block_size;

  std::thread **threads = new std::thread*[n_blocks];
  for (long i = 0; i < n_blocks; ++i) {
    long block_beg = i * block_size;
    long block_end = std::min(block_beg + block_size, length);
    long this_block_size = block_end - block_beg;
    threads[i] = new std::thread(move_bwt_from_sa_to_dest_aux<T>,
        sa + block_beg, dest + block_beg, this_block_size);
  }
  for (long i = 0; i < n_blocks; ++i) threads[i]->join();
  for (long i = 0; i < n_blocks; ++i) delete threads[i];
  delete[] threads;
}

template<typename T>
long bwt_from_sa_replace_text(T *SA, unsigned char *text, long length,
    long max_threads) {
  //----------------------------------------------------------------------------
  // STEP 1: for all i replace SA[i] with B[SA[i] - 1]. All i except i0
  //         such that SA[i] = 0. Such i is returned via index_0 from the
  //         thread that finds it.
  //----------------------------------------------------------------------------
  long block_size = (length + max_threads - 1) / max_threads;
  long n_blocks = (length + block_size - 1) / block_size;

  long *index_0 = new long[n_blocks];
  std::thread **threads = new std::thread*[n_blocks];
  for (long i = 0; i < n_blocks; ++i) {
    long block_beg = i * block_size;
    long block_end = std::min(block_beg + block_size, length);

    threads[i] = new std::thread(replace_sa_with_bwt<T>,
        SA, text, block_beg, block_end, index_0 + i);
  }
  for (long i = 0; i < n_blocks; ++i) threads[i]->join();
  for (long i = 0; i < n_blocks; ++i) delete threads[i];
  delete[] threads;

  // Find i0.
  long i0 = 0;
  for (long i = 0; i < n_blocks; ++i)
    if (index_0[i] != -1) i0 = index_0[i];
  delete[] index_0;

  //----------------------------------------------------------------------------
  // STEP 2: overwrite the text with BWT which now is temporarily stored
  //         inside SA.
  //----------------------------------------------------------------------------
  move_bwt_from_sa_to_dest(SA, length, text, max_threads);
  return i0;
}

#endif  // __AUX_PARALLEL_H_INCLUDED
