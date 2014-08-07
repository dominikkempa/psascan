#ifndef __NAIVE_COMPUTE_PARTIAL_SA_H
#define __NAIVE_COMPUTE_PARTIAL_SA_H

#include "divsufsort.h"
#include "utils.h"

void naive_compute_partial_sa(unsigned char *text, long length,
    long block_beg, long block_end, int* &partial_sa) {

  long block_size = block_end - block_beg;
  fprintf(stderr, "    Allocating: ");
  long double start = utils::wclock();
  partial_sa = new int[block_size];
  int *sa = new int[length];
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
  
  fprintf(stderr, "    Divsufsort: ");
  start = utils::wclock();
  divsufsort(text, sa, (int)length);
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

  fprintf(stderr, "    Selecting suffixes: ");
  start = utils::wclock();
  for (long i = 0, ptr = 0; i < length; ++i)
    if (block_beg <= sa[i] && sa[i] < block_end)
      partial_sa[ptr++] = sa[i] - block_beg;
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

  delete[] sa;
}

void naive_select_partial_sa(long length,
    long block_beg, long block_end, int* &partial_sa,
    const char *sa_filename) {

  long block_size = block_end - block_beg;
  fprintf(stderr, "    Allocating: ");
  long double start = utils::wclock();
  partial_sa = new int[block_size];
  int *sa = new int[length];
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
  
  fprintf(stderr, "    Reading SA from file: ");
  start = utils::wclock();
  utils::read_objects_from_file(sa, length, sa_filename);
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

  fprintf(stderr, "    Selecting suffixes: ");
  start = utils::wclock();
  for (long i = 0, ptr = 0; i < length; ++i)
    if (block_beg <= sa[i] && sa[i] < block_end)
      partial_sa[ptr++] = sa[i] - block_beg;
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

  delete[] sa;
}


#endif  // __NAIVE_COMPUTE_PARTIAL_SA_H
