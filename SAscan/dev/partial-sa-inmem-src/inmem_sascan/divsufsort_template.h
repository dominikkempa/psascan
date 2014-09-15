#ifndef __INMEM_SASCAN_DIVSUFSORT_TEMPLATE_H
#define __INMEM_SASCAN_DIVSUFSORT_TEMPLATE_H

#include <cstdio>
#include <cstdlib>

#include "divsufsort.h"
#include "divsufsort64.h"

namespace inmem_sascan_private {


template<typename T>
void run_divsufsort(unsigned char *, T*, T) {
  fprintf(stderr, "\ndivsufsort: non-standard call. Use either"
      "int or long for second and third argument.\n");
  std::exit(EXIT_FAILURE);
}

template<>
void run_divsufsort(unsigned char *text, int *sa, int length) {
  divsufsort(text, sa, length);
}

template<>
void run_divsufsort(unsigned char *text, long *sa, long length) {
  divsufsort64(text, sa, length);
}

}  // namespace inmem_sascan

#endif  // __DIVSUFSORT_TEMPLATE_H
