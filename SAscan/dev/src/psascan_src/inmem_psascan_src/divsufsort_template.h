#ifndef __PSASCAN_SRC_INMEM_PSASCAN_SRC_DIVSUFSORT_TEMPLATE_H_INCLUDED
#define __PSASCAN_SRC_INMEM_PSASCAN_SRC_DIVSUFSORT_TEMPLATE_H_INCLUDED

#include <cstdio>
#include <cstdlib>

#include "divsufsort.h"
#include "divsufsort64.h"


namespace psascan_private {
namespace inmem_psascan_private {

template<typename T>
void run_divsufsort(const unsigned char *, T*, T) {
  fprintf(stderr, "\ndivsufsort: non-standard call. Use either"
      "int or long for second and third argument.\n");
  std::exit(EXIT_FAILURE);
}

template<>
void run_divsufsort(const unsigned char *text, int *sa, int length) {
  divsufsort(text, sa, length);
}

template<>
void run_divsufsort(const unsigned char *text, long *sa, long length) {
  divsufsort64(text, sa, length);
}

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __PSASCAN_SRC_INMEM_PSASCAN_SRC_DIVSUFSORT_TEMPLATE_H_INCLUDED
