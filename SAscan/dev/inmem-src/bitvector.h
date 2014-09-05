// Simple bitvector class.
#ifndef __BITVECTOR_H_INCLUDED
#define __BITVECTOR_H_INCLUDED

#include <cstdio>
#include <cstdlib>

#include "utils.h"

struct bitvector {
  bitvector(long length, long) {
    if (length <= 0) {
      fprintf(stderr, "Error: attempint to construct "
          "empty bitvector.\n");
      std::exit(EXIT_FAILURE);
    }

    m_alloc_bytes = (length + 7) / 8;
    m_data = (unsigned char *)calloc(m_alloc_bytes, sizeof(unsigned char));
  }

  inline bool get(long i) const {
    return m_data[i >> 3] & (1 << (i & 7));
  }

  inline void set(long i) {
    m_data[i >> 3] |= (1 << (i & 7));
  }

  inline void reset(long i) {
    m_data[i >> 3] &= (~(1 << (i & 7)));
  }

  inline void flip(long i) {
    if (get(i)) reset(i);
    else set(i);
  }

  ~bitvector() {
    free(m_data);
  }

  long m_alloc_bytes;
  unsigned char *m_data;
};

#endif // __BITVECTOR_H_INCLUDED

