// Simple bitvector class.
#ifndef __BITVECTOR_H_INCLUDED
#define __BITVECTOR_H_INCLUDED

#include <cstdio>
#include <cstdlib>

#include "parallel_utils.h"
#include "utils.h"

struct bitvector {
  bitvector(long length, long max_threads) {
    if (length <= 0) {
      fprintf(stderr, "Error: attempint to construct "
          "empty bitvector.\n");
      std::exit(EXIT_FAILURE);
    }

    m_alloc_bytes = (length + 7) / 8;
    m_data = new unsigned char[m_alloc_bytes];
    parallel_utils::fill(m_data, m_alloc_bytes, (unsigned char)0,
        max_threads);
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

  ~bitvector() {
    delete[] m_data;
  }

  long m_alloc_bytes;
  unsigned char *m_data;
};

#endif // __BITVECTOR_H_INCLUDED

