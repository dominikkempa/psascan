#ifndef __BITVECTOR_H_INCLUDED
#define __BITVECTOR_H_INCLUDED

#include <cstdio>
#include <cstdlib>

#include "parallel_utils.h"
#include "utils.h"


struct bitvector {
<<<<<<< HEAD:pSAscan/dev/experimental/expalining-differences-between-random-for-random-seq/final-inmem-sufsort/bench/bitvector.h
  bitvector(long length, long max_threads = 1) {
    if (length <= 0) {
      fprintf(stderr, "Error: attempint to construct "
          "empty bitvector.\n");
=======
  bitvector(std::string filename) {
    utils::read_objects_from_file<unsigned char>(m_data, m_alloc_bytes, filename);
  }

  bitvector(long length) : m_alloc_bytes((length + 7) / 8) {
    if (length <= 0) {
      fprintf(stderr, "Error: constructing a bitvector of length 0.\n");
>>>>>>> master:SAscan/dev/src/bitvector.h
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

