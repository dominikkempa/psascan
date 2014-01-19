#ifndef __BITVECTOR_H
#define __BITVECTOR_H

#include <cstdio>
#include <cstdlib>

#include "utils.h"

struct bitvector {
  bitvector(std::string filename) {
    utils::read_file(data, alloc_bytes, filename);
  }

  bitvector(long n) : length(n), alloc_bytes((length + 7) / 8) {
    if (length <= 0) {
      fprintf(stderr, "Error: constructing a bitvector of length 0.\n");
      std::exit(EXIT_FAILURE);
    }
    data = new unsigned char[alloc_bytes];
    std::fill(data, data + alloc_bytes, 0);
  }

  inline bool get(long i) const {
    return data[i >> 3] & (1 << (i & 7));
  }

  inline void set(long i) {
    data[i >> 3] |= (1 << (i & 7));
  }

  inline void reset(long i) {
    data[i >> 3] &= (~(1 << (i & 7)));
  }

  void save(std::string filename) const {
    utils::write_objects_to_file<unsigned char>(data, alloc_bytes, filename);
  }

  void reverse() {
    for (long i = 0, j = length - 1; i < j; ++i, --j) {
      bool left = get(i), right = get(j);
      reset(i); if (right) set(i);
      reset(j); if (left) set(j);
    }
  }

  ~bitvector() {
    delete[] data;
  }

  long length, alloc_bytes;
  unsigned char *data;
};

#endif // __BITVECTOR_H

