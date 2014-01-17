#ifndef __BITVECTOR_H
#define __BITVECTOR_H

#include "utils.h"

struct bitvector {
  bitvector(std::string filename) {
    utils::read_file(data, length, filename);
  }

  inline bool get(long i) const {
    return data[i >> 3] & (1 << (i & 7));
  }

  inline void set(long i) {
    data[i >> 3] |= (1 << (i & 7));
  }

  void save(std::string filename) {
    utils::write_objects_to_file<unsigned char>(data, length, filename);
  }

  bitvector(long n) {
    length = (n + 7) / 8;
    data = new unsigned char[length];
    std::fill(data, data + length, 0);
  }

  ~bitvector() {
    if (data)
      delete[] data;
  }

  long length;
  unsigned char *data;
};

#endif // __BITVECTOR_H

