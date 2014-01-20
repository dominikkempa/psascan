#ifndef __GAP_ARRAY_H
#define __GAP_ARRAY_H

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>

#include "utils.h"

struct buffered_gap_array {
  buffered_gap_array(long n) {
    if (n <= 0) {
      fprintf(stderr, "Error: constructing a gap array of length 0.\n");
      std::exit(EXIT_FAILURE);
    }
    length = n;
    count = new unsigned char[length];
    std::fill(count, count + length, 0);
    buf = new long[bufsize];
    filled = 0;
  }

  inline void flush() {
    for (int i = 0; i < filled; ++i) {
      long pos = buf[i];
      ++count[pos];
      if (!count[pos])
        excess.push_back(pos);
    }
    filled = 0;
  }

  inline void increment(long i) {
    buf[filled++] = i;
    if (filled == bufsize) flush();
  }

  ~buffered_gap_array() {
    delete[] count;
    delete[] buf;
  }
  
  // Store to file using v-byte encoding.  
  void save_to_file(std::string fname) {
    flush();
    
    fprintf(stderr, "  gap->excess.size() = %lu\n", excess.size());
    std::sort(excess.begin(), excess.end());
    FILE *f = utils::open_file(fname.c_str(), "w");
    unsigned char *buffer = (unsigned char *)buf;
    filled = 0;
    
    long max_pos = excess.size();
    for (long j = 0, pos = 0; j < length; ++j) {
      long c = 0;
      while (pos < max_pos && excess[pos] == j) ++pos, ++c;
      long gap_j = count[j] + (c << 8);
      while (gap_j > 127) {
        buffer[filled++] = ((gap_j & 0x7f) | 0x80);
        gap_j >>= 7;
        if (filled == bufsize) {
          utils::add_objects_to_file<unsigned char>(buffer, filled, f);
          filled = 0;
        }
      }
      buffer[filled++] = gap_j;
      if (filled == bufsize) {
        utils::add_objects_to_file<unsigned char>(buffer, filled, f);
        filled = 0;
      }
    }
    if (filled)
      utils::add_objects_to_file<unsigned char>(buffer, filled, f);

    fclose(f);
  }

  unsigned char *count;
  std::vector<long> excess;

  static const int bufsize = (1 << 19); // 4MB
  int filled;
  long length, *buf;
};

#endif // __GAP_ARRAY_H

