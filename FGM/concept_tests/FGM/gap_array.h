#ifndef __GAP_ARRAY_H
#define __GAP_ARRAY_H

#include <cstdio>
#include <vector>
#include <algorithm>

#include "utils.h"

struct buffered_gap_array {
  buffered_gap_array(int n) {
    length = n;
    count = new unsigned char[length];
    std::fill(count, count + length, 0);
    buf = new int[bufsize];
    filled = 0;
  }

  inline void flush() {
    std::sort(buf, buf + filled);
    for (int i = 0; i < filled; ++i) {
      int pos = buf[i];
      int old_count = count[pos];
      ++count[pos];
      if (old_count == 255)
        excess.push_back(pos);
    }
    filled = 0;
  }

  inline void increment(int i) {
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
    std::sort(excess.begin(), excess.end());
    FILE *f = utils::open_file(fname.c_str(), "w");
    unsigned char *buffer = (unsigned char *)buf;
    filled = 0;
    
    int max_pos = (int)excess.size();
    for (int j = 0, pos = 0; j < length; ++j) {
      int c = 0;
      while (pos < max_pos && excess[pos] == j) ++pos, ++c;
      int gap_j = count[j] + (c << 8);
      while (gap_j > 127) {
        buffer[filled++] = ((gap_j & 0x7f) | 0x80);
        gap_j >>= 7;
        if (filled == bufsize)
          utils::add_objects_to_file<unsigned char>(buffer, filled, f);
      }
      buffer[filled++] = gap_j;
      if (filled == bufsize)
        utils::add_objects_to_file<unsigned char>(buffer, filled, f);
    }
    if (filled)
      utils::add_objects_to_file<unsigned char>(buffer, filled, f);

    fclose(f);
  }

  unsigned char *count;
  std::vector<int> excess;

  static const int bufsize = (1 << 19); // 2MB
  int *buf, filled, length;
};

#endif // __GAP_ARRAY_H

