#ifndef __BUFFERED_GAP_ARRAY
#define __BUFFERED_GAP_ARRAY

#include <cstdio>
#include <vector>
#include <algorithm>

struct buffered_gap_array {
  buffered_gap_array(int length) {
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

  unsigned char *count;
  std::vector<int> excess;

  static const int bufsize = (1 << 19); // 2MB
  int *buf, filled;
};

#endif // __BUFFERED_GAP_ARRAY

