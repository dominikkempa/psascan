#ifndef __BIT_STREAM_WRITER
#define __BIT_STREAM_WRITER

#include <cstdio>
#include <string>

#include "utils.h"

struct bit_stream_writer {
  bit_stream_writer(std::string filename) {
    f = utils::open_file(filename, "w");
    buf = new unsigned char[bufsize];
    std::fill(buf, buf + bufsize, 0);
    filled = pos_bit = 0;
  }

  ~bit_stream_writer() {
    flush();
    fclose(f);
    delete[] buf;
  }

  inline void flush() {
    if (pos_bit) ++filled; // in case this is the final flush.
    utils::add_objects_to_file<unsigned char>(buf, filled, f);
    filled = pos_bit = 0;
    std::fill(buf, buf + bufsize, 0);
  }

  void write(int bit) {
    buf[filled] |= (bit << pos_bit);
    ++pos_bit;
    if (pos_bit == 8) {
      pos_bit = 0;
      ++filled;
      if (filled == bufsize)
        flush();
    }
  }

  std::FILE *f;
  static const int bufsize = (2 << 20); // 2MB
  unsigned char *buf;
  int filled, pos_bit;
};

#endif // __BIT_STREAM_WRITER

