#ifndef __BIT_STREAM_READER
#define __BIT_STREAM_READER

#include <cstdio>
#include <string>

#include "utils.h"

struct bit_stream_reader {
  bit_stream_reader(std::string filename) {
    f = utils::open_file(filename.c_str(), "r");
    buf = new unsigned char[bufsize];
    refill();
  }

  inline void refill() {
    filled = fread(buf, 1, bufsize, f);
    pos_byte = pos_bit = 0;
  }

  inline bool read() {
    bool ret = buf[pos_byte] & (1 << pos_bit);
    ++pos_bit;
    if (pos_bit == 8) {
      pos_bit = 0;
      ++pos_byte;
      if (pos_byte == filled)
        refill();
    }

    return ret;
  }

  ~bit_stream_reader() {
    fclose(f);
    delete[] buf;
  }

  std::FILE *f;
  static const int bufsize = (2 << 20); // 1MB
  unsigned char *buf;
  int filled, pos_byte, pos_bit;
};

#endif // __BIT_STREAM_READER

