// Various types of streamers.
#ifndef __IO_STREAMER_H_INCLUDED
#define __IO_STREAMER_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <string>

#include "utils.h"


struct bit_stream_writer {
  bit_stream_writer(std::string filename) {
    f = utils::open_file(filename, "w");
    buf = new unsigned char[bufsize];
    if (!buf) {
      fprintf(stderr, "\nError: allocation error in bit_stream_writer\n");
      std::exit(EXIT_FAILURE);
    }
    std::fill(buf, buf + bufsize, 0);
    filled = pos_bit = 0;
  }

  inline void flush() {
    if (pos_bit) ++filled; // final flush?
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
  
  ~bit_stream_writer() {
    flush();
    std::fclose(f);
    delete[] buf;
  }

private:
  static const int bufsize = (1 << 20); // 1MB
  
  unsigned char *buf;
  int filled, pos_bit;

  std::FILE *f;
};

#endif  // __IO_STREAMER_H_INCLUDED
