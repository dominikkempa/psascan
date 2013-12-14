// Various types of streamers.
#ifndef __STREAM
#define __STREAM

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <string>

#include "utils.h"

/********************************* usage ***************************************
stream_reader<int> *sr = new stream_reader<int>("input.txt", 1 << 21);
while (!sr->empty()) {
  int next = sr->read();
}
delete sr;
*******************************************************************************/

template<typename data_type>
struct stream_reader {
  stream_reader(std::string fname, int bufsize)
      : m_bufelems(bufsize / sizeof(data_type)) {
    f = utils::open_file(fname, "r");
    buffer = new data_type[m_bufelems];
    if (!buffer) {
      fprintf(stderr, "Error: not enough memory for stream reader\n");
      std::exit(EXIT_FAILURE);
    }
    refill();
  }

  ~stream_reader() {
    delete[] buffer;
    std::fclose(f);
  }

  inline data_type read() {
    data_type ret = buffer[pos++];
    if (pos == filled) refill();
    
    return ret;
  }

  inline bool empty() {
    return (!filled && !refill());
  }
  
private:
  int m_bufelems, filled, pos;
  data_type *buffer;
  std::FILE *f;

  int refill() {
    filled = std::fread(buffer, sizeof(data_type), m_bufelems, f);
    pos = 0;
    return filled;
  }
};

/********************************* usage ***************************************
stream_writer<int> *sw = new stream_writer<int>("output.txt", 1 << 21);
for (int i = 0; i < n; ++i)
  sw->write(SA[i]);
delete sw;
*******************************************************************************/

template<typename data_type>
struct stream_writer {
  stream_writer(std::string fname, int bufsize)
      : m_bufelems(bufsize / sizeof(data_type)) {
    f = utils::open_file(fname.c_str(), "w");
    buffer = new data_type[m_bufelems];
    if (!buffer) {
      fprintf(stderr, "Error: allocation error in stream_writer\n");
      std::exit(EXIT_FAILURE);
    }
    filled = 0;
  }

  ~stream_writer() {
    if (filled) flush();
    delete[] buffer;
    fclose(f);
  }

  void write(data_type x) {
    buffer[filled++] = x;
    if (filled == m_bufelems) flush();
  }

private:
  void flush() {
    utils::add_objects_to_file(buffer, filled, f);
    filled = 0;
  }

  std::FILE *f;

  int m_bufelems, filled;
  data_type *buffer;
};

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

struct vbyte_stream_reader {
  vbyte_stream_reader(std::string fname, int bufsize)
      : m_bufsize(bufsize) {
    f = utils::open_file(fname, "r");
    buf = new unsigned char[m_bufsize];
    if (!buf) {
      fprintf(stderr, "Error: allocation failed in vbyte_stream_reader\n");
      std::exit(EXIT_FAILURE);
    }
    refill();
  }

  ~vbyte_stream_reader() {
    delete[] buf;
    std::fclose(f);
  }

  inline int read() {
    int ret = 0, offset = 0;
    while (buf[pos] & 0x80) {
      ret |= ((buf[pos++] & 0x7f) << offset);
      if (pos == filled) refill();
      offset += 7;
    }
    ret |= (buf[pos++] << offset);
    if (pos == filled) refill();

    return ret;
  }
  
private:
  int m_bufsize, filled, pos;
  unsigned char *buf;
  std::FILE *f;

  inline void refill() {
    filled = fread(buf, 1, m_bufsize, f);
    pos = 0;
  }
};

struct text_reader {
  text_reader(std::string filename) {
    f = utils::open_file(filename, "r");
    buffersize = (1 << 21);
    buf = new unsigned char[buffersize];
  }

  ~text_reader() {
    fclose(f);
    delete[] buf;
  }

  void read_block(int beg, int length, unsigned char *b) {
    fseek(f, beg, SEEK_SET);
    utils::read_objects_from_file<unsigned char>(b, length, f);
  }

  void refill() {
    int curpos = ftell(f);
    int left = std::min(curpos - filled, buffersize);
    fseek(f, -(filled + left), SEEK_CUR);
    filled = fread(buf, 1, left, f);
    pos = filled - 1;
  }

  void init_backward_streaming() {
    fseek(f, 0, SEEK_END);
    filled = 0;
    refill();
  }

  inline unsigned char read_next() {
    unsigned char ret = buf[pos--];
    if (pos < 0) refill();

    return ret;
  }

  std::FILE *f;

//  static const int buffersize = (1 << 21); // ???
  int buffersize;
  unsigned char *buf;
  int filled, pos;
};


#endif // __STREAM
