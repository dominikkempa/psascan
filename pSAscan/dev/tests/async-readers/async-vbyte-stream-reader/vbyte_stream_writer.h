#ifndef __VBYTE_STREAM_WRITER_H_INCLUDED
#define __VBYTE_STREAM_WRITER_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <string>

#include "utils.h"


template<typename value_type>
struct vbyte_stream_writer {
  vbyte_stream_writer(std::string fname, long bufsize = (4L << 20))
      : m_bufsize(bufsize) {
    m_file = utils::open_file(fname, "w");
    m_buf = new unsigned char[m_bufsize + 512];
    m_filled = 0L;
  }

  inline void write(value_type x) {
    if (m_filled > m_bufsize)
      flush();

    while (x > 127) {
      m_buf[m_filled++] = ((x & 0x7f) | 0x80);
      x >>= 7;
    }
    m_buf[m_filled++] = x;
  }

  ~vbyte_stream_writer() {
    if (m_filled)
      flush();

    delete[] m_buf;
    std::fclose(m_file);
  }

  private:
    inline void flush() {
      utils::add_objects_to_file(m_buf, m_filled, m_file);
      m_filled = 0;
    }

  long m_bufsize, m_filled;
  unsigned char *m_buf;
  
  std::FILE *m_file;
};


#endif  // __VBYTE_STREAM_WRITER_H_INCLUDED
