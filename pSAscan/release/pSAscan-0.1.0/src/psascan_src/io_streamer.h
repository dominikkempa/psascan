/**
 * @file    src/psascan_src/io_streamer.h
 * @author  Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * @section LICENCE
 *
 * This file is part of pSAscan v0.1.0
 * See: http://www.cs.helsinki.fi/group/pads/
 *
 * Copyright (C) 2014-2015
 *   Juha Karkkainen <juha.karkkainen (at) cs.helsinki.fi>
 *   Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/

#ifndef __PSASCAN_SRC_IO_STREAMER_H_INCLUDED
#define __PSASCAN_SRC_IO_STREAMER_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <string>

#include "utils.h"


namespace psascan_private {

template<typename value_type>
struct stream_reader {
  stream_reader(std::string fname, long buf_bytes = (4L << 20))
      : m_bufelems((buf_bytes + sizeof(value_type) - 1) / sizeof(value_type)) {
    m_file = utils::open_file(fname, "r");
    m_buffer = new value_type[m_bufelems];
    refill();
  }

  ~stream_reader() {
    delete[] m_buffer;
    std::fclose(m_file);
  }

  inline value_type read() {
    if (m_pos == m_filled)
      refill();

    return m_buffer[m_pos++];
  }

  inline bool empty() {
    return (!m_filled && !refill());
  }
  
private:
  long refill() {
    m_filled = std::fread(m_buffer, sizeof(value_type), m_bufelems, m_file);
    m_pos = 0;

    return m_filled;
  }

  long m_bufelems, m_filled, m_pos;
  value_type *m_buffer;

  std::FILE *m_file;
};


template<typename value_type>
struct backward_stream_reader {
  backward_stream_reader(std::string fname, long buf_bytes = (4L << 20))
      : m_bufelems((buf_bytes + sizeof(value_type) - 1) / sizeof(value_type)), m_filled(0L) {
    m_buffer = new value_type[m_bufelems];
    m_file = utils::open_file(fname, "r");    
    std::fseek(m_file, 0L, SEEK_END);
    refill();
  }

  ~backward_stream_reader() {
    delete[] m_buffer;
    std::fclose(m_file);
  }

  inline value_type read() {
    value_type ret = m_buffer[m_pos--];
    if (m_pos < 0L) refill();
    
    return ret;
  }
  
private:
  void refill() {
    long curpos = std::ftell(m_file) / sizeof(value_type);
    long toread = std::min(m_bufelems, curpos - m_filled);

    std::fseek(m_file, -((m_filled + toread) * sizeof(value_type)), SEEK_CUR);
    m_filled = std::fread(m_buffer, sizeof(value_type), toread, m_file);
    m_pos = m_filled - 1L;
  }

  long m_bufelems, m_filled, m_pos;
  value_type *m_buffer;

  std::FILE *m_file;
};


//==============================================================================
// Special version of backward stream reader that allows
// skipping some number of elements at the end of the file.
//==============================================================================
template<typename value_type>
struct backward_skip_stream_reader {
  backward_skip_stream_reader(std::string fname, long skip_elems, long buf_bytes = (4L << 20))
      : m_bufelems((buf_bytes + sizeof(value_type) - 1) / sizeof(value_type)), m_filled(0L) {
    m_buffer = new value_type[m_bufelems];
    m_file = utils::open_file(fname, "r");    
    std::fseek(m_file, -(skip_elems * sizeof(value_type)), SEEK_END);
    refill();
  }

  ~backward_skip_stream_reader() {
    delete[] m_buffer;
    std::fclose(m_file);
  }

  inline value_type read() {
    value_type ret = m_buffer[m_pos--];
    if (m_pos < 0L) refill();

    return ret;
  }

private:
  void refill() {
    long curpos = std::ftell(m_file) / sizeof(value_type);
    long toread = std::min(m_bufelems, curpos - m_filled);

    std::fseek(m_file, -((m_filled + toread) * sizeof(value_type)), SEEK_CUR);
    m_filled = std::fread(m_buffer, sizeof(value_type), toread, m_file);
    m_pos = m_filled - 1L;
  }

  long m_bufelems, m_filled, m_pos;
  value_type *m_buffer;

  std::FILE *m_file;
};

template<typename value_type>
struct stream_writer {
  stream_writer(std::string fname, long bufsize = (4L << 20))
      : m_bufelems((bufsize + sizeof(value_type) - 1) / sizeof(value_type)) {
    m_file = utils::open_file(fname.c_str(), "w");
    m_buffer = new value_type[m_bufelems];
    m_filled = 0;
  }

  inline void flush() {
    utils::add_objects_to_file(m_buffer, m_filled, m_file);
    m_filled = 0;
  }

  void write(value_type x) {
    m_buffer[m_filled++] = x;

    if (m_filled == m_bufelems)
      flush();
  }

  ~stream_writer() {
    if (m_filled)
      flush();

    delete[] m_buffer;
    std::fclose(m_file);
  }

private:
  long m_bufelems, m_filled;
  value_type *m_buffer;

  std::FILE *m_file;
};


struct bit_stream_reader {
  bit_stream_reader(std::string filename) {
    m_file = utils::open_file(filename.c_str(), "r");
    m_buf = new unsigned char[k_bufsize];
    refill();
  }

  inline bool read() {
    bool ret = m_buf[m_pos_byte] & (1 << m_pos_bit);
    ++m_pos_bit;
    if (m_pos_bit == 8) {
      m_pos_bit = 0;
      ++m_pos_byte;
      if (m_pos_byte == m_filled)
        refill();
    }

    return ret;
  }

  ~bit_stream_reader() {
    delete[] m_buf;
    std::fclose(m_file);
  }

private:
  inline void refill() {
    m_filled = std::fread(m_buf, 1, k_bufsize, m_file);
    m_pos_byte = m_pos_bit = 0;
  }

  static const long k_bufsize = (2L << 20);

  std::FILE *m_file;

  unsigned char *m_buf;
  long m_filled, m_pos_byte;
  int m_pos_bit;
};


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
    if (pos_bit) ++filled;  // in case this is the final flush
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
  static const long bufsize = (1L << 20);
  
  unsigned char *buf;
  long filled;
  int pos_bit;

  std::FILE *f;
};


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


struct vbyte_stream_reader {
  vbyte_stream_reader(std::string fname, long bufsize = (4L << 20))
      : m_bufsize(bufsize) {
    m_file = utils::open_file(fname, "r");
    m_buf = new unsigned char[m_bufsize];
    refill();
  }

  inline long read() {
    long ret = 0, offset = 0;
    while (m_buf[m_pos] & 0x80) {
      ret |= (((long)m_buf[m_pos++] & 0x7f) << offset);
      if (m_pos == m_filled)
        refill();
      offset += 7;
    }
    ret |= ((long)m_buf[m_pos++] << offset);
    if (m_pos == m_filled)
      refill();

    return ret;
  }
  
  ~vbyte_stream_reader() {
    delete[] m_buf;
    std::fclose(m_file);
  }
  
private:
  inline void refill() {
    m_filled = std::fread(m_buf, 1, m_bufsize, m_file);
    m_pos = 0;
  }

  long m_bufsize, m_filled, m_pos;
  unsigned char *m_buf;
  
  std::FILE *m_file;
};

}  // namespace psascan_private

#endif  // __PSASCAN_SRC_IO_STREAMER_H_INCLUDED
