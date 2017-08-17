/**
 * @file    src/psascan_src/io/async_stream_vbyte_reader.hpp
 * @section LICENCE
 *
 * This file is part of pSAscan v0.2.0
 * See: http://www.cs.helsinki.fi/group/pads/
 *
 * Copyright (C) 2014-2017
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

#ifndef __SRC_PSASCAN_SRC_IO_ASYNC_STREAM_VBYTE_READER_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_IO_ASYNC_STREAM_VBYTE_READER_HPP_INCLUDED

#include <cstdint>
#include <string>

#include "async_stream_reader.hpp"


namespace psascan_private {

class async_stream_vbyte_reader {
  private:
    typedef async_stream_reader<std::uint8_t>
      internal_reader_type;
    internal_reader_type *m_internal_reader;

  public:
    async_stream_vbyte_reader(
        std::string filename) {
      init(filename, (8UL << 20), 4);
    }

    async_stream_vbyte_reader(
        std::string filename,
        std::uint64_t total_buf_size_bytes,
        std::uint64_t n_buffers) {
      init(filename, total_buf_size_bytes, n_buffers);
    }

    void init(std::string filename,
        std::uint64_t total_buf_size_bytes,
        std::uint64_t n_buffers) {
      m_internal_reader = new internal_reader_type(
          filename, total_buf_size_bytes, n_buffers);
    }

    inline std::uint64_t read() {
      std::uint64_t result = 0;
      std::uint64_t offset = 0;
      std::uint64_t next_char = m_internal_reader->read();
      while (next_char & 0x80) {
        result |= ((next_char & 0x7f) << offset);
        offset += 7;
        next_char = m_internal_reader->read();
      }
      result |= (next_char << offset);
      return result;
    }

    inline std::uint64_t bytes_read() {
      return m_internal_reader->bytes_read();
    }

    inline void stop_reading() {
      m_internal_reader->stop_reading();
    }

    ~async_stream_vbyte_reader() {
      delete m_internal_reader;
    }
};

}  // namespace psascan_private

#endif  // __SRC_PSASCAN_SRC_IO_ASYNC_STREAM_VBYTE_READER_HPP_INCLUDED
