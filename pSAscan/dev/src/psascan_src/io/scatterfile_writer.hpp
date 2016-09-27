/**
 * @file    psascan_src/io/scatterfile_writer.hpp
 * @section LICENCE
 *
 * This file is part of pSAscan v0.2.0
 * See: http://www.cs.helsinki.fi/group/pads/
 *
 * Copyright (C) 2014-2016
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

#ifndef __PSASCAN_SRC_IO_SCATTERFILE_WRITER_HPP_INCLUDED
#define __PSASCAN_SRC_IO_SCATTERFILE_WRITER_HPP_INCLUDED

#include <cstdio>
#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>

#include "../utils.hpp"
#include "scatterfile.hpp"


namespace psascan_private {

template<typename value_type>
class scatterfile_writer {
  public:
    typedef scatterfile<value_type> scatterfile_type;

  public:
    scatterfile_writer(scatterfile_type *sfile, std::string filename,
        std::uint64_t buf_size_bytes = (1UL << 20)) {
      m_file = NULL;
      m_sfile = sfile;
      m_filename = filename;
      m_sfile->m_items_written = 0;

      // Initialize buffer.
      m_filled = 0;
      m_buf_size_items = std::max(1UL, buf_size_bytes / sizeof(value_type));
      m_buffer = (value_type *)malloc(m_buf_size_items * sizeof(value_type));
    }

    void write(value_type value) {
      // Open new file if necessary.
      if (m_file == NULL)
        make_new_file();

      // Write value to buffer.
      m_buffer[m_filled++] = value;
      m_cur_file_written++;
      m_sfile->m_items_written++;

      // Flush buffer if necessary.
      if (is_buffer_full())
        flush();
    }

    void write(value_type *tab, std::uint64_t length) {
      // Open new file if necessary.
      if (m_file == NULL)
        make_new_file();

      // Compute how many items can still be written to buffers
      // without exceeding its size of max items per file.
      std::uint64_t buf_left = std::min(m_buf_size_items - m_filled,
          m_sfile->m_max_items_per_file - m_cur_file_written);
      if (length <= buf_left) {
        // Case I: there is enough room, add items to the buffer.
        std::copy(tab, tab + length, m_buffer + m_filled);
        m_filled += length;
        m_cur_file_written += length;
        m_sfile->m_items_written += length;

        // Flush buffer if necessary.
        if (is_buffer_full())
          flush();
      } else {
        // Case II: not enough room in the buffer. Write tab
        // directly to file(s), bypassing the buffer.
        if (m_filled > 0)
          flush();

        while (length > 0) {
          if (m_file == NULL)
            make_new_file();

          // Write as many items as possible to the current file.
          std::uint64_t towrite = std::min(length,
              m_sfile->m_max_items_per_file - m_cur_file_written);
          utils::write_to_file(tab, towrite, m_file);
          m_cur_file_written += towrite;
          m_sfile->m_items_written += towrite;
          tab += towrite;
          length -= towrite;

          // Close the file if it's full.
          if (m_cur_file_written == m_sfile->m_max_items_per_file) {
            std::fclose(m_file);
            m_file = NULL;
          }
        }
      }
    }

    inline void write(value_type *begin, value_type *end) {
      if (begin <= end)
        write(begin, (std::uint64_t)(end - begin));
    }

    ~scatterfile_writer() {
      if (m_filled > 0)
        utils::write_to_file(m_buffer, m_filled, m_file);
      if (m_file != NULL)
        std::fclose(m_file);
      free(m_buffer);
    }

  private:
    inline bool is_buffer_full() const {
      return m_filled == m_buf_size_items ||
        m_cur_file_written == m_sfile->m_max_items_per_file;
    }

    inline void flush() {
      utils::write_to_file(m_buffer, m_filled, m_file);
      m_filled = 0;
      if (m_cur_file_written == m_sfile->m_max_items_per_file) {
        std::fclose(m_file);
        m_file = NULL;
      }
    }

    void make_new_file() {
      m_cur_file_written = 0;
      std::string next_filename = m_filename + ".sfile." + utils::random_string_hash();
      m_sfile->m_filenames.push_back(next_filename);
      m_file = utils::file_open(next_filename, "w");
    }

  private:
    std::uint64_t m_filled;
    std::uint64_t m_buf_size_items;
    std::uint64_t m_cur_file_written;
    value_type *m_buffer;

    std::FILE *m_file;
    scatterfile_type *m_sfile;
    std::string m_filename;
};

}  // psascan_private

#endif // __PSASCAN_SRC_IO_SCATTERFILE_WRITER_HPP_INCLUDED
