/**
 * @file    src/psascan_src/io/async_scatterfile_reader.hpp
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

#ifndef __SRC_PSASCAN_SRC_IO_ASYNC_SCATTERFILE_READER_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_IO_ASYNC_SCATTERFILE_READER_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <string>
#include <vector>
#include <thread>
#include <mutex>
#include <algorithm>
#include <condition_variable>

#include "scatterfile.hpp"
#include "../utils.hpp"


namespace psascan_private {

template<typename value_type>
class async_scatterfile_reader {
  public:
    typedef scatterfile<value_type> scatterfile_type;

  public:
    async_scatterfile_reader(const scatterfile_type *sfile,
        std::uint64_t buf_size_bytes = (1UL << 20)) {
      m_sfile = sfile;
      m_total_items_read_by_user = 0;
      m_active_buf_filled = 0;
      m_passive_buf_filled = 0;
      m_active_buf_pos = 0;

      // Initialize buffers.
      m_buf_size_items = std::max(1UL, (buf_size_bytes / sizeof(value_type)) / 2UL);
      m_active_buf = (value_type *)malloc(m_buf_size_items * sizeof(value_type));
      m_passive_buf = (value_type *)malloc(m_buf_size_items * sizeof(value_type));

      // Start the I/O thread and immediatelly start reading.
      m_avail = true;
      m_thread = new std::thread(async_io_code<value_type>, this);
    }

    inline value_type read() {
      if (m_active_buf_pos == m_active_buf_filled)
        receive_new_buffer();
      ++m_total_items_read_by_user;
      return m_active_buf[m_active_buf_pos++];
    }

    ~async_scatterfile_reader() {

      // Check if all items have been read from scatterfile.
      if (m_total_items_read_by_user != m_sfile->m_items_written) {
        fprintf(stderr, "\nError: not all items were read from scatterfile.\n");
        std::exit(EXIT_FAILURE);
      }

      // Wait for the I/O thread to finish. If the condition above
      // is satisfied, it should have already terminated itself.
      m_thread->join();

      // Clean up.
      delete m_thread;
      free(m_active_buf);
      free(m_passive_buf);
    }

  private:
    template<typename T>
    static void async_io_code(async_scatterfile_reader<T> *caller) {
      std::FILE *file_handler = NULL;
      std::uint64_t cur_file_read = 0;
      std::uint64_t total_items_read = 0;

      while (true) {

        // Wait until the passive buffer is available.
        std::unique_lock<std::mutex> lk(caller->m_mutex);
        while (!caller->m_avail)
          caller->m_cv.wait(lk);
        lk.unlock();

        // Open next file if necessary.
        if (file_handler == NULL) {
          std::uint64_t next_file_id = total_items_read / caller->m_sfile->m_max_items_per_file;
          file_handler = utils::file_open(caller->m_sfile->m_filenames[next_file_id], "r");
          cur_file_read = 0;
        }

        // Read the data from disk.
        std::uint64_t max_file_left = caller->m_sfile->m_max_items_per_file - cur_file_read;
        std::uint64_t sequence_left = caller->m_sfile->m_items_written - total_items_read;
        std::uint64_t this_file_left = std::min(max_file_left, sequence_left);
        caller->m_passive_buf_filled = std::min(caller->m_buf_size_items, this_file_left);
        utils::read_from_file(caller->m_passive_buf, caller->m_passive_buf_filled, file_handler);
        total_items_read += caller->m_passive_buf_filled;
        cur_file_read += caller->m_passive_buf_filled;

        // Delete the current file if we've read all the items.
        if (cur_file_read == caller->m_sfile->m_max_items_per_file ||
            total_items_read == caller->m_sfile->m_items_written) {
          std::fclose(file_handler);
          file_handler = NULL;
          std::uint64_t this_file_id = (total_items_read - 1) / caller->m_sfile->m_max_items_per_file;
          utils::file_delete(caller->m_sfile->m_filenames[this_file_id]);
        }

        // Let the caller know that the I/O thread finished reading.
        lk.lock();
        caller->m_avail = false;
        lk.unlock();
        caller->m_cv.notify_one();

        // Terminate the thread if there is no more data to read.
        if (total_items_read == caller->m_sfile->m_items_written)
          break;
      }
    }

    void receive_new_buffer() {

      // Wait until the I/O thread finishes reading the revious
      // buffer. Most of the time this step is instantaneous.
      std::unique_lock<std::mutex> lk(m_mutex);
      while (m_avail == true)
        m_cv.wait(lk);

      // Set the new active buffer.
      std::swap(m_active_buf, m_passive_buf);
      m_active_buf_filled = m_passive_buf_filled;
      m_active_buf_pos = 0;

      // Let the I/O thead know that it can now
      // prefetch another buffer.
      m_avail = true;
      lk.unlock();
      m_cv.notify_one();
    }

  private:
    const scatterfile_type *m_sfile;
    std::uint64_t m_total_items_read_by_user;

    // Buffers used for asynchronous reading.
    value_type *m_active_buf;
    value_type *m_passive_buf;
    std::uint64_t m_buf_size_items;
    std::uint64_t m_active_buf_pos;
    std::uint64_t m_active_buf_filled;
    std::uint64_t m_passive_buf_filled;

    // For synchronization with thread doing asynchronous reading.
    std::thread *m_thread;
    std::mutex m_mutex;
    std::condition_variable m_cv;
    bool m_finished;
    bool m_avail;
};

}  // psascan_private

#endif // __SRC_PSASCAN_SRC_IO_ASYNC_SCATTERFILE_READER_HPP_INCLUDED
