/**
 * @file    src/psascan_src/distributed_file.h
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

#ifndef __PSASCAN_SRC_DISTRIBUTED_FILE_H_INCLUDED
#define __PSASCAN_SRC_DISTRIBUTED_FILE_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <string>
#include <thread>
#include <mutex>
#include <algorithm>
#include <condition_variable>

#include "utils.h"


namespace psascan_private {

template<typename value_type>
class distributed_file {
  public:
    distributed_file(std::string filename_base, std::uint64_t max_bytes) {
      m_max_items = std::max(1UL, max_bytes / sizeof(value_type));
      m_filename = filename_base + ".distrfile." + utils::random_string_hash();
    }

    distributed_file(std::string filename_base, std::uint64_t max_bytes,
        const value_type *begin, const value_type *end) {
      m_max_items = std::max(1UL, max_bytes / sizeof(value_type));
      m_filename = filename_base + ".distrfile." + utils::random_string_hash();

      initialize_writing();
      write(begin, end);
      finish_writing();
    }

    void initialize_writing() {
      m_total_write = 0;
      m_files_cnt = 0;
      make_new_file();
    }

    void write(const value_type *begin, const value_type *end) {
      // Fill the current file.
      std::uint64_t seq_left = end - begin;
      if (m_cur_file_write != m_max_items) {
        std::uint64_t file_left = m_max_items - m_cur_file_write;
        std::uint64_t towrite = std::min(file_left, seq_left);
        utils::write_to_file(begin, towrite, m_file);
        m_cur_file_write += towrite;
        m_total_write += towrite;
        begin += towrite;
        seq_left -= towrite;
      }

      // Write remaining items.
      while (seq_left > 0) {
        std::fclose(m_file);
        make_new_file();
        std::uint64_t towrite = std::min(m_max_items, seq_left);
        utils::write_to_file(begin, towrite, m_file);
        m_cur_file_write += towrite;
        m_total_write += towrite;
        begin += towrite;
        seq_left -= towrite;
      }
    }

    void finish_writing() {
      if (m_cur_file_write == 0) {
        fprintf(stderr, "\nError: nothing was ever written to %s\n", m_filename.c_str());
        std::exit(EXIT_FAILURE);
      }
      std::fclose(m_file);
    }

    void initialize_reading(std::uint64_t bufsize = (4UL << 20)) {
      // Initialize counters.
      m_active_buf_filled = 0;
      m_passive_buf_filled = 0;
      m_active_buf_pos = 0;
      m_total_read_buf = 0;
      m_total_read_user = 0;
      m_cur_file = -1;

      // Initialize buffers.
      m_buf_size = std::max(1UL, bufsize / (2UL * sizeof(value_type)));
      m_active_buf = (value_type *)malloc(m_buf_size * sizeof(value_type));
      m_passive_buf = (value_type *)malloc(m_buf_size * sizeof(value_type));

      // Start the I/O thread and immediatelly start reading.
      m_avail = true;
      m_finished = false;
      m_thread = new std::thread(async_io_code<value_type>, this);
    }

    inline value_type read() {
      if (m_active_buf_pos == m_active_buf_filled)
        receive_new_buffer();
      m_total_read_user++;
      return m_active_buf[m_active_buf_pos++];
    }

    void finish_reading() {
      if (m_total_read_buf != m_total_read_user || m_total_read_user != m_total_write) {
        fprintf(stderr, "\nError: not all elems were read from distributed file %s\n", m_filename.c_str());
        std::exit(EXIT_FAILURE);
      }

      // Let the I/O thread know that we are done.
      std::unique_lock<std::mutex> lk(m_mutex);
      m_finished = true;
      lk.unlock();
      m_cv.notify_one();

      // Wait for the thread to finish.
      m_thread->join();

      // Clean up.
      delete m_thread;
      close_and_destroy_cur_file();
      free(m_active_buf);
      free(m_passive_buf);
    }

  private:
    void close_and_destroy_cur_file() {
      if (!m_file) {
        fprintf(stderr, "\nError: deleting a NULL file\n");
        std::exit(EXIT_FAILURE);
      }

      std::fclose(m_file);
      std::string cur_fname = m_filename + ".part" + utils::intToStr(m_cur_file);
      utils::file_delete(cur_fname);
    }

    template<typename T>
    static void async_io_code(distributed_file<T> *file) {
      while (true) {
        // Wait until the passive buffer is available.
        std::unique_lock<std::mutex> lk(file->m_mutex);
        while (!(file->m_avail) && !(file->m_finished))
          file->m_cv.wait(lk);

        if (!(file->m_avail) && (file->m_finished)) {
          // We're done, terminate the thread.
          lk.unlock();
          return;
        }
        lk.unlock();

        // This should never happen.
        if (file->m_total_read_buf == file->m_total_write) {
          fprintf(stderr, "\nError: trying to read past the end of file\n");
          std::exit(EXIT_FAILURE);
        }

        // Safely process the passive buffer.
        // Check if we need to open next file.
        if (file->m_cur_file == -1 || file->m_cur_file_read == file->m_max_items) {
          if (file->m_cur_file != -1)
            file->close_and_destroy_cur_file();
          file->open_next_file();
        }

        // Read the data from disk.
        std::uint64_t file_left = file->m_max_items - file->m_cur_file_read;
        std::uint64_t items_left = file->m_total_write - file->m_total_read_buf;
        std::uint64_t left = std::min(file_left, items_left);
        file->m_passive_buf_filled = std::min(left, file->m_buf_size);
        file->m_cur_file_read += file->m_passive_buf_filled;
        file->m_total_read_buf += file->m_passive_buf_filled;
        utils::read_from_file(file->m_passive_buf, file->m_passive_buf_filled, file->m_file);

        // Let the caller know that the I/O thread finished reading.
        lk.lock();
        file->m_avail = false;
        lk.unlock();
        file->m_cv.notify_one();
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
      m_avail = (m_total_read_buf < m_total_write);
      lk.unlock();
      m_cv.notify_one();
    }

    void open_next_file() {
      ++m_cur_file;
      m_file = utils::file_open(m_filename + ".part" + utils::intToStr(m_cur_file), "r");
      m_cur_file_read = 0;
    }

    void make_new_file() {
      m_file = utils::file_open(m_filename + ".part" + utils::intToStr(m_files_cnt), "w");
      ++m_files_cnt;
      m_cur_file_write = 0;
    }

  private:
    std::FILE *m_file;         // file handler
    std::string m_filename;    // file name base
    std::uint64_t m_max_items; // max items per file

    // Buffers used for asynchronous reading.
    value_type *m_active_buf;
    value_type *m_passive_buf;
    std::uint64_t m_buf_size;
    std::uint64_t m_active_buf_pos;
    std::uint64_t m_active_buf_filled;
    std::uint64_t m_passive_buf_filled;

    // Various housekeeping statistics about the number of items.
    std::uint64_t m_cur_file_write;   // number of items written to a current file
    std::uint64_t m_total_write;      // total number of written items
    std::uint64_t m_cur_file_read;    // number of items read from the current file
    std::uint64_t m_total_read_buf;   // total number of items read from files into buffers
    std::uint64_t m_total_read_user;  // total number of items read by the user

    // Used to keep track of file count.
    std::uint64_t m_files_cnt; // counts the files during writing
    std::int64_t m_cur_file;  // iterates through [0..m_files_cnt) during reading

    // For synchronization with thread doing asynchronous reading.
    std::thread *m_thread;
    std::mutex m_mutex;
    std::condition_variable m_cv;
    bool m_finished;
    bool m_avail;
};

}  // psascan_private

#endif // __PSASCAN_SRC_DISTRIBUTED_FILE_H_INCLUDED
