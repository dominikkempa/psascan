/**
 * @file    async_multifile_bit_writer.h
 * @section LICENCE
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

#ifndef __ASYNC_MULTIFILE_BIT_WRITER_H_INCLUDED
#define __ASYNC_MULTIFILE_BIT_WRITER_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <queue>
#include <string>
#include <algorithm>
#include <condition_variable>
#include <mutex>
#include <thread>

#include "../utils.h"


namespace psascan_private {

class async_multifile_bit_writer {
  private:
    struct buffer {
      buffer(std::uint64_t size_bytes) {
        m_size = std::max(1UL, size_bytes);
        m_content = (std::uint8_t *)malloc(m_size + 1);
        m_filled = 0;
        m_bit_pos = 0;
        m_content[m_filled] = 0;
      }

      void flush_to_file(std::FILE *f) {
        if (m_bit_pos != 0) ++m_filled;
        utils::write_to_file(m_content, m_filled, f);
        m_filled = 0;
        m_bit_pos = 0;
        m_content[m_filled] = 0;
      }

      inline void write(std::uint8_t bit) {
        m_content[m_filled] |= (bit << m_bit_pos);
        ++m_bit_pos;
        if (m_bit_pos == 8) {
          m_bit_pos = 0;
          ++m_filled;
          m_content[m_filled] = 0;
        }
      }

      ~buffer() {
        free(m_content);
      }

      inline bool full() const { return m_filled == m_size; }
      inline bool empty() const { return m_filled == 0 && m_bit_pos == 0; }
      inline std::uint64_t size_in_bytes() const { return m_filled; }

      std::uint8_t *m_content;
      std::uint64_t m_bit_pos;
      std::uint64_t m_filled;
      std::uint64_t m_size;
    };

    struct request {
      request(buffer *buffer, std::uint64_t file_id) {
        m_buffer = buffer;
        m_file_id = file_id;
      }

      buffer *m_buffer;
      std::uint64_t m_file_id;
    };

    struct request_queue {
      request_queue()
        : m_no_more_requests(false) {}

      request get() {
        request ret = m_requests.front();
        m_requests.pop();
        return ret;
      }

      inline void add(request req) {
        std::lock_guard<std::mutex> lk(m_mutex);
        m_requests.push(req);
      }

      inline bool empty() const { return m_requests.empty(); }

      std::queue<request> m_requests;  // Must be queue
      std::condition_variable m_cv;
      std::mutex m_mutex;
      bool m_no_more_requests;
    };

    struct buffer_collection {
      // Separate method to (1) hide the implementation of
      // the collection (std::vector) and (2) allow locking.
      inline void add(buffer *buffer) {
        std::lock_guard<std::mutex> lk(m_mutex);
        m_buffers.push_back(buffer);
      }

      buffer* get() {
        buffer *ret = m_buffers.back();
        m_buffers.pop_back();
        return ret;
      }

      inline bool empty() const { return m_buffers.empty(); }

      std::vector<buffer*> m_buffers;
      std::condition_variable m_cv;
      std::mutex m_mutex;
    };

  private:
    static void async_io_thread_code(async_multifile_bit_writer *caller) {
      while (true) {
        // Wait for request or until 'no more requests' flag is set.
        std::unique_lock<std::mutex> lk(caller->m_write_requests.m_mutex);
        while (caller->m_write_requests.empty() &&
            !(caller->m_write_requests.m_no_more_requests))
          caller->m_write_requests.m_cv.wait(lk);

        if (caller->m_write_requests.empty() &&
            caller->m_write_requests.m_no_more_requests) {
          // No more requests -- exit.
          lk.unlock();
          break;
        }

        // Extract the buffer from the collection.
        request req = caller->m_write_requests.get();
        lk.unlock();

        // Process the request.
        caller->m_bytes_written += req.m_buffer->size_in_bytes();
        req.m_buffer->flush_to_file(caller->m_files[req.m_file_id]);

        // Add the (now empty) buffer to the collection
        // of empty buffers and notify the waiting thread.
        caller->m_free_buffers.add(req.m_buffer);
        caller->m_free_buffers.m_cv.notify_one();
      }
    }

  private:
    std::uint64_t m_bytes_written;
    std::uint64_t m_buf_size_bytes;
    std::uint64_t m_items_per_buf;

    std::vector<std::FILE*> m_files;
    std::vector<buffer*> m_buffers;
    buffer_collection m_free_buffers;
    request_queue m_write_requests;
    std::thread *m_io_thread;

    // Issue a request to write to buffer.
    void issue_write_request(std::uint64_t file_id) {
      request req(m_buffers[file_id], file_id);
      m_buffers[file_id] = NULL;
      m_write_requests.add(req);
      m_write_requests.m_cv.notify_one();
    }

    // Get a free buffer from the collection of free buffers.
    buffer* get_free_buffer() {
      std::unique_lock<std::mutex> lk(m_free_buffers.m_mutex);
      while (m_free_buffers.empty())
        m_free_buffers.m_cv.wait(lk);
      buffer *ret = m_free_buffers.get();
      lk.unlock();
      return ret;
    }

  public:
    async_multifile_bit_writer(std::uint64_t buf_size_bytes = (1UL << 20),
        std::uint64_t n_free_buffers = 4UL) {
      // Initialize basic parameters.
      // Works even with n_free_buffers == 0.
      m_bytes_written = 0;
      m_buf_size_bytes = buf_size_bytes;
      m_items_per_buf = std::max(1UL, buf_size_bytes);
      for (std::uint64_t j = 0; j < n_free_buffers; ++j)
        m_free_buffers.add(new buffer(m_buf_size_bytes));
      m_io_thread = new std::thread(async_io_thread_code, this);
    }

    // The added file gets the next available ID, i.e., it becomes
    // m_n_buffers'th file (recall that file IDs start from 0).
    void add_file(std::string filename, std::string write_mode =
        std::string("w")) {
      m_buffers.push_back(new buffer(m_buf_size_bytes));
      m_files.push_back(utils::file_open(filename, write_mode));
    }

    // Write value to i-th file. Files are numbered according
    // to the order in which there we added and numbered from 0.
    void write_to_ith_file(std::uint64_t i, std::uint8_t bit) {
      m_buffers[i]->write(bit);
      if (m_buffers[i]->full()) {
        issue_write_request(i);
        m_buffers[i] = get_free_buffer();
      }
    }

    inline std::uint64_t bytes_written() const {
      return m_bytes_written;
    }

    ~async_multifile_bit_writer() {
      // Flush all buffers.
      std::uint64_t n_buffers = m_buffers.size();
      for (std::uint64_t file_id = 0; file_id < n_buffers; ++file_id) {
        if (!(m_buffers[file_id]->empty()))
          issue_write_request(file_id);
      }

      // Let the I/O thread know that there
      // won't be any more requests.
      std::unique_lock<std::mutex> lk(m_write_requests.m_mutex);
      m_write_requests.m_no_more_requests = true;
      lk.unlock();
      m_write_requests.m_cv.notify_one();

      // Wait for the I/O to finish.
      m_io_thread->join();
      delete m_io_thread;

      // Delete buffers and close files.
      for (std::uint64_t file_id = 0; file_id < n_buffers; ++file_id) {
        delete m_buffers[file_id];  // Can be NULL
        std::fclose(m_files[file_id]);
      }

      // Delete free buffers.
      while (!(m_free_buffers.empty())) {
        buffer *buf = m_free_buffers.get();
        delete buf;
      }
    }
};

}  // namespace psascan_private

#endif  // __ASYNC_MULTIFILE_BIT_WRITER_H_INCLUDED
