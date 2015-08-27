/**
 * @file    async_stream_reader.h
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

#ifndef __ASYNC_STREAM_READER_H_INCLUDED
#define __ASYNC_STREAM_READER_H_INCLUDED

#include <cstdio>
#include <thread>
#include <future>
#include <mutex>
#include <condition_variable>
#include <algorithm>

#include "utils.h"


template<typename value_type>
class async_stream_reader {
  private:
    template<typename T>
    static void io_thread_code(async_stream_reader<T> *reader) {
      while (true) {
        // Wait until the passive buffer is available.
        std::unique_lock<std::mutex> lk(reader->m_mutex);
        while (!(reader->m_avail) && !(reader->m_finished))
          reader->m_cv.wait(lk);

        if (!(reader->m_avail) && (reader->m_finished)) {
          // We're done, terminate the thread.
          lk.unlock();
          return;
        }
        lk.unlock();

        // Safely read the data to disk.
        reader->m_passive_buf_filled = std::fread(reader->m_passive_buf,
            sizeof(T), reader->m_items_per_buf, reader->m_file);
        reader->m_bytes_read += reader->m_passive_buf_filled * sizeof(T);

        // Let the caller know what the I/O thread finished reading.
        lk.lock();
        reader->m_avail = false;
        lk.unlock();
        reader->m_cv.notify_one();
      }
    }

    // Check if the reading thread has already prefetched the next
    // buffer (the request should have been done before), and wait
    // if the prefetching was not completed yet.
    void receive_new_buffer() {
      // Wait until the I/O thread finishes reading the previous
      // buffer. Most of the time this step is instantaneous.
      std::unique_lock<std::mutex> lk(m_mutex);
      while (m_avail == true)
        m_cv.wait(lk);

      // Set the new active buffer.
      std::swap(m_active_buf, m_passive_buf);
      m_active_buf_filled = m_passive_buf_filled;
      m_active_buf_pos = 0;

      // Let the I/O thread know that it can now prefetch
      // another buffer.
      m_avail = true;
      lk.unlock();
      m_cv.notify_one();
    }

  public:
    async_stream_reader(std::string filename = std::string(""),
        std::size_t buf_size_bytes = (2UL << 20)) {
      if (filename.empty()) m_file = stdin;
      else m_file = utils::open_file(filename.c_str(), "r");

      // Initialize buffers.
      m_bytes_read = 0;
      m_active_buf_filled = 0;
      m_passive_buf_filled = 0;
      m_active_buf_pos = 0;
      m_items_per_buf = std::max(1UL, buf_size_bytes / (2 * sizeof(value_type)));
      m_active_buf = (value_type *)malloc(m_items_per_buf * sizeof(value_type));
      m_passive_buf = (value_type *)malloc(m_items_per_buf * sizeof(value_type));

      // Start the I/O thread and immediatelly start reading.
      m_finished = false;
      m_avail = true;
      m_thread = new std::thread(io_thread_code<value_type>, this);
    }
  
    ~async_stream_reader() {
      // Let the I/O thread know that we're done.
      std::unique_lock<std::mutex> lk(m_mutex);
      m_finished = true;
      lk.unlock();
      m_cv.notify_one();

      // Wait for the thread to actually finish.
      m_thread->join();
    
      // Clean up.
      delete m_thread;
      free(m_active_buf);
      free(m_passive_buf);
      if (m_file != stdin)
        std::fclose(m_file);
    }

    inline value_type read() {
      if (m_active_buf_pos == m_active_buf_filled) {
        // The active buffer run out of data.
        // At this point we need to swap it with the passive
        // buffer. The request to read that passive buffer should
        // have been scheduled long time ago, so hopefully the
        // buffer is now available. We check for that, but we
        // also might wait a little, if the reading has not yet
        // been finished. At this point we also already schedule
        // the next read.
        receive_new_buffer();
      }

      return m_active_buf[m_active_buf_pos++];
    }

    inline value_type peek() {
      if (m_active_buf_pos == m_active_buf_filled)
        receive_new_buffer();

      return m_active_buf[m_active_buf_pos];
    }

    inline bool empty() {
      if (m_active_buf_pos == m_active_buf_filled)
        receive_new_buffer();
      return (m_active_buf_pos == m_active_buf_filled);
    }

    inline std::size_t bytes_read() const {
      return m_bytes_read;
    }

  private:
    value_type *m_active_buf;
    value_type *m_passive_buf;

    std::size_t m_items_per_buf;
    std::size_t m_active_buf_pos;
    std::size_t m_active_buf_filled;
    std::size_t m_passive_buf_filled;
    std::size_t m_bytes_read;

    // Used for synchronization with the I/O thread.
    std::mutex m_mutex;
    std::condition_variable m_cv;
    bool m_avail;
    bool m_finished;

    std::FILE *m_file;
    std::thread *m_thread;
};

#endif  // __ASYNC_STREAM_READER_H_INCLUDED
