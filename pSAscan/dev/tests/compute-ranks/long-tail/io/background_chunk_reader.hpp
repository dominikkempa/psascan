/**
 * @file    src/psascan_src/io/background_chunk_reader.hpp
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

#ifndef __SRC_PSASCAN_SRC_IO_BACKGROUND_CHUNK_READER_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_IO_BACKGROUND_CHUNK_READER_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <string>
#include <algorithm>
#include <thread>
#include <mutex>
#include <condition_variable>

#include "../utils.hpp"


namespace psascan_private {

class background_chunk_reader {
  private:
    std::FILE *m_file;
    const std::uint64_t m_end;
    const std::uint64_t m_chunk_length;
    
    std::condition_variable m_cv;
    std::mutex m_mutex;
    std::thread *m_thread;
    
    bool m_signal_read_next_chunk;
    bool m_io_thread_stopped;
    bool m_signal_stop;

    std::uint64_t m_file_pos;
    std::uint64_t m_passive_chunk_end;
    std::uint8_t *m_passive_chunk;
    std::uint8_t *m_chunk;

  private:

    //=========================================================================
    // Main function of I/O thread.
    //=========================================================================
    static void async_io_code(background_chunk_reader &r) {

      // Keep reading until explicitly told to stop
      // or reach a specified point in file.
      while (true) {

        // Wait until a message to read
        // next chunk or to terminate is set.
        std::unique_lock<std::mutex> lk(r.m_mutex);
        while (!r.m_signal_read_next_chunk && !r.m_signal_stop)
          r.m_cv.wait(lk);
          
        // Depending on the type of event
        // set the protected flags.
        const bool sig_stop = r.m_signal_stop;
        r.m_signal_read_next_chunk = false;
        if (sig_stop)
          r.m_io_thread_stopped = true;
        lk.unlock();
        
        // Exit if we were told so.
        if (sig_stop)
          break;
        
        // Compute length of next chunk.
        const std::uint64_t next_chunk_length =
          std::min(r.m_chunk_length, r.m_end - r.m_file_pos);

        // If we reached the specified
        // point in file, exit.
        if (next_chunk_length == 0) {
          lk.lock();
          r.m_io_thread_stopped = true;
          lk.unlock();
          r.m_cv.notify_all();
          break;
        }

        // Read next chunk.
        utils::read_from_file(r.m_passive_chunk,
            next_chunk_length, r.m_file);
        
        // Update the end of passive chunk
        // and notify the caller.
        lk.lock();
        r.m_file_pos += next_chunk_length;
        lk.unlock();
        r.m_cv.notify_all();
      }
    }

  public:

    //=========================================================================
    // Constructor.
    //=========================================================================
    background_chunk_reader(
        const std::string filename,
        const std::uint64_t beg,
        const std::uint64_t end,
        const std::uint64_t chunk_length = (1 << 20))
      : m_end(end),
        m_chunk_length(chunk_length),
        m_signal_read_next_chunk(true),
        m_io_thread_stopped(false),
        m_signal_stop(false),
        m_file_pos(beg),
        m_passive_chunk_end(beg) {

      // Sanity check.
      if (beg > end) {
        fprintf(stderr, "Error: beg > end in "
            "background_chunk_reader.\n");
        std::exit(EXIT_FAILURE);
      }

      // Handle special case.
      if (beg == end)
        return;

      // Allocate chunks.
      m_chunk = (std::uint8_t *)malloc(m_chunk_length);
      m_passive_chunk = (std::uint8_t *)malloc(m_chunk_length);
      
      // Open file and move position for reading.
      m_file = utils::file_open(filename, "r");
      std::fseek(m_file, m_file_pos, SEEK_SET);

      // Start the I/O thread.
      m_thread =
        new std::thread(async_io_code, std::ref(*this));
    }

    //=========================================================================
    // Read the next chunk.
    //=========================================================================
    inline void read_next_chunk() {

      // Wait until the I/O thread modifies passive_chunk_beg.
      std::unique_lock<std::mutex> lk(m_mutex);
      while (!m_io_thread_stopped && m_passive_chunk_end == m_file_pos)
        m_cv.wait(lk);

      // Sanity check.
      if (m_io_thread_stopped == true) {
        fprintf(stderr, "Error: I/O thread already stopped "
            "in background_chunk_reader.\n");
        std::exit(EXIT_FAILURE);
      }
        
      // Sanity check.
      if (m_signal_read_next_chunk) {
        fprintf(stderr, "Error: m_signal_read_next_chunk "
            "in the wrong state.\n");
        std::exit(EXIT_FAILURE);
      }

      // Swap active and passive chunk.
      std::swap(m_chunk, m_passive_chunk);
      m_passive_chunk_end = m_file_pos;
      m_signal_read_next_chunk = true;
      lk.unlock();

      // Notify the I/O thread that it
      // can now read another chunk.
      m_cv.notify_all();
    }

    //=========================================================================
    // Return const pointer to chunk.
    //=========================================================================
    inline const std::uint8_t *get_chunk_ptr() const {
      return m_chunk;
    }

    //=========================================================================
    // Destructor.
    //=========================================================================
    ~background_chunk_reader() {

      // Notify the I/O thread it should stop.
      std::unique_lock<std::mutex> lk(m_mutex);
      m_signal_stop = true;
      lk.unlock();
      m_cv.notify_all();

      // Wait until the thread notices the flag and exits.
      // Possibly the thread is already not running, but
      // in this case this call will do nothing.
      m_thread->join();
      std::fclose(m_file);

      // Clean up.  
      delete m_thread;
      free(m_chunk);
      free(m_passive_chunk);
    }

    //=========================================================================
    // Return the chunk size.
    //=========================================================================
    inline std::uint64_t get_chunk_size() const {
      return m_chunk_length;
    }
};

}  // namespace psascan_private

#endif  // __SRC_PSASCAN_SRC_IO_BACKGROUND_CHUNK_READER_HPP_INCLUDED
