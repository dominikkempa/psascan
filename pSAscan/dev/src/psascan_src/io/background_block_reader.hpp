/**
 * @file    psascan_src/io/background_block_reader.hpp
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

#ifndef __PSASCAN_SRC_IO_BACKGROUND_BLOCK_READER_HPP_INCLUDED
#define __PSASCAN_SRC_IO_BACKGROUND_BLOCK_READER_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <thread>
#include <mutex>
#include <condition_variable>

#include "../utils.hpp"


namespace psascan_private {

class background_block_reader {
  public:
    std::uint8_t *m_data;
    std::uint64_t m_start;
    std::uint64_t m_size;

  private:
    static const std::uint64_t k_chunk_size;

    // These variables are protected by m_mutex.
    std::uint64_t m_fetched;
    bool m_signal_stop;
    bool m_joined;

    std::mutex m_mutex;

    // This condition variable is used by the I/O thread to notify
    // the waiting threads when the next chunk is read.
    std::condition_variable m_cv;

    std::thread *m_thread;
    std::FILE *m_file;

  private:
    static void io_thread_main(background_block_reader &reader) {
      while (true) {
        std::unique_lock<std::mutex> lk(reader.m_mutex);
        std::uint64_t fetched = reader.m_fetched;
        bool signal_stop = reader.m_signal_stop;
        lk.unlock();

        if (fetched == reader.m_size || signal_stop) break;

        std::uint64_t toread = std::min(reader.m_size - fetched, reader.k_chunk_size);
        std::uint8_t *dest = reader.m_data + fetched;
        utils::read_from_file(dest, toread, reader.m_file);

        lk.lock();
        reader.m_fetched += toread;
        lk.unlock();
        reader.m_cv.notify_all();
      }
      
      // Close the file and exit.
      std::fclose(reader.m_file);
    }

  public:
    background_block_reader(std::string filename, std::uint64_t start, std::uint64_t size) {
      m_start = start;
      m_size = size;
         
      // Initialize file and buffer.
      m_data = (std::uint8_t *)malloc(m_size);
      m_file = utils::file_open(filename, "r");
      std::fseek(m_file, m_start, SEEK_SET);
      m_fetched = 0;

      // Start the I/O thread.
      m_signal_stop = false;
      m_joined = false;
      m_thread = new std::thread(io_thread_main, std::ref(*this));
    }

    ~background_block_reader() {
      if (!m_joined) {
        fprintf(stderr, "\nError: the I/O thread is still not joined when "
          "destroying an object of backgroud_block_reader.\n");
        std::exit(EXIT_FAILURE);
      }
      
      // Note: m_file is already closed.
      delete m_thread;
      free(m_data);
    }

    inline void stop() {
      // Set the flag for the thread to stop.
      std::unique_lock<std::mutex> lk(m_mutex);
      m_signal_stop = true;
      lk.unlock();

      // Wait until the thread notices the flag and exits. Possibly the thread
      // is already not running, but in this case this call will do nothing.
      m_thread->join();
      
      // To detect (in the destructor) if stop() was called.
      lk.lock();
      m_joined = true;
      lk.unlock();
    }

    inline void wait(std::uint64_t target_fetched) {
      std::unique_lock<std::mutex> lk(m_mutex);
      while (m_fetched < target_fetched)
        m_cv.wait(lk);
      lk.unlock();
    }
};

const std::uint64_t background_block_reader::k_chunk_size = (1UL << 20);

}  // namespace psascan_private

#endif  // __PSASCAN_SRC_IO_BACKGROUND_BLOCK_READER_HPP_INCLUDED
