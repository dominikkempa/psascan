#ifndef __NOCACHE_ASYNC_STREAM_READER_H_INCLUDED
#define __NOCACHE_ASYNC_STREAM_READER_H_INCLUDED
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cerrno>

#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <algorithm>

#include "utils.h"
#include "linux_io_utils.h"
#include "aligned_alloc.h"


template<typename value_type,
  unsigned alignment = 4096 /*modify with extreme caution*/>
struct nocache_async_stream_reader {
  template<typename T, unsigned align>
  static void io_thread_code(nocache_async_stream_reader<T, align> *reader) {
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

      // Safely read the data from disk.
      reader->m_passive_buf_filled_bytes = ::read(reader->m_file_descriptor,
          reader->m_passive_buf, reader->m_buf_size_bytes);

      // Let the caller know what the I/O thread finished reading.
      lk.lock();
      reader->m_avail = false;
      lk.unlock();
      reader->m_cv.notify_all();
    }
  }

  nocache_async_stream_reader(std::string filename, long bufsize = (4L << 20)) {
    static_assert(alignment >= sizeof(value_type), "Error: alignment is too small.");

    m_file_descriptor = utils::open_file_direct(filename.c_str(), O_RDONLY);

    // Set m_buf_size_bytes so that is it a multiple of alignment
    // and m_buf_size_bytes >= max(butsize, sizeof(value_type)). 
    bufsize /= 2;
    bufsize = std::max(bufsize, (long)sizeof(value_type));
    m_buf_size_bytes = (bufsize + alignment - 1) / alignment;
    m_buf_size_bytes *= alignment;

    // Initialize counters.
    m_active_buf_filled_bytes = 0L;
    m_active_buf_filled_items = 0L;
    m_active_buf_extra_bytes = 0L;
    m_passive_buf_filled_bytes = 0L;
    m_active_buf_pos = 0L;

    // Allocate buffers.
    m_active_buf = (unsigned char *)aligned_alloc<alignment>(m_buf_size_bytes + alignment);
    m_passive_buf = (unsigned char *)aligned_alloc<alignment>(m_buf_size_bytes + alignment);
    m_active_buf += alignment;
    m_passive_buf += alignment;

    m_finished = false;

    // Start the I/O thread and immediatelly start reading.
    m_avail = true;
    m_thread = new std::thread(io_thread_code<value_type, alignment>, this);
  }

  ~nocache_async_stream_reader() {
    // Let the I/O thread know that we're done.
    std::unique_lock<std::mutex> lk(m_mutex);
    m_finished = true;
    lk.unlock();
    m_cv.notify_all();

    // Wait for the thread to actually finish.
    m_thread->join();
    
    // Clean up.
    delete m_thread;
    m_active_buf -= alignment;
    m_passive_buf -= alignment;
    aligned_dealloc(m_active_buf);
    aligned_dealloc(m_passive_buf);
    close(m_file_descriptor);
  }

  // This function checks if the reading thread has already
  // prefetched the next buffer (the request should have been
  // done before), and waits if the prefetching was not
  // completed yet.
  void receive_new_buffer() {
    // Wait until the I/O thread finishes reading the previous
    // buffer. Most of the time this step is instantaneous.
    std::unique_lock<std::mutex> lk(m_mutex);
    while (m_avail == true)
      m_cv.wait(lk);

    // Set the new active buffer.
    long extra_bytes = (m_active_buf_extra_bytes + m_active_buf_filled_bytes) % sizeof(value_type);
    if (extra_bytes > 0) {
      unsigned char *src = (m_active_buf + m_active_buf_filled_bytes) - extra_bytes;
      unsigned char *dest = m_passive_buf - extra_bytes;
      std::copy(src, src + extra_bytes, dest);
    }
    std::swap(m_active_buf, m_passive_buf);
    m_active_buf_extra_bytes = extra_bytes;
    m_active_buf_filled_bytes = m_passive_buf_filled_bytes;
    m_active_buf_filled_items = (m_active_buf_filled_bytes + m_active_buf_extra_bytes) / sizeof(value_type);
    m_active_buf_item_ptr = (value_type *)(m_active_buf - m_active_buf_extra_bytes);
    m_active_buf_pos = 0L;

    // Let the I/O thread know that it can now prefetch
    // another buffer.
    m_avail = true;
    lk.unlock();
    m_cv.notify_all();
  }

  inline value_type read() {
    if (m_active_buf_pos == m_active_buf_filled_items) {
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

    return m_active_buf_item_ptr[m_active_buf_pos++];
  }

private:
  unsigned char *m_active_buf;
  unsigned char *m_passive_buf;
  value_type *m_active_buf_item_ptr;

  long m_buf_size_bytes;
  long m_active_buf_pos;
  long m_active_buf_filled_bytes;
  long m_active_buf_filled_items;
  long m_active_buf_extra_bytes;
  long m_passive_buf_filled_bytes;

  // Used for synchronization with the I/O thread.
  std::mutex m_mutex;
  std::condition_variable m_cv;
  bool m_avail;
  bool m_finished;

  int m_file_descriptor;
  std::thread *m_thread;
};

#endif  // __ASYNC_STREAM_READER_H_INCLUDED
