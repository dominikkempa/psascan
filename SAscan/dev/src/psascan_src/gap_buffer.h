#ifndef __PSASCAN_SRC_GAP_BUFFER_H_INCLUDED
#define __PSASCAN_SRC_GAP_BUFFER_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <queue>
#include <mutex>
#include <condition_variable>


namespace psascan_private {

template<typename value_type>
struct gap_buffer {  
  gap_buffer(long size_bytes, long n_increasers)
      : m_filled(0L),
        m_size(size_bytes / sizeof(value_type)) {
    m_content = new value_type[m_size];

    sblock_size = new long[n_increasers];
    sblock_beg = new long[n_increasers];
  }
  
  ~gap_buffer() {
    delete[] m_content;
    delete[] sblock_size;
    delete[] sblock_beg;
  }

  long m_filled, m_size;
  value_type *m_content;

  long *sblock_size;
  long *sblock_beg;
};

// Same class for the poll of empty and full gap buffers.  
template<typename value_type>
struct gap_buffer_poll {
  typedef gap_buffer<value_type> gap_buffer_type;

  gap_buffer_poll(long worker_threads = 0L) {
    m_worker_threads = worker_threads; // unused for the poll of empty buffers.
    m_worker_threads_finished = 0L;
  }
  
  void add(gap_buffer_type *b) {
    m_queue.push(b);
  }
  
  bool available() const {
    return m_queue.size() > 0;
  }

  gap_buffer_type *get() {
    if (m_queue.empty()) {
      fprintf(stderr, "\nError: requesting a gap buffer from empty poll!\n");
      std::exit(EXIT_FAILURE);
    }

    gap_buffer_type *ret = m_queue.front();
    m_queue.pop();

    return ret;
  }

  bool finished() const {
    return m_worker_threads_finished == m_worker_threads;
  }

  void increment_finished_workers() {
    ++m_worker_threads_finished;
  }

  std::condition_variable m_cv;
  std::mutex m_mutex;

private:
  long m_worker_threads;
  long m_worker_threads_finished;  // to detect when all threads finished

  std::queue<gap_buffer_type*> m_queue;
};

}  // psascan_private

#endif  // __PSASCAN_SRC_GAP_BUFFER_H_INCLUDED
