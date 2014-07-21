#ifndef __BUFFER_H_INCLUDED
#define __BUFFER_H_INCLUDED

#include <cstdio> // fprintf
#include <cstdlib> // exit

#include <queue>
#include <condition_variable>
#include <mutex>

const int n_increasers = 24;

template<typename T>
struct buffer {  
  buffer(long size_bytes)
      : m_filled(0L),
        m_size(size_bytes / sizeof(T)) {
    m_content = new T[m_size];

    sblock_size = new long[n_increasers];
    sblock_beg = new long[n_increasers];
  }
  
  ~buffer() {
    delete[] m_content;
    delete[] sblock_size;
    delete[] sblock_beg;
  }

  long m_filled, m_size;
  T *m_content;

  long *sblock_size;
  long *sblock_beg;
};

// Same class for the poll of empty and full buffers.  
template<typename T>
struct buffer_poll {
  buffer_poll(long worker_threads = 0L) {
    m_worker_threads = worker_threads; // unused for the poll of empty buffers.
    m_worker_threads_finished = 0L;
  }
  
  void add(buffer<T> *b) {
    m_queue.push(b);
  }
  
  bool available() const {
    return m_queue.size() > 0;
  }

  buffer<T> *get() {
    if (m_queue.empty()) {
      fprintf(stderr, "Error: requesting a buffer from empty poll!\n");
      std::exit(EXIT_FAILURE);
    }

    buffer<T> *ret = m_queue.front();
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
  long m_worker_threads; // used to detect that worker threads are done
  long m_worker_threads_finished;

  std::queue<buffer<T>* > m_queue;
};

#endif
