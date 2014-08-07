#ifndef __STREAM_INFO_H_INCLUDED
#define __STREAM_INFO_H_INCLUDED

#include <mutex>
#include <algorithm>

#include "utils.h"

//=============================================================================
// Used to store progress information for different threads during streaming.
//=============================================================================
struct stream_info {
  stream_info(long thread_count, long tostream)
    : m_update_count(0L),
      m_thread_count(thread_count),
      m_tostream(tostream) {
    m_streamed = new long[thread_count];
    std::fill(m_streamed, m_streamed + thread_count, 0L);

    m_idle_update = new long double[thread_count];
    m_idle_work  = new long double[thread_count];
    std::fill(m_idle_update, m_idle_update + thread_count, 0.L);
    std::fill(m_idle_work, m_idle_work + thread_count, 0.L);

    m_timestamp = utils::wclock();
  }

  ~stream_info() {
    delete[] m_streamed;
    delete[] m_idle_work;
    delete[] m_idle_update;
  }

  long m_update_count;     // number of updates
  long m_thread_count;     // number of threads
  long m_tostream;         // total text length to stream
  long double m_timestamp; // when the streaming started
  long *m_streamed;        // how many bytes streamed by each thread
  long double *m_idle_update;
  long double *m_idle_work;

  std::mutex m_mutex;
};

#endif // __STREAM_INFO_H_INCLUDED
