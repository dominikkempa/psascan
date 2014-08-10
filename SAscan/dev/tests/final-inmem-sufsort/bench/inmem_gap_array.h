#ifndef __INMEM_GAP_ARRAY_H_INCLUDED
#define __INMEM_GAP_ARRAY_H_INCLUDED

#include <vector>
#include <algorithm>
#include <mutex>

#include "parallel_utils.h"


struct inmem_gap_array {
  inmem_gap_array(long length, long max_threads)
    : m_length(length) {
    m_count = new unsigned char[m_length];
    parallel_utils::fill(m_count, m_length, (unsigned char)0, max_threads);
  }

  ~inmem_gap_array() {
    delete[] m_count;
  }
  
  unsigned char *m_count;
  long m_length;

  std::vector<long> m_excess;
  std::mutex m_excess_mutex;
};

#endif // __INMEM_GAP_ARRAY_H_INCLUDED
