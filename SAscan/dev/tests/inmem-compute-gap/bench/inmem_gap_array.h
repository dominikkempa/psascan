#ifndef __INMEM_GAP_ARRAY_H_INCLUDED
#define __INMEM_GAP_ARRAY_H_INCLUDED

#include <vector>
#include <algorithm>
#include <mutex>

struct inmem_gap_array {
  inmem_gap_array(long length)
    : m_length(length) {
    m_count = new unsigned char[m_length];
    std::fill(m_count, m_count + m_length, 0);
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
