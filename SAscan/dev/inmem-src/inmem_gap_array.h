#ifndef __INMEM_GAP_ARRAY_H_INCLUDED
#define __INMEM_GAP_ARRAY_H_INCLUDED

#include <vector>
#include <algorithm>
#include <mutex>

#include "parallel_utils.h"


struct inmem_gap_array {
  inmem_gap_array(long length, long)
    : m_length(length) {
    m_count = (unsigned char *)calloc(m_length, sizeof(unsigned char));
  }

  ~inmem_gap_array() {
    free(m_count);
  }
  
  unsigned char *m_count;
  long m_length;

  std::vector<long> m_excess;
  std::mutex m_excess_mutex;
};

#endif // __INMEM_GAP_ARRAY_H_INCLUDED
