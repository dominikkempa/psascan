#ifndef __GAP_ARRAY_H_INCLUDED
#define __GAP_ARRAY_H_INCLUDED

#include <vector>
#include <algorithm>

const long kExcessThreshold = 256;

struct gap_array {
  gap_array(long length)
    : m_length(length) {
    m_count = new unsigned char[m_length];
    std::fill(m_count, m_count + m_length, 0);
  }

  void increment(long i) {
    ++m_count[i];
    // Efficient, for final version.
    if (m_count[i] == 0)
      m_excess.push_back(i);

    // For testing purposes.
    // if (m_count[i] == kExcessThreshold) {
    //   m_count[i] = 0;
    //   m_excess.push_back(i);
    // }
  }

  ~gap_array() {
    delete[] m_count;
  }
  
  unsigned char *m_count;
  long m_length;
  std::vector<long> m_excess;
};

#endif // __GAP_ARRAY_H_INCLUDED
