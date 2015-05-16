#ifndef __BITVECTOR_H_INCLUDED
#define __BITVECTOR_H_INCLUDED

#include <cstdio>
#include <cstdlib>

#include <stdint.h>

#include "utils.h"
#include "io_streamer.h"


struct bitvector {
  private:
    long m_alloc_bytes;
    unsigned char *m_data;

  public:
    bitvector(std::string filename) {
      utils::read_objects_from_file<unsigned char>(m_data, m_alloc_bytes, filename);
    }

    bitvector(long length) {
      m_alloc_bytes = (length + 7) / 8;
      m_data = (unsigned char *)calloc(m_alloc_bytes, sizeof(unsigned char));
    }

    inline bool get(long i) const {
      return m_data[i >> 3] & (1 << (i & 7));
    }

    inline void set(long i) {
      m_data[i >> 3] |= (1 << (i & 7));
    }

    inline void reset(long i) {
      m_data[i >> 3] &= (~(1 << (i & 7)));
    }

    inline void flip(long i) {
      if (get(i)) reset(i);
      else set(i);
    }

    inline void save(std::string filename) const {
      utils::write_objects_to_file<unsigned char>(m_data, m_alloc_bytes, filename);
    }

    // Number of 1 bits in the range [beg..end).
    long range_sum(long beg, long end) const {
      long result = 0L;
    
      long j = beg;
      while (j < end && (j & 63))
        result += get(j++);

      uint64_t *ptr64 = (uint64_t *)(m_data + (j >> 3));
      while (j + 64 <= end) {
        result += __builtin_popcountll(*ptr64++);
        j += 64;
      }

      while (j < end)
        result += get(j++);

      return result;
    }

    ~bitvector() {
      free(m_data);
    }
};

#endif  // __BITVECTOR_H_INCLUDED

