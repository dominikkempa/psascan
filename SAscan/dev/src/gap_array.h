#ifndef __GAP_ARRAY_H_INCLUDED
#define __GAP_ARRAY_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>

#include "utils.h"
#include "stream.h"

struct buffered_gap_array {
  buffered_gap_array(long n) {
    if (n <= 0L) {
      fprintf(stderr, "Error: attempting to construct empty gap array.\n");
      std::exit(EXIT_FAILURE);
    }

    m_length = n;
    m_count = new unsigned char[m_length];
    std::fill(m_count, m_count + m_length, 0);

    m_buf = new long[k_buf_limit]; // 4MiB
    m_buf_filled = 0;

    m_excess = new long[k_excess_limit]; // 4MiB
    m_excess_filled = 0;
    m_total_excess = 0L;

    // File used to store excess values.
    storage_filename = "excess." + utils::random_string_hash();
  }

  inline void flush() {
    for (int i = 0; i < m_buf_filled; ++i) {
      long pos = m_buf[i];
      ++m_count[pos];

      if (!m_count[pos]) {
        m_excess[m_excess_filled++] = pos;
        ++m_total_excess;

        if (m_excess_filled == k_excess_limit) {
          utils::add_objects_to_file(m_excess, m_excess_filled, storage_filename);
          m_excess_filled = 0;
        }
      }
    }

    m_buf_filled = 0;
  }

  inline void increment(long i) {
    m_buf[m_buf_filled++] = i;
    
    if (m_buf_filled == k_buf_limit)
      flush();
  }

  void increment(long *buf, long elems) {
    for (int i = 0; i < elems; ++i) {
      long pos = buf[i];
      ++m_count[pos];

      if (!m_count[pos]) {
        m_excess[m_excess_filled++] = pos;
        ++m_total_excess;

        if (m_excess_filled == k_excess_limit) {
          utils::add_objects_to_file(m_excess, m_excess_filled, storage_filename);
          m_excess_filled = 0;
        }
      }
    }
  }

  ~buffered_gap_array() {
    delete[] m_count;
    delete[] m_buf;
    delete[] m_excess;

    if (utils::file_exists(storage_filename))
      utils::file_delete(storage_filename);
  }
  
  // Store to file using v-byte encoding.
  void save_to_file(std::string fname) {
    fprintf(stderr, "  Gap excess (inmem): %d\n", m_excess_filled);
    fprintf(stderr, "  Gap excess (disk): %ld\n", m_total_excess - m_excess_filled);
    fprintf(stderr, "  Saving gap to file: ");
    long double gap_save_start = utils::wclock();
    flush();

    // Gather all excess values together,
    // including ones stored on disk (if any).
    long *sorted_excess = new long[m_total_excess];
    std::copy(m_excess, m_excess + m_excess_filled, sorted_excess);
    if (m_total_excess != m_excess_filled) {
      long *dest = sorted_excess + m_excess_filled;
      long toread = m_total_excess - m_excess_filled;
      utils::read_n_objects_from_file(dest, toread, storage_filename.c_str());
    }

    // Sort the excess values.
    std::sort(sorted_excess, sorted_excess + m_total_excess);

    // Write gap values to file.
    stream_writer<unsigned char> *writer = new stream_writer<unsigned char>(fname);
    for (long j = 0, pos = 0; j < m_length; ++j) {
      long c = 0;
      while (pos < m_total_excess && sorted_excess[pos] == j)
        ++pos, ++c;

      long gap_j = m_count[j] + (c << 8);
      while (gap_j > 127) {
        writer->write((gap_j & 0x7f) | 0x80);
        gap_j >>= 7;
      }

      writer->write(gap_j);
    }
    delete writer;
    delete[] sorted_excess;
    fprintf(stderr, "%.2Lf\n", utils::wclock() - gap_save_start);
  }

  unsigned char *m_count;

  long m_length;

  static const int k_excess_limit = (1 << 19);
  long *m_excess;
  int m_excess_filled;
  long m_total_excess;
  
  static const int k_buf_limit = (1 << 19);
  long *m_buf;
  int m_buf_filled;


  std::string storage_filename;
};

#endif // __GAP_ARRAY_H_INCLUDED

