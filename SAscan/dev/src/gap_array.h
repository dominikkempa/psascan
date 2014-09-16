#ifndef __GAP_ARRAY_H_INCLUDED
#define __GAP_ARRAY_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <mutex>
#include <algorithm>

#include "utils.h"
#include "io_streamer.h"

struct buffered_gap_array {
  buffered_gap_array(long n, unsigned char *count) {
    if (n <= 0L) {
      fprintf(stderr, "Error: attempting to construct empty gap array.\n");
      std::exit(EXIT_FAILURE);
    }

    m_length = n;

    // m_count = new unsigned char[m_length];
    m_count = count;
    std::fill(m_count, m_count + m_length, 0);

    m_excess = new long[k_excess_limit];
    m_excess_filled = 0;

    // File used to store excess values.
    m_storage_filename = "excess." + utils::random_string_hash();
    m_excess_disk = 0L;
  }

  void add_excess(long x) {
    m_excess[m_excess_filled++] = x;
    if (m_excess_filled == k_excess_limit) {
      m_gap_writing_mutex.lock();
      m_excess_disk += m_excess_filled;
      utils::add_objects_to_file(m_excess, m_excess_filled, m_storage_filename);
      m_excess_filled = 0L;
      m_gap_writing_mutex.unlock();
    }
  }

  std::mutex m_excess_mutex;
  std::mutex m_gap_writing_mutex;

  ~buffered_gap_array() {
    // delete[] m_count;
    delete[] m_excess;

    if (utils::file_exists(m_storage_filename))
      utils::file_delete(m_storage_filename);
  }
  
  // Store to file using v-byte encoding.
  void save_to_file(std::string fname) {
    long total_excess = m_excess_filled + m_excess_disk;
    // fprintf(stderr, "  Gap excess (inmem): %ld\n", m_excess_filled);
    // fprintf(stderr, "  Gap excess (disk): %ld\n", m_excess_disk);
    fprintf(stderr, "  Saving gap to file: ");
    long double gap_save_start = utils::wclock();
    
    // Gather all excess values together,
    // including ones stored on disk (if any).
    long *sorted_excess = new long[total_excess];
    long filled = m_excess_filled;
    std::copy(m_excess, m_excess + m_excess_filled, sorted_excess);
    if (m_excess_disk > 0L) {
      long *dest = sorted_excess + filled;
      long toread = m_excess_disk;
      utils::read_n_objects_from_file(dest, toread, m_storage_filename.c_str());
    }

    // Sort the excess values.
    std::sort(sorted_excess, sorted_excess + total_excess);

    // Write gap values to file.
    stream_writer<unsigned char> *writer = new stream_writer<unsigned char>(fname);
    for (long j = 0, pos = 0; j < m_length; ++j) {
      long c = 0;
      while (pos < total_excess && sorted_excess[pos] == j)
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

  static const int k_excess_limit = (1 << /*19*/25);
  long *m_excess, m_excess_filled;

  std::string m_storage_filename;
  long m_excess_disk;
};

#endif // __GAP_ARRAY_H_INCLUDED

