#ifndef __GAP_ARRAY_H_INCLUDED
#define __GAP_ARRAY_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <mutex>
#include <algorithm>
#include <string>

#include "utils.h"
#include "io_streamer.h"

struct buffered_gap_array {
  buffered_gap_array(long n, unsigned char *count,
      std::string storage_fname = std::string("")) {
    if (n <= 0L) {
      fprintf(stderr, "\nError: attempting to construct empty gap array.\n");
      std::exit(EXIT_FAILURE);
    }

    m_length = n;
    m_count = count;
    std::fill(m_count, m_count + m_length, 0);

    m_excess = new long[k_excess_limit];

    // File used to store excess values.
    m_storage_filename = storage_fname;
    if (!m_storage_filename.length())
      m_storage_filename = "excess." + utils::random_string_hash();

    m_excess_filled = 0L;
    m_excess_disk = 0L;
    m_sorted_excess = NULL;
    m_sequential_read_initialized = false;
  }

  void add_excess(long x) {
    m_excess[m_excess_filled++] = x;
    if (m_excess_filled == k_excess_limit) {
      m_gap_writing_mutex.lock();  // XXX necessary?
      m_excess_disk += m_excess_filled;
      utils::add_objects_to_file(m_excess, m_excess_filled, m_storage_filename);
      m_excess_filled = 0L;
      m_gap_writing_mutex.unlock();
    }
  }

  void start_sequential_access() {
    if (!m_sequential_read_initialized) {
      m_sequential_read_initialized = true;
      m_total_excess = m_excess_filled + m_excess_disk;
      m_sorted_excess = new long[m_total_excess];
      std::copy(m_excess, m_excess + m_excess_filled, m_sorted_excess);
      if (m_excess_disk > 0L) {
        long *dest = m_sorted_excess + m_excess_filled;
        long toread = m_excess_disk;
        utils::read_n_objects_from_file(dest, toread, m_storage_filename.c_str());
      }
      std::sort(m_sorted_excess, m_sorted_excess + m_total_excess);
    }

    m_excess_ptr = 0;
    m_current_pos = 0;
  }

  inline long get_next() {
    long c = 0;
    while (m_excess_ptr < m_total_excess && m_sorted_excess[m_excess_ptr] == m_current_pos)
      ++m_excess_ptr, ++c;
    long result = c * 256L + m_count[m_current_pos];

    ++m_current_pos;
    return result;
  }

  void stop_sequential_access() {
    if (m_sequential_read_initialized) {
      delete[] m_sorted_excess;
      m_sequential_read_initialized = false;
    } else {
      fprintf(stderr, "\nError: attempting to stop sequential "
          "access to the gap array before it was initialized.\n");
      std::exit(EXIT_FAILURE);
    }
  }

  std::mutex m_excess_mutex;
  std::mutex m_gap_writing_mutex;

  ~buffered_gap_array() {
    if (m_sequential_read_initialized) {
      fprintf(stderr, "\nError: sequential access to gap was not terminated.");
      std::exit(EXIT_FAILURE);
    }

    delete[] m_excess;
    if (utils::file_exists(m_storage_filename))
      utils::file_delete(m_storage_filename);
  }
  
  // Write to a given file using v-byte encoding.
  void save_to_file(std::string fname) {
    fprintf(stderr, "  Saving gap to file: ");
    long double gap_save_start = utils::wclock();
    start_sequential_access();
    stream_writer<unsigned char> *writer = new stream_writer<unsigned char>(fname);
    for (long j = 0; j < m_length; ++j) {
      long val = get_next();
      while (val > 127) {
        writer->write((val & 0x7f) | 0x80);
        val >>= 7;
      }
      writer->write(val);
    }
    delete writer;
    stop_sequential_access();
    fprintf(stderr, "%.2Lf\n", utils::wclock() - gap_save_start);
  }

  static const int k_excess_limit = (1 << 25); // XXX: isn't that too big?

  unsigned char *m_count;
  long m_length;
  long m_excess_filled;
  long m_excess_disk;
  long *m_excess;

  std::string m_storage_filename;

  bool m_sequential_read_initialized;
  long *m_sorted_excess;
  long m_total_excess;
  long m_excess_ptr;
  long m_current_pos;
};

#endif // __GAP_ARRAY_H_INCLUDED

