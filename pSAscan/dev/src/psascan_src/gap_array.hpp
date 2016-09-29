/**
 * @file    psascan_src/gap_array.hpp
 * @section LICENCE
 *
 * This file is part of pSAscan v0.2.0
 * See: http://www.cs.helsinki.fi/group/pads/
 *
 * Copyright (C) 2014-2016
 *   Juha Karkkainen <juha.karkkainen (at) cs.helsinki.fi>
 *   Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/

#ifndef __PSASCAN_SRC_GAP_ARRAY_HPP_INCLUDED
#define __PSASCAN_SRC_GAP_ARRAY_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <mutex>
#include <string>
#include <thread>
#include <algorithm>
#include <parallel/algorithm>
#include <omp.h>

#include "utils.hpp"
#include "bitvector.hpp"
#include "io/async_stream_writer.hpp"


namespace psascan_private {

struct buffered_gap_array {
  buffered_gap_array(std::uint64_t length, std::string storage_fname = std::string("")) {
    if (length <= 0L) {
      fprintf(stderr, "\nError: attempting to construct empty gap array.\n");
      std::exit(EXIT_FAILURE);
    }

    m_length = length;
    m_count = (std::uint8_t *)malloc(m_length);  // XXX calloc?
    std::fill(m_count, m_count + m_length, 0);

    m_excess = new std::uint64_t[k_excess_limit];

    // File used to store excess values.
    m_storage_filename = storage_fname;
    if (!m_storage_filename.length())
      m_storage_filename = ".excess." + utils::random_string_hash();

    m_excess_filled = 0;
    m_excess_disk = 0;
    m_sorted_excess = NULL;
    m_sequential_read_initialized = false;
  }

  void add_excess(std::uint64_t x) {
    m_excess[m_excess_filled++] = x;
    if (m_excess_filled == k_excess_limit) {
      m_gap_writing_mutex.lock();  // XXX necessary?
      m_excess_disk += m_excess_filled;
      std::FILE *f = utils::file_open(m_storage_filename, "a");
      utils::write_to_file(m_excess, m_excess_filled, f);
      std::fclose(f);
      m_excess_filled = 0;
      m_gap_writing_mutex.unlock();
    }
  }

  void flush_excess_to_disk() {
    if (m_excess_filled > 0) {
      std::FILE *f = utils::file_open(m_storage_filename, "a");
      utils::write_to_file(m_excess, m_excess_filled, f);
      std::fclose(f);
      m_excess_disk += m_excess_filled;
      m_excess_filled = 0;
    }
  }

  void start_sequential_access() {
    if (!m_sequential_read_initialized) {
      m_sequential_read_initialized = true;
      m_total_excess = m_excess_filled + m_excess_disk;
      m_sorted_excess = (std::uint64_t *)malloc(m_total_excess * sizeof(std::uint64_t));  // XXX what if m_total_excess is 0?
      std::copy(m_excess, m_excess + m_excess_filled, m_sorted_excess);
      if (m_excess_disk > 0) {
        std::uint64_t *dest = m_sorted_excess + m_excess_filled;
        std::uint64_t toread = m_excess_disk;
        utils::read_from_file(dest, toread, m_storage_filename.c_str());
      }
      std::sort(m_sorted_excess, m_sorted_excess + m_total_excess);
    }

    m_excess_ptr = 0;
    m_current_pos = 0;
  }

  inline std::uint64_t get_next() {
    std::uint64_t c = 0;
    while (m_excess_ptr < m_total_excess && m_sorted_excess[m_excess_ptr] == m_current_pos)
      ++m_excess_ptr, ++c;
    std::uint64_t result = c * 256L + m_count[m_current_pos];

    ++m_current_pos;
    return result;
  }

  void stop_sequential_access() {
    if (m_sequential_read_initialized) {
      free(m_sorted_excess);
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

    free(m_count);
    delete[] m_excess;
//    if (utils::file_exists(m_storage_filename))
//      utils::file_delete(m_storage_filename);
  }

  void erase_disk_excess() {
    if (utils::file_exists(m_storage_filename))
      utils::file_delete(m_storage_filename);
  }
  
  // Write to a given file using v-byte encoding.
  void save_to_file(std::string fname) {
    fprintf(stderr, "    Write gap to file: ");
    long double gap_write_start = utils::wclock();
    std::uint64_t bytes_written = 0;

    start_sequential_access();
    typedef async_stream_writer<std::uint8_t> stream_writer_type;
    stream_writer_type *writer = new stream_writer_type(fname);

    for (std::uint64_t j = 0; j < m_length; ++j) {
      std::uint64_t val = get_next();
      while (val > 127) {
        writer->write((val & 0x7F) | 0x80);
        val >>= 7;
        ++bytes_written;
      }
      writer->write(val);
    }

    bytes_written += m_length;
    stop_sequential_access();
    delete writer;

    long double gap_write_time = utils::wclock() - gap_write_start;
    long double io_speed = (bytes_written / (1024.L * 1024)) / gap_write_time;
    fprintf(stderr, "%.2Lfs (%.2LfMiB/s)\n", gap_write_time, io_speed);
  }
  
  
  //==============================================================================
  // Note about the input:
  // - j is the maximal integer such that gapsum[j] + j <= beg.
  // - S contains value gapsum[j] + j.
  //==============================================================================
  static void convert_gap_to_bitvector_aux(std::uint64_t beg, std::uint64_t end,
      std::uint64_t j, std::uint64_t S, buffered_gap_array *gap, bitvector *bv) {
    // Initialize pointer to sorted excess values.
    std::uint64_t excess_pointer = std::lower_bound(gap->m_sorted_excess,
        gap->m_sorted_excess + gap->m_total_excess, j) - gap->m_sorted_excess;

    // Compute gap[j].
    std::uint64_t gap_j = gap->m_count[j];
    while (excess_pointer < gap->m_total_excess && gap->m_sorted_excess[excess_pointer] == j) {
      gap_j += 256;
      ++excess_pointer;
    }

    std::uint64_t p = beg;
    std::uint64_t ones = std::min(end - p, (gap_j + S) - beg);
    for (std::uint64_t k = 0; k < ones; ++k) bv->set(p++);
    ++j;

    while (p < end) {
      ++p;

      // Compute gap[j].
      gap_j = gap->m_count[j];
      while (excess_pointer < gap->m_total_excess && gap->m_sorted_excess[excess_pointer] == j) {
        gap_j += 256;
        ++excess_pointer;
      }

      ones = std::min(end - p, gap_j);

      for (std::uint64_t k = 0; k < ones; ++k) bv->set(p++);
      ++j;
    }
  }

  static void compute_j_aux(std::uint64_t range_beg, std::uint64_t n_chunks,
      std::uint64_t max_chunk_size, const std::uint64_t *sparse_gapsum,
      std::uint64_t &initial_gap_ptr, std::uint64_t &initial_gapsum_value,
      const buffered_gap_array *gap) {
    // Fast forward through as many chunks as possible.
    std::uint64_t j = 0L;
    std::uint64_t gapsum_j = 0L;  // At any time gapsum_j = gap[0] + .. + gap[j - 1].
    while (j + 1 < n_chunks && sparse_gapsum[j + 1] + (max_chunk_size * (j + 1)) <= range_beg) ++j;
    gapsum_j = sparse_gapsum[j];
    j = (j * max_chunk_size);

    // Slowly find the right place in a single chunk.
    std::uint64_t excess_ptr = std::lower_bound(gap->m_sorted_excess,
        gap->m_sorted_excess + gap->m_total_excess, j) - gap->m_sorted_excess;
    while (j < gap->m_length) {
      std::uint64_t gap_j = gap->m_count[j];
      while (excess_ptr < gap->m_total_excess && gap->m_sorted_excess[excess_ptr] == j) {
        gap_j += 256;
        ++excess_ptr;
      }

      if (gapsum_j + gap_j + j + 1 <= range_beg) {
        gapsum_j += gap_j;
        ++j;
      } else break;
    }

    // Store the answer.
    initial_gap_ptr = j;
    initial_gapsum_value = gapsum_j + j;
  }


  static void compute_gapsum_for_chunk_group(std::uint64_t group_beg, std::uint64_t group_end,
      std::uint64_t max_chunk_size, std::uint64_t *sparse_gapsum, const buffered_gap_array *gap) {
    for (std::uint64_t chunk_id = group_beg; chunk_id < group_end; ++chunk_id) {
      std::uint64_t chunk_beg = chunk_id * max_chunk_size;
      std::uint64_t chunk_end = std::min(chunk_beg + max_chunk_size, gap->m_length);

      // Compute sum of gap values inside chunk. We assume that
      // the excess values are in RAM and were sorted.
      std::int64_t occ = std::upper_bound(gap->m_sorted_excess, gap->m_sorted_excess + gap->m_total_excess, chunk_end - 1)
        - std::lower_bound(gap->m_sorted_excess, gap->m_sorted_excess + gap->m_total_excess, chunk_beg);
      std::uint64_t gap_sum_inside_chunk = (std::uint64_t)256 * (std::uint64_t)std::max((std::int64_t)0, occ);
      for (std::uint64_t j = chunk_beg; j < chunk_end; ++j)
        gap_sum_inside_chunk += gap->m_count[j];

      // Store the result.
      sparse_gapsum[chunk_id] = gap_sum_inside_chunk;
    }
  }


  //==============================================================================
  //
  // Given a gap array computed for substring of length left_block_size wrt
  // to a range of suffixes of length right_block_size, compute a bitvector of
  // length left_block_size + right_block_size that has 0 if the suffix starts in
  // the left block and 1 if it start in the right block.
  //
  //==============================================================================
  //
  // This should have the following description. Compute the bitvector
  // representation of gap array. Note that this is well defined for any gap
  // array, regardless of whether it was indeed computed wrt tail or some subset
  // of suffixes.
  //
  // Nevertheless, to optimize the computation, we provie the sum of values
  // in the gap array. This allows to immediately allocate the bitvector of
  // correct size.
  //
  // "Compute the bitvector representation of the gap array in parallel".
  // 
  //==============================================================================
  bitvector* convert_to_bitvector(std::uint64_t max_threads) {
    // 1
    //
    // The term chunks is used to compute sparse gapsum array.
    // Chunk is a length such that
    // gapsum[k] = gap[0] + gap[1] + .. + gap[k * max_chunk_size - 1]
    std::uint64_t max_chunk_size = std::min((std::uint64_t)(4 << 20), (m_length + max_threads - 1) / max_threads);
    std::uint64_t n_chunks = (m_length + max_chunk_size - 1) / max_chunk_size;
    std::uint64_t *sparse_gapsum = (std::uint64_t *)malloc(n_chunks * sizeof(std::uint64_t));


    // 2
    //
    // Compute the sum of gap value inside each chunk. Since there can be
    // more chunks than threads, we split chunks into groups and let each
    // thread compute the sum of gap values inside the group of chunks.
    std::uint64_t chunk_group_size = (n_chunks + max_threads - 1) / max_threads;
    std::uint64_t n_chunk_groups = (n_chunks + chunk_group_size - 1) / chunk_group_size;

    start_sequential_access();
    std::thread **threads = new std::thread*[n_chunk_groups];
    for (std::uint64_t t = 0; t < n_chunk_groups; ++t) {
      std::uint64_t chunk_group_beg = t * chunk_group_size;
      std::uint64_t chunk_group_end = std::min(chunk_group_beg + chunk_group_size, n_chunks);

      threads[t] = new std::thread(compute_gapsum_for_chunk_group, chunk_group_beg,
          chunk_group_end, max_chunk_size, sparse_gapsum, this);
    }

    for (std::uint64_t t = 0; t < n_chunk_groups; ++t) threads[t]->join();
    for (std::uint64_t t = 0; t < n_chunk_groups; ++t) delete threads[t];
    delete[] threads;


    // 3
    //
    // Compute comulative sum over sparse_gapsum array.
    std::uint64_t gap_total_sum = 0;
    for (std::uint64_t i = 0; i < n_chunks; ++i) {
      std::uint64_t temp = sparse_gapsum[i];
      sparse_gapsum[i] = gap_total_sum;
      gap_total_sum += temp;
    }


    // 4
    //
    // Compute all initial gap pointers. For a thread handling range [beg..end), the
    // initial_gap_ptr values is the largest j, such that gapsum[j] + j <= beg.
    // After we find j, we store the value of gapsum[j] + j in initial_gapsum_value.
    std::uint64_t result_length = (m_length + gap_total_sum) - 1;
    bitvector *result = new bitvector(result_length + 1);  // +1 is to make room for sentinel

    std::uint64_t max_range_size = (result_length + max_threads - 1) / max_threads;
    while (max_range_size & 7) ++max_range_size;
    std::uint64_t n_ranges = (result_length + max_range_size - 1) / max_range_size;

    std::uint64_t *initial_gap_ptr = new std::uint64_t[n_ranges];  // std::vector XXX
    std::uint64_t *initial_gapsum_value = new std::uint64_t[n_ranges]; // std::vector XXX

    threads = new std::thread*[n_ranges];
    for (std::uint64_t t = 0; t < n_ranges; ++t) {
      std::uint64_t range_beg = t * max_range_size;
      threads[t] = new std::thread(compute_j_aux, range_beg, n_chunks, max_chunk_size,
          sparse_gapsum, std::ref(initial_gap_ptr[t]), std::ref(initial_gapsum_value[t]), this);
    }
    for (std::uint64_t t = 0; t < n_ranges; ++t) threads[t]->join();
    for (std::uint64_t t = 0; t < n_ranges; ++t) delete threads[t];


    // 5
    //
    // Compute the bitvector. Each thread fills in the range of bits.
    for (std::uint64_t t = 0; t < n_ranges; ++t) {
      std::uint64_t range_beg = t * max_range_size;
      std::uint64_t range_end = std::min(range_beg + max_range_size, result_length);

      threads[t] = new std::thread(convert_gap_to_bitvector_aux, range_beg,
          range_end, initial_gap_ptr[t], initial_gapsum_value[t], this, result);
    }

    for (std::uint64_t t = 0; t < n_ranges; ++t) threads[t]->join();
    for (std::uint64_t t = 0; t < n_ranges; ++t) delete threads[t];
    delete[] threads;

    delete[] initial_gap_ptr;
    delete[] initial_gapsum_value;
    stop_sequential_access();
    free(sparse_gapsum);

    return result;
  }
  
  static const std::uint64_t k_excess_limit = (std::uint64_t)(1 << 22); // XXX: isn't that too big?
                                                                        // that surely causes the problems with swapping.

  std::uint8_t *m_count;
  std::uint64_t m_length;
  std::uint64_t m_excess_filled;
  std::uint64_t m_excess_disk;
  std::uint64_t *m_excess;

  std::string m_storage_filename;

  bool m_sequential_read_initialized;
  std::uint64_t m_excess_ptr;
  std::uint64_t m_current_pos;

public:
  std::uint64_t *m_sorted_excess;
  std::uint64_t m_total_excess; 
};


struct gap_array_2n {
  gap_array_2n(const buffered_gap_array *gap) {
    m_length = gap->m_length;
    m_count = (std::uint16_t *)malloc(m_length * sizeof(std::uint16_t));
#ifdef _OPENMP
    #pragma omp parallel for
    for (std::uint64_t i = 0; i < m_length; ++i)
      m_count[i] = gap->m_count[i];
#else
    for (std::uint64_t i = 0; i < m_length; ++i)
      m_count[i] = gap->m_count[i];
#endif
    m_storage_filename = gap->m_storage_filename;
    m_excess_disk = gap->m_excess_disk;
  }

  gap_array_2n(std::uint64_t length) {
    m_length = length;
    m_count = (std::uint16_t *)malloc(m_length * sizeof(std::uint16_t));
  }

  ~gap_array_2n() {
    if (m_count)
      free(m_count);
  }

  static void apply_excess_aux(gap_array_2n *gap, const std::uint64_t *tab,
      std::uint64_t block_beg, std::uint64_t block_end,
      uint64_t &initial_run_length) {
    std::uint64_t block_size = block_end - block_beg;

    // Each thread gathers excess values in a buffer and at the end
    // copies then to the gap array's mutex-protected m_excess vector.
    std::vector<std::uint64_t> excess_buffer;

    // Compute the length of initial run.
    initial_run_length = 1UL;
    while (initial_run_length < (uint64_t)block_size && tab[block_beg] ==
        tab[block_beg + initial_run_length]) ++initial_run_length;

    // Update count values.
    for (std::uint64_t i = block_beg + initial_run_length; i < block_end; ++i) {
      std::uint64_t x = tab[i];
      uint64_t value = (uint64_t)gap->m_count[x] + (std::uint64_t)256;
      if (value >= (std::uint64_t)(1 << 16)) {
        value -= (std::uint64_t)(1 << 16);
        excess_buffer.push_back(x);
      }
      gap->m_count[x] = value;
    }

    // Copy the excess values to the gap array's mutex-protected vector.
    std::unique_lock<std::mutex> lk(gap->m_excess_mutex);
    for (std::uint64_t i = 0; i < excess_buffer.size(); ++i)
      gap->m_excess.push_back(excess_buffer[i]);
    lk.unlock();
  }

  void apply_excess_from_disk(std::uint64_t ram_budget,
      std::uint64_t max_threads) {
    if (!m_excess_disk) return;

    // We only use half of the RAM for buffer, because we will use parallel
    // merge sort for sorting the buffer (which requires double the space
    // for the input).
    std::uint64_t elems = std::max((std::uint64_t)1, ram_budget / 16);
    std::uint64_t *buffer = (std::uint64_t *)malloc(elems * sizeof(std::uint64_t));

    std::FILE *f = utils::file_open(m_storage_filename.c_str(), "r");
    std::thread **threads = new std::thread*[max_threads];

    // After sorting the buffer, when we split it equally between threads
    // we obey the rule, the every thread only counts the number of 
    // elements equal to the first element in the handled range, but does
    // not do any updates for these elements. This prevents two threads
    // trying to update the same elements in the m_count array. The
    // length of the first run is computed and returned by each thread.
    // It is then updated sequentially.
    uint64_t *first_run_length = new uint64_t[max_threads];
 
    while (m_excess_disk > 0) {
      // Read a portion of excess values from disk.
      std::uint64_t toread = std::min(m_excess_disk, elems);
      utils::read_from_file(buffer, toread, f);

      // Sort excess values in parallel.
      __gnu_parallel::sort(buffer, buffer + toread);

      // Update m_count and m_excess with elements from the buffer.
      // The buffer is dividied into blocks, each blocks handles one
      // block. Each thread updates the values except the first run
      // of the block, which is handled separatelly (sequentially).
      std::uint64_t max_block_size = (toread + max_threads - 1) / max_threads;
      std::uint64_t n_blocks = (toread + max_block_size - 1) / max_block_size;

      for (std::uint64_t t = 0; t < n_blocks; ++t) {
        std::uint64_t block_beg = t * max_block_size;
        std::uint64_t block_end = std::min(block_beg + max_block_size, toread);

        threads[t] = new std::thread(apply_excess_aux, this, buffer,
            block_beg, block_end, std::ref(first_run_length[t]));
      }

      for (std::uint64_t t = 0; t < n_blocks; ++t) threads[t]->join();
      for (std::uint64_t t = 0; t < n_blocks; ++t) delete threads[t];

      // Sequentially handle the elements in the first run of each block.
      for (std::uint64_t t = 0; t < n_blocks; ++t) {
        std::uint64_t block_beg = t * max_block_size;
        std::uint64_t first = buffer[block_beg];  // first elements in the block

        uint64_t freq = m_count[first] + first_run_length[t] * (std::uint64_t)256;
        while (freq >= (std::uint64_t)(1 << 16)) {
          freq -= (std::uint64_t)(1 << 16);
          m_excess.push_back(first);
        }
        m_count[first] = freq;
      }

      m_excess_disk -= toread;
    }

    __gnu_parallel::sort(m_excess.begin(), m_excess.end());

    delete[] threads;
    delete[] first_run_length;

    std::fclose(f);
    free(buffer);
  }

  void set_count(std::uint64_t pos, std::uint64_t value) {
    while (value >= (std::uint64_t)(1 << 16)) {
      m_excess.push_back(pos);
      value -= (std::uint64_t)(1 << 16);
    }
    m_count[pos] = value;
  }

  void erase_disk_excess() {
    if (utils::file_exists(m_storage_filename))
      utils::file_delete(m_storage_filename);
  }

  std::uint16_t *m_count;

  std::uint64_t m_length;
  std::uint64_t m_excess_disk;

  std::mutex m_excess_mutex;
  std::string m_storage_filename;
  std::vector<std::uint64_t> m_excess;  // all excess values are in RAM
};

}  // psascan_private

#endif  // __PSASCAN_SRC_GAP_ARRAY_HPP_INCLUDED
