/**
 * @file    src/psascan_src/inmem_psascan_src/pagearray.hpp
 * @section LICENCE
 *
 * This file is part of pSAscan v0.2.0
 * See: http://www.cs.helsinki.fi/group/pads/
 *
 * Copyright (C) 2014-2017
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

#ifndef __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_PAGEARRAY_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_PAGEARRAY_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stack>
#include <algorithm>
#include <thread>
#include <mutex>


namespace psascan_private {
namespace inmem_psascan_private {

/**
 * The paged array representation, as described in Appending B of
 *
 *   Juha Karkkainen, Peter Sanders, Stefan Burkhardt:
 *   Linear Work Suffix Array Construction.
 *   J. ACM 53(6), p. 918-936 (2006).
 **/

template<typename T, std::uint32_t k_pagesize_log = 12U>
class pagearray {
  public:
    static const std::uint32_t pagesize_log;
    static const std::uint32_t pagesize;
    static const std::uint32_t pagesize_mask;

    typedef T value_type;
    typedef pagearray<value_type, k_pagesize_log> pagearray_type;

    std::uint64_t m_length;
    std::uint64_t m_shift;

    value_type *m_origin;
    value_type **m_pageindex;

  public:

    // Initialize empty page array, possibly it will be
    // a result of merging two page arrays.
    pagearray(value_type *origin, std::uint64_t length) {
      m_length = length;
      m_origin = origin;
      m_shift = (pagesize - m_length % pagesize) % pagesize;

      std::uint64_t n_pages = (m_length + pagesize - 1) / pagesize;
      m_pageindex = new value_type*[n_pages + 1];
    }

    // Build page array from plain array.
    pagearray(value_type *begin, value_type *end) {
      m_length = end - begin;
      m_origin = begin;
      m_shift = (pagesize - m_length % pagesize) % pagesize;

      std::uint64_t n_pages = (m_length + pagesize - 1) / pagesize;
      m_pageindex = new value_type*[n_pages + 1];
      for (std::uint64_t i = 0; i < n_pages; ++i)
        m_pageindex[i] = begin + i * pagesize - m_shift;
    }

    inline value_type &operator[] (std::uint64_t i) const {
      i += m_shift;
      return m_pageindex[i >> k_pagesize_log][i & pagesize_mask];
    }

    inline std::uint64_t get_page_offset(std::uint64_t i) const {
      i += m_shift;
      return (i & pagesize_mask);
    }

    inline std::uint64_t get_page_id(std::uint64_t i) const {
      i += m_shift;
      return (i >> k_pagesize_log);
    }

    inline std::uint64_t get_page_id(value_type *p) const {
      p += m_shift;
      return ((p - m_origin) >> k_pagesize_log);
    }

    inline bool owns_page(value_type *p) const {
      p += m_shift;
      return m_origin <= p && p < m_origin + m_length;
    }

    inline value_type *get_page_addr(std::uint64_t id) const {
      return m_origin + (id << k_pagesize_log) - m_shift;
    }

    inline bool fully_contained_page(value_type *p) const {
      p += m_shift;
      return (m_origin <= p && p + pagesize <= m_origin + m_length);
    }

    // Used only for testing.
    void random_shuffle() {
      std::uint64_t trimmed_length = m_length - m_length % pagesize;
      std::uint64_t n_full_pages = (trimmed_length / pagesize);
      for (std::uint64_t t = 0; t < 2 * n_full_pages; ++t) {
        std::uint64_t i = rand() % n_full_pages;
        std::uint64_t j = rand() % n_full_pages;

        // Swap the page content.
        for (std::uint64_t tt = 0; tt < pagesize; ++tt)
          std::swap(m_pageindex[i][tt], m_pageindex[j][tt]);

        // Update page index.
        std::swap(m_pageindex[i], m_pageindex[j]);
      }
    }

    ~pagearray() {
      if (m_pageindex)
        delete[] m_pageindex;
    }

  private:
    static void permute_to_plain_array_aux(pagearray_type &a,
        std::mutex *mutexes, std::uint64_t &selector, std::mutex &selector_mutex) {
      std::uint64_t n_pages = (a.m_length + pagesize - 1) / pagesize;

      // Invariant: at all times, index[i] for any i points
      // to content that should be placed at i-th page of tab.
      while (true) {

        // Find starting point on some cycle.
        std::uint64_t start;
        while (true) {

          // Get the candidate using selector.
          std::unique_lock<std::mutex> lk(selector_mutex);
          while (selector < n_pages && a.m_pageindex[selector] == a.get_page_addr(selector))
            ++selector;

          // Exit if the selector does not give any candidate.
          if (selector == n_pages) {
            lk.unlock();
            return;
          }

          // Unlock selector lock, allow other threads
          // to look for candidates in the meantime.
          start = selector++;
          lk.unlock();

          // Lock a candidate page and check if it's still good.
          // If yes, keep the lock and exit the loop to process the page.
          if (mutexes[start].try_lock() && a.m_pageindex[start] != a.get_page_addr(start)) break;
        }

        // Invariant: we found a good candidate
        // page and have lock on mutexes[start].

        // First, create temporary space for the content
        // of page at index[start] and move the content
        // at index[start] to that temp space.
        value_type *temp = new value_type[pagesize];
        std::copy(a.m_pageindex[start], a.m_pageindex[start] + pagesize, temp);
        std::swap(a.m_pageindex[start], temp);
        mutexes[start].unlock();

        // We now have free space at temp. Keep placing there
        // elements from the cycle and moving temp pointer.
        do {

          // Invariant: temp points to a page inside tab.
          std::uint64_t next = a.get_page_id(temp);
          std::unique_lock<std::mutex> lk(mutexes[next]);
          std::copy(a.m_pageindex[next], a.m_pageindex[next] + pagesize, temp);
          std::swap(a.m_pageindex[next], temp);
          lk.unlock();
        } while (a.owns_page(temp));
        delete[] temp;
      }
    }

  public:
    void permute_to_plain_array(std::uint64_t max_threads) {
      std::uint64_t n_pages = (m_length + pagesize - 1) / pagesize;
      std::uint64_t selector = 0;

      std::mutex selector_mutex;
      std::mutex *mutexes = new std::mutex[n_pages];
      std::thread **threads = new std::thread*[max_threads];
  
      for (std::uint64_t i = 0; i < max_threads; ++i)
        threads[i] = new std::thread(permute_to_plain_array_aux,
            std::ref(*this), mutexes, std::ref(selector), std::ref(selector_mutex));

      for (std::uint64_t i = 0; i < max_threads; ++i) threads[i]->join();
      for (std::uint64_t i = 0; i < max_threads; ++i) delete threads[i];
      delete[] threads;
      delete[] mutexes;
      delete[] m_pageindex;
      m_pageindex = NULL;
    }
};

template<typename T, std::uint32_t k_pagesize_log>
  const std::uint32_t pagearray<T, k_pagesize_log>
  ::pagesize_log = k_pagesize_log;

template<typename T, std::uint32_t k_pagesize_log>
  const std::uint32_t pagearray<T, k_pagesize_log>
  ::pagesize = (1U << k_pagesize_log);

template<typename T, std::uint32_t k_pagesize_log>
  const std::uint32_t pagearray<T, k_pagesize_log>
  ::pagesize_mask = (1U << k_pagesize_log) - 1;

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_PAGEARRAY_HPP_INCLUDED
