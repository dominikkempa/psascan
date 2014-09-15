#ifndef __INMEM_SASCAN_PAGEARRAY_H_INCLUDED
#define __INMEM_SASCAN_PAGEARRAY_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <stack>
#include <algorithm>
#include <thread>
#include <mutex>

namespace inmem_sascan_private {


//==============================================================================
// The pagearray stores the array of objects of given type in a permuted form.
// Accessing element at index i (for i = 0, .., m_length - 1) for a pagearray
// a is done as follows:
//
//  a.pageindex[i >> a.pagesize_log][i & a.pagesize_mask];
//
// All addresses of pages are in the range [a.origin..a.origin + a.m_length).
// Furthermore, if m_length % pagesize is not 0 (that is, there exists only
// partially filled page) then that page is the last page in the index. Also,
// we do not assume anything about its content following the values from the
// array that that pages has to contain. It could be garbage or memory used by
// some other array in the program. You should never read or write anything
// from there.
//==============================================================================
template<typename T, unsigned k_pagesize_log = 12U>
struct pagearray {
  static const unsigned pagesize_log = k_pagesize_log;
  static const unsigned pagesize = (1U << k_pagesize_log);
  static const unsigned pagesize_mask = (1U << k_pagesize_log) - 1;

  typedef T value_type;
  typedef pagearray<value_type, pagesize_log> pagearray_type;

  long m_length;

  value_type *m_origin;
  value_type **m_pageindex;

  // Initialize empty page array, possible it will be
  // a result of merging two page arrays.
  pagearray(value_type *origin, long length) {
    m_length = length;
    m_origin = origin;

    long n_pages = (m_length + pagesize - 1) / pagesize;
    m_pageindex = new value_type*[n_pages + 1];
  }

  // Build page array from plain array.
  pagearray(value_type *begin, value_type *end) {
    m_length = end - begin;
    m_origin = begin;

    long n_pages = (m_length + pagesize - 1) / pagesize;
    m_pageindex = new value_type*[n_pages + 1];
    for (long i = 0; i < n_pages; ++i)
      m_pageindex[i] = begin + i * pagesize;
  }

  inline value_type &operator[] (long i) const {
    return m_pageindex[i >> pagesize_log][i & pagesize_mask];
  }

  // Used only for testing.
  void random_shuffle() {
    long trimmed_length = m_length - m_length % pagesize;
    long n_full_pages = (trimmed_length / pagesize);
    for (long t = 0; t < 2 * n_full_pages; ++t) {
      long i = rand() % n_full_pages;
      long j = rand() % n_full_pages;

      // Swap the page content.
      for (long tt = 0; tt < pagesize; ++tt)
        std::swap(m_pageindex[i][tt], m_pageindex[j][tt]);

      // Update page index.
      std::swap(m_pageindex[i], m_pageindex[j]);
    }
  }

  ~pagearray() {
    if (m_pageindex)
      delete[] m_pageindex;
  }

  static void permute_to_plain_array_aux(pagearray_type &a,
      std::mutex *mutexes, long &selector, std::mutex &selector_mutex) {
    static const unsigned pagesize = (1U << pagesize_log);

    long len = a.m_length - a.m_length % pagesize;
    long n_pages = (len + pagesize - 1) >> pagesize_log;

    // Invariant: at all times, index[i] for any i points
    // to content that should be placed at i-th page of tab.
    while (true) {
      // Find starting point on some cycle.
      long start;
      while (true) {
        // Get the candidate using selector.
        std::unique_lock<std::mutex> lk(selector_mutex);
        while (selector < n_pages &&
            (a.m_pageindex[selector] == a.m_origin + (selector << pagesize_log) ||
             a.m_pageindex[selector] < a.m_origin || a.m_origin + len <= a.m_pageindex[selector]))
            ++selector;

        // Exit, if the selector does not give any candidate.
        if (selector == n_pages) {
          lk.unlock();
          return;
        }

        // Unlock selector lock, allow other threads
        // to look for candidates in the meantime.
        start = selector++;
        lk.unlock();

        // Lock a candidate page and check if it's still good.
        // If yes, keep lock and proceed to process it.
        if (mutexes[start].try_lock() &&
            a.m_pageindex[start] != a.m_origin + (start << pagesize_log) &&
            a.m_origin <= a.m_pageindex[start] && a.m_pageindex[start] < a.m_origin + len) break;
      }

      // Invariant: we have found a good candidate
      // page and have lock on mutexes[start].

      // First, we create temporary space for the
      // content of page at index[start] and move
      // the content at index[start] to that temp space.
      value_type *temp = new value_type[pagesize];
      std::copy(a.m_pageindex[start], a.m_pageindex[start] + pagesize, temp);
      std::swap(a.m_pageindex[start], temp);
      mutexes[start].unlock();

      // We now have free space at temp. Keep placing there
      // elements from the cycle and moving temp pointer.
      do {
        // Invariant: temp points to a page inside tab.
        long next = (temp - a.m_origin) >> pagesize_log;
        std::unique_lock<std::mutex> lk(mutexes[next]);
        std::copy(a.m_pageindex[next], a.m_pageindex[next] + pagesize, temp);
        std::swap(a.m_pageindex[next], temp);
        lk.unlock();
      } while (a.m_origin <= temp && temp < a.m_origin + len);
      delete[] temp;
    }
  }

  void permute_to_plain_array(long max_threads) {
    long len = m_length - m_length % pagesize;
    long n_pages = (len + pagesize - 1) >> pagesize_log;

    long selector = 0;
    std::mutex selector_mutex;
    std::mutex *mutexes = new std::mutex[n_pages];
    std::thread **threads = new std::thread*[max_threads];
  
    for (long i = 0; i < max_threads; ++i)
      threads[i] = new std::thread(permute_to_plain_array_aux,
          std::ref(*this), mutexes, std::ref(selector), std::ref(selector_mutex));

    for (long i = 0; i < max_threads; ++i) threads[i]->join();
    for (long i = 0; i < max_threads; ++i) delete threads[i];
    delete[] threads;
    delete[] mutexes;
    delete[] m_pageindex;
    m_pageindex = NULL;
  }
};

}  // namespace inmem_sascan


#endif  // __PAGEARRAY_H_INCLUDED
