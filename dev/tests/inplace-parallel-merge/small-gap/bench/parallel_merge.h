/*********************************************************************
  This file implements a single method. The function is fully
  parallelized and works almost inplace.

      template<typename T, unsigned pagesize_bits>
      void merge(
                  T*          tab,
                  long        n1,
                  long        n2,
                  gap_array*  gap,
                  long        max_threads
      );

  The function takes the input array tab[0..n1+n2) and treats it as
  two subarrays: tab[0..n1) and tab[n1..n1+n2) that we want to merge.
  The ordering of the elements in the output sequence (after the
  merging) is defined by the gap array. For i = 0, .., n1 we have

      gap[i] = number of elements placed between tab[i] and tab[i - 1]
               in the output (if 0 < i < n1). gap[0] is the number of
               elements placed before tab[0]. gap[n1] is the number
               of elements placed after tab[n1 - 1].

  The function works with the small gap array representation, that is:

      gap[i] = gap->m_count[i] +
               256 * number of occurrences of i in gap->m_excess

  We assume that m_excess is sorted. We only access gap array
  sequentially, so it suffices to binary search in m_excess to find
  the starting point and then scan m_count together with m_excess to
  obtain actual gap values.

  The parameter pagesize_bits is the log_2 of used pagesize. In
  practice a good value is 12. A bigger value does not significantly
  affect speed but increases the space usage.

  The parameter max_threads is the number of threads used during
  computation. For optimal performance a number of physical cores
  should be used here. If the CPU supports hyper threading then the
  number of 2 * #cores should be used.

  Some comments about the extra space usage are in order here.
*********************************************************************/

#ifndef __PARALLEL_MERGE_H_INCLUDED
#define __PARALLEL_MERGE_H_INCLUDED

#include <cstdio>
#include <cstdlib>

#include <vector>
#include <stack>
#include <algorithm>
#include <thread>
#include <mutex>

#include "gap_array.h"


//==============================================================================
// Find the smallest j such that j + gap[0] + .. + gap[j] >= a. Store
// the value of j into b and gap[0] + .. + gap[j] into c. To speed up the
// algorithm, we have array gapsum defined as
//
//    gapsum[i] = gap[0] + .. + gap[i * block_size - 1].
//
//==============================================================================
void answer_single_gap_query(gap_array *gap, long length, long block_size,
    long *gapsum, long a, long &b, long &c) {
  long n_blocks = (length + block_size - 1) / block_size;

  // Find the block containing the correct index. To do that  find the largest
  // j such that gapsum[j] + block_size * j - 1 < a and start searching from
  // j * block_size.
  long j = 0;
  while (j + 1 < n_blocks && gapsum[j + 1] + block_size * (j + 1) - 1 < a) ++j;
  // Invariant: the j we are searching for is > j * block_size - 1.

  long sum = gapsum[j];
  j = block_size * j;
  size_t excess_ptr = std::lower_bound(gap->m_excess.begin(),
      gap->m_excess.end(), j) - gap->m_excess.begin();

  while (true) {
    // Invariant: sum = gap[0] + .. + gap[j - 1].
    // Compute gap[j] using small gap array representation.
    long gap_j = gap->m_count[j];
    while (excess_ptr < gap->m_excess.size() && gap->m_excess[excess_ptr] == j) {
      gap_j += kExcessThreshold;
      ++excess_ptr;
    }

    if (j + sum + gap_j >= a) { b = j; c = sum + gap_j; return; }
    else { sum += gap_j; ++j; }
  }
}


//==============================================================================
// Compute gap[beg] + .. + gap[end - 1] and store into result.
//==============================================================================
void compute_sum(gap_array *gap, long beg, long end, long &result) {
  result = 0;

  size_t excess_ptr = std::lower_bound(gap->m_excess.begin(),
      gap->m_excess.end(), beg) - gap->m_excess.begin();
  for (long i = beg; i < end; ++i) {
    // Compute gap[i].
    long gap_i = gap->m_count[i];
    while (excess_ptr < gap->m_excess.size() && gap->m_excess[excess_ptr] == i) {
      gap_i += kExcessThreshold;
      ++excess_ptr;
    }

    result += gap_i;
  }
}


//==============================================================================
// Parallel computation of answers to n_queries queries of the form:
// What is the smallest j such that j + gap[0] + .. + gap[j] >= a[i]"
//   - the answer to i-th query is stored in b[i]
//   - in addition we also return gap[0] + .. + gap[j] in c[i]
//
// To do that we first split the gap array into blocks and in parallel compute
// sums of gap values inside these blocks. We the accumulate these sums into
// array of prefix sums.
//
// To answer each of the queries we start a separate thread. Each thread uses
// the partial sums of gap array at block boundaries to find a good starting
// point for search and then scans the gap array from there.
//==============================================================================
void answer_gap_queries(gap_array *gap, long length, long n_queries,
    long *a, long *b, long *c, long max_threads) {
  //----------------------------------------------------------------------------
  // STEP 1: split gap array into at most blocks and in parallel compute sum
  // of values inside each block.
  //
  // To speed up STEP 3 of this algorithm we split the sequence into blocks
  // of at most 4 * 2^20 elements. This slightly increases the memory usage,
  // does not slow down the parallel sum computation and greatly speed up
  // STEP 3, since each thread has to scan only at most 4 * 2^20 elements.
  //----------------------------------------------------------------------------
  long double start = utils::wclock();
  fprintf(stderr, "Precomp1: ");
  long block_size = std::min(4L << 20, (length + max_threads - 1) / max_threads);
  long n_blocks = (length + block_size - 1) / block_size;
  long *gapsum = new long[n_blocks];
  std::thread **threads = new std::thread*[max_threads];
  for (long range_beg = 0; range_beg < n_blocks; range_beg += max_threads) {
    long range_end = std::min(range_beg + max_threads, n_blocks);
    long range_size = range_end - range_beg;

    // Compute sum inside blocks range_beg, .., range_end - 1.
    for (long i = range_beg; i < range_end; ++i) {
      long block_beg = i * block_size;
      long block_end = std::min(block_beg + block_size, length);
      threads[i - range_beg] = new std::thread(compute_sum, gap,
          block_beg, block_end, std::ref(gapsum[i]));
    }
    for (long i = 0; i < range_size; ++i) threads[i]->join();
    for (long i = 0; i < range_size; ++i) delete threads[i];
  }
  delete[] threads;
  fprintf(stderr, "%5.2Lf ", utils::wclock() - start);

  //----------------------------------------------------------------------------
  // STEP 2: compute partial sum from block counts.
  //----------------------------------------------------------------------------
  // Change gapsum so that gapsum[i] is the sum of blocks 0, 1, .., i - 1.
  fprintf(stderr, "Precomp2: ");
  start = utils::wclock();
  for (long i = 0, s = 0, t; i < n_blocks; ++i)
    { t = gapsum[i]; gapsum[i] = s; s += t; }
  fprintf(stderr, "%5.2Lf ", utils::wclock() - start);

  //----------------------------------------------------------------------------
  // STEP 3: Answer the queries in parallel.
  //----------------------------------------------------------------------------
  start = utils::wclock();
  fprintf(stderr, "Precomp3: ");
  threads = new std::thread*[n_queries];
  for (long i = 0; i < n_queries; ++i)
    threads[i] = new std::thread(answer_single_gap_query, gap, length,
      block_size, gapsum, a[i], std::ref(b[i]), std::ref(c[i]));
  for (long i = 0; i < n_queries; ++i) threads[i]->join();
  fprintf(stderr, "%5.2Lf ", utils::wclock() - start);
  for (long i = 0; i < n_queries; ++i) delete threads[i];
  delete[] threads;
  delete[] gapsum;
}


//==============================================================================
// Compute the range [res_beg..res_beg+res_size) of the output (i.e., the
// sequence after merging). The range is guaranteed to be aligned with page
// boundaries.
//==============================================================================
template<typename T, unsigned pagesize_bits>
void parallel_merge(T *tab, gap_array *gap, long length, T** pageindex, long left_idx,
    long right_idx, int remaining_gap, long res_beg, long res_size) {
  static const unsigned pagesize = (1U << pagesize_bits);
  static const unsigned pagesize_mask = pagesize - 1;

  std::stack<T*> freepages;
  T *dest = NULL;

  size_t excess_ptr = std::lower_bound(gap->m_excess.begin(),
      gap->m_excess.end(), left_idx + 1) - gap->m_excess.begin();
  long lpage_read = 0L, rpage_read = 0L, filled = 0L;
  for (long i = res_beg; i < res_beg + res_size; ++i) {
    if (!(i & pagesize_mask)) {
      if (freepages.empty()) dest = new T[pagesize];
      else { dest = freepages.top(); freepages.pop(); }
      pageindex[i >> pagesize_bits] = dest;
      filled = 0L;
    }
    if (remaining_gap > 0) {
      --remaining_gap;
      // The next element comes from the right subarray.
      dest[filled++] = tab[right_idx++];
      ++rpage_read;
      if (!(right_idx & pagesize_mask)) {
        // We reached the end of page in the right subarray.
        // We put it into free pages if we read exactly
        // pagesize elements from it. This means the no other
        // thread will attemp to read from it in the future.
        if (rpage_read == pagesize)
          freepages.push(tab + right_idx - pagesize);
        rpage_read = 0;
      }
    } else {
      // Next elem comes from the left subarray.
      dest[filled++] = tab[left_idx++];

      // Compute gap[left_idx].
      long gap_left_idx = gap->m_count[left_idx];
      while (excess_ptr < gap->m_excess.size() &&
          gap->m_excess[excess_ptr] == left_idx) {
        gap_left_idx += kExcessThreshold;
        ++excess_ptr;
      }

      remaining_gap = gap_left_idx;
      ++lpage_read;
      if (!(left_idx & pagesize_mask)) {
        // We reached the end of page in the left
        // subarray, proceed analogously.
        if (lpage_read == pagesize)
          freepages.push(tab + left_idx - pagesize);
        lpage_read = 0;
      }
    }
  }

  // Release the unused auxiliary pages.
  while (!freepages.empty()) {
    T* p = freepages.top();
    freepages.pop();
    if (p < tab || tab + length <= p)
      delete[] p;
  }
}


//==============================================================================
// A function that permutes the pages according to permutation stored in the
// given page index. The function does not deal with particular range of parts
// of cycles in the permutation, but simply chooses and element on any cycle
// and keeps following the cycle until it arrives at a page that was a starting
// point for some other threads.
//
// The function achieves correctness through non-trivial mutual exclusion
// mechanism and a series of invariants.
//==============================================================================
template<typename T, unsigned pagesize_bits>
void parallel_permute(T *tab, T** index, std::mutex *mutexes,
    long length, long n_pages, long &selector, std::mutex &selector_mutex) {
  unsigned pagesize = (1U << pagesize_bits);
  // Invariant: at all times, index[i] for any i points
  // to content that should be placed at i-th page of tab.
  while (true) {
    // Find starting point on some cycle.
    long start;
    while (true) {
      // Get the candidate using selector.
      std::unique_lock<std::mutex> lk(selector_mutex);
      while (selector < n_pages &&
          (index[selector] == tab + (selector << pagesize_bits) ||
           index[selector] < tab || tab + length <= index[selector]))
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
          index[start] != tab + (start << pagesize_bits) &&
          tab <= index[start] && index[start] < tab + length) break;
    }

    // Invariant: we have found a good candidate
    // page and have lock on mutexes[start].

    // First, we create temporary space for the
    // content of page at index[start] and move
    // the content at index[start] to that temp space.
    T *temp = new T[pagesize];
    std::copy(index[start], index[start] + pagesize, temp);
    std::swap(index[start], temp);
    mutexes[start].unlock();

    // We now have free space at temp. Keep placing there
    // elements from the cycle and moving temp pointer.
    do {
      // Invariant: temp points to a page inside tab.
      long next = (temp - tab) >> pagesize_bits;
      std::unique_lock<std::mutex> lk(mutexes[next]);
      std::copy(index[next], index[next] + pagesize, temp);
      std::swap(index[next], temp);
      lk.unlock();
    } while (tab <= temp && temp < tab + length);
    delete[] temp;
  }
}


//==============================================================================
// Merge subarray tab[0..n1) and tab[n1..n1+n2). The ordering of the output
// sequence is defined by the gap array. For i = 0, .., n1 we have
//
//       gap[i] = number of elements placed between tab[i] and tab[i - 1]
//                in the output (if i > 0) or at the beginning (if i = 0).
//
// The function is almost in-place and fully parallelized.
// Add here the exact extra space usage.
//==============================================================================
template<typename T, unsigned pagesize_bits>
void merge(T *tab, long n1, long n2, gap_array *gap, long max_threads) {
  static const unsigned pagesize = (1U << pagesize_bits);
  long length = n1 + n2;

  long n_pages = (length + pagesize - 1) / pagesize;
  T **pageindex = new T*[n_pages];

  //----------------------------------------------------------------------------
  // STEP 1: compute the initial parameters for each thread. Each thread is
  //         assigned a range of elements in the output. This range is
  //         guaranteed to be aligned with page boundaries. To start merging
  //         each thread needs to know where to start the merging in each
  //         of the two input subarrays.
  //----------------------------------------------------------------------------
  long pages_per_thread = (n_pages + max_threads - 1) / max_threads;
  long n_threads = (n_pages + pages_per_thread - 1) / pages_per_thread;

  // Allocate the arrays for initial parameters. Each thread gets:
  //   * long left_idx, right_idx -- indices to first elems from aubarrays
  //   * initial_bckt_size -- how many element from right seq goes first
  //   * res_size --  number of elements to process
  //   * res_beg --  index to the first elements of the output
  // In other words, the thread is responsible for computing the output
  // in the range [res_beg .. res_beg + res_size).
  long *left_idx = new long[n_threads];
  long *right_idx = new long[n_threads];
  int *remaining_gap = new int[n_threads];

  // Prepare gap queries.
  long *gap_query = new long[n_threads];
  long *gap_answer_a = new long[n_threads];
  long *gap_answer_b = new long[n_threads];
  for (long i = 0; i < n_threads; ++i)
    gap_query[i] = i * pages_per_thread * pagesize;

  // Answer these queries in parallel and convert the
  // answers to left_idx, right_idx and remaining_gap values.
  answer_gap_queries(gap, n1 + 1, n_threads, gap_query,
      gap_answer_a, gap_answer_b, max_threads);
  for (long i = 0; i < n_threads; ++i) {
    long res_beg = i * pages_per_thread * pagesize;
    long j = gap_answer_a[i], s = gap_answer_b[i];
    left_idx[i] = j;
    right_idx[i] = n1 + (res_beg - j);
    remaining_gap[i] = (int)(j + s - res_beg);
  }
  delete[] gap_query;
  delete[] gap_answer_a;
  delete[] gap_answer_b;

  //----------------------------------------------------------------------------
  // STEP 2: perform the actual merging.
  //----------------------------------------------------------------------------
  long double start = utils::wclock();
  fprintf(stderr, "Merging: ");
  std::thread **threads = new std::thread*[n_threads];
  for (long i = 0; i < n_threads; ++i) {
    long res_beg = i * pages_per_thread * pagesize;
    long res_end = std::min(res_beg + pages_per_thread * pagesize, length);
    long res_size = res_end - res_beg;

    threads[i] = new std::thread(parallel_merge<T, pagesize_bits>,
      tab, gap, length, pageindex, left_idx[i], right_idx[i],
      remaining_gap[i], res_beg, res_size);
  }
  for (long i = 0; i < n_threads; ++i) threads[i]->join();
  for (long i = 0; i < n_threads; ++i) delete threads[i];
  delete[] threads;
  delete[] left_idx;
  delete[] right_idx;
  delete[] remaining_gap;
  fprintf(stderr, "%5.2Lf ", utils::wclock() - start);

  start = utils::wclock();
  fprintf(stderr, "Erase aux pages: ");
  // If the last input page was incomplete, handle
  // it separatelly and exclude from the computation.
  if (length % pagesize) {
    long size = length % pagesize;
    T *lastpage = pageindex[n_pages - 1];
    T *dest = tab + pagesize * (n_pages - 1);
    std::copy(lastpage, lastpage + size, dest);

    // Release the lastpage if it was temporary.
    if (lastpage < tab || tab + length <= lastpage)
      delete[] lastpage;

    length -= size;
    --n_pages;
  }

  // Find unused input pages.
  bool *usedpage = new bool[n_pages];
  std::fill(usedpage, usedpage + n_pages, false);
  std::vector<std::pair<long, T*> > auxpages;
  for (long i = 0; i < n_pages; ++i) {
    T *p = pageindex[i];
    if (p < tab || tab + length <= p)
      auxpages.push_back(std::make_pair(i, p));
    else usedpage[(p - tab) >> pagesize_bits] = true;
  }

  // Assign aux pages to unused pages in any
  // order and release them (aux pages).
  for (long i = 0, ptr = 0; i < n_pages; ++i) {
    if (!usedpage[i]) {
      long page_id = auxpages[ptr].first;
      T *src = auxpages[ptr++].second;
      T *dest = tab + i * pagesize;
      std::copy(src, src + pagesize, dest);
      pageindex[page_id] = dest;
      delete[] src;
    }
  }
  delete[] usedpage;
  fprintf(stderr, "%5.2Lf ", utils::wclock() - start);

  //----------------------------------------------------------------------------
  // STEP 3: parallel permutation of pages.
  //----------------------------------------------------------------------------
  start = utils::wclock();
  fprintf(stderr, "Permuting: ");
  long selector = 0;
  std::mutex selector_mutex;
  std::mutex *mutexes = new std::mutex[n_pages];
  threads = new std::thread*[max_threads];
  
  for (long i = 0; i < max_threads; ++i)
    threads[i] = new std::thread(parallel_permute<T, pagesize_bits>,
        tab, pageindex, mutexes, length, n_pages,
        std::ref(selector), std::ref(selector_mutex));

  for (long i = 0; i < max_threads; ++i) threads[i]->join();
  fprintf(stderr, "%5.2Lf\n", utils::wclock() - start);
  for (long i = 0; i < max_threads; ++i) delete threads[i];
  delete[] threads;
  delete[] mutexes;
  delete[] pageindex;
}

#endif  // __PARALLEL_MERGE_H_INCLUDED
