/*******************************************************************************
  Generic parallel inplace merging.

  The input is an array tab[0..n) where tab[0..n1) contains elements the
  first subarray nad tab[n1..n) the second. In addition we expect an array
  gap[0..n1] such that gap[i] is the number of elements from the right subarray
  that goes in between tab[i] and tab[i - 1] (or at the beginning if i = 0).

  The procedure computes the resulting ordering of elements and places it
  in tab[0..n). It uses O(1) pages per each thread of extra space.
*******************************************************************************/

#ifndef __PARALLEL_MERGE_H_INCLUDED
#define __PARALLEL_MERGE_H_INCLUDED

#include <cstdio>
#include <cstdlib>

#include <vector>
#include <stack>
#include <algorithm>
#include <thread>
#include <mutex>

template<typename T, unsigned pagesize_bits>
void parallel_merge(T *tab, long *gap, long length,
    T** pageindex, long left_idx, long right_idx,
    long remaining_gap, long res_beg, long res_size) {
  static const unsigned pagesize = (1U << pagesize_bits);
  static const unsigned pagesize_mask = pagesize - 1;

  std::stack<T*> freepages;
  T *dest = NULL;

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
      remaining_gap = gap[left_idx];
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

template<typename T, unsigned pagesize_bits>
void merge(T *tab, long n1, long n2, long *gap, long max_threads) {
  static const unsigned pagesize = (1U << pagesize_bits);
  long length = n1 + n2;

  long n_pages = (length + pagesize - 1) / pagesize;
  T **pageindex = new T*[n_pages];

  //----------------------------------------------------------------------------
  // STEP 1: compute the initial parameters for each thread. For now, we do it
  //         seqentially.
  //
  // Each thread gets:
  //   * long left_idx, right_idx -- indices to first elems from aubarrays
  //   * initial_bckt_size -- how many element from right seq goes first
  //   * res_size --  number of elements to process
  //   * res_beg -- index to the first elements of the output
  //
  // In short, if we did it sequentially, the thread would just produce the
  // elements of the output in the range [res_beg .. rea_beg + res_size).
  //----------------------------------------------------------------------------

  long pages_per_thread = (n_pages + max_threads - 1) / max_threads;
  long n_threads = (n_pages + pages_per_thread - 1) / pages_per_thread;

  // Compute initial parameters for each thread.
  long *left_idx = new long[n_threads];
  long *right_idx = new long[n_threads];
  long *remaining_gap = new long[n_threads];
  for (long i = 0; i < n_threads; ++i) {
    long res_beg = i * pages_per_thread * pagesize;
    long j = 0, jpos = gap[0];
    while (jpos < res_beg) jpos += gap[++j] + 1;
    left_idx[i] = j;
    right_idx[i] = n1 + res_beg - j;
    remaining_gap[i] = jpos - res_beg;
  }
  
  // Okay, we can start the threads.
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
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

  start = utils::wclock();
  fprintf(stderr, "Finalizing: ");
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
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

  // Permute the pages, for now sequentially.
  /*T *output = new T[length];
  T *dest = output;
  for (long i = 0; i < n_pages; ++i, dest += pagesize)
    std::copy(pageindex[i], pageindex[i] + pagesize, dest);
  std::copy(output, output + length, tab);
  delete[] output;*/
  delete[] pageindex;
}

#endif  // __PARALLEL_MERGE_H_INCLUDED
