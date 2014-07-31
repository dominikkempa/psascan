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
void parallel_merge(T *tab, int *gap, long length,
    T** pageindex, long left_idx, long right_idx,
    int remaining_gap, long res_beg, long res_size) {
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
void parallel_permute(T *tab, T** index, std::mutex *mutexes,
    long length, long n_pages, long &selector, std::mutex &selector_mutex) {
  unsigned pagesize = (1U << pagesize_bits);
  long attempts = 0;
  long double waited_total = 0.L;
  // Invariant: at all times, index[i] for any i points
  // to content that should be placed at i-th page of tab.
  while (true) {
    // Find starting point on some cycle.
    long start;
    ++attempts;
    while (true) {
      // Get the candidate using selector.
      long double st = utils::wclock();
      std::unique_lock<std::mutex> lk(selector_mutex);
      waited_total += utils::wclock() - st;
      while (selector < n_pages &&
          (index[selector] == tab + (selector << pagesize_bits) ||
           index[selector] < tab || tab + length <= index[selector]))
          ++selector;

      // Exit, if the selector does not give any candidate.
      if (selector == n_pages) {
        lk.unlock();
        fprintf(stderr, "thread waited %.2Lf sec for lock\n", waited_total);
        fprintf(stderr, "attempted begins: %ld\n", attempts);
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

// Permute pages in index[beg..end).
template<typename T, unsigned pagesize_bits>
void parallel_permute_out_of_place(T* output, T** index, long beg, long end) {
  static const unsigned pagesize = (1U << pagesize_bits);

  for (long i = beg; i < end; ++i) {
    T *src = index[i];
    T *dest = output + (i << pagesize_bits);
    std::copy(src, src + pagesize, dest);
  }
}

template<typename T, unsigned pagesize_bits>
void merge(T *tab, long n1, long n2, int *gap, long max_threads) {
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
  int *remaining_gap = new int[n_threads];
  for (long i = 0; i < n_threads; ++i) {
    long res_beg = i * pages_per_thread * pagesize;
    long j = 0, jpos = gap[0];
    while (jpos < res_beg) jpos += gap[++j] + 1;
    left_idx[i] = j;
    right_idx[i] = n1 + res_beg - j;
    remaining_gap[i] = (int)(jpos - res_beg);
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
  fprintf(stderr, "Detect and erase aux pages: ");
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



  // Compute cycle-length distribution.
  std::vector<long> cycle_lens;
  bool *visited = new bool[n_pages];
  std::fill(visited, visited + n_pages, false);
  start = utils::wclock();
  fprintf(stderr, "Traversing cycles: ");
  for (long i = 0; i < n_pages; ++i) {
    if (!visited[i]) {
      visited[i] = true;
      long cycle_len = 1L;
      long j = ((pageindex[i] - tab) >> pagesize_bits);
      while (!visited[j]) {
        ++cycle_len;
        visited[j] = true;
        j = ((pageindex[j] - tab) >> pagesize_bits);
      }
      cycle_lens.push_back(cycle_len);
    }
  }
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
  delete[] visited;
  std::sort(cycle_lens.begin(), cycle_lens.end());
  fprintf(stderr, "cycle distribution:\n");
  for (size_t beg = 0, end; beg < cycle_lens.size(); beg = end) {
    end = beg;
    while (end < cycle_lens.size() && cycle_lens[end] == cycle_lens[beg]) ++end;
    fprintf(stderr, "cycles[%ld] = %ld\n", cycle_lens[beg], end - beg);
  }


  // Parallel paermutation of pages, about 50% slower than
  // above in practice, does not require decomposition of
  // permutation into cycles.
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
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
  for (long i = 0; i < max_threads; ++i) delete threads[i];
  delete[] threads;
  delete[] mutexes;

  // Sequential-out-of-place-permutation
  /*T *output = new T[length];
  std::fill(output, output + length, (T)0);
  start = utils::wclock();
  fprintf(stderr, "Permuting (out-of-place): ");
  threads = new std::thread*[n_threads];
  for (long i = 0; i < n_threads; ++i) {
    long page_range_beg = i * pages_per_thread;
    long page_range_end = std::min(page_range_beg + pages_per_thread, n_pages);
    threads[i] = new std::thread(parallel_permute_out_of_place<T, pagesize_bits>,
        output, pageindex, page_range_beg, page_range_end);
  }
  for (long i = 0; i < n_threads; ++i) threads[i]->join();
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
  for (long i = 0; i < n_threads; ++i) delete threads[i];
  delete[] threads;
  std::copy(output, output + length, tab);
  delete[] output;*/


  delete[] pageindex;
}

#endif  // __PARALLEL_MERGE_H_INCLUDED
