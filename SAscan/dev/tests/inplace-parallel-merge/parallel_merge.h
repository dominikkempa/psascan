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

template<typename T>
void parallel_merge(T *tab, long *gap, long length, long pagesize,
    T** pagedir, long left_idx, long right_idx,
    long remaining_gap, long res_beg, long res_size,
    long &allocated_pages, long &released_pages) {

  std::stack<T*> free_pages;
  T *dest = NULL;
  allocated_pages = 0L;
  released_pages = 0L;

  long lpage_read = 0L, rpage_read = 0L, page_filled = 0L;
  for (long i = res_beg; i < res_beg + res_size; ++i) {
    if ((i % pagesize) == 0) {
      if (free_pages.empty()) {
        dest = new T[pagesize];
        ++allocated_pages;
      } else {
        dest = free_pages.top(); free_pages.pop();
      }
      pagedir[i / pagesize] = dest;
      page_filled = 0L;
    }
    if (remaining_gap > 0) {
      --remaining_gap;
      // The next element comes from the right subarray.
      dest[page_filled++] = tab[right_idx++];
      ++rpage_read;
      if ((right_idx % pagesize) == 0) {
        // We reached the end of page in the right subarray.
        // We put it into free pages if we read exactly
        // pagesize elements from it. This means the no other
        // thread will attemp to read from it in the future.
        if (rpage_read == pagesize)
          free_pages.push(tab + right_idx - pagesize);
        rpage_read = 0;
      }
    } else {
      // Next elem comes from the left subarray.
      dest[page_filled++] = tab[left_idx++];
      remaining_gap = gap[left_idx];
      ++lpage_read;
      if ((left_idx % pagesize) == 0) {
        // We reached the end of page in the left
        // subarray, proceed analogously.
        if (lpage_read == pagesize)
          free_pages.push(tab + left_idx - pagesize);
        lpage_read = 0;
      }
    }
  }

  // Release the unused temporary pages.
  while (!free_pages.empty()) {
    T* page = free_pages.top();
    free_pages.pop();
    if (page < tab || tab + length <= page) {
      delete[] page;
      ++released_pages;
    }
  }

  if (allocated_pages > 4) {
    fprintf(stderr, "ERROR\n");
    std::exit(EXIT_FAILURE);
  }
}

template<typename T>
void merge(T *tab, long n1, long n2, long *gap, long pagesize, long max_thrd) {
  if (n1 <= 0 || n2 <= 0 || pagesize <= 0) {
    fprintf(stderr, "\nparallel merge: incorrect parameters\n");
    std::exit(EXIT_FAILURE);
  }

  long length = n1 + n2;
  long n_pages = (length + pagesize - 1) / pagesize;
  T **pagedir = new T*[n_pages];  // pages with output

  //----------------------------------------------------------------------------
  // STEP 1: compute the initial parameters for each thread. For now, we do it
  //         seqentially.
  //
  // Each thread gets:
  //   * pages_manager -- pointer to free pages manager
  //   * long left_idx, right_idx -- indices to first elems from aubarrays
  //   * initial_bckt_size -- how many element from right seq goes first
  //   * res_size --  number of elements to process
  //   * res_beg -- index to the first elements of the output
  //
  // In short, if we did it sequentially, the thread would just produce the
  // elements of the output in the range [res_beg .. rea_beg + res_size).
  //----------------------------------------------------------------------------

  long pages_per_thread = (n_pages + max_thrd - 1) / max_thrd;
  long n_threads = (n_pages + pages_per_thread - 1) / pages_per_thread;

  // Compute initial parameters for each thread.
  long *left_idx = new long[n_threads];
  long *right_idx = new long[n_threads];
  long *remaining_gap = new long[n_threads];

  long *allocated_pages = new long[n_threads];
  long *released_pages = new long[n_threads];

  for (long i = 0; i < n_threads; ++i) {
    long res_beg = i * pages_per_thread * pagesize;

    // Position of j-th element from left subsequence in the output
    // is j + gap[0] + gap[1] + .. + gap[j]. Thus 'left_idx' is the
    // smallest element in the left subsequence such that its pos in
    // the output is >= res_beg.
    long j = 0, j_pos = j + gap[j];
    while (j_pos < res_beg) { ++j; j_pos += gap[j] + 1; }

    // In parallel it will be sufficient to only look at
    // gap, to compute j and j_pos.
    left_idx[i] = j;
    right_idx[i] = n1 + res_beg - j;

    remaining_gap[i] = j_pos - res_beg;
    allocated_pages[i] = released_pages[i] = 0L;

  }
  
  // Okay, we can start the threads.
  std::thread **threads = new std::thread*[n_threads];
  for (long i = 0; i < n_threads; ++i) {
    long res_beg = i * pages_per_thread * pagesize;
    long res_end = std::min(res_beg + pages_per_thread * pagesize, length);
    long res_size = res_end - res_beg;

    threads[i] = new std::thread(parallel_merge<T>,
      tab, gap, length, pagesize, pagedir, left_idx[i],
      right_idx[i], remaining_gap[i], res_beg, res_size,
      std::ref(allocated_pages[i]), std::ref(released_pages[i]));
  }
  for (long i = 0; i < n_threads; ++i) threads[i]->join();
  for (long i = 0; i < n_threads; ++i) delete threads[i];
  delete[] threads;
  delete[] left_idx;
  delete[] right_idx;
  delete[] remaining_gap;

  long total_allocated_pages = 0L, total_released_pages = 0L;
  for (long i = 0; i < n_threads; ++i) total_allocated_pages += allocated_pages[i];
  for (long i = 0; i < n_threads; ++i) total_released_pages += released_pages[i];
  delete[] allocated_pages;
  delete[] released_pages;

  // If the last input page was incomplete, handle
  // it separatelly and exclude from the computation.
  if (length % pagesize) {
    long size = length % pagesize;
    T *lastpage = pagedir[n_pages - 1];
    T *dest = tab + pagesize * (n_pages - 1);
    std::copy(lastpage, lastpage + size, dest);

    // Release the lastpage if it was temporary.
    if (lastpage < tab || tab + length <= lastpage) {
      delete[] lastpage;
      ++total_released_pages;
    }

    length -= size;
    --n_pages;
  }

  // Find unused input pages.
  bool *usedpage = new bool[n_pages];
  std::fill(usedpage, usedpage + n_pages, false);
  std::vector<std::pair<long, T*> > temp_pages;
  for (long i = 0; i < n_pages; ++i) {
    T *page = pagedir[i];
    if (page < tab || tab + length <= page)
      temp_pages.push_back(std::make_pair(i, page));
    else usedpage[(page - tab) / pagesize] = true;
  }

  // Assign temp pages to unused pages in any
  // order and release the memory for temp pages.
  for (long i = 0, ptr = 0; i < n_pages; ++i) {
    if (!usedpage[i]) {
      long page_id = temp_pages[ptr].first;
      T *src = temp_pages[ptr++].second;
      T *dest = tab + i * pagesize;
      std::copy(src, src + pagesize, dest);
      pagedir[page_id] = dest;
      delete[] src;
      ++total_released_pages;
    }
  }
  delete[] usedpage;

  if (total_allocated_pages != total_released_pages) {
    fprintf(stderr, "\nError\n");
    fprintf(stderr, "\ttotal_allocated_pages = %ld\n", total_allocated_pages);
    fprintf(stderr, "\ttotal_released_pages = %ld\n", total_released_pages);
    std::exit(EXIT_FAILURE);
  }

  // Now all memory should be correctly released.
  // What's left is to check is to permute the pages
  // in parallel. For now, we do it sequentially.
  T *output = new T[length];
  T *dest = output;
  for (long i = 0; i < n_pages; ++i, dest += pagesize)
    std::copy(pagedir[i], pagedir[i] + pagesize, dest);
  std::copy(output, output + length, tab);
  delete[] output;

  delete[] pagedir;
}


#endif  // __PARALLEL_MERGE_H_INCLUDED
