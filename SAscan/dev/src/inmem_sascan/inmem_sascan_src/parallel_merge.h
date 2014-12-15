#ifndef __INMEM_SASCAN_PARALLEL_MERGE_H_INCLUDED
#define __INMEM_SASCAN_PARALLEL_MERGE_H_INCLUDED

#include <cstdio>
#include <cstdlib>

#include <vector>
#include <stack>
#include <algorithm>
#include <thread>
#include <mutex>

#include "pagearray.h"
#include "inmem_gap_array.h"

namespace inmem_sascan_private {


//==============================================================================
// Compute the range [res_beg..res_beg+res_size) of the output (i.e., the
// sequence after merging). The range is guaranteed to be aligned with page
// boundaries.
//==============================================================================
template<typename pagearray_type>
void parallel_merge_aux(
    pagearray_type *l_pagearray,
    pagearray_type *r_pagearray,
    pagearray_type *output,
    inmem_gap_array *gap,
    long left_idx, long right_idx,
    long remaining_gap,
    long res_beg, long res_size,
    long what_to_add) {

  typedef typename pagearray_type::value_type value_type;
  static const unsigned pagesize = pagearray_type::pagesize;

  long lpage_read = 0L;
  long lpage_id = l_pagearray->get_page_id(left_idx);
  long lpage_offset = l_pagearray->get_page_offset(left_idx);
  value_type *lpage = l_pagearray->m_pageindex[lpage_id++];

  long rpage_read = 0L;
  long rpage_id = r_pagearray->get_page_id(right_idx);
  long rpage_offset = r_pagearray->get_page_offset(right_idx);
  value_type *rpage = r_pagearray->m_pageindex[rpage_id++];

  long pageid = output->get_page_id(res_beg);
  long filled = output->get_page_offset(res_beg);
  value_type *dest = new value_type[pagesize];
  output->m_pageindex[pageid++] = dest;

  std::stack<value_type*> freepages;

  size_t excess_ptr = std::lower_bound(gap->m_excess.begin(),
      gap->m_excess.end(), left_idx + 1) - gap->m_excess.begin();

  for (long i = 0; i < res_size; ++i) {
    if (filled == pagesize) {
      if (freepages.empty()) dest = new value_type[pagesize];
      else { dest = freepages.top(); freepages.pop(); }
      output->m_pageindex[pageid++] = dest;
      filled = 0L;
    }
    if (remaining_gap > 0) {
      --remaining_gap;
      // The next element comes from the right subarray.
      dest[filled] = rpage[rpage_offset++];
      dest[filled++].sa += what_to_add;
      rpage_read++;
      if (rpage_offset == pagesize) {
        // We reached the end of page in the right subarray.
        // We put it into free pages if we read exactly
        // pagesize elements from it. This means the no other
        // thread will attemp to read from it in the future.
        if (rpage_read == pagesize) freepages.push(r_pagearray->m_pageindex[rpage_id - 1]);

        // Note: we don't have to check, if the page below exists, because we have
        // a sentinel page in the page index of every pagearray.
        rpage = r_pagearray->m_pageindex[rpage_id++];
        rpage_offset = 0L;
        rpage_read = 0L;
      }
    } else {
      // Next elem comes from the left subarray.
      dest[filled++] = lpage[lpage_offset++];
      left_idx++;
      lpage_read++;

      // Compute gap[left_idx].
      long gap_left_idx = gap->m_count[left_idx];
      while (excess_ptr < gap->m_excess.size() &&
          gap->m_excess[excess_ptr] == left_idx) {
        gap_left_idx += 256L;
        ++excess_ptr;
      }

      remaining_gap = gap_left_idx;
      if (lpage_offset == pagesize) {
        // We reached the end of page in the left
        // subarray, proceed analogously.
        if (lpage_read == pagesize) freepages.push(l_pagearray->m_pageindex[lpage_id - 1]);

        // Note: we don't have to check, if the page below exists, because we have
        // a sentinel page in the page index of every pagearray.
        lpage = l_pagearray->m_pageindex[lpage_id++];
        lpage_offset = 0L;
        lpage_read = 0L;
      }
    }
  }

  // Release the unused auxiliary pages.
  while (!freepages.empty()) {
    value_type* p = freepages.top();
    freepages.pop();
    if (!output->owns_page(p))
      delete[] p;
  }
}

template<typename pagearray_type>
pagearray_type *parallel_merge(pagearray_type *l_pagearray, pagearray_type *r_pagearray,
    inmem_gap_array *gap, long max_threads, long i0, long &aux_result, long what_to_add) {
  typedef typename pagearray_type::value_type value_type;
  static const unsigned pagesize_log = pagearray_type::pagesize_log;
  static const unsigned pagesize = pagearray_type::pagesize;

  //----------------------------------------------------------------------------
  // STEP 1: compute the initial parameters for each thread. For now, we do it
  //         seqentially. Each thread gets:
  //
  //  - long left_idx, right_idx -- indices to first elems from aubarrays
  //  - initial_bckt_size -- how many element from right seq goes first
  //  - res_size --  number of elements to process
  //  - res_beg -- index to the first elements of the output
  //
  // In short, if we did it sequentially, the thread would just produce the
  // elements of the output in the range [res_beg .. rea_beg + res_size).
  //----------------------------------------------------------------------------
  fprintf(stderr, "queries: ");
  long double start = utils::wclock();
  long length = l_pagearray->m_length + r_pagearray->m_length;
  long n_pages = (length + pagesize - 1) / pagesize;
  long pages_per_thread = (n_pages + max_threads - 1) / max_threads;
  long n_threads = (n_pages + pages_per_thread - 1) / pages_per_thread;

  // Compute initial parameters for each thread.
  long *left_idx = new long[n_threads];
  long *right_idx = new long[n_threads];
  long *remaining_gap = new long[n_threads];

  // Prepare gap queries.
  long *gap_query = new long[n_threads];
  long *gap_answer_a = new long[n_threads];
  long *gap_answer_b = new long[n_threads];
  for (long i = 0; i < n_threads; ++i)
    gap_query[i] = i * pages_per_thread * pagesize;

  // Answer these queries in parallel and convert the answers
  // to left_idx, right_idx and remaining_gap values.
  aux_result = gap->answer_queries(n_threads, gap_query, gap_answer_a, gap_answer_b, max_threads, i0);
  for (long i = 0; i < n_threads; ++i) {
    long res_beg = i * pages_per_thread * pagesize;
    long j = gap_answer_a[i], s = gap_answer_b[i];
    left_idx[i] = j;
    right_idx[i] = res_beg - j;
    remaining_gap[i] = j + s - res_beg;
  }
  delete[] gap_query;
  delete[] gap_answer_a;
  delete[] gap_answer_b;
  fprintf(stderr, "%.2Lf ", utils::wclock() - start);


  //----------------------------------------------------------------------------
  // STEP 2: merging the arrays.
  //----------------------------------------------------------------------------
  fprintf(stderr, "merge: ");
  start = utils::wclock();
  typedef pagearray<value_type, pagesize_log> output_type;
  output_type *result = new output_type(l_pagearray->m_origin, length);

  std::thread **threads = new std::thread*[n_threads];
  for (long t = 0; t < n_threads; ++t) {
    long res_beg = t * pages_per_thread * pagesize;
    long res_end = std::min(res_beg + pages_per_thread * pagesize, length);
    long res_size = res_end - res_beg;

    threads[t] = new std::thread(parallel_merge_aux<pagearray_type>,
        l_pagearray, r_pagearray, result, gap,  left_idx[t], right_idx[t],
        remaining_gap[t], res_beg, res_size, what_to_add);
  }
  for (long t = 0; t < n_threads; ++t) threads[t]->join();
  for (long t = 0; t < n_threads; ++t) delete threads[t];
  delete[] threads;
  delete[] left_idx;
  delete[] right_idx;
  delete[] remaining_gap;

  // If the last input page was incomplete, handle
  // it separatelly and exclude from the computation.
  if (length % pagesize) {
    long size = length % pagesize;
    value_type *lastpage = result->m_pageindex[n_pages - 1];
    value_type *dest = result->m_origin + pagesize * (n_pages - 1);
    std::copy(lastpage, lastpage + size, dest);
    result->m_pageindex[n_pages - 1] = dest;

    // Release the lastpage if it was temporary.
    if (!result->owns_page(lastpage))
      delete[] lastpage;

    length -= size;
    --n_pages;
  }

  // Find unused input pages.
  bool *usedpage = new bool[n_pages];
  std::fill(usedpage, usedpage + n_pages, false);
  std::vector<std::pair<long, value_type*> > auxpages;
  for (long i = 0; i < n_pages; ++i) {
    value_type *p = result->m_pageindex[i];
    if (result->owns_page(p)) usedpage[result->get_page_id(p)] = true;
    else auxpages.push_back(std::make_pair(i, p));
  }

  // Assign aux pages to unused pages in any
  // order and release them (aux pages).
  for (long i = 0, ptr = 0; i < n_pages; ++i) {
    if (!usedpage[i]) {
      long page_id = auxpages[ptr].first;
      value_type *src = auxpages[ptr++].second;
      value_type *dest = result->m_origin + i * pagesize;
      std::copy(src, src + pagesize, dest);
      result->m_pageindex[page_id] = dest;
      delete[] src;
    }
  }
  delete[] usedpage;
  fprintf(stderr, "%.2Lf ", utils::wclock() - start);

  return result;
}

}  // namespace inmem_sascan


#endif  // __PARALLEL_MERGE_H_INCLUDED
