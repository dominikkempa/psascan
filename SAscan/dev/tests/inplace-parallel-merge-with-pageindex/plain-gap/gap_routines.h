#ifndef __GAP_ROUTINES_H_INCLUDED
#define __GAP_ROUTINES_H_INCLUDED

#include <cstdio>
#include <cstdlib>

#include <vector>
#include <stack>
#include <algorithm>
#include <thread>
#include <mutex>


//==============================================================================
// Find and smallest j such that j + gap[0] + .. + gap[j] >= a. Store
// the value of j into b and gap[0] + .. + gap[j] into c. To speed up the
// algorithm, we have array gapsum defined as
//
//    gapsum[i] = gap[0] + .. + gap[i * block_size - 1].
//
//==============================================================================
void answer_single_gap_query(int *gap, long length, long block_size,
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
  while (true) {
    // Invariant: sum = gap[0] + .. + gap[j - 1].
    long gap_j = gap[j]; // compute gap[j]

    if (j + sum + gap_j >= a) { b = j; c = sum + gap_j; return; }
    else { sum += gap_j; ++j; }
  }
}


//==============================================================================
// Compute gap[beg] + .. + gap[end - 1] and store into result.
//==============================================================================
void compute_sum(int *gap, long beg, long end, long &result) {
  result = 0;
  for (long i = beg; i < end; ++i)
    result += gap[i];
}


//==============================================================================
// Parallel computaton of answers to n_queries queries of the form:
// What is the smallest j such that j + gap[0] + .. + gap[j] >= a[i]"
//   - the answer to i-th query is stored in b[i]
//   - in addition we also return gap[0] + .. + gap[j] in c[i]
//
// To do that we first split the gap array into blocks of size of about
// length / max_threads and (in parallel) compute sums of gap values inside
// these blocks. We the accumulate these sums into array of prefix sums.
//
// To answer each of the queries we start a separate thread. Each thread uses
// the partial sums of gap array at block boundaries to find a good starting
// point for search and then scans the gap array from there.
//==============================================================================
void answer_gap_queries(int *gap, long length, long n_queries,
    long *a, long *b, long *c, long max_threads) {
  //----------------------------------------------------------------------------
  // STEP 1: split gap array into at most max_threads blocks
  // and in parallel compute sum of values inside each block.
  //----------------------------------------------------------------------------
  fprintf(stderr, "Precomp1: ");
  long double start = utils::wclock();
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
  fprintf(stderr, "Precomp3: ");
  start = utils::wclock();
  threads = new std::thread*[n_queries];
  for (long i = 0; i < n_queries; ++i)
    threads[i] = new std::thread(answer_single_gap_query, gap, length,
      block_size, gapsum, a[i], std::ref(b[i]), std::ref(c[i]));
  for (long i = 0; i < n_queries; ++i) threads[i]->join();
  for (long i = 0; i < n_queries; ++i) delete threads[i];
  delete[] threads;
  delete[] gapsum;
  fprintf(stderr, "%5.2Lf ", utils::wclock() - start);
}


#endif  // __GAP_ROUTINES_H_INCLUDED
