#include <cstdio>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <thread>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "utils.h"
#include "parallel_merge.h"
#include "pagearray.h"
#include "inmem_gap_array.h"

template<typename T, unsigned pagesize_log>
void test(T *tab, long n1, long n2, inmem_gap_array *gap, long max_threads) {
  static const unsigned pagesize = (1U << pagesize_log);
  static const unsigned pagesize_mask = pagesize - 1;
  long length = n1 + n2;

  T *tab_copy = new T[length];
  std::copy(tab, tab + length, tab_copy);

  // Compute plain gap, to simplify the computation of correct values.
  long *plain_gap = new long[n1 + 1];
  size_t excess_ptr = 0L;
  for (long i = 0; i <= n1; ++i) {
    plain_gap[i] = gap->m_count[i];
    while (excess_ptr < gap->m_excess.size() && gap->m_excess[excess_ptr] == i) {
      plain_gap[i] += 256L;
      ++excess_ptr;
    }
  }

  // First, compute the answer naively.
  T *correct_answer = new T[length];
  long right_ptr = n1, out_ptr = 0;
  for (long j = 0; j < plain_gap[0]; ++j) correct_answer[out_ptr++] = tab[right_ptr++];
  for (long j = 0; j < n1; ++j) {
    correct_answer[out_ptr++] = tab[j];
    for (long k = 0; k < plain_gap[j + 1]; ++k)
      correct_answer[out_ptr++] = tab[right_ptr++];
  }

  pagearray<T, pagesize_log> *l_pagearray = new pagearray<T, pagesize_log>(tab, tab + n1);
  pagearray<T, pagesize_log> *r_pagearray = new pagearray<T, pagesize_log>(tab + n1, tab + n1 + n2);

  l_pagearray->random_shuffle();
  r_pagearray->random_shuffle();
  typedef pagearray<T, pagesize_log> pagearray_type;
  long temp_result = 0;
  pagearray_type *result = parallel_merge(l_pagearray, r_pagearray, gap, max_threads, -1, 0, temp_result);
  delete l_pagearray;
  delete r_pagearray;

  bool eq = true;
  for (long i = 0; i < length; ++i)
    if (correct_answer[i] != result->m_pageindex[i >> pagesize_log][i & pagesize_mask])
      eq = false;

  if (!eq) {
    fprintf(stdout, "\nError:\n");
    fprintf(stdout, "\ttab: ");
    for (long j = 0; j < length; ++j)
      fprintf(stdout, "%ld ", (long)tab_copy[j]);
    fprintf(stdout, "\n");
    fprintf(stdout, "\tn1 = %ld, n2 = %ld\n", n1, n2);
    fprintf(stdout, "\tpagesize = %u\n", pagesize);
    fprintf(stdout, "\tmax_threads = %ld\n", max_threads);
    fprintf(stdout, "\tgap: ");
    for (long j = 0; j <= n1; ++j)
      fprintf(stdout, "%ld ", plain_gap[j]);
    fprintf(stdout, "\n");
    fprintf(stdout, "Correct: ");
    for (long j = 0; j < length; ++j)
      fprintf(stdout, "%ld ", (long)correct_answer[j]);
    fprintf(stdout, "\n");
    fprintf(stdout, "Computed: ");
    for (long j = 0; j < length; ++j)
      fprintf(stdout, "%ld ", (long)result->m_pageindex[j >> pagesize_log][j & pagesize_mask]);
    fprintf(stdout, "\n");
    std::fflush(stdout);
    std::exit(EXIT_FAILURE);
  }

  // Test the correctness of permuting the page array into plain array.
  result->permute_to_plain_array(max_threads);
  delete result;

  if (!std::equal(tab, tab + length, correct_answer)) {
    fprintf(stdout, "\nError:\n");
    fprintf(stdout, "\ttab: ");
    for (long j = 0; j < length; ++j)
      fprintf(stdout, "%ld ", (long)tab_copy[j]);
    fprintf(stdout, "\n");
    fprintf(stdout, "\tn1 = %ld, n2 = %ld\n", n1, n2);
    fprintf(stdout, "\tpagesize = %u\n", pagesize);
    fprintf(stdout, "\tmax_threads = %ld\n", max_threads);
    fprintf(stdout, "\tgap: ");
    for (long j = 0; j <= n1; ++j)
      fprintf(stdout, "%ld ", plain_gap[j]);
    fprintf(stdout, "\n");
    fprintf(stdout, "Correct: ");
    for (long j = 0; j < length; ++j)
      fprintf(stdout, "%ld ", (long)correct_answer[j]);
    fprintf(stdout, "\n");
    fprintf(stdout, "Computed: ");
    for (long j = 0; j < length; ++j)
      fprintf(stdout, "%ld ", (long)tab[j]);
    fprintf(stdout, "\n");
    std::fflush(stdout);
    std::exit(EXIT_FAILURE);
  }


  delete[] plain_gap;
  delete[] correct_answer;
  delete[] tab_copy;
}

template<unsigned pagesize_log>
void test_random(int testcases, long max_length, long maxval) {
  static const unsigned pagesize = (1U << pagesize_log);
  fprintf(stdout,"TEST, pagesize = %4u, testcases = %6d, max_n = %8ld\n", pagesize, testcases, max_length);
  std::fflush(stdout);

  int *tab = new int[max_length * 2 + 1000000];
  for (int tc = 0; tc < testcases; ++tc) {
    // Print progress information.
    fprintf(stdout,"%d (%.2Lf%%)\r", tc, (tc * 100.L) / testcases);
    std::fflush(stdout);

    long max_threads = utils::random_long(1, 50);

    // Generate input.
    long n1 = utils::random_long(1, max_length);
    while (n1 % pagesize) ++n1; // make sure that size of the left array is a multiple of pagesize
    long n2 = utils::random_long(1, max_length);

    for (long j = 0; j < n1 + n2; ++j)
      tab[j] = utils::random_int(0, maxval);

    inmem_gap_array *g = new inmem_gap_array(n1 + 1);
    for (long j = 0; j < n2; ++j) {
      long pos = utils::random_long(0, n1);
      g->m_count[pos] += 1;
      if (g->m_count[pos] == 0)
        g->m_excess.push_back(pos);
    }
    std::sort(g->m_excess.begin(), g->m_excess.end());

    /*long n1 = 3;
    long n2 = 4;
    long pagesize = 6;
    long max_threads = 26;
    tab[0] = 3; tab[1] = 5; tab[2] = 3; tab[3] = 0; tab[4] = 3; tab[5] = 4; tab[6] = 1;
    gap[0] = 2; gap[1] = 0; gap[2] = 1; gap[3] = 1;*/

    // Run the test on generated input.
    test<int, pagesize_log>(tab, n1, n2, g, max_threads);
    delete g;
  }

  // Clean up.
  delete[] tab;
}

int main() {
  static const long maxval = 1000000000L;
  std::srand(std::time(0) + getpid());

  // Redirect stdout to /dev/null
  int redir = open("/dev/null", O_WRONLY);
  dup2(redir, 2);
  close(redir);

  test_random<0> (2000,  10,        maxval);
  test_random<4> (2000,  10,        maxval);
  test_random<7> (2000,  10,        maxval);
  test_random<10>(2000,  10,        maxval);
  test_random<12>(2000,  10,        maxval);

  test_random<0> (2000,  100,       maxval);
  test_random<4> (2000,  100,       maxval);
  test_random<7> (2000,  100,       maxval);
  test_random<10>(2000,  100,       maxval);
  test_random<12>(2000,  100,       maxval);

  test_random<0> (2000,   1000,      maxval);
  test_random<4> (2000,   1000,      maxval);
  test_random<7> (2000,   1000,      maxval);
  test_random<10>(2000,   1000,      maxval);
  test_random<12>(2000,   1000,      maxval);

  test_random<0> (2000,    10000,     maxval);
  test_random<4> (2000,    10000,     maxval);
  test_random<7> (2000,    10000,     maxval);
  test_random<10>(2000,    10000,     maxval);
  test_random<12>(2000,    10000,     maxval);

  test_random<0> (200,     100000,    maxval);
  test_random<4> (200,     100000,    maxval);
  test_random<7> (200,     100000,    maxval);
  test_random<10>(200,     100000,    maxval);
  test_random<12>(200,     100000,    maxval);

  test_random<0> (100,      1000000,   maxval);
  test_random<4> (100,      1000000,   maxval);
  test_random<7> (100,      1000000,   maxval);
  test_random<10>(100,      1000000,   maxval);
  test_random<12>(100,      1000000,   maxval);

  test_random<0> (10,       10000000,  maxval);
  test_random<4> (10,       10000000,  maxval);
  test_random<7> (10,       10000000,  maxval);
  test_random<10>(10,       10000000,  maxval);
  test_random<12>(10,       10000000,  maxval);

  fprintf(stdout,"All tests passed.\n");
}
