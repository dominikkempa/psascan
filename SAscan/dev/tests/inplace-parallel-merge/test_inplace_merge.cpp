#include <cstdio>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <thread>
#include <unistd.h>

#include "utils.h"
#include "parallel_merge.h"

template<typename T>
void test(T *tab, long n1, long n2, long *gap, unsigned pagesize_bits, long max_threads) {
  static const unsigned pagesize = (1U << pagesize_bits);
  long length = n1 + n2;

  T *tab_copy = new T[length];
  std::copy(tab, tab + length, tab_copy);

  long *gap_copy = new long[n1 + 1];
  std::copy(gap, gap + n1 + 1, gap_copy);

  // First, compute the answer naively.
  T *correct_answer = new T[length];
  long right_ptr = n1, out_ptr = 0;
  for (long j = 0; j < gap[0]; ++j) correct_answer[out_ptr++] = tab[right_ptr++];
  for (long j = 0; j < n1; ++j) {
    correct_answer[out_ptr++] = tab[j];
    for (long k = 0; k < gap[j + 1]; ++k)
      correct_answer[out_ptr++] = tab[right_ptr++];
  }

  // Merge elements using parallel in-place procedure.
  merge<T>(tab, n1, n2, gap, pagesize_bits,/*pagesize,*/ max_threads);
  if (!std::equal(tab, tab + length, correct_answer)) {
    fprintf(stderr, "\nError:\n");
    fprintf(stderr, "\ttab: ");
    for (long j = 0; j < length; ++j)
      fprintf(stderr, "%ld ", (long)tab_copy[j]);
    fprintf(stderr, "\n");
    fprintf(stderr, "\tn1 = %ld, n2 = %ld\n", n1, n2);
    fprintf(stderr, "\tpagesize = %u\n", pagesize);
    fprintf(stderr, "\tmax_threads = %ld\n", max_threads);
    fprintf(stderr, "\tgap: ");
    for (long j = 0; j <= n1; ++j)
      fprintf(stderr, "%ld ", gap_copy[j]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Correct: ");
    for (long j = 0; j < length; ++j)
      fprintf(stderr, "%ld ", (long)correct_answer[j]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Computed: ");
    for (long j = 0; j < length; ++j)
      fprintf(stderr, "%ld ", (long)tab[j]);
    fprintf(stderr, "\n");
    std::exit(EXIT_FAILURE);
  }

  delete[] correct_answer;
  delete[] tab_copy;
  delete[] gap_copy;
}

void test_random(int testcases, long max_length, long maxval) {
  fprintf(stderr,"TEST, testcases = %d, max_n = %ld\n", testcases, max_length);

  int *tab = new int[max_length * 2];
  long *gap = new long[max_length + 1];

  for (int tc = 0; tc < testcases; ++tc) {
    // Print progress information.
    if (tc % 1000 == 0)
    fprintf(stderr,"%d (%.2Lf%%)\r", tc, (tc * 100.L) / testcases);

    // Generate input.
    long n1 = utils::random_long(1, max_length);
    long n2 = utils::random_long(1, max_length);

    for (long j = 0; j < n1 + n2; ++j)
      tab[j] = utils::random_int(0, maxval);

    std::fill(gap, gap + n1 + 1, 0L);
    for (long j = 0; j < n2; ++j)
      ++gap[utils::random_long(0, n1)];

    long max_threads = utils::random_long(1, 50);
    unsigned pagesize_bits = utils::random_long(0U, 7U);

    /*long n1 = 3;
    long n2 = 4;
    long pagesize = 6;
    long max_threads = 26;
    tab[0] = 3; tab[1] = 5; tab[2] = 3; tab[3] = 0; tab[4] = 3; tab[5] = 4; tab[6] = 1;
    gap[0] = 2; gap[1] = 0; gap[2] = 1; gap[3] = 1;*/

    // Run the test on generated input.
    test<int>(tab, n1, n2, gap, pagesize_bits, max_threads);
  }

  // Clean up.
  delete[] tab;
  delete[] gap;
}

int main() {
  static const long maxval = 1000000000L;
  std::srand(std::time(0) + getpid());

  test_random(10000000, 10,        maxval);
  test_random(1000000,  100,       maxval);
  test_random(100000,   1000,      maxval);
  test_random(10000,    10000,     maxval);
  test_random(1000,     100000,    maxval);
  test_random(100,      1000000,   maxval);
  test_random(10,       10000000,  maxval);
  fprintf(stderr,"All tests passed.\n");
}
