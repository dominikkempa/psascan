#include <cstdio>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <thread>
#include <unistd.h>

#include "utils.h"
#include "gap_queries.h"

void test(int *gap, long length, long n_queries, long *queries, long max_threads) {
  // First, compute the answer naively.
  long *correct_b = new long[n_queries];
  long *correct_c = new long[n_queries];
  for (long i = 0; i < n_queries; ++i) {
    long x = queries[i];
    
    long j = 0, sum = 0;
    // invariant: sum = gap[0] + .. + gap[j - 1]
    while (j + sum + gap[j] < x) sum += gap[j++];
    
    correct_b[i] = j;
    correct_c[i] = sum + gap[j];
  }

  // Second, compute answers using parallel algorithm.
  long *computed_b = new long[n_queries];
  long *computed_c = new long[n_queries];
  answer_gap_queries(gap, length, n_queries, queries, computed_b, computed_c, max_threads);
  
  // Compare.
  for (long i = 0; i < n_queries; ++i) {
    if (correct_b[i] != computed_b[i] || correct_c[i] != computed_c[i]) {
      fprintf(stderr, "Error\n");
      fprintf(stderr, "\tlength = %ld\n", length);
      fprintf(stderr, "\tgap: ");
      for (long j = 0; j < length; ++j)
        fprintf(stderr, "%d ", gap[j]);
      fprintf(stderr, "\n");
      fprintf(stderr, "\tquery = %ld\n", queries[i]);
      fprintf(stderr, "\tcorrect answers: b = %ld, c = %ld\n",
        correct_b[i], correct_c[i]);
      fprintf(stderr, "\tcomputed answers: b = %ld, c = %ld\n",
        computed_b[i], computed_c[i]);
      std::exit(EXIT_FAILURE);
    }
  }

  delete[] computed_b;
  delete[] computed_c;
  delete[] correct_b;
  delete[] correct_c;
}

void test_random(int testcases, long max_length, long max_queries) {
  fprintf(stderr,"TEST, testcases = %d, max_n = %ld, max_queries = %ld\n",
    testcases, max_length, max_queries);

  int *gap = new int[max_length];
  long *queries = new long[max_queries];

  for (int tc = 0; tc < testcases; ++tc) {
    // Print progress information.
    if (tc % 100 == 0)
      fprintf(stderr,"%d (%.2Lf%%)\r", tc, (tc * 100.L) / testcases);

    // Generate input.
    static const int INF = 1000000;
    long length = utils::random_long(1, max_length);
    long gapsum = 0;
    for (long i = 0; i < length; ++i) {
      gap[i] = utils::random_int(0, INF);
      gapsum += gap[i];
    }

    long max_threads = utils::random_long(1, 50);
    long n_queries = utils::random_long(1, max_queries);
    for (long i = 0; i < n_queries; ++i)
      queries[i] = utils::random_long(0, gapsum + length - 1);

    // Run the test.
    test(gap, length, n_queries, queries, max_threads);
  }

  // Clean up.
  delete[] gap;
  delete[] queries;
}

int main() {
  std::srand(std::time(0) + getpid());

  test_random(100000,  10,       50);
  test_random(100000,  100,      50);
  test_random(10000,   1000,     50);
  test_random(10000,   10000,    50);
  test_random(10000,   100000,   50);
  test_random(1000,    1000000,  50);
  test_random(100,     10000000, 50);
  fprintf(stderr,"All tests passed.\n");
}
