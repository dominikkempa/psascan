#include <cstdio>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <thread>
#include <unistd.h>

#include "uint40.h"
#include "utils.h"
#include "parallel_shrink.h"


template<typename T, typename S>
void test(T *tab, long length, long max_threads) {
  S *correct = new S[length];
  for (long i = 0; i < length; ++i)
    correct[i] = (S)tab[i];

  T *orig_tab = new T[length];
  std::copy(tab, tab + length, orig_tab);

  S *computed = parallel_shrink<T, S>(tab, length, max_threads);

  if (!std::equal(correct, correct + length, computed)) {
    fprintf(stderr, "Error:\n");
    fprintf(stderr, "\tlength = %ld\n", length);
    fprintf(stderr, "\tsizeof(T) = %ld\n", (long)sizeof(T));
    fprintf(stderr, "\tsizeof(S) = %ld\n", (long)sizeof(T));
    if (length < 1000) {
      fprintf(stderr, "\ttab: ");
      for (long j = 0; j < length; ++j)
        fprintf(stderr, "%ld ", (long)orig_tab[j]);
      fprintf(stderr, "\n");
      fprintf(stderr, "correct: ");
      for (long j = 0; j < length; ++j)
        fprintf(stderr, "%ld ", (long)correct[j]);
      fprintf(stderr, "\n");
      fprintf(stderr, "computed: ");
      for (long j = 0; j < length; ++j)
        fprintf(stderr, "%ld ", (long)computed[j]);
      fprintf(stderr, "\n");
    }
    std::exit(EXIT_FAILURE);
  }

  delete[] correct;
  delete[] orig_tab;
}


void test_random(int testcases, long max_length) {
  fprintf(stderr,"TEST, testcases = %d, max_n = %ld\n", testcases, max_length);

  long *tab = new long[max_length];

  for (int tc = 0; tc < testcases; ++tc) {
    // Print progress information.
    if (tc % 10 == 0)
      fprintf(stderr,"%d (%.2Lf%%)\r", tc, (tc * 100.L) / testcases);

    // Generate input.
    long length = utils::random_long(1L, max_length);

    for (long i = 0; i < length; ++i)
      tab[i] = utils::random_long(1L, 137438953472L);

    long max_threads = utils::random_long(1, 50);

    // Run the test on generated input.
    test<long, uint40>(tab, length, max_threads);
  }

  // Clean up.
  delete[] tab;
}

int main() {
  std::srand(std::time(0) + getpid());

  test_random(100000,  10);
  test_random(10000,   100);
  test_random(1000,    1000);
  test_random(1000,    10000);
  test_random(100,     100000);
  test_random(10,      1000000);
  test_random(10,       10000000);
  test_random(1,        100000000);
  fprintf(stderr,"All tests passed.\n");
}

