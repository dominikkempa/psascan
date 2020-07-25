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
  fprintf(stderr, "\nmax_threads = %ld\n", max_threads);
  S *correct = (S *)malloc(length * sizeof(S));
  fprintf(stderr, "Computing correct output: ");
  long double start = utils::wclock();
  for (long i = 0; i < length; ++i)
    correct[i] = (S)tab[i];
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

  T *orig_tab = (T *)malloc(length * sizeof(T));
  std::copy(tab, tab + length, orig_tab);

  fprintf(stderr, "Running the algorithm: ");
  start = utils::wclock();
  S *computed = parallel_shrink<T, S>(tab, length, max_threads);
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

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

  free(correct);
  free(orig_tab);
}


void test_random(long length) {
  fprintf(stderr, "Length = %ld\n", length);

  long *tab = (long *)malloc(length * sizeof(long));

  fprintf(stderr, "Generating input: ");
  long double start = utils::wclock();
  long initial_chunk = std::min(1L << 20, length);
  for (long i = 0; i < initial_chunk; ++i)
    tab[i] = utils::random_long(1L, 137438953472L);
  for (long i = initial_chunk; i < length; ++i)
    tab[i] = tab[i - initial_chunk];
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

  // Run the test on generated input.
  test<long, uint40>(tab, length, 24);
  test<long, uint40>(tab, length, 12);
  test<long, uint40>(tab, length, 6);
  test<long, uint40>(tab, length, 4);
  test<long, uint40>(tab, length, 1);

  // Clean up.
  free(tab);
}

int main() {
  std::srand(std::time(0) + getpid());

  test_random(5L << 30);
}

