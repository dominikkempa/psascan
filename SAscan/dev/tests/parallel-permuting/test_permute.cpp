#include <cstdio>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <thread>
#include <unistd.h>

#include "utils.h"
#include "parallel_permute.h"

template<typename T>
void test(T *tab, T **index, long length, long n_threads) {
  T *tab_copy = new T[length];
  std::copy(tab, tab + length, tab_copy);
  T **index_copy = new T*[length];
  std::copy(index, index + length, index_copy);

  T *correct = new T[length];
  for (long i = 0; i < length; ++i)
    correct[i] = *index[i];

  permute<T>(tab, index, length, n_threads);

  if (!std::equal(tab, tab + length, correct)) {
    fprintf(stderr, "Error: ");
    fprintf(stderr, "\tlength = %ld\n", length);
    fprintf(stderr, "\ttab: ");
    for (long i = 0; i < length; ++i)
      fprintf(stderr, "%ld ", (long)tab_copy[i]);
    fprintf(stderr, "\n");
    fprintf(stderr, "\tindex: ");
    for (long i = 0; i < length; ++i)
      fprintf(stderr, "%ld ", (long)(index_copy[i] - tab));
    fprintf(stderr, "\n");
    fprintf(stderr, "\tcorrect: ");
    for (long i = 0; i < length; ++i)
      fprintf(stderr, "%ld ", (long)correct[i]);
    fprintf(stderr, "\n");
    fprintf(stderr, "\tcomputed: ");
    for (long i = 0; i < length; ++i)
      fprintf(stderr, "%ld ", (long)tab[i]);
    fprintf(stderr, "\n");
    std::exit(EXIT_FAILURE);
  }

  delete[] correct;
  delete[] tab_copy;
  delete[] index_copy;
}

void test_random(int testcases, long max_length) {
  fprintf(stderr,"TEST, testcases = %d, max_n = %ld\n", testcases, max_length);

  int *tab = new int[max_length];
  int **index = new int*[max_length];
  int *perm = new int[max_length];

  for (int tc = 0; tc < testcases; ++tc) {
    // Print progress information.
    if (tc % 10 == 0)
    fprintf(stderr,"%d (%.2Lf%%)\r", tc, (tc * 100.L) / testcases);

    // Generate input.
    long length = utils::random_long(1L, max_length);

    for (long i = 0; i < length; ++i)
      perm[i] = tab[i] = i;
    std::random_shuffle(perm, perm + length);
    for (long i = 0; i < length; ++i)
      index[i] = tab + perm[i];

    long n_threads = utils::random_long(1, 50);

    // Run the test on generated input.
    test<int>(tab, index, length, n_threads);
  }

  // Clean up.
  delete[] tab;
  delete[] index;
  delete[] perm;
}

int main() {
  std::srand(std::time(0) + getpid());

  test_random(100000,   10);
  test_random(100000,   100);
  test_random(10000,    1000);
  test_random(10000,    10000);
  test_random(1000,     100000);
  test_random(100,      1000000);
  fprintf(stderr,"All tests passed.\n");
}
