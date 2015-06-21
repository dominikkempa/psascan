#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <string>
#include <unistd.h>

#include "parallel_utils.h"
#include "utils.h"

void test(long *tab, long length, long max_threads) {
  unsigned char *slab = new unsigned char[sizeof(long) * length];
  parallel_utils::convert_array_to_vbyte_slab(tab, length, slab, max_threads);

  long *tab2 = new long[length];
  long ptr = 0;
  for (long t = 0; t < length; ++t) {
    tab2[t] = 0;
    long offset = 0;
    while (slab[ptr] & 0x80) {
      tab2[t] |= ((((long)slab[ptr++]) & 0x7f) << offset);
      offset += 7;
    }
    tab2[t] |= (((long)slab[ptr++]) << offset);
  }

  if (!std::equal(tab, tab + length, tab2)) {
    fprintf(stdout, "\n\033[22;31mFAILED\033[0m\n");
    fprintf(stderr, "\ttab: ");
    for (long j = 0; j < length; ++j)
      fprintf(stderr, "%ld ", tab[j]);
    fprintf(stderr, "\n");
    fprintf(stderr, "\ttab2: ");
    for (long j = 0; j < length; ++j)
      fprintf(stderr, "%ld ", tab2[j]);
    fprintf(stderr, "\n");

    std::exit(EXIT_FAILURE);
  }

  delete[] slab;
  delete[] tab2;
}

// Test many string chosen according to given paranters.
void test_random(long testcases, long max_length, long max_value) {
  fprintf(stderr, "TEST, testcases = %ld, max_n = %ld, max_value = %ld\r",
      testcases, max_length, max_value);

  long *tab = new long[max_length];

  for (long tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {
    // Print progress information.
    if (dbg == 1000) {
      fprintf(stderr, "TEST, testcases = %ld, max_n = %ld, max_value = %ld: "
          "%ld (%.0Lf%%)\r", testcases, max_length, max_value, tc,
          (tc * 100.L) / testcases);
      dbg = 0;
    }

    // Generate string.
    long length = utils::random_long(1L, max_length);
    long max_threads = utils::random_long(1L, 50L);

    for (long j = 0; j < length; ++j)
      tab[j] = utils::random_long(0, max_value);

    test(tab, length, max_threads);
  }

  // Clean up.
  delete[] tab;
  
  fprintf(stderr, "TEST, testcases = %ld, max_n = %ld, max_value = %ld: "
      "\033[22;32mPASSED\033[0m%10s\n", testcases, max_length, max_value, "");
}

int main(int, char **) {
  std::srand(std::time(0) + getpid());
  
  // Run tests
  test_random(500000, 10,            1000);
  test_random(500000, 10,         1000000);
  test_random(500000, 10,      1000000000);
  test_random(500000, 10,  1000000000000L);
  test_random(50000, 100,            1000);
  test_random(50000, 100,         1000000);
  test_random(50000, 100,      1000000000);
  test_random(50000, 100,  1000000000000L);
  test_random(5000, 1000,            1000);
  test_random(5000, 1000,         1000000);
  test_random(5000, 1000,      1000000000);
  test_random(5000, 1000,  1000000000000L);
  test_random(500, 10000,            1000);
  test_random(500, 10000,         1000000);
  test_random(500, 10000,      1000000000);
  test_random(500, 10000,  1000000000000L);

  fprintf(stderr, "All tests passed.\n");
}

