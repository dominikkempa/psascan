#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <string>
#include <unistd.h>

#include "parallel_utils.hpp"
#include "utils.hpp"


using namespace psascan_private;

void test(
    std::uint64_t *tab,
    long length) {

  unsigned char *slab = new unsigned char[sizeof(long) * length];
  parallel_utils::convert_array_to_vbyte_slab(tab, length, slab);

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
void test_random(
    long testcases,
    long max_length,
    long max_value) {

  fprintf(stderr, "TEST, testcases = %ld, "
      "max_n = %ld, max_value = %ld\r",
      testcases, max_length, max_value);

  std::uint64_t * const tab = new std::uint64_t[max_length];

  for (long tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {
    // Print progress information.
    if (dbg == 10) {
      fprintf(stderr, "TEST, testcases = %ld, max_n = %ld, max_value = %ld: "
          "%ld (%.0Lf%%)\r", testcases, max_length, max_value, tc,
          (tc * 100.L) / testcases);
      dbg = 0;
    }

    // Generate string.
    long length = utils::random_int64(1L, max_length);

    for (long j = 0; j < length; ++j)
      tab[j] = utils::random_int64(0, max_value);

    test(tab, length);
  }

  // Clean up.
  delete[] tab;
  
  fprintf(stderr, "TEST, testcases = %ld, max_n = %ld, max_value = %ld: "
      "\033[22;32mPASSED\033[0m%10s\n", testcases, max_length, max_value, "");
}

int main(int, char **) {
  std::srand(std::time(0) + getpid());
  
  // Run tests
#ifdef NDEBUG
  test_random(100000, 10,            1000);
  test_random(100000, 10,         1000000);
  test_random(100000, 10,      1000000000);
  test_random(100000, 10,  1000000000000L);
  test_random(10000, 100,            1000);
  test_random(10000, 100,         1000000);
  test_random(10000, 100,      1000000000);
  test_random(10000, 100,  1000000000000L);
  test_random(3000, 1000,            1000);
  test_random(3000, 1000,         1000000);
  test_random(3000, 1000,      1000000000);
  test_random(3000, 1000,  1000000000000L);
  test_random(300, 10000,            1000);
  test_random(300, 10000,         1000000);
  test_random(300, 10000,      1000000000);
  test_random(300, 10000,  1000000000000L);
#else
  test_random(2000, 10,            1000);
  test_random(2000, 10,         1000000);
  test_random(2000, 10,      1000000000);
  test_random(2000, 10,  1000000000000L);
  test_random(200, 100,            1000);
  test_random(200, 100,         1000000);
  test_random(200, 100,      1000000000);
  test_random(200, 100,  1000000000000L);
  test_random(70, 1000,            1000);
  test_random(70, 1000,         1000000);
  test_random(70, 1000,      1000000000);
  test_random(70, 1000,  1000000000000L);
  test_random(5, 10000,            1000);
  test_random(5, 10000,         1000000);
  test_random(5, 10000,      1000000000);
  test_random(5, 10000,  1000000000000L);
#endif

  fprintf(stderr, "All tests passed.\n");
}

