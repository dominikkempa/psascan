#include <ctime>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <string>
#include <algorithm>

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "divsufsort64.h"

#include "../uint40.h"
#include "../utils.h"
#include "../sascan.h"

// Test many string chosen according to given paranters.
void test_random(long testcases, long max_length, long max_sigma) {
  printf("TEST, testcases = %ld, max_n = %ld, max_sigma = %ld\r",
      testcases, max_length, max_sigma);
  std::fflush(stdout);
  unsigned char *text = new unsigned char[max_length + 1];
  long *SA = new long[max_length];

  for (long tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {
    // Print progress information.
    if (dbg == 10) {
      printf("TEST, testcases = %ld, max_n = %ld, max_sigma = %ld: "
          "%ld (%.0Lf%%)\r", testcases, max_length, max_sigma,
          tc, (tc * 100.L) / testcases);
      std::fflush(stdout);
      dbg = 0;
    }

    // Generate string.
    long length = utils::random_long(1, max_length);
    long sigma = utils::random_long(1, max_sigma);
    long ram_use = 0;
    do ram_use = utils::random_long(5L, 5L * length);
    while (4L * ram_use < length);
    if (max_sigma <= 26) utils::fill_random_letters(text, length, sigma);
    else utils::fill_random_string(text, length, sigma);

    // debug //
    // long length = strlen("aaaaaaaaa");
    // strcat((char *)text, "aaaaaaaaa");
    // long ram_use = 39;
    ///////////

    text[length] = 0;
    std::string filename = "/tmp/in" + utils::random_string_hash();
    utils::write_objects_to_file<unsigned char>(text, length, filename);

    // debug //
    // printf("text = %s, length = %ld, ram_use = %ld\n", text, length, ram_use);
    // fflush(stdout);
    ///////////

    // Run the test on generated string.
    SAscan(filename, ram_use);
    
    // Compare the result to correct SA.
    divsufsort64(text, SA, length); // recall that SAscan computes the SA of *reversed* text.
    uint40 *computed_SA = new uint40[length];
    utils::read_n_objects_from_file(computed_SA, length, filename + ".sa5");
    utils::file_delete(filename + ".sa5");
    bool eq = true;
    for (long i = 0; i < length; ++i)
      if ((unsigned long)SA[i] != computed_SA[i].ull()) { eq = false; break; }
    if (!eq) {
    printf("\n\033[22;31mFAILED\033[0m\n");
      if (length < 10000) {
        printf("  text = %s\n", text);
        printf("  computed SA: ");
        for (long k = 0; k < length; ++k) printf("%lu ", (unsigned long)computed_SA[k]);
        printf("\n");
        printf("  correct SA:  ");
        for (long k = 0; k < length; ++k) printf("%ld ", SA[k]);
        printf("\n");
      }
      std::fflush(stdout);
      std::exit(EXIT_FAILURE);
    }
    delete[] computed_SA;
    utils::file_delete(filename);
  }

  // Clean up.
  delete[] text;
  delete[] SA;
  
  printf("TEST, testcases = %ld, max_n = %ld, max_sigma = %ld: "
      "\033[22;32mPASSED\033[0m%10s\n", testcases, max_length, max_sigma, "");
  std::fflush(stdout);
}

int main(int, char **) {
  std::srand(std::time(0) + getpid());

  // Redirect stderr to /dev/null
  int redir = open("/dev/null", O_WRONLY);
  dup2(redir, 2);
  close(redir);

  printf("Testing SAscan.\n");
  std::fflush(stdout);
  test_random(500, 10,      5);
  test_random(500, 10,     20);
  test_random(500, 10,    128);
  test_random(500, 10,    254);
  test_random(500, 100,      5);
  test_random(500, 100,     20);
  test_random(500, 100,    128);
  test_random(500, 100,    254);
  test_random(50, 1000,      5);
  test_random(50, 1000,     20);
  test_random(50, 1000,    128);
  test_random(50, 1000,    254);
  test_random(50, 10000,     5);
  test_random(50, 10000,    20);
  test_random(50, 10000,   128);
  test_random(50, 10000,   254);
  test_random(5, 100000,     5);
  test_random(5, 100000,    20);
  test_random(5, 100000,   128);
  test_random(5, 100000,   254);
  std::fflush(stdout);
}

