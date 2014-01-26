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

#include "uint40.h"
#include "divsufsort64.h"
#include "utils.h"
#include "sascan.h"

// Test many string chosen according to given paranters.
void test_random(long testcases, long max_length, long max_sigma) {
  printf("TEST, testcases = %ld, max_n = %ld, max_sigma = %ld\n",
      testcases, max_length, max_sigma);
  fflush(stdout);
  unsigned char *text = new unsigned char[max_length + 1];
  long *SA = new long[max_length];

  for (long tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {
    // Print progress information.
    if (dbg == 10) {
      printf("%ld (%.2Lf%%)\r", tc, (tc * 100.L) / testcases);
      fflush(stdout);
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
    utils::write_text_to_file(text, length, filename);

    // debug //
    // printf("text = %s, length = %ld, ram_use = %ld\n", text, length, ram_use);
    // fflush(stdout);
    ///////////

    // Run the test on generated string.
    SAscan(filename, ram_use);
    
    // Compare the result to correct SA.
    std::reverse(text, text + length);
    divsufsort64(text, SA, length); // recall that SAscan computes the SA of *reversed* text.
    std::reverse(text, text + length);
    uint40 *computed_SA;
    utils::read_n_objects_from_file(computed_SA, length, filename + ".sa5");
    utils::file_delete(filename + ".sa5");
    bool eq = true;
    for (long i = 0; i < length; ++i)
      if ((unsigned long)SA[i] != computed_SA[i].ull()) { eq = false; break; }
    if (!eq) {
      printf("Error!\n");
      if (length < 10000) {
        printf("  text = %s\n", text);
        printf("  computed SA: ");
        for (long k = 0; k < length; ++k) printf("%lu ", (unsigned long)computed_SA[k]);
        printf("\n");
        printf("  correct SA:  ");
        for (long k = 0; k < length; ++k) printf("%ld ", SA[k]);
        printf("\n");
      }
      fflush(stdout);
      std::exit(EXIT_FAILURE);
    }
    delete[] computed_SA;
    utils::file_delete(filename);
  }

  // Clean up.
  delete[] text;
  delete[] SA;
}

int main(int, char **) {
  srand(time(0) + getpid());

  // Redirect stderr to /dev/null
  int redir = open("/dev/null", O_WRONLY);
  dup2(redir, 2);
  close(redir);

  printf("Testing SAscan.\n");
  fflush(stdout);
  test_random(5000, 10,      5);
  test_random(5000, 10,     20);
  test_random(5000, 10,    128);
  test_random(5000, 10,    254);
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
  printf("All tests passed.\n");
  fflush(stdout);

}

