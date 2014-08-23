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
//#include "utils.h"
#include "../sascan.h"

extern long stream_buffer_size;

int random_int(int p, int r) {
  return p + rand() % (r - p + 1);
}

long random_long(long p, long r) {
  long x = random_int(0, 1000000000);
  long y = random_int(0, 1000000000);
  long z = x * 1000000000L + y;
  return p + z % (r - p + 1);
}

void fill_random_string(unsigned char* &s, long length, int sigma) {
  for (long i = 0; i < length; ++i)
    s[i] = random_int(0, sigma - 1);
}

void fill_random_letters(unsigned char* &s, long n, int sigma) {
  fill_random_string(s, n, sigma);
  for (long i = 0; i < n; ++i)
    s[i] += 'a';
}

std::string random_string_hash() {
  uint64_t hash = (uint64_t)rand() * RAND_MAX + rand();
  std::stringstream ss;
  ss << hash;
  return ss.str();
}

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
    long length = random_long(1, max_length);
    long sigma = random_long(1, max_sigma);

    long n_blocks =random_int(1, 50);
    long max_block_size = (length + n_blocks - 1) / n_blocks;
    long ram_use = std::max(6L, (long)(max_block_size * 5.L));

    if (max_sigma <= 26) fill_random_letters(text, length, sigma);
    else fill_random_string(text, length, sigma);

    // debug //
    // long length = strlen("aaaaaaaaa");
    // strcat((char *)text, "aaaaaaaaa");
    // long ram_use = 39;
    ///////////

    text[length] = 0;
    std::string filename = "/tmp/in" + random_string_hash();
    utils::write_objects_to_file<unsigned char>(text, length, filename);

    // debug //
    // printf("text = %s, length = %ld, ram_use = %ld\n", text, length, ram_use);
    // fflush(stdout);
    ///////////

    // Run the test on generated string.
    SAscan(filename, ram_use);
    
    // Compare the result to correct SA.
    divsufsort64(text, SA, length);
    uint40 *computed_SA;
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
  test_random(5000, 10,       5);
  test_random(5000, 10,      20);
  test_random(5000, 10,     128);
  test_random(5000, 10,     254);
  test_random(500, 100,      5);
  test_random(500, 100,     20);
  test_random(500, 100,    128);
  test_random(500, 100,    254);
  test_random(500, 1000,      5);
  test_random(500, 1000,     20);
  test_random(500, 1000,    128);
  test_random(500, 1000,    254);
  test_random(50, 10000,     5);
  test_random(50, 10000,    20);
  test_random(50, 10000,   128);
  test_random(50, 10000,   254);
  test_random(50, 100000,     5);
  test_random(50, 100000,    20);
  test_random(50, 100000,   128);
  test_random(50, 100000,   254);
  std::fflush(stdout);
}

