#include <cstdio>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <thread>
#include <unistd.h>

#include "utils.h"
#include "inmem_smaller_suffixes.h"
#include "bitvector.h"
#include "compute_partial_sa.h"
#include "brute.h"

void test(unsigned char *text, long length, long block_beg, long block_end,
    long suf_start) {

  // First, compute the answer naively.
  long correct_answer;
  inmem_smaller_suffixes_brute(text, length, block_beg, block_end,
      suf_start, correct_answer);

  // Second, compute answers using string range matching.
  long computed_answer;
  int *partial_sa;
  compute_partial_sa(text, length, block_beg, block_end, partial_sa);
  inmem_smaller_suffixes(text, length, block_beg, block_end, suf_start,
      partial_sa, computed_answer);

  if (computed_answer != correct_answer) {
    fprintf(stderr, "\nError:\n");
    fprintf(stderr, "\tlength = %ld\n", length);
    fprintf(stderr, "\ttext: ");
    for (long i = 0; i < length; ++i)
      fprintf(stderr, "%c", text[i]);
    fprintf(stderr, "\n");
    fprintf(stderr, "\tblock_beg = %ld\n", block_beg);
    fprintf(stderr, "\tblock_end = %ld\n", block_end);
    fprintf(stderr, "\tsuf_start = %ld\n", suf_start);
    long block_size = block_end - block_beg;
    fprintf(stderr, "\tpartial sa: ");
    for (long j = 0; j < block_size; ++j)
      fprintf(stderr, "%ld ", (long)partial_sa[j]);
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\tcorrect answer = %ld\n", correct_answer);
    fprintf(stderr, "\tcomputed answer = %ld\n", computed_answer);
    std::exit(EXIT_FAILURE);
  }

  delete partial_sa;
}

void test_random(int testcases, long max_length, long max_sigma) {
  fprintf(stderr,"TEST, testcases = %d, max_n = %ld, max_sigma = %ld\n",
    testcases, max_length, max_sigma);

  unsigned char *text = new unsigned char[max_length];
  for (int tc = 0; tc < testcases; ++tc) {
    // Print progress information.
    if (tc % 100 == 0)
      fprintf(stderr,"%d (%.2Lf%%)\r", tc, (tc * 100.L) / testcases);

    // Generate input.
    long length = utils::random_long(1, max_length);
    long block_size = utils::random_long(1, length);
    long block_beg = utils::random_long(0, length - block_size);
    long block_end = block_beg + block_size;
    long suf_start = utils::random_long(block_end, length);

    long sigma;
    if (length >= 10000) sigma = utils::random_long(2, max_sigma);
    else sigma = utils::random_long(1, max_sigma);

    if (sigma <= 26) utils::fill_random_letters(text, length, sigma);
    else utils::fill_random_string(text, length, sigma);

    // Run the test.
    test(text, length, block_beg, block_end, suf_start);
  }

  // Clean up.
  delete[] text;
}

int main() {
  std::srand(std::time(0) + getpid());

  test_random(1000000,  10,    5);
  test_random(1000000,  10,    20);
  test_random(1000000,  10,    128);
  test_random(1000000,  10,    256);
  test_random(1000000,  100,   5);
  test_random(1000000,  100,   20);
  test_random(1000000,  100,   128);
  test_random(1000000,  100,   256);
  test_random(100000,  1000,   5);
  test_random(100000,  1000,   20);
  test_random(100000,  1000,   128);
  test_random(100000,  1000,   256);
  test_random(1000,  100000,  5);
  test_random(1000,  100000,  20);
  test_random(1000,  100000,  128);
  test_random(1000,  100000,  256);
  test_random(1000,  1000000, 5);
  test_random(1000,  1000000, 20);
  test_random(1000,  1000000, 128);
  test_random(1000,  1000000, 256);

  fprintf(stderr,"All tests passed.\n");
}
