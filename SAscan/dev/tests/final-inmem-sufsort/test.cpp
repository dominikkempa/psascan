#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <unistd.h>

#include "divsufsort.h"
#include "divsufsort64.h"
#include "utils.h"
#include "final_inmem_sufsort.h"


void test(unsigned char *text, long text_length, long max_threads) {
  //----------------------------------------------------------------------------
  // STEP 1: compute correct answer.
  //----------------------------------------------------------------------------
  int *correct_sa = new int[text_length];
  divsufsort(text, correct_sa, (long)text_length);


  //----------------------------------------------------------------------------
  // STEP 2: run the tested algorithm.
  //----------------------------------------------------------------------------
  // 1) compute gt_in of size block_size, where gt_in[i] == 1 iff
  //    text[block_beg + i..) > text[block_end..).
  int *computed_sa = new int[text_length];
  inmem_sascan(text, text_length, computed_sa, max_threads);


  //----------------------------------------------------------------------------
  // STEP 3: compare answers.
  //----------------------------------------------------------------------------
  if (!std::equal(correct_sa, correct_sa + text_length, computed_sa)) {
    fprintf(stderr, "\nError:\n");
    fprintf(stderr, "\tlength = %ld\n", text_length);
    if (text_length <= 1000) {
      fprintf(stderr, "\ttext: ");
      for (long j = 0; j < text_length; ++j)
        fprintf(stderr, "%c", text[j]);
      fprintf(stderr, "\n");
    }
    fprintf(stderr, "\tmax threads = %ld\n", max_threads);
    fprintf(stderr, "\tcorrect sa: ");
    for (long i = 0; i < text_length; ++i)
      fprintf(stderr, "%d ", correct_sa[i]);
    fprintf(stderr, "\n");
    fprintf(stderr, "\tcomputed sa: ");
    for (long i = 0; i < text_length; ++i)
      fprintf(stderr, "%d ", computed_sa[i]);
    fprintf(stderr, "\n");
    std::exit(EXIT_FAILURE);
  }

  //----------------------------------------------------------------------------
  // STEP 4: clean up.
  //----------------------------------------------------------------------------
  delete[] correct_sa;
  delete[] computed_sa;
}


void test_random(long testcases, long max_length, long max_sigma) {
  fprintf(stderr, "TEST, testcases = %ld, max_length = %ld, max_sigma = %ld\n",
      testcases, max_length, max_sigma);
  unsigned char *text = new unsigned char[max_length];

  // Run tests.
  for (long tc = 0; tc < testcases; ++tc) {
    fprintf(stderr, "Progress: %.2Lf%%\r", (100.L * tc) / testcases);

    // Generate input.
    long length = utils::random_long(1L, max_length);
    long max_threads = utils::random_long(1L, 50L);

    long sigma;
    if (length <= 10000) sigma = utils::random_long(1L, max_sigma);
    else sigma = utils::random_long(2L, max_sigma);
    if (sigma <= 26) utils::fill_random_letters(text, length, sigma);
    else utils::fill_random_string(text, length, sigma);

    /*strcpy((char *)text, "ccbbba");
    long length = 6;
    long max_threads = 24;*/

    /*fprintf(stderr, "length = %ld\n", length);
    fprintf(stderr, "text: ");
    for (long j = 0; j < length; ++j)
      fprintf(stderr, "%c", text[j]);
    fprintf(stderr, "\n");
    fprintf(stderr, "max_threads = %ld\n", max_threads);*/

    test(text, length, max_threads);
  }

  // Clean up.
  delete[] text;
}


int main() {
  std::srand(std::time(0) + getpid());

  test_random(10000, 10, 5);
  test_random(10000, 10, 20);
  test_random(10000, 10, 128);
  test_random(10000, 10, 255);

  test_random(1000, 1000, 5);
  test_random(1000, 1000, 20);
  test_random(1000, 1000, 128);
  test_random(1000, 1000, 255);

  test_random(100, 10000, 5);
  test_random(100, 10000, 20);
  test_random(100, 10000, 128);
  test_random(100, 10000, 255);

  test_random(100, 1000000, 5);
  test_random(100, 1000000, 20);
  test_random(100, 1000000, 128);
  test_random(100, 1000000, 255);

  fprintf(stderr, "All tests passed.\n");
}
