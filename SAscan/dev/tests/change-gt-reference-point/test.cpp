#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <unistd.h>

#include "utils.h"
#include "change_gt_reference_point.h"
#include "naive_change_gt_reference_point.h"
#include "naive_compute_gt_in.h"


void test(unsigned char *text, long text_length, long block_beg,
    long block_end, long max_threads) {
  long block_size = block_end - block_beg;


  //----------------------------------------------------------------------------
  // STEP 1: compute correct answer.
  //----------------------------------------------------------------------------
  bitvector *correct_gt_out;
  naive_change_gt_reference_point(text, text_length, block_beg, block_end,
      correct_gt_out);


  //----------------------------------------------------------------------------
  // STEP 2: run the tested algorithm.
  //----------------------------------------------------------------------------
  // 1) compute gt_in of size block_size, where gt_in[i] == 1 iff
  //    text[block_beg + i..) > text[block_end..).
  bitvector *gt_in;
  naive_compute_gt_in(text, text_length, block_beg, block_end, gt_in);

  // 2) run the tested algorithm.
  bitvector *computed_gt_out;
  change_gt_reference_point(text, text_length, block_beg, block_end, gt_in,
      computed_gt_out, max_threads);


  //----------------------------------------------------------------------------
  // STEP 3: compare answers.
  //----------------------------------------------------------------------------
  bool eq = true;
  for (long j = 0; j < block_size + 1; ++j)
    if (correct_gt_out->get(j) != computed_gt_out->get(j))
      eq = false;

  if (!eq) {
    fprintf(stderr, "\nError:\n");
    fprintf(stderr, "\tlength = %ld\n", text_length);
    if (text_length <= 1000) {
      fprintf(stderr, "\ttext: ");
      for (long j = 0; j < text_length; ++j)
        fprintf(stderr, "%c", text[j]);
      fprintf(stderr, "\n");
    }
    fprintf(stderr, "\tblock beg = %ld\n", block_beg);
    fprintf(stderr, "\tblock size = %ld\n", block_size);
    fprintf(stderr, "\tmax threads = %ld\n", max_threads);
    fprintf(stderr, "\tgt_in: ");
    for (long i = 0; i < block_size; ++i)
      if (gt_in->get(i)) fprintf(stderr, "1 ");
      else fprintf(stderr, "0 ");
    fprintf(stderr, "\n");
    fprintf(stderr, "\tcorrect gt out: ");
    for (long i = 0; i < block_size + 1; ++i)
      if (correct_gt_out->get(i)) fprintf(stderr, "1 ");
      else fprintf(stderr, "0 ");
    fprintf(stderr, "\n");
    fprintf(stderr, "\tcomputed gt out: ");
    for (long i = 0; i < block_size + 1; ++i)
      if (computed_gt_out->get(i)) fprintf(stderr, "1 ");
      else fprintf(stderr, "0 ");
    fprintf(stderr, "\n");
    std::exit(EXIT_FAILURE);
  }

  //----------------------------------------------------------------------------
  // STEP 4: clean up.
  //----------------------------------------------------------------------------
  delete gt_in;
  delete correct_gt_out;
  delete computed_gt_out;
}


void test_random(long testcases, long max_length, long max_sigma) {
  fprintf(stderr, "TEST, testcases = %ld, max_length = %ld, max_sigma = %ld\n",
      testcases, max_length, max_sigma);
  unsigned char *text = new unsigned char[max_length];

  // Run tests.
  for (long tc = 0; tc < testcases; ++tc) {
    if (tc % 10 == 0)
      fprintf(stderr, "Progress: %.2Lf%%\r", (100.L * tc) / testcases);

    // Generate input.
    long length = utils::random_long(1L, max_length);
    long block_size = utils::random_long(1L, length);
    long block_beg = utils::random_long(0L, length - block_size);
    long block_end = block_beg + block_size;
    long max_threads = utils::random_long(1L, 50L);

    long sigma;
    if (length <= 10000) sigma = utils::random_long(1L, max_sigma);
    else sigma = utils::random_long(2L, max_sigma);
    if (sigma <= 26) utils::fill_random_letters(text, length, sigma);
    else utils::fill_random_string(text, length, sigma);

    /*strcpy((char *)text, "babbabababbababaaaaabab");
    long length = 23;
    long left_block_beg = 4;
    long left_block_size = 15;
    long max_threads = 12;*/

    test(text, length, block_beg, block_end, max_threads);
  }

  // Clean up.
  delete[] text;
}


int main() {
  std::srand(std::time(0) + getpid());

  test_random(100000, 10, 5);
  test_random(100000, 10, 20);
  test_random(100000, 10, 128);
  test_random(100000, 10, 256);

  test_random(100000, 1000, 5);
  test_random(100000, 1000, 20);
  test_random(100000, 1000, 128);
  test_random(100000, 1000, 256);

  test_random(10000, 10000, 5);
  test_random(10000, 10000, 20);
  test_random(10000, 10000, 128);
  test_random(10000, 10000, 256);

  test_random(1000, 1000000, 5);
  test_random(1000, 1000000, 20);
  test_random(1000, 1000000, 128);
  test_random(1000, 1000000, 256);

  fprintf(stderr, "All tests passed.\n");
}
