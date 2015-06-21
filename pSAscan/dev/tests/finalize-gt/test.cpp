#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <unistd.h>

#include "utils.h"
#include "finalize_gt.h"
#include "naive_finalize_gt.h"

void test(unsigned char *text, long text_length, long left_block_beg,
    long left_block_size, long max_threads) {

  //----------------------------------------------------------------------------
  // STEP 1: compute correct answer.
  //----------------------------------------------------------------------------
  bitvector *correct_gt = new bitvector(left_block_size);
  naive_finalize_gt(text, text_length, left_block_beg, left_block_size, correct_gt);

  //----------------------------------------------------------------------------
  // STEP 2: run the tested algorithm.
  //----------------------------------------------------------------------------

  // 1) compute enougg bits of the gt_in bitvector. We need at most
  //    left_block_size bits, so the size should be min(left_block_size,
  //    text_length - left_block_end), where left_block_end = left_block_beg +
  //    left_block_size.
  long left_block_end = left_block_beg + left_block_size;
  bitvector *gt_in = new bitvector(std::min(left_block_size, text_length - left_block_end));
  for (long i = 0; i < std::min(left_block_size, text_length - left_block_end); ++i) {
    // Compute lcp(text[left_block_end..), text[left_block_end+i..)).
    long lcp = 0L;
    while (left_block_end + i + lcp < text_length && text[left_block_end + lcp] == text[left_block_end + i + lcp])
      ++lcp;

    if (left_block_end + i + lcp != text_length && text[left_block_end + i + lcp] > text[left_block_end + lcp])
      gt_in->set(i);
  }

  // 2) run the tested algorithm
  bitvector *computed_gt = new bitvector(left_block_size);
  finalize_gt(text, text_length, left_block_beg, left_block_size, gt_in, computed_gt, max_threads);

  //----------------------------------------------------------------------------
  // STEP 3: compare answers.
  //----------------------------------------------------------------------------
  bool eq = true;
  for (long j = 0; j < left_block_size; ++j)
    if (correct_gt->get(j) != computed_gt->get(j)) eq = false;
  if (!eq) {
    fprintf(stderr, "\nError:\n");
    fprintf(stderr, "\tlength = %ld\n", text_length);
    if (text_length <= 1000) {
      fprintf(stderr, "\ttext: ");
      for (long j = 0; j < text_length; ++j)
        fprintf(stderr, "%c", text[j]);
      fprintf(stderr, "\n");
    }
    fprintf(stderr, "\tleft block beg = %ld\n", left_block_beg);
    fprintf(stderr, "\tleft block size = %ld\n", left_block_size);
    fprintf(stderr, "\tmax threads = %ld\n", max_threads);
    fprintf(stderr, "\tgt_in: ");
    for (long j = 0; j < std::min(left_block_size, text_length - left_block_end); ++j)
      if (gt_in->get(j)) fprintf(stderr, "1 ");
      else fprintf(stderr, "0 ");
    fprintf(stderr, "\n");
    fprintf(stderr, "\tcorrect answer: ");
    for (long j = 0; j < left_block_size; ++j)
      if (correct_gt->get(j)) fprintf(stderr, "1 ");
      else fprintf(stderr, "0 ");
    fprintf(stderr, "\n");
    fprintf(stderr, "\tcomputed answer: ");
    for (long j = 0; j < left_block_size; ++j)
      if (computed_gt->get(j)) fprintf(stderr, "1 ");
      else fprintf(stderr, "0 ");
    fprintf(stderr, "\n");
    std::exit(EXIT_FAILURE);
  }

  //----------------------------------------------------------------------------
  // STEP 4: clean up.
  //----------------------------------------------------------------------------

  delete gt_in;
  delete correct_gt;
  delete computed_gt;
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
    long left_block_size = utils::random_long(1L, length);
    long left_block_beg = utils::random_long(0L, length - left_block_size);
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

    test(text, length, left_block_beg, left_block_size, max_threads);
  }

  // Clean up.
  delete[] text;
}


int main() {
  std::srand(std::time(0) + getpid());

  test_random(1000000, 10, 5);
  test_random(1000000, 10, 20);
  test_random(1000000, 10, 128);
  test_random(1000000, 10, 256);

  test_random(1000000, 1000, 5);
  test_random(1000000, 1000, 20);
  test_random(1000000, 1000, 128);
  test_random(1000000, 1000, 256);

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
