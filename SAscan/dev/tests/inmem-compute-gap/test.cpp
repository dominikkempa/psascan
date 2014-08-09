#include <cstdio>
#include <cstring>
#include <ctime>
#include <unistd.h>

#include "utils.h"
#include "bitvector.h"
#include "inmem_compute_gap.h"
#include "inmem_gap_array.h"
#include "naive_compute_gt.h"
#include "naive_compute_partial_sa.h"
#include "naive_compute_gap.h"


void test(unsigned char *text, long text_length, long left_block_beg,
    long left_block_size, long right_block_size, long max_threads,
    long stream_buffer_size) {

  //
  /*fprintf(stderr, "\n\ntext: ");
  for (long i = 0; i < text_length; ++i)
    fprintf(stderr, "%c", text[i]);
  fprintf(stderr, "\n");
  fprintf(stderr, "length = %ld\n", text_length);
  fprintf(stderr, "left_block_beg = %ld\n", left_block_beg);
  fprintf(stderr, "left_block_size = %ld\n", left_block_size);
  fprintf(stderr, "right_block_size = %ld\n", right_block_size);*/
  //

  //----------------------------------------------------------------------------
  // STEP 1: compute correct answer.
  //----------------------------------------------------------------------------

  // substep 1
  // Correct gap array using naive method.
  long *correct_gap = NULL;
  naive_compute_gap(text, text_length, left_block_beg, left_block_size,
      right_block_size, correct_gap);

  // substep 2
  // Compute corret gt out using naive method.
  bitvector *correct_gt_out = new bitvector(left_block_size + right_block_size + 1);
  naive_compute_gt(text, text_length, left_block_beg,
      left_block_size + right_block_size + 1, correct_gt_out);


  //----------------------------------------------------------------------------
  // Second, compute gap array using parallel streaming.
  //----------------------------------------------------------------------------

  // Prerequisite: gt bitvector.
  bitvector *gt_in = NULL;
  long left_block_end = left_block_beg + left_block_size;
  naive_compute_gt(text, text_length, left_block_end, right_block_size + 1, gt_in);

  
  //
  /*fprintf(stderr, "gt: ");
  for (long i = 0; i < right_block_size + 1; ++i)
    if (gt->get(i)) fprintf(stderr, "1 ");
    else fprintf(stderr, "0 ");
  fprintf(stderr, "\n");*/
  //


  // Prerequisite: partial suffix array.
  int *partial_sa = NULL;
  naive_compute_partial_sa(text, text_length, left_block_beg,
      left_block_end, std::ref(partial_sa));

  //
  /*fprintf(stderr, "partial sa: ");
  for (long i = 0; i < left_block_size; ++i)
    fprintf(stderr, "%ld ", (long)partial_sa[i]);
  fprintf(stderr, "\n");*/
  //

  // Compute the gap array, this is the method we are testing.
  inmem_gap_array *computed_gap = NULL;
  bitvector *computed_gt_out = NULL;
  inmem_compute_gap(text, text_length, left_block_beg, left_block_size,
      right_block_size, partial_sa, gt_in, computed_gt_out, true, computed_gap,
      max_threads, stream_buffer_size);


  //----------------------------------------------------------------------------
  // Compare the answers.
  //----------------------------------------------------------------------------
  long *computed_gap_plain = new long[left_block_size + 1];
  long excess_ptr = 0L;
  for (long i = 0; i < left_block_size + 1; ++i) {
    // Compute gap[i] from the small-gap-array representation.
    long gap_i = computed_gap->m_count[i];
    while (excess_ptr < (long)computed_gap->m_excess.size() &&
        computed_gap->m_excess[excess_ptr] == i) {
      gap_i += 256;
      ++excess_ptr;
    }
    computed_gap_plain[i] = gap_i;
  }

  bool eq = true;
  for (long i = 0; i < left_block_size + right_block_size + 1; ++i)
    if (computed_gt_out->get(i) != correct_gt_out->get(i)) eq = false;

  if (!std::equal(correct_gap, correct_gap + left_block_size + 1,
        computed_gap_plain) || (!eq)) {
    fprintf(stderr, "\nError:\n");
    fprintf(stderr, "\tlength = %ld\n", text_length);
    if (text_length <= 1000) {
      fprintf(stderr, "\ttext: ");
      for (long j = 0; j < text_length; ++j)
        fprintf(stderr, "%c", text[j]);
      fprintf(stderr, "\n");
    }
    fprintf(stderr, "\tleft_block_beg = %ld\n", left_block_beg);
    fprintf(stderr, "\tleft_block_size = %ld\n", left_block_size);
    fprintf(stderr, "\tright_block_size = %ld\n", right_block_size);
    fprintf(stderr, "\tmax_threads = %ld\n", max_threads);
    fprintf(stderr, "\tstream buffer size = %ld\n", stream_buffer_size);
    fprintf(stderr, "\tcorrect gap: ");
    for (long i = 0; i < left_block_size + 1; ++i)
      fprintf(stderr, "%ld ", correct_gap[i]);
    fprintf(stderr, "\n");
    fprintf(stderr, "\tcomputed gap (plain): ");
    for (long i = 0; i < left_block_size + 1; ++i)
      fprintf(stderr, "%ld ", computed_gap_plain[i]);
    fprintf(stderr, "\n");
    fprintf(stderr, "\tcorrect gt out: ");
    for (long i = 0; i < left_block_size + right_block_size + 1; ++i)
      if (correct_gt_out->get(i)) fprintf(stderr, "1 ");
      else fprintf(stderr, "0 ");
    fprintf(stderr, "\n");
    fprintf(stderr, "\tcomputed gt out: ");
    for (long i = 0; i < left_block_size + right_block_size + 1; ++i)
      if (computed_gt_out->get(i)) fprintf(stderr, "1 ");
      else fprintf(stderr, "0 ");
    fprintf(stderr, "\n");
    std::exit(EXIT_FAILURE);
  }


  //----------------------------------------------------------------------------
  // Clean up.
  //----------------------------------------------------------------------------
  delete gt_in;
  delete computed_gt_out;
  delete[] partial_sa;
  delete[] correct_gap;
  delete[] computed_gap_plain;
  delete computed_gap;
  delete correct_gt_out;
}


void test_random(long testcases, long max_length, long max_sigma) {
  fprintf(stderr, "TEST, testcases = %ld, max_n = %ld, max_sigma = %ld\n",
      testcases, max_length, max_sigma);

  unsigned char *text = new unsigned char[max_length];
  for (long tc = 0; tc < testcases; ++tc) {
    // Print progress.
    if (tc % 10 == 0)
      fprintf(stderr, "Progress: %.2Lf%%\r", (100.L * tc) / testcases);

    // Generate input.
    long length = utils::random_long(2L, std::max(2L, max_length));
    long right_block_size = utils::random_long(1L, length / 2);
    long left_block_size = utils::random_long(1L, right_block_size);
    // long left_block_size = utils::random_long(1L, length - 1);
    // long right_block_size = utils::random_long(1L, length - left_block_size);
    long left_block_beg = utils::random_long(0L, length - left_block_size - right_block_size);
    long max_threads = utils::random_long(1L, 50L);
    long stream_buffer_size = utils::random_long(4L, 100L);
    
    long sigma;
    if (length <= 10000) sigma = utils::random_long(1L, max_sigma);
    else sigma = utils::random_long(2L, std::max(2L, max_sigma));
    if (sigma <= 26) utils::fill_random_letters(text, length, sigma);
    else utils::fill_random_string(text, length, sigma);

    /*long length = 6;
    long left_block_beg = 0;
    long left_block_size = 3;
    long right_block_size = 1;
    long max_threads = 1;
    long stream_buffer_size = 34;
    strcpy((char *)text, "cacacb");*/

    // Run the test.
    test(text, length, left_block_beg, left_block_size,
        right_block_size, max_threads, stream_buffer_size);
  }
  delete[] text;
}

int main() {
  std::srand(std::time(0) + getpid());

  test_random(10000,   10,    5);
  test_random(10000,   10,   20);
  test_random(10000,   10,  128);
  test_random(10000,   10,  256);
  test_random(1000,   1000,    5);
  test_random(1000,   1000,   20);
  test_random(1000,   1000,  128);
  test_random(1000,   1000,  256);
  test_random(1000,   10000,    5);
  test_random(1000,   10000,   20);
  test_random(1000,   10000,  128);
  test_random(1000,   10000,  256);
  test_random(1000,   100000,    5);
  test_random(1000,   100000,   20);
  test_random(1000,   100000,  128);
  test_random(1000,   100000,  256);
  test_random(100,   1000000,    5);
  test_random(100,   1000000,   20);
  test_random(100,   1000000,  128);
  test_random(100,   1000000,  256);

  fprintf(stderr, "All tests passed.\n");
}
