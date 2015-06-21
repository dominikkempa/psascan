#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <string>

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "sais.hxx"
#include "utils.h"
#include "gap_array.h"
#include "io_streamer.h"

void test(unsigned char *text, long text_length, long left_block_beg, long left_block_size,
    long right_block_size, long max_threads) {
  // Compute the suffix array of text.
  int *text_sa = new int[text_length];
  saisxx(text, text_sa, (int)text_length);

  // Compute the gap array of the block (wrt tail).
  long block_beg = left_block_beg;
  long left_block_end = left_block_beg + left_block_size;
  long right_block_beg = left_block_end;
  long right_block_end = right_block_beg + right_block_size;
  long block_end = right_block_end;
  long block_size = block_end - block_beg;

  // Compute the gap array of the left half block wrt to the
  // right half-block (gap array representation).
  long gap_val = 0L;
  long gap_ptr = 0L;
  buffered_gap_array *left_block_gap = new buffered_gap_array(left_block_size + 1);
  for (long j = 0; j < text_length; ++j) {
    if (block_beg <= text_sa[j] && text_sa[j] < block_end) {
      if (right_block_beg <= text_sa[j]) ++gap_val;
      else {
        left_block_gap->set_count(gap_ptr++, gap_val);
        gap_val = 0;
      }
    }
  }
  left_block_gap->set_count(gap_ptr++, gap_val);

  // Compute the gap array of the left half block wrt to the
  // right half-block (bitvector representation).
  bitvector *left_gap_bv_correct = new bitvector(block_size);
  gap_ptr = 0L;
  for (long j = 0; j < text_length; ++j) {
    if (block_beg <= text_sa[j] && text_sa[j] < block_end) {
      if (right_block_beg <= text_sa[j]) left_gap_bv_correct->set(gap_ptr);
      ++gap_ptr;
    }
  }
  delete text_sa;

  bitvector *left_gap_bv_computed = left_block_gap->convert_to_bitvector(max_threads);
  left_block_gap->erase_disk_excess();
  delete left_block_gap;

  // compare.
  bool ok = true;
  for (long j = 0; j < block_size; ++j)
    if (left_gap_bv_computed->get(j) != left_gap_bv_correct->get(j)) ok = false;

  if (!ok) {
    fprintf(stdout, "\n\033[22;31mFAILED\033[0m\n");
    fprintf(stdout, "\ttext = %s\n", text);
    fprintf(stdout, "\ttext_length = %ld\n", text_length);
    fprintf(stdout, "\tleft_block = [%ld..%ld)\n", left_block_beg, left_block_end);
    fprintf(stdout, "\tright_block = [%ld..%ld)\n", right_block_beg, right_block_end);
    fprintf(stdout, "\tcorrect bitvector: ");
    for (long j = 0; j < block_size; ++j)
      fprintf(stdout, "%ld ", (long)(left_gap_bv_correct->get(j)));
    fprintf(stdout, "\n");
    fprintf(stdout, "\tcomputed bitvector: ");
    for (long j = 0; j < block_size; ++j)
      fprintf(stdout, "%ld ", (long)(left_gap_bv_computed->get(j)));
    fprintf(stdout, "\n");
    fflush(stdout);
    std::exit(EXIT_FAILURE);
  }

  delete left_gap_bv_computed;
  delete left_gap_bv_correct;
}

// Test many string chosen according to given paranters.
void test_random(int testcases, int max_length, int max_sigma) {
  fprintf(stdout, "TEST, testcases = %d, max_n = %d, max_sigma = %d\r",
      testcases, max_length, max_sigma);
  fflush(stdout);
  unsigned char *text = new unsigned char[max_length + 1];

  for (int tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {
    // Print progress information.
    if (dbg == 1000) {
      fprintf(stdout, "TEST, testcases = %d, max_n = %d, max_sigma = %d: "
          "%d (%.0Lf%%)\r", testcases, max_length, max_sigma, tc,
          (tc * 100.L) / testcases);
      fflush(stdout);
      dbg = 0;
    }

    // Generate string.
    long length = utils::random_long(2L, max_length);
    long sigma = utils::random_long(1L, max_sigma);
    long block_size = utils::random_long(2L, length);
    long left_block_size = utils::random_long(1L, block_size - 1);
    long right_block_size = block_size - left_block_size;
    long left_block_beg = utils::random_long(0L, length - block_size);
    long max_threads = utils::random_long(1L, 50L);
    if (max_sigma <= 26) utils::fill_random_letters(text, length, sigma);
    else utils::fill_random_string(text, length, sigma);

    // Run the test on generated string.
    test(text, length, left_block_beg, left_block_size, right_block_size, max_threads);
  }

  // Clean up.
  delete[] text;
  
  fprintf(stdout, "TEST, testcases = %d, max_n = %d, max_sigma = %d: "
      "\033[22;32mPASSED\033[0m%10s\n", testcases, max_length, max_sigma, "");
  fflush(stdout);
}

int main(int, char **) {
  std::srand(std::time(0) + getpid());
  
  // Redirect stderr to /dev/null
  int redir = open("/dev/null", O_WRONLY);
  dup2(redir, 2);
  close(redir);

  // Run tests
  test_random(50000, 10,       5);
  test_random(50000, 10,      20);
  test_random(50000, 10,     128);
  test_random(50000, 10,     256);
  test_random(50000, 100,      5);
  test_random(50000, 100,     20);
  test_random(50000, 100,    128);
  test_random(50000, 100,    256);
  test_random(5000, 1000,      5);
  test_random(5000, 1000,     20);
  test_random(5000, 1000,    128);
  test_random(5000, 1000,    256);
  test_random(500, 10000,      5);
  test_random(500, 10000,     20);
  test_random(500, 10000,    128);
  test_random(500, 10000,    256);

  fprintf(stdout, "All tests passed.\n");
  fflush(stdout);
}

