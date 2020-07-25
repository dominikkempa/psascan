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
#include "compute_right_gap.h"
#include "io_streamer.h"

void test(unsigned char *text, long text_length, long left_block_beg, long left_block_size,
    long right_block_size, long max_threads, long ram_budget) {
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
  long *block_gap_array = new long[block_size + 1];
  long gap_val = 0L;
  long gap_ptr = 0L;
  for (long j = 0; j < text_length; ++j) {
    if (block_beg <= text_sa[j] && text_sa[j] < block_end) {
      block_gap_array[gap_ptr++] = gap_val;
      gap_val = 0L;
    } else if (block_end <= text_sa[j])
      ++gap_val;
  }
  block_gap_array[gap_ptr++] = gap_val;
  gap_array_2n *block_gap_2n = new gap_array_2n(block_size + 1);
  for (long j = 0; j < block_size + 1; ++j)
    block_gap_2n->set_count(j, block_gap_array[j]);
  delete block_gap_array;

  // Compute the gap array of the left half block wrt to the
  // right half-block (bitvector representation).
  bitvector *left_gap_bv = new bitvector(block_size + 1);  // +1 is for the sentinel (required by the tested method)
  gap_ptr = 0L;
  for (long j = 0; j < text_length; ++j) {
    if (block_beg <= text_sa[j] && text_sa[j] < block_end) {
      if (right_block_beg <= text_sa[j]) left_gap_bv->set(gap_ptr);
      ++gap_ptr;
    }
  }

  // Compute the gap array of the right half-block wrt the tail.
  long *right_block_correct_gap = new long[right_block_size + 1];
  gap_ptr = 0L;
  gap_val = 0L;
  for (long j = 0; j < text_length; ++j) {
    if (right_block_beg <= text_sa[j] && text_sa[j] < right_block_end) {
      right_block_correct_gap[gap_ptr++] = gap_val;
      gap_val = 0L;
    } else if (right_block_end <= text_sa[j])
      ++gap_val;
  }
  right_block_correct_gap[gap_ptr++] = gap_val;

  // Compute the gap array of the right half-block wrt the tail using the tested method.
  std::string out_filename = "temp" + utils::random_string_hash();
  compute_right_gap(left_block_size, right_block_size, block_gap_2n, left_gap_bv, out_filename, max_threads, ram_budget);

  // Now, read compute gap array from out_filename (encoded using v-byte encoding)
  // and compare to right_block_gap.
  long *right_block_computed_gap = new long[right_block_size + 1];
  vbyte_stream_reader *reader = new vbyte_stream_reader(out_filename);
  for (long j = 0; j < right_block_size + 1; ++j)
    right_block_computed_gap[j] = reader->read();
  delete reader;
  delete left_gap_bv;
  delete block_gap_2n;
  delete text_sa;
  utils::file_delete(out_filename);

  if (!std::equal(right_block_computed_gap, right_block_computed_gap + right_block_size + 1, right_block_correct_gap)) {
    fprintf(stdout, "\n\033[22;31mFAILED\033[0m\n");
    fprintf(stdout, "\ttext = %s\n", text);
    fprintf(stdout, "\ttext_length = %ld\n", text_length);
    fprintf(stdout, "\tleft_block = [%ld..%ld)\n", left_block_beg, left_block_end);
    fprintf(stdout, "\tright_block = [%ld..%ld)\n", right_block_beg, right_block_end);
    fprintf(stdout, "\tcorrect right block gap: ");
    for (long j = 0; j < right_block_size + 1; ++j)
      fprintf(stdout, "%ld ", right_block_correct_gap[j]);
    fprintf(stdout, "\n");
    fprintf(stdout, "\tcomputed right block gap: ");
    for (long j = 0; j < right_block_size + 1; ++j)
      fprintf(stdout, "%ld ", right_block_computed_gap[j]);
    fprintf(stdout, "\n");
    fflush(stdout);
    std::exit(EXIT_FAILURE);
  }

  delete[] right_block_computed_gap;
  delete[] right_block_correct_gap;
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
    long ram_budget = utils::random_long(1L, 24L * (right_block_size + 1));
    if (max_sigma <= 26) utils::fill_random_letters(text, length, sigma);
    else utils::fill_random_string(text, length, sigma);

    // Run the test on generated string.
    test(text, length, left_block_beg, left_block_size, right_block_size, max_threads, ram_budget);
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
  test_random(500000, 10,      5);
  test_random(500000, 10,     20);
  test_random(500000, 10,    128);
  test_random(500000, 10,    256);
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

