#include <cstdio>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <thread>
#include <unistd.h>

#include "compute_initial_gt_bitvectors.h"
#include "initial_partial_sufsort.h"
#include "divsufsort.h"
#include "divsufsort64.h"
#include "bitvector.h"
#include "utils.h"


//==============================================================================
// Fully parallel computation of partial suffix arrays.
//==============================================================================
void partial_sufsort(unsigned char *text, long text_length, long max_threads) {
  long max_block_size = (text_length + max_threads - 1) / max_threads;
  long n_blocks = (text_length + max_block_size - 1) / max_block_size;

  bitvector **gt;
  compute_initial_gt_bitvectors(text, text_length, gt, max_threads);

  // Check the correctness of gt bitvectors.
  for (long i = n_blocks - 1; i >= 0; --i) {
    long block_beg = i * max_block_size;
    long block_end = std::min(block_beg + max_block_size, text_length);
    long block_size = block_end - block_beg;

    for (long j = 0; j < block_size; ++j) {
      // Compute correct gt[j] for block i.
      long lcp = 0L;
      while (block_end + lcp < text_length && text[block_beg + j + lcp] == text[block_end + lcp]) ++lcp;
      bool correct_gt = (block_end + lcp == text_length || text[block_beg + j + lcp] > text[block_end + lcp]);

      if (gt[i]->get(j) != correct_gt) {
        fprintf(stderr, "Error!\n");
        if (text_length <= 1000) {
          fprintf(stderr, "text: ");
          for (long k = 0; k < text_length; ++k)
            fprintf(stderr, "%c", text[k]);
          fprintf(stderr, "\n");
        }
        fprintf(stderr, "max_threads = %ld\n", max_threads);
        fprintf(stderr, "n_blocks = %ld\n", n_blocks);
        fprintf(stderr, "max_block_size = %ld\n",  max_block_size);
        fprintf(stderr, "block id = %ld, j = %ld\n", i, j);
        fprintf(stderr, "correct_gt = %d, computed_gt = %d\n",
          correct_gt, gt[i]->get(j));
        std::exit(EXIT_FAILURE);
      }
    }
  }

  int *full_sa = new int[text_length];
  divsufsort(text, full_sa, (int)text_length);

  unsigned char *original_text = new unsigned char[text_length];
  std::copy(text, text + text_length, original_text);
  int *partial_sa;
  initial_partial_sufsort(text, text_length, gt, partial_sa, max_threads);

  // The the correctness of partial suffix arrays.
  int *temp = new int[max_block_size];
  for (long i = 0; i < n_blocks; ++i) {
    long block_beg = i * max_block_size;
    long block_end = std::min(block_beg + max_block_size, text_length);
    long block_size = block_end - block_beg;
    
    // To verify correctness of partial_sa[i] we collect suffixes starting
    // inside [beg..end) in a temp array and compare to partial_sa[i].
    for (long j = 0, ptr = 0; j < text_length; ++j)
      if (block_beg <= full_sa[j] && full_sa[j] < block_end)
        temp[ptr++] = full_sa[j];

    int *partial_sa_i = partial_sa + block_beg;
    // Compare temp and partial_sa[i].
    for (long j = 0; j < block_size; ++j) {
      if (temp[j] != block_beg + partial_sa_i[j]) {
        fprintf(stderr, "Error!\n");
        if (text_length <= 1000) {
          fprintf(stderr, "text: ");
          for (long jj = 0; jj < text_length; ++jj)
            fprintf(stderr, "%c", text[jj]);
          fprintf(stderr, "\n");
        }
        fprintf(stderr, "block = %ld, beg = %ld, end = %ld\n", i, block_beg, block_end);
        fprintf(stderr, "  temp[%ld] == %d\n", j, temp[j]);
        fprintf(stderr, "  beg + partial_sa[%ld][%ld] == %ld\n", i, j,
            block_beg + partial_sa_i[j]);
        std::exit(EXIT_FAILURE);
      }
    }
  }

  if (!std::equal(original_text, original_text + text_length, text)) {
    fprintf(stderr, "Error!\n");
    fprintf(stderr, "Original text != Text after computation.\n");
    if (text_length <= 1000) {
      fprintf(stderr, "\toriginal text: ");
      for (long jj = 0; jj < text_length; ++jj)
        fprintf(stderr, "%c", original_text[jj]);
      fprintf(stderr, "\n");
      fprintf(stderr, "\ttext after computation: ");
      for (long jj = 0; jj < text_length; ++jj)
        fprintf(stderr, "%c", text[jj]);
      fprintf(stderr, "\n");
    }
    std::exit(EXIT_FAILURE);
  }

  delete[] original_text;
  delete[] full_sa;
  delete[] temp;
  delete[] partial_sa;
}


void test_random(int testcases, long max_length, int max_sigma) {
  fprintf(stderr,"TEST, testcases = %d, max_n = %ld, max_sigma = %d\n",
      testcases, max_length, max_sigma);
  unsigned char *text = new unsigned char[max_length + 1];

  for (int tc = 0; tc < testcases; ++tc) {
    // Print progress information.
    fprintf(stderr,"%d (%.2Lf%%)\r", tc, (tc * 100.L) / testcases);

    // Generate string.
    long length = utils::random_long(1, max_length);
    int sigma = utils::random_int(2, max_sigma);
    if (max_sigma <= 26) utils::fill_random_letters(text, length, sigma);
    else utils::fill_random_string(text, length, sigma);
    
    long max_threads = utils::random_long(1, 50);

    // Run the test on generated string.
    partial_sufsort(text, length, max_threads);
  }

  // Clean up.
  delete[] text;
}


int main() {
  std::srand(std::time(0) + getpid());

  test_random(10000,   10,      5);
  test_random(10000,   10,      255);
  test_random(10000,   1000,    5);
  test_random(10000,   1000,    255);
  test_random(1000,    100000,  5);
  test_random(1000,    100000,  255);
  test_random(1000,     1000000, 5);
  test_random(1000,     1000000, 255);

  fprintf(stderr,"All tests passed.\n");
}
