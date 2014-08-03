#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <string>

#include "utils.h"
#include "rank.h"
#include "new_rank.h"

// Test rank queries on a given strings.
void test(unsigned char *text, int length, int queries) {
  unsigned char *text_copy = new unsigned char[length + 1];
  std::copy(text, text + length, text_copy);
  text_copy[length] = 0;
  long max_threads = utils::random_long(1, 24);
  context_rank_4n *correct_rank = new context_rank_4n(text, length);
  rank4n *rank = new rank4n(text, length, max_threads);


  if (correct_rank->n_block != rank->n_blocks) {
    fprintf(stderr, "\nError: different number of blocks\n");
    std::exit(EXIT_FAILURE);
  }

  // compare, something simple first:
  if (!std::equal(correct_rank->m_mapping, correct_rank->m_mapping + 512 * correct_rank->n_block,
        rank->m_mapping)) {
    fprintf(stderr, "\nError: different mappings\n");
    std::exit(EXIT_FAILURE);
  }
  
  if (!std::equal(correct_rank->block_header, correct_rank->block_header + correct_rank->n_block,
        rank->block_header)) {
    fprintf(stderr, "\nError: different block headers\n");
    std::exit(EXIT_FAILURE);
  }

  if (correct_rank->n_sblock != rank->n_sblock) {
    fprintf(stderr, "\nError: different numbers of sblocks\n");
    std::exit(EXIT_FAILURE);
  }

  if (!std::equal(correct_rank->sb_rank, correct_rank->sb_rank + (256 * correct_rank->n_sblock), rank->sb_rank)) {
    fprintf(stderr, "\nError: different sblocks ranks\n");
    std::exit(EXIT_FAILURE);
  }

  if (!std::equal(correct_rank->c_rank, correct_rank->c_rank + 256, rank->c_rank)) {
    fprintf(stderr, "\nError: different c_rank\n");
    std::exit(EXIT_FAILURE);
  }

  if (correct_rank->rare_trunk.size() != rank->rare_trunk.size()) {
    fprintf(stderr, "\nError: incorrect trunk size: %ld != %ld\n",
        correct_rank->rare_trunk.size(), rank->rare_trunk.size());
    std::exit(EXIT_FAILURE);
  }

  if (!std::equal(correct_rank->rare_trunk.begin(), correct_rank->rare_trunk.end(), rank->rare_trunk.begin())) {
    fprintf(stderr, "\nError: rare trunks are different, compared %ld elems\n", (long)correct_rank->rare_trunk.size());
    std::exit(EXIT_FAILURE);
  }

  if (!std::equal(correct_rank->freq_trunk, correct_rank->freq_trunk + correct_rank->n_block * correct_rank->k_block_size,
        rank->freq_trunk)) {
    fprintf(stderr, "\nError: freq trunks are different, was comparing %ld symbols\n",
        correct_rank->n_block * correct_rank->k_block_size);
    std::exit(EXIT_FAILURE);
  }


  for (int q = 0; q < queries; ++q) {
    int i = utils::random_int(-2 * length, 2 * length);
    unsigned char c = utils::random_int(0, 255);

    // First, compute the correct answer.
    int ans = 0;
    for (int j = 0; j < std::min(length, i); ++j)
      if (text_copy[j] == c) ++ans;
      
    // Now, query the tested data structured.
    int rank_ans = rank->rank(i, c);
    if (rank_ans != ans) {
      fprintf(stderr, "\n\033[22;31mFAILED\033[0m\n");
      fprintf(stderr, "  length = %d\n", length);
      if (length <= 1000) {
        text[length] = 0;
        fprintf(stderr, "  text = ");
        for (int j = 0; j < length; ++j)
          fprintf(stderr, "%d ", text_copy[j]);
        fprintf(stderr, "\n");
      }
      fprintf(stderr, "  max_threads = %ld\n", max_threads);
      fprintf(stderr, "  i = %d, c = %d\n", i, c);
      fprintf(stderr, "  correct_ans = %d, rank_ans = %d\n",
          ans, rank_ans);
      std::exit(EXIT_FAILURE);
    }
  }

  delete[] text_copy;
  delete rank;
  delete correct_rank;
}

// Test many string chosen according to given paranters.
void test_random(int testcases, int max_length, int max_sigma, int queries) {
  fprintf(stderr,"TEST, testcases = %d, max_n = %d, max_sigma = %d, queries = %d\r",
      testcases, max_length, max_sigma, queries);
  unsigned char *text = new unsigned char[max_length + 1];

  for (int tc = 0; tc < testcases; ++tc) {
    // Print progress information.
    if (tc % 10 == 0)
      fprintf(stderr,"TEST, testcases = %d, max_n = %d, max_sigma = %d, "
          "queries = %d: %d (%.0Lf%%)\r", testcases, max_length, max_sigma,
          queries, tc, (tc * 100.L) / testcases);

    // Generate string.
    int length = utils::random_int(1, max_length);
    int sigma = utils::random_int(1, max_sigma);
    int freq_sigma = utils::random_int(1, sigma);
    int rare_sigma = sigma - freq_sigma;
    for (int j = 0; j < length; ++j)
      text[j] = utils::random_int(0, freq_sigma - 1);
    for (int j = 0; j < rare_sigma; ++j) {
      int tries = utils::random_int(1, (length + 9) / 10);
      for (int k = 0; k < tries; ++k)
        text[utils::random_int(0, length - 1)] = freq_sigma - 1 + j;
    }

    // Run the test on generated string.
    test(text, length, queries);
  }

  // Clean up.
  delete[] text;

  fprintf(stderr,"TEST, testcases = %d, max_n = %d, max_sigma = %d, "
      "queries = %d: \033[22;32mPASSED\033[0m%10s\n", testcases, max_length,
      max_sigma, queries, "");
}

int main(int, char **) {
  std::srand(std::time(0) + getpid());

  // Run tests.
  test_random(100000,   10,      5,    10000);
  test_random(100000,   10,      20,  10000);
  test_random(100000,   10,      128,  10000);
  test_random(100000,   10,      256,  10000);
  test_random(100000,   1000,    5,    1000);
  test_random(100000,   1000,    20,  1000);
  test_random(100000,   1000,    128,  1000);
  test_random(100000,   1000,    256,  1000);
  test_random(10000,    100000,  5,    100);
  test_random(10000,    100000,  20,  100);
  test_random(10000,    100000,  128,  100);
  test_random(10000,    100000,  256,  100);
  test_random(1000,     1000000, 5,    100);
  test_random(1000,     1000000, 20,  100);
  test_random(1000,     1000000, 128,  100);
  test_random(1000,     1000000, 256,  100);
}
