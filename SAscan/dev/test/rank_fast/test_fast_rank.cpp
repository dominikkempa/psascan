#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <string>

#include "utils.h"
#include "fast_rank.h"

// Test rank queries on a given strings.
void test(unsigned char *text, int length, int queries) {
  unsigned char *text_copy = new unsigned char[length + 1];
  std::copy(text, text + length, text_copy);
  text_copy[length] = 0;
  fast_rank_4n *rank = new fast_rank_4n(text, length);
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
      fprintf(stderr, "  i = %d, c = %d\n", i, c);
      fprintf(stderr, "  correct_ans = %d, rank_ans = %d\n",
          ans, rank_ans);
      std::exit(EXIT_FAILURE);
    }
  }

  delete[] text_copy;
  delete rank;
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
    utils::fill_random_string(text, length, sigma);

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
  fprintf(stderr, "Testing fast rank.\n");
  test_random(10000,   10,      5,    10000);
  test_random(10000,   10,      256,  10000);
  test_random(10000,   1000,    5,    1000);
  test_random(10000,   1000,    256,  1000);
  test_random(1000,    100000,  5,    100);
  test_random(1000,    100000,  256,  100);
  test_random(100,     1000000, 5,    100);
  test_random(100,     1000000, 256,  100);
}

