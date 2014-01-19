#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <string>

#include "utils.h"
#include "general_rank.h"
#include "rank_4n.h"

// Test rank queries on a given strings.
void test(unsigned char *text, int length, int queries) {
  //general_rank *rank = new general_rank(text, length);
  rank_4n *rank = new rank_4n(text, length);
  for (int q = 0; q < queries; ++q) {
    int i = utils::random_int(-2 * length, 2 * length);
    unsigned char c = utils::random_int(0, 255);
    
    // First, compute the correct answer.
    int ans = 0;
    for (int j = 0; j < std::min(length, i); ++j)
      if (text[j] == c) ++ans;
      
    // Now, query the tested data structured.
    int rank_ans = rank->rank(i, c);
    if (rank_ans != ans) {
      fprintf(stderr, "Test failed!\n");
      if (length <= 1000) {
        text[length] = 0;
        fprintf(stderr, "  text = %s\n", text);
      }
      fprintf(stderr, "  i = %d, c = %d\n", i, c);
      fprintf(stderr, "  correct_ans = %d, rank_ans = %d\n",
          ans, rank_ans);
      std::exit(EXIT_FAILURE);
    }
  }

  delete rank;
}

// Test many string chosen according to given paranters.
void test_random(int testcases, int max_length, int max_sigma, int queries) {
  fprintf(stderr,"TEST, testcases = %d, max_n = %d, max_sigma = %d, queries = %d\n",
      testcases, max_length, max_sigma, queries);
  unsigned char *text = new unsigned char[max_length + 1];

  for (int tc = 0; tc < testcases; ++tc) {
    // Print progress information.
    fprintf(stderr,"%d (%.2Lf%%)\r", tc, (tc * 100.L) / testcases);

    // Generate string.
    int length = utils::random_int(1, max_length);
    int sigma = utils::random_int(2, max_sigma);
    if (max_sigma <= 26) utils::fill_random_letters(text, length, sigma);
    else utils::fill_random_string(text, length, sigma);

    // Run the test on generated string.
    test(text, length, queries);
  }

  // Clean up.
  delete[] text;
}

int main(int, char **) {
  srand(time(0) + getpid());

  // Run tests.
  fprintf(stderr, "Testing general rank.\n");
  test_random(1000,   10,      5,    10000);
  test_random(1000,   10,      256,  10000);
  test_random(1000,   1000,    5,    10000);
  test_random(1000,   1000,    256,  10000);
  test_random(100,    100000,  5,    1000);
  test_random(100,    100000,  256,  1000);
  test_random(10,     1000000, 5,    1000);
  test_random(10,     1000000, 256,  1000);
  fprintf(stderr,"All tests passed.\n");

  return 0;
}

