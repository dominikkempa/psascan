#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <string>

#include "utils.h"
#include "oldrank.h"
#include "rank.h"

// Test rank queries on a given strings.
void test(unsigned char *text, long length, long queries) {
  context_rank_4n *correct_rank = new context_rank_4n(text, length);
  rank4n<13, 9> *rank = new rank4n<13, 9>(text, length);

  for (long q = 0; q < queries; ++q) {
    long i = utils::random_long(-2 * length, 2 * length);
    unsigned char c = utils::random_int(0, 255);

    // First, compute the correct answer.
    long ans = 0;
    for (long j = 0; j < std::min(length, i); ++j)
      if (text[j] == c) ++ans;
      
    // Now, query the tested data structured.
    long rank_ans = rank->rank(i, c);
    if (rank_ans != ans) {
      fprintf(stderr, "\n\033[22;31mFAILED\033[0m\n");
      fprintf(stderr, "  length = %ld\n", length);
      if (length <= 1000) {
        text[length] = 0;
        fprintf(stderr, "  text = ");
        for (long j = 0; j < length; ++j)
          fprintf(stderr, "%d ", text[j]);
        fprintf(stderr, "\n");
      }
      fprintf(stderr, "  i = %ld, c = %ld\n", i, (long)c);
      fprintf(stderr, "  correct_ans = %ld, rank_and = %ld\n",
          ans, rank_ans);
      fprintf(stderr, "  correct_rank->answer = %ld\n",
          correct_rank->rank(i, c));
      std::exit(EXIT_FAILURE);
    }
  }

  delete rank;
  delete correct_rank;
}

// Test many string chosen according to given paranters.
void test_random(long testcases, long max_length, long max_sigma, long queries) {
  fprintf(stderr,"TEST, testcases = %ld, max_n = %ld, max_sigma = %ld, queries = %ld\r",
      testcases, max_length, max_sigma, queries);
  unsigned char *text = new unsigned char[max_length + 1];

  for (long tc = 0; tc < testcases; ++tc) {
    // Print progress information.
    if (tc % 10 == 0)
      fprintf(stderr,"TEST, testcases = %ld, max_n = %ld, max_sigma = %ld, "
          "queries = %ld: %ld (%.0Lf%%)\r", testcases, max_length, max_sigma,
          queries, tc, (tc * 100.L) / testcases);

    // Generate string.
    long length = utils::random_long(0, max_length);
    long sigma = utils::random_long(1, max_sigma);
    long freq_sigma = utils::random_long(1, sigma);
    long rare_sigma = sigma - freq_sigma;
    if (length > 0) {
      for (long j = 0; j < length; ++j)
        text[j] = utils::random_long(0, freq_sigma - 1);
      for (long j = 0; j < rare_sigma; ++j) {
        long tries = utils::random_long(1, (length + 9) / 10);
        for (long k = 0; k < tries; ++k)
          text[utils::random_long(0, length - 1)] = freq_sigma - 1 + j;
      }
    }

    // Run the test on generated string.
    test(text, length, queries);
  }

  // Clean up.
  delete[] text;

  fprintf(stderr,"TEST, testcases = %ld, max_n = %ld, max_sigma = %ld, "
      "queries = %ld: \033[22;32mPASSED\033[0m%10s\n", testcases, max_length,
      max_sigma, queries, "");
}

void test_big_random(long length, long sigma, long queries) {
  fprintf(stderr, "Generating text: ");
  unsigned char *text = new unsigned char[length];
  for (long i = 0; i < length; ++i)
    text[i] = utils::random_int(0, sigma - 1);
  fprintf(stderr, "DONE\n");

  fprintf(stderr, "Building rank seqentially: ");
  long double start = utils::wclock();
  context_rank_4n *correct_rank = new context_rank_4n(text, length);
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

  fprintf(stderr, "Building rank in parallel: ");
  start = utils::wclock();
  rank4n<> *rank = new rank4n<>(text, length);
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

  fprintf(stderr, "Asking %ld queries:\n", queries);
  for (long q = 0; q < queries; ++q) {
    if (q % 10000 == 0) fprintf(stderr, "\rprogress: %.2Lf",
        (100.L * q) / queries);
    long i = utils::random_long(-2L * length, 2L * length);
    unsigned char c = utils::random_int(0, 255);
    long rank_ans = rank->rank(i, c);
    long ans = correct_rank->rank(i, c);

    if (rank_ans != ans) {
      fprintf(stderr, "\n\033[22;31mFAILED\033[0m\n");
      fprintf(stderr, "  length = %ld\n", length);
      if (length <= 1000) {
        text[length] = 0;
        fprintf(stderr, "  text = ");
        for (long j = 0; j < length; ++j)
          fprintf(stderr, "%d ", text[j]);
        fprintf(stderr, "\n");
      }
      fprintf(stderr, "  i = %ld, c = %ld\n", i, (long)c);
      fprintf(stderr, "  correct_ans = %ld, rank_and = %ld\n",
          ans, rank_ans);
      std::exit(EXIT_FAILURE);
    }
  }
  fprintf(stderr, "\nOK\n");

  delete correct_rank;
  delete rank;
  delete[] text;
}

int main(int, char **) {
  std::srand(std::time(0) + getpid());

  // Run tests.
  test_random(10000,   10,      5,    10000);
  test_random(10000,   10,      20,  10000);
  test_random(10000,   10,      128,  10000);
  test_random(10000,   10,      256,  10000);
  test_random(10000,   1000,    5,    1000);
  test_random(10000,   1000,    20,  1000);
  test_random(10000,   1000,    128,  1000);
  test_random(10000,   1000,    256,  1000);
  test_random(1000,    100000,  5,    100);
  test_random(1000,    100000,  20,  100);
  test_random(1000,    100000,  128,  100);
  test_random(1000,    100000,  256,  100);
  test_random(100,     1000000, 5,    100);
  test_random(100,     1000000, 20,  100);
  test_random(100,     1000000, 128,  100);
  test_random(100,     1000000, 256,  100);
  test_big_random(5L << 30, 256, 100000000);
}
