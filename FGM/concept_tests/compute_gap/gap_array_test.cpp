// Tests the correctness of the streaming in FGM, that is, for every suffix of
// A, where A is long and B short, finds a place in lexicographic range of
// suffixes of BA starting in B.
// This is done mostly to correctly figure out the +/-1 cases and correct
// handling of gt and dollar positions. Also, we wish to plug a standard rank,
// not especially tailored for this. By this I mean, it should return number of
// c's in prefix [0..i). The rank should implement all possible values of i
// (including negative).

#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <string>

#include "sais.hxx"
#include "utils.h"
#include "../rank/general_rank.h"

void test(unsigned char *text, int length, int B_length) {
  int A_length = length - B_length, BA_length = length;
  unsigned char *B = text, *BA = text, *A = text + B_length;

  int *BA_SA = new int[BA_length];
  saisxx(BA, BA_SA, BA_length);

  int *sparseSA = new int[B_length];
  int *gap = new int[B_length + 1];
  std::fill(gap, gap + B_length + 1, 0);
  for (int i = 0, j = 0; i < BA_length; ++i) {
    if (BA_SA[i] < B_length) sparseSA[j++] = BA_SA[i];
    else ++gap[j];
  }

  unsigned char *BWT = new unsigned char[B_length];
  int dollar_pos = 0;
  for (int i = 0, j = 0; i < B_length; ++i)
    if (sparseSA[i] == 0) dollar_pos = i;
    else BWT[j++] = B[sparseSA[i] - 1];

  int count[256] = {0};
  for (int i = 0; i < B_length; ++i) count[B[i] + 1]++;
  for (int i = 1; i < 256; ++i) count[i] += count[i - 1];

  unsigned char *gt = new unsigned char[A_length];
  for (int j = 0; j < A_length; ++j) {
    int lcp = 0;
    while (j + lcp < A_length && A[j + lcp] == A[lcp]) ++lcp;
    gt[j] = (j + lcp < A_length && A[j + lcp] > A[lcp]); // A[j..] > A[0..]?
  }

  int *computed_gap = new int[B_length + 1];               
  std::fill(computed_gap, computed_gap + B_length + 1, 0); 
  int i = 0;
  general_rank *rank = new general_rank(BWT, B_length - 1);
  for (int j = A_length - 1; j >= 0; --j) {
    unsigned char c = A[j];
    i = count[c] + rank->rank(i - (i > dollar_pos), c);    
    if (c == B[B_length - 1] && j + 1 != A_length && gt[j + 1]) ++i;                    
    ++computed_gap[i];

    // Compute (by brute force) the rank (smallest suffix which is larger) of
    // string A[j..A_length) among suffixes of BA starting in B.
    int brute_answer = 0;
    while (brute_answer < B_length) {
      int lcp = 0;
      while (lcp < A_length - j && A[j + lcp] == text[sparseSA[brute_answer] + lcp]) ++lcp;
      if (lcp == A_length - j || text[sparseSA[brute_answer] + lcp] > A[j + lcp]) break;
      ++brute_answer;
    }
    if (i != brute_answer) {
      fprintf(stderr, "Error!\n");
      fprintf(stderr, "  A = ");
      for (int k = 0; k < A_length; ++k) fprintf(stderr, "%c", A[k]);
      fprintf(stderr, "\n");
      fprintf(stderr, "  B = ");
      for (int k = 0; k < B_length; ++k) fprintf(stderr, "%c", B[k]);
      fprintf(stderr, "\n");
      fprintf(stderr, "brute_answer = %d, i = %d\n", brute_answer, i);
      std::exit(EXIT_FAILURE);
    }
  }

  bool eq = true;
  for (int j = 0; j <= B_length; ++j)
    if (gap[j] != computed_gap[j]) eq = false;
  if (!eq) {
      fprintf(stderr, "Error!\n");
      fprintf(stderr, "  A = ");
      for (int k = 0; k < A_length; ++k) fprintf(stderr, "%c", A[k]);
      fprintf(stderr, "\n");
      fprintf(stderr, "  B = ");
      for (int k = 0; k < B_length; ++k) fprintf(stderr, "%c", B[k]);
      fprintf(stderr, "\n");
      fprintf(stderr, "correct gap array: ");
      for (int k = 0; k <= B_length; ++k)
        fprintf(stderr, "%d ", gap[k]);
      fprintf(stderr, "\n");
      fprintf(stderr, "computed gap array: ");
      for (int k = 0; k <= B_length; ++k)
        fprintf(stderr, "%d ", computed_gap[k]);
      fprintf(stderr, "\n");
      std::exit(EXIT_FAILURE);
  }

  delete[] gap;
  delete[] computed_gap;
  delete[] gt;
  delete rank;
  delete[] sparseSA;
  delete[] BA_SA;
}

// Test many string chosen according to given paranters.
void test_random(int testcases, int max_length, int max_sigma) {
  fprintf(stderr,"TEST, testcases = %d, max_n = %d, max_sigma = %d\n",
      testcases, max_length, max_sigma);
  unsigned char *text = new unsigned char[max_length + 1];

  for (int tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {
    // Print progress information.
    if (dbg == 1000) {
      fprintf(stderr,"%d (%.2Lf%%)\r", tc, (tc * 100.L) / testcases);
      dbg = 0;
    }

    // Generate string.
    int length = utils::random_int(2, max_length);
    int sigma = utils::random_int(2, max_sigma);
    int block_size = utils::random_int(1, length - 1);
    if (max_sigma <= 26) utils::fill_random_letters(text, length, sigma);
    else utils::fill_random_string(text, length, sigma);

    // Run the test on generated string.
    test(text, length, block_size);
  }

  // Clean up.
  delete[] text;
}

int main(int, char **) {
  srand(time(0) + getpid());

  // Run tests.
  fprintf(stderr, "Testing gap array computation.\n");
  test_random(500000, 10,      5);
  test_random(500000, 10,      256);
  test_random(500000, 100,     5);
  test_random(500000, 100,    256);
  fprintf(stderr,"All tests passed.\n");

  return 0;
}

