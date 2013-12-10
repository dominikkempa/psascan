// We are given a string BA, where B is short and A long. Suppose we know, for
// every suffix starting in B whether it is larger than A or not. This information
// is stored in gt_eof array: gt_eof[i] == 1 iff B[i..] > A.
//
// Having all that, we wish to compute the ordering of suffixes of BA starting in B.
// We want to remap the symbols of B, such that sorting this modified B is equivalent
// to sorting suffixes of BA starting in B.

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

  unsigned char *Bcopy = new unsigned char[B_length];
  std::copy(B, B + B_length, Bcopy);

  int *BA_SA = new int[BA_length];
  saisxx(BA, BA_SA, BA_length);

  int *sparseSA = new int[B_length];
  for (int i = 0, j = 0; i < BA_length; ++i)
    if (BA_SA[i] < B_length) sparseSA[j++] = BA_SA[i];

  /////
  --A;
  ++A_length;
  /////

  unsigned char *gt_eof = new unsigned char[B_length];
  for (int j = 0; j < B_length; ++j) {
    int lcp = 0;
    while (lcp < A_length && B[j + lcp] == A[lcp]) ++lcp;
    gt_eof[j] = ((lcp == A_length) || B[j + lcp] > A[lcp]); // B[j..] > A?
  }

  // Remap symbols in B. We assume that the maximal symbol is <= 254.
/*  unsigned char last = B[B_length - 1];
  for (int i = 0; i < B_length - 1; ++i)
    if (B[i] > last || (B[i] == last && gt_eof[i + 1])) B[i] += 1;
  ++B[B_length - 1];*/

  for (int i = 0; i < B_length; ++i) B[i] += gt_eof[i]; // new! YAY!


  // Compute the SA for modified B.
  int *B_SA = new int[B_length];
  saisxx(B, B_SA, B_length);

  // Compare B_SA and sparseSA.
  if (!std::equal(B_SA, B_SA + B_length, sparseSA)) {
    fprintf(stderr, "Failure:\n");
    if (length < 100) {
      fprintf(stderr, "  B = ");
      for (int k = 0; k < B_length; ++k) fprintf(stderr, "%c", Bcopy[k]);
      fprintf(stderr, "\n");
      fprintf(stderr, "  B'= ");
      for (int k = 0; k < B_length; ++k) fprintf(stderr, "%c", B[k]);
      fprintf(stderr, "\n");
      fprintf(stderr, "  A = ");
      for (int k = 0; k < A_length; ++k) fprintf(stderr, "%c", A[k]);
      fprintf(stderr, "\n");
      fprintf(stderr, "  correct sparseSA:  ");
      for (int k = 0; k < B_length; ++k) fprintf(stderr, "%d ", sparseSA[k]);
      fprintf(stderr, "\n");
      fprintf(stderr, "  computed sparseSA: ");
      for (int k = 0; k < B_length; ++k) fprintf(stderr, "%d ", B_SA[k]);
      fprintf(stderr, "\n");
    }
    std::exit(EXIT_FAILURE);
  }

  // Clean up.
  delete[] B_SA;
  delete[] gt_eof;
  delete[] sparseSA;
  delete[] BA_SA;
  delete[] Bcopy;
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
  fprintf(stderr, "Testing the SA correction.\n");
  test_random(500000, 10,      5);
  test_random(500000, 10,    255);
  test_random(500000, 100,     5);
  test_random(500000, 100,   255);
  test_random(50000,  1000,    5);
  test_random(50000,  1000,  255);
  test_random(500,    10000,   5);
  test_random(500,    10000, 255);
  fprintf(stderr,"All tests passed.\n");

  return 0;
}

