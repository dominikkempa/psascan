// Given BA where B is the current FGM block and A' is the previous block
// (recall, that we consider blocks from the end of text, thus A' is the prefix
// of a longer string A), we compute the gt_eof array, where gt_eof[i] == 1 iff
// the suffix of BA starting at pos i (in BA) is larger than A.
// 
// We accomplish this by only having access to B and A'. In addition, we assume
// that we have the prefix of the gt array for A, namely, the array gt_A of
// length |A'|, such that, gt_A[i] == 1 iff the suffix of A starting at i if
// larger than A.
// 
// The computation is done using Crochemore's algorithm: we treat A' as a
// pattern that is searched inside B. We use the version of suffix ranking,
// where the output is available for reading. 

#include <ctime>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <string>

#include "utils.h"

// Update ms-decomposition of T[0..n) from T[0..n-1).
void next(unsigned char *T, int n, int &s, int &p) {
  if (n == 1) { s = 0; p = 1; return; }
  int i = n - 1, r = (i - s) % p;
  while (i < n) {
    unsigned char a = T[s + r], b = T[i];
    if (a > b) { p = i - s + 1; r = 0; }
    else if (a < b) { i -= r; s = i; p = 1; r = 0; }
    else { ++r; if (r == p) r = 0; } ++i;
  }
}

void compute_bitmap(unsigned char *A, int A_length,
                    unsigned char *B, int B_length,
                    unsigned char *gt_A, unsigned char *gt_eof) {
  std::fill(gt_eof, gt_eof + B_length, 0);
  int i = 0, el = 0, s = 0, p = 0, i_max = 0, el_max = 0, s_max = 0, p_max = 0;
  while (i < B_length) {
    while (i + el < B_length && el < A_length && B[i + el] == A[el]) next(A, ++el, s, p);
    if (el == A_length || (i + el == B_length && !gt_A[el]) ||
        (i + el < B_length && B[i + el] > A[el])) gt_eof[i] = 1;
    int j = i_max;
    if (el > el_max) {
      std::swap(el, el_max);
      std::swap(s, s_max);
      std::swap(p, p_max);
      i_max = i;
    }
    if (p && 3 * p <= el && !memcmp(A, A + p, s)) {
      std::copy(gt_eof + j + 1, gt_eof + j + p, gt_eof + i + 1);
      i += p; el -= p;
    } else {
      int h = (el / 3) + 1;
      std::copy(gt_eof + j + 1, gt_eof + j + h, gt_eof + i + 1);
      i += h; el = 0; s = 0; p = 0;
    }
  }
}

void test(unsigned char *text, int length, int B_length) {
  unsigned char *B = text, *A = text + B_length;
  int A_length = length - B_length;

  // Compute gt_A.
  unsigned char *gt_A = new unsigned char[A_length];
  for (int j = 0; j < A_length; ++j) {
    int lcp = 0;
    while (j + lcp < A_length && A[j + lcp] == A[lcp]) ++lcp;
    gt_A[j] = (j + lcp < A_length && A[j + lcp] > A[lcp]); // A[j..] > A[0..]?
  }

  // Compute gt eof using Crochemore's algorithm.
  unsigned char *gt_eof = new unsigned char[B_length];
  compute_bitmap(A, A_length, B, B_length, gt_A, gt_eof);

  // Compute gt_eof by brute force.
  unsigned char *correct_gt_eof = new unsigned char[B_length];
  for (int j = 0; j < B_length; ++j) {
    int lcp = 0;
    while (lcp < A_length && B[j + lcp] == A[lcp]) ++lcp;
    correct_gt_eof[j] = ((lcp == A_length) || (B[j + lcp] > A[lcp]));
  }

  // Compare
  if (!std::equal(gt_eof, gt_eof + B_length, correct_gt_eof)) {
    fprintf(stderr, "Error!\n");
    if (length < 100) {
      fprintf(stderr, "  B = ");
      for (int k = 0; k < B_length; ++k) fprintf(stderr, "%c", B[k]);
      fprintf(stderr, "\n");
      fprintf(stderr, "  A = ");
      for (int k = 0; k < A_length; ++k) fprintf(stderr, "%c", A[k]);
      fprintf(stderr, "\n");
      fprintf(stderr, "correct  gt_eof: ");
      for (int k = 0; k < B_length; ++k) fprintf(stderr, "%d ", correct_gt_eof[k]);
      fprintf(stderr, "\n");
      fprintf(stderr, "computed gt_eof: ");
      for (int k = 0; k < B_length; ++k) fprintf(stderr, "%d ", gt_eof[k]);
      fprintf(stderr, "\n");
    }  
    std::exit(EXIT_FAILURE);
  }

  delete[] gt_A;
  delete[] gt_eof;
  delete[] correct_gt_eof;
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
    text[length] = 0;

    // Run the test on generated string.
    test(text, length, block_size);
  }

  // Clean up.
  delete[] text;
}

int main(int, char **) {
  srand(time(0) + getpid());

  // Run tests.
  fprintf(stderr, "Testing computation of gt_eof array.\n");
  test_random(5000000, 10,      5);
  test_random(5000000, 10,    256);
  test_random(500000, 100,      5);
  test_random(500000, 100,    256);
  test_random(50000, 1000,      5);
  test_random(50000, 1000,    256);
  test_random(50000, 10000,     5);
  test_random(50000, 10000,   256);
  test_random(500, 100000,      5);
  test_random(500, 100000,    256);
  fprintf(stderr,"All tests passed.\n");

  return 0;
}

