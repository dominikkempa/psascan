// For any string T, the gt array, gt(T) is a bitvector of length |T|, such
// that gt[i] == 1 (i = 0,1,..,|T|-1) iff T[i..n) > T[0..n).
// In this program we test computing gt(BA), given gt(A) and the array of
// sorted suffixes of BA strating in B.

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
  for (int i = 0, j = 0; i < BA_length; ++i)
    if (BA_SA[i] < B_length) sparseSA[j++] = BA_SA[i];

  unsigned char *BWT = new unsigned char[B_length];
  int dollar_pos = 0;
  for (int i = 0, j = 0; i < B_length; ++i)
    if (sparseSA[i] == 0) dollar_pos = i;
    else BWT[j++] = B[sparseSA[i] - 1];

  // Compute the 'whole_suffix_pos'
  int whole_suffix_pos = 0;
  for (int i = 0; i < B_length; ++i)
    if (sparseSA[i] == 0) whole_suffix_pos = i;

  int count[256] = {0};
  for (int i = 0; i < B_length; ++i) count[B[i] + 1]++;
  for (int i = 1; i < 256; ++i) count[i] += count[i - 1];

  unsigned char *gt_A = new unsigned char[A_length];
  for (int j = 0; j < A_length; ++j) {
    int lcp = 0;
    while (j + lcp < A_length && A[j + lcp] == A[lcp]) ++lcp;
    gt_A[j] = (j + lcp < A_length && A[j + lcp] > A[lcp]); // A[j..] > A[0..]?
  }
 
  // Do the streaming as in the 'gap' test.
  unsigned char *gt_BA = new unsigned char[length];
  int i = 0;
  general_rank *rank = new general_rank(BWT, B_length - 1);
  for (int j = A_length - 1; j >= 0; --j) {
    unsigned char c = A[j];
    i = count[c] + rank->rank(i - (i > dollar_pos), c);
    if (c == B[B_length - 1] && j + 1 != A_length && gt_A[j + 1]) ++i;

    gt_BA[B_length + j] = (i > whole_suffix_pos);
  }
  
  // Compute the correct gt_BA (by brute force).
  // The loop below computes gt_BA[B_length..B_length+A_length).
  unsigned char *correct_gt_BA = new unsigned char[length];
  for (int j = 0; j < length; ++j) {
    int lcp = 0;
    while (j + lcp < length && text[j + lcp] == text[lcp]) ++lcp;
    correct_gt_BA[j] = (j + lcp < length && text[j + lcp] > text[lcp]);
  }
  // Now siply compute gt_BA[0..B_length) -- just scan sparse array
  // starting from 'whole_suffix_pos'.
  for (int k = 0; k <= whole_suffix_pos; ++k) gt_BA[sparseSA[k]] = 0;
  for (int k = whole_suffix_pos + 1; k < B_length; ++k) gt_BA[sparseSA[k]] = 1;

  // Compare
  if (!std::equal(gt_BA, gt_BA + length, correct_gt_BA)) {
    fprintf(stderr, "Error!\n");
    if (length < 100) {
      fprintf(stderr, "  B = ");
      for (int k = 0; k < B_length; ++k) fprintf(stderr, "%c", B[k]);
      fprintf(stderr, "\n");
      fprintf(stderr, "  A = ");
      for (int k = 0; k < A_length; ++k) fprintf(stderr, "%c", A[k]);
      fprintf(stderr, "\n");
      fprintf(stderr, "correct  gt_AB: ");
      for (int k = 0; k < length; ++k) fprintf(stderr, "%d ", correct_gt_BA[k]);
      fprintf(stderr, "\n");
      fprintf(stderr, "computed gt_AB: ");
      for (int k = 0; k < length; ++k) fprintf(stderr, "%d ", gt_BA[k]);
      fprintf(stderr, "\n");
    }  
    std::exit(EXIT_FAILURE);
  }

  delete[] gt_A;
  delete[] gt_BA;
  delete[] correct_gt_BA;
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
  fprintf(stderr, "Testing the updates of gt array.\n");
  test_random(500000, 10,      5);
  test_random(500000, 10,      256);
  test_random(500000, 100,     5);
  test_random(500000, 100,    256);
  fprintf(stderr,"All tests passed.\n");

  return 0;
}

