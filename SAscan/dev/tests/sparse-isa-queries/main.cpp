#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <string>

#include "divsufsort64.h"
#include "utils.h"
#include "rank.h"
#include "sparse_isa.h"


void test(unsigned char *text, long length, long B_length) {
  unsigned char *B = text;

  // 1
  //
  // Compute sparse suffix array.
  long *BA_SA = new long[length];
  divsufsort64(text, BA_SA, length);

  // 2
  //
  // Compute BWT.
  long *sparseSA = new long[B_length];
  for (int i = 0, j = 0; i < length; ++i)
    if (BA_SA[i] < B_length) sparseSA[j++] = BA_SA[i];
  unsigned char *BWT = new unsigned char[B_length];
  long i0 = 0;
  for (long i = 0; i < B_length; ++i) {
    if (sparseSA[i] == 0) { i0 = i; BWT[i] = 0; }
    else BWT[i] = B[sparseSA[i] - 1];
  }

  // 3
  //
  // Build rank over the BWT.
  rank4n<> *rank = new rank4n<>(BWT, B_length, 1);

  long *ISA = new long[B_length];
  for (long j = 0; j < B_length; ++j)
    ISA[sparseSA[j]] = j;

  long max_threads = utils::random_long(1L, 24L);
  sparse_isa<> *sp_isa = new sparse_isa<>(sparseSA, B, B_length, i0, rank, max_threads);

  for (long q = 0; q < B_length; ++q) {
    long j = utils::random_long(0, B_length - 1);
    long isa_j = sp_isa->query(j);
    if (isa_j != ISA[j]) {
      fprintf(stderr, "\nError\n");
      std::exit(EXIT_FAILURE);
    }
  }

  delete sp_isa;
  delete rank;
  delete[] sparseSA;
  delete[] ISA;
  delete[] BWT;
  delete[] BA_SA;
}

// Test many string chosen according to given paranters.
void test_random(int testcases, long max_length, int max_sigma) {
  fprintf(stderr,"TEST, testcases = %d, max_n = %ld, max_sigma = %d\r",
      testcases, max_length, max_sigma);
  unsigned char *text = new unsigned char[max_length + 1];

  for (int tc = 0; tc < testcases; ++tc) {
    // Print progress information.
    fprintf(stderr,"TEST, testcases = %d, max_n = %ld, max_sigma = %d: "
        "%d (%.0Lf%%)\r", testcases, max_length, max_sigma, tc,
        (tc * 100.L) / testcases);

    // Generate string.
    long length = utils::random_long(2L, max_length);
    int sigma = utils::random_int(2, max_sigma);
    long block_size = utils::random_long(1L, length - 1);
    utils::fill_random_string(text, length, sigma);

    // Run the test on generated string.
    test(text, length, block_size);
  }

  // Clean up.
  delete[] text;
  
  fprintf(stderr,"TEST, testcases = %d, max_n = %ld, max_sigma = %d: "
      "\033[22;32mPASSED\033[0m%10s\n", testcases, max_length, max_sigma, "");
}

int main(int, char **) {
  std::srand(std::time(0) + getpid());

  test_random(15000, 10,      5);
  test_random(15000, 10,      256);
  test_random(15000, 100,     5);
  test_random(15000, 100,     256);
  test_random(1500,  1000,    5);
  test_random(1500,  1000,    256);
  test_random(150,   10000,   5);
  test_random(150,   10000,   256);
}

