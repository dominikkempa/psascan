#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <string>

#include "divsufsort64.h"
#include "utils.h"
#include "rank.h"
#include "approx_rank.h"
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
  //typedef rank4n<> rank_type;
  typedef approx_rank<8L> rank_type;
  long max_threads = utils::random_long(1L, 24L);
  rank_type *rank = new rank_type(BWT, B_length, max_threads);

  long *ISA = new long[B_length];
  for (long j = 0; j < B_length; ++j)
    ISA[sparseSA[j]] = j;

  typedef sparse_isa<rank_type, long, 3L> isa_type;
  isa_type *sp_isa = new isa_type(sparseSA, B, B_length, i0, rank, max_threads);

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
void test_random(long testcases, long max_length, long max_sigma) {
  fprintf(stderr,"TEST, testcases = %ld, max_n = %ld, max_sigma = %ld\r",
      testcases, max_length, max_sigma);
  unsigned char *text = new unsigned char[max_length + 1];

  for (long tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {
    // Print progress information.
    if (dbg == 10) {
      fprintf(stderr,"TEST, testcases = %ld, max_n = %ld, max_sigma = %ld: "
        "%ld (%.0Lf%%)\r", testcases, max_length, max_sigma, tc,
        (tc * 100.L) / testcases);
      dbg = 0;
    }

    // Generate string.
    long length = utils::random_long(2L, max_length);
    long sigma = utils::random_long(2L, max_sigma);
    long block_size = utils::random_long(1L, length - 1);

    if (sigma <= 26) utils::fill_random_letters(text, length, sigma);
    else utils::fill_random_string(text, length, sigma);

    // Run the test on generated string.
    test(text, length, block_size);
  }

  // Clean up.
  delete[] text;
  
  fprintf(stderr,"TEST, testcases = %ld, max_n = %ld, max_sigma = %ld: "
      "\033[22;32mPASSED\033[0m%10s\n", testcases, max_length, max_sigma, "");
}

int main(int, char **) {
  std::srand(std::time(0) + getpid());

  test_random(150000,  10,      5);
  test_random(150000,  10,      256);
  test_random(150000,  100,     5);
  test_random(150000,  100,     256);
  test_random(15000,   1000,    5);
  test_random(15000,   1000,    256);
  test_random(1000,    10000,   5);
  test_random(1000,    10000,   256);
  test_random(100,     100000,  5);
  test_random(100,     100000,  256);

  fprintf(stderr, "All tests passed.\n");
}

