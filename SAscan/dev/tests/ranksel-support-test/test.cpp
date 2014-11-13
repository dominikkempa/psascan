#include <cstdio>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <thread>
#include <unistd.h>

#include "utils.h"
#include "bitvector.h"
#include "ranksel_support.h"


// valid arguments: 0 <= i <= bv.length.
long naive_rank(bitvector *bv, long i) {
  long result = 0L;
  for (long j = 0; j < i; ++j)
    if (bv->get(j)) ++result;

  return result;
}

// valid arguments: 0 <= i <= bv.length.
long naive_rank0(bitvector *bv, long i) {
  return i - naive_rank(bv, i);
}

// valid arguments: 0 <= i < number of 1-bits in bv.
long naive_select1(bitvector *bv, long i) {
  long rk = 0L, j = 0L;
  while (rk + bv->get(j) <= i)
    rk += bv->get(j++);

  return j;
}

// valid arguments: 0 <= i < number of -bits in bv.
long naive_select0(bitvector *bv, long i) {
  long rk = 0L, j = 0L;
  while (rk + (1 - bv->get(j)) <= i)
    rk += (1 - bv->get(j++));

  return j;
}

void test(bitvector *bv, long length, long max_threads, long queries) {
  ranksel_support *bv_ranksel = new ranksel_support(bv, length, max_threads);
  long bv_ones = 0L;
  for (long j = 0; j < length; ++j)
    if (bv->get(j)) ++bv_ones;
  long bv_zeros = length - bv_ones;

  // rank1 queries.
  for (long q = 0; q < queries; ++q) {
    long i = utils::random_long(0, length);
    if (naive_rank(bv, i) != bv_ranksel->rank(i)) {
      fprintf(stderr, "Error:\n");
      fprintf(stderr, "\trank1(bv, %ld) = %ld (correct)\n", i,  naive_rank(bv, i));
      fprintf(stderr, "\trank1(bv, %ld) = %ld (computed)\n", i, bv_ranksel->rank(i));
      std::exit(EXIT_FAILURE);
    }
  }

  // rank0 queries.
  for (long q = 0; q < queries; ++q) {
    long i = utils::random_long(0, length);
    if (naive_rank0(bv, i) != bv_ranksel->rank0(i)) {
      fprintf(stderr, "Error:\n");
      fprintf(stderr, "\trank0(bv, %ld) = %ld (correct)\n", i,  naive_rank0(bv, i));
      fprintf(stderr, "\trank0(bv, %ld) = %ld (computed)\n", i, bv_ranksel->rank0(i));
      std::exit(EXIT_FAILURE);
    }
  }

  // select0 queries.
  if (bv_zeros > 0) {
    for (long q = 0; q < queries; ++q) {
      long i = utils::random_long(0, bv_zeros - 1);
      if (naive_select0(bv, i) != bv_ranksel->select0(i)) {
        fprintf(stderr, "Error:\n");
        fprintf(stderr, "\tselect0(bv, %ld) = %ld (correct)\n", i,  naive_select0(bv, i));
        fprintf(stderr, "\tselect0(bv, %ld) = %ld (computed)\n", i, bv_ranksel->select0(i));
        std::exit(EXIT_FAILURE);
      }
    }
  }

  // select1 queries.
  if (bv_ones > 0) {
    for (long q = 0; q < queries; ++q) {
      long i = utils::random_long(0, bv_ones - 1);
      if (naive_select1(bv, i) != bv_ranksel->select1(i)) {
        fprintf(stderr, "Error:\n");
        fprintf(stderr, "\tselect1(bv, %ld) = %ld (correct)\n", i,  naive_select1(bv, i));
        fprintf(stderr, "\tselect1(bv, %ld) = %ld (computed)\n", i, bv_ranksel->select1(i));
        std::exit(EXIT_FAILURE);
      }
    }
  }

  delete bv_ranksel;
}

void test_random(int testcases, long max_length, long queries) {
  fprintf(stderr,"TEST, testcases = %d, max_n = %ld, queries = %ld\n", testcases, max_length, queries);

  for (int tc = 0; tc < testcases; ++tc) {
    // Print progress information.
    if (tc % 100 == 0)
      fprintf(stderr,"%d (%.2Lf%%)\r", tc, (tc * 100.L) / testcases);

    // Generate input.
    long length = utils::random_long(1L, max_length);
    bitvector *bv = new bitvector(length);
    for (long j = 0; j < length; ++j)
      if (utils::random_int(0, 1)) bv->set(j);

    long max_threads = utils::random_long(1L, 50L);

    // Run the test on generated input.
    test(bv, length, max_threads, queries);

    delete bv;
  }
}

int main() {
  std::srand(std::time(0) + getpid());

  test_random(100000,  10,        100);
  test_random(10000,   100,       100);
  test_random(1000,    1000,      1000);
  test_random(1000,    10000,     1000);
  test_random(100,     100000,    1000);
  test_random(10,      1000000,   1000);

  fprintf(stderr,"All tests passed.\n");
}

