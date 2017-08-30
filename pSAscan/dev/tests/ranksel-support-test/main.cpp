#include <cstdio>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <thread>
#include <unistd.h>

#include "utils.hpp"
#include "bitvector.hpp"
#include "ranksel_support.hpp"


using namespace psascan_private;

void print_bitvector(
    const bitvector * const bv,
    const std::uint64_t length) {
  fprintf(stderr, "\tbv: ");
  for (std::uint64_t i = 0; i < length; ++i)
    fprintf(stderr, "%lu", (std::uint64_t)bv->get(i));
  fprintf(stderr, "\n");
}

// Valid arguments: 0 <= i <= bv.length.
std::uint64_t naive_rank(
    const bitvector * const bv,
    const std::uint64_t i) {

  std::uint64_t result = 0;
  for (std::uint64_t j = 0; j < i; ++j)
    if (bv->get(j))
      ++result;

  return result;
}

// Valid arguments: 0 <= i <= bv.length.
std::uint64_t naive_rank0(
    const bitvector * const bv,
    const std::uint64_t i) {

  return i - naive_rank(bv, i);
}

// Valid arguments: 0 <= i < number of 1-bits in bv.
std::uint64_t naive_select1(
    const bitvector * const bv,
    const std::uint64_t i) {

  std::uint64_t rk = 0, j = 0;
  while (rk + bv->get(j) <= i)
    rk += bv->get(j++);

  return j;
}

// Valid arguments: 0 <= i < number of -bits in bv.
std::uint64_t naive_select0(
    const bitvector * const bv,
    const std::uint64_t i) {
  std::uint64_t rk = 0, j = 0;
  while (rk + (1 - bv->get(j)) <= i)
    rk += (1 - bv->get(j++));

  return j;
}

template<typename ranksel_support_type>
void test(
    const bitvector * const bv,
    const std::uint64_t length,
    const std::uint64_t queries) {

  ranksel_support_type *bv_ranksel =
    new ranksel_support_type(bv, length);

  std::uint64_t bv_ones = 0;
  for (std::uint64_t j = 0; j < length; ++j)
    if (bv->get(j)) ++bv_ones;
  const std::uint64_t bv_zeros = length - bv_ones;

  // rank1 queries.
  for (std::uint64_t q = 0; q < queries; ++q) {
    const std::uint64_t i = utils::random_int64(0, length);
    if (naive_rank(bv, i) != bv_ranksel->rank(i)) {
      fprintf(stderr, "Error:\n");
      print_bitvector(bv, length);
      fprintf(stderr, "\trank1(bv, %ld) = %ld (correct)\n",
          i,  naive_rank(bv, i));
      fprintf(stderr, "\trank1(bv, %ld) = %ld (computed)\n",
          i, bv_ranksel->rank(i));
      std::exit(EXIT_FAILURE);
    }
  }

  // rank0 queries.
  for (std::uint64_t q = 0; q < queries; ++q) {
    const std::uint64_t i = utils::random_int64(0, length);
    if (naive_rank0(bv, i) != bv_ranksel->rank0(i)) {
      fprintf(stderr, "Error:\n");
      print_bitvector(bv, length);
      fprintf(stderr, "\trank0(bv, %ld) = %ld (correct)\n",
          i,  naive_rank0(bv, i));
      fprintf(stderr, "\trank0(bv, %ld) = %ld (computed)\n",
          i, bv_ranksel->rank0(i));
      std::exit(EXIT_FAILURE);
    }
  }

  // select0 queries.
  if (bv_zeros > 0) {
    for (std::uint64_t q = 0; q < queries; ++q) {
      const std::uint64_t i = utils::random_int64(0, bv_zeros - 1);
      if (naive_select0(bv, i) != bv_ranksel->select0(i)) {
        fprintf(stderr, "Error:\n");
        print_bitvector(bv, length);
        fprintf(stderr, "\tselect0(bv, %ld) = %ld (correct)\n",
            i,  naive_select0(bv, i));
        fprintf(stderr, "\tselect0(bv, %ld) = %ld (computed)\n",
            i, bv_ranksel->select0(i));
        std::exit(EXIT_FAILURE);
      }
    }
  }

  // select1 queries.
  if (bv_ones > 0) {
    for (std::uint64_t q = 0; q < queries; ++q) {
      const std::uint64_t i = utils::random_int64(0, bv_ones - 1);
      if (naive_select1(bv, i) != bv_ranksel->select1(i)) {
        fprintf(stderr, "Error:\n");
        print_bitvector(bv, length);
        fprintf(stderr, "\tselect1(bv, %ld) = %ld (correct)\n",
            i,  naive_select1(bv, i));
        fprintf(stderr, "\tselect1(bv, %ld) = %ld (computed)\n",
            i, bv_ranksel->select1(i));
        std::exit(EXIT_FAILURE);
      }
    }
  }

  delete bv_ranksel;
}

template<std::uint64_t k_sampling_rate_log>
void test_random(
    const std::uint64_t testcases,
    const std::uint64_t max_length,
    const std::uint64_t queries) {

  fprintf(stderr,"TEST, k_srl = %lu, testcases = %lu, "
      "max_n = %lu, queries = %lu\n",
      k_sampling_rate_log, testcases, max_length, queries);

  for (std::uint64_t tc = 0; tc < testcases; ++tc) {

    // Print progress information.
    if (tc % 100 == 0)
      fprintf(stderr,"%lu (%.2Lf%%)\r", tc, (tc * 100.L) / testcases);

    // Generate input.
    const std::uint64_t length = utils::random_int64(1, max_length);

    bitvector *bv = new bitvector(length);
    for (std::uint64_t j = 0; j < length; ++j)
      if (utils::random_int32(0, 1)) bv->set(j);

    typedef std::uint64_t offset_type;
    typedef ranksel_support<offset_type, k_sampling_rate_log>
      ranksel_support_type;

    // Run the test on generated input.
    test<ranksel_support_type>(bv, length, queries);

    delete bv;
  }
}

int main() {
  std::srand(std::time(0) + getpid());

#ifdef NDEBUG
  test_random<0>(100000,  10,        100);
  test_random<0>(10000,   100,       100);
  test_random<0>(1000,    1000,      1000);
  test_random<0>(1000,    10000,     1000);
  test_random<0>(100,     100000,    1000);
  test_random<0>(10,      1000000,   1000);

  test_random<2>(100000,  10,        100);
  test_random<2>(10000,   100,       100);
  test_random<2>(1000,    1000,      1000);
  test_random<2>(1000,    10000,     1000);
  test_random<2>(100,     100000,    1000);
  test_random<2>(10,      1000000,   1000);

  test_random<5>(100000,  10,        100);
  test_random<5>(10000,   100,       100);
  test_random<5>(1000,    1000,      1000);
  test_random<5>(1000,    10000,     1000);
  test_random<5>(100,     100000,    1000);
  test_random<5>(10,      1000000,   1000);
#else
  test_random<0>(1000,  10,        100);
  test_random<0>(100,   100,       100);
  test_random<0>(10,    1000,      1000);
  test_random<0>(10,    10000,     1000);

  test_random<2>(1000,  10,        100);
  test_random<2>(100,   100,       100);
  test_random<2>(10,    1000,      1000);
  test_random<2>(10,    10000,     1000);
  
  test_random<5>(1000,  10,        100);
  test_random<5>(100,   100,       100);
  test_random<5>(10,    1000,      1000);
  test_random<5>(10,    10000,     1000);
#endif

  fprintf(stderr,"All tests passed.\n");
}

