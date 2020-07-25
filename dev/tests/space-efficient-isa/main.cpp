#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <string>

#include "divsufsort64.h"
#include "utils.hpp"
#include "rank.h"
#include "approx_rank.hpp"
#include "space_efficient_isa.hpp"

using namespace psascan_private;


void test(
    const std::uint8_t *text,
    std::uint64_t text_length,
    std::uint64_t block_length) {

  const std::uint8_t *block = text;

  // Compute sparse suffix array.
  long *text_sa = new long[text_length];
  divsufsort64(text, text_sa, text_length);

  // Compute BWT.
  std::uint64_t *block_psa = new std::uint64_t[block_length];
  for (std::uint64_t i = 0, j = 0; i < text_length; ++i)
    if ((std::uint64_t)text_sa[i] < block_length)
      block_psa[j++] = text_sa[i];

  std::uint8_t *block_bwt = new std::uint8_t[block_length];
  std::uint64_t i0 = 0;
  for (std::uint64_t i = 0; i < block_length; ++i) {
    if (block_psa[i] == 0) {
      i0 = i;
      block_bwt[i] = 0;
    } else block_bwt[i] = block[block_psa[i] - 1];
  }

  // Build rank over the BWT.
  //typedef rank4n<> rank_type;
  typedef approx_rank<8L> rank_type;
  rank_type *rank =
    new rank_type(block_bwt, block_length);

  std::uint64_t *block_isa = new std::uint64_t[block_length];
  for (std::uint64_t j = 0; j < block_length; ++j)
    block_isa[block_psa[j]] = j;

  typedef space_efficient_isa<rank_type, std::uint64_t, 3>
    isa_type;
  isa_type *sp_isa =
    new isa_type(block_psa, block, rank, block_length, i0);

  for (std::uint64_t q = 0; q < block_length; ++q) {
    std::uint64_t j = utils::random_int64(0, block_length - 1);
    std::uint64_t isa_j = sp_isa->query(j);
    if (isa_j != block_isa[j]) {
      fprintf(stderr, "\nError\n");
      std::exit(EXIT_FAILURE);
    }
  }

  delete sp_isa;
  delete rank;
  delete[] block_psa;
  delete[] block_isa;
  delete[] block_bwt;
  delete[] text_sa;
}

// Test many string chosen according to given paranters.
void test_random(
    std::uint64_t testcases,
    std::uint64_t max_length,
    std::uint64_t max_sigma) {

  fprintf(stderr,"TEST, testcases = %lu, max_n = %lu, max_sigma = %lu\r",
      testcases, max_length, max_sigma);
  std::uint8_t *text = new std::uint8_t[max_length + 1];

  for (std::uint64_t tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {

    // Print progress information.
    if (dbg == 10) {
      fprintf(stderr,"TEST, testcases = %lu, max_n = %lu, max_sigma = %lu: "
        "%lu (%.0Lf%%)\r", testcases, max_length, max_sigma, tc,
        (tc * 100.L) / testcases);
      dbg = 0;
    }

    // Generate string.
    std::uint64_t length = utils::random_int64(2L, max_length);
    std::uint64_t sigma = utils::random_int64(2L, max_sigma);
    std::uint64_t block_size = utils::random_int64(1L, length - 1);

    if (sigma <= 26) utils::fill_random_letters(text, length, sigma);
    else utils::fill_random_string(text, length, sigma);

    // Run the test on generated string.
    test(text, length, block_size);
  }

  // Clean up.
  delete[] text;
  
  fprintf(stderr,"TEST, testcases = %lu, max_n = %lu, max_sigma = %lu: "
      "\033[22;32mPASSED\033[0m%10s\n", testcases, max_length, max_sigma, "");
}

int main(int, char **) {
  std::srand(std::time(0) + getpid());

#ifdef NDEBUG
  test_random(1500,  10,      5);
  test_random(1500,  10,      256);
  test_random(1500,  100,     5);
  test_random(1500,  100,     256);
  test_random(150,   1000,    5);
  test_random(150,   1000,    256);
  test_random(100,   10000,   5);
  test_random(100,   10000,   256);
  test_random(10,    100000,  5);
  test_random(10,    100000,  256);
#else
  test_random(100, 10,      5);
  test_random(100, 10,      256);
  test_random(30,  100,     5);
  test_random(30,  100,     256);
  test_random(10,  1000,    5);
  test_random(10,  1000,    256);
  test_random(5,   10000,   5);
  test_random(5,   10000,   256);
#endif

  fprintf(stderr, "All tests passed.\n");
}

