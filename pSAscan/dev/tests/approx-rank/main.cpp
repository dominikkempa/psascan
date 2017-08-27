#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <string>

#include "utils.hpp"
#include "simple_rank.hpp"
#include "approx_rank.hpp"


using namespace psascan_private;


// Test rank queries on a given strings.
template<std::uint64_t k_sampling_rate_log>
void test(
    const std::uint8_t *text,
    std::uint64_t text_length,
    std::uint64_t queries) {

  static const std::uint64_t k_sampling_rate =
    (std::uint64_t)1 << k_sampling_rate_log;

  simple_rank *rank = new simple_rank(text, (int)text_length);
  typedef approx_rank<k_sampling_rate_log> approx_rank_type;
  approx_rank_type *approx_rank =
    new approx_rank_type(text, text_length);

  for (std::uint64_t q = 0; q < queries; ++q) {
    std::uint64_t i = utils::random_int64(0, 2 * text_length);
    std::uint8_t c = utils::random_int64(0, 255);

    std::uint64_t correct_answer = rank->rank(i, c);
    std::uint64_t computed_answer = approx_rank->rank(i, c);
    std::uint64_t correct_count = rank->rank(text_length, c);
    std::uint64_t computed_count = approx_rank->count(c);

    if (computed_answer > correct_answer ||
        computed_answer + k_sampling_rate <= correct_answer ||
        correct_count != computed_count) {

      fprintf(stderr, "\n\033[22;31mFAILED\033[0m\n");
      if (text_length <= 1000) {
        fprintf(stderr, "  text = ");
        for (std::uint64_t j = 0; j < text_length; ++j)
          fprintf(stderr, "%c", text[j]);
        fprintf(stderr, "\n");
      }

      fprintf(stderr, "  i = %lu, c = %lu\n", i, (std::uint64_t)c);
      fprintf(stderr, "  correct answer = %lu\n", correct_answer);
      fprintf(stderr, "  computed answer = %lu\n", computed_answer);
      fprintf(stderr, "  k_sampling_rate = %lu\n", k_sampling_rate);
      fprintf(stderr, "  correct count = %lu\n", correct_count);
      fprintf(stderr, "  computed count = %lu\n", computed_count);
      std::exit(EXIT_FAILURE);
    }
  }

  delete rank;
  delete approx_rank;
}

// Test many string chosen according to given paranters.
template<std::uint64_t k_sampling_rate_log>
void test_random(
    std::uint64_t testcases,
    std::uint64_t max_length,
    std::uint64_t max_sigma,
    std::uint64_t queries) {

  fprintf(stderr,"TEST, k_srl = %lu, testcases = %lu, max_n = %lu, "
      "max_sigma = %lu, queries = %lu\r", k_sampling_rate_log,
      testcases, max_length, max_sigma, queries);

  std::uint8_t *text = new std::uint8_t[max_length + 1];
  for (std::uint64_t tc = 0; tc < testcases; ++tc) {

    // Print progress information.
    fprintf(stderr,"TEST, k_srl = %lu, testcases = %lu, max_n = %lu, "
        "max_sigma = %lu, queries = %lu: %lu (%.0Lf%%)\r",
        k_sampling_rate_log, testcases, max_length, max_sigma,
        queries, tc, (tc * 100.L) / testcases);

    // Generate string.
    std::uint64_t length = utils::random_int64(1, (std::int64_t)max_length);
    std::uint64_t sigma = utils::random_int64(2, (std::int64_t)max_sigma);
    if (max_sigma <= 26) utils::fill_random_letters(text, length, sigma);
    else utils::fill_random_string(text, length, sigma);

    // Run the test on generated string.
    test<k_sampling_rate_log>(text, length, queries);
  }

  // Clean up.
  delete[] text;

  fprintf(stderr,"TEST, k_srl = %lu, testcases = %lu, max_n = %lu, "
      "max_sigma = %lu, queries = %lu: \033[22;32mPASSED\033[0m%10s\n",
      k_sampling_rate_log, testcases, max_length, max_sigma, queries, "");
}

int main(int, char **) {
  std::srand(std::time(0) + getpid());

  // Run tests.
#ifdef NDEBUG
  fprintf(stderr, "Testing approx rank.\n");
  test_random<0>(10000,        10,      5,  1000);
  test_random<0>(10000,        10,    256,  1000);
  test_random<0>(10000,       100,      5,  1000);
  test_random<0>(10000,       100,    256,  1000);
  test_random<0>(1000,       1000,      5,  1000);
  test_random<0>(1000,       1000,    256,  1000);
  test_random<0>(100,      100000,      5,  1000);
  test_random<0>(100,      100000,    256,  1000);
  test_random<0>(10,      1000000,      5,  1000);
  test_random<0>(10,      1000000,    256,  1000);

  test_random<3>(10000,        10,      5,  1000);
  test_random<3>(10000,        10,    256,  1000);
  test_random<3>(10000,       100,      5,  1000);
  test_random<3>(10000,       100,    256,  1000);
  test_random<3>(1000,       1000,      5,  1000);
  test_random<3>(1000,       1000,    256,  1000);
  test_random<3>(100,      100000,      5,  1000);
  test_random<3>(100,      100000,    256,  1000);
  test_random<3>(10,      1000000,      5,  1000);
  test_random<3>(10,      1000000,    256,  1000);

  test_random<8>(10000,        10,      5,  1000);
  test_random<8>(10000,        10,    256,  1000);
  test_random<8>(10000,       100,      5,  1000);
  test_random<8>(10000,       100,    256,  1000);
  test_random<8>(1000,       1000,      5,  1000);
  test_random<8>(1000,       1000,    256,  1000);
  test_random<8>(100,      100000,      5,  1000);
  test_random<8>(100,      100000,    256,  1000);
  test_random<8>(10,      1000000,      5,  1000);
  test_random<8>(10,      1000000,    256,  1000);

#else
  fprintf(stderr, "Testing approx rank.\n");
  test_random<0>(100,        10,      5,  1000);
  test_random<0>(100,        10,    256,  1000);
  test_random<0>(100,       100,      5,  1000);
  test_random<0>(100,       100,    256,  1000);
  test_random<0>(10,       1000,      5,  1000);
  test_random<0>(10,       1000,    256,  1000);
  test_random<0>(5,      100000,      5,  1000);
  test_random<0>(5,      100000,    256,  1000);
  test_random<0>(5,      1000000,     5,  1000);
  test_random<0>(5,      1000000,   256,  1000);

  test_random<3>(100,        10,      5,  1000);
  test_random<3>(100,        10,    256,  1000);
  test_random<3>(100,       100,      5,  1000);
  test_random<3>(100,       100,    256,  1000);
  test_random<3>(10,       1000,      5,  1000);
  test_random<3>(10,       1000,    256,  1000);
  test_random<3>(5,      100000,      5,  1000);
  test_random<3>(5,      100000,    256,  1000);
  test_random<3>(5,      1000000,     5,  1000);
  test_random<3>(5,      1000000,   256,  1000);

  test_random<8>(100,        10,      5,  1000);
  test_random<8>(100,        10,    256,  1000);
  test_random<8>(100,       100,      5,  1000);
  test_random<8>(100,       100,    256,  1000);
  test_random<8>(10,       1000,      5,  1000);
  test_random<8>(10,       1000,    256,  1000);
  test_random<8>(5,      100000,      5,  1000);
  test_random<8>(5,      100000,    256,  1000);
  test_random<8>(5,      1000000,     5,  1000);
  test_random<8>(5,      1000000,   256,  1000);

#endif


  fprintf(stderr, "All tests passed.\n");
}

