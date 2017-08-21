#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <ctime>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "../../src/psascan_src/bitvector.hpp"
#include "inmem_psascan.hpp"
#include "utils.hpp"


template<typename text_offset_type>
void next(std::uint8_t *text,
    text_offset_type length,
    text_offset_type &s,
    text_offset_type &p,
    text_offset_type &r) {

  if (length == 1) {
    s = 0;
    p = 1;
    r = 0;
    return;
  }

  text_offset_type i = length - 1;
  while (i < length) {
    std::uint8_t a = text[s + r], b = text[i];
    if (a > b) { p = i - s + 1; r = 0; }
    else if (a < b) { i -= r; s = i; p = 1; r = 0; }
    else { ++r; if (r == p) r = 0; } ++i;
  }
}


//=============================================================================
// Compute gt_begin for text.
//=============================================================================
std::uint64_t compute_gt_begin(
    std::uint8_t *text,
    std::uint64_t length,
    psascan_private::bitvector *gt_begin) {

  std::uint64_t whole_suffix_rank = length - 1;
  std::uint64_t i = 1, el = 0, s = 0, p = 0, r = 0;
  std::uint64_t i_max = 0, el_max = 0, s_max = 0, p_max = 0, r_max = 0;
  while (i < length) {
    while (i + el < length && el < length && text[i + el] == text[el])
      next(text, ++el, s, p, r);
 
    if (i + el < length && (el == length || text[i + el] > text[el])) {
      gt_begin->set(i - 1);
      --whole_suffix_rank;
    }

    std::uint64_t j = i_max;
    if (el > el_max) {
      std::swap(el, el_max);
      std::swap(s, s_max);
      std::swap(p, p_max);
      std::swap(r, r_max);
      i_max = i;
    }

    if (p && 3 * p <= el && !memcmp(text, text + p, s)) {
      for (std::uint64_t k = 1; k < std::min(length - i, p); ++k) { 
        if (gt_begin->get(j + k - 1)) {
          gt_begin->set(i + k - 1);
          --whole_suffix_rank;
        }
      }

      i += p;
      el -= p;
    } else {
      std::uint64_t h = (el / 3) + 1;
      for (std::uint64_t k = 1; k < std::min(length - i, h); ++k) { 
        if (gt_begin->get(j + k - 1)) {
          gt_begin->set(i + k - 1);
          --whole_suffix_rank;
        }
      }

      i += h;
      el = 0;
      s = 0;
      p = 0;
    }
  }

  return whole_suffix_rank;
}


template<unsigned pagesize_log>
void test(
    std::uint8_t *text,
    std::uint64_t text_length,
    std::uint64_t max_threads) {

  //---------------------------------------------------------------------------
  // STEP 1: compute correct answer.
  //---------------------------------------------------------------------------
  psascan_private::bitvector correct_gt_begin(text_length);
  compute_gt_begin(text, text_length, &correct_gt_begin);


  //---------------------------------------------------------------------------
  // STEP 2: run the tested algorithm.
  //---------------------------------------------------------------------------
  // 1) compute gt_in of size block_size, where gt_in[i] == 1 iff
  //    text[block_beg + i..) > text[block_end..).
  std::uint8_t *computed_sa_temp =
    new std::uint8_t[text_length * (sizeof(int) + 1)];
  psascan_private::bitvector computed_gt_begin(text_length);
  inmem_psascan<int, pagesize_log>(text, text_length, computed_sa_temp,
      max_threads, false, true, &computed_gt_begin);

  //---------------------------------------------------------------------------
  // STEP 3: compare answers.
  //---------------------------------------------------------------------------
  bool eq = true;
  for (std::uint64_t i = 0; i < text_length; ++i) {
    if (correct_gt_begin.get(i) !=
        computed_gt_begin.get(text_length - 1 - i)) {
      eq = false;
      break;
    }
  }
  
  if (!eq) {
    fprintf(stdout, "\nError:\n");
    fprintf(stdout, "\tlength = %lu\n", text_length);
    if (text_length <= 1000) {
      fprintf(stdout, "\ttext: ");
      for (std::uint64_t j = 0; j < text_length; ++j)
        fprintf(stdout, "%c", text[j]);
      fprintf(stdout, "\n");
    }
    fprintf(stdout, "\tmax threads = %lu\n", max_threads);
    fprintf(stdout, "\tcorrect gt begin: ");
    for (std::uint64_t i = 0; i < text_length; ++i)
      fprintf(stdout, "%lu", (std::uint64_t)correct_gt_begin.get(i));
    fprintf(stdout, "\n");
    fprintf(stdout, "\tcomputed gt begin: ");
    for (std::uint64_t i = 0; i < text_length; ++i)
      fprintf(stdout, "%lu", (std::uint64_t)computed_gt_begin.get(i));
    fprintf(stdout, "\n");
    std::fflush(stdout);
    std::exit(EXIT_FAILURE);
  }

  //---------------------------------------------------------------------------
  // STEP 4: clean up.
  //---------------------------------------------------------------------------
  delete[] computed_sa_temp;
}

template<unsigned pagesize_log>
void test_random(
    std::uint64_t testcases,
    std::uint64_t max_length,
    std::uint64_t max_sigma) {

  static const unsigned pagesize = (1U << pagesize_log);
  fprintf(stdout, "TEST, testcases = %4ld, pagesize = %5u, "
      "max_length = %8ld, max_sigma = %3ld\n",
      testcases, pagesize, max_length, max_sigma);
  std::fflush(stdout);
  std::uint8_t *text = new std::uint8_t[max_length];

  // Run tests.
  for (std::uint64_t tc = 0; tc < testcases; ++tc) {
    fprintf(stdout, "Progress: %.2Lf%%\r", (100.L * tc) / testcases);
    std::fflush(stdout);

    // Generate input.
    std::uint64_t length = utils::random_int64(1L, max_length);
    std::uint64_t max_threads = utils::random_int64(1L, 50L);

    std::uint64_t sigma;
    if (length <= 10000) sigma = utils::random_int64(1L, max_sigma);
    else sigma = utils::random_int64(2L, max_sigma);
    if (sigma <= 26) utils::fill_random_letters(text, length, sigma);
    else utils::fill_random_string(text, length, sigma);

    /*strcpy((char *)text, "ccbbba");
    std::uint64_t length = 6;
    std::uint64_t max_threads = 24;*/

    /*fprintf(stdout, "length = %lu\n", length);
    fprintf(stdout, "text: ");
    for (std::uint64_t j = 0; j < length; ++j)
      fprintf(stdout, "%c", text[j]);
    fprintf(stdout, "\n");
    fprintf(stdout, "max_threads = %lu\n", max_threads);*/

    test<pagesize_log>(text, length, max_threads);
  }

  // Clean up.
  delete[] text;
}


int main() {
  std::srand(std::time(0) + getpid());

  // Redirect stderr to /dev/null
  int redir = open("/dev/null", O_WRONLY);
  dup2(redir, 2);
  close(redir);

  test_random<7> (50, 1000, 5);
  test_random<12>(50, 1000, 5);
  test_random<7> (50, 1000, 20);
  test_random<12>(50, 1000, 20);
  test_random<7> (50, 1000, 128);
  test_random<12>(50, 1000, 128);
  test_random<7> (50, 1000, 255);
  test_random<12>(50, 1000, 255);
  
  test_random<7> (20, 10000, 5);
  test_random<12>(20, 10000, 5);
  test_random<7> (20, 10000, 20);
  test_random<12>(20, 10000, 20);
  test_random<7> (20, 10000, 128);
  test_random<12>(20, 10000, 128);
  test_random<7> (20, 10000, 255);
  test_random<12>(20, 10000, 255);

  test_random<7> (20, 1000000, 5);
  test_random<12>(20, 1000000, 5);
  test_random<7> (20, 1000000, 20);
  test_random<12>(20, 1000000, 20);
  test_random<7> (20, 1000000, 128);
  test_random<12>(20, 1000000, 128);
  test_random<7> (20, 1000000, 255);
  test_random<12>(20, 1000000, 255);

  fprintf(stdout, "All tests passed.\n");
}
