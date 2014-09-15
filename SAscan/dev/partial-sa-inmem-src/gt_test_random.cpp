#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "utils.h"
#include "inmem_sascan.h"
#include "bitvector.h"

template<typename T>
void next(unsigned char *text, T length, T &s, T &p, T &r) {
  if (length == 1) { s = 0; p = 1; r = 0; return; }
  T i = length - 1;
  while (i < length) {
    unsigned char a = text[s + r], b = text[i];
    if (a > b) { p = i - s + 1; r = 0; }
    else if (a < b) { i -= r; s = i; p = 1; r = 0; }
    else { ++r; if (r == p) r = 0; } ++i;
  }
}


//==============================================================================
// Compute gt_begin for text.
//==============================================================================
long compute_gt_begin(unsigned char *text, long length, bitvector *gt_begin) {
  long whole_suffix_rank = length - 1;
  long i = 1, el = 0, s = 0, p = 0, r = 0;
  long i_max = 0, el_max = 0, s_max = 0, p_max = 0, r_max = 0;
  while (i < length) {
    while (i + el < length && el < length && text[i + el] == text[el])
      next(text, ++el, s, p, r);
 
    if (i + el < length && (el == length || text[i + el] > text[el])) {
      gt_begin->set(i - 1);
      --whole_suffix_rank;
    }

    long j = i_max;
    if (el > el_max) {
      std::swap(el, el_max);
      std::swap(s, s_max);
      std::swap(p, p_max);
      std::swap(r, r_max);
      i_max = i;
    }

    if (p && 3 * p <= el && !memcmp(text, text + p, s)) {
      for (long k = 1; k < std::min(length - i, p); ++k) { 
        if (gt_begin->get(j + k - 1)) {
          gt_begin->set(i + k - 1);
          --whole_suffix_rank;
        }
      }

      i += p;
      el -= p;
    } else {
      long h = (el / 3) + 1;
      for (long k = 1; k < std::min(length - i, h); ++k) { 
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
void test(unsigned char *text, long text_length, long max_threads) {
  //----------------------------------------------------------------------------
  // STEP 1: compute correct answer.
  //----------------------------------------------------------------------------
  bitvector correct_gt_begin(text_length, max_threads);
  compute_gt_begin(text, text_length, &correct_gt_begin);


  //----------------------------------------------------------------------------
  // STEP 2: run the tested algorithm.
  //----------------------------------------------------------------------------
  // 1) compute gt_in of size block_size, where gt_in[i] == 1 iff
  //    text[block_beg + i..) > text[block_end..).
  unsigned char *computed_sa_temp = new unsigned char[text_length * (sizeof(int) + 1)];
  bitvector computed_gt_begin(text_length, max_threads);
  inmem_sascan<int, pagesize_log>(text, text_length, computed_sa_temp, max_threads, false, true, &computed_gt_begin);

  //----------------------------------------------------------------------------
  // STEP 3: compare answers.
  //----------------------------------------------------------------------------
  bool eq = true;
  for (long i = 0; i < text_length; ++i) {
    if (correct_gt_begin.get(i) != computed_gt_begin.get(i)) { eq = false; break; }
  }
  
  if (!eq) {
    fprintf(stdout, "\nError:\n");
    fprintf(stdout, "\tlength = %ld\n", text_length);
    if (text_length <= 1000) {
      fprintf(stdout, "\ttext: ");
      for (long j = 0; j < text_length; ++j)
        fprintf(stdout, "%c", text[j]);
      fprintf(stdout, "\n");
    }
    fprintf(stdout, "\tmax threads = %ld\n", max_threads);
    fprintf(stdout, "\tcorrect gt begin: ");
    for (long i = 0; i < text_length; ++i)
      fprintf(stdout, "%ld", (long)correct_gt_begin.get(i));
    fprintf(stdout, "\n");
    fprintf(stdout, "\tcomputed gt begin: ");
    for (long i = 0; i < text_length; ++i)
      fprintf(stdout, "%ld", (long)computed_gt_begin.get(i));
    fprintf(stdout, "\n");
    std::fflush(stdout);
    std::exit(EXIT_FAILURE);
  }

  //----------------------------------------------------------------------------
  // STEP 4: clean up.
  //----------------------------------------------------------------------------
  delete[] computed_sa_temp;
}

template<unsigned pagesize_log>
void test_random(long testcases, long max_length, long max_sigma) {
  static const unsigned pagesize = (1U << pagesize_log);
  fprintf(stdout, "TEST, testcases = %4ld, pagesize = %5u, max_length = %8ld, max_sigma = %3ld\n",
      testcases, pagesize, max_length, max_sigma);
  std::fflush(stdout);
  unsigned char *text = new unsigned char[max_length];

  // Run tests.
  for (long tc = 0; tc < testcases; ++tc) {
    fprintf(stdout, "Progress: %.2Lf%%\r", (100.L * tc) / testcases);
    std::fflush(stdout);

    // Generate input.
    long length = utils::random_long(1L, max_length);
    long max_threads = utils::random_long(1L, 50L);

    long sigma;
    if (length <= 10000) sigma = utils::random_long(1L, max_sigma);
    else sigma = utils::random_long(2L, max_sigma);
    if (sigma <= 26) utils::fill_random_letters(text, length, sigma);
    else utils::fill_random_string(text, length, sigma);

    /*strcpy((char *)text, "ccbbba");
    long length = 6;
    long max_threads = 24;*/

    /*fprintf(stdout, "length = %ld\n", length);
    fprintf(stdout, "text: ");
    for (long j = 0; j < length; ++j)
      fprintf(stdout, "%c", text[j]);
    fprintf(stdout, "\n");
    fprintf(stdout, "max_threads = %ld\n", max_threads);*/

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
