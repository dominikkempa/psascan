#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <ctime>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "inmem_psascan.h"
#include "divsufsort.h"
#include "divsufsort64.h"
#include "utils.h"


template<unsigned pagesize_log>
void test(unsigned char *text, long text_length, std::uint64_t max_threads) {
  //----------------------------------------------------------------------------
  // STEP 1: compute correct answer.
  //----------------------------------------------------------------------------
  int *correct_sa = new int[text_length];
  divsufsort(text, correct_sa, (long)text_length);


  //----------------------------------------------------------------------------
  // STEP 2: run the tested algorithm.
  //----------------------------------------------------------------------------
  // 1) compute gt_in of size block_size, where gt_in[i] == 1 iff
  //    text[block_beg + i..) > text[block_end..).
  unsigned char *computed_sa_temp = new unsigned char[text_length * (sizeof(int) + 1)];
  int *computed_sa = (int *)computed_sa_temp;
  inmem_psascan<int, pagesize_log>(text, text_length, computed_sa_temp, max_threads);

  //----------------------------------------------------------------------------
  // STEP 3: compare answers.
  //----------------------------------------------------------------------------
  if (!std::equal(correct_sa, correct_sa + text_length, computed_sa)) {
    fprintf(stdout, "\nError:\n");
    fprintf(stdout, "\tlength = %ld\n", text_length);
    fprintf(stdout, "\tmax threads = %lu\n", max_threads);
    if (text_length <= 1000) {
      fprintf(stdout, "\ttext: ");
      for (long j = 0; j < text_length; ++j)
        fprintf(stdout, "%c", text[j]);
      fprintf(stdout, "\n");
      fprintf(stdout, "\tcorrect sa: ");
      for (long i = 0; i < text_length; ++i)
        fprintf(stdout, "%d ", correct_sa[i]);
      fprintf(stdout, "\n");
      fprintf(stdout, "\tcomputed sa: ");
      for (long i = 0; i < text_length; ++i)
        fprintf(stdout, "%d ", computed_sa[i]);
      fprintf(stdout, "\n");
    }
    std::fflush(stdout);
    std::exit(EXIT_FAILURE);
  }

  //----------------------------------------------------------------------------
  // STEP 4: clean up.
  //----------------------------------------------------------------------------
  delete[] correct_sa;
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
    std::uint64_t max_threads = utils::random_long(1L, 50L);

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

  test_random<1> (100, 300, 5);
  test_random<1> (100, 300, 128);
  test_random<3> (100, 300, 5);
  test_random<3> (100, 300, 128);
  test_random<5> (100, 300, 5);
  test_random<5> (100, 300, 128);
  test_random<12>(100, 300, 5);
  test_random<12>(100, 300, 128);
  test_random<1> (100, 300, 255);
  test_random<3> (100, 300, 255);
  test_random<5> (100, 300, 255);
  test_random<12>(100, 300, 255);

  test_random<1> (100, 1000, 5);
  test_random<1> (100, 1000, 128);
  test_random<3> (100, 1000, 5);
  test_random<3> (100, 1000, 128);
  test_random<5> (100, 1000, 5);
  test_random<5> (100, 1000, 128);
  test_random<12>(100, 1000, 5);
  test_random<12>(100, 1000, 128);
  test_random<1> (100, 1000, 255);
  test_random<3> (100, 1000, 255);
  test_random<5> (100, 1000, 255);
  test_random<12>(100, 1000, 255);

  test_random<7> (1000, 1000,    5);
  test_random<12>(1000, 1000,    5);
  test_random<7> (1000, 1000,    20);
  test_random<12>(1000, 1000,    20);
  test_random<7> (1000, 1000,    128);
  test_random<12>(1000, 1000,    128);
  test_random<7> (1000, 1000,    255);
  test_random<12>(1000, 1000,    255);
  test_random<7> (500,  10000,   5);
  test_random<12>(500,  10000,   5);
  test_random<7> (500,  10000,   20);
  test_random<12>(500,  10000,   20);
  test_random<7> (500,  10000,   128);
  test_random<12>(500,  10000,   128);
  test_random<7> (500,  10000,   255);
  test_random<12>(500,  10000,   255);
  test_random<7> (100,  1000000, 5);
  test_random<12>(100,  1000000, 5);
  test_random<7> (100,  1000000, 20);
  test_random<12>(100,  1000000, 20);
  test_random<7> (100,  1000000, 128);
  test_random<12>(100,  1000000, 128);
  test_random<7> (100,  1000000, 255);
  test_random<12>(100,  1000000, 255);

  fprintf(stdout, "All tests passed.\n");
}
