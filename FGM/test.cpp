#include <ctime>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <string>
#include <algorithm>

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "uint40.h"
#include "divsufsort.h"
#include "utils.h"
#include "fgm.h"

// Test many string chosen according to given paranters.
void test_random(int testcases, int max_length, int max_sigma) {
  printf("TEST, testcases = %d, max_n = %d, max_sigma = %d\n",
      testcases, max_length, max_sigma);
  fflush(stdout);
  unsigned char *text = new unsigned char[max_length + 1];
  int *SA = new int[max_length];

  for (int tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {
    // Print progress information.
    if (dbg == 10) {
      printf("%d (%.2Lf%%)\r", tc, (tc * 100.L) / testcases);
      fflush(stdout);
      dbg = 0;
    }

    // Generate string.
    int length = utils::random_int(1, max_length);
    int sigma = utils::random_int(1, max_sigma);
    int block_size = 0;
    do block_size = utils::random_int(1, length);
    while (20 * block_size <= length);
    if (max_sigma <= 26) utils::fill_random_letters(text, length, sigma);
    else utils::fill_random_string(text, length, sigma);
    // int length = strlen("cbcabccaac");
    // strcat((char *)text, "cbcabccaac");
    // int block_size = 6;
    text[length] = 0;
    std::string filename = "/tmp/in" + utils::random_string_hash();
    utils::write_text_to_file(text, length, filename);

    // Run the test on generated string.
    FGM(filename, block_size);
    
    // Compare the result to correct SA.
    std::reverse(text, text + length);
    divsufsort(text, SA, length); // recall that FGM computes the SA of *reversed* text.
    std::reverse(text, text + length);
    uint40 *computed_SA;
    utils::read_n_objects_from_file(computed_SA, length, filename + ".sa5");
    utils::file_delete(filename + ".sa5");
    bool eq = true;
    for (int i = 0; i < length; ++i)
      if ((unsigned long)SA[i] != computed_SA[i].ull()) { eq = false; break; }
    if (!eq) {
      printf("Error!\n");
      if (length < 1000) {
        printf("  text = %s\n", text);
        printf("  computed SA: ");
        for (int k = 0; k < length; ++k) printf("%lu ", (unsigned long)computed_SA[k]);
        printf("\n");
        printf("  correct SA:  ");
        for (int k = 0; k < length; ++k) printf("%d ", SA[k]);
        printf("\n");
      }
      fflush(stdout);
      std::exit(EXIT_FAILURE);
    }
    delete[] computed_SA;
    utils::file_delete(filename);
  }

  // Clean up.
  delete[] text;
  delete[] SA;
}

int main(int, char **) {
  srand(time(0) + getpid());

  // Redirect stderr >> /dev/null
  int redir = open("/dev/null", O_WRONLY);
  dup2(redir, 2);
  close(redir);

  printf("Testing FGM.\n");
  fflush(stdout);
  test_random(5000, 10,      5);
  test_random(5000, 10,     20);
  test_random(5000, 10,    128);
  test_random(5000, 10,    255);
  test_random(500, 100,      5);
  test_random(500, 100,     20);
  test_random(500, 100,    128);
  test_random(500, 100,    255);
  test_random(50, 1000,      5);
  test_random(50, 1000,     20);
  test_random(50, 1000,    128);
  test_random(50, 1000,    255);
  test_random(50, 10000,     5);
  test_random(50, 10000,    20);
  test_random(50, 10000,   128);
  test_random(50, 10000,   255);
  test_random(5, 100000,     5);
  test_random(5, 100000,    20);
  test_random(5, 100000,   128);
  test_random(5, 100000,   255);
  printf("All tests passed.\n");
  fflush(stdout);

}

