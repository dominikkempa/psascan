#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <string>

#include "utils.h"
#include "new_rank.h"

// Test many string chosen according to given paranters.
void test_random(int length, int sigma) {
  unsigned char *text = new unsigned char[length + 1];
  // Generate string.
  fprintf(stderr, "Generating text: \n");
  long double start = utils::wclock();
  int freq_sigma = utils::random_int(1, sigma);
  int rare_sigma = sigma - freq_sigma;
  for (int j = 0; j < length; ++j) {
    if (j % 10000 == 0) fprintf(stderr, "\rprogress: %.2Lf%%", (100.L * j) / length);
    text[j] = utils::random_int(0, freq_sigma - 1);
  }
  fprintf(stderr, "\n");
  for (int j = 0; j < rare_sigma; ++j) {
    fprintf(stderr, "\rprogress: %.2Lf%%", (100.L * j) / rare_sigma);
    int tries = utils::random_int(1, (length + 90) / 100);
    for (int k = 0; k < tries; ++k)
      text[utils::random_int(0, length - 1)] = freq_sigma - 1 + j;
  }
  fprintf(stderr, "\nGenerated text in %.2Lf\n", utils::wclock() - start);

  for (long threads = 1; threads <= 32; threads *= 2) {
      long th = std::min(threads, 24L);
      fprintf(stderr, "\n===== threads = %ld =====\n", th);

      fprintf(stderr, "Building rank:\n");
      start = utils::wclock();
      rank4n *rank = new rank4n(text, length, threads);
      fprintf(stderr, "Built rank in %.2Lf\n", utils::wclock() - start);
      delete rank;
  }
  delete[] text;
}

int main(int, char **) {
  std::srand(std::time(0) + getpid());
  test_random(1L << 30, 256);
}
