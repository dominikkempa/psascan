#include <cstdio>
#include <cstdlib>

#include "utils.h"


int main() {
  long length = 2800L << 20; // 2800MiB
  long sigma = 4;
  unsigned char range_start = 0;
  long double *prob = new long double[sigma];
  prob[0] = 0.95;
  long double other_prob = (long double)(1 - prob[0]) / (long double)(sigma - 1);
  for (long i = 1; i < sigma; ++i)
    prob[i] = other_prob;

  fprintf(stderr, "length = %ld\n", length);
  fprintf(stderr, "sigma = %ld\n", sigma);

  for (long i = 0; i < length; ++i) {
    if (i % 1000000 == 0)
      fprintf(stderr, "progress: %.3Lf%%\r", (100.L * i) / length);
    long double p = (long double)rand() / RAND_MAX;

    unsigned char c = 0;
    long double sum = 0.L;
    while (c + 1 < sigma && sum + prob[c] <= p)
      sum += prob[c++];
    std::fputc(range_start + c, stdout);
  }
  fprintf(stderr, "\nFinished.\n");
  delete[] prob;
}

