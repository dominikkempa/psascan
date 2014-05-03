#include <cstdio>
#include <cstring>
#include <string>
#include <unistd.h>
#include <algorithm>

#include "srank.h"

using namespace std;

inline int random_int(int p, int r) {
  return p + rand() % (r - p + 1);
}

int main() {
  srand(time(0) + getpid());

  fprintf(stderr, "strrank test\n");

  // int testcases = (1 << 25);
  static const int max_n = 500;
  static const int max_sigma = 5;

  unsigned char *X = new unsigned char[max_n];
  unsigned char *Y = new unsigned char[max_n];

  for (int tc = 0; /*tc < testcases*/; ++tc) {
    if ((tc % 1000) == 0)
      fprintf(stderr, "tested = %d, so far so good...\r", tc);

    int sigma = random_int(2, max_sigma);

    int n = random_int(2, max_n);
    for (int i = 0; i < n; ++i)
      X[i] = 'a' + random_int(0, sigma - 1);

    int m = random_int(2, max_n);
    for (int i = 0; i < m; ++i)
      Y[i] = 'a' + random_int(0, sigma - 1);
   
    long computed_rank = (int)string_range_matching::strrank(X, (long)n, Y, (long)m);
    long computed_rank2 = (int)string_range_matching::strrank_ondemand(X, (long)n, Y, (long)m);

    std::string Xstr = std::string(X, X + n);
    std::string Ystr = std::string(Y, Y + m);

    long correct_rank = 0;
    for (int i = 0; i < n; ++i) {
      std::string Xi = Xstr.substr(i);
      if (Xi < Ystr) ++correct_rank;
    }

    if (correct_rank != computed_rank || computed_rank != computed_rank2) {
      fprintf(stderr, "\nError:\n");
      fprintf(stderr, "\tX = %s\n", Xstr.c_str());
      fprintf(stderr, "\tY = %s\n", Ystr.c_str());
      fprintf(stderr, "correct rank = %ld\n", correct_rank);
      fprintf(stderr, "computed rank = %ld\n", computed_rank);
      fprintf(stderr, "computed_rank2 = %ld\n", computed_rank2);

      std::exit(EXIT_FAILURE);
    }
  }

  delete[] Y;
  delete[] X;
}

