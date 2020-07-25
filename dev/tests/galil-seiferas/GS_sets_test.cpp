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

  fprintf(stderr, "GS_sets test\n");

  static const int max_n = 50;
  static const int max_sigma = 5;

  unsigned char *X = new unsigned char[max_n];

  for (long tc = 0, dbg = 0; ; ++tc, ++dbg) {
    if (dbg == 100000) {
      fprintf(stderr, "tested = %ld, so far so good...\r", tc);
      dbg = 0;
    }

    int sigma = random_int(2, max_sigma);

    long n = (long)random_int(2, max_n);
    for (long i = 0; i < n; ++i)
      X[i] = 'a' + random_int(0, sigma - 1);

    long pat_length = 0L;
    string_range_matching::GS_sets GS;
    while (pat_length < n) {
      long old_pat_length = pat_length;
      int increase = random_int(1, n - pat_length);
      pat_length += (long)increase;

      // Update sets using GS_sets.
      GS.extend_pattern(X, pat_length);

      // Compute corret answer.
      std::pair<std::vector<string_range_matching::triple>,
                std::vector<string_range_matching::pair> > S
        = string_range_matching::precompute(X, pat_length);

      if (S.first != GS.S_p || S.second != GS.S_n) {
        std::string Xstr = std::string(X, X + (int)n);
        fprintf(stderr, "Error:\n");
        fprintf(stderr, "\tX = %s\n", Xstr.c_str());
        fprintf(stderr, "\tn = %ld\n", n);
        fprintf(stderr, "\told_pat_length = %ld\n", old_pat_length);
        fprintf(stderr, "\tnew_pat_length = %ld\n\n", pat_length);
        fprintf(stderr, "\t[correct] S_p = {");
        for (long i = 0; i < (long)S.first.size(); ++i)
          fprintf(stderr, "(%ld,%ld,%ld) ", S.first[i].b, S.first[i].e, S.first[i].c);
        fprintf(stderr, "}\n");
        fprintf(stderr, "\t[computed]S_p = {");
        for (long i = 0; i < (long)GS.S_p.size(); ++i)
          fprintf(stderr, "(%ld,%ld,%ld) ", GS.S_p[i].b, GS.S_p[i].e, GS.S_p[i].c);
        fprintf(stderr, "}\n\n");
        fprintf(stderr, "\t[correct] S_n = {");
        for (long i = 0; i < (long)S.second.size(); ++i)
          fprintf(stderr, "(%ld,%ld) ", S.second[i].b, S.second[i].c);
        fprintf(stderr, "}\n");
        fprintf(stderr, "\t[computed]S_n = {");
        for (long i = 0; i < (long)GS.S_n.size(); ++i)
          fprintf(stderr, "(%ld,%ld) ", GS.S_n[i].b, GS.S_n[i].c);
        fprintf(stderr, "}\n");
        std::exit(EXIT_FAILURE);
      }
    }
  }

  delete[] X;
}

