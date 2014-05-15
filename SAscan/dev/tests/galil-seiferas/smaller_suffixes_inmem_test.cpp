#include <cstdio>
#include <cstring>
#include <string>
#include <unistd.h>
#include <algorithm>

#include "srank.h"

inline int random_int(int p, int r) {
  return p + rand() % (r - p + 1);
}

int main() {
  std::srand(time(0) + getpid());

  fprintf(stderr, "smaller suffixes test\n");

  static const int max_n = 50;
  static const int max_sigma = 5;

  unsigned char *X = new unsigned char[max_n];
  unsigned char *Y = new unsigned char[max_n];
  unsigned char *A = new unsigned char[max_n];

  unsigned char *str = new unsigned char[3 * max_n];

  unsigned char *gt = new unsigned char[2 * max_n];

  // string we compose is XAY, where A can be empty.

  for (long tc = 0, dbg = 0; ; ++tc, ++dbg) {
    if (dbg == 100000) {
      fprintf(stderr, "tested = %ld, so far so good...\r", tc);
      dbg = 0;
    }

    int sigma = random_int(2, max_sigma);
    long Xlen = (long)random_int(1, max_n);
    long Alen = (long)random_int(0, max_n);
    long Ylen = (long)random_int(1, max_n);

    for (long i = 0; i < Xlen; ++i) X[i] = 'a' + random_int(0, sigma - 1);
    for (long i = 0; i < Ylen; ++i) Y[i] = 'a' + random_int(0, sigma - 1);
    for (long i = 0; i < Alen; ++i) A[i] = 'a' + random_int(0, sigma - 1);

    long total_len = 0;
    for (long i = 0; i < Xlen; ++i) str[total_len++] = X[i];
    for (long i = 0; i < Alen; ++i) str[total_len++] = A[i];
    for (long i = 0; i < Ylen; ++i) str[total_len++] = Y[i];

    //-------------------------------------------------------------------------
    // Compute gt bitvector of length Ylen where
    // gt[i] == 1 iff suffix of str of length
    // i + 1 is lexicographically larger than AY.
    //
    // NOTE: The indexing of such defined gt bitvector is right-to-left
    //       is reversed compared to the gt bitvector defined in the paper.
    //       The reversed order is more convenient for the implementation,
    //       since gt bitvectors are always computed and accessed in reverse.
    //-------------------------------------------------------------------------

    unsigned char *AY = str + Xlen;
    for (long i = 0; i < Ylen; ++i) {
      long suf_len = i + 1;
      unsigned char *suf = str + total_len - suf_len;
      long lcp = 0;
      while (lcp < suf_len && suf[lcp] == AY[lcp]) ++lcp;
      gt[i] = (lcp < suf_len && suf[lcp] > AY[lcp]);
    }

    // Compute answer using pattern matching.
    long computed_answer = string_range_matching::compute_smaller_suffixes(X, Xlen, Y, Ylen, gt);

    // Compute the correct answer.
    long correct_answer = 0;
    for (long i = 0; i < Xlen; ++i) {
      unsigned char *suf = str + i;
      long lcp = 0;
      while (lcp < Ylen && Y[lcp] == suf[lcp]) ++lcp;
      if (lcp < Ylen && suf[lcp] < Y[lcp]) ++correct_answer;
    }

    if (correct_answer != computed_answer) {
      std::string Xstr = std::string(X, X + (int)Xlen);
      std::string Astr = std::string(A, A + (int)Alen);
      std::string Ystr = std::string(Y, Y + (int)Ylen);

      fprintf(stderr, "Error:\n");
      fprintf(stderr, "\tXAY = %s ", Xstr.c_str());
      fprintf(stderr, "%s ", Astr.c_str());
      fprintf(stderr, "%s\n", Ystr.c_str());

      fprintf(stderr, "\tcomputed answer = %ld\n", computed_answer);
      fprintf(stderr, "\tcorrect answer = %ld\n", correct_answer);

      std::exit(EXIT_FAILURE);
    }
  }

  delete[] X;
  delete[] A;
  delete[] Y;
  delete[] str;
}

