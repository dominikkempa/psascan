#include <cstdio>
#include <cstring>
#include <string>
#include <unistd.h>
#include <algorithm>

#include "smaller_suffixes.h"
#include "utils.h"
#include "stream.h"

inline int random_int(int p, int r) {
  return p + rand() % (r - p + 1);
}

int main() {
  std::srand(time(0) + getpid());

  fprintf(stderr, "smaller suffixes test\n");

  static const int max_n = 20;
  static const int max_sigma = 5;

  unsigned char *X = new unsigned char[max_n];
  unsigned char *Y = new unsigned char[max_n];
  unsigned char *str = new unsigned char[2 * max_n];

  // string we compose is str = XY, where both X and Y are nonempty.

  for (long tc = 0, dbg = 0; ; ++tc, ++dbg) {
    if (dbg == 10) {
      fprintf(stderr, "tested = %ld, so far so good...\r", tc);
      dbg = 0;
    }

    int sigma = random_int(2, max_sigma);
    long Xlen = (long)random_int(1, max_n);
    long Ylen = (long)random_int(1, max_n);

    for (long i = 0; i < Xlen; ++i) X[i] = 'a' + random_int(0, sigma - 1);
    for (long i = 0; i < Ylen; ++i) Y[i] = 'a' + random_int(0, sigma - 1);

    long total_len = 0;
    for (long i = 0; i < Xlen; ++i) str[total_len++] = X[i];
    for (long i = 0; i < Ylen; ++i) str[total_len++] = Y[i];

    // Write str to file.
    utils::write_objects_to_file<unsigned char>(str, total_len, "text");

    long pat_start = Xlen + utils::random_long(0L, Ylen - 1);
    unsigned char *pat = str + pat_start;
    long pat_length = total_len - pat_start;

    //--------------------------------------------------------------------------
    // Compute and write gt bitvector to file. gt is a bitvector of length Ylen
    // where gt[i] == 1 iff suffix of str of length i + 1 is lexicographically
    // larger than Y.
    //
    // NOTE: The indexing of such defined gt bitvector is right-to-left
    //       is reversed compared to the gt bitvector defined in the paper.
    //       The reversed order is more convenient for the implementation,
    //       since gt bitvectors are always computed and accessed in reverse.
    //--------------------------------------------------------------------------
    bit_stream_writer *gt = new bit_stream_writer("text.gt_tail");
    for (long i = 0; i < Ylen; ++i) {
      long suf_len = i + 1;
      unsigned char *suf = str + total_len - suf_len;
      long lcp = 0;
      while (lcp < suf_len && suf[lcp] == Y[lcp]) ++lcp;
      gt->write(lcp < suf_len && suf[lcp] > Y[lcp]);
    }
    delete gt;

    // Compute the correct answer.
    long correct_answer = 0L;
    for (long i = 0; i < Xlen; ++i) {
      unsigned char *suf = str + i;      
      long lcp = 0;
      while (lcp < pat_length && pat[lcp] == suf[lcp]) ++lcp;
      if (lcp < pat_length && suf[lcp] < pat[lcp]) ++correct_answer;
    }

    // Compute the answer using string range matching.
    long computed_answer = parallel_smaller_suffixes(X, Xlen, "text", pat_start);

    // Compate answers.
    if (computed_answer != correct_answer) {
      std::string S = std::string(str, str + total_len);
      fprintf(stderr, "Error:\n");
      fprintf(stderr, "\ttext = %s\n", S.c_str());
      fprintf(stderr, "\tpat_start = %ld\n", pat_start);
      fprintf(stderr, "\tcomputed answer = %ld\n", computed_answer);
      fprintf(stderr, "\tcorrect  answer = %ld\n", correct_answer);
      
      std::exit(EXIT_FAILURE);
    }

    // Clean up.
    utils::file_delete("text");
    utils::file_delete("text.gt_tail");
  }

  delete[] X;
  delete[] Y;
  delete[] str;
}

