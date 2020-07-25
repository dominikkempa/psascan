#include <cstdio>
#include <algorithm>
#include <map>

#include "utils.h"
#include "uint40.h"
#include "stream.h"
#include "disk_pattern.h"


inline bool compare(pattern &text, long length, long x, long y) {
  long lcp = 0;
  while (x + lcp < length && y + lcp < length && text[x + lcp] == text[y + lcp]) ++lcp;
  return (x + lcp == length) || (y + lcp < length && text[x + lcp] < text[y + lcp]);
}

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "usage: %s <file>\n"
                    "Check the correctness of <file>.sa5\n",
                    argv[0]);
    std::exit(EXIT_FAILURE);
  }
  
  long length = utils::file_size(argv[1]);
  pattern text(argv[1], 0L);

  fprintf(stderr, "Filename = %s\n", argv[1]);
  fprintf(stderr, "Length = %ld\n", length);

  std::string sa_fname = std::string(argv[1]) + ".sa5";
  std::FILE *f_sa = utils::open_file(sa_fname.c_str(), "r");
  std::map<long, long> sufs;
  for (long offset = 0; offset < length; offset += (1024L << 20)) {
    fprintf(stderr, "Offset = %ld (%.2Lf%%)\n", offset, (100.L * offset) / length);

    long prev;
    if (offset == 0) prev = length;
    else prev = (long)(utils::read_object_at_offset<uint40>(f_sa, offset - 1));

    sufs[prev] += 1;
    if (sufs[prev] > 1) {
      fprintf(stderr, "\nDetected duplicate suffix!\n");
      std::exit(EXIT_FAILURE);
    }

    long totest = std::min(length - offset, 100L);
    for (long j = 0; j < totest; ++j) {
      long cur = (long)(utils::read_object_at_offset<uint40>(f_sa, offset + j));
      fprintf(stderr, "\r  Testing suffix %ld vs %ld (%ld / %ld)                 ",
          prev, cur, j, totest);
      sufs[cur] += 1;
      if (sufs[cur] > 1) {
        fprintf(stderr, "\nDetected duplicate suffix!\n");
        std::exit(EXIT_FAILURE);
      }
      if (!compare(text, length, prev, cur)) {
        fprintf(stderr, "\nError.\n");
        std::exit(EXIT_FAILURE);
      }
      prev = cur;
    }
    fprintf(stderr, "\r  All tested suffixes were correct.                                          \n");
  }

  std::fclose(f_sa);

  fprintf(stderr, "\nNo errors have been found.\n");
}

