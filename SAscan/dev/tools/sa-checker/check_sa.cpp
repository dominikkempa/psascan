#include <cstdio>
#include <algorithm>

#include "utils.h"
#include "uint40.h"
#include "stream.h"

inline bool compare(unsigned char *text, long length, long x, long y) {
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
  
  long length;
  unsigned char *text;
  utils::read_file(text, length, argv[1]);
  fprintf(stderr, "Filename = %s\n", argv[1]);
  fprintf(stderr, "Length = %ld\n", length);

  stream_reader<uint40> *sr = new stream_reader<uint40>(std::string(argv[1]) + ".sa5", 1 << 20);
  bool ok = true;
  int wrong_suffixes = 0;
  uint40 prev = sr->read();
  for (long i = 1; i < length; ++i) {
    if (i % 1000000 == 0) fprintf(stderr, "\rTesting: %ld", i);
    uint40 cur = sr->read();
    if (!compare(text, length, prev.ull(), cur.ull())) {
      ok = false;
      ++wrong_suffixes;
      fprintf(stderr, "Error.\n");
      std::exit(EXIT_FAILURE);
    }
    prev = cur;
  }
  fprintf(stderr, "\n");

  if (!ok) {
    fprintf(stderr, "Does not match!\n");
    fprintf(stderr, "wrong_suffixes = %d\n", wrong_suffixes);
  } else fprintf(stderr, "OK\n");
  delete sr;
}

