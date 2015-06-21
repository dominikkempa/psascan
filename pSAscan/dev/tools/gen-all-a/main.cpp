#include <cstdio>
#include <cstdlib>

#include "utils.h"


int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "%s <size-in-MiB>\n"
        "Generate string aaaa.. of given size\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }

  long length = atol(argv[1]) << 20;
  fprintf(stderr, "length = %ld\n", length);

  for (long i = 0; i < length; ++i) {
    if (i % (1 << 22) == 0)
      fprintf(stderr, "%.2Lf%%\r", (100.L * i) / length);
    unsigned char c = 'a';
    std::fputc(c, stdout);
  }
  fprintf(stderr, "Finished.\n");
}

