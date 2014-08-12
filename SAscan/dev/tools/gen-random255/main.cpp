#include <cstdio>
#include <cstdlib>

#include "utils.h"


int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "%s <size-in-MiB> <alphabet-size>\n"
        "Generate random file of given size of alphabet"
        "0, 1, .., alphabet-size - 1\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }

  long length = atol(argv[1]) << 20;
  long sigma = atol(argv[2]);

  fprintf(stderr, "length = %ld\n", length);
  fprintf(stderr, "sigma = %ld\n", sigma);

  for (long i = 0; i < length; ++i) {
    if (i % (1 << 22) == 0)
      fprintf(stderr, "%.2Lf%%\r", (100.L * i) / length);
    unsigned char c = utils::random_int(0, sigma - 1);
    std::fputc(c, stdout);
  }
  fprintf(stderr, "Finished.\n");
}

