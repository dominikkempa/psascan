#include <cstdio>
#include <cstdlib>

#include "fgm.h"

int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "usage: %s <file> <mem>\n"
                    "Compute SA of <file> using <mem> MiB of RAM and write to <file>.sa.\n",
                    argv[0]);
    std::exit(EXIT_FAILURE);
  }

  long ram_use = (long)atoi(argv[2]) << 20;
  long block_size = ram_use / 5;

  fprintf(stderr, "Filename = %s\n", argv[1]);
  fprintf(stderr, "RAM use = %ld\n", ram_use);
  fprintf(stderr, "block size = %ld\n", block_size);

  FGM(argv[1], block_size);
}

