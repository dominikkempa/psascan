#include <cstdio>
#include <cstdlib>

#include "utils.h"
#include "stream.h"

int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "%s FILE SIZE\n"
        "Read SIZE MiB of FILE into RAM and measure reading speed.\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }

  long max_size = atol(argv[2]) << 20;
  long size = utils::file_size(argv[1]);
  long toread = std::min(size, max_size);

  fprintf(stderr, "Alloc and clear memory: ");
  unsigned char *block = new unsigned char[toread];
  std::fill(block, block + toread, 0);
  fprintf(stderr, "DONE\n");

  fprintf(stderr, "Reading: ");
  long double start = utils::wclock();
  utils::read_block(argv[1], 0, toread, block);
  long double total_time = utils::wclock() - start;
  long double speed = (toread / (1024.L * 1024)) / total_time;
  fprintf(stderr, "%.2Lfs (%.2LfMiB/s)\n", total_time, speed);
  delete[] block;
}

