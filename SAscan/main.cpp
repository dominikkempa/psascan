#include <cstdio>
#include <cstdlib>

#include "sascan.h"

int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "usage: %s <file> <mem>\n\n"
                    "Compute SA of <file> using <mem> MiB of RAM and write to <file>.sa5\n",
                    argv[0]);
    std::exit(EXIT_FAILURE);
  }

  long ram_use = (long)atoi(argv[2]) << 20;
  SAscan(argv[1], ram_use);
}
