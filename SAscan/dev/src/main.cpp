#include <cstdio>
#include <cstdlib>

#include "sascan.h"

extern long max_threads;

int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "usage: %s <file> <mem>\n\n"
                    "Compute SA of <file> using <mem> MiB of RAM and write to <file>.sa5\n",
                    argv[0]);
    std::exit(EXIT_FAILURE);
  }

  max_threads = 32L;
  // NOTE: the number of threads can (?) be obtained using STL method:
  // http://en.cppreference.com/w/cpp/thread/thread/hardware_concurrency

  long ram_use = (long)atoi(argv[2]) << 20;
  // ram_use -= (11L << 20) * max_threads; // subtract the RAM for threads.
                                  // Disabled for now for testing purposes.
  SAscan(argv[1], ram_use);
}
