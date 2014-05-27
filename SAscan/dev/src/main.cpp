#include <cstdio>
#include <cstdlib>

#include "sascan.h"

extern long max_threads;
extern long n_updaters;
extern long stream_buffer_size;
extern long n_stream_buffers;
extern long max_gap_sections;

int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "usage: %s <file> <mem>\n\n"
                    "Compute SA of <file> using <mem> MiB of RAM and write to <file>.sa5\n",
                    argv[0]);
    std::exit(EXIT_FAILURE);
  }

  // NOTE: the number of threads can (?) be obtained using STL method:
  // http://en.cppreference.com/w/cpp/thread/thread/hardware_concurrency
  max_threads = 24;
  n_updaters = 24;
  stream_buffer_size = (2 << 20);
  n_stream_buffers = 48;
  max_gap_sections = 6;

  long ram_use = (long)atoi(argv[2]) << 20;
  ram_use -= stream_buffer_size * n_stream_buffers; // RAM for buffers

  SAscan(argv[1], ram_use);
}
