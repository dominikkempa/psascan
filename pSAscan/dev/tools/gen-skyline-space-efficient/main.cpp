#include <cstdio>
#include <cstdlib>

#include "async_stream_writer.hpp"


int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "%s ORDER\n"
      "Write Skyline word of given ORDER (and length 2^ORDER-1) on standard output\n",
      argv[0]);
    std::exit(EXIT_FAILURE);
  }

  std::uint64_t order = atol(argv[1]);

  {
    typedef async_stream_writer<std::uint8_t> writer_type;
    writer_type *writer = new writer_type();
    for (std::uint64_t i = 1; i < (1UL << order); ++i) {
      if (i % 10000000 == 0)
        fprintf(stderr, "%.2Lf%%\r", (100.L * i) / (1UL << order));
      writer->write((std::uint8_t)__builtin_ctzll(i));
    }
    delete writer;
  }

  fprintf(stderr, "Finished.\n");
}

