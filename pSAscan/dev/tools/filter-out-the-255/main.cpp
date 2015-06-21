#include <cstdio>
#include <cstdlib>
#include <string>
#include <algorithm>

#include "utils.h"
#include "stream.h"


int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "%s FILE\n"
        "Display all bytes that occur in FILE.\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }

  std::string text_filename = argv[1];
  std::string out_filename = text_filename + ".filetered255";

  long size = utils::file_size(argv[1]);
  stream_reader<unsigned char> *reader = new stream_reader<unsigned char>(text_filename, 2 << 20);
  stream_writer<unsigned char> *writer = new stream_writer<unsigned char>(out_filename, 2 << 20);

  long double start = utils::wclock();
  for (long i = 0, dbg = 0; i < size; ++i, ++dbg) {
    if (dbg == (64 << 20)) {
      long double elapsed = utils::wclock() - start;
      long double processed_MiB = (2L * i) / (1024.L * 1024);
      long double speed = processed_MiB / elapsed;
      fprintf(stderr, "Processed %.0LfMiB (%.1Lf%%). Speed: %.2LfMiB/s.\r", processed_MiB, (100.L * i) / size, speed);
      dbg = 0;
    }

    unsigned char c = reader->read();
    if (c != 255) writer->write(c);
  }

  long double elapsed = utils::wclock() - start;
  long double processed_MiB = (2L * size) / (1024.L * 1024);
  long double speed = processed_MiB / elapsed;
  fprintf(stderr, "Processed %.0LfMiB (100.0%%). Speed: %.2LfMiB/s.\n", processed_MiB, speed);

  delete reader;
  delete writer;
}

