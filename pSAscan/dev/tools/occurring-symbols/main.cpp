#include <cstdio>
#include <cstdlib>

#include "utils.h"
#include "stream.h"


int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "%s FILE\n"
        "Display all bytes that occur in FILE.\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }

  long size = utils::file_size(argv[1]);
  stream_reader<unsigned char> *reader =
    new stream_reader<unsigned char>(argv[1], 2 << 20);

  long double start = utils::wclock();
  bool occurs[256] = {false};
  for (long i = 0, dbg = 0; i < size; ++i, ++dbg) {
    if (dbg == (64 << 20)) {
      long double elapsed = utils::wclock() - start;
      long double processed_MiB = i / (1024.L * 1024);
      long double speed = processed_MiB / elapsed;
      long sigma = 0;
      for (long j = 0; j < 256; ++j) if (occurs[j]) ++sigma;
      fprintf(stderr, "Processed %.0LfMiB (%.1Lf%%). Speed: %.2LfMiB/s. Current sigma: %ld\r",
          processed_MiB, (100.L * i) / size, speed, sigma);
      dbg = 0;
    }

    unsigned char c = reader->read();
    occurs[c] = true;
  }

  long sigma = 0;
  for (long i = 0; i < 256; ++i) if (occurs[i]) ++sigma;
  long double elapsed = utils::wclock() - start;
  long double processed_MiB = size / (1024.L * 1024);
  long double speed = processed_MiB / elapsed;
  fprintf(stderr, "Processed %.0LfMiB (100.0%%). Speed: %.2LfMiB/s. Computed alphabet size: %ld\n",
      processed_MiB, speed, sigma);
  fprintf(stderr, "Occurring symbols: ");
  for (int i = 0; i < 256; ++i)
    if (occurs[i]) fprintf(stderr, "%d ", i);
  fprintf(stderr, "\n");

  delete reader;
}

