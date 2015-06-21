#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <string.h>
#include <errno.h>
#include <stdint.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <algorithm>

#include "nocache_async_stream_reader.h"
#include "utils.h"

#define NOCACHE_ALIGNMENT 4096


int main(int argc, char **argv) {
  if (argc != 2)
    std::exit(EXIT_FAILURE);

  // Experiment:
  //   1. perform some reading and display speed
  //   2. close
  //   3. open again
  //   4. perform some reading and display speed
  //
  // With O_DIRECT the speed should be the same in both cases.

  typedef nocache_async_stream_reader<unsigned char, NOCACHE_ALIGNMENT> reader_type;
  long length = utils::file_size(argv[1]);

  // 1. Read file.
  {
    reader_type *reader = new reader_type(argv[1]);
    long double start = utils::wclock();
    long sum = 0;
    for (long j = 0, dbg = 0; j < length; ++j, ++dbg) {
      sum += reader->read();
      if (dbg > (100L << 20)) {
        long double elapsed = utils::wclock() - start;
        long double io_speed = (1.L * j / (1 << 20)) / elapsed;
        fprintf(stderr, "\r%.2Lf%%, I/O: %.2LfMiB/s", (100.L * j) / length, io_speed);
        dbg = 0;
      }
    }
    delete reader;
    long double elapsed = utils::wclock() - start;
    long double io_speed = (1.L * length / (1 << 20)) / elapsed;
    fprintf(stderr, "\r100.00%%, I/O: %.2LfMiB/s\n", io_speed);
    fprintf(stderr, "sum = %ld\n", sum);
  }

  // 2. Read file again.
  {
    reader_type *reader = new reader_type(argv[1]);
    long double start = utils::wclock();
    long sum = 0;
    for (long j = 0, dbg = 0; j < length; ++j, ++dbg) {
      sum += reader->read();
      if (dbg > (100L << 20)) {
        long double elapsed = utils::wclock() - start;
        long double io_speed = (1.L * j / (1 << 20)) / elapsed;
        fprintf(stderr, "\r%.2Lf%%, I/O: %.2LfMiB/s", (100.L * j) / length, io_speed);
        dbg = 0;
      }
    }
    delete reader;
    long double elapsed = utils::wclock() - start;
    long double io_speed = (1.L * length / (1 << 20)) / elapsed;
    fprintf(stderr, "\r100.00%%, I/O: %.2LfMiB/s\n", io_speed);
    fprintf(stderr, "sum = %ld\n", sum);
  }
}
