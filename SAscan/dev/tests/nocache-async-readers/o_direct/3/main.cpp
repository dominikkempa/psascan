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

#include "aligned_alloc.h"
#include "utils.h"


int open_file_direct(std::string filename, int mode) {
  int fd = open(filename.c_str(), mode);
  if (fd == -1) {
    fprintf(stderr, "\n open() failed with error [%s]\n", strerror(errno));
    std::exit(EXIT_FAILURE);
  }

  return fd;
}

void drop_pages(std::string filename) {
  long double start = utils::wclock();
  int fd = open_file_direct(filename.c_str(), O_RDWR);
  long length = lseek(fd, 0L, SEEK_END);
  lseek(fd, 0L, SEEK_SET);
  posix_fadvise(fd, 0, length, POSIX_FADV_DONTNEED);
  close(fd);
  long double elapsed = utils::wclock() - start;
  fprintf(stderr, "drop-cache: %.4Lfs\n", elapsed);
}

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

  static const long bufsize = (1L << 20);
  unsigned char *buffer = (unsigned char *)aligned_alloc<DIRECT_ALIGNMENT>(bufsize);

  // 1. Read file.
  {
    int fd = open_file_direct(argv[1], O_RDONLY);
    long length = lseek(fd, 0L, SEEK_END);
    lseek(fd, 0L, SEEK_SET);
    long double start = utils::wclock();
    long read_ret;
    long total_read = 0L;
    long dbg = 0L;
    while((read_ret = read(fd, buffer, bufsize)) > 0) {
      if (dbg > (100L << 20)) {
        long double elapsed = utils::wclock() - start;
        long double io_speed = (1.L * total_read / (1 << 20)) / elapsed;
        fprintf(stderr, "\r%.2Lf%%, I/O: %.2LfMiB/s", (100.L * total_read) / length, io_speed);
        dbg = 0L;
      }
      total_read += read_ret;
      dbg += read_ret;
    }
    long double elapsed = utils::wclock() - start;
    long double io_speed = (1.L * total_read / (1 << 20)) / elapsed;
    fprintf(stderr, "\r100.00%%, I/O: %.2LfMiB/s\n", io_speed);
    fprintf(stderr, "Read %ld bytes.\n", total_read);
    close(fd);
  }

  // 2. Read file again.
  drop_pages(argv[1]);
  {
    int fd = open_file_direct(argv[1], O_RDONLY);
    long length = lseek(fd, 0L, SEEK_END);
    lseek(fd, 0L, SEEK_SET);
    long double start = utils::wclock();
    long read_ret;
    long total_read = 0L;
    long dbg = 0L;
    while((read_ret = read(fd, buffer, bufsize)) > 0) {
      if (dbg > (100L << 20)) {
        long double elapsed = utils::wclock() - start;
        long double io_speed = (1.L * total_read / (1 << 20)) / elapsed;
        fprintf(stderr, "\r%.2Lf%%, I/O: %.2LfMiB/s", (100.L * total_read) / length, io_speed);
        dbg = 0L;
      }
      total_read += read_ret;
      dbg += read_ret;
    }
    long double elapsed = utils::wclock() - start;
    long double io_speed = (1.L * total_read / (1 << 20)) / elapsed;
    fprintf(stderr, "\r100.00%%, I/O: %.2LfMiB/s\n", io_speed);
    fprintf(stderr, "Read %ld bytes.\n", total_read);
    close(fd);
  }
  drop_pages(argv[1]);

  aligned_dealloc(buffer);
}
