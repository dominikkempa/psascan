#include <cstdio>
#include <cstdlib>

#include "utils.h"


int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "%s <size-in-MiB> <alphabet-size>\n"
        "Generate random file of given size of alphabet"
        "0, 1, .., alphabet-size - 1\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }

  long length = atol(argv[1]) << 20;
  long sigma = atol(argv[2]);

  fprintf(stderr, "length = %ld\n", length);
  fprintf(stderr, "sigma = %ld\n", sigma);

  static const long bufsize = (4L << 20);
  unsigned char *buf = new unsigned char[bufsize];

  long written = 0;
  while (written < length) {
    fprintf(stderr, "\r%.2Lf%%", (100.L * written) / length);
    long towrite = std::min(length - written, bufsize);
    for (long j = 0; j < towrite; ++j)
      buf[j] = utils::random_int(0, sigma - 1);
    fwrite(buf, 1, towrite, stdout);
    written += towrite;
  }

  fprintf(stderr, "\rFinished.\n");
}

