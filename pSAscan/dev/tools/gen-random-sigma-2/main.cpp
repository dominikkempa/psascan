#include <cstdio>
#include <cstdlib>
#include <string>
#include <ctime>
#include <unistd.h>

#include "utils.h"


int main(int argc, char **argv) {
  srand(time(0) + getpid());

  if (argc != 3) {
    fprintf(stderr, "%s SIZE OUTFILE\n"
        "Generate random file of SIZE GiB over binary alphabet and write to OUTFILE\n",
        argv[0]);
    std::exit(EXIT_FAILURE);
  }

  std::size_t length = atol(argv[1]) << 30;
  std::size_t sigma = 255;
  std::string output_filename = argv[2];
  std::FILE *f = utils::open_file(output_filename, "w");

  fprintf(stderr, "Length = %ld\n", length);
  fprintf(stderr, "Output file = %s\n", output_filename.c_str());

  static const std::size_t bufsize = (4UL << 20);
  unsigned char *buf = new unsigned char[bufsize];

  std::size_t written = 0;
  while (written < length) {
    fprintf(stderr, "\r%.2Lf%%", (100.L * written) / length);
    std::size_t towrite = std::min(length - written, bufsize);
    std::size_t filled = 0;
    while (filled < towrite) {
      unsigned char x = utils::random_int(0, sigma - 1);
      std::size_t totake = std::min(towrite - filled, 8UL);
      for (std::size_t j = 0; j < totake; ++j)
        buf[filled++] = ((x >> (totake - 1 - j)) & 0x1);
    }
    fwrite(buf, 1, towrite, f);
    written += towrite;
  }

  std::fclose(f);
  fprintf(stderr, "\rFinished.\n");
}

