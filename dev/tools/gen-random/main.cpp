#include <cstdio>
#include <cstdlib>
#include <string>

#include "utils.h"


int main(int argc, char **argv) {
  if (argc != 4) {
    fprintf(stderr, "%s SIZE SIGMA OUTFILE\n"
        "Generate random file of SIZE GiB over alphabet of size SIGMA and write to OUTFILE\n",
        argv[0]);
    std::exit(EXIT_FAILURE);
  }

  std::size_t length = atol(argv[1]) << 30;
  std::size_t sigma = atol(argv[2]);
  std::string output_filename = argv[3];
  std::FILE *f = utils::open_file(output_filename, "w");

  fprintf(stderr, "Length = %ld\n", length);
  fprintf(stderr, "Sigma = %ld\n", sigma);
  fprintf(stderr, "Output file = %s\n", output_filename.c_str());

  static const std::size_t bufsize = (4UL << 20);
  unsigned char *buf = new unsigned char[bufsize];

  std::size_t written = 0;
  while (written < length) {
    fprintf(stderr, "\r%.2Lf%%", (100.L * written) / length);
    std::size_t towrite = std::min(length - written, bufsize);
    for (std::size_t j = 0; j < towrite; ++j)
      buf[j] = utils::random_int(0, sigma - 1);
    fwrite(buf, 1, towrite, f);
    written += towrite;
  }

  std::fclose(f);
  fprintf(stderr, "\rFinished.\n");
}

