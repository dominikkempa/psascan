#include <cstdio>
#include <cstring>
#include <string>
#include <algorithm>

#include "sais.hxx"
#include "utils.h"
#include "uint40.h"
#include "stream.h"

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "usage: %s <file>\n\n"
                    "Compute the SA of <file>\n",
                    argv[0]);
    std::exit(EXIT_FAILURE);
  }

  long length;
  unsigned char *text;
  utils::read_file(text, length, argv[1]);
  std::reverse(text, text + length);

  int *SA = new int[length];
  fprintf(stderr, "Sorting: ");
  long double sais_start = utils::wclock();
  saisxx(text, SA, (int)length); // peak 9n, right?
  fprintf(stderr, "%.2Lf\n", utils::wclock() - sais_start);

  fprintf(stderr, "Writing: ");
  long double writing_start = utils::wclock();
  std::string out_fname = std::string(argv[1]) + ".sa5";
  stream_writer<uint40> *writer =
    new stream_writer<uint40>(out_fname, 1 << 20);
  for (long j = 0, dbg = 0; j < length; ++j, ++dbg) {
    if (dbg == (1 << 23)) {
      fprintf(stderr, "\r%ld/%ld (%.1Lf%%)",
          j, length, (100.L * j) / length);
      dbg = 0;
    }
    writer->write(uint40((unsigned long)SA[j]));
  }
  delete writer;
  fprintf(stderr, "\nWriting time: %.2Lf\n", utils::wclock() - writing_start);

  delete[] text;
  delete[] SA;
}

