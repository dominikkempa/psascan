#include <cstdio>
#include <cstdlib>

#include "utils.h"
#include "divsufsort64.h"
#include "uint40.h"


int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "%s FILE\n"
        "Compute suffix array of text stored in FILE and write to FILE.sa8 (using 64-bit integers).\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }

  fprintf(stderr, "Filename = %s\n", argv[1]);

  long length;
  unsigned char *text;
  fprintf(stderr, "Reading text... ");
  long double start = utils::wclock();
  utils::read_objects_from_file(text, length, argv[1]);
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

  long *sa = new long[length];
  fprintf(stderr, "Running divsufsort64... ");
  start = utils::wclock();
  divsufsort64(text, sa, length);
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
  std::string sa_filename = std::string(argv[1]) + ".sa5";

  fprintf(stderr, "Writing SA to file... ");
  start = utils::wclock();
  uint40 *sa40 = new uint40[length];
  for (long j = 0; j < length; ++j)
    sa40[j] = uint40(sa[j]);
  utils::write_objects_to_file(sa40, length, sa_filename);
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

  delete[] text;
  delete[] sa;
  delete[] sa40;
}

