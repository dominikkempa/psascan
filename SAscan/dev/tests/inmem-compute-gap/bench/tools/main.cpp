#include <cstdio>
#include <cstring>

#include "divsufsort.h"
#include "utils.h"

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s <file>\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }

  long double start;

  fprintf(stderr, "Reading text: ");
  start = utils::wclock();
  long length;
  unsigned char *text;
  utils::read_objects_from_file(text, length, argv[1]);
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

  fprintf(stderr, "Computing sa: ");
  start = utils::wclock();
  int *sa = new int[length];
  divsufsort(text, sa, (int)length);
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

  fprintf(stderr, "Storing sa to file: ");
  start = utils::wclock();
  std::string out_fname = std::string(argv[1]) + ".sa";
  utils::write_objects_to_file(sa, length, out_fname.c_str());
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

  delete[] sa;
  delete[] text;
}
