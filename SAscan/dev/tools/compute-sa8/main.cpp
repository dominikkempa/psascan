#include <cstdio>
#include <cstdlib>

#include "utils.h"
#include "divsufsort64.h"

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "%s FILE\n"
        "Display all bytes that occur in FILE.\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }

  long length;
  unsigned char *text;
  utils::read_objects_from_file(text, length, argv[1]);

  long *sa = new long[length];
  divsufsort64(text, sa, length);
  std::string sa_filename = std::string(argv[1]) + ".sa8";
  utils::write_objects_to_file(sa, length, sa_filename);

  delete[] text;
  delete[] sa;
}

