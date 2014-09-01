#include <cstdio>
#include <cstdlib>

#include "utils.h"


int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "%s FILE\n"
        "Display all bytes that occur in FILE.\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }

  long length;
  unsigned char *text;
  utils::read_objects_from_file(text, length, argv[1]);
  
  bool occurs[256] = {false};
  for (long i = 0; i < length; ++i)
    occurs[text[i]] = true;
  delete[] text;

  for (long i = 0; i < 256; ++i)
    if (occurs[i]) fprintf(stderr, "%ld ", i);
  fprintf(stderr, "\n");
}

