#include <cstdio>
#include <cstdlib>

#include "utils.h"

void gen(unsigned char *text, long length, long start) {
  if (!length) return;

  long mid = length / 2;
  text[mid] = start;
  gen(text, mid, start + 1);
  gen(text + mid + 1, length - mid - 1, start + 1);
}

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "%s <size-in-MiB>\n"
        "Generate Skyline word file of given length (in MiB)\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }

  long length = atol(argv[1]) << 20;
  fprintf(stderr, "length = %ld\n", length);

  unsigned char *text = new unsigned char[length];
  gen(text, length, 'a');

  for (long i = 0; i < length; ++i) {
    if (i % (1 << 22) == 0)
      fprintf(stderr, "%.2Lf%%\r", (100.L * i) / length);
    std::fputc(text[i], stdout);
  }

  delete[] text;
  fprintf(stderr, "Finished.\n");
}

