#include <cstdio>
#include <cstdlib>

#include "utils.h"


int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "%s <size-in-MiB>\n"
        "Generate Fibonacci word file of given length (in MiB)\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }

  long length = atol(argv[1]) << 20;
  fprintf(stderr, "length = %ld\n", length);

  unsigned char *text = new unsigned char[length];
  text[0] = 'a';
  text[1] = 'b';
  long prev = 1, cur = 2;
  while (cur < length) {
    long chunk = std::min(prev, length - cur);
    prev = cur;
    for (long j = 0; j < chunk; ++j)
      text[cur++] = text[j];
  }

  for (long i = 0; i < length; ++i) {
    if (i % (1 << 22) == 0)
      fprintf(stderr, "%.2Lf%%\r", (100.L * i) / length);
    std::fputc(text[i], stdout);
  }

  delete[] text;
  fprintf(stderr, "Finished.\n");
}

