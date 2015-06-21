#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <unistd.h>
#include <algorithm>

int main(int argc, char **argv) {
  if (argc != 2) std::exit(EXIT_FAILURE);
  srand(time(0) + getpid());

  FILE *f = fopen(argv[1], "w");

  static const long bufsize = (4 << 20);
  unsigned char *buf = new unsigned char[bufsize];
  for (int i = 0; i < bufsize; ++i) buf[i] = rand() % 50;

  static const long n = 10000000000L; // 10^10 bytes = (approx) 10GiB
  long i = 0;
  while (i < n) {
    int towrite = std::min(bufsize, n - i);
    fwrite(buf, 1, towrite, f);
    i += towrite;
  }

  fclose(f);
  delete[] buf;
}
