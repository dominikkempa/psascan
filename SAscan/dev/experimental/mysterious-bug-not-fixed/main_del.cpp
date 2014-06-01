#include <cstdio>
#include <cstdlib>
#include <errno.h>
#include <cstring>
#include <string>

#include "uint40.h"

void work(std::string fname, long elems, long bufelems) {
  unsigned char *tmp = new unsigned char[39L << 30]; // 39 GiB
  for (long i = 0L; i < (39L << 30); ++i) tmp[i] = 0;
  uint40 *buf = new uint40[bufelems];

  FILE *f = fopen(fname.c_str(), "r");
  if (f == NULL) exit(1);
  long r = fread(buf, sizeof(uint40), elems, f);
  if (r != elems) exit(1);
  fclose(f);

  // This does not work. But try altering something seemingly
  // unrelated in the code (e.g. delete the initialization of
  // tmp or change 39L to 29L) and it will work :|
  std::string cmd = "rm " + fname;
  int res = system(cmd.c_str());
  printf("%s\n", res ? "FAILED" : "OK");

  // This always works
  // if (remove(fname.c_str()))
  //   fprintf(stderr, "FAILED: %s\n", strerror(errno));
  // else fprintf(stderr, "OK\n");

  delete[] buf;
  delete[] tmp;
}

int main(int argc, char **argv) {
  if (argc != 2) exit(1);
  work(argv[1], 2000000000L, 8L << 30);
}
