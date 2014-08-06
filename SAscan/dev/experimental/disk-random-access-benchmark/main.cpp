#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <unistd.h>

#include "utils.h"

void test(const char *filename) {
  long size = utils::file_size(filename);

  static const long n_queries = 100;
  long *queries = new long[n_queries];
  for (long i = 0; i < n_queries; ++i)
    queries[i] = utils::random_long(0, size - 1);

  std::FILE *f = utils::open_file(filename, "r");
  long double start = utils::wclock();
  unsigned char dest;
  for (long i = 0; i < n_queries; ++i) {
    std::fseek(f, queries[i], SEEK_SET);
   std::fread(&dest, 1, 1, f);
  }
  long double elapsed = utils::wclock() - start;
  std::fclose(f);

  fprintf(stderr, "Time for %ld queries: %.2Lfs\n",
      n_queries, elapsed);
  fprintf(stderr, "Speed: %.2Lf queries/sec\n",
      (long double)n_queries / elapsed);
  fprintf(stderr, "Access time: %.2Lfms\n",
      (1000.L * elapsed) / n_queries);
}

int main(int argc, char **argv) {
  std::srand(std::time(0) + getpid());
  if (argc != 2) {
    fprintf(stderr, "Usage: %s <file>\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }
  test(argv[1]);
}
