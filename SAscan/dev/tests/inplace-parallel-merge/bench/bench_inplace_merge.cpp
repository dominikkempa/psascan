#include <cstdio>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <thread>
#include <unistd.h>

#include "utils.h"
#include "parallel_merge.h"

template<typename T>
void sequential(T *tab, long n1, long n2, long *gap) {
  long length = n1 + n2;

  long double start = utils::wclock();
  T *correct = new T[length];
  fprintf(stderr, "Merging sequentially: ");
  long right_ptr = n1, out_ptr = 0;
  for (long j = 0; j < gap[0]; ++j) correct[out_ptr++] = tab[right_ptr++];
  for (long j = 0; j < n1; ++j) {
    correct[out_ptr++] = tab[j];
    for (long k = 0; k < gap[j + 1]; ++k)
      correct[out_ptr++] = tab[right_ptr++];
  }
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);
  delete[] correct;
}


void test(long n1, long n2) {
  fprintf(stderr, "n1 = %ld\n", n1);
  fprintf(stderr, "n2 = %ld\n", n2);

  long length = n1 + n2;
  int *tab = new int[length];
  long *gap = new long[n1 + 1];
  long *gap2 = new long[n1 + 1];

  for (long i = 0; i < length; ++i) {
    if ((i % (1 << 20)) == 0)
      fprintf(stderr, "generating tab: %.2Lf%%\r", (100.L * i) / (n1 + n2));
    tab[i] = i * i;
  }
  fprintf(stderr, "\n");

  std::fill(gap2, gap2 + n1 + 1, 0L);
  static const long bufsize = (1 << 20);
  long *buf = new long[bufsize];
  for (long j = 0; j < n2; j += bufsize) {
      fprintf(stderr, "generating gap: %.2Lf%%\r", (100.L * j) / n2);
    long this_buf_size = std::min(n2 - j, bufsize);
    for (long jj = 0; jj < this_buf_size; ++jj)
      buf[jj] = utils::random_long(0, n1);
    for (long jj = 0; jj < this_buf_size; ++jj)
      ++gap2[buf[jj]];
  }
  delete[] buf;
  fprintf(stderr, "\n");

  // The same as inplace using 1 thread -- quite nice :)
  // sequential<int>(tab, n1, n2, gap2);

  for (long t = 1; t <= 32; t *= 2) {
    long thr = std::min(24L, t);
    fprintf(stderr, "======= threads = %ld =======\n", thr);

    std::copy(gap2, gap2 + n1 + 1, gap);
    merge<int, 11>(tab, n1, n2, gap, thr);
    // Any pagesize between 2^8 - 2^13 was equally good.
    // Of course the larger, the better: smaller
    // overhead of the sequential scans over page index.
  }

  delete[] tab;
  delete[] gap;
  delete[] gap2;
}

int main() {
  std::srand(std::time(0) + getpid());
  test(1L << 32, 1L << 32);
}
