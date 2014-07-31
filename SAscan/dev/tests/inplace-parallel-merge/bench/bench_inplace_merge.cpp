#include <cstdio>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <thread>
#include <unistd.h>

#include "utils.h"
#include "parallel_merge.h"

// Unused, times were the same as using 1 thread
// Means that sequential inplace merging works ok.
template<typename T>
void sequential(T *tab, long n1, long n2, int *gap) {
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
  int *gap = new int[n1 + 1];

  for (long i = 0; i < length; ++i) {
    if ((i % (1 << 20)) == 0)
      fprintf(stderr, "generating tab: %.2Lf%%\r", (100.L * i) / (n1 + n2));
    tab[i] = i * i;
  }
  fprintf(stderr, "\n");

  std::fill(gap, gap + n1 + 1, 0L);
  static const long bufsize = (1 << 20);
  long *buf = new long[bufsize];
  int *buf2 = new int[bufsize];
  long left = n2;
  while (left) {
    fprintf(stderr, "generating gap: %.2Lf%%\r", (100.L * (n2 - left)) / n2);
    for (long j = 0; j < bufsize; ++j) {
      buf[j] = utils::random_long(0, n1);
      buf2[j] = (int)utils::random_long(0L, std::min(5L, left));
      left -= buf2[j];
    }
    for (long j = 0; j < bufsize; ++j)
      gap[buf[j]] += buf2[j];
  }
  delete[] buf;
  delete[] buf2;
  fprintf(stderr, "\n");

  for (long t = 1; t <= 32; t *= 2) {
    long thr = std::min(24L, t);
    fprintf(stderr, "======= threads = %ld =======\n", thr);

    /*merge<int, 8 >(tab, n1, n2, gap, thr);
    merge<int, 9 >(tab, n1, n2, gap, thr);
    merge<int, 10>(tab, n1, n2, gap, thr);*/
    merge<int, 11>(tab, n1, n2, gap, thr);
    /*merge<int, 12>(tab, n1, n2, gap, thr);
    merge<int, 13>(tab, n1, n2, gap, thr);
    merge<int, 14>(tab, n1, n2, gap, thr);
    merge<int, 15>(tab, n1, n2, gap, thr);*/
  }

  delete[] tab;
  delete[] gap;
}

int main() {
  std::srand(std::time(0) + getpid());
  test(1L << 33, 1L << 33);
}
