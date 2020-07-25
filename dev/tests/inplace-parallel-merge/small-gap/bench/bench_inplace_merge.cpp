#include <cstdio>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <thread>
#include <unistd.h>

#include "utils.h"
#include "parallel_merge.h"
#include "gap_array.h"

void test(long n1, long n2) {
  fprintf(stderr, "n1 = %ld\n", n1);
  fprintf(stderr, "n2 = %ld\n", n2);

  long length = n1 + n2;
  int *tab = new int[length];

  for (long i = 0; i < length; ++i) {
    if ((i % (1 << 20)) == 0)
      fprintf(stderr, "\rgenerating tab: %.2Lf%%", (100.L * i) / (n1 + n2));
    tab[i] = i * i;
  }
  fprintf(stderr, "\n");

  gap_array *gap = new gap_array(n1 + 1);
  static const long bufsize = (1 << 20);
  long *buf = new long[bufsize];
  int *buf2 = new int[bufsize];
  long left = n2;
  while (left) {
    fprintf(stderr, "\rgenerating gap: %.2Lf%%", (100.L * (n2 - left)) / n2);
    for (long j = 0; j < bufsize; ++j) {
      buf[j] = utils::random_long(0, n1);
      buf2[j] = (int)utils::random_long(0L, std::min(5L, left));
      left -= buf2[j];
    }
    for (long j = 0; j < bufsize; ++j)
      for (long jj = 0; jj < buf2[j]; ++jj)
        gap->increment(buf[j]);
  }
  delete[] buf;
  delete[] buf2;
  fprintf(stderr, "\n");
  fprintf(stderr, "excess size = %ld\n", (long)gap->m_excess.size());
  std::sort(gap->m_excess.begin(), gap->m_excess.end());

  for (long t = 1; t <= 32; t *= 2) {
    long thr = std::min(24L, t);
    fprintf(stderr, "======= threads = %ld =======\n", thr);

    merge<int, 9 >(tab, n1, n2, gap, thr);
    merge<int, 10>(tab, n1, n2, gap, thr);
    merge<int, 11>(tab, n1, n2, gap, thr);
    merge<int, 12>(tab, n1, n2, gap, thr);
    merge<int, 13>(tab, n1, n2, gap, thr);
  }

  delete[] tab;
  delete gap;
}

int main() {
  std::srand(std::time(0) + getpid());
  test(1L << 30, 1L << 30);
}
