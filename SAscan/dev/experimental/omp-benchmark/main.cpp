#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <omp.h>

#include "utils.h"

int main() {
  long length = 100000000; // 100M
  int *tab = new int[length];
  int *phi = new int[length];
  for (int i = 0; i < length; ++i) tab[i] = i;
  std::random_shuffle(tab, tab + length);

  clock_t start = std::clock();
  long double wstart = utils::wclock();

  int i;
  omp_set_dynamic(0);
  for (int j = 0; j < 10; ++j) {
    #pragma omp parallel for shared(tab, phi, length) private(i) num_threads(2)
    for (i = 1; i < length; ++i) {
      int val = tab[i], val_prev = tab[i - 1];
      phi[val] = val_prev;
    } 
  }

  clock_t end = std::clock();
  long double wend = utils::wclock();
  long double wtime = wend - wstart;
  long double cputime = (long double)(end - start) / CLOCKS_PER_SEC;

  fprintf(stderr, "CPU:  %.2Lf sec\n", cputime);
  fprintf(stderr, "Real: %.2Lf sec\n", wtime);

  delete[] tab;
  delete[] phi;
}

