#include <cstdio>
#include <cstdlib>

#include <mutex>
#include <thread>

#include "utils.h"

std::mutex sum_mutex;

void add_random(int *tab, long *sum, int i) {
//  std::unique_lock<std::mutex> lk(sum_mutex);
  *sum += tab[i];
//  lk.unlock();
}

int main() {
  static const int n = (10 << 20);
  int *tab = new int[n];
  for (int i = 0; i < n; ++i)
    tab[i] = rand() % n;

  static const long n_threads = 20000;
  std::thread **threads = new std::thread*[n_threads];

  long sum = 0L;
  long double start = utils::wclock();
  for (int i = 0; i < n_threads; ++i) threads[i] = new std::thread(add_random, tab, &sum, i);
  for (int i = 0; i < n_threads; ++i) threads[i]->join();
  long double elapsed = utils::wclock() - start;

  fprintf(stderr, "Created %ld threads in %.4Lfsec\n", n_threads, elapsed);
  fprintf(stderr, "Single thread creation: %.8Lfsec\n", elapsed / n_threads);
  fprintf(stderr, "sum = %ld\n", sum);

  for (int i = 0; i < n_threads; ++i) delete threads[i];
  delete[] threads;
  delete[] tab;
}

