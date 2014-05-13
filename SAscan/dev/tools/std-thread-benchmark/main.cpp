#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <algorithm>
#include <thread>
#include <mutex>
#include <unistd.h>

#include "utils.h"

std::mutex stdout_mutex;

void invert(int *tab, int *phi, int beg, int end) {
  long double start = utils::wclock();
  stdout_mutex.lock();
  std::cerr << "Thread " << std::this_thread::get_id() << " begins.\n";
  std::cerr << "Thread " << std::this_thread::get_id() << " range: [" << beg << ", " << end << "]\n";
  stdout_mutex.unlock();

  for (int i = beg; i < end; ++i) {
    int val = tab[i];
    int val_prev = tab[i - 1];
    phi[val] = val_prev;
  }

  long double elapsed = utils::wclock() - start;
  stdout_mutex.lock();
  std::cerr << "Thread " << std::this_thread::get_id() << " finished in "
    << elapsed << "secs.\n";
  stdout_mutex.unlock();
}

int main(int argc, char **argv) {
  std::srand(std::time(0) + getpid());

  if (argc != 2) {
    fprintf(stderr, "usage: %s <n_thread>\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }

  int n_thread = std::atoi(argv[1]);
  fprintf(stderr, "Using %d thread%s.\n", n_thread,
      (n_thread > 1) ? "s" : "");

  long length = 1000000000; // 100M
  fprintf(stderr, "Allocating...");
  int *tab = new int[length];
  int *phi = new int[length];
  fprintf(stderr, "DONE\n");

  fprintf(stderr, "Assigning...");
  for (int i = 0; i < length; ++i)
    tab[i] = i;
  fprintf(stderr, "DONE\n");

  fprintf(stderr, "Random shuffling..");
  std::random_shuffle(tab, tab + length);
  fprintf(stderr, "DONE\n");

  clock_t start = std::clock();
  long double wstart = utils::wclock();

  static const int repeats = 20;

  std::vector<long double> cputimes;
  std::vector<long double> wtimes;

  for (int i = 0; i < repeats; ++i) {
    long double single_start = std::clock();
    long double single_wstart = utils::wclock();

    std::thread **threads = new std::thread*[n_thread];
    for (int j = 0; j < n_thread; ++j)
      threads[j] = new std::thread(invert, tab, phi,
          std::max(1L,           j * (length / n_thread)),
          std::min(length, (j + 1) * (length / n_thread))
      );

    for (int j = 0; j < n_thread; ++j)
      threads[j]->join();

    for (int j = 0; j < n_thread; ++j)
      delete threads[j];

    delete[] threads;
    long double single_wtime = utils::wclock() - single_wstart;
    long double single_end = std::clock();
    long double single_cputime = (long double)(single_end - single_start) / CLOCKS_PER_SEC;
    cputimes.push_back(single_cputime);
    wtimes.push_back(single_wtime);
  }

  clock_t end = std::clock();
  long double wend = utils::wclock();
  long double wtime = wend - wstart;
  long double cputime = (long double)(end - start) / CLOCKS_PER_SEC;

  cputime /= repeats;
  wtime /= repeats;

  std::sort(cputimes.begin(), cputimes.end());
  std::sort(wtimes.begin(), wtimes.end());
  /*fprintf(stderr, "cputimes: ");
  for (unsigned i = 0; i < cputimes.size(); ++i)
    fprintf(stderr, "%.2Lf ", cputimes[i]);
  fprintf(stderr, "\n");*/
  fprintf(stderr, "Real times: ");
  for (unsigned i = 0; i < wtimes.size(); ++i)
    fprintf(stderr, "%.2Lf ", wtimes[i]);
  fprintf(stderr, "\n");

  fprintf(stderr, "Times (average of %d):\n", repeats);
  //fprintf(stderr, "  CPU:  %.2Lf sec\n", cputime);
  fprintf(stderr, "  Real: %.2Lf sec\n", wtime);
  //fprintf(stderr, "  CPU / Real = %.5Lf\n", cputime / wtime);
  fprintf(stderr, "Times (median of %d):\n", repeats);
  //fprintf(stderr, "  CPU:  %.2Lf sec\n", cputimes[repeats / 2]);
  fprintf(stderr, "  Real: %.2Lf sec\n", wtimes[repeats / 2]);
  //fprintf(stderr, "  CPU / Real = %.5Lf\n", cputimes[repeats / 2] / wtimes[repeats / 2]);

  delete[] tab;
  delete[] phi;
}

