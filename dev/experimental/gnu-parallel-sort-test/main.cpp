#include <omp.h>
#include <cstdio>
#include <ctime>
#include <vector>
#include <algorithm>
#include <parallel/algorithm>
#include <sys/time.h>
#include <unistd.h>

// OS-specific timing
long double wclock() {
  timeval now;
  gettimeofday(&now, NULL);

  return now.tv_sec + now.tv_usec / 1000000.0L;
}


int main(int argc, char* argv[]) {
  std::srand(std::time(NULL) + getpid());

  if (argc != 2) {
    fprintf(stderr, "usage: %s <elems>\n"
        "Test gnu parallel sort on an array of <elems> ints\n",
        argv[0]);
    std::exit(EXIT_FAILURE);
  }

  long length = std::atol(argv[1]);

  printf("Sorting %ld elements\n", length);
  printf("Max. OpenMP threads: %d (%d processors)\n",
      omp_get_max_threads(), omp_get_num_procs());

  int *tab = new int[length];
  for (long i = 0L; i < length; i++)
    tab[i] = std::rand();
  fprintf(stderr, "Array takes: %.2LfMiB\n",
      (sizeof(int) * length) / (1024.L * 1024));

  long double start = wclock();

  //std::sort(tab, tab + length);

  // both versions of quicksort are inplace
  // quickosrt achieves only x 1.45 (yes, only 45%) speedup over std::sort
  //__gnu_parallel::sort(tab, tab + length, __gnu_parallel::balanced_quicksort_tag());
  // version above is 11% slower then the quicksort below
  //__gnu_parallel::sort(tab, tab + length, __gnu_parallel::quicksort_tag());


  // all versions of mergesort differ by less than 1% in runtime.
  // mergesort uses twice the space required for the input.
  // with 24 threads (12 cores), this is 13.5 times faster than std::sort.
  //__gnu_parallel::sort(tab, tab + length, __gnu_parallel::multiway_mergesort_exact_tag());
  //__gnu_parallel::sort(tab, tab + length, __gnu_parallel::multiway_mergesort_sampling_tag());
  //__gnu_parallel::sort(tab, tab + length);


  long double elapsed = wclock() - start;
  printf("__gnu_parallel::sort: %.3Lfs\n", elapsed);

  delete[] tab;
}
