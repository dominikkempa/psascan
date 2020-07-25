#ifndef __PARALLEL_PERMUTE_H_INCLUDED
#define __PARALLEL_PERMUTE_H_INCLUDED

#include <cstdio>
#include <cstdlib>

#include <vector>
#include <stack>
#include <algorithm>
#include <thread>
#include <mutex>

template<typename T>
void parallel_permute(T *tab, T **index, std::mutex *mutexes,
    long length, long& selector, std::mutex &selector_mutex) {
  while (true) {
    // Find the starting point on some cycle.
    long start;
    while (true) {
      // Get the candidate using selector.
      std::unique_lock<std::mutex> lk(selector_mutex);
      while (selector < length && (index[selector] == tab + selector ||
            index[selector] < tab || tab + length <= index[selector]))
        ++selector;

      // Exit, if the selector does not give any candidate.
      if (selector == length) {
        lk.unlock();
        return;
      }

      // Unlock selector lock, allow other thread to look
      // for candidates in the meantime.
      start = selector++;
      lk.unlock();

      // Lock a candidate and check if still good.
      // If yes, keep lock and break.
      // NOTE: looks like its correct only with the
      // first check.
      if (mutexes[start].try_lock() && index[start] != tab + start &&
          tab <= index[start] && index[start] < tab + length) break;
    }

    // Invariant: we have found a good candidate.
    // and have lock on mutexes[start].

    // First, we create temporary space for *index[start]
    // and move *index[start] into that temp space.
    T *temp = new T;
    *temp = *index[start];
    std::swap(index[start], temp);
    mutexes[start].unlock();

    // We now have free space at temp,
    // we can write something there.
    //
    // Invariant: at all times, index[i] points to
    // element that should be placed at index i in tab.
    // 
    // Keep following the cycle.
    do {
      long next = temp - tab;
      std::unique_lock<std::mutex> lk(mutexes[next]);
      *temp = *index[next];
      std::swap(index[next], temp);
      lk.unlock();
    } while (tab <= temp && temp < tab + length);
    delete temp;
  }
}

//==============================================================================
// Permute elements inside tab as follows: tab[i] := *index[i];
// index[0..length) contain pointers to tab[0], tab[1], .., tab[length - 1].
//==============================================================================
template<typename T>
void permute(T *tab, T **index, long length, long n_threads) {
  long selector = 0;
  std::mutex selector_mutex;
  std::mutex *mutexes = new std::mutex[length];
  std::thread **threads = new std::thread*[n_threads];

  for (long i = 0; i < n_threads; ++i)
    threads[i] = new std::thread(parallel_permute<T>, tab, index,
        mutexes, length, std::ref(selector), std::ref(selector_mutex));

  for (long i = 0; i < n_threads; ++i) threads[i]->join();
  for (long i = 0; i < n_threads; ++i) delete threads[i];
  delete[] threads;
  delete[] mutexes;
}

#endif  // __PARALLEL_PERMUTE_H_INCLUDED
