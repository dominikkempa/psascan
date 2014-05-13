#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>

#include <thread>
#include <mutex>
#include <condition_variable>

#include <unistd.h>

#include "utils.h"

struct buffer {  
  buffer(long size)
      : m_filled(0L),
        m_size(size),
        is_full(false) {
    m_content = new long[size];
  }
  
  ~buffer() {
    delete[] m_content;
  }

  long m_filled;
  long m_size;

  bool is_full;

  long *m_content;
  std::mutex m_mutex;
  std::condition_variable m_cv;
};

void update(buffer *active, buffer *inactive, long *output, long toprocess) {
  long processed = 0L;
  while (processed < toprocess) {
    std::unique_lock<std::mutex> lk(active->m_mutex);
    while (active->is_full == false)
      active->m_cv.wait(lk);
    for (long i = 0; i < active->m_filled; ++i)
      output[processed++] = active->m_content[i];
    active->is_full = false;
    lk.unlock();
    active->m_cv.notify_all();
    std::swap(active, inactive);
  }
}

void create(buffer *active, buffer *inactive, long toprocess) {
  long processed = 0L;
  while (processed < toprocess) {
    std::unique_lock<std::mutex> lk(active->m_mutex);
    while (active->is_full == true)
      active->m_cv.wait(lk);
    active->m_filled = std::min(toprocess - processed, active->m_size);
    for (long i = 0L; i < active->m_filled; ++i)
      active->m_content[i] = processed++;
    active->is_full = true;
    lk.unlock();
    active->m_cv.notify_all();
    std::swap(active, inactive);
  }
}

int main() {
  std::srand(std::time(0) + getpid());
  
  static const int toprocess = 1000000;
  long *output = new long[toprocess];
  
  long tested = 0L;
  while (true) {
    ++tested;
    if (tested % 10 == 0)
      fprintf(stderr, "\rTested %ld", tested);
    buffer *b1 = new buffer((long)utils::random_int(1, 1000));
    buffer *b2 = new buffer((long)utils::random_int(1, 1000));
    std::thread *cr = new std::thread(create, b1, b2, toprocess);
    std::thread *up = new std::thread(update, b1, b2, output, toprocess);
    
    cr->join();
    up->join();

    bool ok = true;
    for (long i = 0L; i < toprocess; ++i)
      if (output[i] != i) { ok = false; break; }
    if (!ok) {
      fprintf(stderr, "ERROR: communication failed.\n");
      std::exit(EXIT_FAILURE);
    }
    
    delete b1;
    delete b2;

    delete cr;
    delete up;
  }
  
  delete[] output;
}
