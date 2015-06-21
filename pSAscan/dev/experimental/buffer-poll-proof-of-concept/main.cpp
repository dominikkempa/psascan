#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <unistd.h>

#include <iostream>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>

#include "utils.h"

std::mutex global_mutex;
std::mutex stdout_mutex;
long global_sum;

struct buffer {  
  buffer(long size)
      : m_filled(0L),
        m_size(size) {
    m_content = new long[size];
  }
  
  ~buffer() {
    delete[] m_content;
  }

  long m_filled, m_size;
  long *m_content;
};

// Same class for the poll of empty and full buffers.    
struct buffer_poll {
  buffer_poll(long worker_threads = 0L) {
    m_worker_threads = worker_threads; // unused for the poll of empty buffers.
    m_worker_threads_finished = 0L;
  }
  
  void add(buffer *b) {
    m_queue.push(b);
  }
  
  bool available() const {
    return m_queue.size() > 0;
  }

  buffer *get() {
    if (m_queue.empty()) {
      fprintf(stderr, "Error: requesting a buffer from empty poll!\n");
      std::exit(EXIT_FAILURE);
    }

    buffer *ret = m_queue.front();
    m_queue.pop();

    return ret;
  }

  bool finished() const {
    return m_worker_threads_finished == m_worker_threads;
  }

  void increment_finished_workers() {
    ++m_worker_threads_finished;
  }

  std::condition_variable m_cv;
  std::mutex m_mutex;

private:
  long m_worker_threads; // used to detect that worker threads are done
  long m_worker_threads_finished;

  std::queue<buffer*> m_queue;
};

void update(buffer_poll *full_buffers, buffer_poll *empty_buffers) {
  while (true) {
    // Get a buffer from the poll of full buffers.
    std::unique_lock<std::mutex> lk(full_buffers->m_mutex);
    while (!full_buffers->available() && !full_buffers->finished())
      full_buffers->m_cv.wait(lk);

    if (!full_buffers->available() && full_buffers->finished()) {
      // All workers finished. We're exiting, but before, we're letting
      // other updating threads know that they also should exit.
      lk.unlock();
      full_buffers->m_cv.notify_one();
      break;
    }
    buffer *b = full_buffers->get();
    lk.unlock();
    full_buffers->m_cv.notify_one(); // let others know they should try

    // Process buffer.
    long sum = 0L;
    for (long i = 0L; i < b->m_filled; ++i)
      sum += b->m_content[i];
    global_mutex.lock();
    global_sum += sum;
    global_mutex.unlock();

    // Add the buffer to the poll of empty buffers and notify waiting thread.
    std::unique_lock<std::mutex> lk2(empty_buffers->m_mutex);
    empty_buffers->add(b);
    lk2.unlock();
    empty_buffers->m_cv.notify_one();
  }
}

void work(buffer_poll *full_buffers, buffer_poll *empty_buffers, long beg, long end) {
  long range_size = end - beg;
  long processed = 0;
  while (processed < range_size) {
    // Get a buffer from the poll of empty buffers.
    std::unique_lock<std::mutex> lk(empty_buffers->m_mutex);
    while (!empty_buffers->available())
      empty_buffers->m_cv.wait(lk);
    buffer *b = empty_buffers->get();
    lk.unlock();
    empty_buffers->m_cv.notify_one(); // let others know they should re-check

    // Process buffer: fill with numbers.
    b->m_filled = std::min(range_size - processed, b->m_size);
    for (long i = 0L; i < b->m_filled; ++i)
      b->m_content[i] = beg + processed++;
    
    // Add the buffer to the poll of full buffers and notify waiting thread.
    std::unique_lock<std::mutex> lk2(full_buffers->m_mutex);
    full_buffers->add(b);
    lk2.unlock();
    full_buffers->m_cv.notify_one();
  }

  // Report that another worker thread has finished.
  std::unique_lock<std::mutex> lk(full_buffers->m_mutex);
  full_buffers->increment_finished_workers();
  lk.unlock();
  
  // Notify waiting update threads in case no more buffers
  // are going to be produces by worker threads.
  full_buffers->m_cv.notify_one();
}

int main() {
  std::srand(std::time(0) + getpid());
  
  long tested = 0L;
  while (true) {
    ++tested;
    if (tested % 100 == 0)
      fprintf(stderr, "\rTested %ld", tested);

    global_sum = 0L;
      
    // Create some number of different buffers.
    long n_buffers = (long)utils::random_int(1, 10);
    buffer **buffers = new buffer*[n_buffers];
    for (long i = 0L; i < n_buffers; ++i)
      buffers[i] = new buffer((long)utils::random_int(1, 10));

    // Set the number of working/updating threads.
    long n_workers = (long)utils::random_int(1, 10);
    long n_updaters = (long)utils::random_int(1, 10);

    // Create polls with buffers.
    buffer_poll *empty_buffers = new buffer_poll();
    buffer_poll *full_buffers = new buffer_poll(n_workers);

    // Add buffers to the empty poll.
    for (long i = 0L; i < n_buffers; ++i)
      empty_buffers->add(buffers[i]);
      
    long range = 0L; // total range of numbers processed by all threads.

    // Create working threads.
    std::thread **workers = new std::thread*[n_workers];
    for (long i = 0L; i < n_workers; ++i) {
      long this_range = (long)utils::random_int(1, 10);
      workers[i] = new std::thread(work, full_buffers, empty_buffers, range, range + this_range);
      range += this_range;
    }

    // Create updating threads.
    std::thread **updaters = new std::thread*[n_updaters];
    for (long i = 0L; i < n_updaters; ++i)
      updaters[i] = new std::thread(update, full_buffers, empty_buffers);

    // Wait for all threads to finish.
    for (long i = 0L; i < n_workers; ++i) workers[i]->join();
    for (long i = 0L; i < n_updaters; ++i) updaters[i]->join();

    // Clean up.
    for (long i = 0L; i < n_workers; ++i)
      delete workers[i];
    delete[] workers;
    for (long i = 0L; i < n_updaters; ++i)
      delete updaters[i];
    delete[] updaters;
    for (long i = 0L; i < n_buffers; ++i)
      delete buffers[i];
    delete[] buffers;
    delete empty_buffers;
    delete full_buffers;

    long sum_correct = 0L;
    for (long i = 0L; i < range; ++i)
      sum_correct += i;
      
    if (global_sum != sum_correct) {
      fprintf(stderr, "ERROR: communication failed.\n");
      fprintf(stderr, "global sum = %ld\n", global_sum);
      fprintf(stderr, "sum_correct = %ld\n", sum_correct);
      std::exit(EXIT_FAILURE);
    }
  }
}
