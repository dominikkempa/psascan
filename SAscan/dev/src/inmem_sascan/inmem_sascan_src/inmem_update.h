#ifndef __INMEM_SASCAN_INMEM_UPDATE_H_INCLUDED
#define __INMEM_SASCAN_INMEM_UPDATE_H_INCLUDED

#include <thread>
#include <mutex>
#include <condition_variable>
#include <algorithm>

#include "utils.h"
#include "buffer.h"
#include "inmem_gap_array.h"

namespace inmem_sascan_private {

//==============================================================================
// This object creates a given number of threads that will perform gap array
// updates. Most of the time all threads are sleeping on a conditional variable.
// Once the buffer is available for processing, they are all woken up and
// perform the update in parallel. The caller waits until all threads are
// finished and puts the buffer in the poll of empty buffers.
//
// Only one object of this class should exist.
//==============================================================================
template<typename block_offset_type>
struct gap_parallel_updater {

  template<typename T>
  static void parallel_update(gap_parallel_updater<T> *updater, int id) {
    while (true) {
      // Wait until there is a buffer available or the
      // notification 'no more buffers' arrives (in that case exit).
      std::unique_lock<std::mutex> lk(updater->m_avail_mutex);
      while (!(updater->m_avail[id]) && !(updater->m_avail_no_more))
        updater->m_avail_cv.wait(lk);

      if (!(updater->m_avail[id]) && updater->m_avail_no_more) {
        // The msg was that there won't be more buffers -- exit.
        lk.unlock();
        return;
      }

      updater->m_avail[id] = false;
      lk.unlock();

      // Safely perform the update.
      buffer<T> *buf = updater->m_buffer;
      inmem_gap_array *gap = updater->m_gap_array;
      int beg = buf->sblock_beg[id];
      int end = beg + buf->sblock_size[id];

      for (int i = beg; i < end; ++i) {
        T x = buf->m_content[i];
        gap->m_count[x]++;

        // Check if values wrapped-around.
        if (gap->m_count[x] == 0) {
          gap->m_excess_mutex.lock();
          gap->m_excess.push_back(x);
          gap->m_excess_mutex.unlock();
        }
      }

      // Update the number of finished threads.
      bool finished_last = false;
      std::unique_lock<std::mutex> lk2(updater->m_finished_mutex);
      updater->m_finished++;
      if (updater->m_finished == updater->m_threads_cnt)
        finished_last = true;
      lk2.unlock();

      // If this was the last thread finishing, let the caller know.
      if (finished_last)
        updater->m_finished_cv.notify_one();
    }
  }

  gap_parallel_updater(inmem_gap_array *gap_array, int threads_cnt)
      : m_gap_array(gap_array),
        m_threads_cnt(threads_cnt),
        m_avail_no_more(false) {
    m_avail = new bool[m_threads_cnt];
    std::fill(m_avail, m_avail + m_threads_cnt, false);
    m_threads = new std::thread*[m_threads_cnt];

    // After this threads immediatelly hang up on m_avail_cv.
    for (int i = 0; i < m_threads_cnt; ++i)
      m_threads[i] = new std::thread(parallel_update<block_offset_type>, this, i);
  }

  ~gap_parallel_updater() {
    // Signal to all threads to finish.
    std::unique_lock<std::mutex> lk(m_avail_mutex);
    m_avail_no_more = true;
    lk.unlock();
    m_avail_cv.notify_all();

    // Wait until they actually finish and release memory.
    for (int i = 0; i < m_threads_cnt; ++i) {
      m_threads[i]->join();
      delete m_threads[i];
    }
    delete[] m_threads;
    delete[] m_avail;
  }

  void update(buffer<block_offset_type> *buffer) {
    // Prepare a message for each thread that new buffer is available.
    std::unique_lock<std::mutex> lk(m_avail_mutex);
    m_finished = 0;
    m_buffer = buffer;
    for (int i = 0; i < m_threads_cnt; ++i)
      m_avail[i] = true;
    lk.unlock();

    // Wake up all threads, they will now perform the update.
    m_avail_cv.notify_all();

    // Wait until all threads report that they are done.
    std::unique_lock<std::mutex> lk2(m_finished_mutex);
    while (m_finished != m_threads_cnt)
      m_finished_cv.wait(lk2);
    lk2.unlock();

    // We are done processing the buffer. The caller of this method
    // can now place the buffer into the poll of empty buffers.
  }

private:
  inmem_gap_array *m_gap_array;

  std::thread **m_threads;
  int m_threads_cnt;

  buffer<block_offset_type> *m_buffer;

  // For notifying threads about available buffer.
  std::mutex m_avail_mutex;
  std::condition_variable m_avail_cv;
  bool *m_avail;
  bool m_avail_no_more;

  // The mutex below is to protect m_finished. The condition
  // variable allows the caller to wait (and be notified when done)
  // until threads finish processing their section of buffer.
  int m_finished;
  std::mutex m_finished_mutex;
  std::condition_variable m_finished_cv;
};

template<typename block_offset_type>
void inmem_gap_updater(buffer_poll<block_offset_type> *full_buffers,
    buffer_poll<block_offset_type> *empty_buffers,
    inmem_gap_array *gap, long n_increasers) {

  gap_parallel_updater<block_offset_type> *updater =
    new gap_parallel_updater<block_offset_type>(gap, n_increasers);

  while (true) {
    // Get a buffer from the poll of full buffers.
    std::unique_lock<std::mutex> lk(full_buffers->m_mutex);
    while (!full_buffers->available() && !full_buffers->finished())
      full_buffers->m_cv.wait(lk);

    if (!full_buffers->available() && full_buffers->finished()) {
      // All workers finished -- exit, but before, let let other
      // updating threads know that they also should check for exit.
      lk.unlock();
      full_buffers->m_cv.notify_one();
      break;
    }

    buffer<block_offset_type> *b = full_buffers->get();
    lk.unlock();
    full_buffers->m_cv.notify_one(); // let others know they should try

    // Process buffer.
    updater->update(b);

    // Add the buffer to the poll of empty buffers and notify waiting thread.
    std::unique_lock<std::mutex> lk2(empty_buffers->m_mutex);
    empty_buffers->add(b);
    lk2.unlock();
    empty_buffers->m_cv.notify_one();
  }

  delete updater;
}

}  // namespace inmem_sascan

#endif // __INMEM_UPDATE_H_INCLUDED
