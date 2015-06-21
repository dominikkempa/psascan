// Class implementing a file that we can read sequentially and simultanously
// release the disk space for the already scanned prefix.
//
// The constructor is given a `filename' (partial files with be named
// filename.part0, filename.part1, ...) and a maximal size (in bytes) of a
// single file.
//
// Example: create, write and read a sequence distributed into 10MiB files.
// ----------------------------------------------------------------------------
// static const int n = (100 << 20);              // prepare some data to write
// int *t = new int[n], i, S = 0;
// for (i=0; i<n; ++i) t[i] = 5*i;
//
// distributed_file<int> *f = new distributed_file<int>("input.sa", 10 << 20);
// f->initialize_writing(2 << 20);                   // 2MiB buffer for writing
// for (i=0; i<n; ++i) f->write(t[i]);                         // single writes
// f->write(t, t + n);                        // block writes are also possible
// for (i=0; i<b; ++i) f->write(t[i]);   // both write types can be interleaved
// f->finish_writing();                   // all files remain closed after this
//
// f->initialize_reading(1 << 20);                   // 1MiB buffer for reading
// for (i=0; i<n; ++i) S += f->read();   // simultanously destroys scanned file
// f->finish_reading();   // NOTE: this fails if you don't read the whole file!
//
// delete f;
// ----------------------------------------------------------------------------

#ifndef __PSASCAN_SRC_DISTRIBUTED_FILE_H_INCLUDED
#define __PSASCAN_SRC_DISTRIBUTED_FILE_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <string>
#include <thread>
#include <mutex>
#include <algorithm>
#include <condition_variable>

#include "utils.h"


namespace psascan_private {

template<typename value_type>
struct distributed_file {
  distributed_file(std::string filename_base, long max_bytes) {
    m_state = STATE_INIT;
    m_max_items = std::max(1UL, max_bytes / sizeof(value_type));
    m_filename = filename_base + ".distrfile." + utils::random_string_hash();
  }

  distributed_file(std::string filename_base, long max_bytes,
      const value_type *begin, const value_type *end) {
    m_state = STATE_INIT;
    m_max_items = std::max(1UL, max_bytes / sizeof(value_type));
    m_filename = filename_base + ".distrfile." + utils::random_string_hash();

    initialize_writing();
    write(begin, end);
    finish_writing();
  }


  void initialize_writing() {
    if (m_state != STATE_INIT) {
      fprintf(stderr, "\nError: initializing writing in state %s\n", state_string().c_str());
      std::exit(EXIT_FAILURE);
    }

    m_state = STATE_WRITING;
    m_total_write = 0;
    m_files_cnt = 0;
    make_new_file();
  }

  void write(const value_type *begin, const value_type *end) {
    if (m_state != STATE_WRITING) {
      fprintf(stderr, "\nError: write in state %s\n", state_string().c_str());
      std::exit(EXIT_FAILURE);
    }

    // Fill the current file.
    if (m_cur_file_write != m_max_items) {
      long left = m_max_items - m_cur_file_write;
      long towrite = std::min(left, end - begin);
      utils::add_objects_to_file(begin, towrite, m_file);
      m_cur_file_write += towrite;
      m_total_write += towrite;
      begin += towrite;
    }

    // Write remaining items.
    while (begin < end) {
      std::fclose(m_file);
      make_new_file();

      long towrite = std::min(m_max_items, end - begin);
      utils::add_objects_to_file(begin, towrite, m_file);
      m_cur_file_write += towrite;
      m_total_write += towrite;
      begin += towrite;
    }
  }

  void finish_writing() {
    if (m_state != STATE_WRITING) {
      fprintf(stderr, "\nError: finishing writing when in state %s\n", state_string().c_str());
      std::exit(EXIT_FAILURE);
    }
    if (m_cur_file_write == 0) {
      fprintf(stderr, "\nError: nothing was ever written to %s\n", m_filename.c_str());
      std::exit(EXIT_FAILURE);
    }

    std::fclose(m_file);
    m_state = STATE_WRITTEN;
  }

  void initialize_reading(long bufsize = (4 << 20)) {
    if (m_state != STATE_WRITTEN) {
      fprintf(stderr, "\nError: initializing reading in state %s\n", state_string().c_str());
      std::exit(EXIT_FAILURE);
    }

    // Compute buffer size.
    m_state = STATE_READING;
    long items = std::max(2UL, (bufsize + sizeof(value_type) - 1) / sizeof(value_type));
    m_buf_size = items / 2L;

    // Reset counters.
    m_active_buf_filled = 0;
    m_passive_buf_filled = 0;
    m_active_buf_pos = 0;
    m_total_read_buf = 0;
    m_total_read_user = 0;
    m_cur_file = -1;

    // Initialize buffers.
    m_active_buf = (value_type *)malloc(m_buf_size * sizeof(value_type));
    m_passive_buf = (value_type *)malloc(m_buf_size * sizeof(value_type));

    // Start the I/O thread and immediatelly start reading.
    m_avail = true;
    m_finished = false;
    m_thread = new std::thread(async_io_code<value_type>, this);
  }

  inline value_type read() {
    if (m_state != STATE_READING) {
      fprintf(stderr, "\nError: reading in state %s\n", state_string().c_str());
      std::exit(EXIT_FAILURE);
    }

    if (m_active_buf_pos == m_active_buf_filled)
      receive_new_buffer();

    m_total_read_user++;
    return m_active_buf[m_active_buf_pos++];
  }

  void finish_reading() {
    if (m_state != STATE_READING) {
      fprintf(stderr, "\nError: finishing reading in state %s\n", state_string().c_str());
      std::exit(EXIT_FAILURE);
    }

    if (m_total_read_buf != m_total_read_user || m_total_read_user != m_total_write) {
      fprintf(stderr, "\nError: not all elems were read from distributed file %s\n", m_filename.c_str());
      std::exit(EXIT_FAILURE);
    }

    // Let the I/O thread know that we are done.
    std::unique_lock<std::mutex> lk(m_mutex);
    m_finished = true;
    lk.unlock();
    m_cv.notify_one();

    // Wait for the thread to finish.
    m_thread->join();

    // Clean up.
    delete m_thread;
    close_and_destroy_cur_file();
    free(m_active_buf);
    free(m_passive_buf);

    // Enter the terminal state.
    m_state = STATE_READ;
  }

  std::string state_string() const {
    switch(m_state) {
      case STATE_INIT:    return "STATE_INIT";
      case STATE_WRITING: return "STATE_WRITING";
      case STATE_WRITTEN: return "STATE_WRITTEN";
      case STATE_READING: return "STATE_READING";
      case STATE_READ:    return "STATE_READ";
      default: return "undefined state";
    }
  }

  void close_and_destroy_cur_file() {
    if (m_state != STATE_READING) {
      fprintf(stderr, "\nError: destroying a file in state %s\n", state_string().c_str());
      std::exit(EXIT_FAILURE);
    }

    if (!m_file) {
      fprintf(stderr, "\nError: deleting a NULL file\n");
      std::exit(EXIT_FAILURE);
    }

    std::fclose(m_file);
    std::string cur_fname = m_filename + ".part" + utils::intToStr(m_cur_file);
    utils::file_delete(cur_fname);
  }

  template<typename T>
  static void async_io_code(distributed_file<T> *file) {
    while (true) {
      // Wait until the passive buffer is available.
      std::unique_lock<std::mutex> lk(file->m_mutex);
      while (!(file->m_avail) && !(file->m_finished))
        file->m_cv.wait(lk);

      if (!(file->m_avail) && (file->m_finished)) {
        // We're done, terminate the thread.
        lk.unlock();
        return;
      }
      lk.unlock();

      // This should never happen.
      if (file->m_total_read_buf == file->m_total_write) {
        fprintf(stderr, "\nError: trying to read past the end of file\n");
        std::exit(EXIT_FAILURE);
      }

      // Safely process the passive buffer.
      // Check if we need to open next file.
      if (file->m_cur_file == -1 || file->m_cur_file_read == file->m_max_items) {
        if (file->m_cur_file != -1)
          file->close_and_destroy_cur_file();
        file->open_next_file();
      }

      // Read the data from disk.
      long file_left = file->m_max_items - file->m_cur_file_read;
      long items_left = file->m_total_write - file->m_total_read_buf;
      long left = std::min(file_left, items_left);
      file->m_passive_buf_filled = std::min(left, file->m_buf_size);
      file->m_cur_file_read += file->m_passive_buf_filled;
      file->m_total_read_buf += file->m_passive_buf_filled;
      utils::read_objects_from_file(file->m_passive_buf,
          file->m_passive_buf_filled, file->m_file);

      // Let the caller know that the I/O thread finished reading.
      lk.lock();
      file->m_avail = false;
      lk.unlock();
      file->m_cv.notify_one();
    }
  }

  void receive_new_buffer() {
    if (m_state != STATE_READING) {
      fprintf(stderr, "\nError: refilling in state %s\n", state_string().c_str());
      std::exit(EXIT_FAILURE);
    }

    // Wait until the I/O thread finishes reading the revious
    // buffer. Most of the time this step is instantaneous.
    std::unique_lock<std::mutex> lk(m_mutex);
    while (m_avail == true)
      m_cv.wait(lk);

    // Set the new active buffer.
    std::swap(m_active_buf, m_passive_buf);
    m_active_buf_filled = m_passive_buf_filled;
    m_active_buf_pos = 0;

    // Let the I/O thead know that it can now
    // prefetch another buffer.
    m_avail = (m_total_read_buf < m_total_write);
    lk.unlock();
    m_cv.notify_one();
  }

  void open_next_file() {
    if (m_state != STATE_READING) {
      fprintf(stderr, "\nError: opening a new file in state %s\n", state_string().c_str());
      std::exit(EXIT_FAILURE);
    }

    ++m_cur_file;
    m_file = utils::open_file(m_filename + ".part" + utils::intToStr(m_cur_file), "r");
    m_cur_file_read = 0;
  }

  void make_new_file() {
    if (m_state != STATE_WRITING) {
      fprintf(stderr, "\nError: making new file in state %s\n", state_string().c_str());
      std::exit(EXIT_FAILURE);
    }

    m_file = utils::open_file(m_filename + ".part" + utils::intToStr(m_files_cnt), "w");
    ++m_files_cnt;
    m_cur_file_write = 0;
  }


  enum { STATE_INIT,    // right after creating (before init_writing)
         STATE_WRITING, // after initialize_writing, writing possible
         STATE_WRITTEN, // after finish_writing, waiting for initialize_reading
         STATE_READING, // after initialize_reading, reading possible
         STATE_READ     // after finish_reading, waiting for death
  } m_state;

  std::FILE *m_file;       // file handler
  std::string m_filename;  // file name base
  long m_max_items;        // max items per file

  // Buffers used for asynchronous reading.
  value_type *m_active_buf;
  value_type *m_passive_buf;
  long m_buf_size;
  long m_active_buf_pos;
  long m_active_buf_filled;
  long m_passive_buf_filled;

  // Various housekeeping statistics about the number of items.
  long m_cur_file_write;   // number of items written to a current file
  long m_total_write;      // total number of written items
  long m_cur_file_read;    // number of items read from the current file
  long m_total_read_buf;   // total number of items read from files into buffers
  long m_total_read_user;  // total number of items read by the user

  // Used to keep track of file count.
  long m_files_cnt; // counts the files during writing
  long m_cur_file;  // iterates through [0..m_files_cnt) during reading

  // For synchronization with thread doing asynchronous reading.
  std::thread *m_thread;
  std::mutex m_mutex;
  std::condition_variable m_cv;
  bool m_finished;
  bool m_avail;
};

}  // psascan_private

#endif // __PSASCAN_SRC_DISTRIBUTED_FILE_H_INCLUDED
