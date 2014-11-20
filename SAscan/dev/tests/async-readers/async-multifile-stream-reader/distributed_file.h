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

#ifndef __DISTRIBUTED_FILE_H_INCLUDED
#define __DISTRIBUTED_FILE_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <string>
#include <thread>
#include <mutex>
#include <vector>
#include <algorithm>
#include <condition_variable>

#include "utils.h"

struct single_file_info {
  long m_beg;
  long m_end;
  std::string m_filename;

  single_file_info(long beg, long end, std::string filename) {
    m_beg = beg;
    m_end = end;
    m_filename = filename;
  }
};

struct multifile {
  std::vector<single_file_info> files_info;

  void add_file(long beg, long end, std::string filename) {
    files_info.push_back(single_file_info(beg, end, filename));
  }

  ~multifile() {
    for (size_t i = 0; i < files_info.size(); ++i)
      utils::file_delete(files_info[i].m_filename);
  }
};

struct distributed_file {
  distributed_file(multifile *m, long bufsize = (4L << 20)) {
    m_files_info = m->files_info;
  
    long items = std::max(2L, bufsize);
    m_buf_size = items / 2L;

    // Reset counters.
    m_active_buf_filled = 0;
    m_passive_buf_filled = 0;
    m_active_buf_pos = 0;    
    m_total_read_buf = 0;
    m_file = NULL;

    // Initialize buffers.
    m_active_buf = (unsigned char *)malloc(m_buf_size);
    m_passive_buf = (unsigned char *)malloc(m_buf_size);

    // Start the I/O thread and immediatelly start reading.
    m_avail = true;
    m_finished = false;
    m_thread = new std::thread(async_io_code, this);
  }

  inline unsigned char read() {
    if (m_active_buf_pos == m_active_buf_filled)
      receive_new_buffer();

    return m_active_buf[m_active_buf_pos++];
  }

  ~distributed_file() {
    // Let the I/O thread know that we are done.
    std::unique_lock<std::mutex> lk(m_mutex);
    m_finished = true;
    lk.unlock();
    m_cv.notify_one();

    // Wait for the thread to finish.
    m_thread->join();

    // Clean up.
    delete m_thread;
    free(m_active_buf);
    free(m_passive_buf);
    if (m_file)
      std::fclose(m_file);
  }

  static void async_io_code(distributed_file *file) {
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
      
      if (file->m_file == NULL) {
        // Find the next file to open.
        for (size_t j = 0; j < file->m_files_info.size(); ++j)
          if (file->m_files_info[j].m_beg == file->m_total_read_buf) {
            file->m_file_id = j;
            file->m_file = utils::open_file(file->m_files_info[j].m_filename, "r");
            break;
          }
      }

      // If file ID was found, we perform the read.
      // Otherwise there is no more data to prefetch.
      if (file->m_file != NULL) {
        long file_left = file->m_files_info[file->m_file_id].m_end - file->m_total_read_buf;
        file->m_passive_buf_filled = std::min(file_left, file->m_buf_size);
        utils::read_objects_from_file(file->m_passive_buf, file->m_passive_buf_filled, file->m_file);
        file->m_total_read_buf += file->m_passive_buf_filled;
        if (file->m_total_read_buf == file->m_files_info[file->m_file_id].m_end) {
          std::fclose(file->m_file);
          file->m_file = NULL;
        }
      }

      // Let the caller know that the I/O thread finished reading.
      lk.lock();
      file->m_avail = false;
      lk.unlock();
      file->m_cv.notify_one();
    }
  }

  void receive_new_buffer() {
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
    m_avail = true;
    lk.unlock();
    m_cv.notify_one();
  }

private:
  std::FILE *m_file;       // file handler  
  long m_total_read_buf;   // total number of items read from files into buffers
  long m_file_id;
  std::vector<single_file_info> m_files_info;

  // Buffers used for asynchronous reading.
  unsigned char *m_active_buf;
  unsigned char *m_passive_buf;
  long m_buf_size;
  long m_active_buf_pos;
  long m_active_buf_filled;
  long m_passive_buf_filled;

  // For synchronization with thread doing asynchronous reading.
  std::thread *m_thread;
  std::mutex m_mutex;
  std::condition_variable m_cv;
  bool m_finished;
  bool m_avail;
};

#endif // __DISTRIBUTED_FILE_H_INCLUDED
