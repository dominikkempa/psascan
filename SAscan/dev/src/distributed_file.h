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

#ifndef __DISTRIBUTED_FILE_INCLUDED
#define __DISTRIBUTED_FILE_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <string>

#include "utils.h"

template<typename data_type>
struct distributed_file {
  // maxsize is the maximal size of a single file (in bytes).
  distributed_file(const char *filename, long maxsize)
    : m_filename(filename),
      m_maxsize(maxsize),
      m_state(STATE_INIT),
      m_elems_per_file(maxsize / sizeof(data_type)) {
    if (!m_elems_per_file) {
      fprintf(stderr, "\nError: maxsize = %ld is too small to hold elems "
          "of size %lu\n", maxsize, sizeof(data_type));
      std::exit(EXIT_FAILURE);
    }
  }

  void initialize_writing(long bufsize = 4 << 20) { // 4MiB default bufsize
    if (m_state != STATE_INIT) {
      fprintf(stderr, "\nError: initializing writing while in state %s\n",
          state_string().c_str());
      std::exit(EXIT_FAILURE);
    }
    m_state = STATE_WRITING;

    m_bufelems = (bufsize + sizeof(data_type) - 1) / sizeof(data_type);
    if (m_bufelems == 0) {
      fprintf(stderr, "\nError: bufsize = %ld is too small to hold elems "
          "of size %lu\n", bufsize, sizeof(data_type));
      std::exit(EXIT_FAILURE);
    }
    m_bufelems = std::min(m_bufelems, m_elems_per_file); // shrink buffer

    m_buf = new data_type[m_bufelems];
    m_filled = 0;

    m_total_written = 0;
    m_files_cnt = 0;
    make_new_file();
  }

  void write(data_type x) {
    if (m_state != STATE_WRITING) {
      fprintf(stderr, "\nError: attempting a single write in state %s\n",
          state_string().c_str());
      std::exit(EXIT_FAILURE);
    }

    if (m_cur_file_written == m_elems_per_file) {
      flush();
      std::fclose(m_file);
      make_new_file();
    } else if (m_filled == m_bufelems)
      flush();

    m_buf[m_filled++] = x;
    ++m_cur_file_written;
    ++m_total_written;
  }

  void write(data_type *begin, data_type *end) {
    if (m_state != STATE_WRITING) {
      fprintf(stderr, "\nError: attempting a block write in state %s\n",
          state_string().c_str());
      std::exit(EXIT_FAILURE);
    }
    m_total_written += (end - begin);

    // First, fill up the current file (if not already full).
    // Optimization: we bypass the buffer.
    flush();
    if (m_cur_file_written != m_elems_per_file) {
      long file_elems_left = m_elems_per_file - m_cur_file_written;
      long elems_to_write = std::min(file_elems_left, end - begin);
      utils::add_objects_to_file(begin, elems_to_write, m_file);
      m_cur_file_written += elems_to_write;
      begin += elems_to_write;
    }

    // Writing remaining elements (possibly into many files).
    while (begin < end) {
      std::fclose(m_file);
      make_new_file();

      long elems_to_write = std::min(m_elems_per_file, end - begin);
      utils::add_objects_to_file(begin, elems_to_write, m_file);
      m_cur_file_written += elems_to_write;
      begin += elems_to_write;
    }
  }

  void finish_writing() {
    if (m_state != STATE_WRITING) {
      fprintf(stderr, "\nError: finishing writing when in state %s\n",
          state_string().c_str());
      std::exit(EXIT_FAILURE);
    }
    if (m_cur_file_written == 0) {
      fprintf(stderr, "\nError: nothing was ever written to %s\n",
          m_filename.c_str());
      std::exit(EXIT_FAILURE);
    }

    flush();
    std::fclose(m_file);
    delete[] m_buf;

    m_state = STATE_WRITTEN;
  }

  void initialize_reading(long bufsize = 4 << 20) { // 4MiB default bufsize
    if (m_state != STATE_WRITTEN) {
      fprintf(stderr, "\nError: initializing reading in state %s\n",
          state_string().c_str());
      std::exit(EXIT_FAILURE);
    }
    m_state = STATE_READING;

    m_bufelems = (bufsize + sizeof(data_type) - 1) / sizeof(data_type);
    if (m_bufelems == 0) {
      fprintf(stderr, "\nError: bufsize %ld is too small to hold elems "
          "of size %lu\n", bufsize, sizeof(data_type));
      std::exit(EXIT_FAILURE);
    }
    m_buf = new data_type[m_bufelems];

    m_total_bufread = m_total_read = 0;
    m_cur_file = -1;
    refill();
  }

  inline data_type read() {
    if (m_state != STATE_READING) {
      fprintf(stderr, "\nError: reading in state %s\n",
          state_string().c_str());
      std::exit(EXIT_FAILURE);
    }
    if (m_bufpos == m_filled)
      refill();

    ++m_total_read;

    return m_buf[m_bufpos++];
  }

  void finish_reading() {
    if (m_state != STATE_READING) {
      fprintf(stderr, "\nError: finishing reading in state %s\n",
          state_string().c_str());
      std::exit(EXIT_FAILURE);
    }

    if (m_total_bufread != m_total_read || m_total_read != m_total_written) {
      fprintf(stderr, "\nError: not all elems were read from "
          "distributed file %s\n", m_filename.c_str());
      std::exit(EXIT_FAILURE);
    }

    close_and_destroy_cur_file();
    delete[] m_buf;
    m_state = STATE_READ;
  }

private:
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
      fprintf(stderr, "\nError: destroying a file in state %s\n",
          state_string().c_str());
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

  void flush() {
    if (m_state != STATE_WRITING) {
      fprintf(stderr, "\nError: flushing in state %s\n",
          state_string().c_str());
      std::exit(EXIT_FAILURE);
    }

    if (m_filled > 0)
      utils::add_objects_to_file(m_buf, m_filled, m_file);

    m_filled = 0;
  }

  void refill() {
    if (m_total_bufread == m_total_written) {
      fprintf(stderr, "\nError: trying to read past the end of file\n");
      std::exit(EXIT_FAILURE);
    }

    if (m_state != STATE_READING) {
      fprintf(stderr, "\nError: refilling in state %s\n",
          state_string().c_str());
      std::exit(EXIT_FAILURE);
    }

    if (m_cur_file == -1 || m_cur_file_read == m_elems_per_file) {
      if (m_cur_file != -1)
        close_and_destroy_cur_file();
      open_next_file();
    }

    long file_left = std::min(m_total_written - m_total_bufread,
        m_elems_per_file - m_cur_file_read);
    m_filled = std::min(file_left, m_bufelems);
    size_t r = std::fread(m_buf, sizeof(data_type), m_filled, m_file);
    if ((long)r != m_filled) {
      fprintf(stderr, "\nError: fread failed during refill()\n");
      std::exit(EXIT_FAILURE);
    }

    m_total_bufread += m_filled;
    m_cur_file_read += m_filled;

    m_bufpos = 0;
  }

  void open_next_file() {
    if (m_state != STATE_READING) {
      fprintf(stderr, "\nError: opening a new file in state %s\n",
          state_string().c_str());
      std::exit(EXIT_FAILURE);
    }

    ++m_cur_file;
    m_file = utils::open_file(m_filename + ".part" + utils::intToStr(m_cur_file), "r");
    m_cur_file_read = 0;
  }

  void make_new_file() {
    if (m_state != STATE_WRITING) {
      fprintf(stderr, "\nError: making new file in state %s\n",
          state_string().c_str());
      std::exit(EXIT_FAILURE);
    }

    m_file = utils::open_file(m_filename + ".part" + utils::intToStr(m_files_cnt), "w");
    ++m_files_cnt;
    m_cur_file_written = 0;
  }

  std::string m_filename; // base of the filename
  long m_maxsize; // maximal size of a single file (in bytes)

  enum { STATE_INIT,    // right after creating (before init_writing)
         STATE_WRITING, // after initialize_writing, writing possible
         STATE_WRITTEN, // after finish_writing, waiting for initialize_reading
         STATE_READING, // after initialize_reading, reading possible
         STATE_READ     // after finish_reading, waiting for death
  } m_state;

  data_type *m_buf;
  long m_bufelems; // maximal number of elemens in a buffer
  long m_filled; // number of elems in a buffer (in reading / writing mode)
  long m_bufpos; // number of elems read from the buffer (in reading mode)

  long m_cur_file_written, m_elems_per_file; // current / max elems written to file
  long m_cur_file_read;                      // current elems read from file
  long m_total_written;     // read / written elems in a sequence

  long m_total_bufread; // total elems read from files (into buffers)
  long m_total_read; // total elems read by the user with `read' function

  long m_files_cnt; // counts the files during writing
  long m_cur_file;  // iterates through [0..m_files_cnt) during reading
  std::FILE *m_file;
};

#endif // __DISTRIBUTED_FILE_INCLUDED

