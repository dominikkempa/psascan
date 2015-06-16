#ifndef __BACKGROUND_CHUNK_READER_H_INCLUDED
#define __BACKGROUND_CHUNK_READER_H_INCLUDED

#include <cstdio>
#include <string>

#include "utils.h"


struct background_chunk_reader {
  private:
    std::FILE *m_file;
    long m_chunk_length;
    long m_last_wait_argument;
    long m_end;

  public:
    unsigned char *m_chunk;

  public:
    background_chunk_reader(std::string filename, long beg,
        long end, long chunk_length = (1L << 20)) {
      m_chunk_length = chunk_length;
      m_file = utils::open_file(filename, "r");
      m_chunk = (unsigned char *)malloc(chunk_length);
      m_last_wait_argument = beg;
      m_end = end;
    }

    void wait(long end) {
      long length = end - m_last_wait_argument;
      utils::read_block(m_file, m_last_wait_argument, length, m_chunk);
      m_last_wait_argument = end;
    }

    ~background_chunk_reader() {
      std::fclose(m_file);
      free(m_chunk);
    }

    inline long get_chunk_size() const {
      return m_chunk_length;
    }
};

#endif  // __BACKGROUND_CHUNK_READER_H_INCLUDED
