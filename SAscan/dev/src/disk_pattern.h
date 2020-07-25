#ifndef __DISK_PATTERN_H_INCLUDED
#define __DISK_PATTERN_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <algorithm>

#include "utils.h"

//=============================================================================
// A class pattern implements I/O-efficient random access to pattern (which is
// a substring of text stored in 'filename' starting at position 'pat_start')
// using 0-based indexing, e.g., pattern[0] is the pat_start-th byte of text.
//=============================================================================

// Debug implementation of pattern class.
/*struct pattern {
  pattern(std::string filename, long pat_start) {
    long filesize = utils::file_size(filename);
    long pat_length = filesize - pat_start;
    data = new unsigned char[pat_length];

    // Read the whole pattern at once and keep in memory.
    utils::read_block(filename, pat_start, pat_length, data);
  }

  inline unsigned char operator[] (int i) const {
    return data[i];
  }
  
  ~pattern() {
    delete[] data;
  }
  
  unsigned char *data;
};*/

#define SMALLER_SUFFIXES_DISK_BLOCK_SIZE 1048576L
#define pattern_bufsize (2L * SMALLER_SUFFIXES_DISK_BLOCK_SIZE)

struct pattern {
  pattern(std::string filename, long pat_start) {
    pattern_start = pat_start;
    length = utils::file_size(filename);
    buf = new unsigned char[pattern_bufsize];
    origin = filled = fpos = 0L;
    f = utils::open_file(filename, "r");
  }
  
  ~pattern() {
    std::fclose(f);
    delete[] buf;
  }

  inline unsigned char operator[] (long i) {
    return get_absolute(pattern_start + i);
  }

private:
  inline unsigned char get_absolute(long i) {
    if (i >= origin && i < origin + filled)  return buf[i - origin];
    else {
      // Refill buffer centered at position i.
      long neworigin = std::max(0L, i - pattern_bufsize / 2L);
      long toread = std::min(length - neworigin, pattern_bufsize);
      long roffset = 0L;
      filled = 0L;

      // Save some disk I/O if possible.
      if (i < origin && neworigin + toread > origin) {
        long overlap = neworigin + toread - origin;
        for (long j = 1L; j <= overlap; ++j)
          buf[toread - j] = buf[overlap - j];
        filled = overlap;
        toread -= overlap;
      } else if (i > origin + filled && neworigin < origin + filled) {
        long overlap = origin + filled - neworigin;
        for (long j = 0L; j < overlap; ++j)
          buf[j] = buf[filled - overlap + j];
        filled = overlap;
        toread -= overlap;
        roffset = overlap;
      } else filled = 0L;

      // Do the disk I/O.
      long offset = neworigin + roffset - fpos;
      if (offset && std::fseek(f, offset, SEEK_CUR)) {
        std::perror("Error: pattern fseek1 failed.\n");
        std::exit(EXIT_FAILURE);
      }
      if (std::ftell(f) != neworigin) {
        std::perror("Error: incorrect pattern neworigin.\n");
        std::exit(EXIT_FAILURE);
      }
      long r = std::fread(buf + roffset, 1, toread, f);
      if (r != toread) {
        std::perror("Error: pattern fread1 failed.\n");
        std::exit(EXIT_FAILURE);
      }
      filled += toread;
      fpos = std::ftell(f);
      origin = neworigin;

      return buf[i - origin];
    }
  }

  unsigned char *buf;
  long filled, origin, length, fpos;
  std::FILE *f;
  
  long pattern_start;
};

#endif  // __DISK_PATTERN_H_INCLUDED

