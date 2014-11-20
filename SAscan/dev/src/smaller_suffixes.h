#ifndef __SMALLER_SUFFIXES_H_INCLUDED
#define __SMALLER_SUFFIXES_H_INCLUDED

#include <cstdio>
#include <string>

#include "utils.h"
#include "multifile.h"
#include "disk_pattern.h"


#if 0
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
#endif


//=============================================================================
// A class gt_accessor implements access to gt bitvector on disk. We only
// access gt bitvector sequentially so the class essentially just streams the
// bitvector one disk block after another.
//
// However, it is optimized to perform well in 2 scenarios that often occur
// in practice:
// * only a very small number of bits in gt are accessed. This is why the
//   class does not actually stream the whole file, but whenever the accessed
//   bit is further apart from the bits inside the buffer than one disk block,
//   it fseeks to that place in file and then refills the buffer.
// * almost all bits are accessed. In this case the access is also efficient
//   because we read the data in disk blocks.
//=============================================================================

// Debug implementation of gt_accessor.
/*struct gt_accessor {
  gt_accessor(std::string filename) {
    long filesize;
    utils::read_objects_from_file<unsigned char>(data, filesize, filename);
  }

  ~gt_accessor() {
    delete[] data;
  }

  inline bool operator[] (int i) const {
    long byte_idx = i >> 3;
    long bit_idx = i & 7;
    
    return data[byte_idx] & (1 << bit_idx);
  }

  unsigned char *data;
};*/

// Efficient implementation of gt_accessor.
struct gt_accessor {
  gt_accessor(std::string filename) {
    f = utils::open_file(filename, "r");
    buf = new unsigned char[bufsize];
    filled = 0;
    bit_offset = 0;

    // NOTE: if gt is not accessed at all, we
    // do not read even a single block.
  }

  ~gt_accessor() {
    std::fclose(f);
    delete[] buf;
  }

  //---------------------------------------------------------------------------
  // Get the i-th bit.
  //
  // Range of bits stored in a buffer is [bit_offset...bit_offset + 8*filled).
  //---------------------------------------------------------------------------
  inline bool operator[] (long i) {
    if (!filled) {

      // First read.
      reposition_buffer(i);
    } else if (i >= bit_offset + 8L * filled) {

      // Requested bit not in a buffer.
      if (i < bit_offset + 8L * (filled + bufsize)) {

        // Requested bit is inside the next block: refill buffer, no fseek.
        bit_offset += 8L * filled;
        filled = std::fread(buf, 1, bufsize, f);
      } else {

        // Requested bit is more than one block apart: fseek.
        reposition_buffer(i);
      }
    }

    i -= bit_offset;
    return buf[i >> 3] & (1 << (i & 7));
  }

  void reposition_buffer(long bit_id) {
    long byte_offset = bit_id >> 3;
    bit_offset = byte_offset << 3;

    std::fseek(f, byte_offset, SEEK_SET);
    filled = std::fread(buf, 1, bufsize, f);
  }

  std::FILE *f;
  unsigned char *buf;
  static const long bufsize = SMALLER_SUFFIXES_DISK_BLOCK_SIZE;
  long filled;     // the buffer holds 8 * filled bits
  long bit_offset; // index of the first bit in the buffer
};


void parallel_smaller_suffixes(unsigned char *block, long block_beg, long block_end,
    long text_length, std::string text_filename, long suffix_start_pos, long &ret, multifile *tail_gt_begin);

// Return true iff text[i..) (but we always stop the comparison at text_length)
// is smaller than pat[0..pat_length).

bool lcp_compare2(unsigned char *text, long text_length, long text_absolute_end, pattern &pat, long pat_length, long pat_absolute_beg,
    long supertext_length, long j, std::string supertext_filename);

template<typename saidx_t>
void parallel_smaller_suffixes2(
    unsigned char *block,
    long block_beg,
    long block_end,
    long text_length,
    saidx_t *block_partial_sa,
    std::string text_filename,
    long suf_start,
    long &ret) {

  long block_size = block_end - block_beg;
  long pat_length = text_length - suf_start;

  if (!pat_length) {
    ret = 0L;
    return;
  }

  pattern pat(text_filename, suf_start);

  long left = -1L;
  long right = block_size;

  while (left + 1 != right) {
    long mid = (left + right) / 2;

    if (lcp_compare2(block, block_size, block_end, pat, pat_length, suf_start, text_length, block_partial_sa[mid], text_filename))
      right = mid;
    else left = mid;
  }

  ret = right;
}


#endif // __SMALLER_SUFFIXES_H_INCLUDED
