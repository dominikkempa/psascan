#ifndef __SMALLER_SUFFIXES_H_INCLUDED
#define __SMALLER_SUFFIXES_H_INCLUDED

#include <cstdio>
#include <string>

#include "utils.h"
#include "multifile_bitvector.h"

#define SMALLER_SUFFIXES_DISK_BLOCK_SIZE 1048576L

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


void parallel_smaller_suffixes(unsigned char *block, long block_size,
    std::string text_filename, long suffix_start_pos, long &ret, multifile &gt_files);

#endif // __SMALLER_SUFFIXES_H_INCLUDED
