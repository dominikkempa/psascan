//=============================================================================
//
// A module implementing the computation of rank necessary to initiate
// parallel streaming in SAscan.
//
// The whole module implements one method, which should be passed to thread
// constructor. The method is:
//
//    long parallel_smaller_suffixes(
//         unsigned char*  block,
//         long            block_size,
//         std::string     text_filename,
//         long            suffix_start_pos
//    );
//
// The method takes the current SAscan block (it is in memory and can be
// accessed sequentially) and computes the number of suffixes of text starting
// inside that block (and extending until the end of text) that are smaller
// than the suffix of text starting at position suffix_starting_pos.
//
// In addition, the method assumes that there is a file named
// text_filename.gt_head on disk containing the gt bitvector from the previous
// iteration of SAscan. This file contains the bitvector defined as follows:
//   
//   gt[i] == 1
//   (i = 0, 1, ...)
//
//     iff
//
//   The suffix of text of length i + 1 is lexicographically larger than
//   the suffix of text starting immediatelly after the 'block' in the text.
//
//-----------------------------------------------------------------------------
//
// The function requires random access to the block (which is in memory). In
// addition, it allocates only 3 disk blocks (currently each block is 1MiB) of
// disk space. These blocks are used to access pattern and gt bitvector. Also,
// some negligible space is required for the list of scopes of the pattern,
// which is computed on demand by the 'GS_sets' class and its 'extend_pattern'
// method.
//
//=============================================================================

#include "smaller_suffixes.h"

#include <string>
#include <vector>

#include "utils.h"

#define SMALLER_SUFFIXES_DISK_BLOCK_SIZE 1048576L

struct triple {
  long b, e, c;
  triple(long b_ = 0L, long e_ = 0L, long c_ = 0L)
    : b(b_), e(e_), c(c_) {}

  bool operator == (const triple &t) const {
    return b == t.b && e == t.e && c == t.c;
  }
};

struct pair {
  long b, c;
  pair(long b_ = 0L, long c_ = 0L)
    : b(b_), c(c_) {}

  bool operator == (const pair &p) const {
    return b == p.b && c == p.c;
  }
};

triple contains(std::vector<triple> &S_p, long j) {
  for (long i = 0L; i < (long)S_p.size(); ++i)
    if (S_p[i].b <= j && j < S_p[i].e) return S_p[i];
    else if (S_p[i].b > j) break;
  return triple(0L, 0L, 0L);
}

pair pred(std::vector<pair> &S_n, long j) {
  long i = 0;
  while (i + 1 < (long)S_n.size() && S_n[i + 1].b <= j) ++i;
  return S_n[i];
}

struct GS_sets {
  GS_sets()
      : i(1L), L(0L), last(1L), count(0L) {
    S_n.push_back(pair(1, 0));
  }

  void extend_pattern(unsigned char *new_pat, long new_length) {
    while (i < new_length) {
      while (i + L < new_length && new_pat[i + L] == new_pat[L]) ++L;

      triple bec = contains(S_p, L);

      if (3L * i <= i + L && (bec.b == 0L || S_p.back().b == 2L * i)) {
        bec = triple(2L * i, i + L + 1L, count);

        // Modification: we check if we are trying to update the most
        // recently added scope. If yes, we drop the old and add the
        // new (extended) one.
        if (!S_p.empty() && S_p.back().b == 2L * i)
          S_p.pop_back();
        S_p.push_back(bec);
      }

      if (2L * last <= i) {
        S_n.push_back(pair(i, count));
        last = i;
      }

      // Next line is crucial to ensure the correctness of this nethod.
      // Normally, it is not necessary in the Galil-Seiferas, but the
      // elements computed without exiting at this point are not used
      // in the GS matching procedure anyway.
      if (i + L == new_length) return;
      else if (new_pat[i + L] < new_pat[L]) ++count;

      if (bec.b != 0L) { count += bec.c; i += bec.b / 2L; L -= bec.b / 2L; }
      else { pair bc = pred(S_n, (L / 3L) + 1L); count += bc.c; i += bc.b; L = 0L; }
    }
  }

  long i, L, last, count;
  std::vector<triple> S_p;
  std::vector<pair> S_n;
};


//=============================================================================
// A class gt_accessor implements access to gt bitvector on disk. We only
// access gt bitvector sequentially so the class essentially just streams the
// bitvector one disk block after another.
//
// However, it is optimized to perform well in 2 scenarions that often occur
// in practice:
// * only a very small number of bits in gt are accessed. This is why the
//   class does not actually stream the whole file, but whenever the accessed
//   bit is further apart from the bits inside the buffer than one disk block
//   it fseek's to that place in file and then refills the buffer.
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
  long filled;  // the buffer holds 8 * filled bits
  long bit_offset; // index of the first bit in the buffer
};

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

struct pattern {
  pattern(std::string filename, long pat_start) {
    pattern_start = pat_start;
    length = utils::file_size(filename);
    buf = new unsigned char[bufsize];
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
      long neworigin = std::max(0L, i - bufsize / 2L);
      long toread = std::min(length - neworigin, bufsize);
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

  static const long bufsize = 2L * SMALLER_SUFFIXES_DISK_BLOCK_SIZE;

  unsigned char *buf;
  long filled, origin, length, fpos;
  std::FILE *f;
  
  long pattern_start;
};

//==============================================================================
// Computes the number of suffixes starting inside the 'block' (which
// is a substring of text stored in 'text_filename') that are smaller than
// the suffix of text starting at position 'suffix_start_pos'.
//==============================================================================
long parallel_smaller_suffixes(unsigned char *block, long block_size,
    std::string text_filename, long suffix_start_pos) {
  long text_length = utils::file_size(text_filename);
  long pat_length = text_length - suffix_start_pos;
  gt_accessor gt(text_filename + ".gt_tail");
  pattern pat(text_filename, suffix_start_pos);

  GS_sets GS;
  long count = 0L, i = 0L, L = 0L, maxL = 0L;

  while (i < block_size) {
    while (i + L < block_size && L < pat_length && block[i + L] == pat[L]) ++L;
    if (L > maxL) { maxL = L; GS.extend_pattern(block + i, L); }
    triple bec = contains(GS.S_p, L);

    if (L < pat_length && (
           (i + L  < block_size && block[i + L] < pat[L]) ||
           (i + L == block_size && gt[pat_length - L - 1])
                          )
       ) ++count;

    if (bec.b != 0L) { count += bec.c; i += bec.b / 2L; L -= bec.b / 2L; }
    else { pair bc = pred(GS.S_n, (L / 3L) + 1L); count += bc.c; i += bc.b; L = 0L; }
  }

  return count;
}

