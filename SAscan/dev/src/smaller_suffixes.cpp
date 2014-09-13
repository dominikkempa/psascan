//=============================================================================
// A module implementing the computation of rank necessary to initiate
// parallel streaming in SAscan.
//
// The whole module implements one method, which should be passed to thread
// constructor. The method is:
//
//    void parallel_smaller_suffixes(
//         unsigned char*  block,
//         long            block_size,
//         std::string     text_filename,
//         long            suffix_start_pos,
//         long*           ret
//    );
//
// The method takes the current SAscan block (it is in memory and can be
// random-accessed) and computes the number of suffixes of text starting
// inside that block (and extending until the end of text) that are smaller
// than the suffix of text starting at position suffix_start_pos. The
// result is returned via ret pointer.
//
// In addition, the method assumes that there is a file named
// text_filename.gt on disk containing the gt bitvector from the previous
// iteration of SAscan. This file contains the bitvector defined as follows:
//   
//   gt[i] == 1
//   (i = 0, 1, ...)
//
//     iff
//
//   The suffix of text of length i + 1 is lexicographically larger than
//   the suffix of text starting immediatelly after the 'block' in the text.
//-----------------------------------------------------------------------------
// The function requires random access to the block (which is in memory). In
// addition, it allocates only 3 disk blocks (currently each block is 1MiB) of
// RAM. These blocks are used to access pattern and gt bitvector. Also,
// some negligible space is required for the list of scopes of the pattern,
// which is computed on demand by the 'extend_pattern' method of 'GS_sets'
// class.
//=============================================================================

#include "smaller_suffixes.h"
#include "multifile_bitvector.h"

#include <string>
#include <vector>

#include "utils.h"

struct triple {
  long b, e, c;
  triple(long b_ = 0L, long e_ = 0L, long c_ = 0L)
    : b(b_), e(e_), c(c_) {}

  inline bool operator == (const triple &t) const {
    return b == t.b && e == t.e && c == t.c;
  }
};

struct pair {
  long b, c;
  pair(long b_ = 0L, long c_ = 0L)
    : b(b_), c(c_) {}

  inline bool operator == (const pair &p) const {
    return b == p.b && c == p.c;
  }
};

inline triple contains(std::vector<triple> &S_p, long j) {
  for (long i = 0L; i < (long)S_p.size(); ++i)
    if (S_p[i].b <= j && j < S_p[i].e) return S_p[i];
    else if (S_p[i].b > j) break;
  return triple(0L, 0L, 0L);
}

inline pair pred(std::vector<pair> &S_n, long j) {
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

//==============================================================================
// Computes the number of suffixes starting inside the 'block' (which
// is a substring of text stored in 'text_filename') that are smaller than
// the suffix of text starting at position 'suffix_start_pos'.
//==============================================================================
void parallel_smaller_suffixes(unsigned char *block, long block_beg, long block_end,
    long text_length, std::string text_filename, long suffix_start_pos, long &ret,
    multifile *tail_gt_begin_reversed) {
  long block_size = block_end - block_beg;
  long pat_length = text_length - suffix_start_pos;

  if (!pat_length) {
    ret = 0L;
    return;
  }

  multifile_bitvector_reader gt_reader(tail_gt_begin_reversed);
  pattern pat(text_filename, suffix_start_pos);

  GS_sets GS;
  long count = 0L, i = 0L, L = 0L, maxL = 0L;

  while (i < block_size) {
    while (i + L < block_size && L < pat_length && block[i + L] == pat[L]) ++L;
    if (L > maxL) { maxL = L; GS.extend_pattern(block + i, L); }
    triple bec = contains(GS.S_p, L);

    if (L < pat_length && (
           (i + L  < block_size && block[i + L] < pat[L]) ||
           (i + L == block_size && gt_reader.access((pat_length - block_size) + i))
                          )
       ) ++count;

    if (bec.b != 0L) { count += bec.c; i += bec.b / 2L; L -= bec.b / 2L; }
    else { pair bc = pred(GS.S_n, (L / 3L) + 1L); count += bc.c; i += bc.b; L = 0L; }
  }

  ret = count;
}

