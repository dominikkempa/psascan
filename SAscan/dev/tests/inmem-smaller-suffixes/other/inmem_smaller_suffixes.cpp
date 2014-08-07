//==============================================================================
// A module implementing the in-memory computation of rank necessary to
// initiate parallel streaming in the in-memory SAscan.
//
// The whole module implements a method:
//
//    void inmem_smaller_suffixes(
//         unsigned char*   text,
//         long             text_length,
//         long             block_beg,
//         long             block_end,
//         long             suf_start
//         bitvector*       gt,
//         long*            res
//    );
//
// By block we mean text[block_beg..block_end). The function computes the
// number of suffixes starting inside the block (and extending until
// the end of the string) that are smaller than the suffix of text
// text[suf_start..text_length). The result is returned via res.
//
// The method also uses bitvector gt where for i = 0, 1, ... gt[i] == 1
// iff the suffix of text starting at position block_end + i is
// lexicographically larger than text[block_size..text_length).
//==============================================================================

#include "inmem_smaller_suffixes.h"

#include <vector>

#include "bitvector.h"

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
      : i(1L), match(0L), last(1L), count(0L) {
    S_n.push_back(pair(1L, 0L));
  }

  void extend_pattern(unsigned char *new_pat, long new_length) {
    while (i < new_length) {
      while (i + match < new_length && new_pat[i + match] == new_pat[match])
        ++match;

      triple bec = contains(S_p, match);

      if (3L * i <= i + match && (bec.b == 0L || S_p.back().b == 2L * i)) {
        bec = triple(2L * i, i + match + 1L, count);

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
      if (i + match == new_length) return;
      else if (new_pat[i + match] < new_pat[match]) ++count;

      if (bec.b) {
        count += bec.c;
        i += bec.b / 2L;
        match -= bec.b / 2L;
      } else {
        pair bc = pred(S_n, match / 3L + 1L);
        count += bc.c;
        i += bc.b;
        match = 0L; 
      }
    }
  }

  long i, match, last, count;
  std::vector<triple> S_p;
  std::vector<pair> S_n;
};

//==============================================================================
// Compute the number of suffixes starting inside a given block of positions
// that are smaller than some other specified suffix.
//
// More precisely: compute the number of suffixes starting inside
// [block_beg..blocke_end) (and extending untill the end of the string)
// that are smaller than the suffix of text starting at position suf_start.
//
// To speed up the computation we are given a bitvector gt. For i = 0, 1, ...
// gt[i] = 1 iff text[range_end + i..text_length) > text[block_end..block_length).
//
// The bitvector is guaranteed to be long enough.
//------------------------------------------------------------------------------
// We assume 0 <= block_beg < block_end <= suf_start <= text_length.
//==============================================================================
void inmem_smaller_suffixes(unsigned char *text, long text_length,
    long block_beg, long block_end, long suf_start, bitvector *gt, long &res) {
  res = 0L;

  unsigned char *block = text + block_beg;
  long block_size = block_end - block_beg;

  unsigned char *pat = text + suf_start;
  long pat_length = text_length - suf_start;
  if (!pat_length) return;

  GS_sets GS;

  long i = 0L;
  long match = 0L, max_match = 0L;
  while (i < block_size) {
    while (i + match < block_size && match < pat_length &&
        block[i + match] == pat[match]) ++match;

    if (match > max_match) {
      max_match = match;
      GS.extend_pattern(block + i, match);
    }

    triple bec = contains(GS.S_p, match);
    if (match < pat_length &&
        ((i + match  < block_size && block[i + match] < pat[match]) ||
         (i + match == block_size && gt->get(suf_start + match))))
      ++res;

    fprintf(stderr, "  pat_length = %ld\n", pat_length);
    fprintf(stderr, "  suf_start = %ld\n", suf_start);
    fprintf(stderr, "  i = %ld, match = %ld\n", i, match);
    fprintf(stderr, "  block_size = %ld, block[i + match = %ld] = %ld, pat[match = %ld] = %ld\n",
        block_size, i + match, (long)block[i + match], match, (long)pat[match]);

    if (bec.b) {
      res += bec.c;
      i += bec.b / 2L;
      match -= bec.b / 2L;
    } else {
      pair bc = pred(GS.S_n, match / 3L + 1L);
      res += bc.c;
      i += bc.b;
      match = 0L;
    }
  }
}

