// Simple, 4n bytes rank answering general rank queries.
#ifndef __SIMPLE_RANK_INCLUDED
#define __SIMPLE_RANK_INCLUDED

#include <algorithm>
#include <string>

#include "utils.h"

struct simple_rank {

  simple_rank(unsigned char *b, int size) {
    n = size;
    int pow8 = 1 << 8;
    int pow24 = 1 << 24;
    int mask24 = pow24 - 1;
    int mask8 = pow8 - 1;

    int b_count = ((size + mask8)/ pow8) + 1;
    int sb_count = ((size + mask24) / pow24) + 1;

    bwt32 = new int[b_count * 256];
    if (!bwt32) {
      fprintf(stderr, "Error: cannot allocate bwt32.\n");
      std::exit(EXIT_FAILURE);
    }
    for (int i = 0; i < n; ++i) bwt32[i] = b[i];
    for (int i = n; i < 256 * b_count; ++i) bwt32[i] = 0;

    int *rank_c = new int[256];
    if (!rank_c) {
      fprintf(stderr, "Error: cannot allocate rank.\n");
      std::exit(EXIT_FAILURE);
    }

    sb_rank = new int[256 * sb_count];
    if (!sb_rank) {
      fprintf(stderr, "Error: cannot allocate sb_rank.\n");
      std::exit(EXIT_FAILURE);
    }

    std::fill(rank_c, rank_c + 256, 0);
    int sb_ptr = 0;
    for (int i = 0; i < b_count * 256; ++i) {
      if (!(i & mask24)) {
        for (int j = 0; j < 256; ++j)
          sb_rank[256 * sb_ptr + j] = rank_c[j];
        ++sb_ptr;
      }
      if (!(i & mask8))
        for (int j = 0; j < 256; ++j) {
          int diff = rank_c[j] - sb_rank[(sb_ptr - 1) * 256 + j];
          bwt32[i + j] += (diff << 8);
        }
      rank_c[bwt32[i] & mask8]++;
    }

    delete[] rank_c;
  }

  inline int rank(int i, unsigned char c) const {
    if (i <= 0) return 0;
    if (i > n) i = n;
    static const int mask16 = (1 << 16) - 1;
    static const int mask8 = (1 << 8) - 1;
    int sb_id = (i >> 24);
    int sb_count = sb_rank[256 * sb_id + c];
    int b_id = (i >> 8);
    int b_count = (bwt32[b_id * 256 + c] >> 8);
    int next_b = b_id + 1;
    int nextb_count = 0;
    if (!(next_b & mask16)) nextb_count = sb_rank[(next_b >> 16) * 256 + c] - sb_count;
    else nextb_count = bwt32[256 * next_b + c] >> 8;
    if (nextb_count == b_count) return sb_count + b_count;
    int extra = 0;
    if (i & 128) {
      for (int j = i; j < next_b * 256; ++j)
        if ((bwt32[j] & mask8) == c) ++extra;
      return sb_count + nextb_count - extra;
    } else {
      for (int j = 256 * b_id; j < i; ++j)
        if ((bwt32[j] & mask8) == c) ++extra;
      return sb_count + b_count + extra;
    }
  }

  ~simple_rank() {
    if (bwt32) delete[] bwt32;
    if (sb_rank) delete[] sb_rank;
  }

  int n;
  int *sb_rank;
  int *bwt32;
};

#endif // __SIMPLE_RANK_INCLUDED
