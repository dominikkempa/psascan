#ifndef __RANK_H
#define __RANK_H

#include <algorithm>
#include <string>

#include "utils.h"

struct rank_4n {

  rank_4n(unsigned char *b, long size) {
    n = size;
    int pow8 = 1 << 8;
    int pow24 = 1 << 24;
    int mask24 = pow24 - 1;
    int mask8 = pow8 - 1;

    long b_count = ((size + mask8)/ pow8) + 1;
    long sb_count = ((size + mask24) / pow24) + 1;

    bwt32 = new unsigned[b_count * 256];
    if (!bwt32) {
      fprintf(stderr, "Error: cannot allocate bwt32.\n");
      std::exit(EXIT_FAILURE);
    }
    for (long i = 0; i < n; ++i) bwt32[i] = b[i];
    for (long i = n; i < 256 * b_count; ++i) bwt32[i] = 0;

    unsigned *rank_c = new unsigned[256];
    if (!rank_c) {
      fprintf(stderr, "Error: cannot allocate rank.\n");
      std::exit(EXIT_FAILURE);
    }

    sb_rank = new unsigned[256 * sb_count];
    if (!sb_rank) {
      fprintf(stderr, "Error: cannot allocate sb_rank.\n");
      std::exit(EXIT_FAILURE);
    }

    std::fill(rank_c, rank_c + 256, 0);
    int sb_ptr = 0;
    for (long i = 0; i < b_count * 256; ++i) {
      if (!(i & mask24)) {
        for (int j = 0; j < 256; ++j)
          sb_rank[256 * sb_ptr + j] = rank_c[j];
        ++sb_ptr;
      }
      if (!(i & mask8))
        for (int j = 0; j < 256; ++j) {
          unsigned diff = rank_c[j] - sb_rank[(sb_ptr - 1) * 256 + j];
          bwt32[i + j] += (diff << 8);
        }
      rank_c[bwt32[i] & mask8]++;
    }

    delete[] rank_c;
  }

  inline long rank(long i, unsigned char c) const {
    if (i <= 0) return 0;
    if (i > n) i = n;
    static const int mask16 = (1 << 16) - 1;
    static const int mask8 = (1 << 8) - 1;
    unsigned sb_id = (i >> 24);
    long sb_count = sb_rank[(sb_id << 8) + c];
    unsigned b_id = (i >> 8);
    long b_count = (bwt32[(b_id << 8) + c] >> 8);
    unsigned next_b = b_id + 1;
    unsigned nextb_count = 0;    
    if (!(next_b & mask16)) nextb_count = sb_rank[((next_b >> 16) << 8) + c] - sb_count;
    else nextb_count = bwt32[(next_b << 8) + c] >> 8;
    if (nextb_count == b_count) return sb_count + b_count;
    long extra = 0;
    if (i & 128) {
      for (long j = i; j < next_b * 256L; ++j)
        if ((bwt32[j] & mask8) == c) ++extra;
      return (sb_count + nextb_count) - extra;
    } else {
      for (long j = 256L * b_id; j < i; ++j)
        if ((bwt32[j] & mask8) == c) ++extra;
      return sb_count + b_count + extra;
    }
  }

  ~rank_4n() {
    if (bwt32) delete[] bwt32;
    if (sb_rank) delete[] sb_rank;
  }

  long n;
  unsigned *sb_rank;
  unsigned *bwt32;
};

#endif // __RANK_H

