#ifndef __FAST_RANK_H
#define __FAST_RANK_H

#include <algorithm>
#include <vector>
#include <string>

#include "utils.h"

struct fast_sigma_rank_4n {

  fast_sigma_rank_4n(unsigned char *text, long length, int sigma)
      : m_length(length), m_sigma(sigma) {
    sb_rate = sigma * (1 << 16);

    long b_count = ((length + sigma - 1) / sigma) + 1;
    long sb_count = ((length + sb_rate - 1) / sb_rate) + 1;

    bwt32 = new unsigned[b_count * sigma];
    if (!bwt32) {
      fprintf(stderr, "Error: cannot allocate bwt32.\n");
      std::exit(EXIT_FAILURE);
    }
    for (long i = 0; i < length; ++i) bwt32[i] = text[i];
    for (long i = length; i < b_count * sigma; ++i) bwt32[i] = 0;

    rank_c = new unsigned[sigma];
    if (!rank_c) {
      fprintf(stderr, "Error: cannot allocate rank.\n");
      std::exit(EXIT_FAILURE);
    }

    sb_rank = new unsigned[sigma * sb_count];
    if (!sb_rank) {
      fprintf(stderr, "Error: cannot allocate sb_rank.\n");
      std::exit(EXIT_FAILURE);
    }

    std::fill(rank_c, rank_c + sigma, 0);
    int sb_ptr = 0;
    for (long i = 0; i < b_count * sigma; ++i) {
      if (!(i % sb_rate)) {
        for (int j = 0; j < sigma; ++j)
          sb_rank[sigma * sb_ptr + j] = rank_c[j];
        ++sb_ptr;
      }
      if (!(i % sigma))
        for (int j = 0; j < sigma; ++j) {
          unsigned diff = rank_c[j] - sb_rank[(sb_ptr - 1) * sigma + j];
          bwt32[i + j] += (diff << 8);
        }
      rank_c[bwt32[i] & 255]++;
    }
    rank_c[0] -= b_count * sigma - length;
  }

  inline long rank(long i, unsigned char c) const {
    if (c >= m_sigma) return 0;
    if (i <= 0) return 0;
    if (i >= m_length) return rank_c[c];
    unsigned sb_id = i / sb_rate;
    long sb_count = sb_rank[m_sigma * sb_id + c];
    unsigned b_id = i / m_sigma;
    long b_count = (bwt32[b_id * m_sigma + c] >> 8);
    unsigned next_b = b_id + 1;
    unsigned nextb_count = 0;

    if (!(next_b & mask16)) nextb_count = sb_rank[(next_b >> 16) * m_sigma + c] - sb_count;
    else nextb_count = bwt32[m_sigma * next_b + c] >> 8;
    if (nextb_count == b_count) return sb_count + b_count;

    long extra = 0;
    if (i & 128) {
      for (unsigned j = i; j < next_b * m_sigma; ++j)
        if ((bwt32[j] & mask8) == c) ++extra;
      return (sb_count + nextb_count) - extra;
    } else {
      for (unsigned j = m_sigma * b_id; j < i; ++j)
        if ((bwt32[j] & mask8) == c) ++extra;
      return sb_count + b_count + extra;
    }
  }

  ~fast_sigma_rank_4n() {
    if (bwt32) delete[] bwt32;
    if (sb_rank) delete[] sb_rank;
    if (rank_c) delete[] rank_c;
  }

  static const int mask8 = (1 << 8) - 1;
  static const int mask16 = (1 << 16) - 1;

  long m_length;
  int m_sigma, sb_rate;
  unsigned *sb_rank, *rank_c;
  unsigned *bwt32;
};

struct fast_rank_4n {
  fast_rank_4n(unsigned char *text, long length)
      : m_length(length), count(NULL), is_freq(NULL),
      mapping(NULL), freq_rank(NULL), rare_rank(NULL) {
    if (!length) return;
    count = new long[256];
    std::fill(count, count + 256, 0);
    for (long i = 0; i < length; ++i) ++count[text[i]];

    std::vector<std::pair<unsigned, int> > sorted;
    for (int i = 0; i < 256; ++i)
      if (count[i] > 0) sorted.push_back(std::make_pair(count[i], i));
    std::sort(sorted.begin(), sorted.end());
    std::reverse(sorted.begin(), sorted.end());
    
    // Separate chars into frequent and rare symbol.
    is_freq = new bool[256];
    mapping = new int[256];

    // First, build the rank for rare symbols.
    std::fill(is_freq, is_freq + 256, true);
    std::fill(mapping, mapping + 256, -1);
    unsigned char rare_count = 0;
    long rare_sum = 0;
    while (32L * (rare_sum + sorted[rare_count].first) <= length) {
      is_freq[sorted[rare_count].second] = false;
      rare_sum += sorted[rare_count].first;
      mapping[sorted[rare_count].second] = rare_count;
      ++rare_count;
    }
    if (rare_count > 0) {    
      unsigned char *rbwt = new unsigned char[rare_sum];
      for (int i = 0, ptr = 0; i < length; ++i)
        if (!is_freq[text[i]]) rbwt[ptr++] = mapping[text[i]];
      rare_rank = new fast_sigma_rank_4n(rbwt, rare_sum, rare_count);
      delete[] rbwt;
    } else rare_rank = NULL;

    // Build the rank over frequent chars.
    freq_count = 0;
    for (unsigned i = rare_count; i < sorted.size(); ++i)
      mapping[sorted[i].second] = freq_count++;
    for (int i = 0; i < length; ++i)
      if (is_freq[text[i]]) text[i] = mapping[text[i]];
      else text[i] = freq_count;
    int new_sigma = freq_count + (rare_count > 0);
    freq_rank = new fast_sigma_rank_4n(text, length, new_sigma);
  }
  
  inline long rank(long i, unsigned char c) const {
    if (!m_length || i <= 0 || mapping[c] == -1) return 0;
    if (i >= m_length) return count[c];
    if (is_freq[c]) return freq_rank->rank(i, mapping[c]);
    else return rare_rank->rank(freq_rank->rank(i, freq_count), mapping[c]);
  }

  ~fast_rank_4n() {
    delete[] is_freq;
    delete[] mapping;
    delete[] count;
    if (freq_rank) delete freq_rank;
    if (rare_rank) delete rare_rank;
  }

  long m_length;
  int freq_count;
  long *count;
  bool *is_freq;
  int *mapping;
  fast_sigma_rank_4n *freq_rank;
  fast_sigma_rank_4n *rare_rank;
};

#endif // __FAST_RANK_H

