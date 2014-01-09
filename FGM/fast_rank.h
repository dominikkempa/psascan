#ifndef __FAST_RANK_H
#define __FAST_RANK_H

#include <algorithm>
#include <vector>
#include <string>

#include "utils.h"

struct fast_sigma_rank_4n {
  static const int k_sb_size_bits = 24; // log_2 of the super block size

  fast_sigma_rank_4n(unsigned char *text, long length, int sigma_bits)
      : m_sigma_bits(sigma_bits), m_sigma(1 << m_sigma_bits),
        m_blocks_in_sb(k_sb_size / m_sigma), m_length(length) {
    long b_count = (length + m_sigma - 1) / m_sigma + 1;
    long sb_count = (b_count + m_blocks_in_sb - 1) / m_blocks_in_sb;
    trunk = new unsigned[b_count * m_sigma];
    std::fill(trunk, trunk + b_count * m_sigma, 0);
    for (long i = 0; i < length; ++i) trunk[i] = text[i];
    sb_rank = new unsigned[m_sigma * sb_count];
    count = new unsigned[m_sigma];
    std::fill(count, count + m_sigma, 0);
    for (long i = 0, sb_ptr = 0; i < b_count * m_sigma; ++i) {
      if (i % k_sb_size == 0) for (int j = 0; j < m_sigma; ++j) sb_rank[sb_ptr++] = count[j];
      if (i % m_sigma == 0) for (int j = 0; j < m_sigma; ++j)
        trunk[i + j] += (count[j] - sb_rank[sb_ptr - m_sigma + j]) << 8;
      count[trunk[i] & 255]++;
    }
    count[0] -= b_count * m_sigma - m_length;
  }

  inline long rank(long i, unsigned char c) const {
    if (c >= m_sigma || i <= 0) return 0;
    else if (i >= m_length) return count[c];

    unsigned b_id = i >> m_sigma_bits, nextb_id = b_id + 1;
    unsigned sb_id = (i >> k_sb_size_bits);
    long sb_count = sb_rank[(sb_id << m_sigma_bits) + c];
    long b_count = (trunk[(b_id << m_sigma_bits) + c] >> 8);
    long nextb_count = 0;
    if (!(nextb_id & (m_blocks_in_sb - 1)))
      nextb_count = sb_rank[((sb_id + 1) << m_sigma_bits) + c] - sb_count;
    else nextb_count = trunk[(nextb_id << m_sigma_bits) + c] >> 8;

    // optimization
    if (!((nextb_count - b_count) & (m_sigma - 1)))
      return sb_count + b_count + ((b_count == nextb_count) ? 0 : (i & (m_sigma - 1)));
 
    long extra = 0;
    if (i & 128) {
      for (unsigned j = i; j < nextb_id * m_sigma; ++j) if ((trunk[j] & 255) == c) ++extra;
      return sb_count + nextb_count - extra;
    } else {
      for (unsigned j = m_sigma * b_id; j < i; ++j) if ((trunk[j] & 255) == c) ++extra;
      return sb_count + b_count + extra;
    }
  }

  ~fast_sigma_rank_4n() {
    delete[] trunk;
    delete[] sb_rank;
    delete[] count;
  }

private:
  static const int k_sb_size = (1 << k_sb_size_bits);
  static const int k_sb_size_mask = k_sb_size - 1;
  const long m_sigma_bits, m_sigma, m_blocks_in_sb;

  long m_length;
  unsigned *sb_rank, *count, *trunk;
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
      rare_rank = new fast_sigma_rank_4n(rbwt, rare_sum, utils::log2ceil(rare_count));
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
    freq_rank = new fast_sigma_rank_4n(text, length, utils::log2ceil(new_sigma));
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

