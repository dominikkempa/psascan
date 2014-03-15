#ifndef __RANK_H
#define __RANK_H

#include <algorithm>
#include <vector>

#include "utils.h"

#define FREQUENT_CHAR 0
#define RARE_CHAR 1
#define NOT_OCCURRING 2

struct context_rank_4n {
  static const int k_sb_size_bits = 24;    // log_2(super block size), 1MiB
  static const int k_block_size_bits = 18; // log 2(block size), 64KiB
                                           // anything 16-20 is a pretty good value.

  context_rank_4n(unsigned char *text, long length)
      : m_length(length),
        n_block((length + k_block_size - 1) / k_block_size),
        n_sblock((n_block + k_blocks_in_sb - 1) / k_blocks_in_sb) {

    c_rank       = new long[256];
    block_header = new long[n_block];
    sb_rank      = new long[256 * n_sblock];

    m_mapping    = new unsigned char[2 * 256 * n_block];
    freq_trunk   = new unsigned[n_block * k_block_size];

    // debug ///////
    fprintf(stderr, "  m_mapping takes %.2Lf MiB\n", (long double)(512 * n_block) / (1 << 20));
    fprintf(stderr, "  k_sb_size_bits = %d, k_block_size_bits = %d\n", k_sb_size_bits, k_block_size_bits);
    freq_cnt_total = rare_cnt_total = 0;
    ////////////////

    std::fill(c_rank, c_rank + 256, 0);
    for (long start = 0, sb_ptr = 0; start < n_block * k_block_size; start += k_block_size) {
      long end = start + k_block_size; // text[start..end) is current block
      long block_id = start / k_block_size;

      // 1. Sort symbols in the current block by frequency.
      unsigned count[256] = {0};
      for (long i = start; i < end; ++i) ++count[(i < length ? text[i] : 0)];
      std::vector<std::pair<unsigned, unsigned> > sorted_chars;
      for (int i = 0; i < 256; ++i) if (count[i])
        sorted_chars.push_back(std::make_pair(count[i], i));
      std::sort(sorted_chars.begin(), sorted_chars.end());

      // 2. Separate (at most, due to rounding of freq_cnt) ~3% of rarest symbols.
      long rare_cnt = 0, rare_sum = 0;
      while (rare_cnt < (int)sorted_chars.size() && 16 * (rare_sum + sorted_chars[rare_cnt].first) <= k_block_size)
        rare_sum += sorted_chars[rare_cnt++].first;
      long freq_cnt = (int)sorted_chars.size() - rare_cnt; // Compute freq_cnt. Then round up freq_cnt + 1 (+1 is
      long freq_cnt_bits = utils::log2ceil(freq_cnt + 1);  // for rare char marker) to the smallest power of two.
      freq_cnt = (1 << freq_cnt_bits);                     // NOTE: this does not work when freq_cnt == 256.
      rare_cnt = std::max(0L, (long)sorted_chars.size() - freq_cnt + 1); // Recompute rare_cnt (note the +1)
      std::vector<unsigned char> freq_chars, rare_chars;
      for (int i = 0; i < rare_cnt; ++i) {
        rare_chars.push_back(sorted_chars[i].second);
        // debug ///////
        rare_cnt_total += sorted_chars[i].first;
        ////////////////
      }
      for (int i = rare_cnt; i < (int)sorted_chars.size(); ++i) {
        freq_chars.push_back(sorted_chars[i].second);
        // debug //////
        freq_cnt_total += sorted_chars[i].first;
        ///////////////
      }
      long rare_cnt_bits = 0;
      if (rare_cnt) {                              // If there are rare symbols,
        rare_cnt_bits = utils::log2ceil(rare_cnt); // round up rare_cnt to the
        rare_cnt = (1 << rare_cnt_bits);           // smallest power of two.
      }

      // 3. Update block header and mapping.
      block_header[block_id] = (rare_trunk.size() << 16) + // 6 bytes
                               (rare_cnt_bits << 8) +      // one byte
                               freq_cnt_bits;              // one byte

      // 4. Compute and store symbols mapping.
      bool is_freq[256] = {false};
      unsigned char freq_map[256], rare_map[256];
      for (unsigned c = 0; c < 256; ++c)
        m_mapping[2 * (c * n_block + block_id)] = NOT_OCCURRING;
      for (unsigned i = 0; i < freq_chars.size(); ++i) {
        unsigned char c = freq_chars[i];
        is_freq[c] = true;
        m_mapping[2 * (c * n_block + block_id) + 1] = freq_map[c] = i;
        m_mapping[2 * (c * n_block + block_id)] = FREQUENT_CHAR;
      }
      for (unsigned i = 0; i < rare_chars.size(); ++i) {
        unsigned char c = rare_chars[i];
        freq_map[c] = freq_cnt - 1; // Rare symbol marker.
        m_mapping[2 * (c * n_block + block_id) + 1] = rare_map[c] = i;
        m_mapping[2 * (c * n_block + block_id)] = RARE_CHAR;
      }

      // 5. Update rank at the super block boundary.
      if (start % k_sb_size == 0)
        for (int j = 0; j < 256; ++j)
          sb_rank[sb_ptr++] = c_rank[j];

      // 6. Compute trunk.
      unsigned rare_trunk_ptr = rare_trunk.size(); 
      for (long i = start, nofreq_cnt = 0; i < end; ++i) {
        unsigned char c = (i < length ? text[i] : 0);
        if (i % freq_cnt == 0) {
          for (int j = 0; j < freq_cnt - 1; ++j) {
            unsigned char freq_ch = (j < (int)freq_chars.size() ? freq_chars[j] : 0);
            freq_trunk[i + j] = (c_rank[freq_ch] - sb_rank[sb_ptr - 256 + freq_ch]) << 8;
          }
          freq_trunk[i + freq_cnt - 1] = (nofreq_cnt << 8);
        }
        freq_trunk[i] |= freq_map[c]; // Mapping of frequent c.
        if (!is_freq[c]) {
          if (nofreq_cnt % rare_cnt == 0)
            for (int j = 0; j < rare_cnt; ++j) {
              unsigned char rare_ch = (j < (int)rare_chars.size() ? rare_chars[j] : 0);
              rare_trunk.push_back((c_rank[rare_ch] - sb_rank[sb_ptr - 256 + rare_ch]) << 8);
            }
          rare_trunk[rare_trunk_ptr++] |= rare_map[c]; // Mapping of rare c.
          ++nofreq_cnt;
        }
        ++c_rank[c];
      }
      for (int j = 0; j < rare_cnt; ++j) {
        unsigned char rare_ch = (j < (int)rare_chars.size() ? rare_chars[j] : 0);
        rare_trunk.push_back((c_rank[rare_ch] - sb_rank[sb_ptr - 256 + rare_ch]) << 8);
      }
    }
    c_rank[0] -= n_block * k_block_size - length;

    // debug ////////
    fprintf(stderr, "  rare_trunk.size() = %lu (%.2Lf millions)\n",
        rare_trunk.size(), (long double)rare_trunk.size() / 1000000.L);
    // freq_hit = rare_hit = no_hit = avg_freq_scanning = avg_rare_scanning = 0;
    fprintf(stderr, "  freq_cnt_total = %ld, rare_cnt_total = %ld\n",
        freq_cnt_total, rare_cnt_total);
    ////////////////
  }

  inline long rank(long i, unsigned char c) {
    if (i <= 0) return 0L;
    else if (i >= m_length) return c_rank[c];

    long block_id = (i >> k_block_size_bits), sb_id = (i >> k_sb_size_bits);
    long sb_count = sb_rank[(sb_id << 8) + c];
    unsigned char type = m_mapping[2 * (c * n_block + block_id)];
    unsigned char c_map = m_mapping[2 * (c * n_block + block_id) + 1];

    long freq_cnt_bits = (block_header[block_id] & 255L);
    long rare_cnt_bits = ((block_header[block_id] >> 8) & 255L);
    long micro_block_id = (i >> freq_cnt_bits);

    if (type == FREQUENT_CHAR) {
      // debug ////
      // ++freq_hit;
      //////////////

      long b_count = freq_trunk[(micro_block_id << freq_cnt_bits) + c_map] >> 8;

      // Full/empty context-block optimization -- not effective.
      // unsigned freq_cnt = (1 << freq_cnt_bits);
      // long nextb_count = freq_trunk[((micro_block_id + 1) << freq_cnt_bits) + c_map] >> 8;
      // if ((micro_block_id + 1) & ((1U << (k_block_size_bits - freq_cnt_bits)) - 1) &&
      //     (!((nextb_count - b_count) & (freq_cnt - 1)))) {
      //   if (b_count + freq_cnt == nextb_count) return sb_count + b_count + (i & (freq_cnt - 1));
      //   else if (b_count == nextb_count) return sb_count + b_count;
      // }

      long extra = 0;
      for (long j = (micro_block_id << freq_cnt_bits); j < i; ++j) {
        if ((freq_trunk[j] & 255) == c_map) ++extra;
        // debug ///////
        // ++avg_freq_scanning;
        ////////////////
      }

      return sb_count + b_count + extra;
    } else if (type == RARE_CHAR) {
      // debug ////////
      // ++rare_hit;
      /////////////////

      // Compute new_i.
      long rare_trunk_ptr = (block_header[block_id] >> 16);
      long new_i = freq_trunk[((micro_block_id + 1) << freq_cnt_bits) - 1] >> 8;

      // Full/empty context-block optimization -- not effective.
      // unsigned freq_cnt = (1 << freq_cnt_bits);
      // unsigned rare_cnt = (1 << rare_cnt_bits);      
      // long nextb_count = freq_trunk[((micro_block_id + 2) << freq_cnt_bits) - 1] >> 8;
      // if (   ((micro_block_id + 1) & ((1U << (k_block_size_bits - freq_cnt_bits)) - 1)) != 0
      //     && ((nextb_count - new_i) & (freq_cnt - 1)) == 0) { // full or empty context block
      //   if (new_i != nextb_count) new_i += (i & (freq_cnt - 1)); // full block optimization
      // } else {
      for (long j = (micro_block_id << freq_cnt_bits); j < i; ++j)
        if ((freq_trunk[j] & 255) + 1 == (1U << freq_cnt_bits)) ++new_i;
      // }
      
      // Answer a query on rare trunk.
      long rare_micro_block_id = (new_i >> rare_cnt_bits);
      long b_count = rare_trunk[rare_trunk_ptr + (rare_micro_block_id << rare_cnt_bits) + c_map] >> 8;
      long extra = 0;

      // Full/empty rare context-block optimization -- not effective.
      // long rare_nextb_count = rare_trunk[rare_trunk_ptr + ((rare_micro_block_id + 1) << rare_cnt_bits) + c_map] >> 8;
      // if (b_count + rare_cnt == rare_nextb_count) return sb_count + b_count + (new_i & (rare_cnt - 1));
      // else if (b_count == rare_nextb_count) return sb_count + b_count;

      for (long j = (rare_micro_block_id << rare_cnt_bits); j < new_i; ++j) {
        if ((rare_trunk[rare_trunk_ptr + j] & 255) == c_map) ++extra;
        // debug ///////
        // ++avg_rare_scanning;
        ////////////////
      }
        
      return sb_count + b_count + extra;
    } else {
      // debug ///////
      // ++no_hit;
      ////////////////
      
      while (block_id < n_block && (block_id & k_blocks_in_sb_mask) && m_mapping[2 * (c * n_block + block_id)] == NOT_OCCURRING)
        ++block_id;
      if (block_id == n_block) return c_rank[c];
      else if (!(block_id & k_blocks_in_sb_mask)) return sb_rank[256 * (block_id >> k_blocks_in_sb_bits) + c];
      else return rank(block_id << k_block_size_bits, c);
    }
  }

  ~context_rank_4n() {
    delete[] sb_rank;
    delete[] freq_trunk;
    delete[] block_header;
    delete[] m_mapping;
    delete[] c_rank;

    // debug //////
    // fprintf(stderr, "  freq_hit=%ld, rare_hit=%ld, no_hit=%ld\n", freq_hit, rare_hit, no_hit);
    // fprintf(stderr, "  avg_freq_scanning_length = %.3Lf\n", (long double)avg_freq_scanning / freq_hit);
    // fprintf(stderr, "  avg_rare_scanning_length = %.3Lf\n", (long double)avg_rare_scanning / rare_hit);
    ///////////////
  }

private:
  static const int k_sb_size = (1 << k_sb_size_bits);
  static const int k_block_size = (1 << k_block_size_bits);
  static const int k_blocks_in_sb_bits = k_sb_size_bits - k_block_size_bits;
  static const int k_blocks_in_sb = (1 << k_blocks_in_sb_bits);
  static const int k_blocks_in_sb_mask = k_blocks_in_sb - 1;

  const long m_length, n_block, n_sblock;
  long *sb_rank, *c_rank, *block_header;
  unsigned char *m_mapping;
  
  unsigned *freq_trunk;
  std::vector<unsigned> rare_trunk;

  // debug /////
  // long rare_hit, freq_hit, no_hit, avg_freq_scanning, avg_rare_scanning;
  long freq_cnt_total, rare_cnt_total;
  //////////////
};

/*struct fast_sigma_rank_4n {
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
};*/

#endif // __RANK_H

