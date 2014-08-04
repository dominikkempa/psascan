// Data structure answering general rank queries over byte alphabet based
// on the rank data structure from the 'bwtdisk' implementation of
// algorithm (http://people.unipmn.it/manzini/bwtdisk/) from
//
//     Paolo Ferragina, Travis Gagie, Giovanni Manzini:
//     Lightweight Data Indexing and Compression in External Memory.
//     Algorithmica 63(3): 707-730 (2012)
//
// Our key modification is applying two new techniques:
//
//   * alphabet partitioning
//
//     Jérémy Barbay, Travis Gagie, Gonzalo Navarro, Yakov Nekrich:
//     Alphabet Partitioning for Compressed Rank/Select and Applications.
//     Proc. ISAAC 2010.
//
//  * fixed block boosting
//
//    Juha Kärkkäinen, Simon J. Puglisi:
//    Fixed Block Compression Boosting in FM-Indexes.
//    Proc. SPIRE 2011.

#ifndef __RANK_H_INCLUDED
#define __RANK_H_INCLUDED

#include <algorithm>
#include <vector>

#include "utils.h"

#define FREQUENT_CHAR 0
#define RARE_CHAR 1
#define NOT_OCCURRING 2

struct context_rank_4n {
  static const int k_sb_size_bits = 13;    // log_2(super block size), 1MiB
  static const int k_block_size_bits = 9;  // log 2(block size), 64KiB
                                           // anything 16-20 is a pretty good value.

  context_rank_4n(unsigned char *text, long length)
      : m_length(length),
        n_block((length + k_block_size - 1) / k_block_size),
        n_sblock((n_block + k_blocks_in_sb - 1) / k_blocks_in_sb) {

    c_rank       = new long[256];
    std::fill(c_rank, c_rank + 256, 0L);
    block_header = new long[n_block];
    sb_rank      = new long[256 * n_sblock];

    m_mapping    = new unsigned char[2 * 256 * n_block];
    freq_trunk   = new unsigned[n_block * k_block_size];

    if (!length) return;

    // debug ///////
    // fprintf(stderr, "  m_mapping takes %.2Lf MiB\n", (long double)(512 * n_block) / (1 << 20));
    // fprintf(stderr, "  k_sb_size_bits = %d, k_block_size_bits = %d\n", k_sb_size_bits, k_block_size_bits);
    // freq_cnt_total = rare_cnt_total = 0;
    ////////////////

    // Compute rare trunk size.
    long rare_trunk_size = 0;
    for (long start = 0; start < n_block * k_block_size; start += k_block_size) {
      long end = start + k_block_size; // text[start..end) is current block

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
        rare_cnt_total += sorted_chars[i].first;
      }
      for (int i = rare_cnt; i < (int)sorted_chars.size(); ++i) {
        freq_chars.push_back(sorted_chars[i].second);
        freq_cnt_total += sorted_chars[i].first;
      }
      long rare_cnt_bits = 0;
      if (rare_cnt) {                              // If there are rare symbols,
        rare_cnt_bits = utils::log2ceil(rare_cnt); // round up rare_cnt to the
        rare_cnt = (1 << rare_cnt_bits);           // smallest power of two.
      }

      // 3. Compute and store symbols mapping.
      bool is_freq[256] = {false};
      for (unsigned i = 0; i < freq_chars.size(); ++i) {
        unsigned char c = freq_chars[i];
        is_freq[c] = true;
      }

      for (long i = start, nofreq_cnt = 0; i < end; ++i) {
        unsigned char c = (i < length ? text[i] : 0);
        if (!is_freq[c]) {
          if (nofreq_cnt % rare_cnt == 0)
            rare_trunk_size += rare_cnt;
          ++nofreq_cnt;
        }
      }
      rare_trunk_size += rare_cnt;
    }
    
    rare_trunk.reserve(rare_trunk_size);

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
        rare_cnt_total += sorted_chars[i].first;
      }
      for (int i = rare_cnt; i < (int)sorted_chars.size(); ++i) {
        freq_chars.push_back(sorted_chars[i].second);
        freq_cnt_total += sorted_chars[i].first;
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

      // Not necessary -- easier to debug.
      std::sort(freq_chars.begin(), freq_chars.end());
      std::sort(rare_chars.begin(), rare_chars.end());

      // 4. Compute and store symbols mapping.
      bool is_freq[256] = {false};
      unsigned char freq_map[256], rare_map[256];
      for (unsigned c = 0; c < 256; ++c) {
        m_mapping[2 * (c * n_block + block_id)] = NOT_OCCURRING;
        m_mapping[2 * (c * n_block + block_id) + 1] = 0;
      }
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

    if (rare_trunk_size != (long)rare_trunk.size() ||
        rare_trunk_size != (long)rare_trunk.capacity()) {
      fprintf(stderr, "\n\nError during rank initialization:\n");
      fprintf(stderr, "\trare_trunk_size = %ld\n", rare_trunk_size);
      fprintf(stderr, "\trare_trunk.size() = %ld\n", (long)rare_trunk.size());
      fprintf(stderr, "\trare_trunk.capacity() = %ld\n", (long)rare_trunk.capacity());
      std::exit(EXIT_FAILURE);
    }

    // debug //
    // fprintf(stderr, "  rare_trunk.size() = %lu (%.2Lf millions)\n",
    //     rare_trunk.size(), (long double)rare_trunk.size() / 1000000.L);
    // fprintf(stderr, "  freq_cnt_total = %ld, rare_cnt_total = %ld\n",
    //     freq_cnt_total, rare_cnt_total);
    ///////////
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
      long b_count = freq_trunk[(micro_block_id << freq_cnt_bits) + c_map] >> 8;
      long extra = 0;
      for (long j = (micro_block_id << freq_cnt_bits); j < i; ++j)
        if ((freq_trunk[j] & 255) == c_map) ++extra;

      return sb_count + b_count + extra;
    } else if (type == RARE_CHAR) {
      // Compute new_i.
      long rare_trunk_ptr = (block_header[block_id] >> 16);
      long new_i = freq_trunk[((micro_block_id + 1) << freq_cnt_bits) - 1] >> 8;

      for (long j = (micro_block_id << freq_cnt_bits); j < i; ++j)
        if ((freq_trunk[j] & 255) + 1 == (1U << freq_cnt_bits)) ++new_i;
      
      // Answer a query on rare trunk.
      long rare_micro_block_id = (new_i >> rare_cnt_bits);
      long b_count = rare_trunk[rare_trunk_ptr + (rare_micro_block_id << rare_cnt_bits) + c_map] >> 8;
      long extra = 0;

      for (long j = (rare_micro_block_id << rare_cnt_bits); j < new_i; ++j)
        if ((rare_trunk[rare_trunk_ptr + j] & 255) == c_map) ++extra;
        
      return sb_count + b_count + extra;
    } else {
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
  }

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

  long freq_cnt_total, rare_cnt_total;
};

#endif // __RANK_H_INCLUDED
