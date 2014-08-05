//==============================================================================
// Data structure answering general rank queries over byte alphabet based
// on the rank data structure from the bwtdisk implementation
// (http://people.unipmn.it/manzini/bwtdisk/) of the algorithm from:
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
//   * fixed block boosting
//
//     Juha Kärkkäinen, Simon J. Puglisi:
//     Fixed Block Compression Boosting in FM-Indexes.
//     Proc. SPIRE 2011.
//==============================================================================

#ifndef __RANK4N_H_INCLUDED
#define __RANK4N_H_INCLUDED

#include <algorithm>
#include <vector>

#include "utils.h"

#define FREQUENT_CHAR 0
#define RARE_CHAR 1
#define NOT_OCCURRING 2


template<long k_sb_size_bits = 24, long k_block_size_bits = 18>
struct rank4n {
  rank4n(unsigned char *text, long length)
      : m_length(length),
        n_blocks((length + k_block_size - 1) >> k_block_size_bits),
        n_sblocks((n_blocks + k_blocks_in_sb - 1) / k_blocks_in_sb) {

    c_rank = new long[256];
    std::fill(c_rank, c_rank + 256, 0L);
    if (!length) return;

    //--------------------------------------------------------------------------
    // STEP 1: compute rare trunk size, block_header and m_mapping.
    //--------------------------------------------------------------------------
    block_header = new long[n_blocks];
    m_mapping = new unsigned char[2L * 256L * n_blocks];
    long rare_trunk_size = 0;

    bool isfreq[256];
    std::vector<std::pair<unsigned, unsigned> > sorted_chars;
    std::vector<unsigned char> freq_chars, rare_chars;
    unsigned block_count[256];
    for (long block_id = 0; block_id < n_blocks; ++block_id) {
      long block_beg = block_id << k_block_size_bits;
      long block_end = block_beg + k_block_size;

      // Compute symbol counts inside a block.
      std::fill(block_count, block_count + 256, 0L);
      for (long j = block_beg; j < block_end; ++j) {
        unsigned char c = (j < m_length ? text[j] : 0);
        ++block_count[c];
      }

      // Sort symbol counts by frequencies.
      sorted_chars.clear();
      for (long j = 0; j < 256; ++j) if (block_count[j])
        sorted_chars.push_back(std::make_pair(block_count[j], j));
      std::sort(sorted_chars.begin(), sorted_chars.end());

      // Separate (at most, due to rounding of freq_cnt)
      // about 3% of rarest symbols.
      long rare_cnt = 0L, rare_sum = 0L;
      while (rare_cnt < (long)sorted_chars.size() &&
          16L * (rare_sum + sorted_chars[rare_cnt].first) <= k_block_size)
        rare_sum += sorted_chars[rare_cnt++].first;

      // Compute freq_cnt. Then round up freq_cnt + 1 (+1 is
      // for rare char marker) to the smallest power of two.
      // Note: this does not work when freq_cnt == 256.
      long freq_cnt = (long)sorted_chars.size() - rare_cnt;
      long freq_cnt_bits = utils::log2ceil(freq_cnt + 1);
      freq_cnt = (1L << freq_cnt_bits);

      // Recomputed rare_cnt (note the +1).
      rare_cnt = std::max(0L, (long)sorted_chars.size() - freq_cnt + 1);

      // Compute freq and rare chars.
      rare_chars.clear();
      freq_chars.clear();
      for (int i = 0; i < rare_cnt; ++i)
        rare_chars.push_back(sorted_chars[i].second);
      for (int i = rare_cnt; i < (long)sorted_chars.size(); ++i)
        freq_chars.push_back(sorted_chars[i].second);
 
      // If there are rare symbols, round up
      // rare_cnt to the smallest power of two.
      long rare_cnt_bits = 0L;
      if (rare_cnt) {
        rare_cnt_bits = utils::log2ceil(rare_cnt);
        rare_cnt = (1L << rare_cnt_bits);
      }

      // Store log2 of rare_cnt and freq_cnt into block header.
      block_header[block_id] = (rare_cnt_bits << 8) + freq_cnt_bits;

      // Compute and store symbols mapping.
      std::sort(freq_chars.begin(), freq_chars.end());
      std::sort(rare_chars.begin(), rare_chars.end());
      std::fill(isfreq, isfreq + 256, false);
      for (unsigned c = 0; c < 256; ++c)
        m_mapping[2 * (c * n_blocks + block_id)] = NOT_OCCURRING;
      for (unsigned i = 0; i < freq_chars.size(); ++i) {
        unsigned char c = freq_chars[i];
        isfreq[c] = true;
        m_mapping[2 * (c * n_blocks + block_id) + 1] = i;
        m_mapping[2 * (c * n_blocks + block_id)] = FREQUENT_CHAR;
      }
      for (unsigned i = 0; i < rare_chars.size(); ++i) {
        unsigned char c = rare_chars[i];
        m_mapping[2 * (c * n_blocks + block_id) + 1] = i;
        m_mapping[2 * (c * n_blocks + block_id)] = RARE_CHAR;
      }

      // Update rare_trunk_size.
      if (rare_cnt) {
        long nofreq_cnt = 0L;
        for (long i = block_beg; i < block_end; ++i) {
          unsigned char c = (i < m_length ? text[i] : 0);
          if (!isfreq[c]) ++nofreq_cnt;
        }
        long rare_micro_blocks = 1 + (nofreq_cnt + rare_cnt - 1) / rare_cnt;
        rare_trunk_size += rare_micro_blocks * rare_cnt;
      }
    }
    
    rare_trunk.resize(rare_trunk_size);
    freq_trunk = (unsigned *)malloc(n_blocks * k_block_size * sizeof(unsigned));
    sb_rank = new long[256 * n_sblocks];

    //--------------------------------------------------------------------------
    // STEP 2: compute sb_rank, freq_trunk and rare_trunk.
    //--------------------------------------------------------------------------
    unsigned char freq_map[256], rare_map[256];
    long r_filled = 0L; // for rare trunk
    for (long block_id = 0; block_id < n_blocks; ++block_id) {
      long r_ptr = r_filled; // where to write next
      block_header[block_id] |= (r_filled << 16);

      long block_beg = block_id << k_block_size_bits;
      long block_end = block_beg + k_block_size;

      //------------------------------------------------------------------------
      // Process block block_id.
      //------------------------------------------------------------------------

      long freq_cnt_bits = (block_header[block_id] & 255L);
      long rare_cnt_bits = ((block_header[block_id] >> 8) & 255L);
      long freq_cnt = (1L << freq_cnt_bits);
      long freq_cnt_mask = freq_cnt - 1;
      long rare_cnt = (1L << rare_cnt_bits);
      long rare_cnt_mask = rare_cnt - 1;

      freq_chars.clear();
      rare_chars.clear();
      std::fill(isfreq, isfreq + 256, false);
      for (long j = 0; j < 256; ++j) {
        unsigned char type = m_mapping[2 * (j * n_blocks + block_id)];
        if (type == FREQUENT_CHAR) {
          isfreq[j] = true;
          freq_chars.push_back(j);
          freq_map[j] = m_mapping[2 * (j * n_blocks + block_id) + 1];
        } else if (type == RARE_CHAR) {
          rare_chars.push_back(j);
          rare_map[j] = m_mapping[2 * (j * n_blocks + block_id) + 1];
          freq_map[j] = freq_cnt - 1;
        }
      }

      if (rare_chars.empty()) {
        rare_cnt_bits = 0;
        rare_cnt = 0;
      }

      // Compute ranks at superblock boundary.
      long sb_id = (block_beg >> k_sb_size_bits);
      if (!(block_beg & k_sb_size_mask))
        for (long i = 0; i < 256; ++i)
          sb_rank[(sb_id << 8) + i] = c_rank[i];

      // Compute freq and rare trunk of current block.
      long nofreq_cnt = 0L;
      for (long i = block_beg; i < block_end; ++i) {
        unsigned char c = (i < m_length ? text[i] : 0);

        //----------------------------------------------------------------------
        // Invariant: for any symbol a, c_rank[a] = number of occurrence of
        // symbols a in prefix text[0..i).
        //----------------------------------------------------------------------
        // Invariant: for any symbol a, sb_rank[(sb_id << 8) + a] = the
        // number of occurrence of a in a prefix of text up to the closest
        // superblock boundary. All ranks we store in the trunk are relative
        // to this boundary thus we compute them as
        // c_rank[a] - sb_rank[(sb_id << 8) + a].
        //----------------------------------------------------------------------

        // Compute ranks at the freq microblock boundary.
        if (!(i & freq_cnt_mask)) {
          freq_trunk[i + freq_cnt - 1] = (nofreq_cnt << 8);
          for (long j = 0; j + 1 < freq_cnt; ++j) {
            unsigned char ch = (j < (long)freq_chars.size() ? freq_chars[j] : 0);
            long local_rank = c_rank[ch] - sb_rank[(sb_id << 8) + ch];
            freq_trunk[i + j] = (local_rank << 8);
          }
        }

        // Store freq symbol mapping.
        freq_trunk[i] |= freq_map[c];

        // Handle rare symbol.
        if (!isfreq[c]) {
          // Compute ranks at the rare microblock boundary.
          if (!(nofreq_cnt & rare_cnt_mask)) {
            for (long j = 0; j < rare_cnt; ++j) {
              unsigned char ch = (j < (long)rare_chars.size() ? rare_chars[j] : 0);
              long local_rank = c_rank[ch] - sb_rank[(sb_id << 8) + ch];
              rare_trunk[r_filled++] = (local_rank << 8);
            }
          }

          // Store rare symbol mapping.
          rare_trunk[r_ptr++] |= rare_map[c];
          ++nofreq_cnt;
        }

        ++c_rank[c];
      }

      for (long j = 0; j < rare_cnt; ++j) {
        unsigned char ch = (j < (long)rare_chars.size() ? rare_chars[j] : 0);
        long local_rank = c_rank[ch] - sb_rank[(sb_id << 8) + ch];
        rare_trunk[r_filled++] = (local_rank << 8);
      }
    }

    c_rank[0] -= n_blocks * k_block_size - length;
  }

  inline long rank(long i, unsigned char c) {
    if (i <= 0) return 0L;
    else if (i >= m_length) return c_rank[c];

    long block_id = (i >> k_block_size_bits);
    long sb_id = (i >> k_sb_size_bits);
    long sb_count = sb_rank[(sb_id << 8) + c];

    unsigned char type = m_mapping[2 * (c * n_blocks + block_id)];
    unsigned char c_map = m_mapping[2 * (c * n_blocks + block_id) + 1];

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
      long b_count = rare_trunk[rare_trunk_ptr +
        (rare_micro_block_id << rare_cnt_bits) + c_map] >> 8;
      long extra = 0;

      for (long j = (rare_micro_block_id << rare_cnt_bits); j < new_i; ++j)
        if ((rare_trunk[rare_trunk_ptr + j] & 255) == c_map) ++extra;
        
      return sb_count + b_count + extra;
    } else {
      while (block_id < n_blocks && (block_id & k_blocks_in_sb_mask) &&
          m_mapping[2 * (c * n_blocks + block_id)] == NOT_OCCURRING)
        ++block_id;

      if (block_id == n_blocks) return c_rank[c];
      else if (!(block_id & k_blocks_in_sb_mask))
        return sb_rank[256 * (block_id >> k_blocks_in_sb_bits) + c];
      else return rank(block_id << k_block_size_bits, c);
    }
  }

  ~rank4n() {
    if (m_length) {
      delete[] sb_rank;
      delete[] block_header;
      delete[] m_mapping;
      free(freq_trunk);
    }
    delete[] c_rank;
  }

  static const int k_sb_size = (1 << k_sb_size_bits);
  static const int k_sb_size_mask = k_sb_size - 1;
  static const int k_block_size = (1 << k_block_size_bits);
  static const int k_blocks_in_sb_bits = k_sb_size_bits - k_block_size_bits;
  static const int k_blocks_in_sb = (1 << k_blocks_in_sb_bits);
  static const int k_blocks_in_sb_mask = k_blocks_in_sb - 1;

  const long m_length, n_blocks, n_sblocks;
  long *sb_rank, *c_rank, *block_header;
  unsigned char *m_mapping;
  
  unsigned *freq_trunk;
  std::vector<unsigned> rare_trunk;
};

#endif // __RANK4N_H_INCLUDED
