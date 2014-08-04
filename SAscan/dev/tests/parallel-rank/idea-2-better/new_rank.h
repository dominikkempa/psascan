//==============================================================================
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
//==============================================================================

#ifndef __NEW_RANK_H_INCLUDED
#define __NEW_RANK_H_INCLUDED

#include <algorithm>
#include <vector>
#include <thread>

#include "utils.h"

#define FREQUENT_CHAR 0
#define RARE_CHAR 1
#define NOT_OCCURRING 2

template<long k_sb_size_bits = 24, long k_block_size_bits = 18>
struct rank4n {
  static void construction_step_1(rank4n &r, unsigned char *text, long range_beg,
      long range_end, long *range_count, long &rare_trunk_size) {
    // Initialize output values.
    std::fill(range_count, range_count + 256, 0L);
    rare_trunk_size = 0L;

    std::vector<std::pair<unsigned, unsigned> > sorted_chars;
    for (long block_id = range_beg; block_id < range_end; ++block_id) {
      long block_beg = block_id << k_block_size_bits;
      long block_end = block_beg + k_block_size;

      //--------------------------------------------------------------------------
      // Process block i.
      //--------------------------------------------------------------------------

      // Compute symbols counts inside a block and update
      // symbols counts inside the range of blocks.
      unsigned block_count[256] = {0U};
      for (long j = block_beg; j < block_end; ++j) {
        unsigned char c = (j < r.m_length ? text[j] : 0);
        ++block_count[c];
        ++range_count[c];
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
      std::vector<unsigned char> freq_chars, rare_chars;
      for (long i = 0; i < rare_cnt; ++i)
        rare_chars.push_back(sorted_chars[i].second);
      for (long i = rare_cnt; i < (long)sorted_chars.size(); ++i)
        freq_chars.push_back(sorted_chars[i].second);

      // If there are rare symbols, round up rare_cnt
      // to the smallest power of two.
      long rare_cnt_bits = 0L;
      if (rare_cnt) {
        rare_cnt_bits = utils::log2ceil(rare_cnt);
        rare_cnt = (1L << rare_cnt_bits);
      }

      // Store log2 of rare_cnt and freq_cnt into block header.
      r.block_header[block_id] = (rare_cnt_bits << 8) +  // one byte
                                 freq_cnt_bits;          // one byte

      std::sort(freq_chars.begin(), freq_chars.end());
      std::sort(rare_chars.begin(), rare_chars.end());

      // Compute and store symbols mapping.
      long n_blocks = (r.m_length + k_block_size - 1) / k_block_size;
      bool is_freq[256] = {false};
      for (unsigned c = 0; c < 256; ++c) {
        r.m_mapping[2 * (c * n_blocks + block_id)] = NOT_OCCURRING;
        r.m_mapping[2 * (c * n_blocks + block_id) + 1] = 0; // anything.
      }
      for (unsigned i = 0; i < freq_chars.size(); ++i) {
        unsigned char c = freq_chars[i];
        is_freq[c] = true;
        r.m_mapping[2 * (c * n_blocks + block_id) + 1] = i;
        r.m_mapping[2 * (c * n_blocks + block_id)] = FREQUENT_CHAR;
      }
      for (unsigned i = 0; i < rare_chars.size(); ++i) {
        unsigned char c = rare_chars[i];
        r.m_mapping[2 * (c * n_blocks + block_id) + 1] = i;
        r.m_mapping[2 * (c * n_blocks + block_id)] = RARE_CHAR;
      }

      // Update rare_trunk_size.
      if (rare_cnt) {
        long nofreq_cnt = 0L;
        for (long i = block_beg; i < block_end; ++i) {
          unsigned char c = (i < r.m_length ? text[i] : 0);
          if (!is_freq[c]) ++nofreq_cnt;
        }
        long rare_micro_blocks = 1 + (nofreq_cnt + rare_cnt - 1) / rare_cnt;
        rare_trunk_size += rare_micro_blocks * rare_cnt;
      }
    }
  }


  static void construction_step_3(rank4n &r, unsigned char *text, long range_beg,
      long range_end, long *count) {
    long *cur_count = new long[256];
    std::copy(count, count + 256, cur_count);
    for (long block_id = range_beg; block_id < range_end; ++block_id) {
      long block_beg = block_id << k_block_size_bits;
      long block_end = block_beg + k_block_size;

      if (!(block_beg & k_sb_size_mask)) {
        // Current block starts at the superblock boundary --
        // copy current symbols counts as rank values at
        // superblock boundary.
        long sb_id = (block_beg >> k_sb_size_bits);
        std::copy(cur_count, cur_count + 256, r.sb_rank + (sb_id << 8));
      }

      // Update symbol counts.
      for (long j = block_beg; j < block_end; ++j) {
        unsigned char c = (j < r.m_length ? text[j] : 0);
        ++cur_count[c];
      }
    }

    delete[] cur_count;
  }


  static void construction_step_4(rank4n &r, unsigned char *text, long range_beg,
      long range_end, long *count, long rare_trunk_ptr) {
    bool isfreq[256];
    unsigned char freq_map[256], rare_map[256];
    std::vector<unsigned char> freq_chars, rare_chars;

    for (long block_id = range_beg; block_id < range_end; ++block_id) {
      long rare_trunk_writing_ptr = rare_trunk_ptr;
      r.block_header[block_id] |= (rare_trunk_ptr << 16);

      long block_beg = block_id << k_block_size_bits;
      long block_end = block_beg + k_block_size;

      //------------------------------------------------------------------------
      // Process block i.
      //------------------------------------------------------------------------

      long freq_cnt_bits = (r.block_header[block_id] & 255L);
      long rare_cnt_bits = ((r.block_header[block_id] >> 8) & 255L);
      long freq_cnt = (1L << freq_cnt_bits);
      long freq_cnt_mask = freq_cnt - 1;
      long rare_cnt = (1L << rare_cnt_bits);
      long rare_cnt_mask = rare_cnt - 1;

      freq_chars.clear();
      rare_chars.clear();
      std::fill(isfreq, isfreq + 256, false);
      for (long j = 0; j < 256; ++j) {
        unsigned char type = r.m_mapping[2 * (j * r.n_blocks + block_id)];
        if (type == FREQUENT_CHAR) {
          isfreq[j] = true;
          freq_chars.push_back(j);
          freq_map[j] = r.m_mapping[2 * (j * r.n_blocks + block_id) + 1];
        } else if (type == RARE_CHAR) {
          rare_chars.push_back(j);
          rare_map[j] = r.m_mapping[2 * (j * r.n_blocks + block_id) + 1];
          freq_map[j] = freq_cnt - 1;
        }
      }

      if (rare_chars.empty()) {
        rare_cnt_bits = 0;
        rare_cnt = 0;
      }

      // Compute the freq and rare trunk of block i.
      long nofreq_cnt = 0L;
      long sb_id = block_beg >> k_sb_size_bits;
      for (long i = block_beg; i < block_end; ++i) {
        unsigned char c = (i < r.m_length ? text[i] : 0);

        //----------------------------------------------------------------------
        // Invariant: for any symbol a, count[a] = number of occurrence of
        // symbols a in prefix text[0..i).
        //----------------------------------------------------------------------
        // Invariant: for any symbol a, r.sb_rank[(sb_id << 8) + a] = the
        // number of occurrence of a in a prefix of text up to the closest
        // superblock boundary. All ranks we store in the trunk are relative
        // to this boundary thus we compute them as
        // count[a] - r.sb_rank[(sb_id << 8) + a].
        //----------------------------------------------------------------------

        if (!(i & freq_cnt_mask)) {
          for (long j = 0; j + 1 < freq_cnt; ++j) {
            unsigned char freq_ch = (j < (long)freq_chars.size() ? freq_chars[j] : 0);
            r.freq_trunk[i + j] = (count[freq_ch] - r.sb_rank[(sb_id << 8) + freq_ch]) << 8;
          }
          r.freq_trunk[i + freq_cnt - 1] = (nofreq_cnt << 8);
        }
        r.freq_trunk[i] |= freq_map[c]; // mapping of frequent c
        if (!isfreq[c]) {
          if (!(nofreq_cnt & rare_cnt_mask)) {
            for (long j = 0; j < rare_cnt; ++j) {
              unsigned char rare_ch = (j < (long)rare_chars.size() ? rare_chars[j] : 0);
              r.rare_trunk[rare_trunk_ptr++] = (count[rare_ch] - r.sb_rank[(sb_id << 8) + rare_ch]) << 8;
            }
          }
          r.rare_trunk[rare_trunk_writing_ptr++] |= rare_map[c]; // mapping of rare c
          ++nofreq_cnt;
        }

        ++count[c];
      }

      for (long j = 0; j < rare_cnt; ++j) {
        unsigned char rare_ch = (j < (long)rare_chars.size() ? rare_chars[j] : 0);
        r.rare_trunk[rare_trunk_ptr++] =
          (count[rare_ch] - r.sb_rank[(sb_id << 8) + rare_ch]) << 8;
      }
    }
  }

  rank4n(unsigned char *text, long length, long max_threads)
    : m_length(length),
      n_blocks((length + k_block_size - 1) / k_block_size),
      n_sblock((n_blocks + k_blocks_in_sb - 1) / k_blocks_in_sb) {

    // Makes the structure work if length == 0.
    c_rank = new long[256];
    std::fill(c_rank, c_rank + 256, 0L);
    if (!length) return;

    //--------------------------------------------------------------------------
    // STEP 1: split all blocks into (rougly equal size) ranges of blocks.
    //
    // Each thread will handle one range of blocks. In first round each thread
    // returns:
    //   * total symbol count inside all blocks in the handled range
    //   * the size of rare trunk for that range of blocks
    //
    // During its execution, each thread computes headers for all blocks
    // the handled range. This includes:
    //   * symbols types (rare / freq / non-occurring)
    //   * mapping of any symbol
    //   * values of freq_cnt_bits and rare_cnt_bits for each block
    // All this is stored inside the block_header and m_mapping arrays.
    //--------------------------------------------------------------------------
    long range_size = (n_blocks + max_threads - 1) / max_threads; // measured in blocks
    long n_ranges = (n_blocks + range_size - 1) / range_size;

    long **count = new long*[n_ranges];
    for (long i = 0; i < n_ranges; ++i)
      count[i] = new long[256];
    long *rare_trunk_size = new long[n_ranges];
    std::thread **threads = new std::thread*[n_ranges];

    // These will be filled by threads created below.
    block_header = new long[n_blocks];
    m_mapping    = new unsigned char[2 * 256 * n_blocks];

    for (long i = 0; i < n_ranges; ++i) {
      long range_beg = i * range_size;
      long range_end = std::min(range_beg + range_size, n_blocks);
 
      threads[i] = new std::thread(
          construction_step_1, std::ref(*this), text, range_beg,
          range_end, count[i], std::ref(rare_trunk_size[i]));
    }

    for (long i = 0; i < n_ranges; ++i) threads[i]->join();
    for (long i = 0; i < n_ranges; ++i) delete threads[i];
    delete[] threads;


    //--------------------------------------------------------------------------
    // STEP 2: compute non-inclusive partial sum over computed symbols
    //         counts. Same for rare trunk sizes.
    //
    // In addition compute:
    //   * cumulative counts of all symbols (necessary during queries).
    //   * the total rare trunk size (necessary to allocate enough memory).
    //--------------------------------------------------------------------------
    long temp_count[256];
    for (long i = 0; i < n_ranges; ++i) {
      std::copy(count[i], count[i] + 256, temp_count);
      std::copy(c_rank, c_rank + 256, count[i]);
      for (long j = 0; j < 256; ++j) c_rank[j] += temp_count[j];
    }
    c_rank[0] -= n_blocks * k_block_size - length;
    long total_rare_trunk_size = 0L;
    for (long i = 0, t; i < n_ranges; ++i) {
      t = rare_trunk_size[i];
      rare_trunk_size[i] = total_rare_trunk_size;
      total_rare_trunk_size += t;
    }

    //--------------------------------------------------------------------------
    // STEP 3: compute rank at superblock boundaries.
    //--------------------------------------------------------------------------
    sb_rank = new long[256 * n_sblock];
    threads = new std::thread*[n_ranges];
    for (long i = 0; i < n_ranges; ++i) {
      long range_beg = i * range_size;
      long range_end = std::min(range_beg + range_size, n_blocks);

      threads[i] = new std::thread(
          construction_step_3, std::ref(*this), text,
          range_beg, range_end, count[i]);
    }
    for (long i = 0; i < n_ranges; ++i) threads[i]->join();
    for (long i = 0; i < n_ranges; ++i) delete threads[i];
    delete[] threads;

    //--------------------------------------------------------------------------
    // STEP 4: split the blocks into range as in step 1, but this time each
    //         thread can compute the trunk, since it has symbols counts up to
    //         the beginning of the sequence and it can also fill in the rare
    //         trunk.
    //
    // In this step we compute:
    //   * frequent and rare trunk
    //   * symbol counts at superblock boundaries
    //--------------------------------------------------------------------------

    // These will be filled by threads created below.
    // Note: g++ automatically zero-initializes arrays allocated with 'new'
    // in template-parametrized classes. To avoid this, we use malloc.
    freq_trunk = (unsigned *)malloc(n_blocks * k_block_size * sizeof(unsigned));
    rare_trunk.resize(total_rare_trunk_size);

    threads = new std::thread*[n_ranges];
    for (long i = 0; i < n_ranges; ++i) {
      long range_beg = i * range_size;
      long range_end = std::min(range_beg + range_size, n_blocks);

      threads[i] = new std::thread(
          construction_step_4, std::ref(*this), text, range_beg,
          range_end, count[i], rare_trunk_size[i]);
    }
    for (long i = 0; i < n_ranges; ++i) threads[i]->join();
    for (long i = 0; i < n_ranges; ++i) delete threads[i];
    delete[] threads;

    // Clean up.
    for (long i = 0; i < n_ranges; ++i) delete[] count[i];
    delete[] count;
    delete[] rare_trunk_size;
  }


  inline long rank(long i, unsigned char c) {
    if (i <= 0) return 0L;
    else if (i >= m_length) return c_rank[c];

    long block_id = (i >> k_block_size_bits), sb_id = (i >> k_sb_size_bits);
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
      long b_count = rare_trunk[rare_trunk_ptr + (rare_micro_block_id << rare_cnt_bits) + c_map] >> 8;
      long extra = 0;

      for (long j = (rare_micro_block_id << rare_cnt_bits); j < new_i; ++j)
        if ((rare_trunk[rare_trunk_ptr + j] & 255) == c_map) ++extra;
        
      return sb_count + b_count + extra;
    } else {
      while (block_id < n_blocks && (block_id & k_blocks_in_sb_mask) && m_mapping[2 * (c * n_blocks + block_id)] == NOT_OCCURRING)
        ++block_id;
      if (block_id == n_blocks) return c_rank[c];
      else if (!(block_id & k_blocks_in_sb_mask)) return sb_rank[256 * (block_id >> k_blocks_in_sb_bits) + c];
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

  const long m_length, n_blocks, n_sblock;
  long *sb_rank, *c_rank, *block_header;
  unsigned char *m_mapping;
  
  unsigned *freq_trunk;
  std::vector<unsigned> rare_trunk;

  long freq_cnt_total, rare_cnt_total;
};

#endif // __NEW_RANK_H_INCLUDED
