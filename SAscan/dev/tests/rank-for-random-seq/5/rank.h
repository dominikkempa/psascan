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
//     JÃ©rÃ©my Barbay, Travis Gagie, Gonzalo Navarro, Yakov Nekrich:
//     Alphabet Partitioning for Compressed Rank/Select and Applications.
//     Proc. ISAAC 2010.
//
//   * fixed block boosting
//
//     Juha KÃ¤rkkÃ¤inen, Simon J. Puglisi:
//     Fixed Block Compression Boosting in FM-Indexes.
//     Proc. SPIRE 2011.
//==============================================================================

#ifndef __RANK4N_H_INCLUDED
#define __RANK4N_H_INCLUDED

#include <algorithm>
#include <vector>

#include "utils.h"


template<unsigned k_sblock_size_log = 24, unsigned k_cblock_size_log = 18, unsigned k_sigma_log = 8>
class rank4n {
  private:
    static const unsigned k_cblock_size;
    static const unsigned k_cblock_size_mask;
    static const unsigned k_cblock_size_mask_neg;
    static const unsigned k_cblocks_in_sblock_log;
    static const unsigned k_cblocks_in_sblock;
    static const unsigned k_cblocks_in_sblock_mask;
    static const unsigned k_2cblock_size;
    static const unsigned k_2cblock_size_mask;
    static const unsigned k_sblock_size;
    static const unsigned k_sblock_size_mask;
    static const int k_sigma;
    static const int k_sigma_mask;

    static const unsigned k_char_type_freq =    0x01;
    static const unsigned k_char_type_rare =    0x02;
    static const unsigned k_char_type_missing = 0x03;

    long m_length;   // length of original sequence
    long n_cblocks;  // number of context blocks
    long n_sblocks;  // number of super block

    long *m_cblock_header;
    unsigned long *m_cblock_header2;
    long *m_sblock_header;
    unsigned char *m_cblock_type;
    unsigned char *m_cblock_mapping;

    unsigned *m_freq_trunk;
    unsigned *m_rare_trunk;

  public:
    long *m_count; // symbol counts

  public:
    rank4n(unsigned char *text, long length, long) {
      m_length = length;

      // Compute the number of blocks.
      n_cblocks = (m_length + k_cblock_size - 1) / k_cblock_size;
      n_sblocks = (n_cblocks + k_cblocks_in_sblock - 1) / k_cblocks_in_sblock;

      // Allocate and initialize symbol counts.
      m_count = (long *)malloc(256L * sizeof(long));
      std::fill(m_count, m_count + 256, 0L);

      if (!m_length) return;

      m_sblock_header = (long *)malloc(n_sblocks * 256L * sizeof(long));

      //------------------------------------------------------------------------
      // STEP 1: compute rare trunk size, m_cblock_header and m_cblock_mapping.
      //------------------------------------------------------------------------
      m_cblock_header = (long *)malloc(n_cblocks * sizeof(long));
      m_cblock_header2 = (unsigned long *)malloc(k_sigma * n_cblocks * sizeof(unsigned long));
      m_cblock_mapping = (unsigned char *)malloc(n_cblocks * 512L);
      long rare_trunk_size = 0;

      long *list_beg = new long[k_sigma];      // beginnings of occurrence lists
      std::vector<long> *occ = new std::vector<long>[k_sigma];

      long cblock_type_bytes = (n_cblocks + 7) / 8;
      m_cblock_type = (unsigned char *)malloc(cblock_type_bytes);
      std::fill(m_cblock_type, m_cblock_type + cblock_type_bytes, 0);


      m_freq_trunk = (unsigned *)malloc(n_cblocks * k_cblock_size * sizeof(unsigned));


      bool isfreq[k_sigma];
      std::vector<std::pair<uint32_t, unsigned char> > sorted_chars;
      std::vector<unsigned char> freq_chars;
      std::vector<unsigned char> rare_chars;
      unsigned cblock_count[k_sigma];
      for (long cblock_id = 0; cblock_id < n_cblocks; ++cblock_id) {
        long cblock_beg = cblock_id << k_cblock_size_log;
        long cblock_end = cblock_beg + k_cblock_size;

        // Compute symbol counts inside a block.
        std::fill(cblock_count, cblock_count + k_sigma, 0L);
        for (long j = cblock_beg; j < cblock_end; ++j) {
          unsigned char c = (j < m_length ? text[j] : 0);
          ++cblock_count[c];
        }
        
        // Compute starting positions of occurrences lists.
        for (long j = 0, t, s = 0; j < k_sigma; ++j) {
          t = cblock_count[j];
          list_beg[j] = s;
          s += t;
        }
        
        // Compute first part the cblock header: ranks up to cblock
        // beginning pointers to beginnings of occurrence lists.
        for (long c = 0; c < 256; ++c) {
          m_cblock_header2[(cblock_id << 8) + c] = (list_beg[c] << 5);
          m_cblock_header2[(cblock_id << 8) + c] |= (m_count[c] << (k_cblock_size_log + 6));
        }

        // Sort symbol counts by frequencies.
        sorted_chars.clear();
        for (long j = 0; j < 256; ++j)
          if (cblock_count[j])
            sorted_chars.push_back(std::make_pair(cblock_count[j], j));
        std::sort(sorted_chars.begin(), sorted_chars.end());

        // Separate (at most, due to rounding of freq_cnt)
        // about 3% of rarest symbols.
        unsigned rare_cnt = 0L, rare_sum = 0L;
        while (rare_cnt < sorted_chars.size() &&
            16L * (rare_sum + sorted_chars[rare_cnt].first) <= k_cblock_size)
          rare_sum += sorted_chars[rare_cnt++].first;

        // Compute freq_cnt. Then round up freq_cnt + 1 (+1 is
        // for rare char marker) to the smallest power of two.
        // Note: rare_cnt > 0, so after rounding freq_cnt <= 256.
        unsigned freq_cnt = sorted_chars.size() - rare_cnt;
        unsigned freq_cnt_log = utils::log2ceil(freq_cnt + 1);
        freq_cnt = (1 << freq_cnt_log);

        if (freq_cnt >= 128) {
          // Set that the cblock is encoded using method for random strings.
          m_cblock_type[cblock_id >> 3] |= (1 << (cblock_id & 7));
 
          // Compute lists of occurrences.
          for (long c = 0; c < k_sigma; ++c) occ[c].clear();
          unsigned *cblock_trunk = m_freq_trunk + cblock_beg;
          std::fill(cblock_trunk, cblock_trunk + k_cblock_size, 0);

          for (long i = cblock_beg; i < cblock_end; ++i) {
            unsigned char c = (i < m_length ? text[i] : 0);
            long cblock_i = i - cblock_beg;
            occ[c].push_back(cblock_i);
          }

          for (long c = 0; c < k_sigma; ++c) {
            long freq = occ[c].size();

            // Compute the number of bits necessary to encode lookup table entries.
            // Lookup table can store values in the range [0..freq], thus we need
            // ceil(log2(freq + 1)) bits. This values is called lookup_bits and we
            // store it in the block header using 5 bits.
            long lookup_bits = utils::log2ceil(freq + 1);
            m_cblock_header2[(cblock_id << 8) + c] |= lookup_bits;

            // Compute log2 of the distance between two reference points.
            long refpoint_dist_log = 31 - lookup_bits;
            long refpoint_dist = (1L << refpoint_dist_log);
            long refpoint_dist_mask = refpoint_dist - 1;

            // Used to compute the closest reference point on the left.
            long refpoint_dist_mask_neg = (~refpoint_dist_mask);

            // For each block compute:
            //   * its boundaries (begin and end)
            //   * lookup table at entry corresponding to the block
            //   * store lists of occurrences of all symbols. Each
            //     value in the occurrence list is a distance to the
            //     reference block associated with the block containing
            //     the position. Such reference point is defined as the
            //     closest reference point to the left of the block begin
            long occ_ptr = 0;
            for (long j = 0; j < freq; ++j) {
              // Process j-th block.

              // 1
              //
              // Compute block boundaries.

              // The beginning of j-th block is the smallest i such that
              // floor((i * freq) / k_cblock_size) = j. We first overshoot
              // a little, and then go back; probably there exists exact
              // formula.
              long block_beg = ((j << k_cblock_size_log) + freq - 1) / freq;
              while (block_beg && (((block_beg - 1) * freq) >> k_cblock_size_log) == j)
                --block_beg;

              // Block end = beginning of the next block, or
              // k_cblock_size if this is already last block.
              long block_end = k_cblock_size;
              if (j + 1 != freq) {
                block_end = (((j + 1) << k_cblock_size_log) + freq - 1) / freq;
                while (block_end && (((block_end - 1) * freq) >> k_cblock_size_log) == j + 1)
                  --block_end;
              }

              // 2
              //
              // Find the range of elements from the current block inside occ[c].
              while (occ_ptr < freq && occ[c][occ_ptr] < block_beg) ++occ_ptr;
              long range_beg = occ_ptr;
              while (occ_ptr < freq && occ[c][occ_ptr] < block_end) ++occ_ptr;
              long range_end = occ_ptr;

              // 3
              //
              // Store the value in the lookup table.
              cblock_trunk[list_beg[c] + j] |= range_beg;

              // 4
              //
              // Add the occurrences occ[c][range_beg..range_end) to the
              // list of c's occurrences in the trunk. Encode them with
              // respect to the closes reference point on the left of block_beg.
              long closest_ref_point = (block_beg & refpoint_dist_mask_neg);
              for (long occ_id = range_beg; occ_id < range_end; ++occ_id)
                cblock_trunk[list_beg[c] + occ_id] |= ((occ[c][occ_id] - closest_ref_point) << lookup_bits);
            }
          }
        } else {
          // Recompute rare_cnt (note the +1).
          rare_cnt = 0;
          if (sorted_chars.size() + 1 > freq_cnt)
            rare_cnt = sorted_chars.size() + 1 - freq_cnt;

          // Compute freq and rare chars.
          rare_chars.clear();
          freq_chars.clear();
          for (unsigned i = 0; i < rare_cnt; ++i)
            rare_chars.push_back(sorted_chars[i].second);
          for (unsigned i = rare_cnt; i < sorted_chars.size(); ++i)
            freq_chars.push_back(sorted_chars[i].second);
 
          // If there are rare symbols, round up
          // rare_cnt to the smallest power of two.
          unsigned rare_cnt_log = 0;
          if (rare_cnt) {
            rare_cnt_log = utils::log2ceil(rare_cnt);
            rare_cnt = (1 << rare_cnt_log);
          }

          // Store log2 of rare_cnt and freq_cnt into block header.
          m_cblock_header[cblock_id] = freq_cnt_log;
          m_cblock_header[cblock_id] |= (rare_cnt_log << 8);

          // Compute and store symbols mapping.
          std::sort(freq_chars.begin(), freq_chars.end());
          std::sort(rare_chars.begin(), rare_chars.end());
          std::fill(isfreq, isfreq + 256, false);
          for (unsigned c = 0; c < 256; ++c)
            m_cblock_mapping[2 * (c * n_cblocks + cblock_id)] = k_char_type_missing;
          for (unsigned i = 0; i < freq_chars.size(); ++i) {
            unsigned char c = freq_chars[i];
            isfreq[c] = true;
            m_cblock_mapping[2 * (c * n_cblocks + cblock_id) + 1] = i;
            m_cblock_mapping[2 * (c * n_cblocks + cblock_id)] = k_char_type_freq;
          }
          for (unsigned i = 0; i < rare_chars.size(); ++i) {
            unsigned char c = rare_chars[i];
            m_cblock_mapping[2 * (c * n_cblocks + cblock_id) + 1] = i;
            m_cblock_mapping[2 * (c * n_cblocks + cblock_id)] = k_char_type_rare;
          }

          // Update rare_trunk_size.
          if (rare_cnt) {
            unsigned nofreq_cnt = 0L;
            for (long i = cblock_beg; i < cblock_end; ++i) {
              unsigned char c = (i < m_length ? text[i] : 0);
              if (!isfreq[c]) ++nofreq_cnt;
            }
            long rare_micro_blocks = 1 + (nofreq_cnt + rare_cnt - 1) / rare_cnt;
            rare_trunk_size += rare_micro_blocks * rare_cnt;
          }
        }
        
        // Compute ranks at superblock boundary.
        long sblock_id = (cblock_beg >> k_sblock_size_log);
        if (!(cblock_beg & k_sblock_size_mask)) {
          long sblock_ptr = (sblock_id << 8);
          for (unsigned i = 0; i < 256; ++i)
            m_sblock_header[sblock_ptr++] = m_count[i];
        }

        // Update global counts.        
        for (long i = cblock_beg; i < cblock_end; ++i) {
          unsigned char c = (i < m_length ? text[i] : 0);
          ++m_count[c];
        }
      }
      m_count[0] -= n_cblocks * k_cblock_size - m_length;

      m_rare_trunk = (unsigned *)malloc(rare_trunk_size * sizeof(unsigned));

      //------------------------------------------------------------------------
      // STEP 2: compute sb_rank, freq_trunk and rare_trunk.
      //------------------------------------------------------------------------
      unsigned char freq_map[256];
      unsigned char rare_map[256];
      long r_filled = 0; // for rare trunk
      unsigned long *cur_count = new unsigned long[256];
      for (long cblock_id = 0; cblock_id < n_cblocks; ++cblock_id) {
        // If the cblock was encoded using method for random sequences,
        // skip this cblock and proceed to next. The trunk of this cblock
        // was already encoded in the first stage.
        if (m_cblock_type[cblock_id >> 3] & (1 << (cblock_id & 7)))
          continue;

        // Fill in cur_counts.
        for (long c = 0; c < 256; ++c)
          cur_count[c] = (m_cblock_header2[(cblock_id << 8) + c] >> (k_cblock_size_log + 6));
      
        long r_ptr = r_filled; // where to write next
        m_cblock_header[cblock_id] |= (r_filled << 16);

        long cblock_beg = cblock_id << k_cblock_size_log;
        long cblock_end = cblock_beg + k_cblock_size;

        //----------------------------------------------------------------------
        // Process cblock [cblock_beg..cblock_end).
        //----------------------------------------------------------------------

        long freq_cnt_log = (m_cblock_header[cblock_id] & 255L);
        long rare_cnt_log = ((m_cblock_header[cblock_id] >> 8) & 255L);
        long freq_cnt = (1L << freq_cnt_log);
        long freq_cnt_mask = freq_cnt - 1;
        long rare_cnt = (1L << rare_cnt_log);
        long rare_cnt_mask = rare_cnt - 1;

        freq_chars.clear();
        rare_chars.clear();
        std::fill(isfreq, isfreq + 256, false);
        for (unsigned j = 0; j < 256; ++j) {
          unsigned char type = m_cblock_mapping[2 * (j * n_cblocks + cblock_id)];
          if (type == k_char_type_freq) {
            isfreq[j] = true;
            freq_chars.push_back(j);
            freq_map[j] = m_cblock_mapping[2 * (j * n_cblocks + cblock_id) + 1];
          } else if (type == k_char_type_rare) {
            rare_chars.push_back(j);
            rare_map[j] = m_cblock_mapping[2 * (j * n_cblocks + cblock_id) + 1];
            freq_map[j] = freq_cnt - 1;
          }
        }

        if (rare_chars.empty()) {
          rare_cnt_log = 0;
          rare_cnt = 0;
        }

        long sblock_id = (cblock_beg >> k_sblock_size_log);

        // Compute freq and rare trunk of current block.
        long nofreq_cnt = 0;
        for (long i = cblock_beg; i < cblock_end; ++i) {
          unsigned char c = (i < m_length ? text[i] : 0);

          //--------------------------------------------------------------------
          // Invariant: for any symbol a, c_rank[a] = number of occurrence of
          // symbols a in prefix text[0..i).
          //--------------------------------------------------------------------
          // Invariant: for any symbol a, sb_rank[(sb_id << 8) + a] = the
          // number of occurrence of a in a prefix of text up to the closest
          // superblock boundary. All ranks we store in the trunk are relative
          // to this boundary thus we compute them as
          // c_rank[a] - sb_rank[(sb_id << 8) + a].
          //--------------------------------------------------------------------

          // Compute ranks at the freq microblock boundary.
          if (!(i & freq_cnt_mask)) {
            m_freq_trunk[i + freq_cnt - 1] = (nofreq_cnt << 8);
            for (long j = 0; j + 1 < freq_cnt; ++j) {
              unsigned char ch = (j < (long)freq_chars.size() ? freq_chars[j] : 0);
              long local_rank = cur_count[ch] - m_sblock_header[(sblock_id << 8) + ch];
              m_freq_trunk[i + j] = (local_rank << 8);
            }
          }

          // Store freq symbol mapping.
          m_freq_trunk[i] |= freq_map[c];

          // Handle rare symbol.
          if (!isfreq[c]) {
            // Compute ranks at the rare microblock boundary.
            if (!(nofreq_cnt & rare_cnt_mask)) {
              for (long j = 0; j < rare_cnt; ++j) {
                unsigned char ch = (j < (long)rare_chars.size() ? rare_chars[j] : 0);
                long local_rank = cur_count[ch] - m_sblock_header[(sblock_id << 8) + ch];
                m_rare_trunk[r_filled++] = (local_rank << 8);
              }
            }

            // Store rare symbol mapping.
            m_rare_trunk[r_ptr++] |= rare_map[c];
            ++nofreq_cnt;
          }

          ++cur_count[c];
        }

        for (long j = 0; j < rare_cnt; ++j) {
          unsigned char ch = (j < (long)rare_chars.size() ? rare_chars[j] : 0);
          long local_rank = cur_count[ch] - m_sblock_header[(sblock_id << 8) + ch];
          m_rare_trunk[r_filled++] = (local_rank << 8);
        }
      }

      delete[] cur_count;
      delete[] list_beg;
      delete[] occ;
    }

    inline long rank(long i, unsigned char c) {
      if (i <= 0) return 0L;
      else if (i >= m_length) return m_count[c];

      long cblock_id = (i >> k_cblock_size_log);    
      if (m_cblock_type[cblock_id >> 3] & (1 << (cblock_id & 7))) {    
        long cblock_beg = (i & k_cblock_size_mask_neg);
        long cblock_i = (i & k_cblock_size_mask);     // offset in cblock
      
        //--------------------------------------------------------------------------
        // STEP 1: extract the rank up to the start of cblock.
        //--------------------------------------------------------------------------
        long rank_up_to_cblock = (m_cblock_header2[(cblock_id << k_sigma_log) + c] >> (k_cblock_size_log + 6));

        //--------------------------------------------------------------------------
        // STEP 2: compute the number of occurrences of c inside the cblock.
        //--------------------------------------------------------------------------

        // 1
        //
        // Decode the beginning and end of c's occurrence list.
        long list_beg = ((m_cblock_header2[(cblock_id << k_sigma_log) + c] >> 5) & k_2cblock_size_mask);
        long list_end = ((c == k_sigma - 1) ? k_cblock_size :
            ((m_cblock_header2[(cblock_id << k_sigma_log) + c + 1] >> 5) & k_2cblock_size_mask));
        if (list_beg == list_end) return rank_up_to_cblock;

        // 2
        //
        // Compute the distance from i to the closest reference point on the left.
        long lookup_bits = (m_cblock_header2[(cblock_id << k_sigma_log) + c] & 31);
        long refpoint_dist_log = 31 - lookup_bits;
        long refpoint_disk_mask = (1 << refpoint_dist_log) - 1;
        long i_refpoint_offset = (cblock_i & refpoint_disk_mask);

        // 3
        //
        // Compute threshold of symbol c inside the current cblock.
        // The threshold is a small power of two that is >= max block size,
        // where max block size is the maximal block size for symbol c in
        // the current cblock.
        long threshold = (1 << (k_cblock_size_log - lookup_bits + 1));

        // 4
        //
        // Compute the id if block containing i. Note: we only have id.
        // Computing boundaries would require division.
        long list_size = list_end - list_beg;
        long approx = ((cblock_i * list_size) >> k_cblock_size_log);

        // 5
        //
        // Extract the lookup table entry.
        long lookup_mask = (1 << lookup_bits) - 1;
        long begin = (m_freq_trunk[cblock_beg + list_beg + approx] & lookup_mask);
        long next_block_begin =  (approx + 1 == list_size) ? list_size : (m_freq_trunk[cblock_beg + list_beg + approx + 1] & lookup_mask);

        // 6
        //
        // Correct the value of begin and return the answer. To correct begin
        // we need to know wrt to which reference point have we encoded the
        // occurrences of c inside the block. This would be easy to compute
        // if we knew the block beginning (but we don't). However, we can
        // infer the unknown reference point from i and/or the first element
        // in the occurrence list (or rather whether the first element is big
        // or small).
        if (i_refpoint_offset >= threshold) {
          // This is the easy case which will happen most of the
          // time, so we should get good branch prediction here.
          while (begin < next_block_begin && (m_freq_trunk[cblock_beg + list_beg + begin] >> lookup_bits) < i_refpoint_offset)
            ++begin;

          return rank_up_to_cblock + begin;
        } else {
          // This is executed very rarely, so we move more
          // expensive code (i.e., another branch) here.
          if (begin == next_block_begin || (m_freq_trunk[cblock_beg + list_beg + begin] >> lookup_bits) < (2 * threshold)) {
            // The value in the occ list was small -> the ref
            // point for i and for the block are the same, we
            // proceed as before, without modifying i_refpoint_offset.
            while (begin < next_block_begin && (m_freq_trunk[cblock_beg + list_beg + begin] >> lookup_bits) < i_refpoint_offset)
              ++begin;

            return rank_up_to_cblock + begin;
          } else {
            // Block occurrences were encoded wrt to the previous
            // ref point -> we increase i_refpoint_offset by
            // refpoint_dist and proceed as before.
            i_refpoint_offset += (1 << refpoint_dist_log);
            while (begin < next_block_begin && (m_freq_trunk[cblock_beg + list_beg + begin] >> lookup_bits) < i_refpoint_offset)
              ++begin;

            return rank_up_to_cblock + begin;
          }
        }
      } else {    
        long sblock_id = (i >> k_sblock_size_log);
        long sblock_rank = m_sblock_header[(sblock_id << 8) + c];

        unsigned char type = m_cblock_mapping[2 * (c * n_cblocks + cblock_id)];
        unsigned char c_map = m_cblock_mapping[2 * (c * n_cblocks + cblock_id) + 1];

        long freq_cnt_bits = (m_cblock_header[cblock_id] & 255L);
        long rare_cnt_bits = ((m_cblock_header[cblock_id] >> 8) & 255L);
        long block_id = (i >> freq_cnt_bits);

        if (type == k_char_type_freq) {
          // Case 1 (fastest): symbol c was frequent in the context block.
          // Answer a query using frequent trunk.
          long block_rank = m_freq_trunk[(block_id << freq_cnt_bits) + c_map] >> 8;
          long extra = 0;
          for (long j = (block_id << freq_cnt_bits); j < i; ++j)
            if ((m_freq_trunk[j] & 255) == c_map) ++extra;

          return sblock_rank + block_rank + extra;
        } else if (type == k_char_type_rare) {
          // Case 2: symbol c was rare inside the context block.
          // Compute new_i.
          long rare_trunk_ptr = (m_cblock_header[cblock_id] >> 16);
          long new_i = m_freq_trunk[((block_id + 1) << freq_cnt_bits) - 1] >> 8;
          for (long j = (block_id << freq_cnt_bits); j < i; ++j)
            if ((m_freq_trunk[j] & 255) + 1 == (1U << freq_cnt_bits)) ++new_i;
      
          // Answer a query on rare trunk.
          long rare_block_id = (new_i >> rare_cnt_bits);
          long block_rank = m_rare_trunk[rare_trunk_ptr +
            (rare_block_id << rare_cnt_bits) + c_map] >> 8;
          long extra = 0;
          for (long j = (rare_block_id << rare_cnt_bits); j < new_i; ++j)
            if ((m_rare_trunk[rare_trunk_ptr + j] & 255) == c_map) ++extra;

          return sblock_rank + block_rank + extra;
        } else {
          // Case 3: symbol c does not occur in the context block.
          // Find the first cblock where c occurrs.
          while (cblock_id < n_cblocks && (cblock_id & k_cblocks_in_sblock_mask) &&
              m_cblock_mapping[2 * (c * n_cblocks + cblock_id)] == k_char_type_missing)
            ++cblock_id;

          if (cblock_id == n_cblocks) {
            // We reached the end of encoding, return count[c].
            return m_count[c];
          } else if (!(cblock_id & k_cblocks_in_sblock_mask)) {
            // We reached the boundary of superblock,
            // retreive the answer from superblock header.
            return m_sblock_header[256 * (cblock_id >> k_cblocks_in_sblock_log) + c];
          } else {
            // We found cblock where c occurrs, but it wasn't on the
            // sblock boundary. In the recursive call this will either
            // be case 1 or case 2.
            return rank(cblock_id << k_cblock_size_log, c);
          }
        }
      }
    }

    ~rank4n() {
      if (m_length) {
        free(m_sblock_header);
        free(m_cblock_header);
        free(m_cblock_header2);
        free(m_cblock_mapping);
        free(m_cblock_type);
        free(m_freq_trunk);
        free(m_rare_trunk);
      }
      free(m_count);
    }
};


template<unsigned k_sblock_size_log, unsigned k_cblock_size_log, unsigned k_sigma_log>
  const unsigned rank4n<k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_cblock_size = (1 << k_cblock_size_log);

template<unsigned k_sblock_size_log, unsigned k_cblock_size_log, unsigned k_sigma_log>
  const unsigned rank4n<k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_cblock_size_mask = (1 << k_cblock_size_log) - 1;
  
template<unsigned k_sblock_size_log, unsigned k_cblock_size_log, unsigned k_sigma_log>
  const unsigned rank4n<k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_2cblock_size = (2 << k_cblock_size_log);

template<unsigned k_sblock_size_log, unsigned k_cblock_size_log, unsigned k_sigma_log>
  const unsigned rank4n<k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_2cblock_size_mask = (2 << k_cblock_size_log) - 1;

template<unsigned k_sblock_size_log, unsigned k_cblock_size_log, unsigned k_sigma_log>
  const int rank4n<k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_sigma = (1 << k_sigma_log);

template<unsigned k_sblock_size_log, unsigned k_cblock_size_log, unsigned k_sigma_log>
  const int rank4n<k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_sigma_mask = (1 << k_sigma_log) - 1;

template<unsigned k_sblock_size_log, unsigned k_cblock_size_log, unsigned k_sigma_log>
  const unsigned rank4n<k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_cblock_size_mask_neg = ~((1 << k_cblock_size_log) - 1);

template<unsigned k_sblock_size_log, unsigned k_cblock_size_log, unsigned k_sigma_log>
  const unsigned rank4n<k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_cblocks_in_sblock_log = k_sblock_size_log - k_cblock_size_log;

template<unsigned k_sblock_size_log, unsigned k_cblock_size_log, unsigned k_sigma_log>
  const unsigned rank4n<k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_cblocks_in_sblock = (1 << (k_sblock_size_log - k_cblock_size_log));

template<unsigned k_sblock_size_log, unsigned k_cblock_size_log, unsigned k_sigma_log>
  const unsigned rank4n<k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_cblocks_in_sblock_mask = (1 << (k_sblock_size_log - k_cblock_size_log)) - 1;

template<unsigned k_sblock_size_log, unsigned k_cblock_size_log, unsigned k_sigma_log>
  const unsigned rank4n<k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_sblock_size = (1 << k_sblock_size_log);
    
template<unsigned k_sblock_size_log, unsigned k_cblock_size_log, unsigned k_sigma_log>
  const unsigned rank4n<k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_sblock_size_mask = (1 << k_sblock_size_log) - 1;


#endif // __RANK4N_H_INCLUDED
