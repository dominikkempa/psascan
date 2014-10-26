#ifndef __RANK4N_H_INCLUDED
#define __RANK4N_H_INCLUDED

#include <algorithm>
#include <vector>

#include "utils.h"


//=============================================================================
// Data structure answering rank queries over byte alphabet.
// For an input sequence of length n the data structure occupies ~4.1n bytes.
//=============================================================================
template<unsigned k_sblock_size_log = 24,
  unsigned k_cblock_size_log = 20,
  unsigned k_sigma_log = 8>
class rank4n {
  private:
    static const unsigned long k_cblock_size;
    static const unsigned long k_cblock_size_mask;
    static const unsigned long k_cblock_size_mask_neg;
    static const unsigned k_cblocks_in_sblock_log;
    static const unsigned k_cblocks_in_sblock;
    static const unsigned k_cblocks_in_sblock_mask;
    static const unsigned k_2cblock_size;
    static const unsigned k_2cblock_size_mask;
    static const unsigned k_sblock_size;
    static const unsigned k_sblock_size_mask;
    static const unsigned k_sigma;
    static const unsigned k_sigma_mask;

    static const unsigned k_char_type_freq =    0x01;
    static const unsigned k_char_type_rare =    0x02;
    static const unsigned k_char_type_missing = 0x03;

    unsigned long m_length;   // length of original sequence
    unsigned long n_cblocks;  // number of context blocks
    unsigned long n_sblocks;  // number of super blocks

    unsigned long *m_sblock_header;
    unsigned long *m_cblock_header;
    unsigned long *m_cblock_header2;

    unsigned char *m_cblock_type;
    unsigned char *m_cblock_mapping;

    unsigned *m_freq_trunk;
    unsigned *m_rare_trunk;

  public:
    unsigned long *m_count; // symbol counts

  public:
    rank4n(const unsigned char *text, unsigned long length) {
      m_length = length;
      n_cblocks = (m_length + k_cblock_size - 1) / k_cblock_size;
      n_sblocks = (n_cblocks + k_cblocks_in_sblock - 1) / k_cblocks_in_sblock;

      m_count = (unsigned long *)malloc(256L * sizeof(unsigned long));
      std::fill(m_count, m_count + 256, 0UL);
      if (!m_length) return;

      m_sblock_header = (unsigned long *)malloc(n_sblocks * sizeof(unsigned long) * k_sigma);
      m_cblock_header = (unsigned long *)malloc(n_cblocks * sizeof(unsigned long));
      m_cblock_header2 = (unsigned long *)malloc(n_cblocks * k_sigma * sizeof(unsigned long));
      m_cblock_mapping = (unsigned char *)malloc(n_cblocks * k_sigma * 2);
      m_cblock_type = (unsigned char *)malloc((n_cblocks + 7) / 8);
      m_freq_trunk = (unsigned *)calloc(n_cblocks * k_cblock_size, sizeof(unsigned));
      std::fill(m_cblock_type, m_cblock_type + (n_cblocks + 7) / 8, 0);

      encode_type_I(text);
      encode_type_II(text);

      m_count[0] -= n_cblocks * k_cblock_size - m_length;  // remove extra zeros
    }

    void encode_type_I(const unsigned char *text) {
      std::vector<unsigned char> freq_chars;
      std::vector<unsigned char> rare_chars;
      std::vector<std::pair<uint32_t, unsigned char> > sorted_chars;

      unsigned *occ = (unsigned *)malloc((k_cblock_size + 1) * sizeof(unsigned));
      unsigned *refpoint_precomputed = (unsigned *)malloc(k_cblock_size * sizeof(unsigned));
      unsigned *cblock_count = new unsigned[k_sigma];
      unsigned *list_beg = new unsigned[k_sigma];
      unsigned *list_beg2 = new unsigned[k_sigma];
      unsigned *lookup_bits_precomputed = new unsigned[k_sigma];
      unsigned *min_block_size_precomputed = new unsigned[k_sigma];
      unsigned long *refpoint_mask_precomputed = new unsigned long[k_sigma];
      bool *isfreq = new bool[k_sigma];

      // Process cblocks one by one.
      long type_II_cnt = 0;
      unsigned long rare_trunk_total_size = 0L;
      for (unsigned long cblock_id = 0; cblock_id < n_cblocks; ++cblock_id) {
        unsigned long cblock_beg = cblock_id << k_cblock_size_log;
        unsigned long cblock_end = cblock_beg + k_cblock_size;

        // Store ranks up to cblock boundary.
        unsigned long cblock_header_beg = (cblock_id << k_sigma_log);
        for (unsigned c = 0; c < k_sigma; ++c)
          m_cblock_header2[cblock_header_beg + c] = (m_count[c] << (k_cblock_size_log + 6));

        // Update sblock header
        if (!(cblock_beg & k_sblock_size_mask)) {
          unsigned long sblock_id = (cblock_beg >> k_sblock_size_log);
          unsigned long sblock_header_beg = (sblock_id << k_sigma_log);
          for (unsigned c = 0; c < k_sigma; ++c)
            m_sblock_header[sblock_header_beg + c] = m_count[c];
        }

        // Compute symbol counts inside cblock.
        std::fill(cblock_count, cblock_count + k_sigma, 0);
        unsigned long maxj = std::min(cblock_end, m_length);
        for (unsigned long j = cblock_beg; j < maxj; ++j)
          ++cblock_count[text[j]];
        cblock_count[0] += cblock_end - maxj;

        // Compute starting positions of occurrences lists and update m_count.
        for (unsigned c = 0, t, s = 0; c < k_sigma; ++c) {
          t = cblock_count[c];
          m_count[c] += t;
          list_beg[c] = s;
          list_beg2[c] = s;
          s += t;
        }

        // Store pointers to beginnings of occurrence lists in the type-II
        // cblock header. Note: this implicitly encodes cblock counts.
        for (unsigned c = 0; c < k_sigma; ++c)
          m_cblock_header2[cblock_header_beg + c] |= (list_beg[c] << 5);

        // Sort symbol counts by frequencies.
        sorted_chars.clear();
        for (unsigned j = 0; j < k_sigma; ++j)
          if (cblock_count[j])
            sorted_chars.push_back(std::make_pair(cblock_count[j], j));
        std::sort(sorted_chars.begin(), sorted_chars.end());

        // Separate (at most, due to rounding of freq_cnt) about 3% of
        // rarest symbols.
        unsigned rare_cnt = 0L, rare_sum = 0L;
        while (rare_cnt < sorted_chars.size() &&
            16L * (rare_sum + sorted_chars[rare_cnt].first) <= k_cblock_size)
          rare_sum += sorted_chars[rare_cnt++].first;

        // Compute freq_cnt. Then round up freq_cnt + 1 (+1 is for rare char
        // marker) to the smallest power of two. Note: rare_cnt > 0, so after
        // rounding freq_cnt <= 256.
        unsigned freq_cnt = sorted_chars.size() - rare_cnt;
        unsigned freq_cnt_log = utils::log2ceil(freq_cnt + 1);
        freq_cnt = (1 << freq_cnt_log);

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

        // If there are rare symbols, round up rare_cnt to the smallest power
        // of two.
        unsigned rare_cnt_log = 0;
        if (rare_cnt) {
          rare_cnt_log = utils::log2ceil(rare_cnt);
          rare_cnt = (1 << rare_cnt_log);
        }

        // Update cblock type-I header.
        m_cblock_header[cblock_id] = freq_cnt_log;
        m_cblock_header[cblock_id] |= (rare_cnt_log << 8);
        m_cblock_header[cblock_id] |= (rare_trunk_total_size << 16);

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

        unsigned nofreq_cnt = 0L;
        for (unsigned c = 0; c < k_sigma; ++c)
          if (!isfreq[c]) nofreq_cnt += cblock_count[c];

        if (freq_cnt >= 128) { // Type-I cblock.
          m_cblock_type[cblock_id >> 3] |= (1 << (cblock_id & 7));
 
          // Compute lists of occurrences.
          for (unsigned long i = cblock_beg; i < maxj; ++i)
            occ[list_beg2[text[i]]++] = i - cblock_beg;
          for (unsigned long i = maxj; i < cblock_end; ++i)
            occ[list_beg2[0]++] = i - cblock_beg;

          // Precompute helper arrays and and store lookup bits into the header.
          for (unsigned c = 0; c < k_sigma; ++c) {
            lookup_bits_precomputed[c] = utils::log2ceil(cblock_count[c] + 2);
            m_cblock_header2[(cblock_id << 8) + c] |= lookup_bits_precomputed[c];
            if (cblock_count[c])
              min_block_size_precomputed[c] = k_cblock_size / cblock_count[c];
            else min_block_size_precomputed[c] = 0;

            unsigned refpoint_dist_log = 31 - lookup_bits_precomputed[c];
            unsigned long refpoint_dist = (1UL << refpoint_dist_log);
            unsigned long refpoint_dist_mask = refpoint_dist - 1;
            unsigned long refpoint_dist_mask_neg = (~refpoint_dist_mask);
            refpoint_mask_precomputed[c] = refpoint_dist_mask_neg;
          }

          // Actual encoding follows.
          unsigned *cblock_trunk = m_freq_trunk + cblock_beg;
          for (unsigned c = 0; c < k_sigma; ++c) {
            unsigned freq = cblock_count[c];
            unsigned min_block_size = min_block_size_precomputed[c];
            unsigned lookup_bits = lookup_bits_precomputed[c];
            unsigned refpoint_dist_mask_neg = refpoint_mask_precomputed[c];
            unsigned c_list_beg = list_beg[c];

            for (unsigned j = 0; j < freq; ++j)
              cblock_trunk[c_list_beg + j] = freq + 1;
            if (freq) cblock_trunk[c_list_beg + freq - 1] = freq;

            unsigned block_beg = 0;
            for (unsigned j = 0; j < freq; ++j) {
              refpoint_precomputed[j] = (block_beg & refpoint_dist_mask_neg);
              block_beg += min_block_size;
              if ((((unsigned long)block_beg * freq) >> k_cblock_size_log) == j) ++block_beg;
            }

            unsigned refpoint, block_id;
            unsigned mask = (~((1UL << lookup_bits) - 1));
            if (freq) {
              for (long j = freq - 1; j >= 0; --j) {
                block_id = (((unsigned long)occ[c_list_beg + j] * freq) >> k_cblock_size_log);
                refpoint = refpoint_precomputed[block_id];
                cblock_trunk[c_list_beg + block_id] &= mask;
                cblock_trunk[c_list_beg + block_id] |= (unsigned)j;
                cblock_trunk[c_list_beg + j] |= ((occ[c_list_beg + j] - refpoint) << lookup_bits);
              }
            }
          }
        } else {
          ++type_II_cnt;
          // Update total rare trunk size.
          if (rare_cnt) {
            long rare_blocks = 1 + (nofreq_cnt + rare_cnt - 1) / rare_cnt;
            rare_trunk_total_size += rare_blocks * rare_cnt;
          }
        }
      }

      // Allocate rare trunk.
      m_rare_trunk = NULL;
      if (rare_trunk_total_size)
        m_rare_trunk = (unsigned *)calloc(rare_trunk_total_size, sizeof(unsigned));

      // Clean up.
      delete[] list_beg;
      delete[] list_beg2;
      delete[] isfreq;
      delete[] cblock_count;
      delete[] lookup_bits_precomputed;
      delete[] min_block_size_precomputed;
      delete[] refpoint_mask_precomputed;
      free(refpoint_precomputed);
      free(occ);
    }

    void encode_type_II(const unsigned char *text) {
      unsigned char *freq_map = new unsigned char[k_sigma];
      unsigned char *rare_map = new unsigned char[k_sigma];
      unsigned long *cur_count = new unsigned long[k_sigma];
      unsigned long *off = new unsigned long[k_sigma];

      long *sblock_h = new long[k_sigma];
      int *israre = new int[k_sigma];

      std::vector<unsigned char> freq_chars;
      std::vector<unsigned char> rare_chars;

      for (unsigned long cblock_id = 0; cblock_id < n_cblocks; ++cblock_id) {
        unsigned long cblock_beg = cblock_id << k_cblock_size_log;
        unsigned long cblock_end = cblock_beg + k_cblock_size;

        // Skip the cblock if it was type-I encoded.
        if (m_cblock_type[cblock_id >> 3] & (1 << (cblock_id & 7))) continue;
 
        // Retreive symbol counts up to this cblock begin and
        // pointer to rare trunk size from cblock headers.
        for (unsigned c = 0; c < k_sigma; ++c)
          cur_count[c] = (m_cblock_header2[(cblock_id << 8) + c] >> (k_cblock_size_log + 6));

        long r_filled  = (m_cblock_header[cblock_id] >> 16);
        long r_ptr = r_filled;

        long freq_cnt_log = (m_cblock_header[cblock_id] & 255L);
        long rare_cnt_log = ((m_cblock_header[cblock_id] >> 8) & 255L);
        long freq_cnt = (1L << freq_cnt_log);
        long rare_cnt = (1L << rare_cnt_log);
        long rare_cnt_mask = rare_cnt - 1;

        freq_chars.clear();
        rare_chars.clear();
        std::fill(israre, israre + k_sigma, 1);
        for (unsigned c = 0; c < k_sigma; ++c) {
          unsigned char type = m_cblock_mapping[2 * (c * n_cblocks + cblock_id)];
          if (type == k_char_type_freq) {
            israre[c] = 0;
            freq_chars.push_back(c);
            freq_map[c] = m_cblock_mapping[2 * (c * n_cblocks + cblock_id) + 1];
          } else if (type == k_char_type_rare) {
            rare_chars.push_back(c);
            rare_map[c] = m_cblock_mapping[2 * (c * n_cblocks + cblock_id) + 1];
            freq_map[c] = freq_cnt - 1;
          }
        }

        if (rare_chars.empty()) {
          rare_cnt_log = 0;
          rare_cnt = 0;
        }


        long sblock_id = (cblock_beg >> k_sblock_size_log);
        std::copy(m_sblock_header + (sblock_id << 8), m_sblock_header + (sblock_id << 8) + k_sigma, sblock_h);
        for (long j = 0; j < k_sigma; ++j) off[j] = cur_count[j] - sblock_h[j];

        long nofreq_cnt = 0;
        long freq_chars_size = (long)freq_chars.size();
        long rare_chars_size = (long)rare_chars.size();
        if (cblock_end <= m_length) {
          for (unsigned long i = cblock_beg; i < cblock_end; i += freq_cnt) {
            for (long j = 0; j < freq_chars_size; ++j) {
              unsigned char ch = freq_chars[j];
              m_freq_trunk[i + j] = (off[ch] << 8);
            }
            m_freq_trunk[i + freq_cnt - 1] = (nofreq_cnt << 8);
            for (unsigned long j = i; j < i + freq_cnt; ++j) {
              unsigned char c = text[j];
              m_freq_trunk[j] |= freq_map[c];
              if (israre[c]) {
                if (!(nofreq_cnt & rare_cnt_mask)) {
                  for (long jj = 0; jj < rare_chars_size; ++jj) {
                    unsigned char ch = rare_chars[jj];
                    m_rare_trunk[r_filled++] = (off[ch] << 8);
                  }
                  r_filled += rare_cnt - rare_chars_size;
                }
                m_rare_trunk[r_ptr++] |= rare_map[c];
              }
              ++off[c];
              nofreq_cnt += israre[c];
            }
          }
          for (long i = 0; i < k_sigma; ++i)
            cur_count[i] = sblock_h[i] + off[i];
        } else {
          for (unsigned long i = cblock_beg; i < cblock_end; i += freq_cnt) {
            for (long j = 0; j < freq_chars_size; ++j) {
              unsigned char ch = freq_chars[j];
              m_freq_trunk[i + j] = (off[ch] << 8);
            }
            m_freq_trunk[i + freq_cnt - 1] = (nofreq_cnt << 8);
            for (unsigned long j = i; j < i + freq_cnt; ++j) {
              unsigned char c = (j < m_length ? text[j] : 0);
              m_freq_trunk[j] |= freq_map[c];
              if (israre[c]) {
                if (!(nofreq_cnt & rare_cnt_mask)) {
                  for (long jj = 0; jj < rare_chars_size; ++jj) {
                    unsigned char ch = rare_chars[jj];
                    m_rare_trunk[r_filled++] = (off[ch] << 8);
                  }
                  r_filled += rare_cnt - rare_chars_size;
                }
                m_rare_trunk[r_ptr++] |= rare_map[c];
              }
              ++off[c];
              nofreq_cnt += israre[c];
            }
          }
          for (long i = 0; i < k_sigma; ++i)
            cur_count[i] = sblock_h[i] + off[i];
        }


        for (long j = 0; j < rare_cnt; ++j) {
          unsigned char ch = (j < (long)rare_chars.size() ? rare_chars[j] : 0);
          long local_rank = cur_count[ch] - m_sblock_header[(sblock_id << 8) + ch];
          m_rare_trunk[r_filled++] = (local_rank << 8);
        }
      }

      delete[] cur_count;
      delete[] sblock_h;
      delete[] freq_map;
      delete[] rare_map;
      delete[] israre;
      delete[] off;
    }

    inline long rank(long i, unsigned char c) const {
      if (i <= 0) return 0L;
      else if ((unsigned long)i >= m_length) return m_count[c];

      unsigned long cblock_id = (i >> k_cblock_size_log);    
      if (m_cblock_type[cblock_id >> 3] & (1 << (cblock_id & 7))) {  // type-I cblock
        long cblock_beg = (i & k_cblock_size_mask_neg);
        long cblock_i = (i & k_cblock_size_mask);  // offset in cblock
      
        // Extract the rank up to the start of cblock.
        long rank_up_to_cblock = (m_cblock_header2[(cblock_id << k_sigma_log) + c] >> (k_cblock_size_log + 6));

        // Decode the beginning and end of c's occurrence list.
        long list_beg = ((m_cblock_header2[(cblock_id << k_sigma_log) + c] >> 5) & k_2cblock_size_mask);
        long list_end = ((c == k_sigma - 1) ? k_cblock_size :
            ((m_cblock_header2[(cblock_id << k_sigma_log) + c + 1] >> 5) & k_2cblock_size_mask));
        if (list_beg == list_end) return rank_up_to_cblock;

        // Compute the distance from i to the closest reference point on the left.
        long lookup_bits = (m_cblock_header2[(cblock_id << k_sigma_log) + c] & 31);
        long refpoint_dist_log = 31 - lookup_bits;
        long refpoint_disk_mask = (1 << refpoint_dist_log) - 1;
        long i_refpoint_offset = (cblock_i & refpoint_disk_mask);

        // Compute threshold of symbol c inside the current cblock.
        long threshold = (1 << (k_cblock_size_log - lookup_bits + 1));

        // Compute the id of block containing i.
        long list_size = list_end - list_beg;
        long approx = ((cblock_i * list_size) >> k_cblock_size_log);

        // Extract the lookup table entry.
        long lookup_mask = (1 << lookup_bits) - 1;
        long begin = (m_freq_trunk[cblock_beg + list_beg + approx] & lookup_mask);

        // Empty block optimization.
        if (begin == list_size + 1) {
          // Block containing cblock_i is empty, just find the beginning.
          ++approx;
          while ((m_freq_trunk[cblock_beg + list_beg + approx] & lookup_mask) == list_size + 1) ++approx;
          begin = (m_freq_trunk[cblock_beg + list_beg + approx] & lookup_mask);
          return rank_up_to_cblock + begin;
        }
        
        long next_block_begin =  (approx + 1 == list_size) ? list_size :
          (m_freq_trunk[cblock_beg + list_beg + approx + 1] & lookup_mask);

        // Correct next_block_begin
        if (approx + 1 != list_size && next_block_begin == list_size + 1) {
          ++approx;
          while ((m_freq_trunk[cblock_beg + list_beg + approx + 1] & lookup_mask) == list_size + 1) ++approx;
          next_block_begin = (m_freq_trunk[cblock_beg + list_beg + approx + 1] & lookup_mask);
        }

        // Correct the value of begin and return the answer.
        if (i_refpoint_offset >= threshold) {
          // Case 1: easy case, will happen most of the time.
          while (begin < next_block_begin && (m_freq_trunk[cblock_beg + list_beg + begin] >> lookup_bits) < i_refpoint_offset)
            ++begin;

          return rank_up_to_cblock + begin;
        } else {
          // Case 2: executed very rarely.
          if (begin == next_block_begin || (m_freq_trunk[cblock_beg + list_beg + begin] >> lookup_bits) < (2 * threshold)) {
            // Case 2a: the value in the occ list was small -> the ref
            // point for i and for the block are the same, we
            // proceed as before, without modifying i_refpoint_offset.
            while (begin < next_block_begin && (m_freq_trunk[cblock_beg + list_beg + begin] >> lookup_bits) < i_refpoint_offset)
              ++begin;

            return rank_up_to_cblock + begin;
          } else {
            // Case 2b: block occurrences were encoded wrt to the
            // previous ref point -> we increase i_refpoint_offset
            // by refpoint_dist and proceed as before.
            i_refpoint_offset += (1 << refpoint_dist_log);
            while (begin < next_block_begin && (m_freq_trunk[cblock_beg + list_beg + begin] >> lookup_bits) < i_refpoint_offset)
              ++begin;

            return rank_up_to_cblock + begin;
          }
        }
      } else {  // type-II cblock
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
        if (m_rare_trunk)
          free(m_rare_trunk);
      }

      free(m_count);
    }
};


template<unsigned k_sblock_size_log, unsigned k_cblock_size_log, unsigned k_sigma_log>
  const unsigned long rank4n<k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_cblock_size = (1 << k_cblock_size_log);

template<unsigned k_sblock_size_log, unsigned k_cblock_size_log, unsigned k_sigma_log>
  const unsigned long rank4n<k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_cblock_size_mask = (1 << k_cblock_size_log) - 1;
  
template<unsigned k_sblock_size_log, unsigned k_cblock_size_log, unsigned k_sigma_log>
  const unsigned rank4n<k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_2cblock_size = (2 << k_cblock_size_log);

template<unsigned k_sblock_size_log, unsigned k_cblock_size_log, unsigned k_sigma_log>
  const unsigned rank4n<k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_2cblock_size_mask = (2 << k_cblock_size_log) - 1;

template<unsigned k_sblock_size_log, unsigned k_cblock_size_log, unsigned k_sigma_log>
  const unsigned rank4n<k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_sigma = (1 << k_sigma_log);

template<unsigned k_sblock_size_log, unsigned k_cblock_size_log, unsigned k_sigma_log>
  const unsigned rank4n<k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_sigma_mask = (1 << k_sigma_log) - 1;

template<unsigned k_sblock_size_log, unsigned k_cblock_size_log, unsigned k_sigma_log>
  const unsigned long rank4n<k_sblock_size_log, k_cblock_size_log, k_sigma_log>
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
