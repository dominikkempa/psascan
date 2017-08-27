#ifndef __RANK_H_INCLUDED
#define __RANK_H_INCLUDED

#include <algorithm>
#include <vector>
#include <thread>
#include <omp.h>

#include "utils.hpp"

using namespace psascan_private;


template<unsigned k_sblock_size_log = 24, unsigned k_cblock_size_log = 20, unsigned k_sigma_log = 8>
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
    rank4n(const std::uint8_t *text, std::uint64_t length) {

      std::uint64_t max_threads = omp_get_max_threads();
      m_length = length;
      n_cblocks = (m_length + k_cblock_size - 1) / k_cblock_size;
      n_sblocks = (n_cblocks + k_cblocks_in_sblock - 1) / k_cblocks_in_sblock;

      m_count = (unsigned long *)malloc(256L * sizeof(unsigned long));
      std::fill(m_count, m_count + 256, 0UL);
      if (!m_length) return;


      //------------------------------------------------------------------------
      // STEP 1: split all cblocks into (roughly) equal size ranges. Each
      //         range is processed by one thread.
      //
      // During step 1 we compute:
      //  * type of each cblock (but don't yet encode it in bitvector)
      //  * depending on the type we compute:
      //    - of the cblock is encoded using lookup-table-based encoding, then
      //      we compute the whole trunk
      //    - in other case we compute and store: symbol mapping, symbol type
      //      (frequent / rare / non-occurring) and values of freq_cnt_log
      //      and rare_cnt_log.
      //------------------------------------------------------------------------
      unsigned long range_size = (n_cblocks + max_threads - 1) / max_threads;
      unsigned long n_ranges = (n_cblocks + range_size - 1) / range_size;

      m_sblock_header = (unsigned long *)malloc(n_sblocks * sizeof(unsigned long) * k_sigma);
      m_cblock_header = (unsigned long *)malloc(n_cblocks * sizeof(unsigned long));
      m_cblock_header2 = (unsigned long *)malloc(n_cblocks * k_sigma * sizeof(unsigned long));
      m_cblock_mapping = (unsigned char *)malloc(n_cblocks * k_sigma * 2);
      m_cblock_type = (unsigned char *)malloc((n_cblocks + 7) / 8);
      m_freq_trunk = (unsigned *)calloc(n_cblocks * k_cblock_size, sizeof(unsigned));
      std::fill(m_cblock_type, m_cblock_type + (n_cblocks + 7) / 8, 0);

      unsigned long *rare_trunk_size = new unsigned long[n_cblocks];
      std::fill(rare_trunk_size, rare_trunk_size + n_cblocks, 0);

      bool *cblock_type = new bool[n_cblocks];
      std::fill(cblock_type, cblock_type + n_cblocks, 0);

      unsigned **occ = (unsigned **)malloc(n_ranges * sizeof(unsigned *));
      for (unsigned long i = 0; i < n_ranges; ++i)
        occ[i] = (unsigned *)malloc((k_cblock_size + 1) * sizeof(unsigned));

      std::thread **threads = new std::thread*[n_ranges];
      for (unsigned long i = 0; i < n_ranges; ++i) {
        unsigned long range_beg = i * range_size;
        unsigned long range_end = std::min(range_beg + range_size, n_cblocks);

        threads[i] = new std::thread(construction_step_1, std::ref(*this),
            text, range_beg, range_end, rare_trunk_size, cblock_type, occ[i]);
      }
      for (unsigned long i = 0; i < n_ranges; ++i) threads[i]->join();
      for (unsigned long i = 0; i < n_ranges; ++i) delete threads[i];
      delete[] threads;

      for (unsigned long i = 0; i < n_ranges; ++i) {
        free(occ[i]);
      }
      free(occ);


      //------------------------------------------------------------------------
      // STEP 2: compute global information based on local cblock computation:
      //  * cblock types (we do it now, since two local cblock computations
      //    might try to update bits in the same byte
      //  * total size of rare trunk
      //  * pointers to the beginning of each rare trunk
      //  * cumulative counts of all symbols (necessary during queries)
      //  * non-inclusive partial sum over cblock range counts
      //------------------------------------------------------------------------
      unsigned long rare_trunk_total_size = 0;
      for (unsigned long cblock_id = 0; cblock_id < n_cblocks; ++cblock_id) {
        unsigned long cblock_beg = (cblock_id << k_cblock_size_log);

        // 1
        // Store cblock type
        if (cblock_type[cblock_id])
          m_cblock_type[cblock_id >> 3] |= (1 << (cblock_id & 7));

        // 2
        // Compute the pointer to rare trunk and update total rare trunk size.
        unsigned long this_cblock_rare_trunk_size = rare_trunk_size[cblock_id];
        m_cblock_header[cblock_id] |= (rare_trunk_total_size << 16);
        rare_trunk_total_size += this_cblock_rare_trunk_size;

        // 3
        // Update cblock header
        unsigned long cblock_header_beg = (cblock_id << k_sigma_log);
        for (unsigned c = 0; c < k_sigma; ++c)
          m_cblock_header2[cblock_header_beg + c] |= (m_count[c] << (k_cblock_size_log + 6));

        // 4
        // Update sblock header
        if (!(cblock_beg & k_sblock_size_mask)) {
          unsigned long sblock_id = (cblock_beg >> k_sblock_size_log);
          unsigned long sblock_header_beg = (sblock_id << k_sigma_log);
          for (unsigned c = 0; c < k_sigma; ++c)
            m_sblock_header[sblock_header_beg + c] = m_count[c];
        }

        // 5
        // Update m_count.
        unsigned long ptr = (cblock_id << k_sigma_log);
        for (unsigned c = 0; c + 1 < k_sigma; ++c)
          m_count[c] += ((m_cblock_header2[ptr + c + 1] >> 5) & k_2cblock_size_mask) - 
                        ((m_cblock_header2[ptr + c]     >> 5) & k_2cblock_size_mask);
        m_count[k_sigma - 1] += k_cblock_size - ((m_cblock_header2[ptr + k_sigma - 1] >> 5) & k_2cblock_size_mask);
      }
      delete[] cblock_type;
      delete[] rare_trunk_size;


      //------------------------------------------------------------------------
      // STEP 3: compute freq and rare trunk for cblocks encoded with
      //         alphabet partitioning.
      //------------------------------------------------------------------------
      m_rare_trunk = (unsigned *)calloc(rare_trunk_total_size, sizeof(unsigned));
      threads = new std::thread*[n_ranges];
      for (unsigned long i = 0; i < n_ranges; ++i) {
        unsigned long range_beg = i * range_size;
        unsigned long range_end = std::min(range_beg + range_size, n_cblocks);

        threads[i] = new std::thread(construction_step_3, std::ref(*this),
            text, range_beg, range_end);
      }
      for (unsigned long i = 0; i < n_ranges; ++i) threads[i]->join();
      for (unsigned long i = 0; i < n_ranges; ++i) delete threads[i];
      delete[] threads;

      m_count[0] -= n_cblocks * k_cblock_size - m_length;  // remove extra zeros
    }

    static void construction_step_1(rank4n &r, const unsigned char *text,
        unsigned long cblock_range_beg, unsigned long cblock_range_end,
        unsigned long *rare_trunk_size, bool *cblock_type, unsigned *occ) {
      std::vector<std::pair<uint32_t, unsigned char> > sorted_chars;
      std::vector<unsigned char> freq_chars;
      std::vector<unsigned char> rare_chars;

      unsigned *refpoint_precomputed = (unsigned *)malloc(k_cblock_size * sizeof(unsigned));

      unsigned *cblock_count = new unsigned[k_sigma];
      unsigned *list_beg = new unsigned[k_sigma];
      unsigned *list_beg2 = new unsigned[k_sigma];
      bool *isfreq = new bool[k_sigma];
      unsigned *lookup_bits_precomputed = new unsigned[k_sigma];
      unsigned *min_block_size_precomputed = new unsigned[k_sigma];
      unsigned long *refpoint_mask_precomputed = new unsigned long[k_sigma];

      for (unsigned long cblock_id = cblock_range_beg; cblock_id < cblock_range_end; ++cblock_id) {
        unsigned long cblock_beg = cblock_id << k_cblock_size_log;
        unsigned long cblock_end = cblock_beg + k_cblock_size;

        std::fill(cblock_count, cblock_count + k_sigma, 0);

        // compute cblock symbols into cblock_symbols tab.
        unsigned long maxj = std::min(cblock_end, r.m_length);
        for (unsigned long j = cblock_beg; j < maxj; ++j)
          ++cblock_count[text[j]];
        cblock_count[0] += cblock_end - maxj;

        // Compute starting positions of occurrences lists.
        for (unsigned j = 0, t, s = 0; j < k_sigma; ++j) {
          t = cblock_count[j];
          list_beg[j] = s;
          list_beg2[j] = s;
          s += t;
        }
        
        // Compute first part the cblock header: ranks up to cblock
        // beginning pointers to beginnings of occurrence lists.
        // Note: this implicitly encodes cblock counts.
        for (unsigned c = 0; c < k_sigma; ++c)
          r.m_cblock_header2[(cblock_id << k_sigma_log) + c] = (list_beg[c] << 5);

        // Sort symbol counts by frequencies.
        sorted_chars.clear();
        for (unsigned j = 0; j < k_sigma; ++j)
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
        r.m_cblock_header[cblock_id] = freq_cnt_log;
        r.m_cblock_header[cblock_id] |= (rare_cnt_log << 8);

        // Compute and store symbols mapping.
        std::sort(freq_chars.begin(), freq_chars.end());
        std::sort(rare_chars.begin(), rare_chars.end());
        std::fill(isfreq, isfreq + 256, false);
        for (unsigned c = 0; c < 256; ++c)
          r.m_cblock_mapping[2 * (c * r.n_cblocks + cblock_id)] = k_char_type_missing;
        for (unsigned i = 0; i < freq_chars.size(); ++i) {
          unsigned char c = freq_chars[i];
          isfreq[c] = true;
          r.m_cblock_mapping[2 * (c * r.n_cblocks + cblock_id) + 1] = i;
          r.m_cblock_mapping[2 * (c * r.n_cblocks + cblock_id)] = k_char_type_freq;
        }
        for (unsigned i = 0; i < rare_chars.size(); ++i) {
          unsigned char c = rare_chars[i];
          r.m_cblock_mapping[2 * (c * r.n_cblocks + cblock_id) + 1] = i;
          r.m_cblock_mapping[2 * (c * r.n_cblocks + cblock_id)] = k_char_type_rare;
        }

        unsigned nofreq_cnt = 0L;
        for (unsigned c = 0; c < k_sigma; ++c)
          if (!isfreq[c]) nofreq_cnt += cblock_count[c];


        if (freq_cnt >= 128) {
          // Set that the cblock is encoded using method for random strings.
          cblock_type[cblock_id] = true;
 
          // Compute lists of occurrences.
          for (unsigned long i = cblock_beg; i < maxj; ++i)
            occ[list_beg2[text[i]]++] = i - cblock_beg;
          for (unsigned long i = maxj; i < cblock_end; ++i)
            occ[list_beg2[0]++] = i - cblock_beg;

          // Precompute some values and store lookup bits into the header.
          for (unsigned c = 0; c < k_sigma; ++c) {
            lookup_bits_precomputed[c] = utils::log2ceil(cblock_count[c] + 2);
            r.m_cblock_header2[(cblock_id << 8) + c] |= lookup_bits_precomputed[c];
            if (cblock_count[c])
              min_block_size_precomputed[c] = k_cblock_size / cblock_count[c];
            else min_block_size_precomputed[c] = 0;

            unsigned refpoint_dist_log = 31 - lookup_bits_precomputed[c];
            unsigned long refpoint_dist = (1UL << refpoint_dist_log);
            unsigned long refpoint_dist_mask = refpoint_dist - 1;
            unsigned long refpoint_dist_mask_neg = (~refpoint_dist_mask);
            refpoint_mask_precomputed[c] = refpoint_dist_mask_neg;
          }

          unsigned *cblock_trunk = r.m_freq_trunk + cblock_beg;
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
          // Update rare_trunk_size.
          if (rare_cnt) {
            long rare_blocks = 1 + (nofreq_cnt + rare_cnt - 1) / rare_cnt;
            rare_trunk_size[cblock_id] = rare_blocks * rare_cnt;
          }
        }
      }
      delete[] list_beg;
      delete[] list_beg2;
      delete[] isfreq;
      delete[] cblock_count;
      delete[] lookup_bits_precomputed;
      delete[] min_block_size_precomputed;
      delete[] refpoint_mask_precomputed;
      free(refpoint_precomputed);
    }

    static void construction_step_3(rank4n &r, const unsigned char *text,
        unsigned long cblock_range_beg, unsigned long cblock_range_end) {
      unsigned char *freq_map = new unsigned char[k_sigma];
      unsigned char *rare_map = new unsigned char[k_sigma];
      unsigned long *cur_count = new unsigned long[k_sigma];
      unsigned long *off = new unsigned long[k_sigma];

      long *sblock_h = new long[k_sigma];
      int *israre = new int[k_sigma];

      std::vector<unsigned char> freq_chars;
      std::vector<unsigned char> rare_chars;

      for (unsigned long cblock_id = cblock_range_beg; cblock_id < cblock_range_end; ++cblock_id) {
        unsigned long cblock_beg = cblock_id << k_cblock_size_log;
        unsigned long cblock_end = cblock_beg + k_cblock_size;

        //----------------------------------------------------------------------
        // Process cblock [cblock_beg..cblock_end).
        //----------------------------------------------------------------------

        // Skip this cblock if its trunk was already encoded.
        if (r.m_cblock_type[cblock_id >> 3] & (1 << (cblock_id & 7))) continue;

        //----------------------------------
        // PHASE 1
        //----------------------------------
        // Retreive symbol counts up to this cblock begin and pointer
        // to rare trunk size from cblock headers.
        for (unsigned c = 0; c < k_sigma; ++c)
          cur_count[c] = (r.m_cblock_header2[(cblock_id << 8) + c] >> (k_cblock_size_log + 6));

        long r_filled  = (r.m_cblock_header[cblock_id] >> 16);
        long r_ptr = r_filled;

        long freq_cnt_log = (r.m_cblock_header[cblock_id] & 255L);
        long rare_cnt_log = ((r.m_cblock_header[cblock_id] >> 8) & 255L);
        long freq_cnt = (1L << freq_cnt_log);
        long rare_cnt = (1L << rare_cnt_log);
        long rare_cnt_mask = rare_cnt - 1;

        freq_chars.clear();
        rare_chars.clear();
        std::fill(israre, israre + k_sigma, 1);
        for (unsigned c = 0; c < k_sigma; ++c) {
          unsigned char type = r.m_cblock_mapping[2 * (c * r.n_cblocks + cblock_id)];
          if (type == k_char_type_freq) {
            israre[c] = 0;
            freq_chars.push_back(c);
            freq_map[c] = r.m_cblock_mapping[2 * (c * r.n_cblocks + cblock_id) + 1];
          } else if (type == k_char_type_rare) {
            rare_chars.push_back(c);
            rare_map[c] = r.m_cblock_mapping[2 * (c * r.n_cblocks + cblock_id) + 1];
            freq_map[c] = freq_cnt - 1;
          }
        }

        if (rare_chars.empty()) {
          rare_cnt_log = 0;
          rare_cnt = 0;
        }



        //--------------------
        // PHASE 2
        //--------------------
        long sblock_id = (cblock_beg >> k_sblock_size_log);
        std::copy(r.m_sblock_header + (sblock_id << 8), r.m_sblock_header + (sblock_id << 8) + k_sigma, sblock_h);
        for (long j = 0; j < k_sigma; ++j) off[j] = cur_count[j] - sblock_h[j];

        long nofreq_cnt = 0;
        long freq_chars_size = (long)freq_chars.size();
        long rare_chars_size = (long)rare_chars.size();

        if (cblock_end <= r.m_length) {
          for (unsigned long i = cblock_beg; i < cblock_end; i += freq_cnt) {
            for (long j = 0; j < freq_chars_size; ++j) {
              unsigned char ch = freq_chars[j];
              r.m_freq_trunk[i + j] = (off[ch] << 8);
            }
            r.m_freq_trunk[i + freq_cnt - 1] = (nofreq_cnt << 8);
            for (unsigned long j = i; j < i + freq_cnt; ++j) {
              unsigned char c = text[j];
              r.m_freq_trunk[j] |= freq_map[c];
              if (israre[c]) {
                if (!(nofreq_cnt & rare_cnt_mask)) {
                  for (long jj = 0; jj < rare_chars_size; ++jj) {
                    unsigned char ch = rare_chars[jj];
                    r.m_rare_trunk[r_filled++] = (off[ch] << 8);
                  }
                  r_filled += rare_cnt - rare_chars_size;
                }
                r.m_rare_trunk[r_ptr++] |= rare_map[c];
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
              r.m_freq_trunk[i + j] = (off[ch] << 8);
            }
            r.m_freq_trunk[i + freq_cnt - 1] = (nofreq_cnt << 8);
            for (unsigned long j = i; j < i + freq_cnt; ++j) {
              unsigned char c = (j < r.m_length ? text[j] : 0);
              r.m_freq_trunk[j] |= freq_map[c];
              if (israre[c]) {
                if (!(nofreq_cnt & rare_cnt_mask)) {
                  for (long jj = 0; jj < rare_chars_size; ++jj) {
                    unsigned char ch = rare_chars[jj];
                    r.m_rare_trunk[r_filled++] = (off[ch] << 8);
                  }
                  r_filled += rare_cnt - rare_chars_size;
                }
                r.m_rare_trunk[r_ptr++] |= rare_map[c];
              }
              ++off[c];
              nofreq_cnt += israre[c];
            }
          }
          for (long i = 0; i < k_sigma; ++i)
            cur_count[i] = sblock_h[i] + off[i];
        }


        //---------------------------
        // PHASE 3
        //---------------------------
        for (long j = 0; j < rare_cnt; ++j) {
          unsigned char ch = (j < (long)rare_chars.size() ? rare_chars[j] : 0);
          long local_rank = cur_count[ch] - r.m_sblock_header[(sblock_id << 8) + ch];
          r.m_rare_trunk[r_filled++] = (local_rank << 8);
        }
      }

      delete[] cur_count;
      delete[] sblock_h;
      delete[] freq_map;
      delete[] rare_map;
      delete[] israre;
      delete[] off;
    }

    inline long query(long i, unsigned char c) const {
      if (i <= 0) return 0L;
      else if ((unsigned long)i >= m_length) return m_count[c];

      unsigned long cblock_id = (i >> k_cblock_size_log);    
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

        // Empty block optimization.
        if (begin == list_size + 1) { // block containing cblock_i is empty, just find the beginning.
          ++approx;
          while ((m_freq_trunk[cblock_beg + list_beg + approx] & lookup_mask) == list_size + 1) ++approx;
          begin = (m_freq_trunk[cblock_beg + list_beg + approx] & lookup_mask);
          return rank_up_to_cblock + begin;
        }
        
        long next_block_begin =  (approx + 1 == list_size) ? list_size :
          (m_freq_trunk[cblock_beg + list_beg + approx + 1] & lookup_mask);

        // correct next_block_begin
        if (approx + 1 != list_size && next_block_begin == list_size + 1) {
          ++approx;
          while ((m_freq_trunk[cblock_beg + list_beg + approx + 1] & lookup_mask) == list_size + 1) ++approx;
          next_block_begin = (m_freq_trunk[cblock_beg + list_beg + approx + 1] & lookup_mask);
        }

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
            return query(cblock_id << k_cblock_size_log, c);
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
