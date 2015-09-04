/**
 * @file    src/psascan_src/inmem_psascan_src/rank.h
 * @author  Juha Karkkainen <juha.karkkainen (at) cs.helsinki.fi>
 *          Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * @section DESCRIPTION
 *
 * A general rank data structure. Basic idea of the encoding is from
 * the rank data structure used in the external-memory algorithm for
 * constructing the Burrows-Wheeler transform called bwtdisk (available
 * at: http://people.unipmn.it/manzini/bwtdisk/) described in [1]. We
 * extended the data structure by applying the fixed block boosting [2]
 * and alphabet partitioning [3] techniques. The resulting data structure
 * was described in [4]. This file extends the implementation used in [4]
 * by parallelizing the construction and introducting an alternative
 * encoding (called type-I in the code). Type-I encoding is a novel
 * encoding due to present authors.
 *
 * References:
 * [1] Paolo Ferragina, Travis Gagie, Giovanni Manzini:
 *     Lightweight Data Indexing and Compression in External Memory.
 *     Algorithmica 63(3), p. 707-730 (2012).
 * [2] Juha Karkkainen, Simon J. Puglisi:
 *     Fixed Block Compression Boosting in FM-Indexes.
 *     In Proc. SPIRE 2011, p. 174-184.
 * [3] Jeremy Barbay, Travis Gagie, Gonzalo Navarro, Yakov Nekrich:
 *     Alphabet Partitioning for Compressed Rank/Select and Applications.
 *     In Proc. ISAAC 2010, p. 315-326.
 * [4] Juha Karkkainen, Dominik Kempa:
 *     Engineering a Lightweight External Memory Suffix Array Construction
 *     Algorithm.
 *     In Proc. ICABD 2014, p. 53-60.
 *
 * @section LICENCE
 *
 * This file is part of pSAscan v0.1.0
 * See: http://www.cs.helsinki.fi/group/pads/
 *
 * Copyright (C) 2014-2015
 *   Juha Karkkainen <juha.karkkainen (at) cs.helsinki.fi>
 *   Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/

#ifndef __PSASCAN_SRC_INMEM_PSASCAN_SRC_RANK_H_INCLUDED
#define __PSASCAN_SRC_INMEM_PSASCAN_SRC_RANK_H_INCLUDED

#include <cstdio>
#include <algorithm>
#include <vector>
#include <thread>

#include "../utils.h"
#include "bwtsa.h"
#include "pagearray.h"


namespace psascan_private {
namespace inmem_psascan_private {

template<
  typename saidx_t,
  std::uint32_t pagesize_log,
  std::uint32_t k_sblock_size_log = 24,
  std::uint32_t k_cblock_size_log = 20,
  std::uint32_t k_sigma_log = 8>
class rank4n {
  private:
    typedef pagearray<bwtsa_t<saidx_t>, pagesize_log> pagearray_type;

    static const std::uint64_t k_cblock_size;
    static const std::uint64_t k_cblock_size_mask;
    static const std::uint64_t k_cblock_size_mask_neg;
    static const std::uint32_t k_cblocks_in_sblock_log;
    static const std::uint32_t k_cblocks_in_sblock;
    static const std::uint32_t k_cblocks_in_sblock_mask;
    static const std::uint32_t k_2cblock_size;
    static const std::uint32_t k_2cblock_size_mask;
    static const std::uint32_t k_sblock_size;
    static const std::uint32_t k_sblock_size_mask;
    static const std::uint32_t k_sigma;
    static const std::uint32_t k_sigma_mask;

    static const std::uint32_t pagesize = (1U << pagesize_log);
    static const std::uint32_t pagesize_mask = (1U << pagesize_log) - 1;

    static const std::uint32_t k_char_type_freq =    0x01;
    static const std::uint32_t k_char_type_rare =    0x02;
    static const std::uint32_t k_char_type_missing = 0x03;

    std::uint64_t m_length;   // length of original sequence
    std::uint64_t n_cblocks;  // number of context blocks
    std::uint64_t n_sblocks;  // number of super blocks

    std::uint64_t *m_sblock_header;
    std::uint64_t *m_cblock_header;
    std::uint64_t *m_cblock_header2;

    std::uint8_t *m_cblock_type;
    std::uint8_t *m_cblock_mapping;

    std::uint32_t *m_freq_trunk;
    std::uint32_t *m_rare_trunk;

  public:
    std::uint64_t *m_count;  // symbol counts

  public:
    rank4n(const pagearray_type *ptext, std::uint64_t length, std::uint64_t max_threads) {
      m_length = length;
      n_cblocks = (m_length + k_cblock_size - 1) / k_cblock_size;
      n_sblocks = (n_cblocks + k_cblocks_in_sblock - 1) / k_cblocks_in_sblock;

      m_count = (std::uint64_t *)malloc(256L * sizeof(std::uint64_t));
      std::fill(m_count, m_count + 256, 0UL);
      if (!m_length) return;

      long double start = utils::wclock();
      m_sblock_header = (std::uint64_t *)malloc(n_sblocks * sizeof(std::uint64_t) * k_sigma);
      m_cblock_header = (std::uint64_t *)malloc(n_cblocks * sizeof(std::uint64_t));
      m_cblock_header2 = (std::uint64_t *)malloc(n_cblocks * k_sigma * sizeof(std::uint64_t));
      m_cblock_mapping = (std::uint8_t *)malloc(n_cblocks * k_sigma * 2);
      m_cblock_type = (std::uint8_t *)malloc((n_cblocks + 7) / 8);
      m_freq_trunk = (std::uint32_t *)calloc(n_cblocks * k_cblock_size, sizeof(std::uint32_t));
      std::fill(m_cblock_type, m_cblock_type + (n_cblocks + 7) / 8, 0);
      std::uint8_t *bwt = (std::uint8_t *)malloc(length + k_cblock_size);
      long double alloc_time = utils::wclock() - start;
      if (alloc_time > 0.05L)
        fprintf(stderr, "alloc: %.2Lfs ", alloc_time);

      encode_type_I(ptext, bwt, max_threads);
      encode_type_II(bwt, max_threads);

      m_count[0] -= n_cblocks * k_cblock_size - m_length;  // remove extra zeros
      free(bwt);
    }

    void encode_type_I(const pagearray_type *ptext, std::uint8_t *bwt,
        std::uint64_t max_threads) {
      //------------------------------------------------------------------------
      // STEP 1: split all cblocks into equal size ranges (except possible the
      //         last one). Each range is processed by one thread. During this
      //         step we compute: (i) type of each cblock, (ii) encode all
      //         type-I cblocks and for all type-II cblocks, we compute and
      //         store: symbol mapping, symbol type (freq / rare / non-occurring)
      //         and values of freq_cnt_log and rare_cnt_log.
      //------------------------------------------------------------------------
      std::uint64_t range_size = (n_cblocks + max_threads - 1) / max_threads;
      std::uint64_t n_ranges = (n_cblocks + range_size - 1) / range_size;

      std::uint64_t *rare_trunk_size = new std::uint64_t[n_cblocks];
      std::fill(rare_trunk_size, rare_trunk_size + n_cblocks, 0);

      bool *cblock_type = new bool[n_cblocks];
      std::fill(cblock_type, cblock_type + n_cblocks, 0);

      std::uint32_t **occ = (std::uint32_t **)malloc(n_ranges * sizeof(std::uint32_t *));
      for (std::uint64_t i = 0; i < n_ranges; ++i)
        occ[i] = (std::uint32_t *)malloc((k_cblock_size + 1) * sizeof(std::uint32_t));

      fprintf(stderr, "s1: ");
      long double start = utils::wclock();
      std::thread **threads = new std::thread*[n_ranges];
      for (std::uint64_t i = 0; i < n_ranges; ++i) {
        std::uint64_t range_beg = i * range_size;
        std::uint64_t range_end = std::min(range_beg + range_size, n_cblocks);

        threads[i] = new std::thread(encode_type_I_aux, std::ref(*this),
            ptext, range_beg, range_end, rare_trunk_size, cblock_type, occ[i], bwt);
      }

      for (std::uint64_t i = 0; i < n_ranges; ++i) threads[i]->join();
      for (std::uint64_t i = 0; i < n_ranges; ++i) delete threads[i];
      delete[] threads;

      for (std::uint64_t i = 0; i < n_ranges; ++i)
        free(occ[i]);
      free(occ);

      fprintf(stderr, "%.2Lfs ", utils::wclock() - start);


      //------------------------------------------------------------------------
      // STEP 2: compute global information based on local cblock computation:
      //   * store cblock types,
      //   * total size of rare trunk,
      //   * pointers to the beginning of each rare trunk,
      //   * cumulative counts of all symbols,
      //   * non-inclusive partial sum over cblock range counts.
      //------------------------------------------------------------------------
      fprintf(stderr, "s2: ");
      start = utils::wclock();
      std::uint64_t rare_trunk_total_size = 0;
      for (std::uint64_t cblock_id = 0; cblock_id < n_cblocks; ++cblock_id) {
        std::uint64_t cblock_beg = (cblock_id << k_cblock_size_log);

        // 1
        // Store cblock type.
        if (cblock_type[cblock_id])
          m_cblock_type[cblock_id >> 3] |= (1 << (cblock_id & 7));

        // 2
        // Compute the pointer to rare trunk and update total rare trunk size.
        std::uint64_t this_cblock_rare_trunk_size = rare_trunk_size[cblock_id];
        m_cblock_header[cblock_id] |= (rare_trunk_total_size << 16);
        rare_trunk_total_size += this_cblock_rare_trunk_size;

        // 3
        // Update cblock header.
        std::uint64_t cblock_header_beg = (cblock_id << k_sigma_log);
        for (std::uint32_t c = 0; c < k_sigma; ++c)
          m_cblock_header2[cblock_header_beg + c] |= (m_count[c] << (k_cblock_size_log + 6));

        // 4
        // Update sblock header,
        if (!(cblock_beg & k_sblock_size_mask)) {
          std::uint64_t sblock_id = (cblock_beg >> k_sblock_size_log);
          std::uint64_t sblock_header_beg = (sblock_id << k_sigma_log);
          for (std::uint32_t c = 0; c < k_sigma; ++c)
            m_sblock_header[sblock_header_beg + c] = m_count[c];
        }

        // 5
        // Update m_count.
        std::uint64_t ptr = (cblock_id << k_sigma_log);
        for (std::uint32_t c = 0; c + 1 < k_sigma; ++c)
          m_count[c] += ((m_cblock_header2[ptr + c + 1] >> 5) & k_2cblock_size_mask) -
            ((m_cblock_header2[ptr + c]     >> 5) & k_2cblock_size_mask);
        m_count[k_sigma - 1] += k_cblock_size -
          ((m_cblock_header2[ptr + k_sigma - 1] >> 5) & k_2cblock_size_mask);
      }
      m_rare_trunk = (std::uint32_t *)calloc(rare_trunk_total_size, sizeof(std::uint32_t));

      delete[] cblock_type;
      delete[] rare_trunk_size;

      fprintf(stderr, "%.2Lfs ", utils::wclock() - start);
    }

    static void encode_type_I_aux(rank4n &r, const pagearray_type *ptext,
        std::uint64_t cblock_range_beg, std::uint64_t cblock_range_end,
        std::uint64_t *rare_trunk_size, bool *cblock_type, std::uint32_t *occ, std::uint8_t *bwt) {
      std::vector<std::pair<uint32_t, std::uint8_t> > sorted_chars;
      std::vector<std::uint8_t> freq_chars;
      std::vector<std::uint8_t> rare_chars;

      std::uint32_t *refpoint_precomputed = (std::uint32_t *)malloc(k_cblock_size * sizeof(std::uint32_t));
      std::uint32_t *cblock_count = new std::uint32_t[k_sigma];
      std::uint32_t *list_beg = new std::uint32_t[k_sigma];
      std::uint32_t *list_beg2 = new std::uint32_t[k_sigma];
      bool *isfreq = new bool[k_sigma];
      std::uint32_t *lookup_bits_precomputed = new std::uint32_t[k_sigma];
      std::uint32_t *min_block_size_precomputed = new std::uint32_t[k_sigma];
      std::uint64_t *refpoint_mask_precomputed = new std::uint64_t[k_sigma];

      typedef typename pagearray_type::value_type value_type;

      // Process cblocks one by one.
      for (std::uint64_t cblock_id = cblock_range_beg; cblock_id < cblock_range_end; ++cblock_id) {
        std::uint64_t cblock_beg = cblock_id << k_cblock_size_log;
        std::uint64_t cblock_end = cblock_beg + k_cblock_size;

        // Compute symbol counts inside cblock and store bwt symbols.
        std::fill(cblock_count, cblock_count + k_sigma, 0);
        std::uint64_t maxj = std::min(cblock_end, r.m_length);
        std::uint64_t page_id = (cblock_beg >> pagesize_log);
        value_type *cur_page = ptext->m_pageindex[page_id++];
        std::uint64_t page_offset = ptext->get_page_offset(cblock_beg);
        for (std::uint64_t j = cblock_beg; j < maxj; ++j) {
          std::uint8_t c = cur_page[page_offset].m_bwt;
          bwt[j] = c;
          ++cblock_count[c];
          ++page_offset;
          if (page_offset == pagesize) {
            cur_page = ptext->m_pageindex[page_id];
            ++page_id;
            page_offset = 0;
          }
        }
        for (std::uint64_t j = maxj; j < cblock_end; ++j) {
          bwt[j] = 0;
          ++cblock_count[0];
        }


        // Compute starting positions of occurrences lists.
        for (std::uint32_t j = 0, t, s = 0; j < k_sigma; ++j) {
          t = cblock_count[j];
          list_beg[j] = s;
          list_beg2[j] = s;
          s += t;
        }
        
        // Store pointers to beginnings of occurrence lists in the type-I
        // cblock header. Note: this implicitly encodes cblock counts.
        for (std::uint32_t c = 0; c < k_sigma; ++c)
          r.m_cblock_header2[(cblock_id << k_sigma_log) + c] = (list_beg[c] << 5);

        // Sort symbol counts by frequencies.
        sorted_chars.clear();
        for (std::uint32_t j = 0; j < k_sigma; ++j)
          if (cblock_count[j])
            sorted_chars.push_back(std::make_pair(cblock_count[j], j));
        std::sort(sorted_chars.begin(), sorted_chars.end());

        // Separate (at most, due to rounding of freq_cnt)
        // about 3% of rarest symbols.
        std::uint32_t rare_cnt = 0L, rare_sum = 0L;
        while (rare_cnt < sorted_chars.size() &&
            16L * (rare_sum + sorted_chars[rare_cnt].first) <= k_cblock_size)
          rare_sum += sorted_chars[rare_cnt++].first;

        // Compute freq_cnt. Then round up freq_cnt + 1 (+1 is
        // for rare char marker) to the smallest power of two.
        // Note: rare_cnt > 0, so after rounding freq_cnt <= 256.
        std::uint32_t freq_cnt = sorted_chars.size() - rare_cnt;
        std::uint32_t freq_cnt_log = utils::log2ceil(freq_cnt + 1);
        freq_cnt = (1 << freq_cnt_log);

        // Recompute rare_cnt (note the +1).
        rare_cnt = 0;
        if (sorted_chars.size() + 1 > freq_cnt)
          rare_cnt = sorted_chars.size() + 1 - freq_cnt;

        // Compute freq and rare chars.
        rare_chars.clear();
        freq_chars.clear();
        for (std::uint32_t i = 0; i < rare_cnt; ++i)
          rare_chars.push_back(sorted_chars[i].second);
        for (std::uint32_t i = rare_cnt; i < sorted_chars.size(); ++i)
          freq_chars.push_back(sorted_chars[i].second);

        // If there are rare symbols, round up
        // rare_cnt to the smallest power of two.
        std::uint32_t rare_cnt_log = 0;
        if (rare_cnt) {
          rare_cnt_log = utils::log2ceil(rare_cnt);
          rare_cnt = (1 << rare_cnt_log);
        }

        // Update cblock type-I header.
        r.m_cblock_header[cblock_id] = freq_cnt_log;
        r.m_cblock_header[cblock_id] |= (rare_cnt_log << 8);

        // Compute and store symbols mapping.
        std::sort(freq_chars.begin(), freq_chars.end());
        std::sort(rare_chars.begin(), rare_chars.end());
        std::fill(isfreq, isfreq + 256, false);
        for (std::uint32_t c = 0; c < 256; ++c)
          r.m_cblock_mapping[2 * (c * r.n_cblocks + cblock_id)] = k_char_type_missing;
        for (std::uint32_t i = 0; i < freq_chars.size(); ++i) {
          std::uint8_t c = freq_chars[i];
          isfreq[c] = true;
          r.m_cblock_mapping[2 * (c * r.n_cblocks + cblock_id) + 1] = i;
          r.m_cblock_mapping[2 * (c * r.n_cblocks + cblock_id)] = k_char_type_freq;
        }
        for (std::uint32_t i = 0; i < rare_chars.size(); ++i) {
          std::uint8_t c = rare_chars[i];
          r.m_cblock_mapping[2 * (c * r.n_cblocks + cblock_id) + 1] = i;
          r.m_cblock_mapping[2 * (c * r.n_cblocks + cblock_id)] = k_char_type_rare;
        }

        std::uint32_t nofreq_cnt = 0L;
        for (std::uint32_t c = 0; c < k_sigma; ++c)
          if (!isfreq[c]) nofreq_cnt += cblock_count[c];


        if (freq_cnt >= 128) {  // type-I cblock
          cblock_type[cblock_id] = true;
 
          // Compute lists of occurrences.
          for (std::uint64_t i = cblock_beg; i < cblock_end; ++i)
            occ[list_beg2[bwt[i]]++] = i - cblock_beg;

          // Precompute helper arrays and and store lookup bits into the header.
          for (std::uint32_t c = 0; c < k_sigma; ++c) {
            lookup_bits_precomputed[c] = utils::log2ceil(cblock_count[c] + 2);
            r.m_cblock_header2[(cblock_id << 8) + c] |= lookup_bits_precomputed[c];
            if (cblock_count[c])
              min_block_size_precomputed[c] = k_cblock_size / cblock_count[c];
            else min_block_size_precomputed[c] = 0;

            std::uint32_t refpoint_dist_log = 31 - lookup_bits_precomputed[c];
            std::uint64_t refpoint_dist = (1UL << refpoint_dist_log);
            std::uint64_t refpoint_dist_mask = refpoint_dist - 1;
            std::uint64_t refpoint_dist_mask_neg = (~refpoint_dist_mask);
            refpoint_mask_precomputed[c] = refpoint_dist_mask_neg;
          }

          // Actual encoding follows.
          std::uint32_t *cblock_trunk = r.m_freq_trunk + cblock_beg;
          for (std::uint32_t c = 0; c < k_sigma; ++c) {
            std::uint32_t freq = cblock_count[c];
            std::uint32_t min_block_size = min_block_size_precomputed[c];
            std::uint32_t lookup_bits = lookup_bits_precomputed[c];
            std::uint32_t refpoint_dist_mask_neg = refpoint_mask_precomputed[c];
            std::uint32_t c_list_beg = list_beg[c];

            for (std::uint32_t j = 0; j < freq; ++j)
              cblock_trunk[c_list_beg + j] = freq + 1;
            if (freq) cblock_trunk[c_list_beg + freq - 1] = freq;

            std::uint32_t block_beg = 0;
            for (std::uint32_t j = 0; j < freq; ++j) {
              refpoint_precomputed[j] = (block_beg & refpoint_dist_mask_neg);
              block_beg += min_block_size;
              if ((((std::uint64_t)block_beg * freq) >> k_cblock_size_log) == j) ++block_beg;
            }

            std::uint32_t refpoint, block_id;
            std::uint32_t mask = (~((1UL << lookup_bits) - 1));
            if (freq) {
              for (std::uint64_t j = freq; j > 0; --j) {
                block_id = (((std::uint64_t)occ[c_list_beg + j - 1] * freq) >> k_cblock_size_log);
                refpoint = refpoint_precomputed[block_id];
                cblock_trunk[c_list_beg + block_id] &= mask;
                cblock_trunk[c_list_beg + block_id] |= (j - 1);
                cblock_trunk[c_list_beg + j - 1] |= ((occ[c_list_beg + j - 1] - refpoint) << lookup_bits);
              }
            }
          }
        } else {
          // Update rare_trunk_size.
          if (rare_cnt) {
            std::uint64_t rare_blocks = 1 + (nofreq_cnt + rare_cnt - 1) / rare_cnt;
            rare_trunk_size[cblock_id] = rare_blocks * rare_cnt;
          }
        }
      }

      // Clean up.
      delete[] list_beg;
      delete[] list_beg2;
      delete[] isfreq;
      delete[] cblock_count;
      delete[] lookup_bits_precomputed;
      delete[] min_block_size_precomputed;
      delete[] refpoint_mask_precomputed;
      free(refpoint_precomputed);
    }

    void encode_type_II(const std::uint8_t *bwt, std::uint64_t max_threads) {
      fprintf(stderr, "s3: ");
      long double start = utils::wclock();

      std::uint64_t range_size = (n_cblocks + max_threads - 1) / max_threads;
      std::uint64_t n_ranges = (n_cblocks + range_size - 1) / range_size;

      std::thread **threads = new std::thread*[n_ranges];
      for (std::uint64_t i = 0; i < n_ranges; ++i) {
        std::uint64_t range_beg = i * range_size;
        std::uint64_t range_end = std::min(range_beg + range_size, n_cblocks);

        threads[i] = new std::thread(encode_type_II_aux,
            std::ref(*this), range_beg, range_end, bwt);
      }

      for (std::uint64_t i = 0; i < n_ranges; ++i) threads[i]->join();
      for (std::uint64_t i = 0; i < n_ranges; ++i) delete threads[i];
      delete[] threads;

      fprintf(stderr, "%.2Lfs ", utils::wclock() - start);
    }

    static void encode_type_II_aux(rank4n &r, std::uint64_t cblock_range_beg,
        std::uint64_t cblock_range_end, const std::uint8_t *bwt) {
      std::uint8_t *freq_map = new std::uint8_t[k_sigma];
      std::uint8_t *rare_map = new std::uint8_t[k_sigma];
      std::uint64_t *cur_count = new std::uint64_t[k_sigma];
      std::uint64_t *off = new std::uint64_t[k_sigma];

      std::uint64_t *sblock_h = new std::uint64_t[k_sigma];
      bool *israre = new bool[k_sigma];

      std::vector<std::uint8_t> freq_chars;
      std::vector<std::uint8_t> rare_chars;

      for (std::uint64_t cblock_id = cblock_range_beg; cblock_id < cblock_range_end; ++cblock_id) {
        std::uint64_t cblock_beg = cblock_id << k_cblock_size_log;
        std::uint64_t cblock_end = cblock_beg + k_cblock_size;

        // Skip the cblock if it was type-I encoded.
        if (r.m_cblock_type[cblock_id >> 3] & (1 << (cblock_id & 7))) continue;

        // Retreive symbol counts up to this cblock begin and
        // pointer to rare trunk size from cblock headers.
        for (std::uint32_t c = 0; c < k_sigma; ++c)
          cur_count[c] = (r.m_cblock_header2[(cblock_id << 8) + c] >> (k_cblock_size_log + 6));

        std::uint64_t r_filled  = (r.m_cblock_header[cblock_id] >> 16);
        std::uint64_t r_ptr = r_filled;

        std::uint64_t freq_cnt_log = (r.m_cblock_header[cblock_id] & 255UL);
        std::uint64_t rare_cnt_log = ((r.m_cblock_header[cblock_id] >> 8) & 255UL);
        std::uint64_t freq_cnt = (1UL << freq_cnt_log);
        std::uint64_t rare_cnt = (1UL << rare_cnt_log);
        std::uint64_t rare_cnt_mask = rare_cnt - 1;

        freq_chars.clear();
        rare_chars.clear();
        std::fill(israre, israre + k_sigma, true);
        for (std::uint64_t c = 0; c < k_sigma; ++c) {
          std::uint8_t type = r.m_cblock_mapping[2 * (c * r.n_cblocks + cblock_id)];
          if (type == k_char_type_freq) {
            israre[c] = false;
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

        std::uint64_t sblock_id = (cblock_beg >> k_sblock_size_log);
        std::copy(r.m_sblock_header + (sblock_id << 8), r.m_sblock_header + (sblock_id << 8) + k_sigma, sblock_h);
        for (std::uint64_t j = 0; j < k_sigma; ++j) off[j] = cur_count[j] - sblock_h[j];

        std::uint64_t nofreq_cnt = 0;
        std::uint64_t freq_chars_size = freq_chars.size();
        std::uint64_t rare_chars_size = rare_chars.size();
        for (std::uint64_t i = cblock_beg; i < cblock_end; i += freq_cnt) {
          for (std::uint64_t j = 0; j < freq_chars_size; ++j) {
            std::uint8_t ch = freq_chars[j];
            r.m_freq_trunk[i + j] = (off[ch] << 8);
          }
          r.m_freq_trunk[i + freq_cnt - 1] = (nofreq_cnt << 8);
          for (std::uint64_t j = i; j < i + freq_cnt; ++j) {
            std::uint8_t c = bwt[j];
            r.m_freq_trunk[j] |= freq_map[c];
            if (israre[c]) {
              if (!(nofreq_cnt & rare_cnt_mask)) {
                for (std::uint64_t jj = 0; jj < rare_chars_size; ++jj) {
                  std::uint8_t ch = rare_chars[jj];
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
        for (std::uint32_t i = 0; i < k_sigma; ++i)
          cur_count[i] = sblock_h[i] + off[i];

        for (std::uint64_t j = 0; j < rare_cnt; ++j) {
          std::uint8_t ch = (j < rare_chars.size() ? rare_chars[j] : 0);
          std::uint64_t local_rank = cur_count[ch] - r.m_sblock_header[(sblock_id << 8) + ch];
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

    inline std::uint64_t rank(std::uint64_t i, std::uint8_t c) const {
      if (i == 0) return 0UL;
      else if (i >= m_length) return m_count[c];

      std::uint64_t cblock_id = (i >> k_cblock_size_log);    
      if (m_cblock_type[cblock_id >> 3] & (1 << (cblock_id & 7))) {  // type-I cblock
        std::uint64_t cblock_beg = (i & k_cblock_size_mask_neg);
        std::uint64_t cblock_i = (i & k_cblock_size_mask);     // offset in cblock
      
        // Extract the rank up to the start of cblock.
        std::uint64_t rank_up_to_cblock = (m_cblock_header2[(cblock_id << k_sigma_log) + c] >> (k_cblock_size_log + 6));

        // Now we compute the number of occurrences of c inside the cblock.
        // First, decode the beginning and end of c's occurrence list.
        std::uint64_t list_beg = ((m_cblock_header2[(cblock_id << k_sigma_log) + c] >> 5) & k_2cblock_size_mask);
        std::uint64_t list_end = ((c == k_sigma - 1) ? k_cblock_size :
            ((m_cblock_header2[(cblock_id << k_sigma_log) + c + 1] >> 5) & k_2cblock_size_mask));
        if (list_beg == list_end) return rank_up_to_cblock;

        // Compute the distance from i to the closest reference point on the left.
        std::uint64_t lookup_bits = (m_cblock_header2[(cblock_id << k_sigma_log) + c] & 31);
        std::uint64_t refpoint_dist_log = 31 - lookup_bits;
        std::uint64_t refpoint_disk_mask = (1UL << refpoint_dist_log) - 1;
        std::uint64_t i_refpoint_offset = (cblock_i & refpoint_disk_mask);

        // Compute threshold of symbol c inside the current cblock.
        std::uint64_t threshold = (1UL << (k_cblock_size_log - lookup_bits + 1));

        // Compute the id of block containing i.
        std::uint64_t list_size = list_end - list_beg;
        std::uint64_t approx = ((cblock_i * list_size) >> k_cblock_size_log);

        // Extract the lookup table entry.
        std::uint64_t lookup_mask = (1L << lookup_bits) - 1;
        std::uint64_t begin = (m_freq_trunk[cblock_beg + list_beg + approx] & lookup_mask);

        // Empty block optimization.
        if (begin == list_size + 1) {
          // Block containing cblock_i is empty, just find the beginning.
          ++approx;
          while ((m_freq_trunk[cblock_beg + list_beg + approx] & lookup_mask) == list_size + 1) ++approx;
          begin = (m_freq_trunk[cblock_beg + list_beg + approx] & lookup_mask);
          return rank_up_to_cblock + begin;
        }
        
        std::uint64_t next_block_begin =  (approx + 1 == list_size) ? list_size :
          (m_freq_trunk[cblock_beg + list_beg + approx + 1] & lookup_mask);

        // Correct next_block_begin.
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
            i_refpoint_offset += (1L << refpoint_dist_log);
            while (begin < next_block_begin && (m_freq_trunk[cblock_beg + list_beg + begin] >> lookup_bits) < i_refpoint_offset)
              ++begin;

            return rank_up_to_cblock + begin;
          }
        }
      } else {  // type-II cblock
        std::uint64_t sblock_id = (i >> k_sblock_size_log);
        std::uint64_t sblock_rank = m_sblock_header[(sblock_id << 8) + c];

        std::uint8_t type = m_cblock_mapping[2 * (c * n_cblocks + cblock_id)];
        std::uint8_t c_map = m_cblock_mapping[2 * (c * n_cblocks + cblock_id) + 1];

        std::uint64_t freq_cnt_bits = (m_cblock_header[cblock_id] & 255UL);
        std::uint64_t rare_cnt_bits = ((m_cblock_header[cblock_id] >> 8) & 255UL);
        std::uint64_t block_id = (i >> freq_cnt_bits);

        if (type == k_char_type_freq) {
          // Case 1 (fastest): symbol c was frequent in the context block.
          // Answer a query using frequent trunk.
          std::uint64_t block_rank = (m_freq_trunk[(block_id << freq_cnt_bits) + c_map] >> 8);
          std::uint64_t extra = 0;
          for (std::uint64_t j = (block_id << freq_cnt_bits); j < i; ++j)
            if ((m_freq_trunk[j] & 255) == c_map) ++extra;

          return sblock_rank + block_rank + extra;
        } else if (type == k_char_type_rare) {
          // Case 2: symbol c was rare inside the context block.
          // Compute new_i.
          std::uint64_t rare_trunk_ptr = (m_cblock_header[cblock_id] >> 16);
          std::uint64_t new_i = (m_freq_trunk[((block_id + 1) << freq_cnt_bits) - 1] >> 8);
          for (std::uint64_t j = (block_id << freq_cnt_bits); j < i; ++j)
            if ((m_freq_trunk[j] & 255) + 1 == (1U << freq_cnt_bits)) ++new_i;
      
          // Answer a query on rare trunk.
          std::uint64_t rare_block_id = (new_i >> rare_cnt_bits);
          std::uint64_t block_rank = (m_rare_trunk[rare_trunk_ptr +
            (rare_block_id << rare_cnt_bits) + c_map] >> 8);
          std::uint64_t extra = 0;
          for (std::uint64_t j = (rare_block_id << rare_cnt_bits); j < new_i; ++j)
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


template<typename saidx_t, std::uint32_t pagesize_log, std::uint32_t k_sblock_size_log, std::uint32_t k_cblock_size_log, std::uint32_t k_sigma_log>
  const std::uint64_t rank4n<saidx_t, pagesize_log, k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_cblock_size = (1UL << k_cblock_size_log);

template<typename saidx_t, std::uint32_t pagesize_log, std::uint32_t k_sblock_size_log, std::uint32_t k_cblock_size_log, std::uint32_t k_sigma_log>
  const std::uint64_t rank4n<saidx_t, pagesize_log, k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_cblock_size_mask = (1UL << k_cblock_size_log) - 1;
  
template<typename saidx_t, std::uint32_t pagesize_log, std::uint32_t k_sblock_size_log, std::uint32_t k_cblock_size_log, std::uint32_t k_sigma_log>
  const std::uint32_t rank4n<saidx_t, pagesize_log, k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_2cblock_size = (2U << k_cblock_size_log);

template<typename saidx_t, std::uint32_t pagesize_log, std::uint32_t k_sblock_size_log, std::uint32_t k_cblock_size_log, std::uint32_t k_sigma_log>
  const std::uint32_t rank4n<saidx_t, pagesize_log, k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_2cblock_size_mask = (2U << k_cblock_size_log) - 1;

template<typename saidx_t, std::uint32_t pagesize_log, std::uint32_t k_sblock_size_log, std::uint32_t k_cblock_size_log, std::uint32_t k_sigma_log>
  const std::uint32_t rank4n<saidx_t, pagesize_log, k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_sigma = (1U << k_sigma_log);

template<typename saidx_t, std::uint32_t pagesize_log, std::uint32_t k_sblock_size_log, std::uint32_t k_cblock_size_log, std::uint32_t k_sigma_log>
  const std::uint32_t rank4n<saidx_t, pagesize_log, k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_sigma_mask = (1U << k_sigma_log) - 1;

template<typename saidx_t, std::uint32_t pagesize_log, std::uint32_t k_sblock_size_log, std::uint32_t k_cblock_size_log, std::uint32_t k_sigma_log>
  const std::uint64_t rank4n<saidx_t, pagesize_log, k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_cblock_size_mask_neg = ~((1L << k_cblock_size_log) - 1);

template<typename saidx_t, std::uint32_t pagesize_log, std::uint32_t k_sblock_size_log, std::uint32_t k_cblock_size_log, std::uint32_t k_sigma_log>
  const std::uint32_t rank4n<saidx_t, pagesize_log, k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_cblocks_in_sblock_log = k_sblock_size_log - k_cblock_size_log;

template<typename saidx_t, std::uint32_t pagesize_log, std::uint32_t k_sblock_size_log, std::uint32_t k_cblock_size_log, std::uint32_t k_sigma_log>
  const std::uint32_t rank4n<saidx_t, pagesize_log, k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_cblocks_in_sblock = (1U << (k_sblock_size_log - k_cblock_size_log));

template<typename saidx_t, std::uint32_t pagesize_log, std::uint32_t k_sblock_size_log, std::uint32_t k_cblock_size_log, std::uint32_t k_sigma_log>
  const std::uint32_t rank4n<saidx_t, pagesize_log, k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_cblocks_in_sblock_mask = (1U << (k_sblock_size_log - k_cblock_size_log)) - 1;

template<typename saidx_t, std::uint32_t pagesize_log, std::uint32_t k_sblock_size_log, std::uint32_t k_cblock_size_log, std::uint32_t k_sigma_log>
  const std::uint32_t rank4n<saidx_t, pagesize_log, k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_sblock_size = (1U << k_sblock_size_log);
    
template<typename saidx_t, std::uint32_t pagesize_log, std::uint32_t k_sblock_size_log, std::uint32_t k_cblock_size_log, std::uint32_t k_sigma_log>
  const std::uint32_t rank4n<saidx_t, pagesize_log, k_sblock_size_log, k_cblock_size_log, k_sigma_log>
  ::k_sblock_size_mask = (1U << k_sblock_size_log) - 1;

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __PSASCAN_SRC_INMEM_PSASCAN_SRC_RANK_H_INCLUDED
