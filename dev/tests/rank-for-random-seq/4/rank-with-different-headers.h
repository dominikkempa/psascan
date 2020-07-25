#ifndef __RANK_H_INCLUDED
#define __RANK_H_INCLUDED

#include <algorithm>
#include <vector>

#include "utils.h"

template<long k_cblock_size_log = 18, long k_sigma_log = 8>
struct rank4n {
  rank4n(unsigned char *text, long length, long) {
    m_length = length;

    // Global symbol counts.
    m_count = new unsigned long[k_sigma];
    std::fill(m_count, m_count + k_sigma, 0L);

    // Allocate trunk and block header.
    n_cblocks = (m_length + k_cblock_size - 1) / k_cblock_size;
    m_trunk = (unsigned *)malloc(n_cblocks * k_cblock_size * sizeof(unsigned));
    std::fill(m_trunk, m_trunk + n_cblocks * k_cblock_size, 0U);
    m_cblock_header = new unsigned long[k_sigma * n_cblocks];
  
    // Process cblocks.
    long *list_beg = new long[k_sigma];      // beginnings of occurrence lists
    long *cblock_count = new long[k_sigma];  // symbol counts in the cblock

    std::vector<long> *occ = new std::vector<long>[k_sigma];
    for (long cblock_id = 0; cblock_id < n_cblocks; ++cblock_id) {
      long cblock_beg = k_cblock_size * cblock_id;
      long cblock_end = cblock_beg + k_cblock_size;

      //------------------------------------------------------------------------
      // STEP 1: compute symbol counts inside cblock.
      //------------------------------------------------------------------------
      std::fill(cblock_count, cblock_count + k_sigma, 0L);
      for (long i = cblock_beg; i < cblock_end; ++i) {
        unsigned char c = (i < m_length ? text[i] : 0);
        ++cblock_count[c];
      }

      //------------------------------------------------------------------------
      // STEP 2: compute starting positions of symbol occurrences lists.
      //------------------------------------------------------------------------
      for (long j = 0, t, s = 0; j < k_sigma; ++j) {
        t = cblock_count[j];
        list_beg[j] = s;
        s += t;
      }

      //------------------------------------------------------------------------
      // STEP 3: compute first part the cblock header: ranks up to cblock
      //         beginning pointers to beginnings of occurrence lists.
      //------------------------------------------------------------------------
      for (long c = 0; c < k_sigma; ++c) {
//        m_cblock_header[(cblock_id << k_sigma_log) + c] |= (list_beg[c] << 3);
//        m_cblock_header[(cblock_id << k_sigma_log) + c] |= (m_count[c] << (k_cblock_size_log + 4));
        long c2 = (c & (~1));
        long header_ptr = (c2 * n_cblocks) + (2 * cblock_id);
        if (c == c2) {
          m_cblock_header[header_ptr] |= m_count[c];
          m_cblock_header[header_ptr + 1] |= list_beg[c];
        } else {
          m_cblock_header[header_ptr] |= (m_count[c] << 32);
          m_cblock_header[header_ptr + 1] |= (list_beg[c] << (k_cblock_size_log + 1));
          m_cblock_header[header_ptr + 1] |= (cblock_count[c] << (2 * (k_cblock_size_log + 1)));
        }
      }

      //------------------------------------------------------------------------
      // STEP 4: update global symbol counts.
      //------------------------------------------------------------------------
      for (long i = cblock_beg; i < cblock_end; ++i) {
        unsigned char c = (i < m_length ? text[i] : 0);
        ++m_count[c];
      }

      for (long c = 0; c < k_sigma; ++c)
        occ[c].clear();

      //------------------------------------------------------------------------
      // STEP 5: compute lists of occurrences.
      //------------------------------------------------------------------------
      unsigned *cblock_trunk = m_trunk + cblock_beg;
      for (long i = cblock_beg; i < cblock_end; ++i) {
        unsigned char c = (i < m_length ? text[i] : 0);
        long cblock_i = i - cblock_beg;
        occ[c].push_back(cblock_i);
      }

      //------------------------------------------------------------------------
      // STEP 6: store list of occurrences and lookup tables.
      //------------------------------------------------------------------------
      for (long c = 0; c < k_sigma; ++c) {
        long freq = occ[c].size();

        // Compute the number of bits necessary to encode lookup table entries.
        // Lookup table can store values in the range [0..freq], thus we need
        // ceil(log2(freq + 1)) bits. This values is called lookup_bits and we
        // store it in the block header using 5 bits.

        // long lookup_bits = utils::log2ceil(freq + 1);

        // !!!
        long lookup_bits = 0;
        while ((1 << lookup_bits) < freq + 1) lookup_bits += 3;
        // !!!

        // m_cblock_header[(cblock_id << k_sigma_log) + c] |= (lookup_bits / 3);
        long c2 = (c & (~1));
        long header_ptr = (c2 * n_cblocks) + (2 * cblock_id);
        if (c == c2) m_cblock_header[header_ptr + 1] |= ((lookup_bits / 3) << (3 * (k_cblock_size_log + 1)));
        else m_cblock_header[header_ptr + 1] |= ((lookup_bits / 3) << (3 * (k_cblock_size_log + 2)));

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
    }

    m_count[0] -= n_cblocks * k_cblock_size - m_length; // remove extra zeros
    delete[] occ;
    delete[] cblock_count;
    delete[] list_beg;
  }
  
  inline long rank(long i, unsigned char c) {
    if (i <= 0) return 0L;
    else if (i >= m_length) return m_count[c];

    long cblock_id = (i >> k_cblock_size_log); // which cblock
    long cblock_beg = (i & k_cblock_size_mask_neg);
    long cblock_i = (i & k_cblock_size_mask);     // offset in cblock

    //--------------------------------------------------------------------------
    // STEP 1: extract the rank up to the start of cblock.
    //--------------------------------------------------------------------------
    unsigned char c2 = (c & (~1));
    long header_ptr = (c2 * n_cblocks) + (2 * cblock_id);
    unsigned long rank_up_to_cblock = ((m_cblock_header[header_ptr] >> (32 * (c - c2))) & 0xffffffff);

    long list_beg, list_size, lookup_bits;
    if (c == c2) {
      list_beg = (m_cblock_header[header_ptr + 1] & k_2cblock_size_mask);
      long next_list_beg = ((m_cblock_header[header_ptr + 1] >> (k_cblock_size_log + 1)) & k_2cblock_size_mask);
      list_size = next_list_beg - list_beg;
      lookup_bits = 3 * ((m_cblock_header[header_ptr + 1] >> (3 * (k_cblock_size_log + 1))) & 7);
    } else {
      list_beg = ((m_cblock_header[header_ptr + 1] >> (k_cblock_size_log + 1)) & k_2cblock_size_mask);
      list_size = ((m_cblock_header[header_ptr + 1] >> (2 * (k_cblock_size_log + 1))) & k_2cblock_size_mask);
      lookup_bits = 3 * ((m_cblock_header[header_ptr + 1] >> (3 * (k_cblock_size_log + 2))) & 7);
    }

    /*long temp_list_beg[2], temp_list_size[2], temp_lookup_bits[2];
    unsigned long val = m_cblock_header[header_ptr + 1];
    temp_list_beg[0] = (val & k_2cblock_size_mask); val >>= (k_cblock_size_log + 1);
    temp_list_beg[1] = (val & k_2cblock_size_mask); val >>= (k_cblock_size_log + 1);
    temp_list_size[0] = temp_list_beg[1] - temp_list_beg[0];
    temp_list_size[1] = (val & k_2cblock_size_mask); val >>= (k_cblock_size_log + 1);
    temp_lookup_bits[0] = 3 * (val & 7); val >>= 3;
    temp_lookup_bits[1] = 3 * (val & 7);*/

    /*long temp_list_beg[2];
    temp_list_beg[0] = (m_cblock_header[header_ptr + 1] & k_2cblock_size_mask);
    temp_list_beg[1] = ((m_cblock_header[header_ptr + 1] >> (k_cblock_size_log + 1)) & k_2cblock_size_mask);
    long temp_list_size[2];
    temp_list_size[0] = temp_list_beg[1] - temp_list_beg[0];
    temp_list_size[1] = ((m_cblock_header[header_ptr + 1] >> (2 * (k_cblock_size_log + 1))) & k_2cblock_size_mask);
    long temp_lookup_bits[2];
    temp_lookup_bits[0] = 3 * ((m_cblock_header[header_ptr + 1] >> (3 * (k_cblock_size_log + 1))) & 7);
    temp_lookup_bits[1] = 3 * ((m_cblock_header[header_ptr + 1] >> (3 * (k_cblock_size_log + 2))) & 7);*/

    // XXX consider actually getting the value of m_cblock_header[header_ptr + 1] into variable
    // and shifting it to the right, then anding, then shifting and anding and so on.

//    long list_beg = temp_list_beg[c - c2];
//    long list_size = temp_list_size[c - c2];
    if (!list_size) return rank_up_to_cblock;

//    long lookup_bits = temp_lookup_bits[c - c2];

//    long rank_up_to_cblock = (m_cblock_header[(cblock_id << k_sigma_log) + c] >> (k_cblock_size_log + 4));

    //--------------------------------------------------------------------------
    // STEP 2: compute the number of occurrences of c inside the cblock.
    //--------------------------------------------------------------------------

    // 1
    //
    // Decode the beginning and end of c's occurrence list.
    /*long list_beg = ((m_cblock_header[(cblock_id << k_sigma_log) + c] >> 3) & k_2cblock_size_mask);
    long list_end = ((c == k_sigma - 1) ? k_cblock_size :
        ((m_cblock_header[(cblock_id << k_sigma_log) + c + 1] >> 3) & k_2cblock_size_mask));
    if (list_beg == list_end) return rank_up_to_cblock;*/

    // 2
    //
    // Compute the distance from i to the closest reference point on the left.
    // long lookup_bits = 3 * (m_cblock_header[(cblock_id << k_sigma_log) + c] & 7);

    long refpoint_dist_log = 31 - lookup_bits;
    long refpoint_disk_mask = (1 << refpoint_dist_log) - 1;
    long i_refpoint_offset = (cblock_i & refpoint_disk_mask);

    // 3
    //
    // Compute threshold of symbol c inside the current cblock.
    // The threshold is a small power of two that is >= max block size,
    // where max block size is the maximal block size for symbol c in
    // the current cblock.
    long threshold = (1 << (k_cblock_size_log - lookup_bits + 3));

    // 4
    //
    // Compute the id if block containing i. Note: we only have id.
    // Computing boundaries would require division.
    // long list_size = list_end - list_beg;
    long approx = ((cblock_i * list_size) >> k_cblock_size_log);

    // 5
    //
    // Extract the lookup table entry.
    long lookup_mask = (1 << lookup_bits) - 1;
    long begin = (m_trunk[cblock_beg + list_beg + approx] & lookup_mask);
    long next_block_begin =  (approx + 1 == list_size) ? list_size : (m_trunk[cblock_beg + list_beg + approx + 1] & lookup_mask);

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
      while (begin < next_block_begin && (m_trunk[cblock_beg + list_beg + begin] >> lookup_bits) < i_refpoint_offset)
        ++begin;

      return rank_up_to_cblock + begin;
    } else {
      // This is executed very rarely, so we move more
      // expensive code (i.e., another branch) here.
      if (begin == next_block_begin || (m_trunk[cblock_beg + list_beg + begin] >> lookup_bits) < (2 * threshold)) {
        // The value in the occ list was small -> the ref
        // point for i and for the block are the same, we
        // proceed as before, without modifying i_refpoint_offset.
        while (begin < next_block_begin && (m_trunk[cblock_beg + list_beg + begin] >> lookup_bits) < i_refpoint_offset)
          ++begin;

        return rank_up_to_cblock + begin;
      } else {
        // Block occurrences were encoded wrt to the previous
        // ref point -> we increase i_refpoint_offset by
        // refpoint_dist and proceed as before.
        i_refpoint_offset += (1 << refpoint_dist_log);
        while (begin < next_block_begin && (m_trunk[cblock_beg + list_beg + begin] >> lookup_bits) < i_refpoint_offset)
          ++begin;

        return rank_up_to_cblock + begin;
      }
    }
  }

  ~rank4n() {
    delete[] m_count;
    delete[] m_cblock_header;
    free(m_trunk);
  }

  static const int k_cblock_size = (1 << k_cblock_size_log);
  static const int k_cblock_size_mask = k_cblock_size - 1;
  static const int k_cblock_size_mask_neg = (~k_cblock_size_mask);
  static const int k_2cblock_size = (2 << k_cblock_size_log);
  static const int k_2cblock_size_mask = k_2cblock_size - 1;
  static const int k_sigma = (1 << k_sigma_log);
  static const int k_sigma_mask = k_sigma - 1;

  long m_length;
  long n_cblocks;

  unsigned long *m_count;
  unsigned long *m_cblock_header;
  unsigned *m_trunk;
};

#endif // __RANK_H_INCLUDED
