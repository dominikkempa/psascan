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
    m_count = new long[k_sigma];
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

//      fprintf(stderr, "Processing cblock [%ld..%ld)\n", cblock_beg, cblock_end);

      //------------------------------------------------------------------------
      // STEP 1: compute symbol counts inside cblock.
      //------------------------------------------------------------------------
      std::fill(cblock_count, cblock_count + k_sigma, 0L);
      for (long i = cblock_beg; i < cblock_end; ++i) {
        unsigned char c = (i < m_length ? text[i] : 0);
        ++cblock_count[c];
      }

      //------------------------------------------------------------------------
      // STEP 3: compute starting positions of symbol occurrence lists.
      //------------------------------------------------------------------------
      for (long j = 0, t, s = 0; j < k_sigma; ++j) {
        t = cblock_count[j];
        list_beg[j] = s;
        s += t;
      }

      //------------------------------------------------------------------------
      // STEP 4: compute the rest of cblock header: ranks up to cblock
      //         beginning pointers to beginnings of occurrence lists.
      //------------------------------------------------------------------------
      for (long c = 0; c < k_sigma; ++c) {
        m_cblock_header[(cblock_id << k_sigma_log) + c] |= (list_beg[c] << 5);
        m_cblock_header[(cblock_id << k_sigma_log) + c] |= (m_count[c] << (k_cblock_size_log + 6));
      }

      //------------------------------------------------------------------------
      // STEP 5: update global symbol counts.
      //------------------------------------------------------------------------
      for (long i = cblock_beg; i < cblock_end; ++i) {
        unsigned char c = (i < m_length ? text[i] : 0);
        ++m_count[c];
      }

      for (long c = 0; c < k_sigma; ++c)
        occ[c].clear();

      //------------------------------------------------------------------------
      // STEP 6: build lists of occurrences.
      //------------------------------------------------------------------------
      std::fill(cblock_count, cblock_count + k_sigma, 0L);
      unsigned *cblock_trunk = m_trunk + cblock_beg;
      for (long i = cblock_beg; i < cblock_end; ++i) {
        unsigned char c = (i < m_length ? text[i] : 0);
        long cblock_i = i - cblock_beg;

        cblock_count[c]++;
        occ[c].push_back(cblock_i);
      }

      //------------------------------------------------------------------------
      // STEP 7: compute list of occurrences and lookup tables.
      //------------------------------------------------------------------------
      for (long c = 0; c < k_sigma; ++c) {
        long freq = occ[c].size();

//        fprintf(stderr, "\tProcessing symbol %ld\n", c);
//        fprintf(stderr, "\tocc list: ");
//        for (long jjj = 0; jjj < freq; ++jjj)
//          fprintf(stderr, "%ld ", occ[c][j]);
//        fprintf(stderr, "\n");

        // Compute the number of bits necessary to encode lookup table entries.
        // Lookup table can store values in the range [0..freq], thus we need
        // ceil(log2(freq + 1)) bits. This values is called lookup_bits and we
        // store it in the block header (using 5 bits).
        long lookup_bits = utils::log2ceil(freq + 1);
        m_cblock_header[(cblock_id << k_sigma_log) + c] |= lookup_bits;

//        fprintf(stderr, "\tlookup_bits = %ld\n", lookup_bits);

        // Compute log2 of the distance between two referece points.
        long refpoint_dist_log = 31 - lookup_bits;
        long refpoint_dist = (1L << refpoint_dist_log);
        long refpoint_dist_mask = refpoint_dist - 1;

//        fprintf(stderr, "\trefpoint_dist = %ld\n", refpoint_dist);

        // Used to compute the closest reference point on the left.
        long refpoint_dist_mask_neg = (~refpoint_dist_mask);

        // For each block compute:
        //   * its boundaries (beg, end) and then compute
        //   * the value of lookup table at entry corresponding to
        //     this block
        //   * encode the list of occurrences of that symbol inside
        //     the block with respect to the reference point
        //     associated with the block. The reference point
        //     associated with the block is largest reference
        //     point p satisfying p <= block_beg, i.e., we encode
        //     occurrences inside a block wrt. to the closest
        //     reference point on the left.
        long occ_ptr = 0;
        for (long j = 0; j < freq; ++j) {
          // Process j-th block.

//          fprintf(stderr, "\t\tEncoding %ld-th block\n", j);

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

//          if (c == 0 && j == 43155)
//            fprintf(stderr, "\t\tBlock boundaries: [%ld..%ld)\n", block_beg, block_end);
//          if (block_beg >= 29719 && block_beg <= 29725) {
//            fprintf(stderr, "\nblock_beg = %ld, block_end = %ld\n", block_beg, block_end);
//          }


          // 2
          //
          // Find the range of elements from the current block inside occ[c].
          while (occ_ptr < freq && occ[c][occ_ptr] < block_beg) ++occ_ptr;
          long range_beg = occ_ptr;
          while (occ_ptr < freq && occ[c][occ_ptr] < block_end) ++occ_ptr;
          long range_end = occ_ptr;

//          if (range_beg >= 29719 && range_beg <= 29725) {
//            fprintf(stderr, "block_beg = %ld, block_end = %ld, range_beg = %ld, range_end = %ld\n",
//                block_beg, block_end, range_beg, range_end);
//          }

//          fprintf(stderr, "\t\tBlock range: [%ld..%ld)\n", range_beg, range_end);

          // 3
          //
          // Store the value in the lookup table.
          cblock_trunk[list_beg[c] + j] |= range_beg;
//          fprintf(stderr, "storing %ld in the lookup table at index %ld\n",
//              range_beg, j);

          // 4
          //
          // Add the occurrences occ[c][range_beg..range_end) to the
          // list of c's occurrences in the trunk. Encode them with
          // respect to the closes reference point on the left of block_beg.
          long closest_ref_point = (block_beg & refpoint_dist_mask_neg);
          for (long occ_id = range_beg; occ_id < range_end; ++occ_id) {
            cblock_trunk[list_beg[c] + occ_id] |= ((occ[c][occ_id] - closest_ref_point) << lookup_bits);
//            fprintf(stderr, "\t\tStoring value %ld\n", (occ[c][occ_id] - closest_ref_point));

            /*if (block_beg >= 29719 && block_beg <= 29725) {
              fprintf(stderr, "written offset of occ[c = %ld][occ_id = %ld] == %ld minus closect_ref_point == %ld "
                  " at cblock_trunk[list_beg[c] == %ld + occ_id == %ld]\n", (long)c, occ_id, (long)occ[c][occ_id],
                  closest_ref_point, list_beg[c], occ_id);
              fprintf(stderr, "range_beg = %ld, range_end = %ld\n", range_beg, range_end);
            }*/

          }

//          fprintf(stderr, "\t\tClosest reference point: %ld\n", closest_ref_point);
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


//    fprintf(stderr, "QUERY: i = %ld, c = %ld\n", i, (long)c);
//    fprintf(stderr, "\tcblock_i = %ld, cblock_beg = %ld, cblock_id = %ld\n",
//        cblock_i, cblock_beg, cblock_id);

    //--------------------------------------------------------------------------
    // STEP 1: extract the rank up to the start of cblock.
    //--------------------------------------------------------------------------
    long rank_up_to_cblock = (m_cblock_header[(cblock_id << k_sigma_log) + c] >> (k_cblock_size_log + 6));

//    fprintf(stderr, "\trank up to cblock = %ld\n", rank_up_to_cblock);

    //--------------------------------------------------------------------------
    // STEP 2: compute the number of occurrences of c inside the cblock.
    //--------------------------------------------------------------------------

    // 1
    //
    // Decode the beginning and end of c's occurrence list.
    long list_beg = ((m_cblock_header[(cblock_id << k_sigma_log) + c] >> 5) & k_2cblock_size_mask);
    long list_end = ((c == k_sigma - 1) ? k_cblock_size :
        ((m_cblock_header[(cblock_id << k_sigma_log) + c + 1] >> 5) & k_2cblock_size_mask));

//    fprintf(stderr, "\tlist_beg = %ld, list_end = %ld\n", list_beg, list_end);

    if (list_beg == list_end) return rank_up_to_cblock;

    // 2
    //
    // Compute the distance from i to the closest reference point on the left.
    long lookup_bits = (m_cblock_header[(cblock_id << k_sigma_log) + c] & 31);
    long refpoint_dist_log = 31 - lookup_bits;
    long refpoint_disk_mask = (1 << refpoint_dist_log) - 1;
    long i_refpoint_offset = (cblock_i & refpoint_disk_mask);

//    fprintf(stderr, "\tlookup_bits = %ld, refpoint_dist_log = %ld, i_refpoint_offset = %ld\n",
//        lookup_bits, refpoint_dist_log, i_refpoint_offset);

    // 3
    //
    // Compute threshold of symbol c inside the current cblock.
    // The threshold is the small power of two that is >= max block size,
    // where max block size is the maximal block size for symbol c in the
    // current cblock.
    long threshold = (1 << (k_cblock_size_log - lookup_bits + 1));

//    fprintf(stderr, "\tthreshold = %ld\n", threshold);

    // 4
    //
    // Interpolate i, i.e., compute the index of lookup table entry.
    long list_size = list_end - list_beg;
    long approx = ((cblock_i * list_size) >> k_cblock_size_log);

//    fprintf(stderr, "\tapprox = %ld\n", approx);

    // 5
    //
    // Extract the lookup table entry.
    long lookup_mask = (1 << lookup_bits) - 1;
    long begin = (m_trunk[cblock_beg + list_beg + approx] & lookup_mask);
    long next_block_begin =  (approx + 1 == list_size) ? list_size : (m_trunk[cblock_beg + list_beg + approx + 1] & lookup_mask);

//    fprintf(stderr, "\tbegin = %ld\n", begin);
//    fprintf(stderr, "\tnext_block_begin = %ld\n", next_block_begin);
//    fprintf(stderr, "\tlooked up begin at m_trunk[%ld]\n", cblock_beg + list_beg + approx);
//    fprintf(stderr, "\tobtained begin by shifting m_trunk[..] by %ld bits right\n",
//        lookup_bits);

    // 6
    //
    // Correct the value of begin and return the answer. To correct begin
    // we need to know wrt to which reference point have we encoded the
    // occurrences of c inside the block. This would be easy to compute
    // if we knew the block beginning (but we don't). However, we can
    // infer the unknown reference point from i and/or the first element
    // on the occurrence list (or rather whether the first element is big
    // or small).
    if (i_refpoint_offset >= threshold) {

//      fprintf(stderr, "\tcase I:\n");
//      fprintf(stderr, "\tfirst element in the trunk: %u\n", (m_trunk[cblock_beg + list_beg + begin] >> lookup_bits));
//      fprintf(stderr, "\ti_refpoint_offset = %ld\n", i_refpoint_offset);

      // This is the easy case which will happen most of the time, so
      // we should get good branch prediction here.
      while (begin < next_block_begin && (m_trunk[cblock_beg + list_beg + begin] >> lookup_bits) < i_refpoint_offset)
        ++begin;

//      fprintf(stderr, "\tcomputed begin = %ld\n", begin);
      return rank_up_to_cblock + begin;
    } else {
      // This is executed very rarely, so we move more expensive
      // code (i.e., another branch) here.

//      fprintf(stderr, "\tfirst element in the trunk = %ld\n", (long)(m_trunk[cblock_beg + list_beg + begin] >> lookup_bits));
//      fprintf(stderr, "\tbegin = %ld, next_block_begin = %ld\n", begin, next_block_begin);

      if (begin == next_block_begin || (m_trunk[cblock_beg + list_beg + begin] >> lookup_bits) < (2 * threshold)) {

//        fprintf(stderr, "\tcase II:\n");
//        fprintf(stderr, "\tfirst elem in the list: %u\n", (m_trunk[cblock_beg + list_beg + begin] >> lookup_bits));
//        fprintf(stderr, "\ti_refpoint_offset = %ld\n", i_refpoint_offset);
//        fprintf(stderr, "\ttest = %ld\n",
//          (long)((m_trunk[cblock_beg + list_beg + begin] >> lookup_bits) < i_refpoint_offset));
       

        // The value in the occ list was small -> the ref point for i and
        // for the block are the same, we proceed as before, without modifying i_refpoint_offset.
        while (begin < next_block_begin && (m_trunk[cblock_beg + list_beg + begin] >> lookup_bits) < i_refpoint_offset)
          ++begin;

//        fprintf(stderr, "\tbegin after update: %ld\n", begin);

        return rank_up_to_cblock + begin;
      } else {

//        fprintf(stderr, "\tcase III:\n");

        // Block occurrences were encoded wrt to the previous ref point -> we
        // increase i_refpoint_offset by refpoint_dist and proceed as before.
        i_refpoint_offset += (1 << refpoint_dist_log);
//        fprintf(stderr, "\tincreasing i_refpoint_offset by %ld up to value %ld\n",
//            (long)(1 << refpoint_dist_log), (long)i_refpoint_offset);

        while (begin < next_block_begin && (m_trunk[cblock_beg + list_beg + begin] >> lookup_bits) < i_refpoint_offset) {
//          fprintf(stderr, "\t!!! increasing begin\n");
//          fprintf(stderr, "\tbecause m_trunk[%ld] == %u\n", cblock_beg + list_beg + begin, (m_trunk[cblock_beg + list_beg + begin] >> lookup_bits));
//          fprintf(stderr, "\trecall that i _refpoint_offset = %ld\n", i_refpoint_offset);
          ++begin;
        }
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

  long *m_count;
  unsigned long *m_cblock_header;
  unsigned *m_trunk;
};

#endif // __RANK_H_INCLUDED
