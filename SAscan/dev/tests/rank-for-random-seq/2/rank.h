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
    long *block_size = new long[k_sigma];    // block size for each symbol
    long *block_size_log = new long[k_sigma];

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
      // STEP 2: compute the block size (different for every symbol inside the
      //         current cblock) and store in the cblock header.
      //------------------------------------------------------------------------
      // We require that the block size for any symbol is a power of two, so we
      // reduce the computation to log2 of block size. For a given
      // block_size_log, the lookup table requires
      //
      //   k_cblock_size / 2 ^ block_size_log
      //
      // entries. We want to lookup table entries to be perfectly aligned with
      // the list of c's occurrences, so we only have cblock_count[c] entries
      // available. Therefore, to make sure we can accommodate the whole lookup
      // table we need block_size_log to be the smallest integer satisfying
      //
      //   k_cblock_size / 2 ^ block_size_log <= cblock_count[c]    <=>
      //   2 ^ block_size_log * cblock_count[c] >= k_cblock_size    <=>
      //   (1 << block_size_log) * cblock_count[c] >= k_cblock_size
      //
      // Mathematically block_size_log is equal to
      //
      //   block_size_log = k_cblock_size_log - floor(log2(cblock_count[c])) (1)
      //
      // The is the formula used below. Note that the elements on the c's list
      // of occurrences are offsets with respect to the beginning of the block,
      // so they will be encoded using only block_size_log bits.
      //
      // On the other hand, if we encode the lookup table entries as offsets
      // with respect to the beginning of the list of c's occurrences, the
      // values will be in the range [0 .. cblock_count[c]] (note, that i-th
      // entry of the lookup table tells the position (or rather an ofset with
      // respect to the beginning of the block; as expalined above) of the
      // first element in the i-th block of size 2 ^ block_size_log. Therefore,
      // if we use r bits to encode values in the lookup table, r need to
      // satisfy
      //
      //   2 ^ r - 1 >= cblock_count[c]    <=>
      //   2 ^ r     >= cblock_count[c] + 1
      //
      // Thus mathematically
      //
      //   r = floor(log2(cblock_count[c] + 1))
      //     <= 1 + floor(log2(cblock_count[c]))                             (2)
      //
      // From (1) and (2) we see that the total number of bits per element
      // in the trunk for symbol c is 
      //
      //    r + block_size_log <= 1 + floor(log2(cblock_count[c])) +
      //                       k_cblock_size_log - floor(log2(cblock_count[c]))
      //                       <= 1 + k_cblock_size_log <= 21
      //
      // Note: during rank query, we need to know both r and block_size_log
      // but clearly it suffices to store floor_log2_cblock_count to compute
      // both. floor_log2_cblock_count is a small values between 0 and 20, so
      // we reserve 5 bits in the cblock header to store it. There will be 64
      // bits reserved for every symbol in the header of every cblock, and the
      // 5 bits will be taken from these 64.
      //------------------------------------------------------------------------
      std::fill(block_size, block_size + k_sigma, 0L);
      for (long c = 0; c < k_sigma; ++c) {
        if (!cblock_count[c]) continue;
        long floor_log2_cblock_count = utils::log2floor(cblock_count[c]);    

        block_size_log[c] = k_cblock_size_log - floor_log2_cblock_count;
        block_size[c] = (1 << block_size_log[c]);
        m_cblock_header[(cblock_id << k_sigma_log) + c] = block_size_log[c];
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
        if (!cblock_count[c]) m_cblock_header[(cblock_id << k_sigma_log) + c] = 0;
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
      // STEP 6: compute the trunk.
      //------------------------------------------------------------------------
      std::fill(cblock_count, cblock_count + k_sigma, 0L);
      unsigned *cblock_trunk = m_trunk + cblock_beg;
      for (long i = cblock_beg; i < cblock_end; ++i) {
        long cblock_i = i - cblock_beg;
        unsigned char c = (i < m_length ? text[i] : 0);

        // Add c to its list of occurrences.
        long block_i = i % block_size[c]; // offset in the block
        cblock_trunk[list_beg[c] + cblock_count[c]++] |= block_i;
        occ[c].push_back(cblock_i);
      }

      for (long c = 0; c < 256; ++c) {
        long blocksize = block_size[c];
        if (!blocksize) continue;
        long blocksize_log = block_size_log[c];
        long ptr = 0;
        for (long blockbeg = 0, t = 0; blockbeg < k_cblock_size; blockbeg += blocksize, ++t) {
          while (ptr < (long)occ[c].size() && occ[c][ptr] < blockbeg) ++ptr;
          cblock_trunk[list_beg[c] + t] |= (ptr << blocksize_log);
        }
      }
    }
    m_count[0] -= n_cblocks * k_cblock_size - m_length; // remove extra zeros
    delete[] occ;
    delete[] cblock_count;
    delete[] block_size;
    delete[] block_size_log;
    delete[] list_beg;
  }
  
  inline long rank(long i, unsigned char c) {
    if (i <= 0) return 0L;
    else if (i >= m_length) return m_count[c];

    long cblock_id = (i >> k_cblock_size_log); // which cblock
    long cblock_beg = (i & k_cblock_size_mask_neg);

    long block_size_log = (m_cblock_header[(cblock_id << k_sigma_log) + c] & 31);
    long block_size = (1 << block_size_log);
    long block_size_mask = block_size - 1;

    long block_i =  (i & block_size_mask);        // offset in block
    long cblock_i = (i & k_cblock_size_mask);     // offset in cblock
    long block_id = (cblock_i >> block_size_log); // which block in cblock

    //--------------------------------------------------------------------------
    // STEP 1: extract the rank up to the start of cblock.
    //--------------------------------------------------------------------------
    long rank_up_to_cblock = (m_cblock_header[(cblock_id << k_sigma_log) + c] >> (k_cblock_size_log + 6));

    //--------------------------------------------------------------------------
    // STEP 2: compute the number of occurrences of c inside the cblock.
    //--------------------------------------------------------------------------
    long header_idx = (cblock_id << k_sigma_log) + c;
    long list_beg = ((m_cblock_header[header_idx] >> 5) & k_2cblock_size_mask);
    long list_end = ((c == k_sigma - 1) ? k_cblock_size : ((m_cblock_header[header_idx + 1] >> 5) & k_2cblock_size_mask));
    if (list_beg == list_end) return rank_up_to_cblock;

    // Find the first element in the c's list that is >= cblock_i.
    // Compute indices such that all occurrences of c inside block block_id
    // are within that range [begin..end). The indices are indexed from the
    // beginning of the list, so to access them in the trunk we need to look
    // at m_trunk[cblock_beg + list_beg + begin..cblock_beg + list_beg + end).
    long begin = (m_trunk[cblock_beg + list_beg + block_id] >> block_size_log);
    long end = (cblock_i + block_size >= k_cblock_size) ?
      list_end - list_beg : (m_trunk[cblock_beg + list_beg + block_id + 1] >> block_size_log);
    while (begin < end && (m_trunk[cblock_beg + list_beg + begin] & block_size_mask) < block_i) ++begin;

    return rank_up_to_cblock + begin;
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
