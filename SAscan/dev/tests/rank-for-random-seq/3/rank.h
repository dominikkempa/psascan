#ifndef __RANK_H_INCLUDED
#define __RANK_H_INCLUDED

#include <algorithm>
#include <vector>

#include "utils.h"

template<long k_cblock_size_log = 16, long k_sigma_log = 8>
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
        m_cblock_header[(cblock_id << k_sigma_log) + c] |= list_beg[c];
        m_cblock_header[(cblock_id << k_sigma_log) + c] |= (m_count[c] << (k_cblock_size_log + 1));
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
        unsigned char c = (i < m_length ? text[i] : 0);
        long cblock_i = i - cblock_beg;

        // Add c to its list of occurrences, we store
        // the absolute position of c inside the cblock,
        // using 16 bits.
        cblock_trunk[list_beg[c] + cblock_count[c]++] |= cblock_i;
        occ[c].push_back(cblock_i);
      }

      // Compute lookup for each symbol.
      for (long c = 0; c < k_sigma; ++c) {
        long k = (long)occ[c].size();

        long occ_ptr = 0;
        for (long j = 0; j < k; ++j) {
          // Find the beginning is j-th block, i.e., the smallest
          // smallest i such that floor((i * k) / k_cblock_size) = j.
          // The last formula we evaluate: (i * k) >> k_cblock_size_log.
          long i = ((j << k_cblock_size_log) + k - 1) / k;
          while (i > 0 && (((i - 1) * k) >> k_cblock_size_log) == j) --i;

          // Now find the smallest element in the j-th block, i.e., the
          // smallest occ_ptr, such that occ[c][occ_ptr] >= i.
          while (occ_ptr < k && occ[c][occ_ptr] < i) ++occ_ptr;
          long ans = occ_ptr;

#if 0
          //--------------------------------------------------------------------
          // Experimental optimization, that actually doesn't help at all.
          //--------------------------------------------------------------------

          // Find the beginning of (j+1)-th block (or just k_cblock_size)
          long i2;
          if (j + 1 == k) i2 = k_cblock_size;
          else {
            i2 = (((j + 1) << k_cblock_size_log) + k - 1) / k;
            while (i2 > 0 && (((i2 - 1) * k) >> k_cblock_size_log) == j + 1) --i2;
          }

          // Find the first element in the (j + 1)-th block.
          long next_occ_ptr = occ_ptr;
          while (next_occ_ptr < k && occ[c][next_occ_ptr] < i2) ++next_occ_ptr;
          ans = (occ_ptr + next_occ_ptr) >> 1;
#endif

          // The lookup table at index j gives the
          // position of occ[c][ptr] in occ[c], so
          // we actually store just ptr. For this 16 bits
          // is always sufficient.
          cblock_trunk[list_beg[c] + j] |= (ans << 16);
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
    long rank_up_to_cblock = (m_cblock_header[(cblock_id << k_sigma_log) + c] >> (k_cblock_size_log + 1));

    //--------------------------------------------------------------------------
    // STEP 2: compute the number of occurrences of c inside the cblock.
    //--------------------------------------------------------------------------
    long list_beg = (m_cblock_header[(cblock_id << k_sigma_log) + c] & k_2cblock_size_mask);
    long list_end = ((c == k_sigma - 1) ? k_cblock_size : (m_cblock_header[(cblock_id << k_sigma_log) + c + 1] & k_2cblock_size_mask));
    if (list_beg == list_end) return rank_up_to_cblock;
    long list_size = list_end - list_beg;

    long approx = ((cblock_i * list_size) >> k_cblock_size_log);
    long begin = (m_trunk[cblock_beg + list_beg + approx] >> 16);
    while (begin < list_size && (m_trunk[cblock_beg + list_beg + begin] & 0xffff) < cblock_i)
      ++begin;

#if 0
    // Necessary for the experimental optimization.
    while (begin && (m_trunk[cblock_beg + list_beg + begin - 1] & 0xffff) >= cblock_i) --begin;
#endif

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
