#ifndef __RANK_H_INCLUDED
#define __RANK_H_INCLUDED

#include <algorithm>
#include <vector>

#include "utils.h"


// XXX consider reversing the order in the trunk.
template<long k_cblock_size_bits = 18, long k_sigma_bits = 8>
struct rank4n {
  rank4n(unsigned char *text, long length, long) {
    m_length = length;

    // Global symbol counts.
    m_count = new long[k_sigma];
    std::fill(m_count, m_count + k_sigma, 0L);

    // Compute trunk and block header.
    n_cblocks = (m_length + k_cblock_size - 1) / k_cblock_size;
    m_trunk = (unsigned *)malloc(n_cblocks * k_cblock_size * sizeof(unsigned));
    std::fill(m_trunk, m_trunk + n_cblocks * k_cblock_size, 0U);
    m_cblock_header = new long[k_sigma * n_cblocks];
    long *cblock_count = new long[k_sigma];

    for (long cblock_id = 0; cblock_id < n_cblocks; ++cblock_id) {
      long cblock_beg = k_cblock_size * cblock_id;
      long cblock_end = cblock_beg + k_cblock_size;

      // Compute local (inside block) counts.
      std::fill(cblock_count, cblock_count + k_sigma, 0L);
      for (long i = cblock_beg; i < cblock_end; ++i) {
        unsigned char c = (i < m_length ? text[i] : 0);
        ++cblock_count[c];
      }

      // Partial sum over block counts, gives starting
      // positions of lists with symbol occurrences.
      for (long j = 0, t, s = 0; j < k_sigma; ++j) {
        t = cblock_count[j];
        cblock_count[j] = s;
        s += t;
      }

      // Compute cblock header: ranks up to cblock start and
      // pointers to beginnings of lists with symbols occurrences.
      for (long c = 0; c < k_sigma; ++c) {
        m_cblock_header[(cblock_id << k_sigma_bits) + c] = (m_count[c] << (k_cblock_size_bits + 1));
        m_cblock_header[(cblock_id << k_sigma_bits) + c] |= cblock_count[c];
      }

      // Now we can update global counts.
      for (long i = cblock_beg; i < cblock_end; ++i) {
        unsigned char c = (i < m_length ? text[i] : 0);
        ++m_count[c];
      }

      // Compute the list of occurrences and lookup
      // table for each symbol inside a cblock.
      unsigned *cblock_trunk = m_trunk + cblock_beg;
      for (long i = cblock_beg; i < cblock_end; ++i) {
        unsigned char c = (i < m_length ? text[i] : 0);

        // Update lookup tables.
        long cblock_local_i = (i & k_cblock_size_mask);
        if (!(i & k_sigma_mask))
          for (long j = 0; j < k_sigma; ++j)
            cblock_trunk[(j << k_lookup_table_size_bits) + (cblock_local_i >> k_sigma_bits)] |= cblock_count[j];

        // Add c to its list of occurrences.
        long block_local_i = (cblock_local_i & k_sigma_mask);
        cblock_trunk[cblock_count[c]++] |= (block_local_i << (k_cblock_size_bits + 1));
      }
    }

    m_count[0] -= n_cblocks * k_cblock_size - m_length; // remove extra zeros
    delete[] cblock_count;
  }


  inline long rank(long i, unsigned char c) {
    if (i <= 0) return 0L;
    else if (i >= m_length) return m_count[c];

    long cblock_id = (i >> k_cblock_size_bits);
    long cblock_beg = (i & k_cblock_size_mask_neg);
    long cblock_local_i = (i & k_cblock_size_mask);
    long block_local_i = (cblock_local_i & k_sigma_mask);

    // Extract the rank up to the start of cblock.
    long rank_up_to_cblock =
      (m_cblock_header[(cblock_id << k_sigma_bits) + c] >> (k_cblock_size_bits + 1));

    // Compute the number of occurrences of c inside the cblock.
    // The beginning of c's list.
    long start = (m_cblock_header[(cblock_id << k_sigma_bits) + c] & k_2cblock_size_mask);

    // Find the first element in the c's list that is >= cblock_local_i.
    long idx = (c << k_lookup_table_size_bits) + (cblock_local_i >> k_sigma_bits);
    long end = (m_trunk[cblock_beg + idx] & k_2cblock_size_mask);
    long endmax = (idx + 1 == k_cblock_size) ? k_cblock_size : (m_trunk[cblock_beg + idx + 1] & k_2cblock_size_mask);
    while (end < endmax && (m_trunk[cblock_beg + end] >> (k_cblock_size_bits + 1)) < block_local_i) ++end;

    return rank_up_to_cblock + (end - start);
  }

  ~rank4n() {
    delete[] m_count;
    delete[] m_cblock_header;
    free(m_trunk);
  }

  static const int k_cblock_size = (1 << k_cblock_size_bits);
  static const int k_2cblock_size = (2 << k_cblock_size_bits);
  static const int k_cblock_size_mask = k_cblock_size - 1;
  static const int k_2cblock_size_mask = k_2cblock_size - 1;
  static const int k_cblock_size_mask_neg = (~k_cblock_size_mask);
  static const int k_lookup_table_size_bits = (k_cblock_size_bits - k_sigma_bits);
  static const int k_sigma = (1 << k_sigma_bits);
  static const int k_sigma_mask = k_sigma - 1;

  long m_length;
  long n_cblocks;

  long *m_count;
  long *m_cblock_header;
  unsigned *m_trunk;
};

#endif // __RANK_H_INCLUDED
