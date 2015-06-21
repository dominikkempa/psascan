#ifndef __PARALLEL_RANK_H_INCLUDED
#define __PARALLEL_RANK_H_INCLUDED

#include <algorithm>
#include <thread>

#include "rank.h"

void initialize_rank(unsigned char *text, long length, context_rank_4n* &rank) {
  rank = new context_rank_4n(text, length);
}

struct parallel_rank {
  parallel_rank(unsigned char *text, long length, long max_threads) {
    m_block_size = (length + max_threads - 1) / max_threads;
    m_blocks = (length + m_block_size - 1) / m_block_size;

    // Build rank structures.
    ranks = new context_rank_4n*[m_blocks];
    std::thread **threads = new std::thread*[m_blocks];
    for (long i = 0; i < m_blocks; ++i) {
      long block_beg = i * m_block_size;
      long block_end = std::min(block_beg + m_block_size, length);
      long this_block_size = block_end - block_beg;
      threads[i] = new std::thread(initialize_rank,
          text + block_beg, this_block_size, std::ref(ranks[i]));
    }
    for (long i = 0; i < m_blocks; ++i) threads[i]->join();
    for (long i = 0; i < m_blocks; ++i) delete threads[i];
    delete[] threads;

    // Compute m_block_rank.
    m_block_rank = new long*[m_blocks];
    for (long i = 0; i < m_blocks; ++i) {
      m_block_rank[i] = new long[256];
      for (unsigned j = 0; j < 256; ++j) {
        if (i == 0) m_block_rank[i][j] = 0;
        else m_block_rank[i][j] = m_block_rank[i - 1][j] +
          ranks[i - 1]->rank(m_block_size, j);
      }
    }
  }

  inline long rank(long i, unsigned char c) {
    long block_id = 0, block_beg = 0;

    // Find the right block. Could also be done
    // by i / m_block_beg, but I think this is faster.
    while (block_id + 1 < m_blocks &&
        block_beg + m_block_size <= i) {
      ++block_id;
      block_beg += m_block_size;
    }

    // Ask a query in a block.
    return m_block_rank[block_id][c] +
      ranks[block_id]->rank(i - block_beg, c);
  }

  ~parallel_rank() {
    for (long i = 0; i < m_blocks; ++i) {
      delete[] m_block_rank[i];
      delete ranks[i];
    }
    delete[] ranks;
    delete[] m_block_rank;
  }

  context_rank_4n **ranks;

  long m_block_size;
  long m_blocks;
  long **m_block_rank;
};

#endif  // __ PARALLEL_RANK_H_INCLUDED
