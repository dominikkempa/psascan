#ifndef __MERGE
#define __MERGE

#include <string>
#include <algorithm>

#include "utils.h"
#include "stream.h"

void merge(long length, long max_block_size, std::string out_filename) {
  stream_writer<int> *output = new stream_writer<int>(out_filename, 1 << 21);
  int n_block = (length + max_block_size - 1) / max_block_size;
  
  // Initialize buffers for merging.
  stream_reader<int> **sparseSA = new stream_reader<int>*[n_block];
  vbyte_stream_reader **gap = new vbyte_stream_reader*[n_block];
  for (int i = 0; i < n_block; ++i) {
    sparseSA[i] = new stream_reader<int>("sparseSA." + utils::intToStr(i), 1 << 20);
    gap[i] = new vbyte_stream_reader("gap." + utils::intToStr(i), 1 << 20);
  }
  
  // Merge.
  int *block_rank  = new int[n_block];
  int *suffix_rank = new int[n_block];
  int *suf_ptr = new int[n_block]; // First non-extracted suffix in the block.
  for (int i = 0; i < n_block; ++i) suffix_rank[i] = gap[i]->read();
  std::fill(block_rank, block_rank + n_block, 0);
  std::fill(suf_ptr, suf_ptr + n_block, 0);
  for (long i = 0; i < length; ++i) {
    // Find the leftmost block j with block_rank[j] == suffix_rank[j].
    int j = 0;
    while (j < n_block && block_rank[j] != suffix_rank[j]) ++j;

    // Extract the suffix.
    output->write(sparseSA[j]->read() + j * max_block_size); // SA[i]

    // Update suffix_rank[j].    
    suffix_rank[j]++;
    suffix_rank[j] += gap[j]->read();
    
    // Update block_rank[0..j].
    for (int k = 0; k <= j; ++k) ++block_rank[k];
  }

  // Clean up.
  delete output;
  for (int i = 0; i < n_block; ++i) {
    delete sparseSA[i];
    delete gap[i];
  }
  delete[] gap;
  delete[] sparseSA;
}

#endif // __MERGE

