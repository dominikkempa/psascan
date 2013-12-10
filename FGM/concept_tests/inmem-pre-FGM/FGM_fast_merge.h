// TODO: write SA directly to disk.
#ifndef __FGM_FAST_MERGE
#define __FGM_FAST_MERGE

#include <algorithm>

#include "utils.h"
#include "file_streamer.h"
#include "vbyte_file_streamer.h"

void FGM_fast_merge(int *outputSA, int length, int max_block_size) {
  int n_block = (length + max_block_size - 1) / max_block_size;
  
  // Initialize buffers for merging.
  file_streamer<int> **sparseSA = new file_streamer<int>*[n_block];
  vbyte_file_streamer **gap = new vbyte_file_streamer*[n_block];
//  file_streamer<int> **gap = new file_streamer<int>*[n_block];
  for (int i = 0; i < n_block; ++i) {
    sparseSA[i] = new file_streamer<int>("sparseSA." + utils::intToStr(n_block - 1 - i), 1 << 20);
//    gap[i] = new file_streamer<int>("gap." + utils::intToStr(n_block - 1 - i), 1 << 20);
    gap[i] = new vbyte_file_streamer("gap." + utils::intToStr(n_block - 1 - i), 1 << 20);
  }
  
  // Merge.
  int *block_rank  = new int[n_block];
  int *suffix_rank = new int[n_block];
  int *suf_ptr = new int[n_block]; // First non-extracted suffix in the block.
  for (int i = 0; i < n_block; ++i) suffix_rank[i] = gap[i]->read();
  std::fill(block_rank, block_rank + n_block, 0);
  std::fill(suf_ptr, suf_ptr + n_block, 0);
  for (int i = 0; i < length; ++i) {
    // Find the leftmost block j with block_rank[j] == suffix_rank[j].
    int j = 0;
    while (j < n_block && block_rank[j] != suffix_rank[j]) ++j;

    // Extract the suffix.
    outputSA[i] = sparseSA[j]->read();

    // Update suffix_rank[j].    
    suffix_rank[j]++;
    suffix_rank[j] += gap[j]->read();
    
    // Update block_rank[0..j].
    for (int k = 0; k <= j; ++k) ++block_rank[k];
  }

  // Clean up.
  for (int i = 0; i < n_block; ++i) {
    delete sparseSA[i];
    delete gap[i];
  }
  delete[] gap;
  delete[] sparseSA;
}

#endif // __FGM_FAST_MERGE

