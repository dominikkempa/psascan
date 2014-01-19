#ifndef __MERGE
#define __MERGE

#include <string>
#include <algorithm>

#include "utils.h"
#include "stream.h"
#include "uint40.h"

void merge(long length, long max_block_size, std::string out_filename) {
  long n_block = (length + max_block_size - 1) / max_block_size;
  long pieces = 5 * n_block + 5;
  long ram_use = 5L * max_block_size;
  long buffer_size = (ram_use + pieces - 1) / pieces;

  fprintf(stderr, "Buffer size for merging: %ld\n", buffer_size);
  stream_writer<uint40> *output = new stream_writer<uint40>(out_filename, 5 * buffer_size);

  // Initialize buffers for merging.
  stream_reader<int> **sparseSA = new stream_reader<int>*[n_block];
  vbyte_stream_reader **gap = new vbyte_stream_reader*[n_block];
  for (long i = 0; i < n_block; ++i) {
    sparseSA[i] = new stream_reader<int>("sparseSA." + utils::intToStr(i), 4 * buffer_size);
    gap[i] = new vbyte_stream_reader("gap." + utils::intToStr(i), buffer_size);
  }
  
  // Merge.
  long *block_rank  = new long[n_block];
  long *suffix_rank = new long[n_block];
  long *suf_ptr = new long[n_block]; // First non-extracted suffix in the block.
  for (long i = 0; i < n_block; ++i) suffix_rank[i] = gap[i]->read();
  std::fill(block_rank, block_rank + n_block, 0);
  std::fill(suf_ptr, suf_ptr + n_block, 0);
  fprintf(stderr, "Merging:\r");
  long double merge_start = utils::wclock();
  for (long i = 0, dbg = 0; i < length; ++i, ++dbg) {
    if (dbg == (1 << 23)) {
      long double elapsed = utils::wclock() - merge_start;
      fprintf(stderr, "Merging: %.1Lf%%. Time: %.2Lfs\r",
          (100.L * i) / length, elapsed);
      dbg = 0;
    }
    // Find the leftmost block j with block_rank[j] == suffix_rank[j].
    long j = 0;
    while (j < n_block && block_rank[j] != suffix_rank[j]) ++j;

    // Extract the suffix.
    output->write(uint40((unsigned long)sparseSA[j]->read() + (unsigned long)max_block_size * j)); // SA[i]

    // Update suffix_rank[j].    
    suffix_rank[j]++;
    suffix_rank[j] += gap[j]->read();
    
    // Update block_rank[0..j].
    for (long k = 0; k <= j; ++k) ++block_rank[k];
  }
  long double merge_time = utils::wclock() - merge_start;
  fprintf(stderr, "Merging: 100.0%%. Time: %.2Lfs\n", merge_time);

  // Clean up.
  delete output;

  for (long i = 0; i < n_block; ++i) {
    delete sparseSA[i];
    delete gap[i];
  }
  delete[] gap;
  delete[] sparseSA;
  
  // Delete auxiliary files.
  for (int i = 0; i < n_block; ++i) {
    utils::file_delete("sparseSA." + utils::intToStr(i));
    utils::file_delete("gap." + utils::intToStr(i));
  }
  
  delete[] block_rank;
  delete[] suffix_rank;
  delete[] suf_ptr;
}

#endif // __MERGE

