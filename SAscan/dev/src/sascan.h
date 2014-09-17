#ifndef __SASCAN_H_INCLUDED
#define __SASCAN_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <string>
#include <algorithm>

#include "partial_sufsort.h"
#include "merge.h"
#include "utils.h"
#include "uint40.h"
#include "settings.h"



void SAscan(std::string input_filename, std::string output_filename, long ram_use, long max_threads,
    long stream_buffer_size = (1L << 21)) {
  long n_stream_buffers = 2 * max_threads;
  if (ram_use < 6L) {
    fprintf(stderr, "Error: not enough memory to run SAscan.\n");
    std::exit(EXIT_FAILURE);
  }
  
  // Turn paths absolute.
  input_filename = utils::absolute_path(input_filename);
  output_filename = utils::absolute_path(output_filename);
  long length = utils::file_size(input_filename);
  fprintf(stderr, "Input filename = %s\n", input_filename.c_str());
  fprintf(stderr, "Output filename = %s\n", output_filename.c_str());
  fprintf(stderr, "Input length = %ld (%.1LfMiB)\n", length, 1.L * length / (1L << 20));
  fprintf(stderr, "\n");

  long ram_use_excluding_threads = ram_use - n_stream_buffers * stream_buffer_size;
  if (ram_use_excluding_threads < 6L) {
    long required_bytes = n_stream_buffers * stream_buffer_size;
    long required_MiB = (required_bytes + (1L << 20) - 1) / (1L << 20);
    fprintf(stderr, "Error: not enough memory to start threads. You need "
        "at least %ldMiB\n", required_MiB + 1);
    std::exit(EXIT_FAILURE);
  }
    
  fprintf(stderr, "RAM budget = %ld (%.1LfMiB)\n", ram_use, 1.L * ram_use / (1L << 20));
  fprintf(stderr, "RAM budget (excluding threads) = %ld (%.1LfMiB)\n",
      ram_use_excluding_threads, 1.L * ram_use_excluding_threads / (1L << 20));
  long max_block_size = std::max(2L, (long)(ram_use_excluding_threads / 5.2L));
//  if (max_block_size > MAX_32BIT_DIVSUFSORT_LENGTH) {
//    long cur_n_block = (length + max_block_size - 1) / max_block_size;
//    long block_size_2GiB = MAX_32BIT_DIVSUFSORT_LENGTH;
//    long tmp_n_block = (length + block_size_2GiB - 1) / block_size_2GiB;
//    if (tmp_n_block == cur_n_block)
//      max_block_size = block_size_2GiB;
//  }

  fprintf(stderr, "Max block size = %ld (%.1LfMiB)\n\n", max_block_size, 1.L * max_block_size / (1L << 20));
  fprintf(stderr, "Parallel settings:\n");
  fprintf(stderr, "  streaming threads = %ld\n", max_threads);
  fprintf(stderr, "  stream buffer size = %ld\n", stream_buffer_size);
  fprintf(stderr, "  #stream buffers = %ld\n\n", n_stream_buffers);

  long double start = utils::wclock();
  if (max_block_size <= MAX_32BIT_DIVSUFSORT_LENGTH) {
    distributed_file<int> **sparseSA = partial_sufsort<int>(input_filename, length, max_block_size, ram_use, max_threads, stream_buffer_size);
    merge<int>(input_filename, output_filename, length, max_block_size, ram_use, sparseSA);
  } else {
    distributed_file<uint40> **sparseSA = partial_sufsort<uint40>( input_filename, length, max_block_size, ram_use, max_threads, stream_buffer_size);
    merge<uint40>(input_filename, output_filename, length, max_block_size, ram_use, sparseSA);
  }
  long double total_time = utils::wclock() - start;
  fprintf(stderr, "\nTotal time:\n");
  fprintf(stderr, "\tabsolute: %.2Lf\n", total_time);
  fprintf(stderr, "\trelative: %.4Lfs/MiB\n", total_time / ((1.L * length) / (1L << 20)));
  fprintf(stderr, "Speed: %.2LfMiB/s\n", ((1.L * length) / (1L << 20)) / total_time);
}


#endif // __SASCAN_H_INCLUDED
