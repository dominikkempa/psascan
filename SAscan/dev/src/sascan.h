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


//==============================================================================
// Compute SA of <filename> and write to <filename>.sa5
//==============================================================================
template<typename output_type>
void SAscan(std::string input_filename, std::string output_filename, long ram_use) {
  if (ram_use < 5L) {
    fprintf(stderr, "Error: not enough memory to run SAscan.\n");
    std::exit(EXIT_FAILURE);
  }

  // Turn paths absolute.
  input_filename = utils::absolute_path(input_filename);
  output_filename = utils::absolute_path(output_filename);

  fprintf(stderr, "Input file = %s\n", input_filename.c_str());
  fprintf(stderr, "Output file = %s\n", output_filename.c_str());

  long length = utils::file_size(input_filename);
  fprintf(stderr, "Input length = %ld (%.1LfMiB)\n\n", length, (long double)length / (1 << 20));
  fprintf(stderr, "RAM budget = %ld (%.1LfMiB)\n", ram_use, (long double)ram_use / (1 << 20));
  
  long max_block_size;

  // Compute max_block_size.
  max_block_size = ram_use / 5;
  if (max_block_size >= (1L << 31)) {
    long cur_n_block = (length + max_block_size - 1) / max_block_size;
    long block_size_2GiB = (1L << 31) - 1;
    long tmp_n_block = (length + block_size_2GiB - 1) / block_size_2GiB;
    if (tmp_n_block == cur_n_block)
      max_block_size = block_size_2GiB;
  }

  // Run the algorithm.
  long double alg_start = utils::wclock();

  long n_block = (length + max_block_size - 1) / max_block_size;
  if (n_block == 1) { // Only a single block -- run divsufsort.
    // Read text.
    unsigned char *text = NULL;
    utils::read_objects_from_file(text, length, input_filename);
    stream_writer<uint40> *sa_writer = new stream_writer<uint40>(output_filename);

    // Use appropriate version (32- or 64-bit) of divsufsort and write SA to disk.
    if (length < (1L << 31)) {
      int *SA = new int[length];

      fprintf(stderr, "Running divsufsort: ");
      long double sa_start = utils::wclock();
      divsufsort(text, SA, (int)length);
      fprintf(stderr, "%.2Lfs\n", utils::wclock() - sa_start);

      for (long i = 0; i < length; ++i)
        sa_writer->write(SA[i]);
      delete[] SA;
    } else {
      long *SA = new long[length];

      fprintf(stderr, "Running divsufsort64: ");
      long double sa_start = utils::wclock();
      divsufsort64(text, SA, length);
      fprintf(stderr, "%.2Lfs\n", utils::wclock() - sa_start);

      for (long i = 0; i < length; ++i)
        sa_writer->write(SA[i]);
      delete[] SA;
    }

    delete[] text;
    delete sa_writer;
  } else {
    fprintf(stderr, "Max block size = %ld (%.1LfMiB)\n", max_block_size, (long double)max_block_size / (1 << 20));
    fprintf(stderr, "sizeof(output_type) = %ld\n", sizeof(output_type));

    if (max_block_size < (1L << 31)) {
      distributed_file<int> **sparseSA = partial_sufsort<int>(input_filename,
          length, output_filename, max_block_size, ram_use);
      merge<int>(length, output_filename, max_block_size, ram_use, sparseSA);
    } else {
      distributed_file<uint40> **sparseSA = partial_sufsort<uint40>(input_filename,
          length, output_filename, max_block_size, ram_use);
      merge<uint40>(length, output_filename, max_block_size, ram_use, sparseSA);
    }
  }
  
  long double alg_time_abs = utils::wclock() - alg_start;
  long double alg_time_rel = alg_time_abs / ((1.L * length) / (1 << 20));

  fprintf(stderr, "\nTotal time:\n");
  fprintf(stderr, "  absolute: %.2Lf\n", alg_time_abs);
  fprintf(stderr, "  relative: %.4Lfs/MiB\n", alg_time_rel);
  fprintf(stderr, "Speed: %.2LfMiB/s\n", ((1.L * length) / (1L << 20)) / alg_time_abs);
}

void SAscan(std::string input_filename, std::string output_filename, long ram_use) {
  SAscan<uint40>(input_filename, output_filename, ram_use);
}

#endif // __SASCAN_H_INCLUDED
