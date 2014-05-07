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

extern long max_threads;

template<typename block_offset_type>
distributed_file<block_offset_type> **partial_sufsort(std::string filename, long length, long max_block_size, long ram_use);

// Compute SA of <filename> and write to <filename>.sa5 (as normal file).
template<typename output_type>
void SAscan(std::string input_filename, long ram_use) {
  if (ram_use < 5L) {
    fprintf(stderr, "Error: not enough memory to run SAscan.\n");
    std::exit(EXIT_FAILURE);
  }
  
  fprintf(stderr, "Input file = %s\n", input_filename.c_str());
  fprintf(stderr, "RAM use = %ld (%.1LfMiB)\n", ram_use, (long double)ram_use / (1 << 20));
  
  long length = utils::file_size(input_filename);
  fprintf(stderr, "Input length = %ld (%.1LfMiB)\n", length, (long double)length / (1 << 20));
  
  long max_block_size;

  // Compute max_block_size.
  max_block_size = ram_use / 5;
  if (max_block_size > MAX_32BIT_DIVSUFSORT_LENGTH) {
    long cur_n_block = (length + max_block_size - 1) / max_block_size;
    long block_size_2GiB = MAX_32BIT_DIVSUFSORT_LENGTH;
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
    utils::read_objects_from_file<unsigned char>(text, length, input_filename);
    stream_writer<uint40> *sa_writer = new stream_writer<uint40>(input_filename + ".sa5");

    // Use appropriate version (32- or 64-bit) of divsufsort and write SA to disk.
    if (length <= MAX_32BIT_DIVSUFSORT_LENGTH) {
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
    fprintf(stderr, "Using block size = %ld (%.1LfMiB)\n", max_block_size, (long double)max_block_size / (1 << 20));
    fprintf(stderr, "sizeof(output_type) = %ld\n", sizeof(output_type));

    if (max_block_size <= MAX_32BIT_DIVSUFSORT_LENGTH) {
      fprintf(stderr, "sizeof(block_offset_type) = %ld\n", sizeof(int));
      distributed_file<int> **sparseSA = partial_sufsort<int>(input_filename, length, max_block_size, ram_use);
      merge<int>(input_filename, length, max_block_size, ram_use, sparseSA);
    } else {
      distributed_file<uint40> **sparseSA = partial_sufsort<uint40>(input_filename, length, max_block_size, ram_use);
      fprintf(stderr, "sizeof(block_offset_type) = %ld\n", sizeof(uint40));
      merge<uint40>(input_filename, length, max_block_size, ram_use, sparseSA);
    }
  }
  
  long double alg_time_abs = utils::wclock() - alg_start;
  long double alg_time_rel = alg_time_abs / ((1.L * length) / (1 << 20));

  fprintf(stderr, "Total time: %.2Lfs (%.3Lfs/MiB)\n",
      alg_time_abs, alg_time_rel);
}

// Compute SA of <filename> and write to <filename>.sa5 (as distributed file)
// (returned as a result). In addition also compute the BWT.
template<typename output_type>
distributed_file<output_type> *partial_SAscan(
    std::string      input_filename,
    long             ram_use,
    unsigned char  **BWT,
    std::string      text_filename,
    long             text_offset) {
  if (ram_use < 5L) {
    fprintf(stderr, "Error: not enough memory to run SAscan.\n");
    std::exit(EXIT_FAILURE);
  }
  
  fprintf(stderr, "Input file = %s\n", input_filename.c_str());
  fprintf(stderr, "Text file = %s\n", text_filename.c_str());
  fprintf(stderr, "Text offset = %ld\n", text_offset);
  fprintf(stderr, "RAM use = %ld (%.1LfMiB)\n", ram_use, (long double)ram_use / (1 << 20));

  long length = utils::file_size(input_filename);
  fprintf(stderr, "Input length = %ld (%.1LfMiB)\n", length, (long double)length / (1 << 20));

  long max_block_size = ram_use / 9L;

  // Run the algorithm.
  long double alg_start = utils::wclock();

  long n_block = (length + max_block_size - 1) / max_block_size;
  if (n_block == 1) {
    fprintf(stderr, "A single block in partial_SAscan. This should not happen. "
        "Please report this bug (especially filesize and RAM use)!\n");
    std::exit(EXIT_FAILURE);
  }

  fprintf(stderr, "Using block size = %ld (%.1LfMiB)\n", max_block_size, (long double)max_block_size / (1 << 20));
  fprintf(stderr, "sizeof(output_type) = %ld\n", sizeof(output_type));

  distributed_file<output_type> *result = NULL;

  if (max_block_size <= MAX_32BIT_DIVSUFSORT_LENGTH) {
    fprintf(stderr, "sizeof(block_offset_type) = %ld\n", sizeof(int));
    distributed_file<int> **partialSA = partial_sufsort<int>(input_filename, length, max_block_size, ram_use);
    result = partial_merge<int, output_type>(input_filename, length, max_block_size, ram_use, BWT, text_filename, text_offset, partialSA);
  } else {
    fprintf(stderr, "sizeof(block_offset_type) = %ld\n", sizeof(uint40));
    distributed_file<uint40> **partialSA = partial_sufsort<uint40>(input_filename, length, max_block_size, ram_use);
    result = partial_merge<uint40, output_type>(input_filename, length, max_block_size, ram_use, BWT, text_filename, text_offset, partialSA);
  }

  long double alg_time_abs = utils::wclock() - alg_start;
  long double alg_time_rel = alg_time_abs / ((1.L * length) / (1 << 20));

  fprintf(stderr, "Total time: %.2Lfs (%.3Lfs/MiB)\n",
      alg_time_abs, alg_time_rel);

  return result;
}

void SAscan(std::string input_filename, long ram_use) {
  SAscan<uint40>(input_filename, ram_use);
}

#endif // __SASCAN_H_INCLUDED
