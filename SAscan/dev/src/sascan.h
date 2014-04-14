#ifndef __SASCAN_H_INCLUDED
#define __SASCAN_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <string>
#include <algorithm>

#include "sascan.h"
#include "partial_sufsort.h"
#include "merge.h"
#include "utils.h"
#include "uint40.h"
#include "settings.h"

void partial_sufsort(std::string filename, long length, long max_block_size, long ram_use);

// Compute SA of <filename> and write to <filename>.sa5 using given
// block size. Optionally also compute the BWT during merging.
template<typename block_offset_type, typename output_type>
void SAscan_block_size(std::string input_filename, long length, long max_block_size,
    long ram_use, unsigned char **BWT, bool compute_bwt,
    std::string text_filename, long text_offset) {
  fprintf(stderr, "Using block size = %ld\n", max_block_size);
  fprintf(stderr, "sizeof(output_type) = %ld\n", sizeof(output_type));
  fprintf(stderr, "sizeof(block_offset_type) = %ld\n", sizeof(block_offset_type));

  partial_sufsort(input_filename, length, max_block_size, ram_use);
  merge<block_offset_type, output_type>(input_filename,
      length,
      max_block_size,
      ram_use,
      BWT,
      compute_bwt,
      text_filename,
      text_offset);
}

// Compute SA of <filename> and write to <filename>.sa5 using given
// RAM limit. Optionally also compute the BWT.
template<typename output_type>
void SAscan(std::string      input_filename,
            long             ram_use,
            unsigned char ** BWT,
            bool             in_recursion,
            bool             compute_bwt,
            std::string      text_filename,
            long             text_offset) {
  if (ram_use < 5L) {
    fprintf(stderr, "Error: not enough memory to run SAscan.\n");
    std::exit(EXIT_FAILURE);
  }
  
  fprintf(stderr, "Input file = %s\n", input_filename.c_str());
  if (compute_bwt) {
    fprintf(stderr, "Text file = %s\n", text_filename.c_str());
    fprintf(stderr, "Text offset = %ld\n", text_offset);
  }
  fprintf(stderr, "RAM use = %ld (%.1LfMiB)\n",
      ram_use, (long double)ram_use / (1 << 20));
  fprintf(stderr, "in-recursion = %s\n", in_recursion ? "TRUE" : "FALSE");
  
  long length = utils::file_size(input_filename);
  fprintf(stderr, "Input length = %ld (%.2LfMiB)\n",
      length, (long double)length / (1 << 20));
  
  long max_block_size;

  // Compute max_block_size.
  if (in_recursion) max_block_size = ram_use / 9;
  else {
    max_block_size = ram_use / 5;
    if (max_block_size > MAX_32BIT_DIVSUFSORT_LENGTH) {
      long cur_n_block = (length + max_block_size - 1) / max_block_size;
      long block_size_2GiB = MAX_32BIT_DIVSUFSORT_LENGTH;
      long tmp_n_block = (length + block_size_2GiB - 1) / block_size_2GiB;
      if (tmp_n_block == cur_n_block)
        max_block_size = block_size_2GiB;
    }
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
    if (max_block_size <= MAX_32BIT_DIVSUFSORT_LENGTH)
      SAscan_block_size<int, output_type>(input_filename, length, max_block_size,
        ram_use, BWT, compute_bwt, text_filename, text_offset);
    else SAscan_block_size<uint40, output_type>(input_filename, length, max_block_size,
        ram_use, BWT, compute_bwt, text_filename, text_offset);
  }
  
  long double alg_time_abs = utils::wclock() - alg_start;
  long double alg_time_rel = alg_time_abs / ((1.L * length) / (1 << 20));

  fprintf(stderr, "Total time: %.2Lfs (%.2Lfs/MiB)\n",
      alg_time_abs, alg_time_rel);
}

void SAscan(std::string filename, long ram_use) {
  SAscan<uint40>(filename, ram_use, NULL, false, false, "", 0);
}

#endif // __SASCAN_H_INCLUDED
