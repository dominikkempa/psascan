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

extern long n_streamers;
extern long stream_buffer_size;
extern long n_stream_buffers;

const long kMiB = (1L << 20);

template<typename block_offset_type>
distributed_file<block_offset_type> **partial_sufsort(std::string filename,
    long length, long max_block_size, long ram_use);

// Compute SA of input_filename and write to output_filename (as normal file).
template<typename output_type>
void SAscan(std::string input_filename, std::string output_filename, long ram_use) {
  if (ram_use < 6L) {
    fprintf(stderr, "Error: not enough memory to run SAscan.\n");
    std::exit(EXIT_FAILURE);
  }
  
  // Turn paths absolute.
  input_filename = utils::absolute_path(input_filename);
  output_filename = utils::absolute_path(output_filename);
  fprintf(stderr, "Input filename = %s\n", input_filename.c_str());
  fprintf(stderr, "Output filename = %s\n", output_filename.c_str());

  long length = utils::file_size(input_filename);
  fprintf(stderr, "Input length = %ld (%.1LfMiB)\n", length, 1.L * length / kMiB);
  fprintf(stderr, "RAM use = %ld (%.1LfMiB)\n", ram_use, 1.L * ram_use / kMiB);

  // Start the time.
  long double alg_start = utils::wclock();

  // Check if running divsufsort is enough.
  if ((length <= MAX_32BIT_DIVSUFSORT_LENGTH && length * 5L <= ram_use) ||
      (length > MAX_32BIT_DIVSUFSORT_LENGTH && length * 9L <= ram_use)) {
    // Read text.
    unsigned char *text = NULL;
    utils::read_objects_from_file<unsigned char>(text, length, input_filename);
    stream_writer<uint40> *sa_writer = new stream_writer<uint40>(output_filename);

    // Use appropriate (32- or 64-bit) version of divsufsort and write SA to disk.
    if (length <= MAX_32BIT_DIVSUFSORT_LENGTH) {
      int *SA = new int[length];

      fprintf(stderr, "Running divsufsort: ");
      long double sa_start = utils::wclock();
      divsufsort(text, SA, (int)length);
      fprintf(stderr, "%.2Lfs\n", utils::wclock() - sa_start);

      for (long i = 0L; i < length; ++i)
        sa_writer->write(SA[i]);
      delete[] SA;
    } else {
      long *SA = new long[length];

      fprintf(stderr, "Running divsufsort64: ");
      long double sa_start = utils::wclock();
      divsufsort64(text, SA, length);
      fprintf(stderr, "%.2Lfs\n", utils::wclock() - sa_start);

      for (long i = 0L; i < length; ++i)
        sa_writer->write(SA[i]);
      delete[] SA;
    }

    delete[] text;
    delete sa_writer;
  } else { // We need to use SAscan.

    fprintf(stderr, "Parallel settings:\n");
    fprintf(stderr, "  Streaming threads = %ld\n", n_streamers);
    fprintf(stderr, "  Stream-buffer size = %ld (%.1LfMiB)\n",
        stream_buffer_size, 1.L * stream_buffer_size / kMiB);
    fprintf(stderr, "  Number of stream-buffers = %ld\n", n_stream_buffers);

    long ram_use_excluding_threads = ram_use - n_stream_buffers * stream_buffer_size;
    if (ram_use_excluding_threads < 6L) {
      long required_bytes = n_stream_buffers * stream_buffer_size;
      long required_MiB = (required_bytes + kMiB - 1) / kMiB;
      fprintf(stderr, "Error: not enough memory to start threads. You need "
          "at least %ldMiB\n", required_MiB + 1);
      std::exit(EXIT_FAILURE);
    }
    
    // Compute max_block_size.
    fprintf(stderr, "RAM use (excluding threads) = %ld (%.1LfMiB)\n",
        ram_use_excluding_threads, 1.L * ram_use_excluding_threads / kMiB);
    long max_block_size = ram_use_excluding_threads / 5.2L;
    if (max_block_size > MAX_32BIT_DIVSUFSORT_LENGTH) {
      long cur_n_block = (length + max_block_size - 1) / max_block_size;
      long block_size_2GiB = MAX_32BIT_DIVSUFSORT_LENGTH;
      long tmp_n_block = (length + block_size_2GiB - 1) / block_size_2GiB;
      if (tmp_n_block == cur_n_block)
        max_block_size = block_size_2GiB;
    }
 
    fprintf(stderr, "Max block size = %ld (%.1LfMiB)\n",
        max_block_size, 1.L * max_block_size / kMiB);
    fprintf(stderr, "sizeof(output_type) = %ld\n", sizeof(output_type));

    if (max_block_size <= MAX_32BIT_DIVSUFSORT_LENGTH) {
      fprintf(stderr, "sizeof(block_offset_type) = %ld\n", sizeof(int));
      distributed_file<int> **sparseSA =
        partial_sufsort<int>(input_filename, length, max_block_size, ram_use);
      merge<int>(input_filename, output_filename, length, max_block_size,
          ram_use, sparseSA);
    } else {
      distributed_file<uint40> **sparseSA =
        partial_sufsort<uint40>(input_filename, length, max_block_size, ram_use);
      fprintf(stderr, "sizeof(block_offset_type) = %ld\n", sizeof(uint40));
      merge<uint40>(input_filename, output_filename, length, max_block_size,
          ram_use, sparseSA
      );
    }
  }
  
  long double alg_time_abs = utils::wclock() - alg_start;
  long double alg_time_rel = alg_time_abs / ((1.L * length) / kMiB);

  fprintf(stderr, "Total time: %.2Lfs (%.3Lfs/MiB)\n",
      alg_time_abs, alg_time_rel);
}

// Compute SA of input_filename and write to input_filename.sa5 (as
// distributed file) (returned as a result). If compute_bwt == true
// also compute the BWT associated with the SA.
template<typename output_type>
distributed_file<output_type> *partial_SAscan(
    std::string      input_filename,
    bool             compute_bwt,
    long             ram_use,
    unsigned char**  BWT,
    std::string      text_filename,
    long             text_offset) {
  if (ram_use < 6L) {
    fprintf(stderr, "Error: not enough memory to run SAscan.\n");
    std::exit(EXIT_FAILURE);
  }
  
  std::string output_filename = input_filename + ".sa5";

  fprintf(stderr, "Input filename = %s\n", input_filename.c_str());
  fprintf(stderr, "Output filename = %s\n", output_filename.c_str());

  long length = utils::file_size(input_filename);
  fprintf(stderr, "Input length = %ld (%.1LfMiB)\n", length, 1.L * length / kMiB);

  fprintf(stderr, "Text file = %s\n", text_filename.c_str());
  fprintf(stderr, "Text offset = %ld\n", text_offset);
  
  fprintf(stderr, "RAM use = %ld (%.1LfMiB)\n", ram_use, 1.L * ram_use / kMiB);
  long ram_use_excluding_threads = ram_use - n_stream_buffers * stream_buffer_size;
    fprintf(stderr, "RAM use (excluding threads) = %ld (%.1LfMiB)\n",
        ram_use_excluding_threads, 1.L * ram_use_excluding_threads / kMiB);

  long max_block_size = ram_use_excluding_threads / 9L;

  // Run the algorithm.
  long double alg_start = utils::wclock();

  long n_block = (length + max_block_size - 1) / max_block_size;
  if (n_block == 1) {
    fprintf(stderr, "A single block in partial_SAscan. This should not happen. "
        "Please report this bug (especially filesize and RAM use)!\n");
    std::exit(EXIT_FAILURE);
  }

  fprintf(stderr, "Max block size = %ld (%.1LfMiB)\n",
      max_block_size, 1.L * max_block_size / kMiB);
  fprintf(stderr, "sizeof(output_type) = %ld\n", sizeof(output_type));

  distributed_file<output_type> *result = NULL;

  if (max_block_size <= MAX_32BIT_DIVSUFSORT_LENGTH) {
    fprintf(stderr, "sizeof(block_offset_type) = %ld\n", sizeof(int));
    distributed_file<int> **partialSA =
      partial_sufsort<int>(input_filename, length, max_block_size, ram_use);
    result = partial_merge<int, output_type>(input_filename,
        output_filename, length, compute_bwt, max_block_size,
        ram_use, BWT, text_filename, text_offset, partialSA);
  } else {
    fprintf(stderr, "sizeof(block_offset_type) = %ld\n", sizeof(uint40));
    distributed_file<uint40> **partialSA =
      partial_sufsort<uint40>(input_filename, length, max_block_size, ram_use);
    result = partial_merge<uint40, output_type>(input_filename,
        output_filename, length, compute_bwt, max_block_size,
        ram_use, BWT, text_filename, text_offset, partialSA);
  }

  long double alg_time_abs = utils::wclock() - alg_start;
  long double alg_time_rel = alg_time_abs / ((1.L * length) / kMiB);

  fprintf(stderr, "Total time: %.2Lfs (%.3Lfs/MiB)\n",
      alg_time_abs, alg_time_rel);

  return result;
}

void SAscan(std::string input_filename, std::string output_filename, long ram_use) {
  SAscan<uint40>(input_filename, output_filename, ram_use);
}

#endif // __SASCAN_H_INCLUDED
