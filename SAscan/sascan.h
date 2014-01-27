#ifndef __SASCAN_H
#define __SASCAN_H

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
void SAscan_block_size(std::string input_filename, long max_block_size,
    long ram_use, unsigned char **BWT, bool compute_bwt,
    std::string text_filename, long text_offset) {
  long length = utils::file_size(input_filename);
  if (!length) {
    fprintf(stderr, "Error: input is empty.\n");
    std::exit(EXIT_FAILURE);
  }

  fprintf(stderr, "Input length = %ld\n", length);
  fprintf(stderr, "Using block size = %ld\n", max_block_size);
  fprintf(stderr, "sizeof(output_type) = %ld\n", sizeof(output_type));
  fprintf(stderr, "sizeof(block_offset_type) = %ld\n", sizeof(block_offset_type));

  // Run the algorithm.
  long double start = utils::wclock();
  {
    // Compute partial SA and gap arrays.
    partial_sufsort(input_filename, length, max_block_size, ram_use);

    // Merge partial SA arrays if necessary.
    long double n_block = (length + max_block_size - 1) / max_block_size;
    if (n_block > 1) {
      merge<block_offset_type, output_type>(input_filename, length,
          max_block_size, n_block, ram_use, input_filename + ".sa5",
          BWT, compute_bwt, text_filename, text_offset);
    } else {
      // Invariant: compute_bwt == false (we only need BWT from recursive SAscan
      // calls and these always require merging, i.e., n_block > 1).
      utils::execute("mv " + input_filename + ".partial_sa.0 " + input_filename + ".sa5");
    }
  }
  long double total_time = utils::wclock() - start;
  long double speed = total_time / ((1.L * length) / (1 << 20));

  fprintf(stderr, "Total time: %.2Lfs. Speed: %.2Lfs/MiB\n",
      total_time, speed);
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
  fprintf(stderr, "RAM use = %ld\n", ram_use);
  fprintf(stderr, "in-recursion = %s\n", in_recursion ? "TRUE" : "FALSE");

  long max_block_size;
  if (in_recursion) max_block_size = ram_use / 9;
  else {
    max_block_size = ram_use / 5;
    if (max_block_size > MAX_32BIT_DIVSUFSORT_LENGTH) {
      long length = utils::file_size(input_filename);
      long n_block = (length + max_block_size - 1) / max_block_size;
      long block_size_2GiB = MAX_32BIT_DIVSUFSORT_LENGTH;
      long tmp_n_block = (length + block_size_2GiB - 1) / block_size_2GiB;
      if (tmp_n_block == n_block) max_block_size = block_size_2GiB;
    }
  }

  // Currently we use simple criterion: for recursive calls (compute_bwt == true)
  // we know that it was true that ram_use == 5 * length thus if we set
  // max_block_size := ram_use / 9 there will be (I skip the exact inequalities)
  //   ceil( length / max_block_size )     =
  //   ceil( (length * 9) / ram_use )      =
  //   ceil( (length * 9) / (length * 5) ) =
  //   ceil( 9 / 5 )                       = 2
  // blocks and no further recursion.

  // TODO: test the option of setting the block size to 2^31-1, (even if
  // ram_use/5 is bigger). This is obviously beneficial if results in the same
  // number of blocks, e.g., in the second level of recursion when e.g. length
  // = 3.8GiB and ram_use = 15GiB. Currently, we would set block_size to 3GiB
  // and end up using recursion for the leftmost block. *But*, we could just
  // use block size = 1.8GiB. It also results in two blocks, but avoids the
  // recursion which reduces the I/O and thus also time.

  // The I/O saving gained by avoiding recursion could also save enough time
  // that using block size == 2^31 - 1 could be even beneficial in the top level
  // recursive calls. We could always compute the number of blocks when setting
  // block_size to ram_use / 5 (n1) and to 2^31-1 (n2) and compare them. Some
  // heuristic here could decide at what ratio of n1 and n2 is it good to go
  // with smaller blocks and sacrifice a little streaming time for the reduced
  // I/O.

  if (max_block_size <= MAX_32BIT_DIVSUFSORT_LENGTH)
    SAscan_block_size<int, output_type>(input_filename, max_block_size,
      ram_use, BWT, compute_bwt, text_filename, text_offset);
  else SAscan_block_size<uint40, output_type>(input_filename, max_block_size,
      ram_use, BWT, compute_bwt, text_filename, text_offset);
}

void SAscan(std::string filename, long ram_use) {
  SAscan<uint40>(filename, ram_use, NULL, false, false, "", 0);
}

#endif // __SASCAN_H
