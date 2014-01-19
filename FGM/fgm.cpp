#include <cstdio>
#include <cstdlib>
#include <string>
#include <algorithm>

#include "fgm.h"
#include "partial_sufsort.h"
#include "merge.h"
#include "utils.h"

void FGM(std::string filename, long ram_use) {
  fprintf(stderr, "RAM use = %ld\n", ram_use);

  long max_block_size = (ram_use + 4) / 5;
  FGM_block_size(filename, max_block_size);
}

void FGM_block_size(std::string filename, long max_block_size) {
  long double start = utils::wclock();
  long length = utils::file_size(filename);
  
  if (!length) {
    fprintf(stderr, "Error: input text is empty.\n");
    std::exit(EXIT_FAILURE);
  }

  fprintf(stderr, "Input file = %s\n", filename.c_str());
  fprintf(stderr, "Input length = %ld\n", length);
  fprintf(stderr, "Using block size = %ld\n", max_block_size);

  partial_sufsort(filename, length, max_block_size);
  merge(filename, length, max_block_size, filename + ".sa5");
  
  long double total_time = utils::wclock() - start;
  long double speed = total_time / ((1.L * length) / (1 << 20));
  fprintf(stderr, "Total time: %.2Lfs. Speed: %.2Lfs/MiB\n",
      total_time, speed);
}

