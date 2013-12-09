#ifndef __GENERAL_RANK_H
#define __GENERAL_RANK_H

#include <cstdio>
#include <algorithm>

struct general_rank {
  static const int bits = 18;
  static const int block_size = (1 << bits);

  general_rank(unsigned char *text, int length) {
    n = length;
    nblocks = (n + block_size - 1) / block_size;

    size = new int[256];
    std::fill(size, size + 256, 0);
    for (int i = 0; i < n; ++i)
      ++size[text[i]];

    occ = new int*[256];
    for (int i = 0; i < 256; ++i) {
      if (!size[i]) occ[i] = NULL;
      else {
        occ[i] = new int[size[i]];
        if (!occ[i]) {
          fprintf(stderr,"Error: cannot allocate occ[%d].\n", i);
          exit(1);
        }
      }
    }
    std::fill(size, size + 256, 0);
    for (int i = 0; i < n; ++i) {
      unsigned char c = text[i];
      occ[c][size[c]++] = i;
    }
    
    R = new int[(nblocks + 1) * 512];
    if (!R) {
      fprintf(stderr, "Error: cannot allocate R.\n");
      exit(1);
    }

    for (int j = 0; j < 256; ++j) R[j * 2 * (nblocks + 1)] = 0;
    for (int j = 1; j < nblocks; ++j) {
      for (int k = 0; k < 256; ++k)
        R[k * 2 * (nblocks + 1) + 2 * j] = R[k * 2 * (nblocks + 1) + 2 * (j - 1)];
      for (int k = 0; k < std::min(n - (j - 1) * block_size, block_size); ++k)
        ++R[text[(j - 1) * block_size + k] * 2 * (nblocks + 1) + 2 * j];
    }
    for (int j = 0; j < 256; ++j) {
      int *occ_ptr = occ[j], ptr = 0;
      for (int i = 0; i < nblocks; ++i) {
        while (ptr < size[j] && occ_ptr[ptr] < i * block_size)
          ++ptr;
        R[j * 2 * (nblocks + 1) + 2 * i + 1] = ptr;
      }
      R[j * 2 * (nblocks + 1) + 2 * nblocks + 1] = size[j];
    }

    raw_ptr = new int*[256];
    for (int i = 0; i < 256; ++i)
      raw_ptr[i] = occ[i];
  }

  inline int rank(int i, unsigned char c) {
    if (i <= 0) return 0;
    if (i >= n) return size[c];
    int block = (i >> bits);

    int addr = 2 * c * (nblocks + 1) + 2 * block;
    int partial_answer = R[addr];
    int left_ptr = R[addr + 1];
    int right_ptr = R[addr + 3] - 1;
    int old_left_ptr = left_ptr;

    if (left_ptr > right_ptr || occ[c][left_ptr] >= i)
      return partial_answer;
    while (left_ptr != right_ptr) {
      int sr = (left_ptr + right_ptr + 1) >> 1;
      if (occ[c][sr] < i) left_ptr = sr;
      else right_ptr = sr - 1;
    }
    return partial_answer +
      (left_ptr - old_left_ptr + 1);
  }
  
  ~general_rank() {
    delete[] size;
    for (int i = 0; i < 256; ++i)
      if (occ[i]) delete[] occ[i];
    delete[] occ;
    delete[] raw_ptr;
    delete[] R;
  }

  int nblocks;
  int **occ, *size, *R, **raw_ptr, n;
};


#endif // __GENERAL_RANK_H
