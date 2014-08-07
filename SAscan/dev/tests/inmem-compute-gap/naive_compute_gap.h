#ifndef __NAIVE_COMPUTE_GAP_H
#define __NAIVE_COMPUTE_GAP_H

#include <algorithm>

void naive_compute_gap(unsigned char *text, long text_length,
    long left_block_beg, long left_block_size, long right_block_size,
    long* &gap) {
  gap = new long[left_block_size + 1];
  std::fill(gap, gap + left_block_size + 1, 0);

  int *sa = new int[text_length];
  divsufsort(text, sa, (int)text_length);

  long gap_ptr = 0;
  long gap_val = 0;
  long left_block_end = left_block_beg + left_block_size;
  long right_block_beg = left_block_end;
  long right_block_end = right_block_beg + right_block_size;
  for (long i = 0; i < text_length; ++i) {
    if (right_block_beg <= sa[i] && sa[i] < right_block_end) ++gap_val;
    else if (left_block_beg <= sa[i] && sa[i] < left_block_end) {
      gap[gap_ptr++] = gap_val;
      gap_val = 0;
    }
  }
  gap[gap_ptr] = gap_val;

  delete[] sa;
}

#endif  // __NAIVE_COMPUTE_GAP_H

