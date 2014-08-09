#include <cstdio>
#include <cstring>
#include <ctime>
#include <unistd.h>

#include "utils.h"
#include "bitvector.h"
#include "inmem_compute_gap.h"
#include "inmem_gap_array.h"
#include "naive_compute_gt.h"
#include "naive_compute_partial_sa.h"
#include "naive_compute_gap.h"


void test(unsigned char *text, long text_length, long left_block_beg,
    long left_block_size, long right_block_size, long max_threads,
    long stream_buffer_size) {

  // Prerequisite: gt bitvector.
  fprintf(stderr, "  Computing gt: ");
  long double start = utils::wclock();
  bitvector *gt_in = NULL;
  long left_block_end = left_block_beg + left_block_size;
  naive_compute_gt(text, text_length, left_block_end, right_block_size + 1, gt_in);
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

  fprintf(stderr, "  Computing partial suffix array:\n");
  start = utils::wclock();
  // Prerequisite: partial suffix array.
  int *partial_sa = NULL;
  naive_compute_partial_sa(text, text_length, left_block_beg, left_block_end,
      std::ref(partial_sa));
  fprintf(stderr, "  Total time: %.2Lf\n", utils::wclock() - start);

  // Compute the gap array, this is the method we are testing.
  fprintf(stderr, "  Computing gap array:\n");
  start = utils::wclock();
  bitvector *gt_out = NULL;
  inmem_gap_array *computed_gap = NULL;
  inmem_compute_gap(text, text_length, left_block_beg, left_block_size,
      right_block_size, partial_sa, gt_in, gt_out, true, computed_gap,
      max_threads, stream_buffer_size);
  fprintf(stderr, "  Total: %.2Lf\n", utils::wclock() - start);

  // Clean up.
  delete gt_in;
  delete gt_out;
  delete[] partial_sa;
  delete computed_gap;
}


void test_random(long length) {
  fprintf(stderr, "Generating text: ");
  unsigned char *text = new unsigned char[length];
  long left_block_size = length / 2;
  long right_block_size = length - left_block_size;
  long left_block_beg = 0;
  long max_threads = 24;
  long stream_buffer_size = (1L << 20);
  utils::fill_random_string(text, length, 64);
  fprintf(stderr, "DONE\n");

  test(text, length, left_block_beg, left_block_size,
      right_block_size, max_threads, stream_buffer_size);

  delete[] text;
}

int main() {
  test_random(1500L << 20);
}
