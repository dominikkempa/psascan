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
    long stream_buffer_size, const char *sa_filename) {

  long double start;

  /*fprintf(stderr, "  Compute gap naively: ");
  start = utils::wclock();
  long *correct_gap = NULL;
  naive_compute_gap(text, text_length, left_block_beg, left_block_size,
      right_block_size, correct_gap);
  fprintf(stderr, "\n");*/

  // Prerequisite: gt bitvector.
  fprintf(stderr, "  Computing gt: ");
  start = utils::wclock();
  bitvector *gt_in = NULL;
  long left_block_end = left_block_beg + left_block_size;
  naive_compute_gt(text, text_length, left_block_end, right_block_size + 1, gt_in);
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);

  fprintf(stderr, "  Computing partial suffix array:\n");
  start = utils::wclock();
  // Prerequisite: partial suffix array.
  int *partial_sa = NULL;
  naive_select_partial_sa(text_length, left_block_beg, left_block_end,
      std::ref(partial_sa), sa_filename);
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

  /*fprintf(stderr, "  Computing gap plain: ");
  start = utils::wclock();
  long *computed_gap_plain = new long[left_block_size + 1];
  long excess_ptr = 0L;
  for (long i = 0; i < left_block_size + 1; ++i) {
    // Compute gap[i] from the small-gap-array representation.
    long gap_i = computed_gap->m_count[i];
    while (excess_ptr < (long)computed_gap->m_excess.size() &&
        computed_gap->m_excess[excess_ptr] == i) {
      gap_i += 256;
      ++excess_ptr;
    }
    computed_gap_plain[i] = gap_i;
  }

  if (!std::equal(correct_gap, correct_gap + left_block_size + 1,
        computed_gap_plain)) {
    fprintf(stderr, "FAIL\n");
  } else fprintf(stderr, "OK\n");*/

  // Clean up.
  //delete[] computed_gap_plain;
  //delete[] correct_gap;

  delete gt_in;
  delete gt_out;
  delete[] partial_sa;
  delete computed_gap;
}


void test_file(const char *filename) {
  fprintf(stderr, "Reading text: ");
  long length;
  unsigned char *text;
  utils::read_objects_from_file(text, length, filename);
  long left_block_size = length / 2;
  long right_block_size = length - left_block_size;
  long left_block_beg = 0;
  long max_threads = 24;
  long stream_buffer_size = (1L << 20);
  fprintf(stderr, "DONE\n");

  std::string sa_filename = std::string(filename) + ".sa";
  test(text, length, left_block_beg, left_block_size,
      right_block_size, max_threads, stream_buffer_size, sa_filename.c_str());

  delete[] text;
}

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s <file>\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }
  test_file(argv[1]);
}
