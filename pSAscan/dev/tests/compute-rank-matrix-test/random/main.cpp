#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <string>

#include "bwtsa.h"
#include "sais.hxx"
#include "utils.h"
#include "bitvector.h"
#include "inmem_compute_block_rank_matrix.h"

using namespace inmem_sascan_private;


void compute_gt_begin_reversed(unsigned char *text, long text_length, bitvector *gt_begin_reversed) {
  for (long j = 1; j < text_length; ++j) {
    long lcp = 0L;
    while (j + lcp < text_length && text[lcp] == text[j + lcp]) ++lcp;
    bool gt_j = (j + lcp < text_length && text[lcp] < text[j + lcp]);
    if (gt_j) gt_begin_reversed->set(text_length - j);
  }
}

template<typename saidx_t>
void test(unsigned char *supertext, long supertext_length, long text_beg, long text_length, long max_block_size) {
  int *supertext_sa = new int[supertext_length];
  saisxx(supertext, supertext_sa, (int)supertext_length);

  bwtsa_t<saidx_t> *bwtsa = new bwtsa_t<saidx_t>[text_length];
  long text_end = text_beg + text_length;
  long tail_length = supertext_length - text_end;
  long n_blocks = (text_length + max_block_size - 1) / max_block_size;

  // Compute reversed gt_begin for the tail (tail_gt_begin_reversed);
  bitvector *tail_gt_begin_reversed = new bitvector(tail_length);
  compute_gt_begin_reversed(supertext + text_end, tail_length, tail_gt_begin_reversed);

  // Compute bwtsa.
  for (long block_id = 0; block_id < n_blocks; ++block_id) {
    long block_end = text_length - (n_blocks - 1 - block_id) * max_block_size;  // local in the text
    long block_beg = std::max(0L, block_end - max_block_size);  // local in the text
    for (long j = 0, ptr = block_beg; j < supertext_length; ++j)
      if (block_beg <= (supertext_sa[j] - text_beg) && (supertext_sa[j] - text_beg) < block_end)
        bwtsa[ptr++].sa = (saidx_t)(supertext_sa[j] - text_beg - block_beg);  // local in the block
  }

  // Run the tested algorithm.
  unsigned char *text = supertext + text_beg;
  unsigned char *next_block = supertext + text_end;
  long **computed_matrix = new long*[n_blocks];
  for (long row = 0; row < n_blocks; ++row) computed_matrix[row] = new long[n_blocks];
  compute_block_rank_matrix<saidx_t>(text, text_length, bwtsa, max_block_size,
      text_beg, supertext_length, tail_gt_begin_reversed, NULL, next_block, computed_matrix);

  // Compute the answer by brute force.
  long **correct_matrix = new long*[n_blocks];
  for (long row = 0; row < n_blocks; ++row) correct_matrix[row] = new long[n_blocks];
  for (long row = 0; row < n_blocks; ++row) {
    for (long col = row + 1; col < n_blocks; ++col) {
      long col_block_end = text_length - (n_blocks - 1 - col) * max_block_size;
      long pattern = text_beg + col_block_end;
      long block_end = text_length - (n_blocks - 1 - row) * max_block_size;
      long block_beg = std::max(0L, block_end - max_block_size);
      block_beg += text_beg;
      block_end += text_beg;
      long ans = 0;
      for (long j = 0; j < supertext_length; ++j) {
        if (supertext_sa[j] == pattern || pattern == supertext_length) break;
        if (block_beg <= supertext_sa[j] && supertext_sa[j] < block_end)
          ++ans;
      }
      correct_matrix[row][col] = ans;
    }
  }

  bool equal = true;
  for (long row = 0; row < n_blocks; ++row)
    for (long col = row + 1; col < n_blocks; ++col)
      if (computed_matrix[row][col] != correct_matrix[row][col]) equal = false;

  if (!equal) {
    fprintf(stderr, "\nError:\n");
    fprintf(stderr, "\tsupertext = ");
    for (long j = 0; j < supertext_length; ++j)
      fprintf(stderr, "%c", supertext[j]);
    fprintf(stderr, "\n");
    fprintf(stderr, "\ttext_beg = %ld\n", text_beg);
    fprintf(stderr, "\ttext_length = %ld\n", text_length);
    fprintf(stderr, "\ttail_length = %ld\n", tail_length);
    fprintf(stderr, "\ttail_gt_begin = ");
    for (long j = 0; j < tail_length; ++j)
      fprintf(stderr, "%d ", tail_gt_begin_reversed->get(j));
    fprintf(stderr, "\n");
    fprintf(stderr, "\tsa: ");
    for (long j = 0; j < text_length; ++j)
      fprintf(stderr, "%ld ", (long)bwtsa[j].sa);
    fprintf(stderr, "\n");
    fprintf(stderr, "\tmax_block_size = %ld\n", max_block_size);
    for (long row = 0; row < n_blocks; ++row) {
      for (long col = row + 1; col < n_blocks; ++col) {
        if (computed_matrix[row][col] != correct_matrix[row][col]) {
          fprintf(stderr, "\t\tcomputed[%ld][%ld] = %ld, correct[%ld][%ld] = %ld\n",
              row, col, computed_matrix[row][col], row, col, correct_matrix[row][col]);
        }
      }
    }
    std::exit(EXIT_FAILURE);
  }

  // Clean up.
   delete tail_gt_begin_reversed;
  for (long row = 0; row < n_blocks; ++row) {
    delete[] correct_matrix[row];
    delete[] computed_matrix[row];
  }
  delete[] correct_matrix;
  delete[] computed_matrix;
  delete[] bwtsa;
  delete[] supertext_sa;
}


template<typename saidx_t>
void test_random(int testcases, int max_length) {
  fprintf(stderr,"TEST, testcases = %d, max_n = %d\r", testcases, max_length);
  unsigned char *supertext = new unsigned char[max_length + 1];

  for (int tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {
    // Print progress information.
    if (dbg == 1000) {
      fprintf(stderr,"TEST, testcases = %d, max_n = %d: ""%d (%.0Lf%%)\r",
          testcases, max_length, tc, (tc * 100.L) / testcases);
      dbg = 0;
    }

    long supertext_length = utils::random_long(1L, max_length);
    long text_length = utils::random_long(1L, supertext_length);
    long text_beg = utils::random_long(0L, supertext_length - text_length);
    utils::fill_random_letters(supertext, supertext_length, 2);

    for (long max_block_size = 1; max_block_size <= text_length; ++max_block_size)
      test<saidx_t>(supertext, supertext_length, text_beg, text_length, max_block_size);
  }

  // Clean up.
  delete[] supertext;

  fprintf(stderr,"TEST, testcases = %d, max_n = %d: "
      "\033[22;32mPASSED\033[0m%10s\n", testcases, max_length, "");
}

int main(int, char **) {
  std::srand(std::time(0) + getpid());

  test_random<int>   (2000000,  20);
  test_random<uint40>(2000000,  20);

  fprintf(stderr, "All tests passed.\n");
}

