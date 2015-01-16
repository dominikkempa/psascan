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
void test(unsigned char *supertext, long supertext_length) {
  int *supertext_sa = new int[supertext_length];
  saisxx(supertext, supertext_sa, (int)supertext_length);
  long combinations = 0L;

#if 0
  fprintf(stderr, "\tsupertext = ");
  for (long j = 0; j < supertext_length; ++j)
    fprintf(stderr, "%c", supertext[j]);
  fprintf(stderr, "\n");
#endif

  for (long text_length = 1; text_length <= supertext_length; ++text_length) {
    bwtsa_t<saidx_t> *bwtsa = new bwtsa_t<saidx_t>[text_length];

    for (long text_beg = 0; text_beg <= supertext_length - text_length; ++text_beg) {
      for (long max_block_size = 1; max_block_size <= text_length; ++max_block_size) {
        long text_end = text_beg + text_length;
        long tail_length = supertext_length - text_end;
        long n_blocks = (text_length + max_block_size - 1) / max_block_size;
        ++combinations;

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
      }
    }
    delete[] bwtsa;
  }
  delete[] supertext_sa;
//  fprintf(stderr, "combinations = %ld\n", combinations);
}


template<typename saidx_t>
void test_all(long max_supertext_length) {
  for (long supertext_length = 1; supertext_length <= max_supertext_length; ++supertext_length) {
    fprintf(stderr,"TEST, sizeof(saidx_t) = %ld, length = %ld\r", (long)sizeof(saidx_t), supertext_length);
    unsigned char *supertext = new unsigned char[supertext_length];
    long testcases = (1L << supertext_length);
    for (long supertext_enc = 0L, dbg = 0L; supertext_enc < testcases; ++supertext_enc, ++dbg) {
      if (dbg == 100) {
        fprintf(stderr,"TEST, sizeof(saidx_t) = %ld, length = %ld: %ld (%.0Lf%%)\r", (long)sizeof(saidx_t),
            supertext_length, supertext_enc + 1, ((supertext_enc + 1) * 100.L) / testcases);
        dbg = 0;
      }

      for (long j = 0; j < supertext_length; ++j)
        if (supertext_enc & (1L << j)) supertext[j] = 'b';
        else supertext[j] = 'a';

      test<saidx_t>(supertext, supertext_length);
    }
    delete[] supertext;
    fprintf(stderr,"TEST, sizeof(saidx_t) = %ld, length = %ld: \033[22;32mPASSED\033[0m%10s\n",
        (long)sizeof(saidx_t), supertext_length, "");
  }
}

int main(int, char **) {
  std::srand(std::time(0) + getpid());

  test_all<int>   (17);
  test_all<uint40>(17);

  fprintf(stderr, "All tests passed.\n");
}

