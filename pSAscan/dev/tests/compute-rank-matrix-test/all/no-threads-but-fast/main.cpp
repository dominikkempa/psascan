#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <string>

#include "bwtsa.hpp"
#include "sais.hxx"
#include "bitvector.hpp"
#include "inmem_compute_block_rank_matrix.hpp"

using namespace inmem_psascan_private;


void compute_gt_begin_reversed(std::uint8_t *text, std::uint64_t text_length, bitvector *gt_begin_reversed) {
  for (std::uint64_t j = 1; j < text_length; ++j) {
    std::uint64_t lcp = 0L;
    while (j + lcp < text_length && text[lcp] == text[j + lcp]) ++lcp;
    bool gt_j = (j + lcp < text_length && text[lcp] < text[j + lcp]);
    if (gt_j) gt_begin_reversed->set(text_length - j);
  }
}

template<typename saidx_t>
void test(std::uint8_t *supertext, std::uint64_t supertext_length) {
  int *supertext_sa = new int[supertext_length];
  saisxx(supertext, supertext_sa, (int)supertext_length);
  std::uint64_t combinations = 0L;

#if 0
  fprintf(stderr, "\tsupertext = ");
  for (std::uint64_t j = 0; j < supertext_length; ++j)
    fprintf(stderr, "%c", supertext[j]);
  fprintf(stderr, "\n");
#endif

  for (std::uint64_t text_length = 1; text_length <= supertext_length; ++text_length) {
    bwtsa_t<saidx_t> *bwtsa = new bwtsa_t<saidx_t>[text_length];

    for (std::uint64_t text_beg = 0; text_beg <= supertext_length - text_length; ++text_beg) {
      for (std::uint64_t max_block_size = 1; max_block_size <= text_length; ++max_block_size) {
        std::uint64_t text_end = text_beg + text_length;
        std::uint64_t tail_length = supertext_length - text_end;
        std::uint64_t n_blocks = (text_length + max_block_size - 1) / max_block_size;
        ++combinations;

        // Compute reversed gt_begin for the tail (tail_gt_begin_reversed);
        bitvector *tail_gt_begin_reversed = new bitvector(tail_length);
        compute_gt_begin_reversed(supertext + text_end, tail_length, tail_gt_begin_reversed);

        // Compute bwtsa.
        for (std::uint64_t block_id = 0; block_id < n_blocks; ++block_id) {
          std::uint64_t block_end = text_length - (n_blocks - 1 - block_id) * max_block_size;  // local in the text
          std::uint64_t block_beg = std::max((std::int64_t)0,
              (std::int64_t)block_end - (std::int64_t)max_block_size);  // local in the text
          for (std::uint64_t j = 0, ptr = block_beg; j < supertext_length; ++j)
            if (block_beg <= (std::uint64_t)(supertext_sa[j] - text_beg) &&
                (std::uint64_t)(supertext_sa[j] - text_beg) < block_end)
              bwtsa[ptr++].m_sa = (saidx_t)(supertext_sa[j] - text_beg - block_beg);  // local in the block
        }

        // Run the tested algorithm.
        const std::uint8_t *text = supertext + text_beg;
        const std::uint8_t *next_block = supertext + text_end;
        std::uint64_t **computed_matrix = new std::uint64_t*[n_blocks];
        for (std::uint64_t row = 0; row < n_blocks; ++row)
          computed_matrix[row] = new std::uint64_t[n_blocks];
        compute_block_rank_matrix<saidx_t>(text, text_length, bwtsa,
            max_block_size, text_beg, supertext_length,
            tail_gt_begin_reversed, NULL, next_block, computed_matrix);

        // Compute the answer by brute force.
        std::uint64_t **correct_matrix = new std::uint64_t*[n_blocks];
        for (std::uint64_t row = 0; row < n_blocks; ++row)
          correct_matrix[row] = new std::uint64_t[n_blocks];
        for (std::uint64_t row = 0; row < n_blocks; ++row) {
          for (std::uint64_t col = row + 1; col < n_blocks; ++col) {
            std::uint64_t col_block_end = text_length - (n_blocks - 1 - col) * max_block_size;
            std::uint64_t pattern = text_beg + col_block_end;
            std::uint64_t block_end = text_length - (n_blocks - 1 - row) * max_block_size;
            std::uint64_t block_beg = std::max((std::int64_t)0,
                (std::int64_t)block_end - (std::int64_t)max_block_size);
            block_beg += text_beg;
            block_end += text_beg;
            std::uint64_t ans = 0;
            for (std::uint64_t j = 0; j < supertext_length; ++j) {
              if ((std::uint64_t)supertext_sa[j] == pattern ||
                  pattern == supertext_length) break;
              if (block_beg <= (std::uint64_t)supertext_sa[j] &&
                  (std::uint64_t)supertext_sa[j] < block_end)
                ++ans;
            }
            correct_matrix[row][col] = ans;
          }
        }

        bool equal = true;
        for (std::uint64_t row = 0; row < n_blocks; ++row)
          for (std::uint64_t col = row + 1; col < n_blocks; ++col)
            if (computed_matrix[row][col] != correct_matrix[row][col]) equal = false;

        if (!equal) {
          fprintf(stderr, "\nError:\n");
          fprintf(stderr, "\tsupertext = ");
          for (std::uint64_t j = 0; j < supertext_length; ++j)
            fprintf(stderr, "%c", supertext[j]);
          fprintf(stderr, "\n");
          fprintf(stderr, "\ttext_beg = %lu\n", text_beg);
          fprintf(stderr, "\ttext_length = %lu\n", text_length);
          fprintf(stderr, "\ttail_length = %lu\n", tail_length);
          fprintf(stderr, "\ttail_gt_begin = ");
          for (std::uint64_t j = 0; j < tail_length; ++j)
            fprintf(stderr, "%lu ", (std::uint64_t)tail_gt_begin_reversed->get(j));
          fprintf(stderr, "\n");
          fprintf(stderr, "\tsa: ");
          for (std::uint64_t j = 0; j < text_length; ++j)
            fprintf(stderr, "%lu ", (std::uint64_t)bwtsa[j].m_sa);
          fprintf(stderr, "\n");
          fprintf(stderr, "\tmax_block_size = %lu\n", max_block_size);
          for (std::uint64_t row = 0; row < n_blocks; ++row) {
            for (std::uint64_t col = row + 1; col < n_blocks; ++col) {
              if (computed_matrix[row][col] != correct_matrix[row][col]) {
                fprintf(stderr, "\t\tcomputed[%lu][%lu] = %lu, correct[%lu][%lu] = %lu\n",
                    row, col, computed_matrix[row][col], row, col, correct_matrix[row][col]);
              }
            }
          }
          std::exit(EXIT_FAILURE);
        }
        // Clean up.
        delete tail_gt_begin_reversed;
        for (std::uint64_t row = 0; row < n_blocks; ++row) {
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
//  fprintf(stderr, "combinations = %lu\n", combinations);
}


template<typename saidx_t>
void test_all(std::uint64_t max_supertext_length) {
  for (std::uint64_t supertext_length = 1; supertext_length <= max_supertext_length; ++supertext_length) {
    fprintf(stderr,"TEST, sizeof(saidx_t) = %lu, length = %lu\r", (std::uint64_t)sizeof(saidx_t), supertext_length);
    std::uint8_t *supertext = new std::uint8_t[supertext_length];
    std::uint64_t testcases = (1L << supertext_length);
    for (std::uint64_t supertext_enc = 0L, dbg = 0L; supertext_enc < testcases; ++supertext_enc, ++dbg) {
      if (dbg == 100) {
        fprintf(stderr,"TEST, sizeof(saidx_t) = %lu, length = %lu: %lu (%.0Lf%%)\r", (std::uint64_t)sizeof(saidx_t),
            supertext_length, supertext_enc + 1, ((supertext_enc + 1) * 100.L) / testcases);
        dbg = 0;
      }

      for (std::uint64_t j = 0; j < supertext_length; ++j)
        if (supertext_enc & (1L << j)) supertext[j] = 'b';
        else supertext[j] = 'a';

      test<saidx_t>(supertext, supertext_length);
    }
    delete[] supertext;
    fprintf(stderr,"TEST, sizeof(saidx_t) = %lu, length = %lu: \033[22;32mPASSED\033[0m%10s\n",
        (std::uint64_t)sizeof(saidx_t), supertext_length, "");
  }
}

int main(int, char **) {
  std::srand(std::time(0) + getpid());

  test_all<int>   (15);
  test_all<uint40>(15);

  fprintf(stderr, "All tests passed.\n");
}

