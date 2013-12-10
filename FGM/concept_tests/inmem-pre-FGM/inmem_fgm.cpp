// Initial version of FGM, where gt arrays and the text are hold in memory
// but all other things are final. In particular, gap and partial SA are stored
// on disk and merged in the same way and they will be in the final version.

// TODO better memory management (less allocation).
// TODO: streaming left-to-right, not right-to-left.
// TODO: encoding single gt value on a bit, not byte.
// TODO: make the rightmost block the smallest.

#include <ctime>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <string>

#include "sais.hxx"
#include "utils.h"
#include "rank_4n.h"
#include "suffix_ranking.h"
#include "buffered_gap_array.h"
#include "FGM_fast_merge.h"

// Write the gap array to file using v-byte encoding.
void write_gap_to_file(buffered_gap_array *gap, std::string fname) {
  FILE *f = utils::open_file(fname.c_str(), "w");

  static const int bufsize = (1 << 19);
  unsigned char *buf = new unsigned char[bufsize];
  int filled = 0;

  gap->flush();
  std::sort(gap->excess.begin(), gap->excess.end());
  for (int j = 0, pos = 0; j < gap->length; ++j) {
    int c = 0;
    while (pos < (int)gap->excess.size() && gap->excess[pos] == j) ++pos, ++c;
    int gap_j = gap->count[j] + (c << 8);
    
    // fwrite(&gap_j, sizeof(int), 1, f);
    // v-byte encoding of gap_j:
    while (gap_j > 127) {
      buf[filled++] = ((gap_j & 0x7f) | 0x80);
      gap_j >>= 7;
      if (filled == bufsize)
        utils::add_objects_to_file<unsigned char>(buf, filled, f);
    }
    buf[filled++] = gap_j;
    if (filled == bufsize)
      utils::add_objects_to_file<unsigned char>(buf, filled, f);
  }
  if (filled)
    utils::add_objects_to_file<unsigned char>(buf, filled, f);

  delete[] buf;
  fclose(f);
}

void inmem_FGM(unsigned char *text, int length, int max_block_size) {
  int end = length, block_id = 0;
  unsigned char *gt = NULL;
  while (end > 0) {
    int beg = std::max(end - max_block_size, 0);
    int block_size = end - beg;
    // Current block: B = text[beg..end). If A = text[end..length)
    // then gt[end..length) holds the gt array of A.

    // 1. We start by computing the gt_eof bitvector.
    unsigned char *gt_eof = new unsigned char[block_size];
    std::fill(gt_eof, gt_eof + block_size, 0);
    compute_bitmap(text + end, length - end, text + beg, block_size, gt, gt_eof);

    // 2. Remap symbols of B on a copy of B to compute ordering of
    // suffixes of BA start start in B. Then restore original B back.
    unsigned char *B = new unsigned char[block_size];
    std::copy(text + beg, text + end, B);
    unsigned char last = B[block_size - 1];
    for (int j = 0; j < block_size - 1; ++j)
      if (B[j] > last || (B[j] == last && gt_eof[j + 1])) B[j] += 2;
    ++B[block_size - 1];
    int *SA = new int[block_size];
    saisxx(B, SA, block_size);
    std::copy(text + beg, text + end, B);

    // 3. Store partial SA to disk and compute the prefix of new_gt array.
    unsigned char *new_gt = new unsigned char[length - beg];
    int whole_suffix_pos = 0;
    for (int k = 0; k < block_size; ++k) if (!SA[k]) whole_suffix_pos = k;
    for (int k = 0; k < block_size; ++k) new_gt[SA[k]] = (k > whole_suffix_pos);

    // 4. Store partial SA on disk
    for (int k = 0; k < block_size; ++k) SA[k] += beg;
    utils::write_ints_to_file(SA, block_size, "sparseSA." + utils::intToStr(block_id));
    for (int k = 0; k < block_size; ++k) SA[k] -= beg;

    // ************************************************************************/
    /*                     EXTRA TESTS, NOW ALSO PASSING                      */
    /**************************************************************************/
    /*int *dbg_SA = new int[length - beg];
    saisxx(text + beg, dbg_SA, length - beg);
    int *dbg_sparseSA = new int[block_size];
    int *dbg_gap = new int[block_size + 1];
    std::fill(dbg_gap, dbg_gap + block_size + 1, 0);
    for (int dbg_j = 0, dbg_jj = 0; dbg_j < length - beg; ++dbg_j)
      if (dbg_SA[dbg_j] < block_size) dbg_sparseSA[dbg_jj++] = dbg_SA[dbg_j];
      else ++dbg_gap[dbg_jj];
    for (int dbg_j = 0; dbg_j < block_size; ++dbg_j)
      if (SA[dbg_j] != dbg_sparseSA[dbg_j]) {
        fprintf(stderr, "Error!\n");
        if (length < 100) {
          fprintf(stderr, "  text = %s\n", text);
          fprintf(stderr, "  max_block_size = %d\n", max_block_size);
          fprintf(stderr, "when processing block [%d..%d):\n", beg, end);
          fprintf(stderr, "  correct sparse SA:  ");
          for (int k = 0; k < block_size; ++k)
            fprintf(stderr, "%d ", dbg_sparseSA[k]);
          fprintf(stderr, "\n");
          fprintf(stderr, "  computed sparse SA: ");
          for (int k = 0; k < block_size; ++k)
            fprintf(stderr, "%d ", SA[k]);
          fprintf(stderr, "\n");
        }
        std::exit(EXIT_FAILURE);
      }
    delete[] dbg_SA;
    delete[] dbg_sparseSA;*/
    /**************************************************************************/
    
    // 4. Compute the BWT from SA and build rank on top of it.
    int count[256] = {0}, dollar_pos = 0;
    for (int j = 0; j < block_size; ++j) count[B[j] + 1]++;
    for (int j = 1; j < 256; ++j) count[j] += count[j - 1];
    unsigned char *tmpBWT = (unsigned char *)SA, *BWT = B;
    for (int j = 0, jj = 0; j < block_size; ++j)
      if (SA[j] == 0) dollar_pos = j;
      else tmpBWT[jj++] = B[SA[j] - 1];
    std::copy(tmpBWT, tmpBWT + block_size, BWT);
    delete[] SA;
    rank_4n *rank = new rank_4n(BWT, block_size - 1);
    delete[] BWT;

    // 5. Allocate the buffered n-bytes gap array and do the streaming.
    buffered_gap_array *gap = new buffered_gap_array(block_size + 1);
    int i = 0, scan_length = length - end;
    for (int j = scan_length - 1, idx = length - 1; j >= 0; --j, --idx) {
      unsigned char c = text[idx];
      i = count[c] + rank->rank(i - (i > dollar_pos), c);
      if (c == last && j + 1 != scan_length && gt[j + 1]) ++i;
      new_gt[block_size + j] = (i > whole_suffix_pos);
      gap->increment(i);
    }

    // 6. Store gap array to disk.
    write_gap_to_file(gap, "gap." + utils::intToStr(block_id));

    /**************************************************************************/
    /*                    EXTRA TESTS, NOW ALSO PASSING                       */
    /**************************************************************************/
    /*gap->flush();
    std::sort(gap->excess.begin(), gap->excess.end());
    for (int j = 0, pos = 0; j <= block_size; ++j) {
      int c = 0;
      while (pos < (int)gap->excess.size() && gap->excess[pos] == j) ++pos, ++c;
      int gap_j = gap->count[j] + (c << 8);
      if (dbg_gap[j] != gap_j) {
        fprintf(stderr, "Error!\n");
        if (length < 100) {
          fprintf(stderr, "  text = %s\n", text);
          fprintf(stderr, "  max_block_size = %d\n", max_block_size);
          fprintf(stderr, "  processing block [%d..%d):\n", beg, end);
          fprintf(stderr, "  correct  gap[%d] == %d\n", j, dbg_gap[j]);
          fprintf(stderr, "  computed gap[%d] == %d\n", j, gap_j);
        }
        std::exit(EXIT_FAILURE);
      }
    }
    delete[] dbg_gap;*/
    /**************************************************************************/
    
    // 7. Clean up.
    delete[] gt;
    gt = new_gt;
    delete gap;
    delete rank;
    end = beg;
    ++block_id;
  }
  
  if (gt)
    delete[] gt;

  // Merge the resulting 'gap' and 'sparseSA' arrays
  // and delete the files.
  int *computed_SA = new int[length];
  FGM_fast_merge(computed_SA, length, max_block_size);
  /*for (int i = 0; i < block_id; ++i) {
    utils::file_delete("sparseSA." + utils::intToStr(i));
    utils::file_delete("gap." + utils::intToStr(i));
  }*/

  // Compare the resulting SA to correct SA.
  int *correct_SA = new int[length];
  saisxx(text, correct_SA, length);
  if (!std::equal(correct_SA, correct_SA + length, computed_SA)) {
    fprintf(stderr, "Error!\n");
    if (length < 10000) {
      fprintf(stderr, "  text = %s\n", text);
      fprintf(stderr, "  computed SA: ");
      for (int k = 0; k < length; ++k) fprintf(stderr, "%d ", computed_SA[k]);
      fprintf(stderr, "\n");
      fprintf(stderr, "  correct SA:  ");
      for (int k = 0; k < length; ++k) fprintf(stderr, "%d ", correct_SA[k]);
      fprintf(stderr, "\n");
    }
    std::exit(EXIT_FAILURE);
  }

  delete[] correct_SA;
  delete[] computed_SA;
}

// Test many string chosen according to given paranters.
void test_random(int testcases, int max_length, int max_sigma) {
  fprintf(stderr,"TEST, testcases = %d, max_n = %d, max_sigma = %d\n",
      testcases, max_length, max_sigma);
  unsigned char *text = new unsigned char[max_length + 1];

  for (int tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {
    // Print progress information.
    if (dbg == 1000) {
      fprintf(stderr,"%d (%.2Lf%%)\r", tc, (tc * 100.L) / testcases);
      dbg = 0;
    }

    // Generate string.
    int length = utils::random_int(2, max_length);
    int sigma = utils::random_int(2, max_sigma);
    int block_size = 0;
    do block_size = utils::random_int(1, length);
    while (20 * block_size <= length);
    if (max_sigma <= 26) utils::fill_random_letters(text, length, sigma);
    else utils::fill_random_string(text, length, sigma);
    text[length] = 0;

    // Run the test on generated string.
    inmem_FGM(text, length, block_size);
  }

  // Clean up.
  delete[] text;
}

int main(int, char **) {
  srand(time(0) + getpid());

  // Run tests.
  fprintf(stderr, "Testing inmem pre-FGM.\n");
  test_random(50000, 10,      5);
  test_random(50000, 10,    254);
  test_random(5000, 100,      5);
  test_random(5000, 100,    254);
  test_random(500, 1000,      5);
  test_random(500, 1000,    254);
  test_random(50, 10000,     5);
  test_random(50, 10000,   254);
  test_random(5, 100000,      5);
  test_random(5, 100000,    254);
  fprintf(stderr,"All tests passed.\n");

  return 0;
}

