// Handling gt bitvectors
//=======================
//
// At the beginning of phase, the gt bitvector for A is split into two pieces:
//  1. the gt bitvector for the previous block (file: 'gt_head')
//  2. the gt bitvector for the rest of A (file: 'gt_tail')
// 
// We use these in the current phase as follows:
//  a) first, for the computation of gt_eof, we need a random access to
//     the gt of the prev block. To achieve this, we simply load the
//     bitvector into 'bitvector' class from the 'gt_head' file and use
//     it during the Crochemore's algorithm. Of course the bitvector class
//     encodes 8bits/byte and provides random access. After we're done with
//     this, we drop this bitvector, although, it could be kept in memory for
//     a little longer, because we need it later.
//  b) next, as we are about to drop (after saving to disk) the sparseSA
//     we save the head of the new gt (resulting from this phase) into
//     'new_gt_head' (which will be renamed later). For the current phase
//     we will not need the 'new_gt_head' any more.
//  c) next follows the scanning part, at this point we need the stream gt of
//     A. To achieve this, we first stream the 'gt_tail' and then the 'gt_head'
//     (note: we stream it now, whereas at the beginning we have accessed it
//     randomly). Simultanously, we use bit_stream_writer to write the
//     tail of the new gt array into 'new_gt_tail'.
//  d) finally, we rename 'new_gt_head' into 'gt_head' and 'new_gt_tail'
//     into 'gt_tail' and are now ready to process next block.

// TODO: better memory management (less allocation).
// TODO: streaming left-to-right, not right-to-left.
// TODO: make the rightmost block the smallest.

#include <ctime>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <string>

#include "sais.hxx"
#include "utils.h"
#include "rank.h"
#include "suffix_ranking.h"
#include "buffered_gap_array.h"
#include "merge.h"
#include "bitvector.h"
#include "bit_stream_reader.h"
#include "bit_stream_writer.h"

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

void FGM(unsigned char *text, int length, int max_block_size) {
  int end = length, block_id = 0, prev_end = length;
  while (end > 0) {
    int beg = std::max(end - max_block_size, 0);
    int block_size = end - beg;
    // Current block: B = text[beg..end). If A = text[end..length)
    // then gt[end..length) holds the gt array of A.

    // 1. We start by computing the gt_eof bitvector.
    bitvector *gt_eof = new bitvector(block_size);
    unsigned char *B = new unsigned char[block_size];
    std::copy(text + beg, text + end, B);
    unsigned char *prevB = new unsigned char[block_size];
    std::copy(text + end, text + prev_end, prevB);
    int prevB_size = prev_end - end;
    bitvector *gt_head_bv = (end < length) ? (new bitvector("gt_head")) : NULL;
    compute_gt_eof(prevB, prevB_size, B, block_size, gt_head_bv, gt_eof);
    delete gt_head_bv;
    delete[] prevB;

    // 2. Remap symbols of B on a copy of B to compute ordering of
    // suffixes of BA starting in B. Then restore original.
    unsigned char last = B[block_size - 1];
    for (int j = 0; j < block_size - 1; ++j)
      if (B[j] > last || (B[j] == last && gt_eof->get(j + 1))) B[j] += 2;
    delete gt_eof;
    ++B[block_size - 1];
    int *SA = new int[block_size];
    saisxx(B, SA, block_size);
    std::copy(text + beg, text + end, B);

    // 3. Store partial SA to disk and compute the prefix of new_gt array.
    int whole_suffix_pos = 0;
    for (int k = 0; k < block_size; ++k)
      if (!SA[k]) { whole_suffix_pos = k; break; }
    bitvector *new_gt_head = new bitvector(block_size);
    for (int k = whole_suffix_pos + 1; k < block_size; ++k)
      new_gt_head->set(block_size - 1 - SA[k]);
    new_gt_head->save("new_gt_head");
    delete new_gt_head;

    // 4. Store partial SA on disk
    for (int k = 0; k < block_size; ++k) SA[k] += beg;
    utils::write_ints_to_file(SA, block_size, "sparseSA." + utils::intToStr(block_id));
    for (int k = 0; k < block_size; ++k) SA[k] -= beg;

    // 5. Compute the BWT from SA and build rank on top of it.
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

    // 6. Allocate the buffered n-bytes gap array and do the streaming.
    buffered_gap_array *gap = new buffered_gap_array(block_size + 1);
    bit_stream_writer *new_gt_tail = new bit_stream_writer("new_gt_tail");
    bit_stream_reader *gt_tail = (prev_end < length) ? (new bit_stream_reader("gt_tail")) : NULL;
    int i = 0;
    unsigned char next_gt = 0;
    for (long j = length - 1; j >= prev_end; --j) { // stream the tail, except prev block
      // Invariant: next_gt = gt[j + 1]
      unsigned char c = text[j];
      i = count[c] + rank->rank(i - (i > dollar_pos), c);
      if (c == last && next_gt) ++i;
      new_gt_tail->write(i > whole_suffix_pos);
      gap->increment(i);
      next_gt = gt_tail->read();
    }
    delete gt_tail;
    bit_stream_reader *gt_head = (end < length) ? (new bit_stream_reader("gt_head")) : NULL;
    for (int j = prev_end - 1; j >= end; --j) { // stream previous block
      // Invariant: next_gt = gt[j + 1]
      unsigned char c = text[j];
      i = count[c] + rank->rank(i - (i > dollar_pos), c);
      if (c == last && next_gt) ++i;
      new_gt_tail->write(i > whole_suffix_pos);
      gap->increment(i);
      next_gt = gt_head->read();
    }
    delete gt_head;
    delete new_gt_tail;

    // 7. Store gap array to disk.
    write_gap_to_file(gap, "gap." + utils::intToStr(block_id));

    // 8. Clean up.
    delete gap;
    delete rank;
    prev_end = end;
    end = beg;
    ++block_id;
    utils::execute("mv new_gt_head gt_head");
    utils::execute("mv new_gt_tail gt_tail");
  }

  // Merge the resulting 'gap' and 'sparseSA' arrays
  // and delete the files.
  int *computed_SA = new int[length];
  merge(computed_SA, length, max_block_size);
  for (int i = 0; i < block_id; ++i) {
    utils::file_delete("sparseSA." + utils::intToStr(i));
    utils::file_delete("gap." + utils::intToStr(i));
  }

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
    if (dbg == 10) {
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
    // int length = strlen("cccedbacba");
    // strcat((char *)text, "cccedbacba");
    // int block_size = 5;
    text[length] = 0;

    // Run the test on generated string.
    FGM(text, length, block_size);
  }

  // Clean up.
  delete[] text;
}

int main(int, char **) {
  srand(time(0) + getpid());

  // Run tests.
  fprintf(stderr, "Testing inmem pre-FGM.\n");
  test_random(5000, 10,      5);
  test_random(5000, 10,    254);
  test_random(500, 100,      5);
  test_random(500, 100,    254);
  test_random(50, 1000,      5);
  test_random(50, 1000,    254);
  test_random(50, 10000,     5);
  test_random(50, 10000,   254);
  test_random(5, 100000,      5);
  test_random(5, 100000,    254);
  fprintf(stderr,"All tests passed.\n");

  return 0;
}

