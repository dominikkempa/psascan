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
#include "srank.h"
#include "gap_array.h"
#include "merge.h"
#include "bitvector.h"
#include "stream.h"

void FGM(std::string filename, int max_block_size) {
  long length = utils::file_size(filename);
  int end = length, block_id = 0, prev_end = length;
  text_reader *reader = new text_reader(filename);
  while (end > 0) {
    int beg = std::max(end - max_block_size, 0);
    int block_size = end - beg; // B = text[beg..end), current block

    // 1. Read current and previously processed block.
    unsigned char *B = new unsigned char[block_size];
    reader->read_block(beg, block_size, B);
    unsigned char *prevB = new unsigned char[max_block_size];
    int prev_block_size = prev_end - end;
    reader->read_block(end, prev_block_size, prevB);

    // 2. Compute gt_eof.
    bitvector *gt_eof = new bitvector(block_size);
    bitvector *gt_head_bv = end < length ? new bitvector("gt_head") : NULL;
    compute_gt_eof(prevB, prev_block_size, B, block_size, gt_head_bv, gt_eof);
    delete gt_head_bv;
    delete[] prevB;

    // 2. Remap symbols of B to compute ordering of suffixes of BA
    // starting in B and then restore original B.
    unsigned char last = B[block_size - 1];
    for (int j = 0; j < block_size - 1; ++j)
      if (B[j] > last || (B[j] == last && gt_eof->get(j + 1))) B[j] += 2;
    delete gt_eof;
    ++B[block_size - 1];
    int *SA = new int[block_size];
    saisxx(B, SA, block_size);
    reader->read_block(beg, end - beg, B);

    // 3. Store the head of new gt bitvector to disk.
    int whole_suffix_pos = 0;
    for (int k = 0; k < block_size; ++k)
      if (!SA[k]) { whole_suffix_pos = k; break; }
    bitvector *new_gt_head = new bitvector(block_size);
    for (int k = whole_suffix_pos + 1; k < block_size; ++k)
      new_gt_head->set(block_size - 1 - SA[k]);
    new_gt_head->save("new_gt_head");
    delete new_gt_head;

    // 4. Store partial SA on disk.
    for (int k = 0; k < block_size; ++k) SA[k] += beg;
    utils::write_ints_to_file(SA, block_size, "sparseSA." + utils::intToStr(block_id));
    for (int k = 0; k < block_size; ++k) SA[k] -= beg;

    // 5. Compute the BWT from SA and build rank on top of it.
    int count[256] = {0}, dollar_pos = 0;
    for (int j = 0; j < block_size; ++j) count[(int)B[j] + 1]++;
    for (int j = 1; j < 256; ++j) count[j] += count[j - 1];
    unsigned char *tmpBWT = (unsigned char *)SA, *BWT = B;
    for (int j = 0, jj = 0; j < block_size; ++j)
      if (SA[j] == 0) dollar_pos = j;
      else tmpBWT[jj++] = B[SA[j] - 1];
    std::copy(tmpBWT, tmpBWT + block_size - 1, BWT);
    delete[] SA;
    rank_4n *rank = new rank_4n(BWT, block_size - 1);
    delete[] BWT;

    // 6. Allocate the gap array, do the streaming and store gap to disk.
    buffered_gap_array *gap = new buffered_gap_array(block_size + 1);
    bit_stream_writer *new_gt_tail = new bit_stream_writer("new_gt_tail");
    bit_stream_reader *gt_tail = prev_end < length ? new bit_stream_reader("gt_tail") : NULL;
    int i = 0;
    unsigned char next_gt = 0;
    reader->init_backward_streaming();
    for (long j = length - 1; j >= prev_end; --j) { // stream the tail
      unsigned char c = reader->read_next();        // c = text[j]
      i = count[c] + rank->rank(i - (i > dollar_pos), c);
      if (c == last && next_gt) ++i;               // next_gt = gt[j + 1]
      new_gt_tail->write(i > whole_suffix_pos);
      gap->increment(i);
      next_gt = gt_tail->read();
    }
    delete gt_tail;
    bit_stream_reader *gt_head = end < length ? new bit_stream_reader("gt_head") : NULL;
    for (int j = prev_end - 1; j >= end; --j) { // stream the head
      unsigned char c = reader->read_next();    // c = text[j]
      i = count[c] + rank->rank(i - (i > dollar_pos), c);
      if (c == last && next_gt) ++i;            // next_gt = gt[j + 1]
      new_gt_tail->write(i > whole_suffix_pos);
      gap->increment(i);
      next_gt = gt_head->read();
    }
    delete gt_head;
    delete new_gt_tail;
    gap->save_to_file("gap." + utils::intToStr(block_id));

    // 7. Clean up.
    delete gap;
    delete rank;
    prev_end = end;
    end = beg;
    ++block_id;
    utils::execute("mv new_gt_head gt_head");
    utils::execute("mv new_gt_tail gt_tail");
  }

  delete reader;

  // Merge gap and sparseSA arrays into final SA.
  std::string out_filename = filename + ".sa";
  merge(length, max_block_size, out_filename);

  // Delete auxiliary files.
  for (int i = 0; i < block_id; ++i) {
    utils::file_delete("sparseSA." + utils::intToStr(i));
    utils::file_delete("gap." + utils::intToStr(i));
  }
  if (utils::file_exists("gt_head")) utils::file_delete("gt_head");
  if (utils::file_exists("gt_tail")) utils::file_delete("gt_tail");
}

// Test many string chosen according to given paranters.
void test_random(int testcases, int max_length, int max_sigma) {
  fprintf(stderr,"TEST, testcases = %d, max_n = %d, max_sigma = %d\n",
      testcases, max_length, max_sigma);
  unsigned char *text = new unsigned char[max_length + 1];
  int *SA = new int[max_length];
  int *computed_SA = new int[max_length];

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
    std::string filename = "/tmp/in" + utils::random_string_hash();
    utils::write_text_to_file(text, length, filename);

    // Run the test on generated string.
    FGM(filename, block_size);
    utils::file_delete(filename);
    
    // Compare the result to correct SA.
    saisxx(text, SA, length);
    FILE *f = utils::open_file(filename + ".sa", "r");
    utils::read_objects_from_file<int>(computed_SA, length, f);
    fclose(f);
    utils::file_delete(filename + ".sa");
    if (!std::equal(SA, SA + length, computed_SA)) {
      fprintf(stderr, "Error!\n");
      if (length < 1000) {
        fprintf(stderr, "  text = %s\n", text);
        fprintf(stderr, "  computed SA: ");
        for (int k = 0; k < length; ++k) fprintf(stderr, "%d ", computed_SA[k]);
        fprintf(stderr, "\n");
        fprintf(stderr, "  correct SA:  ");
        for (int k = 0; k < length; ++k) fprintf(stderr, "%d ", SA[k]);
        fprintf(stderr, "\n");
      }
      std::exit(EXIT_FAILURE);
    }
  }

  // Clean up.
  delete[] text;
  delete[] SA;
  delete[] computed_SA;
}

int main(int, char **) {
  srand(time(0) + getpid());
  fprintf(stderr, "Testing FGM.\n");
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
}

