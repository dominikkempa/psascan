// For a given string, we generate the random partition into blocks. Then for
// each block we compute:
//   * the sparse suffix array (containing only suffixes starting inside)
//   * the gap array, as computed in the FGM
// Then we perform the fast merging (as suggested by Juha) and compare the
// result to the complete suffix array of the text.

#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <vector>
#include <string>

#include "sais.hxx"
#include "utils.h"


// Represents the output of a single FGM phase: gap and sparse SA.
struct fgm_phase_output {
  fgm_phase_output() :gap(NULL), sparse_sa(NULL) {}
  
  // The block is text[beg..end). SA is the suf array of text[0..length).
  void construct(int length, int *SA, int beg, int end) {
    block_length = end - beg;
    gap = new int[block_length + 1];
    sparse_sa = new int[block_length];
    
    std::fill(gap, gap + block_length + 1, 0);
    for (int i = 0, j = 0; i < length; ++i) {
      if (SA[i] < beg) continue;
      else if (SA[i] < end) sparse_sa[j++] = SA[i];
      else ++gap[j];
    }
  }

  ~fgm_phase_output() {
    if (gap) delete[] gap;
    if (sparse_sa) delete[] sparse_sa;
  }

  int block_length, text_length;
  int *gap, *sparse_sa;
};

void test(unsigned char *text, int length) {  
  // Generate the random partition of text, if possible, at least
  // min_blocks blocks.

  int min_blocks = utils::random_int(1, 10);
  std::vector<int> ends;
  for (int endpos = length, begpos; endpos > 0; endpos = begpos) {
    begpos = utils::random_int(0,
      std::min(endpos - 1, (length + min_blocks - 1) / min_blocks));
    ends.push_back(endpos);
  }
  std::reverse(ends.begin(), ends.end());
  int n_block = (int)ends.size();
  

  // Compute the suffix array of the text.
  int *SA = new int[length];
  saisxx(text, SA, length);

  // Compute the FGM output.
  fgm_phase_output *output = new fgm_phase_output[n_block];
  for (int i = 0; i < n_block; ++i) {
    int end = ends[i];
    int beg = (i == 0) ? 0 : ends[i - 1];
    output[i].construct(length, SA, beg, end);
  }

  // Do the fast merging.
  // For each block we keep two numbers:
  //  - # of already extracted suffixes that are lex-smaller than but occur
  //    later in the text. We call this the 'block rank'.
  //  - the positions of the smallest non-extracted suffix of the block, among
  //    the suffixes starting in and after the block, in other words: this is
  //    simply sum of gaps + number of suffixes in the block up to the next
  //    non-extracted suffix. We call this the 'suffix_rank' of a block.
  int *computed_sa = new int[length]; // this is the output of merging

  int *block_rank  = new int[n_block];
  int *suffix_rank = new int[n_block];
  int *suf_ptr = new int[n_block]; // first non-extracted suffix in the block.
  for (int i = 0; i < n_block; ++i) suffix_rank[i] = output[i].gap[0];
  std::fill(block_rank, block_rank + n_block, 0);
  std::fill(suf_ptr, suf_ptr + n_block, 0);
  for (int i = 0; i < length; ++i) {  // merge
    // Find the leftmost block j with block_rank[j] == suffix_rank[j].
    int j = 0;
    while (j < n_block && block_rank[j] != suffix_rank[j]) ++j;

    // Extract the suffix.
    computed_sa[i] = output[j].sparse_sa[suf_ptr[j]];

    // Update suffix_rank[j].
    suffix_rank[j]++;
    suffix_rank[j] += output[j].gap[++suf_ptr[j]];
    
    // Update block_rank[0..j].
    for (int k = 0; k <= j; ++k) ++block_rank[k];
  }

  // Compare the computed suffix array to the correct suffix array.
  if (!std::equal(SA, SA + length, computed_sa)) {
    fprintf(stderr, "\n\033[22;31mFAILED\033[0m\n");
    if (length < 100) {
      fprintf(stderr, "  text = %s\n", text);
      fprintf(stderr, "  SA = ");
      for (int j = 0; j < length; ++j) fprintf(stderr, "%d ", SA[j]);
      fprintf(stderr, "\n");
      fprintf(stderr, "  computed sa = ");
      for (int j = 0; j < length; ++j) fprintf(stderr, "%d ", computed_sa[j]);
      fprintf(stderr, "\n");
    }
    std::exit(EXIT_FAILURE);
  }

  // Clean up.  
  delete[] output;
  delete[] SA;
  delete[] suf_ptr;
  delete[] computed_sa;
  delete[] block_rank;
  delete[] suffix_rank;
}

// Test many string chosen according to given paranters.
void test_random(int testcases, int max_length, int max_sigma) {
  fprintf(stderr,"TEST, testcases = %d, max_n = %d, max_sigma = %d\r",
      testcases, max_length, max_sigma);
  unsigned char *text = new unsigned char[max_length + 1];

  for (int tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {
    // Print progress information.
    if (dbg == 100) {
      fprintf(stderr,"TEST, testcases = %d, max_n = %d, max_sigma = %d: "
          "%d (%.0Lf%%)\r", testcases, max_length, max_sigma, tc,
          (tc * 100.L) / testcases);
      dbg = 0;
    }

    // Generate string.
    int length = utils::random_int(2, max_length);
    int sigma = utils::random_int(2, max_sigma);
    if (max_sigma <= 26) utils::fill_random_letters(text, length, sigma);
    else utils::fill_random_string(text, length, sigma);
    text[length] = 0;

    // Run the test on generated string.
    test(text, length);
  }

  // Clean up.
  delete[] text;

  fprintf(stderr,"TEST, testcases = %d, max_n = %d, max_sigma = %d: "
      "\033[22;32mPASSED\033[0m%10s\n", testcases, max_length, max_sigma, "");
}

int main(int, char **) {
  std::srand(std::time(0) + getpid());

  // Run tests.
  fprintf(stderr, "Testing fast merging.\n");
  test_random(500000, 10,       5);
  test_random(500000, 10,     256);
  test_random(100000, 100,      5);
  test_random(100000, 100,    256);
  test_random(50000,  1000,     5);
  test_random(50000,  1000,   256);
  test_random(10000,  10000,    5);
  test_random(10000,  10000,  256);
  test_random(1000,   100000,   5);
  test_random(1000,   100000, 256);
}

