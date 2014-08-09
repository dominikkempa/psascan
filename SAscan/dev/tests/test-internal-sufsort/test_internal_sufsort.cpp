#include <cstdio>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <thread>
#include <unistd.h>

#include "divsufsort.h"
#include "divsufsort64.h"
#include "bitvector.h"
#include "utils.h"


//==============================================================================
// Compute ms-decomposition of T[0..n) from ms-decomposition of T[0..n-1).
// The result is returned via updated values s, p, r.
//==============================================================================
void next(unsigned char *T, long n, long &s, long &p, long &r) {
  if (n == 1) { s = 0; p = 1; r = 0; return; }
  long i = n - 1;
  while (i < n) {
    unsigned char a = T[s + r], b = T[i];
    if (a > b) { p = i - s + 1; r = 0; }
    else if (a < b) { i -= r; s = i; p = 1; r = 0; }
    else { ++r; if (r == p) r = 0; } ++i;
  }
}


//==============================================================================
// Compute bitvectors bv[0..ref_pos) and decided[0..ref_pos), where:
//   decided[i] == 1 iff lcp(text[0..txtlen), text[ref_pos..txtlen)) < max_lcp
//   gt[i] == 1 iff decided[i] == 1 and text[0..txtlen) > text[ref_pos..txtlen)
//==============================================================================
void compute_partial_gt(unsigned char *text, long ref_pos, long max_lcp,
    long txtlen, bitvector *gt, bitvector *decided, bool &all_decided) {
  long i = 0, el = 0, s = 0, p = 0, r = 0;
  long i_max = 0, el_max = 0, s_max = 0, p_max = 0, r_max = 0;

  unsigned char *pat = text + ref_pos;

  all_decided = true;
  while (i < ref_pos) {
    while (el < max_lcp && text[i + el] == pat[el])
      next(pat, ++el, s, p, r); 
     
    if (el < max_lcp) {
      decided->set(i);
      if (text[i + el] > pat[el]) gt->set(i);
    } else if (ref_pos + el == txtlen) {
      decided->set(i);
      gt->set(i);
    } else all_decided = false;

    long j = i_max;
    if (el > el_max) {
      std::swap(el, el_max);
      std::swap(s, s_max);
      std::swap(p, p_max);
      std::swap(r, r_max);
      i_max = i;
    }

    if (p > 0 && 3 * p <= el && !memcmp(pat, pat + p, s)) {
      for (long k = 1; k < std::min(p, ref_pos - i); ++k) {
        if (decided->get(j + k)) decided->set(i + k);
        if (gt->get(j + k)) gt->set(i + k);
      }

      i += p;
      el -= p;
    } else {
      long h = (el / 3) + 1;
      for (long k = 1; k < std::min(h, ref_pos - i); ++k) {
        if (decided->get(j + k)) decided->set(i + k);
        if (gt->get(j + k)) gt->set(i + k);
      }

      i += h;
      el = 0;
      s = 0;
      p = 0;
    }
  }
}


//==============================================================================
// Set all undecided bits inside the given microblock (that is, the range
// [mb_beg..mb_end)) of all gt bitvectors to their correct values.
//==============================================================================
void compute_final_gt(long length, long max_block_size,
    long mb_beg, long mb_end, bitvector** &gt, bitvector** &decided,
    bool *all_decided) {

  // Go through blocks right-to-left.
  long n_blocks = (length + max_block_size - 1) / max_block_size;
  for (long i = n_blocks - 1; i >= 0; --i) {
    if (all_decided[i]) continue; // optimization
    long beg = i * max_block_size;
    long end = std::min(beg + max_block_size, length);
    long block_size = end - beg;
    long this_mb_end = std::min(block_size, mb_end);

    // Scan the bits inside the microblock of block i.
    for (long j = mb_beg; j < this_mb_end; ++j) {
      if (decided[i]->get(j) == false) {
        // j-th bit of gt[i] was undecided -> copy it from the right.
        if (gt[i + 1]->get(j)) gt[i]->set(j);
      }
    }
  }
}


//==============================================================================
// Fully parallel computation of gt bitvectors.
//==============================================================================
void compute_gt(unsigned char *text, long length, bitvector** &gt,
    long max_threads) {
  long max_block_size = (length + max_threads - 1) / max_threads;
  long n_blocks = (length + max_block_size - 1) / max_block_size;


  //----------------------------------------------------------------------------
  // STEP 1: compute gt bitvectors, some bits may still be undecided after this.
  //----------------------------------------------------------------------------

  // Allocate ane zero-initialize (in parallel) bitvectors.
  bitvector **decided = new bitvector*[n_blocks];
  gt = new bitvector*[n_blocks];
  for (long i = 0; i < n_blocks; ++i) {
    long beg = i * max_block_size;
    long end = std::min(beg + max_block_size, length);
    long this_block_size = end - beg;

    decided[i] = new bitvector(this_block_size, max_threads);
    gt[i] = new bitvector(this_block_size, max_threads);
  }

  // all_decided[i] == true, if all bits inside block i were
  // decided in the first state. This can be used by threads in the
  // second stage to completely skip inspecting some blocks.
  bool *all_decided = new bool[n_blocks];

  // Process blocks right-to-left.
  std::thread **threads = new std::thread*[n_blocks];
  for (long i = n_blocks - 1, next_block_size = 0L; i >= 0; --i) {
    long beg = i * max_block_size;
    long end = std::min(beg + max_block_size, length);
    long this_block_size = end - beg;

    // Compute bitvectors 'gt' and 'decided' for block i.
    threads[i] = new std::thread(compute_partial_gt, text + beg,
      this_block_size, next_block_size, length - beg, gt[i],
      decided[i], std::ref(all_decided[i]));

    next_block_size = this_block_size;
  }

  // Wait for the threads to finish and clean up.
  for (long i = 0; i < n_blocks; ++i) threads[i]->join();
  for (long i = 0; i < n_blocks; ++i) delete threads[i];
  delete[] threads;


  //----------------------------------------------------------------------------
  // STEP 2: compute the undecided bits in the gt bitvectors.
  //----------------------------------------------------------------------------
  
  // The size of micro block has to be a multiple of 8, otherwise two
  // threads might try to update the same char inside bitvector.
  long max_microblock_size = (max_block_size + max_threads - 1) / max_threads;
  while (max_microblock_size & 7) ++max_microblock_size;
  long n_microblocks = (max_block_size + max_microblock_size - 1) / max_microblock_size;

  threads = new std::thread*[n_microblocks];
  for (long i = 0; i < n_microblocks; ++i) {
    long mb_beg = i * max_microblock_size;
    long mb_end = (i + 1) * max_microblock_size;
    threads[i] = new std::thread(compute_final_gt, length,
        max_block_size, mb_beg, mb_end, std::ref(gt),
        std::ref(decided), all_decided);
  }

  // Wait for the threads to finish and clean up.
  for (long i = 0; i < n_microblocks; ++i) threads[i]->join();
  for (long i = 0; i < n_microblocks; ++i) delete threads[i];
  for (long i = 0; i < n_blocks; ++i) delete decided[i];
  delete[] threads;
  delete[] decided;
  delete[] all_decided;
}


//==============================================================================
// Rename the given block using its gt bitvector.
//==============================================================================
void rename_block(unsigned char *block, long block_length, bitvector *gt) {
  unsigned char last = block[block_length - 1];
  for (long i = 0; i + 1 < block_length; ++i)
    if (block[i] > last || (block[i] == last && gt->get(i + 1))) ++block[i];
  ++block[block_length - 1];
}


//==============================================================================
// Given gt bitvectors, compute partial suffix arrays of blocks.
// To do this, in parallel:
//   1) rename the blocks
//   2) run divsufsort on each block
//==============================================================================
void compute_partial_sa(unsigned char *text, long text_length,
    bitvector** &gt, int* &partial_sa, long max_threads) {
  long max_block_size = (text_length + max_threads - 1) / max_threads;
  long n_blocks = (text_length + max_block_size - 1) / max_block_size;

  // Rename the blocks in parallel.
  std::thread **threads = new std::thread*[n_blocks];
  for (long i = 0; i < n_blocks; ++i) {
    long block_beg = i * max_block_size;
    long block_end = std::min(block_beg + max_block_size, text_length);
    long block_size = block_end - block_beg;
    threads[i] = new std::thread(rename_block,
        text + block_beg, block_size, gt[i]);
  }
  for (long i = 0; i < n_blocks; ++i) threads[i]->join();
  for (long i = 0; i < n_blocks; ++i) delete threads[i];

  // We can now delete gt bitvectors. Re-renaming back to
  // original blocks can be done without them.
  for (long i = 0; i < n_blocks; ++i) delete gt[i];
  delete[] gt;

  // Compute the partial suffix arrays in parallel.
  // 
  // First, allocate the space.
  partial_sa = new int[text_length];

  // Now run the threads.
  for (long i = 0; i < n_blocks; ++i) {
    long block_beg = i * max_block_size;
    long block_end = std::min(block_beg + max_block_size, text_length);
    long block_size = block_end - block_beg;
    threads[i] = new std::thread(divsufsort, text + block_beg,
        partial_sa + block_beg, (int)block_size);
  }
  for (long i = 0; i < n_blocks; ++i) threads[i]->join();
  for (long i = 0; i < n_blocks; ++i) delete threads[i];
  delete[] threads;
}


//==============================================================================
// Fully parallel computation of partial suffix arrays.
//==============================================================================
void partial_sufsort(unsigned char *text, long text_length, long max_threads) {
  long max_block_size = (text_length + max_threads - 1) / max_threads;
  long n_blocks = (text_length + max_block_size - 1) / max_block_size;

  bitvector **gt;
  compute_gt(text, text_length, gt, max_threads);

  // Check the correctness of gt bitvectors.
  for (long i = n_blocks - 1; i >= 0; --i) {
    long block_beg = i * max_block_size;
    long block_end = std::min(block_beg + max_block_size, text_length);
    long block_size = block_end - block_beg;

    for (long j = 0; j < block_size; ++j) {
      // Compute correct gt[j] for block i.
      long lcp = 0L;
      while (block_end + lcp < text_length && text[block_beg + j + lcp] == text[block_end + lcp]) ++lcp;
      bool correct_gt = (block_end + lcp == text_length || text[block_beg + j + lcp] > text[block_end + lcp]);

      if (gt[i]->get(j) != correct_gt) {
        fprintf(stderr, "Error!\n");
        if (text_length <= 1000) {
          fprintf(stderr, "text: ");
          for (long k = 0; k < text_length; ++k)
            fprintf(stderr, "%c", text[k]);
          fprintf(stderr, "\n");
        }
        fprintf(stderr, "max_threads = %ld\n", max_threads);
        fprintf(stderr, "n_blocks = %ld\n", n_blocks);
        fprintf(stderr, "max_block_size = %ld\n",  max_block_size);
        fprintf(stderr, "block id = %ld, j = %ld\n", i, j);
        fprintf(stderr, "correct_gt = %d, computed_gt = %d\n",
          correct_gt, gt[i]->get(j));
        std::exit(EXIT_FAILURE);
      }
    }
  }

  int *full_sa = new int[text_length];
  divsufsort(text, full_sa, (int)text_length);

  int *partial_sa;
  compute_partial_sa(text, text_length, gt, partial_sa, max_threads);

  // The the correctness of partial suffix arrays.
  int *temp = new int[max_block_size];
  for (long i = 0; i < n_blocks; ++i) {
    long block_beg = i * max_block_size;
    long block_end = std::min(block_beg + max_block_size, text_length);
    long block_size = block_end - block_beg;
    
    // To verify correctness of partial_sa[i] we collect suffixes starting
    // inside [beg..end) in a temp array and compare to partial_sa[i].
    for (long j = 0, ptr = 0; j < text_length; ++j)
      if (block_beg <= full_sa[j] && full_sa[j] < block_end)
        temp[ptr++] = full_sa[j];

    int *partial_sa_i = partial_sa + block_beg;
    // Compare temp and partial_sa[i].
    for (long j = 0; j < block_size; ++j) {
      if (temp[j] != block_beg + partial_sa_i[j]) {
        fprintf(stderr, "Error!\n");
        if (text_length <= 1000) {
          fprintf(stderr, "text: ");
          for (long jj = 0; jj < text_length; ++jj)
            fprintf(stderr, "%c", text[jj]);
          fprintf(stderr, "\n");
        }
        fprintf(stderr, "block = %ld, beg = %ld, end = %ld\n", i, block_beg, block_end);
        fprintf(stderr, "  temp[%ld] == %d\n", j, temp[j]);
        fprintf(stderr, "  beg + partial_sa[%ld][%ld] == %ld\n", i, j,
            block_beg + partial_sa_i[j]);
        std::exit(EXIT_FAILURE);
      }
    }
  }
  delete[] full_sa;
  delete[] temp;
  delete[] partial_sa;
}


void test_random(int testcases, long max_length, int max_sigma) {
  fprintf(stderr,"TEST, testcases = %d, max_n = %ld, max_sigma = %d\n",
      testcases, max_length, max_sigma);
  unsigned char *text = new unsigned char[max_length + 1];

  for (int tc = 0; tc < testcases; ++tc) {
    // Print progress information.
    fprintf(stderr,"%d (%.2Lf%%)\r", tc, (tc * 100.L) / testcases);

    // Generate string.
    long length = utils::random_long(1, max_length);
    int sigma = utils::random_int(2, max_sigma);
    if (max_sigma <= 26) utils::fill_random_letters(text, length, sigma);
    else utils::fill_random_string(text, length, sigma);
    
    long max_threads = utils::random_long(1, 50);

    // Run the test on generated string.
    partial_sufsort(text, length, max_threads);
  }

  // Clean up.
  delete[] text;
}


int main() {
  std::srand(std::time(0) + getpid());

  test_random(100000,   10,      5);
  test_random(100000,   10,      255);
  test_random(100000,   1000,    5);
  test_random(100000,   1000,    255);
  test_random(10000,    100000,  5);
  test_random(10000,    100000,  255);
  test_random(1000,     1000000, 5);
  test_random(1000,     1000000, 255);
  fprintf(stderr,"All tests passed.\n");
}
