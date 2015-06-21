// This program verifies that using gt_begin in range refining
// may produce smaller intermediate ranges than without gt_begin,
// but the `left' value of the final range is always correct.
// Possibly the `right' is different.

#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <string>

#include "sais.hxx"
#include "utils.h"
#include "bitvector.h"


// Return true iff text[i..) (but we always stop the comparison
// at text_length) is smaller than pat[0..pat_length).
int lcp_compare_2(unsigned char *text, long text_length, unsigned char *pat, long pat_length,
    long tail_gt_begin_reversed_length, long j, bitvector &tail_gt_begin_reversed, long &lcp) {
  while (lcp < pat_length && j + lcp < text_length && pat[lcp] == text[j + lcp]) ++lcp;

  if (lcp == pat_length) return 0;
  if ((j + lcp < text_length && pat[lcp] < text[j + lcp]) ||
      (j + lcp >= text_length && !(tail_gt_begin_reversed.get(tail_gt_begin_reversed_length - (text_length - j)))))
    return -1;
  else return 1;
}

void refine_range(unsigned char *text, long text_length,
    long tail_gt_begin_reversed_length,
    long block_beg, long *block_psa,                             
    long left, long right,                                       
    bitvector &tail_gt_begin_reversed,                           
    long old_pat_length, long pat_length, unsigned char *pat,    
    long &newleft, long &newright) {                             // OUTPUT: new range

  long low = left - 1;
  long high = right;
  long llcp = old_pat_length;
  long rlcp = old_pat_length;

//  static const long min_discrepancy = /*(1L << 16)*/1L;
//  static const long balancing_constant = /*64L*/1L;

  long min_discrepancy = utils::random_long(0L, 10L);
  long balancing_constant = utils::random_long(1L, 10L);

  while (low + 1 != high) {
    // Invariant: newleft is in the range (low, high].

    // Compute mid.
    // Valid values for mid are: low + 1, .., high - 1.
    long mid = 0L;
    if (llcp + min_discrepancy < rlcp) {
      // Choose the pivot that split the range into two parts of sizes
      // with ratio equal to logd / d.
      long d = rlcp - llcp;
      long logd = utils::log2ceil(d);
      mid = low + 1 + ((high - low - 1) * balancing_constant * logd) / (d + balancing_constant * logd);
    } else if (rlcp + min_discrepancy < llcp) {
      // Choose the pivot that split the range into two parts of sizes
      // with ratio equal to logd / d.
      long d = llcp - rlcp;
      long logd = utils::log2ceil(d);
      mid = high - 1 - ((high - low - 1) * balancing_constant * logd) / (d + balancing_constant * logd);
    } else {
      // Discrepancy between lcp values is small, use standard binary search.
      mid = (low + high) / 2;
    }

    long lcp = std::min(llcp, rlcp);

    if (lcp_compare_2(text, text_length, pat, pat_length, tail_gt_begin_reversed_length,
          block_beg + block_psa[mid], tail_gt_begin_reversed, lcp) <= 0) {
      high = mid;
      rlcp = lcp;
    } else {
      low = mid;
      llcp = lcp;
    }
  }

  newleft = high;

  if (rlcp < pat_length) newright = newleft;
  else {
    high = right;
    rlcp = old_pat_length;

    while (low + 1 != high) {
      // Invariant: newright is in the range (low, high].

      // Compute mid.
      // Valid values for mid are: low + 1, .., high - 1.
      long mid = 0L;
      if (llcp + min_discrepancy < rlcp) {
        // Choose the pivot that split the range into two parts of sizes
        // with ratio equal to logd / d.
        long d = rlcp - llcp;
        long logd = utils::log2ceil(d);
        mid = low + 1 + ((high - low - 1) * balancing_constant * logd) / (d + balancing_constant * logd);
      } else if (rlcp + min_discrepancy < llcp) {
        // Choose the pivot that split the range into two parts of sizes
        // with ratio equal to logd / d.
        long d = llcp - rlcp;
        long logd = utils::log2ceil(d);
        mid = high - 1 - ((high - low - 1) * balancing_constant * logd) / (d + balancing_constant * logd);
      } else {
        // Discrepancy between lcp values is small, use standard binary search.
        mid = (low + high) / 2;
      }

      long lcp = std::min(llcp, rlcp);

      if (lcp_compare_2(text, text_length, pat, pat_length, tail_gt_begin_reversed_length,
            block_beg + block_psa[mid], tail_gt_begin_reversed, lcp) < 0) {
        high = mid;
        rlcp = lcp;
      } else {
        low = mid;
        llcp = lcp;
      }
    }

    newright = high;
  }
}

void compute_gt_begin_reversed(unsigned char *text, long text_length, bitvector *gt_begin_reversed) {
  for (long j = 1; j < text_length; ++j) {
    long lcp = 0L;
    while (j + lcp < text_length && text[lcp] == text[j + lcp]) ++lcp;
    bool gt_j = (j + lcp < text_length && text[lcp] < text[j + lcp]);

    if (gt_j) gt_begin_reversed->set(text_length - j);
  }
}

void test(unsigned char *text, long text_length, long tail_length,
    long block_beg, long block_end) {
  // Compute reversed gt_begin for the tail (tail_gt_begin_reversed);
  bitvector *tail_gt_begin_reversed = new bitvector(tail_length);
  compute_gt_begin_reversed(text + text_length, tail_length, tail_gt_begin_reversed);

  // Compute block partial suffix array (block_psa);
  long block_length = block_end - block_beg;
  int *text_sa = new int[text_length + tail_length];
  long *block_psa = new long[block_length];
  saisxx(text, text_sa, (int)text_length + (int)tail_length);
  for (long j = 0, ptr = 0; j < text_length + tail_length; ++j) {
    if (block_beg <= text_sa[j] && text_sa[j] < block_end)
      block_psa[ptr++] = text_sa[j] - block_beg;
  }
  delete[] text_sa;

  long max_prefix_length = std::min(tail_length, text_length - block_end + 1);

  // Test all possible prefix of the tail:
  unsigned char *tail = text + text_length;
  for (long prefix_length = 0L; prefix_length <= max_prefix_length; ++prefix_length) {
    // Find the range of suffixes sarting inside the block
    // that is prefixed with tail[0..prefix_length).
    // The range will be [left..right).
    long smaller = 0L;
    long equal = 0L;
    for (long j = block_beg; j < block_end; ++j) {
      long lcp = 0L;
      while (lcp < prefix_length && tail[lcp] == text[j + lcp]) ++lcp;
      if (lcp < prefix_length && tail[lcp] > text[j + lcp]) ++smaller;
      else if (lcp == prefix_length) ++equal;
    }
    long left = smaller;
    long right = smaller + equal;

    // Now check all possible extensions of the prefix_length (including extension by zero symbols).
    for (long new_prefix_length = prefix_length; new_prefix_length <= max_prefix_length; ++new_prefix_length) {
      // Compute the correct range for the new_prefix_length.
      smaller = 0L;
      equal = 0L;
      for (long j = block_beg; j < block_end; ++j) {
        long lcp = 0L;
        while (lcp < new_prefix_length && tail[lcp] == text[j + lcp]) ++lcp;
        if (lcp < new_prefix_length && tail[lcp] > text[j + lcp]) ++smaller;
        else if (lcp == new_prefix_length) ++equal;
      }
      long newleft = smaller;
      long newright = smaller + equal;

      // Now test if the same range [newleft..newright) is obtained from the tested function.

      long computed_newleft;
      long computed_newright;

      refine_range(text, text_length, tail_length, block_beg, block_psa, left, right,
        *tail_gt_begin_reversed, prefix_length, new_prefix_length, tail, computed_newleft, computed_newright);

      if (newleft != computed_newleft || newright != computed_newright) {
        fprintf(stderr, "\nError!\n");
        fprintf(stderr, "\ttext = ");
        for (long j = 0; j < text_length; ++j)
          fprintf(stderr, "%c", text[j]);
        fprintf(stderr, "\n");
        fprintf(stderr, "\ttail = ");
        for (long j = 0; j < tail_length; ++j)
          fprintf(stderr, "%c", tail[j]);
        fprintf(stderr, "\n");
        fprintf(stderr, "\tblock = [%ld..%ld)\n", block_beg, block_end);
        fprintf(stderr, "\ttail_gt_begin_reversed: ");
        for (long j = 0; j < tail_length; ++j)
          fprintf(stderr, "%ld ", (long)tail_gt_begin_reversed->get(j));
        fprintf(stderr, "\n");
        fprintf(stderr, "\tblock_psa = ");
        for (long j = 0; j < block_length; ++j)
          fprintf(stderr, "%ld ", block_psa[j]);
        fprintf(stderr, "\n");
        fprintf(stderr, "\tprefix_length = %ld, left = %ld, right = %ld\n", prefix_length, left, right);
        fprintf(stderr, "\tnew_prefix_length = %ld, newleft = %ld, newright = %ld\n", new_prefix_length, newleft, newright);
        fprintf(stderr, "\tcomputed_newleft = %ld, computed_newright = %ld\n", computed_newleft, computed_newright);
        std::exit(EXIT_FAILURE);
      }
    }
  }

  delete[] block_psa;
  delete tail_gt_begin_reversed;
}


void test_random(int testcases, int max_length, int max_sigma) {
  fprintf(stderr,"TEST, testcases = %d, max_n = %d, max_sigma = %d\r",
      testcases, max_length, max_sigma);
  unsigned char *text = new unsigned char[max_length + 1];

  for (int tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {
    // Print progress information.
    if (dbg == 1000) {
      fprintf(stderr,"TEST, testcases = %d, max_n = %d, max_sigma = %d: "
          "%d (%.0Lf%%)\r", testcases, max_length, max_sigma,
          tc, (tc * 100.L) / testcases);
      dbg = 0;
    }

    // XXX can the tail be empty?
    long text_tail_length = utils::random_long(1L, max_length);
    long text_length = utils::random_long(1L, text_tail_length);
    long tail_length = text_tail_length - text_length;
    long block_length = utils::random_long(1L, text_length);
    long block_beg = utils::random_long(0L, text_length - block_length);
    long block_end = block_beg + block_length;
    long sigma = utils::random_long(1L, max_sigma);

    if (sigma <= 26) utils::fill_random_letters(text, text_length, sigma);
    else utils::fill_random_string(text, text_length, sigma);

    // Run the test on generated string.
    test(text, text_length, tail_length, block_beg, block_end);
  }

  // Clean up.
  delete[] text;

  fprintf(stderr,"TEST, testcases = %d, max_n = %d, max_sigma = %d: "
      "\033[22;32mPASSED\033[0m%10s\n", testcases, max_length, max_sigma, "");
}

int main(int, char **) {
  std::srand(std::time(0) + getpid());

  test_random(500000000,  10,      5);
  test_random(500000000,  10,     20);
  test_random(500000000,  10,    256);

  test_random(50000000,  100,     5);
  test_random(50000000,  100,    20);
  test_random(50000000,  100,   256);

  test_random(500000,   1000,    5);
  test_random(500000,   1000,   20);
  test_random(500000,   1000,  256);

  fprintf(stderr, "All tests passed.\n");
}

