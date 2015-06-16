#include <cstdio>
#include <cstdlib>
#include <ctime>

#include <unistd.h>

#include <string>
#include <algorithm>

#include "divsufsort.h"
#include "utils.h"
#include "multifile.h"
#include "em_compute_initial_ranks.h"
#include "io_streamer.h"
#include "multifile_bit_stream_reader.h"


void test(unsigned char *text, long block_beg, long block_end, long text_length) {
  // Write text to disk.
  std::string text_filename = "tempfile" + utils::random_string_hash();
  utils::write_objects_to_file(text, text_length, text_filename);

  // Compute SA of text.
  int *text_sa = new int[text_length];
  divsufsort(text, text_sa, (int)text_length);

  // Compute partial SA of block.
  long block_length = block_end - block_beg;
  int *block_psa = new int[block_length];
  long counter = 0;
  for (long j = 0; j < text_length; ++j)
    if (block_beg <= text_sa[j] && text_sa[j] < block_end)
      block_psa[counter++] = text_sa[j] - block_beg;
  delete[] text_sa;

  long i0 = 0;
  unsigned char *block = text + block_beg;
  unsigned char *block_pbwt = new unsigned char[block_length];
  for (long j = 0; j < block_length; ++j) {
    if (block_psa[j]) block_pbwt[j] = block[block_psa[j] - 1];
    else { i0 = j; block_pbwt[j] = 0; }
  }

  // Compute gt_begin_reversed for the tail.
  long tail_length = text_length - block_end;
  unsigned char *tail = text + block_end;

  unsigned char *gt_begin_reversed = new unsigned char[tail_length];
  for (long j = 1; j <= tail_length; ++j) {
    long lcp = 0;
    while (j + lcp < tail_length && tail[lcp] == tail[j + lcp]) ++lcp;
    gt_begin_reversed[tail_length - j] = (j + lcp < tail_length && tail[lcp] < tail[j + lcp]);
  }
  
  // Write gt_begin_reversed to several files.
  multifile gt_begin_rev_multifile;
  long bits_left = tail_length;
  long bits_beg = 0;
  while (bits_left > 0) {
    long size = utils::random_long(1L, bits_left);
    std::string filename = "tempfile_gt" + utils::random_string_hash();
    bit_stream_writer *writer = new bit_stream_writer(filename);
    for (long j = 0; j < size; ++j)
      writer->write(gt_begin_reversed[bits_beg + j]);
    delete writer;
    
    gt_begin_rev_multifile.add_file(bits_beg, bits_beg + size, filename);
    bits_beg += size;
    bits_left -= size; 
  }
  delete[] gt_begin_reversed;
  
  // Run the tested algorithm.
  long max_threads = utils::random_long(1L, 20L);
  long stream_max_block_size = (tail_length + max_threads - 1) / max_threads;
  long n_threads = (tail_length + stream_max_block_size - 1) / stream_max_block_size;
  std::vector<long> result;
  em_compute_initial_ranks(block, block_psa, block_pbwt, i0, block_beg, block_end,
      text_length, text_filename, &gt_begin_rev_multifile, result, n_threads);

  // Compare computed answers to correct answers.
  for (long t = 0; t < n_threads; ++t) {
    long stream_block_beg = block_end + t * stream_max_block_size;
    long stream_block_end = std::min(stream_block_beg + stream_max_block_size, text_length);

    unsigned char *pat = text + stream_block_end;
    long pat_length = text_length - stream_block_end;

    long srank = 0;
    for (long j = block_beg; j < block_end; ++j) {
      long lcp = 0;
      while (lcp < pat_length && text[j + lcp] == pat[lcp]) ++lcp;
      if (lcp < pat_length && text[j + lcp] < pat[lcp]) ++srank;
    }

    if (srank != result[t]) {
      fprintf(stderr, "\nError!\n");
      fprintf(stderr, "\ttext = ");
      for (long j = 0; j < text_length; ++j)
        fprintf(stderr, "%c", text[j]);
      fprintf(stderr, "\n");
      fprintf(stderr, "\tblock_beg = %ld\n", block_beg);
      fprintf(stderr, "\tblock_end = %ld\n", block_end);
      fprintf(stderr, "\tstream_block_max_size = %ld\n", stream_max_block_size);
      fprintf(stderr, "\tn_threads = %ld\n", n_threads);
      fprintf(stderr, "\tt = %ld\n", t);
      fprintf(stderr, "\tcorrect srank = %ld\n", srank);
      fprintf(stderr, "\tcomputed srank = %ld\n", result[t]);
      std::exit(EXIT_FAILURE);
    }
  }

  delete[] block_psa;
  delete[] block_pbwt;
  utils::file_delete(text_filename);
}

void test_random(int testcases, int max_length, int max_sigma) {
  fprintf(stderr,"TEST, testcases = %d, max_n = %d, max_sigma = %d\r",
      testcases, max_length, max_sigma);

  unsigned char *text = new unsigned char[3 * max_length];
  for (int tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {
    // Print progress information.
    if (dbg == 100) {
      fprintf(stderr,"TEST, testcases = %d, max_n = %d, max_sigma = %d: "
          "%d (%.0Lf%%)\r", testcases, max_length, max_sigma,
          tc, (tc * 100.L) / testcases);
      dbg = 0;
    }

    long tail_length = utils::random_long(1L, max_length);
    long block_length = utils::random_long(1L, max_length);
    long head_length = utils::random_long(0L, max_length);
    
    long block_beg   = head_length;
    long block_end   = head_length + block_length;
    long text_length = head_length + block_length + tail_length;

    long sigma = utils::random_long(1L, max_sigma);

    if (sigma <= 26) utils::fill_random_letters(text, text_length, sigma);
    else utils::fill_random_string(text, text_length, sigma);

    // Run the test on generated string.
    test(text, block_beg, block_end, text_length);
  }

  // Clean up.
  delete[] text;

  fprintf(stderr,"TEST, testcases = %d, max_n = %d, max_sigma = %d: "
      "\033[22;32mPASSED\033[0m%10s\n", testcases, max_length, max_sigma, "");
}


int main() {
  std::srand(std::time(0) + getpid());

  test_random(50000,  10,      5);
  test_random(50000,  10,     20);
  test_random(50000,  10,    256);

  test_random(50000,  100,     5);
  test_random(50000,  100,    20);
  test_random(50000,  100,   256);

  test_random(30000,  300,     5);
  test_random(30000,  300,    20);
  test_random(30000,  300,   256);

  test_random(50000,  1000,    5);
  test_random(50000,  1000,   20);
  test_random(50000,  1000,  256);

  test_random(5000,  10000,    5);
  test_random(5000,  10000,   20);
  test_random(5000,  10000,  256);

  test_random(500,  100000,    5);
  test_random(500,  100000,   20);
  test_random(500,  100000,  256);


  fprintf(stderr, "All tests passed.\n");
}

