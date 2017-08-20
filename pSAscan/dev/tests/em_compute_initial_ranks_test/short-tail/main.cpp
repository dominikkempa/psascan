#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <unistd.h>
#include <string>
#include <algorithm>

#include "divsufsort.h"
#include "utils.hpp"
#include "io/multifile.hpp"
#include "em_compute_initial_ranks.hpp"
#include "io/io_streamer.hpp"
#include "io/multifile_bit_stream_reader.hpp"


using namespace psascan_private;


void test(
   const std::uint8_t *text,
   std::uint64_t block_beg,
   std::uint64_t block_end,
   std::uint64_t tail_beg,
   std::uint64_t text_length) {

  // Write text to disk.
  std::string text_filename = "tempfile" + utils::random_string_hash();
  utils::write_to_file(text, text_length, text_filename);

  // Compute SA of text.
  int *text_sa = new int[text_length];
  divsufsort(text, text_sa, (int)text_length);

  // Compute partial SA of block.
  std::uint64_t block_length = block_end - block_beg;
  int *block_psa = new int[block_length];
  std::uint64_t counter = 0;
  for (std::uint64_t j = 0; j < text_length; ++j)
    if (block_beg <= (std::uint64_t)text_sa[j] &&
        (std::uint64_t)text_sa[j] < block_end)
      block_psa[counter++] = text_sa[j] - block_beg;
  delete[] text_sa;

  const std::uint8_t *block = text + block_beg;

  // Compute gt_begin_reversed for the tail.
  std::uint64_t tail_length = text_length - tail_beg;
  const std::uint8_t *tail = text + tail_beg;

  std::uint8_t *gt_begin_reversed = new std::uint8_t[tail_length];
  for (std::uint64_t j = 1; j <= tail_length; ++j) {
    std::uint64_t lcp = 0;
    while (j + lcp < tail_length && tail[lcp] == tail[j + lcp]) ++lcp;
    gt_begin_reversed[tail_length - j] = (j + lcp < tail_length && tail[lcp] < tail[j + lcp]);
  }
  
  // Write gt_begin_reversed to several files.
  multifile gt_begin_rev_multifile;
  std::uint64_t bits_left = tail_length;
  std::uint64_t bits_beg = 0;
  while (bits_left > 0) {
    std::uint64_t size = utils::random_int64(1L, bits_left);
    std::string filename = "tempfile_gt" + utils::random_string_hash();
    bit_stream_writer *writer = new bit_stream_writer(filename);
    for (std::uint64_t j = 0; j < size; ++j)
      writer->write(gt_begin_reversed[bits_beg + j]);
    delete writer;
    
    gt_begin_rev_multifile.add_file(bits_beg, bits_beg + size, filename);
    bits_beg += size;
    bits_left -= size; 
  }
  delete[] gt_begin_reversed;
  
  // Run the tested algorithm.
  std::uint64_t max_threads = utils::random_int64(1L, 20L);
  std::uint64_t stream_max_block_size = (tail_length + max_threads - 1) / max_threads;
  std::uint64_t n_threads = (tail_length + stream_max_block_size - 1) / stream_max_block_size;
  std::vector<std::uint64_t> result;
  em_compute_initial_ranks(block, block_psa, block_beg, block_end, text_length,
      text_filename, &gt_begin_rev_multifile, result, n_threads, tail_beg);

  // Compare computed answers to correct answers.
  for (std::uint64_t t = 0; t < n_threads; ++t) {
    std::uint64_t stream_block_beg = tail_beg + t * stream_max_block_size;

    const std::uint8_t *pat = text + stream_block_beg;
    std::uint64_t pat_length = text_length - stream_block_beg;

    std::uint64_t srank = 0;
    for (std::uint64_t j = block_beg; j < block_end; ++j) {
      std::uint64_t lcp = 0;
      while (lcp < pat_length && text[j + lcp] == pat[lcp]) ++lcp;
      if (lcp < pat_length && text[j + lcp] < pat[lcp]) ++srank;
    }

    if (srank != result[t]) {
      fprintf(stderr, "\nError!\n");
      fprintf(stderr, "\ttext = ");
      for (std::uint64_t j = 0; j < text_length; ++j)
        fprintf(stderr, "%c", text[j]);
      fprintf(stderr, "\n");
      fprintf(stderr, "\tblock_beg = %lu\n", block_beg);
      fprintf(stderr, "\tblock_end = %lu\n", block_end);
      fprintf(stderr, "\ttail_beg = %lu\n", tail_beg);
      fprintf(stderr, "\tstream_block_max_size = %lu\n", stream_max_block_size);
      fprintf(stderr, "\tn_threads = %lu\n", n_threads);
      fprintf(stderr, "\tt = %lu\n", t);
      fprintf(stderr, "\tcorrect srank = %lu\n", srank);
      fprintf(stderr, "\tcomputed srank = %lu\n", result[t]);
      std::exit(EXIT_FAILURE);
    }
  }

  delete[] block_psa;
  utils::file_delete(text_filename);
}

void test_random(
    std::uint64_t testcases,
    std::uint64_t max_length,
    std::uint64_t max_sigma) {

  fprintf(stderr,"TEST, testcases = %lu, max_n = %lu, max_sigma = %lu\r",
      testcases, max_length, max_sigma);

  std::uint8_t *text = new std::uint8_t[4 * max_length];
  for (std::uint64_t tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {

    // Print progress information.
    if (dbg == 100) {
      fprintf(stderr,"TEST, testcases = %lu, max_n = %lu, max_sigma = %lu: "
          "%lu (%.0Lf%%)\r", testcases, max_length, max_sigma,
          tc, (tc * 100.L) / testcases);
      dbg = 0;
    }

    std::uint64_t mid_length = utils::random_int64(0L, max_length);
    std::uint64_t tail_length = utils::random_int64(1L, max_length);
    std::uint64_t block_length = utils::random_int64(1L, max_length);
    std::uint64_t head_length = utils::random_int64(0L, max_length);
    
    std::uint64_t block_beg   = head_length;
    std::uint64_t block_end   = head_length + block_length;
    std::uint64_t tail_beg    = block_end + mid_length;
    std::uint64_t text_length = head_length + block_length + mid_length + tail_length;

    std::uint64_t sigma = utils::random_int64(1L, max_sigma);

    if (sigma <= 26) utils::fill_random_letters(text, text_length, sigma);
    else utils::fill_random_string(text, text_length, sigma);

    // Run the test on generated string.
    test(text, block_beg, block_end, tail_beg, text_length);
  }

  // Clean up.
  delete[] text;

  fprintf(stderr,"TEST, testcases = %lu, max_n = %lu, max_sigma = %lu: "
      "\033[22;32mPASSED\033[0m%10s\n", testcases, max_length, max_sigma, "");
}


int main() {
  std::srand(std::time(0) + getpid());

  test_random(10000,  10,      5);
  test_random(10000,  10,     20);
  test_random(10000,  10,    256);

  test_random(8000,  100,     5);
  test_random(8000,  100,    20);
  test_random(8000,  100,   256);

  test_random(3000,  300,     5);
  test_random(3000,  300,    20);
  test_random(3000,  300,   256);

  test_random(1000,  1000,    5);
  test_random(1000,  1000,   20);
  test_random(1000,  1000,  256);

  test_random(300,  10000,    5);
  test_random(300,  10000,   20);
  test_random(300,  10000,  256);

  fprintf(stderr, "All tests passed.\n");
}

