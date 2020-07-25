#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <unistd.h>
#include <string>
#include <algorithm>

#include "divsufsort.h"
#include "utils.hpp"
#include "io/multifile.hpp"
#include "compute_ranks.hpp"
#include "io/io_streamer.hpp"
#include "io/multifile_bit_stream_reader.hpp"


using namespace psascan_private;


void test(
    const std::uint8_t * const text,
    const std::uint64_t block_beg,
    const std::uint64_t block_end,
    const std::uint64_t tail_end,
    const std::uint64_t text_length) {

  // Write text to disk.
  std::string text_filename = "tempfile" + utils::random_string_hash();
  utils::write_to_file(text, text_length, text_filename);

  // Compute SA of text.
  int * const text_sa = new int[text_length];
  divsufsort(text, text_sa, (int)text_length);

  // Compute partial SA of block.
  const std::uint64_t block_length = block_end - block_beg;
  int *const block_psa = new int[block_length];
  std::uint64_t counter = 0;
  for (std::uint64_t j = 0; j < text_length; ++j)
    if (block_beg <= (std::uint64_t)text_sa[j] &&
        (std::uint64_t)text_sa[j] < block_end)
      block_psa[counter++] = text_sa[j] - block_beg;
  delete[] text_sa;

  std::uint64_t i0 = 0;
  const std::uint8_t * const block = text + block_beg;
  std::uint8_t * const block_pbwt = new std::uint8_t[block_length];
  for (std::uint64_t j = 0; j < block_length; ++j) {
    if (block_psa[j]) block_pbwt[j] = block[block_psa[j] - 1];
    else { i0 = j; block_pbwt[j] = 0; }
  }

  // Compute gt_begin_reversed for the tail.
  const std::uint64_t tail_length = tail_end - block_end;
  const std::uint8_t *pat = text + block_end;
  std::uint64_t pat_length = text_length - block_end;

  std::uint8_t * const gt_begin_reversed =
    new std::uint8_t[text_length - block_end];
  for (std::uint64_t j = 1; j <= tail_length; ++j) {
    std::uint64_t lcp = 0;
    while (j + lcp < pat_length && pat[lcp] == pat[j + lcp]) ++lcp;
    gt_begin_reversed[(text_length - tail_end) + (tail_length - j)] =
      (j + lcp < pat_length && pat[lcp] < pat[j + lcp]);
  }
  
  // Write gt_begin_reversed to several files.
  // We only write range of bits
  // [text_length - tail_end..text_length - block_end).
  multifile gt_begin_rev_multifile;
  std::uint64_t bits_left = tail_length;
  std::uint64_t bits_beg = text_length - tail_end;
  while (bits_left > 0) {
    const std::uint64_t size = utils::random_int64(1L, bits_left);
    const std::string filename = "tempfile_gt" + utils::random_string_hash();
    bit_stream_writer * const writer = new bit_stream_writer(filename);
    for (std::uint64_t j = 0; j < size; ++j)
      writer->write(gt_begin_reversed[bits_beg + j]);
    delete writer;
    
    gt_begin_rev_multifile.add_file(bits_beg, bits_beg + size, filename);
    bits_beg += size;
    bits_left -= size; 
  }
  delete[] gt_begin_reversed;

  // Compute srank for suffix starting after tail.
  pat = text + tail_end;
  pat_length = text_length - tail_end;
  std::uint64_t srank_after_tail = 0;
  for (std::uint64_t j = block_beg; j < block_end; ++j) {
    std::uint64_t lcp = 0;
    while (lcp < pat_length && text[j + lcp] == pat[lcp]) ++lcp;
    if (lcp < pat_length && text[j + lcp] < pat[lcp]) ++srank_after_tail;
  }
  
  // Run the tested algorithm.
  const std::uint64_t max_threads = utils::random_int64(1L, 20L);
  const std::uint64_t stream_max_block_size =
    (tail_length + max_threads - 1) / max_threads;
  const std::uint64_t n_threads =
    (tail_length + stream_max_block_size - 1) / stream_max_block_size;
  std::vector<std::uint64_t> result;
  compute_ranks(
      block_beg, block_end, text_length,
      stream_max_block_size, tail_end,
      srank_after_tail, i0, block, block_pbwt,
      block_psa, &gt_begin_rev_multifile,
      text_filename, result);

  // Compare computed answers to correct answers.
  for (std::uint64_t t = 0; t < n_threads; ++t) {
    const std::uint64_t stream_block_beg =
      block_end + t * stream_max_block_size;

    pat = text + stream_block_beg;
    pat_length = text_length - stream_block_beg;

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
      fprintf(stderr, "\tblock_beg = %ld\n", block_beg);
      fprintf(stderr, "\tblock_end = %ld\n", block_end);
      fprintf(stderr, "\tstream_block_max_size = %ld\n",
          stream_max_block_size);
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

void test_random(
    const std::uint64_t testcases,
    const std::uint64_t max_length,
    const std::uint64_t max_sigma) {

  fprintf(stderr,"TEST, testcases = %lu, max_n = %lu, max_sigma = %lu\r",
      testcases, max_length, max_sigma);

  std::uint8_t * const text = new std::uint8_t[4 * max_length];
  for (std::uint64_t tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {

    // Print progress information.
    if (dbg == 100) {
      fprintf(stderr,"TEST, testcases = %lu, max_n = %lu, max_sigma = %lu: "
          "%lu (%.0Lf%%)\r", testcases, max_length, max_sigma,
          tc, (tc * 100.L) / testcases);
      dbg = 0;
    }

    const std::uint64_t end_block_length = utils::random_int64(0L, max_length);
    const std::uint64_t tail_length = utils::random_int64(1L, max_length);
    const std::uint64_t block_length = utils::random_int64(1L, max_length);
    const std::uint64_t head_length = utils::random_int64(0L, max_length);
    
    const std::uint64_t block_beg   = head_length;
    const std::uint64_t block_end   = head_length + block_length;
    const std::uint64_t tail_end    = block_end + tail_length;
    const std::uint64_t text_length =
      head_length + block_length + tail_length + end_block_length;

    const std::uint64_t sigma = utils::random_int64(1L, max_sigma);

    if (sigma <= 26) utils::fill_random_letters(text, text_length, sigma);
    else utils::fill_random_string(text, text_length, sigma);

    // Run the test on generated string.
    test(text, block_beg, block_end, tail_end, text_length);
  }

  // Clean up.
  delete[] text;

  fprintf(stderr,"TEST, testcases = %lu, max_n = %lu, max_sigma = %lu: "
      "\033[22;32mPASSED\033[0m%10s\n", testcases, max_length, max_sigma, "");
}


int main() {
  std::srand(std::time(0) + getpid());

#ifdef NDEBUG
  test_random(1000,  10,      5);
  test_random(1000,  10,     20);
  test_random(1000,  10,    256);

  test_random(800,  100,     5);
  test_random(800,  100,    20);
  test_random(800,  100,   256);

  test_random(300,  300,     5);
  test_random(300,  300,    20);
  test_random(300,  300,   256);

  test_random(100,  1000,    5);
  test_random(100,  1000,   20);
  test_random(100,  1000,  256);

  test_random(30,  10000,    5);
  test_random(30,  10000,   20);
  test_random(30,  10000,  256);
#else
  test_random(100,  10,      5);
  test_random(100,  10,     20);
  test_random(100,  10,    256);

  test_random(80,  100,     5);
  test_random(80,  100,    20);
  test_random(80,  100,   256);

  test_random(30,  300,     5);
  test_random(30,  300,    20);
  test_random(30,  300,   256);

  test_random(10,  1000,    5);
  test_random(10,  1000,   20);
  test_random(10,  1000,  256);

  test_random(5,  10000,    5);
  test_random(5,  10000,   20);
  test_random(5,  10000,  256);
#endif

  fprintf(stderr, "All tests passed.\n");
}

