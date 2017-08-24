#include <cstdio>
#include <cstring>
#include <cstdint>
#include <ctime>
#include <algorithm>
#include <thread>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "../../src/psascan_src/bitvector.hpp"
#include "../../src/psascan_src/io/multifile.hpp"
#include "inmem_psascan.hpp"
#include "divsufsort.h"
#include "divsufsort64.h"
#include "utils.hpp"
#include "io_streamer.hpp"


namespace psa_test_random_private {

void compute_gt_begin_reversed(
    const std::uint8_t *text,
    std::uint64_t text_length,
    psascan_private::bitvector *gt_begin_reversed) {

  std::uint64_t i = 1, el = 0;
  while (i < text_length) {
    while (i + el < text_length && text[i + el] == text[el]) ++el;
    if (i + el < text_length && text[i + el] > text[el])
      gt_begin_reversed->set(text_length - i);

    el = 0;
    ++i;
  }
}

template<typename saidx_t, unsigned pagesize_log>
void test(
    std::uint8_t *supertext,
    std::uint64_t supertext_length,
    std::uint64_t text_beg,
    std::uint64_t text_end,
    std::uint64_t max_threads) {

  // Sort supertext using divsufsort.
  long *supertext_sa = (long *)malloc(supertext_length * sizeof(long));
  divsufsort64(supertext, supertext_sa, (long)supertext_length);

  // Separate the ordering of the suffixes of text (correct answer).
  std::uint64_t text_length = text_end - text_beg;
  std::uint64_t *correct_answer =
    (std::uint64_t *)malloc(text_length * sizeof(std::uint64_t));
  std::uint64_t ptr = 0;
  for (std::uint64_t i = 0; i < supertext_length; ++i)
    if (text_beg <= (std::uint64_t)supertext_sa[i] &&
        (std::uint64_t)supertext_sa[i] < text_end)
      correct_answer[ptr++] = supertext_sa[i] - text_beg;

  // Compute tail_gt_begin_reversed.
  const std::uint8_t *tail = supertext + text_end;
  std::uint64_t tail_length = supertext_length - text_end;
  psascan_private::bitvector tail_gt_begin_reversed_bv(tail_length);
  compute_gt_begin_reversed(tail, tail_length, &tail_gt_begin_reversed_bv);

  // Store tail_gt_begin_reversed on disk as a multifile bitvector.
  psascan_private::multifile *tail_gt_begin_reversed_multifile =
    new psascan_private::multifile();
  ptr = 0;
  while (ptr < tail_length) {
    std::uint64_t left = tail_length - ptr;
    std::uint64_t chunk = utils::random_int64(1L, left);
   
    // Store bits [ptr..ptr + chunk) from
    // tail_gt_begin_reversed_bv into one file.
    std::string chunk_filename =
      "gt_begin_reversed_bv" + utils::random_string_hash();
    bit_stream_writer *writer = new bit_stream_writer(chunk_filename);
    for (std::uint64_t j = ptr; j < ptr + chunk; ++j)
      writer->write(tail_gt_begin_reversed_bv.get(j));
    delete writer;

    // Add this file to tail_gt_begin_reversed_multifile.
    tail_gt_begin_reversed_multifile->add_file(ptr,
        ptr + chunk, chunk_filename);
    
    ptr += chunk;
  }

  // Write supertext to file.
  std::string supertext_filename = "supertext.txt";
  stream_writer<std::uint8_t> *supertext_writer =
    new stream_writer<std::uint8_t>(supertext_filename);
  for (std::uint64_t i = 0; i < supertext_length; ++i)
    supertext_writer->write(supertext[i]);
  delete supertext_writer;






  // Run the tested algorithm.
  std::uint8_t *text = supertext + text_beg;
  std::uint8_t *bwtsa =
    (std::uint8_t *)malloc(text_length * (1 + sizeof(saidx_t)));
  saidx_t *computed_sa = (saidx_t *)bwtsa;
  std::uint64_t max_blocks = 0; // 0 is setting max_blocks := max_threads
  if (utils::random_int64(0, 1)) max_blocks = utils::random_int64(1L, 50L);
  bool compute_bwt = (bool)utils::random_int64(0L, 1L);
  std::uint64_t tail_prefix_length = std::min(text_length, tail_length);
  std::uint8_t *tail_prefix = NULL;
  if (utils::random_int64(0L, 1L) && tail_length > 0) {
    tail_prefix = (std::uint8_t *)malloc(tail_prefix_length);
    std::copy(tail, tail + tail_prefix_length, tail_prefix);
  }

  inmem_psascan<saidx_t, pagesize_log>(text, text_length, bwtsa,
      max_threads, compute_bwt, false, NULL, max_blocks, text_beg,
      text_end, supertext_length, supertext_filename,
      tail_gt_begin_reversed_multifile, NULL, tail_prefix);

  // Compare answers.
  bool eq = true;
  for (std::uint64_t i = 0; i < text_length; ++i)
    if ((std::uint64_t)computed_sa[i] != correct_answer[i]) eq = false;
  if (!eq) {
    fprintf(stdout, "Error:\n");
    fprintf(stdout, "\tsupertext_length = %lu\n", supertext_length);
    fprintf(stdout, "\tmax threads = %lu\n", max_threads);
    fprintf(stdout, "\tsupertext = ");
    for (std::uint64_t j = 0; j < supertext_length; ++j)
      fprintf(stdout, "%c", supertext[j]);
    fprintf(stdout, "\n");
    fprintf(stdout, "\ttext_beg = %lu\n", text_beg);
    fprintf(stdout, "\ttext_end = %lu\n", text_end);
    fprintf(stdout, "\ttail_gt_begin_reversed_bv = ");
    for (std::uint64_t j = 0; j < tail_length; ++j)
      fprintf(stdout, "%lu", (std::uint64_t)tail_gt_begin_reversed_bv.get(j));
    fprintf(stdout, "\n");
    fprintf(stdout, "\tsupertext sa = ");
    for (std::uint64_t j = 0; j < supertext_length; ++j)
      fprintf(stdout, "%lu ", (std::uint64_t)supertext_sa[j]);
    fprintf(stdout, "\n");
    fprintf(stdout, "\tcorrect answer = ");
    for (std::uint64_t j = 0; j < text_length; ++j)
      fprintf(stdout, "%lu ", (std::uint64_t)correct_answer[j]);
    fprintf(stdout, "\n");
    fprintf(stdout, "\tcomputed answer = ");
    for (std::uint64_t j = 0; j < text_length; ++j)
      fprintf(stdout, "%lu ", (std::uint64_t)computed_sa[j]);
    fprintf(stdout, "\n");
    std::fflush(stdout);

    std::exit(EXIT_FAILURE);
  }


  delete tail_gt_begin_reversed_multifile; // also deletes files
  free(bwtsa);
  free(correct_answer);
  free(supertext_sa);
}


template<typename saidx_t, unsigned pagesize_log>
void test_random(
    int testcases,
    std::uint64_t max_length,
    int max_sigma) {

  fprintf(stdout,"TEST, testcases = %d, max_n = %lu, "
      "max_sigma = %d, sizeof(saidx_t) = %lu, pagesize_log = %lu\n",
      testcases, max_length, max_sigma,
      (std::uint64_t)sizeof(saidx_t), (std::uint64_t)pagesize_log);
  std::uint8_t *supertext = new std::uint8_t[max_length + 1];

  for (int tc = 0; tc < testcases; ++tc) {
    // Print progress information.
    fprintf(stdout,"%d (%.2Lf%%)\r", tc, (tc * 100.L) / testcases);
    std::fflush(stdout);

    // Generate string.
    std::uint64_t supertext_length = utils::random_int64(1, max_length);
    int sigma = utils::random_int32(1, max_sigma);
    if (max_sigma <= 26)
      utils::fill_random_letters(supertext, supertext_length, sigma);
    else utils::fill_random_string(supertext, supertext_length, sigma);
    std::uint64_t max_threads = utils::random_int64(1, 50);
    std::uint64_t text_beg = utils::random_int64(0, supertext_length - 1);
    std::uint64_t text_end =
      utils::random_int64(text_beg + 1, supertext_length);



    /*std::uint64_t supertext_length = 325;
    std::uint64_t max_threads = 16;
    supertext[0] = 0;
    strcpy((char *)supertext, "adbbddacdcdacacccddcbdccbaaaabcaadcadcccabccbaaddabadadbadbaadabdcbcadaaacdbcdbbdccdcbacabcaadbdbcbcbbcbdbdbbadacbdbcddcaccbbaacccaaddbdaaabadcdabacbdabbdccddbbbbbaaddbdacadadacdcdcdaacdcbcdcdbadbddccdbccbbcdabcdaddccbdabbdcbcdabcdadbdadbadccccbbaddddabcccbbdcdcdcdcdcddccbaadcbcbacbbadabadaabdcabbdaabccbdbdadaabbccacdbbdbcc");
    std::uint64_t text_beg = 125;
    std::uint64_t text_end = 269;*/


    // Run the test on generated string.
    test<saidx_t, pagesize_log>(supertext, supertext_length,
        text_beg, text_end, max_threads);
  }

  // Clean up.
  delete[] supertext;
}

}  // namespace psa_test_random_private
