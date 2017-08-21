#include <cstdio>
#include <cstring>
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
#include "utils.hpp"
#include "io_streamer.hpp"
#include "divsufsort_template.hpp"



void compute_gt_begin_for_text(
    const std::uint8_t *supertext,
    std::uint64_t supertext_length,
    std::uint64_t text_beg,
    std::uint64_t text_end,
    psascan_private::bitvector *text_gt_begin) {

  for (std::uint64_t i = text_beg + 1; i <= text_end; ++i) {
    std::uint64_t lcp = 0;
    while (i + lcp < supertext_length &&
        supertext[i + lcp] == supertext[text_beg + lcp]) ++lcp;
    
    if (i + lcp < supertext_length &&
        supertext[i + lcp] > supertext[text_beg + lcp])
      text_gt_begin->set(i - text_beg - 1);
  }
}

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
    std::string supertext_filename,
    std::uint64_t text_length,
    std::uint64_t max_threads) {

  long double start;

  fprintf(stderr, "Input filename: %s\n", supertext_filename.c_str());
  fprintf(stderr, "Reading text: ");
  std::uint64_t supertext_length = utils::file_size(supertext_filename);
  std::uint8_t *supertext = new std::uint8_t[supertext_length];
  utils::read_from_file(supertext, supertext_length, supertext_filename);
  fprintf(stderr, "DONE\n");

  std::string sa_filename = supertext_filename +
    ".sa" + utils::intToStr(sizeof(std::uint64_t));
  if (!utils::file_exists(sa_filename)) {
    fprintf(stderr, "Running divsufsort\n");
    start = utils::wclock();
    long *correct_sa = new long[supertext_length];
    run_divsufsort(supertext, correct_sa, (std::int64_t)supertext_length);
    utils::write_to_file(correct_sa, supertext_length, sa_filename);
    delete[] correct_sa;
    fprintf(stderr, "Total time: %.2Lf\n", utils::wclock() - start);
  }

  text_length = std::min(text_length, supertext_length);
  std::uint64_t text_beg = utils::random_int64((std::int64_t)0,
      supertext_length - text_length);
  std::uint64_t text_end = text_beg + text_length;

  psascan_private::bitvector *text_gt_begin_correct =
    new psascan_private::bitvector(text_length);
  compute_gt_begin_for_text(supertext, supertext_length,
      text_beg, text_end, text_gt_begin_correct);

  // Compute tail_gt_begin_reversed.
  const std::uint8_t *tail = supertext + text_end;
  std::uint64_t tail_length = supertext_length - text_end;
  psascan_private::bitvector *tail_gt_begin_reversed_bv =
    new psascan_private::bitvector(tail_length);
  compute_gt_begin_reversed(tail, tail_length, tail_gt_begin_reversed_bv);

  // Store tail_gt_begin_reversed on disk as a multifile bitvector.
  psascan_private::multifile *tail_gt_begin_reversed_multifile =
    new psascan_private::multifile();
  std::uint64_t ptr = 0;
  while (ptr < tail_length) {
    std::uint64_t left = tail_length - ptr;
    std::uint64_t chunk = utils::random_int64(1L, left);
   
    // Store bits [ptr..ptr+chunk) from
    // tail_gt_begin_reversed_bv into one file.
    std::string chunk_filename = "gt_begin_reversed_bv" +
      utils::random_string_hash();
    bit_stream_writer *writer = new bit_stream_writer(chunk_filename);
    for (std::uint64_t j = ptr; j < ptr + chunk; ++j)
      writer->write(tail_gt_begin_reversed_bv->get(j));
    delete writer;

    // Add this file to tail_gt_begin_reversed_multifile.
    tail_gt_begin_reversed_multifile->add_file(ptr,
        ptr + chunk, chunk_filename);
    
    ptr += chunk;
  }
  delete tail_gt_begin_reversed_bv;

  std::uint8_t *text = (std::uint8_t *)malloc(text_length);
  std::copy(supertext + text_beg, supertext + text_end, text);
  delete[] supertext;

  // Run the tested algorithm.
  fprintf(stderr, "Running inmem sascan\n\n");
  std::uint8_t *bwtsa =
    (std::uint8_t *)malloc(text_length * (1 + sizeof(saidx_t)));
  psascan_private::bitvector *text_gt_begin_computed =
    new psascan_private::bitvector(text_length);
  inmem_psascan<saidx_t, pagesize_log>(text, text_length, bwtsa,
      max_threads, false, // try also true here, to see if it works for both.
      true, text_gt_begin_computed, -1, text_beg, text_end, supertext_length,
      supertext_filename, tail_gt_begin_reversed_multifile);


  // Compute correct answer.
  ptr = 0;
  fprintf(stderr, "\nComparing:\n");
  bool eq = true;
  std::uint64_t compared = 0;
  for (std::uint64_t i = 0, dbg = 0; i < text_length; ++i) {
    ++dbg;
    ++compared;
    if (dbg == 10000000) {
      fprintf(stderr, "progress: %.3Lf%%\r", (100.L * i) / text_length);
      dbg = 0;
    }
    
    bool next_text_gt_begin_correct = text_gt_begin_correct->get(i);
    bool next_text_gt_begin_computed =
      text_gt_begin_computed->get(text_length - 1 - i);

    if (next_text_gt_begin_correct != next_text_gt_begin_computed) {
      eq = false;
      break;
    }
  }
  fprintf(stderr, "Compared %ld values", compared);
  fprintf(stderr, "\nResult: %s\n", eq ? "OK" : "FAIL");

  free(bwtsa);
  free(text);
  delete tail_gt_begin_reversed_multifile; // also deletes files
  delete text_gt_begin_computed;
  delete text_gt_begin_correct;
}


int main(int argc, char **argv) {
  std::srand(std::time(0) + getpid());
  if (argc == 1) {
    fprintf(stderr, "%s <file> <min-text-length-in-MiB>\n",
        argv[0]);
    std::exit(EXIT_FAILURE);
  }

  std::uint64_t min_text_length = atol(argv[2]) << 20;
  test<uint40/*int*/, 12>(argv[1], min_text_length, 24);
}

