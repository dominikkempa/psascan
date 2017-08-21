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
#include "divsufsort.h"
#include "divsufsort64.h"
#include "utils.hpp"
#include "io_streamer.hpp"
#include "divsufsort_template.hpp"


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
    run_divsufsort(supertext, correct_sa, (long)supertext_length);
    utils::write_to_file(correct_sa, supertext_length, sa_filename);
    delete[] correct_sa;
    fprintf(stderr, "Total time: %.2Lf\n", utils::wclock() - start);
  }

  text_length = std::min(text_length, supertext_length);
  std::uint64_t text_beg =
    utils::random_int64(0L, supertext_length - text_length);
  std::uint64_t text_end = text_beg + text_length;


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
   
    // Store bits [ptr..ptr + chunk) from
    // tail_gt_begin_reversed_bv into one file.
    std::string chunk_filename =
      "gt_begin_reversed_bv" + utils::random_string_hash();
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
  saidx_t *computed_sa = (saidx_t *)bwtsa;
  inmem_psascan<saidx_t, pagesize_log>(text, text_length,
      bwtsa, max_threads, false, false, NULL, 0, text_beg,
      text_end, supertext_length, supertext_filename,
      tail_gt_begin_reversed_multifile);

  ptr = 0;
  fprintf(stderr, "\nComparing:\n");
  stream_reader<std::uint64_t> *sa_reader =
    new stream_reader<std::uint64_t>(sa_filename);
  bool eq = true;
  std::uint64_t compared = 0;
  for (std::uint64_t i = 0, dbg = 0; i < supertext_length; ++i) {
    ++dbg;
    ++compared;
    if (dbg == 10000000) {
      fprintf(stderr, "progress: %.3Lf%%\r", (100.L * i) / supertext_length);
      dbg = 0;
    }

    std::uint64_t next_correct_sa = sa_reader->read();
    if (text_beg <= next_correct_sa && next_correct_sa < text_end) {
      if (next_correct_sa - text_beg != (std::uint64_t)computed_sa[ptr++]) {
        eq = false;
        break;
      }
    }
  }
  fprintf(stderr, "Compared %lu values", compared);
  fprintf(stderr, "\nResult: %s\n", eq ? "OK" : "FAIL");

  free(bwtsa);
  free(text);
  delete sa_reader;
  delete tail_gt_begin_reversed_multifile; // also deletes files
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

