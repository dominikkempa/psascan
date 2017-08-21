#include <cstdio>
#include <string>
#include <cstring>
#include <ctime>
#include <unistd.h>

#include "inmem_psascan.hpp"
#include "utils.hpp"
#include "divsufsort_template.hpp"
#include "uint40.hpp"
#include "io_streamer.hpp"


template<typename saidx_t>
void test(
    std::uint8_t *text,
    std::uint64_t text_length,
    std::uint64_t max_threads,
    std::uint64_t max_blocks,
    std::string filename) {

  long double start;

  std::string sa_filename = filename + ".sa" +
    utils::intToStr(sizeof(std::uint64_t));
  if (!utils::file_exists(sa_filename)) {
    fprintf(stderr, "Running divsufsort\n");
    start = utils::wclock();
    long *correct_sa = new long[text_length];
    run_divsufsort(text, correct_sa, (long)text_length);
    utils::write_to_file(correct_sa, text_length, sa_filename);
    delete[] correct_sa;
    fprintf(stderr, "Total time: %.2Lf\n", utils::wclock() - start);
  }

  fprintf(stderr, "Running inmem sascan\n\n");
  std::uint8_t *computed_sa_temp =
    (std::uint8_t *)malloc(text_length * (sizeof(saidx_t) + 1));
  saidx_t *computed_sa = (saidx_t *)computed_sa_temp;
  std::uint8_t *computed_bwt = (std::uint8_t *)(computed_sa + text_length);
  std::uint64_t computed_i0;
  inmem_psascan<saidx_t>(text, text_length, computed_sa_temp,
      max_threads, true, false, NULL, max_blocks, 0, 0, 0,
      "", NULL, &computed_i0);

  fprintf(stderr, "\nComparing:\n");
  stream_reader<std::uint64_t> *sa_reader =
    new stream_reader<std::uint64_t>(sa_filename);
  bool eq = true;
  std::uint64_t compared = 0;
  std::uint64_t correct_i0 = 0;
  for (std::uint64_t i = 0, dbg = 0; i < text_length; ++i) {
    ++dbg;
    ++compared;
    if (dbg == 10000000) {
      fprintf(stderr, "progress: %.3Lf%%\r", (100.L * i) / text_length);
      dbg = 0;
    }

    std::uint64_t next_correct_sa = sa_reader->read();
    std::uint8_t next_correct_bwt =
      ((next_correct_sa == 0) ? 0 : text[next_correct_sa - 1]);
    if (next_correct_bwt != computed_bwt[i]) {
      eq = false;
      break;
    }
    if (next_correct_sa == 0) correct_i0 = i;
  }
  if (correct_i0 != computed_i0) eq = false;
  fprintf(stderr, "Compared %lu values", compared);
  fprintf(stderr, "\nResult: %s\n", eq ? "OK" : "FAIL");

  free(computed_sa_temp);
  delete sa_reader;
}


void test_file(const char *text_filename) {
  fprintf(stderr, "Input filename: %s\n", text_filename);
  fprintf(stderr, "Reading text: ");
  std::uint64_t text_length = utils::file_size(text_filename);
  std::uint8_t *text = new std::uint8_t[text_length];
  utils::read_from_file(text, text_length, text_filename);
  fprintf(stderr, "DONE\n");

  // test<uint40>(text, length, 24, 8, filename);
  // test<uint40>(text, length, 24, 12, filename);
  test<uint40>(text, text_length, 24, 16, text_filename);
  // test<uint40>(text, length, 24, 24, filename);
  // test<uint40>(text, length, 24, 32, filename);

  delete[] text;
}

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s <file>\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }

  fprintf(stderr, "Command line:");
  for (int i = 0; i < argc; ++i)
    fprintf(stderr, " %s", argv[i]);
  fprintf(stderr, "\n");

  test_file(argv[1]);
}
