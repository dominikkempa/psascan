#include <cstdio>
#include <string>
#include <cstring>
#include <ctime>
#include <unistd.h>

#include "utils.h"
#include "inmem_sascan.h"
#include "divsufsort_template.h"
#include "uint40.h"
#include "io_streamer.h"


template<typename saidx_t>
void read_sa(saidx_t* &sa, std::string filename) {
  long length;
  utils::read_objects_from_file(sa, length, filename);
}


template<typename saidx_t>
void test(unsigned char *text, long text_length, long max_threads,
    long max_blocks, std::string filename) {
  long double start;

  std::string sa_filename = filename + ".sa" + utils::intToStr(sizeof(long));
  if (!utils::file_exists(sa_filename)) {
    fprintf(stderr, "Running divsufsort\n");
    start = utils::wclock();
    long *correct_sa = new long[text_length];
    run_divsufsort(text, correct_sa, text_length);
    utils::write_objects_to_file(correct_sa, text_length, sa_filename);
    delete[] correct_sa;
    fprintf(stderr, "Total time: %.2Lf\n", utils::wclock() - start);
  }

  fprintf(stderr, "Running inmem sascan\n\n");
  unsigned char *computed_sa_temp = (unsigned char *)malloc(text_length * (sizeof(saidx_t) + 1));
  saidx_t *computed_sa = (saidx_t *)computed_sa_temp;
  unsigned char *computed_bwt = (unsigned char *)(computed_sa + text_length);
  start = utils::wclock();
  inmem_sascan<saidx_t>(text, text_length, computed_sa_temp, max_threads, true, max_blocks);
  long double total_time = utils::wclock() - start;
  fprintf(stderr, "\nTotal time:\n");
  fprintf(stderr, "\tabsolute: %.2Lf\n", total_time);
  fprintf(stderr, "\trelative: %.4Lfs/MiB\n", total_time / ((long double)text_length / (1 << 20)));
  fprintf(stderr, "Speed: %.2LfMiB/s\n", ((long double)text_length / (1 << 20)) / total_time);

  fprintf(stderr, "\nComparing:\n");
  stream_reader<long> *sa_reader = new stream_reader<long>(sa_filename);
  bool eq = true;
  long compared = 0;
  for (long i = 0, dbg = 0; i < text_length; ++i) {
    ++dbg;
    ++compared;
    if (dbg == 10000000) {
      fprintf(stderr, "progress: %.3Lf%%\r", (100.L * i) / text_length);
      dbg = 0;
    }

    long next_correct_sa = sa_reader->read();
    unsigned char next_correct_bwt = ((next_correct_sa == 0) ? 0 : text[next_correct_sa - 1]);
    if (next_correct_bwt != computed_bwt[i]) {
      eq = false;
      break;
    }
  }
  fprintf(stderr, "Compared %ld values", compared);
  fprintf(stderr, "\nResult: %s\n", eq ? "OK" : "FAIL");

  free(computed_sa_temp);
  delete sa_reader;
}


void test_file(const char *filename) {
  fprintf(stderr, "Input filename: %s\n", filename);
  fprintf(stderr, "Reading text: ");
  long length;
  unsigned char *text;
  utils::read_objects_from_file(text, length, filename);
  fprintf(stderr, "DONE\n");

  // test<uint40>(text, length, 24, 8, filename);
  // test<uint40>(text, length, 24, 12, filename);
  test<uint40>(text, length, 24, 16, filename);
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
  for (long i = 0; i < argc; ++i)
    fprintf(stderr, " %s", argv[i]);
  fprintf(stderr, "\n");

  test_file(argv[1]);
}
