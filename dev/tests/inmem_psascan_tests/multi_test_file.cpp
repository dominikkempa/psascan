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
long double test(
    std::uint8_t *text,
    std::uint64_t text_length,
    std::uint64_t max_threads,
    std::uint64_t max_blocks,
    std::string filename,
    bool compare = true) {

  long double start;
  std::string sa_filename = filename + ".sa" +
    utils::intToStr(sizeof(std::uint64_t));
  if (compare) {
    if (!utils::file_exists(sa_filename)) {
      fprintf(stderr, "Running divsufsort\n");
      start = utils::wclock();
      long *correct_sa = new long[text_length];
      run_divsufsort(text, correct_sa, (long)text_length);
      utils::write_to_file(correct_sa, text_length, sa_filename);
      delete[] correct_sa;
      fprintf(stderr, "Total time: %.2Lf\n", utils::wclock() - start);
    }
  }

  fprintf(stderr, "Running inmem sascan\n\n");
  std::uint8_t *computed_sa_temp =
    (std::uint8_t *)malloc(text_length * (sizeof(saidx_t) + 1));
  saidx_t *computed_sa = (saidx_t *)computed_sa_temp;
  start = utils::wclock();
  inmem_psascan<saidx_t>(text, text_length, computed_sa_temp,
      max_threads, false, false, NULL, max_blocks);
  long double total_time = utils::wclock() - start;

  if (compare) {
    fprintf(stderr, "\nComparing:\n");
    stream_reader<std::uint64_t> *sa_reader =
      new stream_reader<std::uint64_t>(sa_filename);
    bool eq = true;
    for (std::uint64_t i = 0; i < text_length; ++i) {
      std::uint64_t next_correct_sa = sa_reader->read();
      if (next_correct_sa != (std::uint64_t)computed_sa[i]) {
        eq = false;
        break;
      }
    }
    fprintf(stderr, "\nResult: %s\n",
        eq ? "\033[22;32mPASSED\033[0m" : "\033[22;31mFAILED\033[0m");
    if (!eq) std::exit(EXIT_FAILURE);
    delete sa_reader;
  }

  free(computed_sa_temp);
  return total_time;
}

template<typename saidx_t>
void test_file(
    const char *text_filename,
    std::uint64_t max_threads,
    std::uint64_t max_blocks,
    std::uint64_t runs) {

  fprintf(stderr, "Input filename: %s\n", text_filename);
  fprintf(stderr, "Reading text: ");
  std::uint64_t text_length = utils::file_size(text_filename);
  std::uint8_t *text = new std::uint8_t[text_length];
  utils::read_from_file(text, text_length, text_filename);
  fprintf(stderr, "DONE\n");

  std::vector<long double> times;
  for (std::uint64_t i = 0; i < runs; ++i)
    times.push_back(test<saidx_t>(text, text_length,
          max_threads, max_blocks, text_filename, (i == 0)));
  std::sort(times.begin(), times.end());

  long double mode_time = times[runs / 2];
  fprintf(stderr, "SUMMARY (mode of %lu): filename = %s, "
      "sizeof(saidx_t) = %lu, max_threads = %lu, "
      "time = %.2Lf (%.4Lfs/MiB), speed = %.2LfMiB/s\n",
      runs,
      text_filename,
      (std::uint64_t)sizeof(saidx_t),
      max_threads,
      mode_time,
      mode_time / ((long double)text_length / (1 << 20)),
      ((long double)text_length / (1 << 20)) / mode_time);

  delete[] text;
}

int main(int argc, char **argv) {
  if (argc == 1) {
    fprintf(stderr, "Usage: %s <max-threads> <max-blocks> "
        "<file1> <file2> ...\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }

  fprintf(stderr, "Command line:");
  for (int i = 0; i < argc; ++i)
    fprintf(stderr, " %s", argv[i]);
  fprintf(stderr, "\n");

  std::uint64_t max_threads = std::atol(argv[1]);
  std::uint64_t max_blocks = std::atol(argv[2]);

  for (int i = 3; i < argc; ++i) {
    test_file<uint40>(argv[i], max_threads, max_blocks, 3);
    test_file<int>(argv[i], max_threads, max_blocks, 3);
  }
}
