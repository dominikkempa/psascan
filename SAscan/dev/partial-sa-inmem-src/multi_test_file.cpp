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
long double test(unsigned char *text, long text_length, long max_threads,
    long max_blocks, std::string filename, bool compare = true) {
  long double start;
  std::string sa_filename = filename + ".sa" + utils::intToStr(sizeof(long));
  if (compare) {
    if (!utils::file_exists(sa_filename)) {
      fprintf(stderr, "Running divsufsort\n");
      start = utils::wclock();
      long *correct_sa = new long[text_length];
      run_divsufsort(text, correct_sa, text_length);
      utils::write_objects_to_file(correct_sa, text_length, sa_filename);
      delete[] correct_sa;
      fprintf(stderr, "Total time: %.2Lf\n", utils::wclock() - start);
    }
  }

  fprintf(stderr, "Running inmem sascan\n\n");
  unsigned char *computed_sa_temp = (unsigned char *)malloc(text_length * (sizeof(saidx_t) + 1));
  saidx_t *computed_sa = (saidx_t *)computed_sa_temp;
  start = utils::wclock();
  inmem_sascan<saidx_t>(text, text_length, computed_sa_temp, max_threads, false, false, NULL, max_blocks);
  long double total_time = utils::wclock() - start;

  if (compare) {
    fprintf(stderr, "\nComparing:\n");
    stream_reader<long> *sa_reader = new stream_reader<long>(sa_filename);
    bool eq = true;
    for (long i = 0; i < text_length; ++i) {
      long next_correct_sa = sa_reader->read();
      if (next_correct_sa != (long)computed_sa[i]) {
        eq = false;
        break;
      }
    }
    fprintf(stderr, "\nResult: %s\n", eq ? "OK" : "FAIL");
    free(computed_sa_temp);
    delete sa_reader;
  }

  return total_time;
}

template<typename saidx_t>
void test_file(const char *filename, long max_threads, long max_blocks, long runs) {
  fprintf(stderr, "Input filename: %s\n", filename);
  fprintf(stderr, "Reading text: ");
  long text_length;
  unsigned char *text;
  utils::read_objects_from_file(text, text_length, filename);
  fprintf(stderr, "DONE\n");

  std::vector<long double> times;
  for (long i = 0; i < runs; ++i)
    times.push_back(test<saidx_t>(text, text_length, max_threads, max_blocks, filename, (i == 0)));
  std::sort(times.begin(), times.end());

  long double mode_time = times[runs / 2];
  fprintf(stderr, "SUMMARY (mode of %ld): filename = %s, sizeof(saidx_t) = %ld, time = %.2Lf (%.4Lfs/MiB), speed = %.2LfMiB/s\n",
      runs,
      filename,
      (long)sizeof(saidx_t),
      mode_time,
      mode_time / ((long double)text_length / (1 << 20)),
      ((long double)text_length / (1 << 20)) / mode_time);

  delete[] text;
}

int main(int argc, char **argv) {
  if (argc == 1) {
    fprintf(stderr, "Usage: %s <max-threads> <max-blocks> <file1> <file2> ...\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }

  fprintf(stderr, "Command line:");
  for (long i = 0; i < argc; ++i)
    fprintf(stderr, " %s", argv[i]);
  fprintf(stderr, "\n");

  long max_threads = std::atol(argv[1]);
  long max_blocks = std::atol(argv[2]);

  for (long i = 3; i < argc; ++i) {
    test_file<uint40>(argv[i], max_threads, max_blocks, 3);
    test_file<int>(argv[i], max_threads, max_blocks, 3);
  }
}
