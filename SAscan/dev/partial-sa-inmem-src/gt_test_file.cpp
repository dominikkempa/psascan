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
#include "bitvector.h"
#include "io_streamer.h"

template<typename T>
void next(unsigned char *text, T length, T &s, T &p, T &r) {
  if (length == 1) { s = 0; p = 1; r = 0; return; }
  T i = length - 1;
  while (i < length) {
    unsigned char a = text[s + r], b = text[i];
    if (a > b) { p = i - s + 1; r = 0; }
    else if (a < b) { i -= r; s = i; p = 1; r = 0; }
    else { ++r; if (r == p) r = 0; } ++i;
  }
}

//==============================================================================
// Compute gt_begin for text.
//==============================================================================
long compute_gt_begin(unsigned char *text, long length, bitvector *gt_begin) {
  long whole_suffix_rank = length - 1;
  long i = 1, el = 0, s = 0, p = 0, r = 0;
  long i_max = 0, el_max = 0, s_max = 0, p_max = 0, r_max = 0;
  while (i < length) {
    while (i + el < length && el < length && text[i + el] == text[el])
      next(text, ++el, s, p, r);
 
    if (i + el < length && (el == length || text[i + el] > text[el])) {
      gt_begin->set(i - 1);
      --whole_suffix_rank;
    }

    long j = i_max;
    if (el > el_max) {
      std::swap(el, el_max);
      std::swap(s, s_max);
      std::swap(p, p_max);
      std::swap(r, r_max);
      i_max = i;
    }

    if (p && 3 * p <= el && !memcmp(text, text + p, s)) {
      for (long k = 1; k < std::min(length - i, p); ++k) { 
        if (gt_begin->get(j + k - 1)) {
          gt_begin->set(i + k - 1);
          --whole_suffix_rank;
        }
      }

      i += p;
      el -= p;
    } else {
      long h = (el / 3) + 1;
      for (long k = 1; k < std::min(length - i, h); ++k) { 
        if (gt_begin->get(j + k - 1)) {
          gt_begin->set(i + k - 1);
          --whole_suffix_rank;
        }
      }

      i += h;
      el = 0;
      s = 0;
      p = 0;
    }
  }

  return whole_suffix_rank;
}

template<typename saidx_t>
void read_sa(saidx_t* &sa, std::string filename) {
  long length;
  utils::read_objects_from_file(sa, length, filename);
}


template<typename saidx_t>
void test(unsigned char *text, long text_length, long max_threads,
    long max_blocks, std::string filename) {
  long double start;

  std::string gt_begin_filename = filename + ".gt_begin";

  if (!utils::file_exists(gt_begin_filename)) {
    fprintf(stderr, "Computing gt begin\n");
    start = utils::wclock();
    bitvector gt_begin(text_length, max_threads);
    bit_stream_writer *writer = new bit_stream_writer(gt_begin_filename);
    compute_gt_begin(text, text_length, &gt_begin);
    for (long i = 0; i < text_length; ++i)
      writer->write(gt_begin.get(i));
    delete writer;
    fprintf(stderr, "Total time: %.2Lf\n", utils::wclock() - start);
  }

  fprintf(stderr, "Running inmem sascan\n\n");
  unsigned char *computed_sa_temp = (unsigned char *)malloc(text_length * (sizeof(saidx_t) + 1));
  start = utils::wclock();
  bitvector *gt_begin = new bitvector(text_length, max_threads);
  inmem_sascan<saidx_t>(text, text_length, computed_sa_temp, max_threads, false, true, gt_begin, max_blocks);
  long double total_time = utils::wclock() - start;
  fprintf(stderr, "\nTotal time:\n");
  fprintf(stderr, "\tabsolute: %.2Lf\n", total_time);
  fprintf(stderr, "\trelative: %.4Lfs/MiB\n", total_time / ((long double)text_length / (1 << 20)));
  fprintf(stderr, "Speed: %.2LfMiB/s\n", ((long double)text_length / (1 << 20)) / total_time);

  fprintf(stderr, "\nComparing:\n");
  bit_stream_reader *gt_reader = new bit_stream_reader(gt_begin_filename);
  bool eq = true;
  long compared = 0;
  for (long i = 0, dbg = 0; i < text_length; ++i) {
    ++dbg;
    ++compared;
    if (dbg == 10000000) {
      fprintf(stderr, "progress: %.3Lf%%\r", (100.L * i) / text_length);
      dbg = 0;
    }

    bool next_correct_gt = gt_reader->read();
    bool next_computed_gt = gt_begin->get(i);
    if (next_correct_gt != next_computed_gt) {
      eq = false;
      break;
    }
  }
  fprintf(stderr, "Compared %ld values", compared);
  fprintf(stderr, "\nResult: %s\n", eq ? "OK" : "FAIL");

  free(computed_sa_temp);
  delete gt_reader;
  delete gt_begin;
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
