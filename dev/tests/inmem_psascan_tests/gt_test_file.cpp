#include <cstdio>
#include <string>
#include <cstring>
#include <ctime>
#include <unistd.h>

#include "../../src/psascan_src/bitvector.hpp"
#include "inmem_psascan.hpp"
#include "utils.hpp"
#include "divsufsort_template.hpp"
#include "uint40.hpp"
#include "io_streamer.hpp"


template<typename text_offset_type>
void next(
    std::uint8_t *text,
    text_offset_type length,
    text_offset_type &s,
    text_offset_type &p,
    text_offset_type &r) {

  if (length == 1) {
    s = 0;
    p = 1;
    r = 0;
    return;
  }

  text_offset_type i = length - 1;
  while (i < length) {
    std::uint8_t a = text[s + r], b = text[i];
    if (a > b) { p = i - s + 1; r = 0; }
    else if (a < b) { i -= r; s = i; p = 1; r = 0; }
    else { ++r; if (r == p) r = 0; } ++i;
  }
}

//==============================================================================
// Compute gt_begin for text.
//==============================================================================
std::uint64_t compute_gt_begin(
    std::uint8_t *text,
    std::uint64_t length,
    psascan_private::bitvector *gt_begin) {

  std::uint64_t whole_suffix_rank = length - 1;
  std::uint64_t i = 1, el = 0, s = 0, p = 0, r = 0;
  std::uint64_t i_max = 0, el_max = 0, s_max = 0, p_max = 0, r_max = 0;
  while (i < length) {
    while (i + el < length && el < length && text[i + el] == text[el])
      next(text, ++el, s, p, r);
 
    if (i + el < length && (el == length || text[i + el] > text[el])) {
      gt_begin->set(i - 1);
      --whole_suffix_rank;
    }

    std::uint64_t j = i_max;
    if (el > el_max) {
      std::swap(el, el_max);
      std::swap(s, s_max);
      std::swap(p, p_max);
      std::swap(r, r_max);
      i_max = i;
    }

    if (p && 3 * p <= el && !memcmp(text, text + p, s)) {
      for (std::uint64_t k = 1; k < std::min(length - i, p); ++k) { 
        if (gt_begin->get(j + k - 1)) {
          gt_begin->set(i + k - 1);
          --whole_suffix_rank;
        }
      }

      i += p;
      el -= p;
    } else {
      std::uint64_t h = (el / 3) + 1;
      for (std::uint64_t k = 1; k < std::min(length - i, h); ++k) { 
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
void test(
    std::uint8_t *text,
    std::uint64_t text_length,
    std::uint64_t max_threads,
    std::uint64_t max_blocks,
    std::string filename) {

  long double start;

  std::string gt_begin_filename = filename + ".gt_begin";

  if (!utils::file_exists(gt_begin_filename)) {
    fprintf(stderr, "Computing gt begin\n");
    start = utils::wclock();
    psascan_private::bitvector gt_begin(text_length);
    bit_stream_writer *writer = new bit_stream_writer(gt_begin_filename);
    compute_gt_begin(text, text_length, &gt_begin);
    for (std::uint64_t i = 0; i < text_length; ++i)
      writer->write(gt_begin.get(i));
    delete writer;
    fprintf(stderr, "Total time: %.2Lf\n", utils::wclock() - start);
  }

  fprintf(stderr, "Running inmem sascan\n\n");
  std::uint8_t *computed_sa_temp =
    (std::uint8_t *)malloc(text_length * (sizeof(saidx_t) + 1));
  psascan_private::bitvector *gt_begin =
    new psascan_private::bitvector(text_length);
  inmem_psascan<saidx_t>(text, text_length, computed_sa_temp,
      max_threads, false, true, gt_begin, max_blocks);

  fprintf(stderr, "\nComparing:\n");
  bit_stream_reader *gt_reader = new bit_stream_reader(gt_begin_filename);
  bool eq = true;
  std::uint64_t compared = 0;
  for (std::uint64_t i = 0, dbg = 0; i < text_length; ++i) {
    ++dbg;
    ++compared;
    if (dbg == 10000000) {
      fprintf(stderr, "progress: %.3Lf%%\r", (100.L * i) / text_length);
      dbg = 0;
    }

    bool next_correct_gt = gt_reader->read();
    bool next_computed_gt = gt_begin->get(text_length - 1 - i);
    if (next_correct_gt != next_computed_gt) {
      eq = false;
      break;
    }
  }
  fprintf(stderr, "Compared %lu values", compared);
  fprintf(stderr, "\nResult: %s\n", eq ? "OK" : "FAIL");

  free(computed_sa_temp);
  delete gt_reader;
  delete gt_begin;
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
