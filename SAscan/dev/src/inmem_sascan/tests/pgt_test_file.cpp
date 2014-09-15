#include <cstdio>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <thread>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "../../bitvector.h"
#include "../../multifile_bitvector.h"
#include "utils.h"
#include "io_streamer.h"
#include "../inmem_sascan.h"
#include "divsufsort_template.h"

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

void compute_gt_begin_for_text(unsigned char *supertext, long supertext_length, long text_beg, long text_end, bitvector *text_gt_begin) {
  for (long i = text_beg + 1; i <= text_end; ++i) {
    long lcp = 0;
    while (i + lcp < supertext_length && supertext[i + lcp] == supertext[text_beg + lcp]) ++lcp;
    
    if (i + lcp < supertext_length && supertext[i + lcp] > supertext[text_beg + lcp])
      text_gt_begin->set(i - text_beg - 1);
  }
}

void compute_gt_begin_reversed(unsigned char *text, long text_length, bitvector *gt_begin_reversed) {
  long i = 1, el = 0;
  while (i < text_length) {
    while (i + el < text_length && text[i + el] == text[el]) ++el;
    if (i + el < text_length && text[i + el] > text[el])
      gt_begin_reversed->set(text_length - i);

    el = 0;
    ++i;
  }
}

template<typename saidx_t, unsigned pagesize_log>
void test(std::string supertext_filename, long text_length, long max_threads) {
  long double start;

  fprintf(stderr, "Input filename: %s\n", supertext_filename.c_str());
  fprintf(stderr, "Reading text: ");
  long supertext_length;
  unsigned char *supertext;
  utils::read_objects_from_file(supertext, supertext_length, supertext_filename);
  fprintf(stderr, "DONE\n");


  std::string sa_filename = supertext_filename + ".sa" + utils::intToStr(sizeof(long));
  if (!utils::file_exists(sa_filename)) {
    fprintf(stderr, "Running divsufsort\n");
    start = utils::wclock();
    long *correct_sa = new long[supertext_length];
    run_divsufsort(supertext, correct_sa, supertext_length);
    utils::write_objects_to_file(correct_sa, supertext_length, sa_filename);
    delete[] correct_sa;
    fprintf(stderr, "Total time: %.2Lf\n", utils::wclock() - start);
  }


  text_length = std::min(text_length, supertext_length);
  long text_beg = utils::random_long(0L, supertext_length - text_length);
  long text_end = text_beg + text_length;


  bitvector *text_gt_begin_correct = new bitvector(text_length, max_threads);
  compute_gt_begin_for_text(supertext, supertext_length, text_beg, text_end, text_gt_begin_correct);


  // Compute tail_gt_begin_reversed.
  unsigned char *tail = supertext + text_end;
  long tail_length = supertext_length - text_end;
  bitvector *tail_gt_begin_reversed_bv = new bitvector(tail_length, max_threads);
  compute_gt_begin_reversed(tail, tail_length, tail_gt_begin_reversed_bv);

  // Store tail_gt_begin_reversed on disk as a multifile bitvector.
  multifile *tail_gt_begin_reversed_multifile = new multifile();
  long ptr = 0;
  while (ptr < tail_length) {
    long left = tail_length - ptr;
    long chunk = utils::random_long(1L, left);
   
    // Store bits [ptr..ptr+chunk) from tail_gt_begin_reversed_bv into one file.
    std::string chunk_filename = "gt_begin_reversed_bv" + utils::random_string_hash();
    bit_stream_writer *writer = new bit_stream_writer(chunk_filename);
    for (long j = ptr; j < ptr + chunk; ++j)
      writer->write(tail_gt_begin_reversed_bv->get(j));
    delete writer;

    // Add this file to tail_gt_begin_reversed_multifile.
    tail_gt_begin_reversed_multifile->add_file(ptr, ptr + chunk, chunk_filename);
    
    ptr += chunk;
  }
  delete tail_gt_begin_reversed_bv;


  unsigned char *text = (unsigned char *)malloc(text_length);
  std::copy(supertext + text_beg, supertext + text_end, text);
  delete[] supertext;



  // Run the tested algorithm.
  fprintf(stderr, "Running inmem sascan\n\n");
  unsigned char *bwtsa = (unsigned char *)malloc(text_length * (1 + sizeof(saidx_t)));
  bitvector *text_gt_begin_computed = new bitvector(text_length, max_threads);
  inmem_sascan<saidx_t, pagesize_log>(text, text_length, bwtsa, max_threads, false, // try also true here, to see if it works for both.
      true, text_gt_begin_computed, -1, text_beg, text_end, supertext_length, supertext_filename,
      tail_gt_begin_reversed_multifile);


  // Compute correct answer.

  ptr = 0;
  fprintf(stderr, "\nComparing:\n");
  bool eq = true;
  long compared = 0;
  for (long i = 0, dbg = 0; i < text_length; ++i) {
    ++dbg;
    ++compared;
    if (dbg == 10000000) {
      fprintf(stderr, "progress: %.3Lf%%\r", (100.L * i) / text_length);
      dbg = 0;
    }
    
    bool next_text_gt_begin_correct = text_gt_begin_correct->get(i);
    bool next_text_gt_begin_computed = text_gt_begin_computed->get(i);
    if (next_text_gt_begin_correct != next_text_gt_begin_computed) { eq = false; break; }
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

  long min_text_length = atol(argv[2]) << 20;
  test<uint40/*int*/, 12>(argv[1], min_text_length, 24);
}

