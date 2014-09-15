#include <cstdio>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <thread>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "divsufsort.h"
#include "divsufsort64.h"
#include "../../bitvector.h"
#include "../../multifile_bitvector.h"
#include "utils.h"
#include "io_streamer.h"
#include "../inmem_sascan.h"


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
void test(unsigned char *supertext, long supertext_length,
    long text_beg, long text_end, long max_threads) {

  // Sort supertext using divsufsort.
  long *supertext_sa = (long *)malloc(supertext_length * sizeof(long));
  divsufsort64(supertext, supertext_sa, supertext_length);

  // Separate the ordering of the suffixes of text (correct answer).
  long text_length = text_end - text_beg;
  long *correct_answer = (long *)malloc(text_length * sizeof(long));
  long ptr = 0;
  for (long i = 0; i < supertext_length; ++i)
    if (text_beg <= supertext_sa[i] && supertext_sa[i] < text_end)
      correct_answer[ptr++] = supertext_sa[i] - text_beg;

  // Compute tail_gt_begin_reversed.
  unsigned char *tail = supertext + text_end;
  long tail_length = supertext_length - text_end;
  bitvector tail_gt_begin_reversed_bv(tail_length, max_threads);
  compute_gt_begin_reversed(tail, tail_length, &tail_gt_begin_reversed_bv);

  // Store tail_gt_begin_reversed on disk as a multifile bitvector.
  multifile *tail_gt_begin_reversed_multifile = new multifile();
  ptr = 0;
  while (ptr < tail_length) {
    long left = tail_length - ptr;
    long chunk = utils::random_long(1L, left);
   
    // Store bits [ptr..ptr+chunk) from tail_gt_begin_reversed_bv into one file.
    std::string chunk_filename = "gt_begin_reversed_bv" + utils::random_string_hash();
    bit_stream_writer *writer = new bit_stream_writer(chunk_filename);
    for (long j = ptr; j < ptr + chunk; ++j)
      writer->write(tail_gt_begin_reversed_bv.get(j));
    delete writer;

    // Add this file to tail_gt_begin_reversed_multifile.
    tail_gt_begin_reversed_multifile->add_file(ptr, ptr + chunk, chunk_filename);
    
    ptr += chunk;
  }

  // Write supertext to file.
  std::string supertext_filename = "supertext.txt";
  stream_writer<unsigned char> *supertext_writer = new stream_writer<unsigned char>(supertext_filename);
  for (long i = 0; i < supertext_length; ++i)
    supertext_writer->write(supertext[i]);
  delete supertext_writer;






  // Run the tested algorithm.
  unsigned char *text = supertext + text_beg;
  unsigned char *bwtsa = (unsigned char *)malloc(text_length * (1 + sizeof(saidx_t)));
  saidx_t *computed_sa = (saidx_t *)bwtsa;
  unsigned char *computed_bwt = (unsigned char *)(computed_sa + text_length);
  long max_blocks = -1;
  if (utils::random_long(0, 1)) max_blocks = utils::random_long(1L, 50L);
  long computed_i0;
  inmem_sascan<saidx_t, pagesize_log>(text, text_length, bwtsa, max_threads, true,
      false, NULL, max_blocks, text_beg, text_end, supertext_length, supertext_filename,
      tail_gt_begin_reversed_multifile, &computed_i0);



  // Compare answers.
  bool eq = true;
  long correct_i0;
  for (long i = 0; i < text_length; ++i) {
    unsigned char correct_bwt = ((correct_answer[i] == 0) ? 0 : text[correct_answer[i] - 1]);
    if (computed_bwt[i] != correct_bwt) { eq = false; break; }
    if (correct_answer[i] == 0) correct_i0 = i;
  }
  if (correct_i0 != computed_i0) eq = false;

  if (!eq) {
    fprintf(stdout, "Error:\n");
    fprintf(stdout, "\tsupertext_length = %ld\n", supertext_length);
    fprintf(stdout, "\tmax threads = %ld\n", max_threads);
    fprintf(stdout, "\tsupertext = ");
    for (long j = 0; j < supertext_length; ++j)
      fprintf(stdout, "%c", supertext[j]);
    fprintf(stdout, "\n");
    fprintf(stdout, "\ttext_beg = %ld\n", text_beg);
    fprintf(stdout, "\ttext_end = %ld\n", text_end);
    fprintf(stderr, "\tcorrect i0 = %ld\n", correct_i0);
    fprintf(stderr, "\tcomputed i0 = %ld\n", computed_i0);
    fprintf(stdout, "\ttail_gt_begin_reversed_bv = ");
    for (long j = 0; j < tail_length; ++j)
      fprintf(stdout, "%ld", (long)tail_gt_begin_reversed_bv.get(j));
    fprintf(stdout, "\n");
    fprintf(stdout, "\tsupertext sa = ");
    for (long j = 0; j < supertext_length; ++j)
      fprintf(stdout, "%ld ", (long)supertext_sa[j]);
    fprintf(stdout, "\n");
    fprintf(stdout, "\tcorrect bwt: ");
    for (long i = 0; i < text_length; ++i)
      fprintf(stdout, "%c", ((correct_answer[i] == 0) ? 0 : text[correct_answer[i] - 1]));
    fprintf(stdout, "\n");
    fprintf(stdout, "\tcomputed bwt: ");
    for (long i = 0; i < text_length; ++i)
      fprintf(stdout, "%c", computed_bwt[i]);
    fprintf(stdout, "\n");
    std::fflush(stdout);
    std::fflush(stdout);

    std::exit(EXIT_FAILURE);
  }


  delete tail_gt_begin_reversed_multifile; // also deletes files
  free(bwtsa);
  free(correct_answer);
  free(supertext_sa);
}


template<typename saidx_t, unsigned pagesize_log>
void test_random(int testcases, long max_length, int max_sigma) {
  fprintf(stdout,"TEST, testcases = %d, max_n = %ld, max_sigma = %d, sizeof(saidx_t) = %ld, pagesize_log = %ld\n",
      testcases, max_length, max_sigma, (long)sizeof(saidx_t), (long)pagesize_log);
  unsigned char *supertext = new unsigned char[max_length + 1];

  for (int tc = 0; tc < testcases; ++tc) {
    // Print progress information.
    fprintf(stdout,"%d (%.2Lf%%)\r", tc, (tc * 100.L) / testcases);
    std::fflush(stdout);




    // Generate string.
    long supertext_length = utils::random_long(1, max_length);
    int sigma = utils::random_int(2, max_sigma);
    if (max_sigma <= 26) utils::fill_random_letters(supertext, supertext_length, sigma);
    else utils::fill_random_string(supertext, supertext_length, sigma);
    long max_threads = utils::random_long(1, 50);
    long text_beg = utils::random_long(0, supertext_length - 1);
    long text_end = utils::random_long(text_beg + 1, supertext_length);



    /*long supertext_length = 325;
    long max_threads = 16;
    supertext[0] = 0;
    strcpy((char *)supertext, "adbbddacdcdacacccddcbdccbaaaabcaadcadcccabccbaaddabadadbadbaadabdcbcadaaacdbcdbbdccdcbacabcaadbdbcbcbbcbdbdbbadacbdbcddcaccbbaacccaaddbdaaabadcdabacbdabbdccddbbbbbaaddbdacadadacdcdcdaacdcbcdcdbadbddccdbccbbcdabcdaddccbdabbdcbcdabcdadbdadbadccccbbaddddabcccbbdcdcdcdcdcddccbaadcbcbacbbadabadaabdcabbdaabccbdbdadaabbccacdbbdbcc");
    long text_beg = 125;
    long text_end = 269;*/


    // Run the test on generated string.
    test<saidx_t, pagesize_log>(supertext, supertext_length, text_beg, text_end, max_threads);
  }

  // Clean up.
  delete[] supertext;
}


int main() {
  std::srand(std::time(0) + getpid());

  // Redirect stdout to /dev/null
  int redir = open("/dev/null", O_WRONLY);
  dup2(redir, 2);
  close(redir);

  test_random<uint40, 2>(10000,   10,      5);
  test_random<uint40, 5>(10000,   10,      5);
  test_random<uint40, 8>(10000,   10,      5);
  test_random<uint40, 2>(10000,   10,    255);
  test_random<uint40, 5>(10000,   10,    255);
  test_random<uint40, 8>(10000,   10,    255);
  test_random<int,    2>(10000,   10,      5);
  test_random<int,    5>(10000,   10,      5);
  test_random<int,    8>(10000,   10,      5);
  test_random<int,    2>(10000,   10,    255);
  test_random<int,    5>(10000,   10,    255);
  test_random<int,    8>(10000,   10,    255);

  test_random<uint40, 2>(1000,   100,      5);
  test_random<uint40, 5>(1000,   100,      5);
  test_random<uint40, 8>(1000,   100,      5);
  test_random<uint40, 2>(1000,   100,    255);
  test_random<uint40, 5>(1000,   100,    255);
  test_random<uint40, 8>(1000,   100,    255);
  test_random<int,    2>(1000,   100,      5);
  test_random<int,    5>(1000,   100,      5);
  test_random<int,    8>(1000,   100,      5);
  test_random<int,    2>(1000,   100,    255);
  test_random<int,    5>(1000,   100,    255);
  test_random<int,    8>(1000,   100,    255);

  test_random<uint40, 2>(200,   1000,      5);
  test_random<uint40, 5>(200,   1000,      5);
  test_random<uint40, 8>(200,   1000,      5);
  test_random<uint40, 2>(200,   1000,    255);
  test_random<uint40, 5>(200,   1000,    255);
  test_random<uint40, 8>(200,   1000,    255);
  test_random<int,    2>(200,   1000,      5);
  test_random<int,    5>(200,   1000,      5);
  test_random<int,    8>(200,   1000,      5);
  test_random<int,    2>(200,   1000,    255);
  test_random<int,    5>(200,   1000,    255);
  test_random<int,    8>(200,   1000,    255);

  test_random<uint40, 2>(20,   1000000,      5);
  test_random<uint40, 5>(20,   1000000,      5);
  test_random<uint40, 8>(20,   1000000,      5);
  test_random<uint40, 2>(20,   1000000,    255);
  test_random<uint40, 5>(20,   1000000,    255);
  test_random<uint40, 8>(20,   1000000,    255);
  test_random<int,    2>(20,   1000000,      5);
  test_random<int,    5>(20,   1000000,      5);
  test_random<int,    8>(20,   1000000,      5);
  test_random<int,    2>(20,   1000000,    255);
  test_random<int,    5>(20,   1000000,    255);
  test_random<int,    8>(20,   1000000,    255);

  fprintf(stdout,"All tests passed.\n");
  std::fflush(stdout);
}

