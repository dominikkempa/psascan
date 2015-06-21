#include <cstdio>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <thread>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "../../src/psascan_src/bitvector.h"
#include "../../src/psascan_src/multifile.h"
#include "inmem_psascan.h"
#include "divsufsort.h"
#include "divsufsort64.h"
#include "utils.h"
#include "io_streamer.h"


void compute_gt_begin_for_text(unsigned char *supertext, long supertext_length, long text_beg, long text_end, psascan_private::bitvector *text_gt_begin) {
  for (long i = text_beg + 1; i <= text_end; ++i) {
    long lcp = 0;
    while (i + lcp < supertext_length && supertext[i + lcp] == supertext[text_beg + lcp]) ++lcp;
    
    if (i + lcp < supertext_length && supertext[i + lcp] > supertext[text_beg + lcp])
      text_gt_begin->set(i - text_beg - 1);
  }
}

void compute_gt_begin_reversed(unsigned char *text, long text_length, psascan_private::bitvector *gt_begin_reversed) {
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

  long text_length = text_end - text_beg;
  psascan_private::bitvector *text_gt_begin_correct = new psascan_private::bitvector(text_length);
  compute_gt_begin_for_text(supertext, supertext_length, text_beg, text_end, text_gt_begin_correct);


  // Compute tail_gt_begin_reversed.
  unsigned char *tail = supertext + text_end;
  long tail_length = supertext_length - text_end;
  psascan_private::bitvector tail_gt_begin_reversed_bv(tail_length);
  compute_gt_begin_reversed(tail, tail_length, &tail_gt_begin_reversed_bv);

  // Store tail_gt_begin_reversed on disk as a multifile bitvector.
  psascan_private::multifile *tail_gt_begin_reversed_multifile = new psascan_private::multifile();
  long ptr = 0;
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
  psascan_private::bitvector *text_gt_begin_computed = new psascan_private::bitvector(text_length);
  long max_blocks = -1;
  if (utils::random_long(0, 1)) max_blocks = utils::random_long(1, 50L);
  bool compute_bwt = (bool)utils::random_long(0, 1);
  inmem_psascan<saidx_t, pagesize_log>(text, text_length, bwtsa, max_threads, compute_bwt,
      true, text_gt_begin_computed, max_blocks, text_beg, text_end, supertext_length, supertext_filename,
      tail_gt_begin_reversed_multifile);



  // Compare answers.
  bool eq = true;
  for (long i = 0; i < text_length; ++i)
    if (text_gt_begin_correct->get(i) != text_gt_begin_computed->get(text_length - 1 - i)) { eq = false; break; }

  if (!eq) {
    fprintf(stdout, "Error:\n");
    fprintf(stdout, "\tsupertext_length = %ld\n", supertext_length);
    fprintf(stdout, "\tmax threads = %ld\n", max_threads);
    fprintf(stderr, "\tmax blocks = %ld\n", max_blocks);
    fprintf(stdout, "\tsupertext = ");
    for (long j = 0; j < supertext_length; ++j)
      fprintf(stdout, "%c", supertext[j]);
    fprintf(stdout, "\n");
    fprintf(stdout, "\ttext_beg = %ld\n", text_beg);
    fprintf(stdout, "\ttext_end = %ld\n", text_end);
    fprintf(stdout, "\ttail_gt_begin_reversed_bv = ");
    for (long j = 0; j < tail_length; ++j)
      fprintf(stdout, "%ld", (long)tail_gt_begin_reversed_bv.get(j));
    fprintf(stdout, "\n");
    fprintf(stdout, "\tcorrect text_gt_begin: ");
    for (long i = 0; i < text_length; ++i)
      fprintf(stdout, "%ld", (long)text_gt_begin_correct->get(i));
    fprintf(stdout, "\n");
        fprintf(stdout, "\tcomputed text_gt_begin: ");
    for (long i = 0; i < text_length; ++i)
      fprintf(stdout, "%ld", (long)text_gt_begin_computed->get(i));
    fprintf(stdout, "\n");
    std::fflush(stdout);

    std::exit(EXIT_FAILURE);
  }


  free(bwtsa);
  delete tail_gt_begin_reversed_multifile; // also deletes files
  delete text_gt_begin_correct;
  delete text_gt_begin_computed;  
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

