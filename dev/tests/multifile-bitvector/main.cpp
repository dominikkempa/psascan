#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <vector>
#include <string>

#include "utils.h"
#include "multifile_bitvector.h"
#include "io_streamer.h"


void test(long length, bool *seq, multifile_bitvector &bv) {
  // Test access operation.
  for (long i = 0; i < 1000; ++i) {
    // Choose starting position.
    long pos = utils::random_long(0L, length - 1);
    long left = length - pos;

    while (pos < length) {
      bool extracted_bit = bv.access(pos);
      if (extracted_bit != seq[pos]) {
        fprintf(stderr, "\nError in access:\n");
        fprintf(stderr, "\tlength = %ld\n", length);
        fprintf(stderr, "\tseq: ");
        for (long j = 0; j < length; ++j) fprintf(stderr, "%ld", (long)seq[j]);
        fprintf(stderr, "\n");
        fprintf(stderr, "\tpos = %ld\n", pos);
        fprintf(stderr, "\tcorrect bit = %ld\n", (long)seq[pos]);
        fprintf(stderr, "\tcomputed bit = %ld\n", (long)extracted_bit);
        std::exit(EXIT_FAILURE);
      }

      long skip = utils::random_long(1L, left);
      pos += skip;
    }
  }

  // Test sequential reading.
  for (long i = 0; i < 1000; ++i) {
    // Choose starting position.
    long pos = utils::random_long(0L, length - 1);
    bv.initialize_sequential_reading(pos);

    for (long j = pos; j < length; ++j) {
      bool extracted_bit = bv.read();
      if (extracted_bit != seq[j]) {
        fprintf(stderr, "\nError in read:\n");
        fprintf(stderr, "\tseq: ");
        for (long jj = 0; jj < length; ++jj)
          fprintf(stderr, "%ld", (long)seq[jj]);
        fprintf(stderr, "\n");
        fprintf(stderr, "\tj = %ld\n", j);
        fprintf(stderr, "\tcorrect bit = %ld\n", (long)seq[j]);
        fprintf(stderr, "\tcomputed bit = %ld\n", (long)extracted_bit);
        std::exit(EXIT_FAILURE);
      }
    }
  }
}

void test(int testcases, long max_length) {
  fprintf(stderr, "TEST, testcases = %d, max_n = %ld\r", testcases, max_length);

  bool *seq = new bool[max_length];
  for (long tc = 0; tc < testcases; ++tc) {
    fprintf(stderr, "TEST, testcases = %d, max_n = %ld: "
        "%ld (%.0Lf%%)\r", testcases, max_length, tc, (100.L * tc) / testcases);

    long length = utils::random_long(1L, max_length);
    for (long i = 0; i < length; ++i)
      seq[i] = (bool)utils::random_long(0, 1);

    // splite seq into many parts, write to files and
    // create the multifile bitvector object.
    multifile_bitvector bv;
    std::vector<single_file_info> infos;

    long pos = 0;
    long left = length;
    long max_chunk = 1000;
    while (pos < length) {
      long chunk = utils::random_long(1L, std::min(max_chunk, left));
      std::string filename = "bv" + utils::random_string_hash();
      bit_stream_writer *writer = new bit_stream_writer(filename);
      for (long j = 0; j < chunk; ++j)
        writer->write(seq[pos + j]);
      delete writer;
 
      infos.push_back(single_file_info(pos, pos + chunk, filename));
      pos += chunk;
      left -= chunk;
    }

    std::random_shuffle(infos.begin(), infos.end());
    for (long i = 0; i < (long)infos.size(); ++i)
      bv.add_file(infos[i].m_beg, infos[i].m_end, infos[i].m_filename);

    test(length, seq, bv);
  }

  delete[] seq;

  fprintf(stderr, "TEST, testcases = %d, max_n = %ld: "
      "\033[22;32mPASSED\033[0m%10s\n", testcases, max_length, "");
}


int main(int, char **) {
   std::srand(std::time(0) + getpid());

  test(1000, 50);
  test(1000, 500);
  test(1000, 5000);
  test(1000, 50000);
  test(1000, 500000);

  fprintf(stderr, "All tests passed.\n");
}
