#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <vector>
#include <string>

#include "utils.h"
#include "distributed_file.h"

void test(unsigned char *seq, long length) {

  // Write data to disk (to several files).
  multifile *m = new multifile();
  long left = length, beg = 0;
  while (left > 0) {
    long towrite = utils::random_int(1, left);

    std::string filename = "tempfile." + utils::random_string_hash();
    utils::write_objects_to_file(seq + beg, towrite, filename);
    m->add_file(beg, beg + towrite, filename);
 
    beg += towrite;
    left -= towrite;
  }

  // Read the data back to memory using distributed file.
  long bufsize = utils::random_long(1L, 100L);
  distributed_file *f = new distributed_file(m, bufsize);

  unsigned char *result = new unsigned char[length];
  for (long i = 0; i < length; ++i)
    result[i] = f->read();

  delete f;
  delete m;  // this also delete all files from multifile

  if (!std::equal(seq, seq + length, result)) {
    fprintf(stderr, "\n\033[22;31mFAILED\033[0m\n");
    std::exit(EXIT_FAILURE);
  }

  delete[] result;
}

void test(int testcases, long max_length) {
  fprintf(stderr, "TEST, testcases = %d, max_n = %ld\r", testcases, max_length);
  unsigned char *seq = new unsigned char[max_length];

  for (int tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {
    if (dbg == 100) {
      fprintf(stderr, "TEST, testcases = %d, max_n = %ld: "
          "%d (%.0Lf%%)\r", testcases, max_length, tc, (100.L * tc) / testcases);
      dbg = 0;
    }

    long length = utils::random_long(1L, max_length);
    for (long i = 0; i < length; ++i)
      seq[i] = utils::random_int(0, 255);

    test(seq, length);
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
  test(100, 500000);
}
