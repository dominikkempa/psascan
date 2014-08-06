#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <vector>
#include <string>

#include "utils.h"
#include "uint40.h"
#include "distributed_file.h"

template<typename data_type>
void test(data_type *seq, long length, long file_maxsize) {

  distributed_file<data_type> *f =
    new distributed_file<data_type>("temp", file_maxsize);

  long write_buf_size = utils::random_int((long)sizeof(data_type),
      (long)sizeof(data_type) + 100); // in bytes

  // debug //
  // fprintf(stderr, "\nWrite buffer size = %ld\n", write_buf_size);
  // fprintf(stderr, "Anticipated number of files = %ld\n",
  //     ((long)sizeof(data_type) * length + file_maxsize - 1) / file_maxsize);
  ///////////

  f->initialize_writing(write_buf_size);

  long prob = utils::random_int(0, 19);
  long left = length, beg = 0;

  // debug //
  // fprintf(stderr, "A total of %ld elems (%ld bytes) will now be written to"
  //     " distr-file (with %ld bytes limit on a single file)\n",
  //   length, (long)sizeof(data_type) * length, file_maxsize);
  ///////////

  while (left > 0) {
    // Flip a coin to decide whether we write a single element
    // or a block.
    if (utils::random_int(0, 19) < prob) {
      long block_size = utils::random_int(1, left);

      // debug //
      // fprintf(stderr, "  writing a block of %ld elems (%ld bytes)\n",
      //     block_size, (long)sizeof(data_type) * block_size);
      ///////////

      f->write(seq + beg, seq + beg + block_size);
      beg += block_size;
      left -= block_size;
    } else {
      // debug //
      // fprintf(stderr, "  writing a single elem (%ld bytes)\n",
      //     (long)sizeof(data_type));
      ///////////

      f->write(seq[beg]);
      ++beg;
      --left;
    }
  }
  // debug //
  // fprintf(stderr, "Writing done\n");
  ///////////

  f->finish_writing();

  long read_buf_size = utils::random_int((long)sizeof(data_type),
      (long)sizeof(data_type) + 100);
  f->initialize_reading(read_buf_size);

  // debug //
  // fprintf(stderr, "Now %ld elems (%ld bytes) will be read\n", length,
  //     (long)sizeof(long) * length);
  // fprintf(stderr, "Reading buffer size = %ld\n", read_buf_size);
  ///////////

  data_type *result = new data_type[length];
  for (long i = 0; i < length; ++i) {
    // debug //
    // fprintf(stderr, "  reading single elem\n");
    ///////////

    result[i] = f->read();
  }

  // debug //
  // fprintf(stderr, "Reading done\n");
  ///////////

  f->finish_reading();
  delete f;

  if (!std::equal(seq, seq + length, result)) {
    fprintf(stderr, "\n\033[22;31mFAILED\033[0m\n");
    std::exit(EXIT_FAILURE);
  }

  delete[] result;
}

void test_long(int testcases, long max_length) {
  fprintf(stderr, "TEST(long), testcases = %d, max_n = %ld\r", testcases, max_length);
  long *seq = new long[max_length];

  static const long MAX_LONG = 1000000000000000000L;
  for (int tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {
    if (dbg == 100) {
      fprintf(stderr, "TEST(long), testcases = %d, max_n = %ld: "
          "%d (%.0Lf%%)\r", testcases, max_length, tc, (100.L * tc) / testcases);
      dbg = 0;
    }

    long length = utils::random_long(1L, max_length);
    for (long i = 0; i < length; ++i)
      seq[i] = utils::random_long(-MAX_LONG, MAX_LONG);

    // Ensure that the file size is random, but not more than 100 files are created.
    int max_files = utils::random_int(1, 100);
    long seq_bytes = (long)sizeof(long) * length;
    long min_file_size = std::max((long)sizeof(long), (seq_bytes + max_files - 1) / max_files);
    long max_file_size = utils::random_long(min_file_size, seq_bytes);

    test(seq, length, max_file_size);
  }

  delete[] seq;

  fprintf(stderr, "TEST(long), testcases = %d, max_n = %ld: "
      "\033[22;32mPASSED\033[0m%10s\n", testcases, max_length, "");
}

void test_char(int testcases, long max_length) {
  fprintf(stderr, "TEST(char), testcases = %d, max_n = %ld\r", testcases, max_length);
  char *seq = new char[max_length];

  static const int MAX_CHAR = 127;
  for (int tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {
    if (dbg == 100) {
      fprintf(stderr, "TEST(char), testcases = %d, max_n = %ld: "
          "%d (%.0Lf%%)\r", testcases, max_length, tc, (100.L * tc) / testcases);
      dbg = 0;
    }

    long length = utils::random_long(1L, max_length);
    for (long i = 0; i < length; ++i)
      seq[i] = (char)utils::random_int(-MAX_CHAR, MAX_CHAR);

    // Ensure that the file size is random, but not more than 100 files are created.
    int max_files = utils::random_int(1, 100);
    long seq_bytes = (long)sizeof(char) * length;
    long min_file_size = std::max((long)sizeof(char), (seq_bytes + max_files - 1) / max_files);
    long max_file_size = utils::random_long(min_file_size, seq_bytes);

    test(seq, length, max_file_size);
  }

  delete[] seq;

  fprintf(stderr, "TEST(char), testcases = %d, max_n = %ld: "
      "\033[22;32mPASSED\033[0m%10s\n", testcases, max_length, "");
}

void test_short(int testcases, long max_length) {
  fprintf(stderr, "TEST(short), testcases = %d, max_n = %ld\r", testcases, max_length);
  short *seq = new short[max_length];

  static const int MAX_SHORT = 32000;
  for (int tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {
    if (dbg == 100) {
      fprintf(stderr, "TEST(short), testcases = %d, max_n = %ld: "
          "%d (%.0Lf%%)\r", testcases, max_length, tc, (100.L * tc) / testcases);
      dbg = 0;
    }

    long length = utils::random_long(1L, max_length);
    for (long i = 0; i < length; ++i)
      seq[i] = utils::random_int(-MAX_SHORT, MAX_SHORT);

    // Ensure that the file size is random, but not more than 100 files are created.
    int max_files = utils::random_int(1, 100);
    long seq_bytes = (long)sizeof(short) * length;
    long min_file_size = std::max((long)sizeof(short), (seq_bytes + max_files - 1) / max_files);
    long max_file_size = utils::random_long(min_file_size, seq_bytes);

    test(seq, length, max_file_size);
  }

  delete[] seq;

  fprintf(stderr, "TEST(short), testcases = %d, max_n = %ld: "
      "\033[22;32mPASSED\033[0m%10s\n", testcases, max_length, "");
}

void test_int(int testcases, long max_length) {
  fprintf(stderr, "TEST(int), testcases = %d, max_n = %ld\r", testcases, max_length);
  int *seq = new int[max_length];

  static const int MAX_INT = 1000000000;
  for (int tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {
    if (dbg == 100) {
      fprintf(stderr, "TEST(int), testcases = %d, max_n = %ld: "
          "%d (%.0Lf%%)\r", testcases, max_length, tc, (100.L * tc) / testcases);
      dbg = 0;
    }

    long length = utils::random_long(1L, max_length);
    for (long i = 0; i < length; ++i)
      seq[i] = utils::random_int(-MAX_INT, MAX_INT);

    // Ensure that the file size is random, but not more than 100 files are created.
    int max_files = utils::random_int(1, 100);
    long seq_bytes = (long)sizeof(int) * length;
    long min_file_size = std::max((long)sizeof(int), (seq_bytes + max_files - 1) / max_files);
    long max_file_size = utils::random_long(min_file_size, seq_bytes);

    test(seq, length, max_file_size);
  }

  delete[] seq;

  fprintf(stderr, "TEST(int), testcases = %d, max_n = %ld: "
      "\033[22;32mPASSED\033[0m%10s\n", testcases, max_length, "");
}


void test_double(int testcases, long max_length) {
  fprintf(stderr, "TEST(double), testcases = %d, max_n = %ld\r", testcases, max_length);
  double *seq = new double[max_length];

  static const double max_value = 1000000.0;
  for (int tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {
    if (dbg == 100) {
      fprintf(stderr, "TEST(double), testcases = %d, max_n = %ld: "
          "%d (%.0Lf%%)\r", testcases, max_length, tc, (100.L * tc) / testcases);
      dbg = 0;
    }

    long length = utils::random_long(1L, max_length);
    for (long i = 0; i < length; ++i)
      seq[i] = utils::random_double(-max_value, max_value);

    // Ensure that the file size is random, but not more than 100 files are created.
    int max_files = utils::random_int(1, 100);
    long seq_bytes = (long)sizeof(double) * length;
    long min_file_size = std::max((long)sizeof(double), (seq_bytes + max_files - 1) / max_files);
    long max_file_size = utils::random_long(min_file_size, seq_bytes);

    test(seq, length, max_file_size);
  }

  delete[] seq;

  fprintf(stderr, "TEST(double), testcases = %d, max_n = %ld: "
      "\033[22;32mPASSED\033[0m%10s\n", testcases, max_length, "");
}

void test_long_double(int testcases, long max_length) {
  fprintf(stderr, "TEST(long double), testcases = %d, max_n = %ld\r", testcases, max_length);
  long double *seq = new long double[max_length];

  static const double max_value = 1000000.0;
  for (int tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {
    if (dbg == 100) {
      fprintf(stderr, "TEST(long double), testcases = %d, max_n = %ld: "
          "%d (%.0Lf%%)\r", testcases, max_length, tc, (100.L * tc) / testcases);
      dbg = 0;
    }

    long length = utils::random_long(1L, max_length);
    for (long i = 0; i < length; ++i)
      seq[i] = (long double)utils::random_double(-max_value, max_value);

    // Ensure that the file size is random, but not more than 100 files are created.
    int max_files = utils::random_int(1, 100);
    long seq_bytes = (long)sizeof(long double) * length;
    long min_file_size = std::max((long)sizeof(long double), (seq_bytes + max_files - 1) / max_files);
    long max_file_size = utils::random_long(min_file_size, seq_bytes);

    test(seq, length, max_file_size);
  }

  delete[] seq;

  fprintf(stderr, "TEST(long double), testcases = %d, max_n = %ld: "
      "\033[22;32mPASSED\033[0m%10s\n", testcases, max_length, "");
}

void test_uint40(int testcases, long max_length) {
  fprintf(stderr, "TEST(uint40), testcases = %d, max_n = %ld\r", testcases, max_length);
  uint40 *seq = new uint40[max_length];

  static const long MAX_LONG = 100000000000000000L;
  for (int tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {
    if (dbg == 100) {
      fprintf(stderr, "TEST(uint40), testcases = %d, max_n = %ld: "
          "%d (%.0Lf%%)\r", testcases, max_length, tc, (100.L * tc) / testcases);
      dbg = 0;
    }

    long length = utils::random_long(1L, max_length);
    for (long i = 0; i < length; ++i)
      seq[i] = uint40((unsigned long)utils::random_long(0L, MAX_LONG));

    // Ensure that the file size is random, but not more than 100 files are created.
    int max_files = utils::random_int(1, 100);
    long seq_bytes = (long)sizeof(uint40) * length;
    long min_file_size = std::max((long)sizeof(uint40), (seq_bytes + max_files - 1) / max_files);
    long max_file_size = utils::random_long(min_file_size, seq_bytes);

    test(seq, length, max_file_size);
  }

  delete[] seq;

  fprintf(stderr, "TEST(uint40), testcases = %d, max_n = %ld: "
      "\033[22;32mPASSED\033[0m%10s\n", testcases, max_length, "");
}

int main(int, char **) {
  std::srand(std::time(0) + getpid());

  // Run tests.
  fprintf(stderr, "Testing distributed file.\n");

  test_char(1000, 50);
  test_char(1000, 500);
  test_char(1000, 5000);
  test_char(1000, 50000);
  test_char(1000, 500000);

  test_short(1000, 50);
  test_short(1000, 500);
  test_short(1000, 5000);
  test_short(1000, 50000);
  test_short(1000, 500000);

  test_int(1000, 50);
  test_int(1000, 500);
  test_int(1000, 5000);
  test_int(1000, 50000);
  test_int(1000, 500000);
  
  test_long(1000, 50);
  test_long(1000, 500);
  test_long(1000, 5000);
  test_long(1000, 50000);
  test_long(1000, 500000);

  test_uint40(1000, 50);
  test_uint40(1000, 500);
  test_uint40(1000, 5000);
  test_uint40(1000, 50000);
  test_uint40(1000, 500000);
  
  test_double(1000, 50);
  test_double(1000, 500);
  test_double(1000, 5000);
  test_double(1000, 50000);
  test_double(1000, 500000);
  
  test_long_double(1000, 50);
  test_long_double(1000, 500);
  test_long_double(1000, 5000);
  test_long_double(1000, 50000);
  test_long_double(1000, 500000);
}
