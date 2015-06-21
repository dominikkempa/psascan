#include <cstdio>
#include <algorithm>
#include <ctime>

#include <unistd.h>

#include "utils.h"
#include "async_bit_stream_writer.h"
#include "io_streamer.h"


void test(unsigned char *tab, long length) {
  // Choose the filename as random.
  std::string filename = "testfile." + utils::random_string_hash();

  // Create the async stream writer.
  long bufsize = utils::random_long(1, 100);  // in bytes
  typedef async_bit_stream_writer writer_type;
  writer_type *writer = new writer_type(filename, bufsize);

  // Write tab into filename using writer.
  for (long i = 0; i < length; ++i)
    writer->write(tab[i]);

  // Delete the writer.
  delete writer;

  // Now check if the writing was correct. Read the content of the file.
  int *tab_written = new int[length];
  bit_stream_reader *reader = new bit_stream_reader(filename);
  for (long i = 0; i < length; ++i)
    tab_written[i] = reader->read();
  delete reader;

  if (!std::equal(tab, tab + length, tab_written)) {
    fprintf(stderr, "Error!\n");
    if (length < 100) {
      fprintf(stderr, "\tlength = %ld\n", length);
      fprintf(stderr, "\tOriginal tab: ");
      for (long i = 0; i < length; ++i)
        fprintf(stderr, "%d ", tab[i]);
      fprintf(stderr, "\n");

      fprintf(stderr, "Read tab: ");
      for (long i = 0; i < length; ++i)
        fprintf(stderr, "%d ", tab_written[i]);
      fprintf(stderr, "\n");
    }

    if (utils::file_exists(filename))
      utils::file_delete(filename);
    std::exit(EXIT_FAILURE);
  }

  if (utils::file_exists(filename))
    utils::file_delete(filename);
  delete[] tab_written;
}


void test_random(int testcases, long max_length) {
  fprintf(stderr,"TEST, testcases = %d, max_n = %ld\n", testcases, max_length);

  unsigned char *tab = new unsigned char[max_length];

  for (int tc = 0; tc < testcases; ++tc) {
    // Print progress information.
    if (tc % 10 == 0)
      fprintf(stderr,"%d (%.2Lf%%)\r", tc, (tc * 100.L) / testcases);

    // Generate input.
    long length = utils::random_long(1L, max_length);

    for (long i = 0; i < length; ++i)
      tab[i] = utils::random_int(0, 1);

    // Run the test on generated input.
    test(tab, length);
  }

  // Clean up.
  delete[] tab;
}

int main() {
  std::srand(std::time(0) + getpid());

  test_random(1000000,  10);
  test_random(1000000,  100);
  test_random(100000,   1000);
  test_random(100000,   10000);
  test_random(1000,    100000);
  test_random(100,     1000000);

  fprintf(stderr,"All tests passed.\n");
}

