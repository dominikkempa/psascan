#include <cstdio>
#include <algorithm>
#include <ctime>

#include <unistd.h>

#include "utils.h"
#include "async_backward_stream_reader.h"


void test(int *tab, long length) {
  // Choose the filename as random.
  std::string filename = "testfile." + utils::random_string_hash();

  // Write tab into file.
  utils::write_objects_to_file(tab, length, filename);

  // Create the async stream reader.
  long bufsize = utils::random_long(1, 100);  // in bytes
  async_backward_stream_reader<int> *reader = new async_backward_stream_reader<int>(filename, bufsize);

  // Read tab into filename using writer.
  int *tab_read = new int[length];
  for (long i = 0; i < length; ++i)
    tab_read[i] = reader->read();
  std::reverse(tab_read, tab_read + length);

  // Delete the reader.
  delete reader;

  if (!std::equal(tab, tab + length, tab_read)) {
    fprintf(stderr, "Error!\n");
    if (length < 100) {
      fprintf(stderr, "\tlength = %ld\n", length);
      fprintf(stderr, "\tOriginal tab: ");
      for (long i = 0; i < length; ++i)
        fprintf(stderr, "%d ", tab[i]);
      fprintf(stderr, "\n");

      fprintf(stderr, "Read tab: ");
      for (long i = 0; i < length; ++i)
        fprintf(stderr, "%d ", tab_read[i]);
      fprintf(stderr, "\n");
    }

    if (utils::file_exists(filename))
      utils::file_delete(filename);
    std::exit(EXIT_FAILURE);
  }

  if (utils::file_exists(filename))
    utils::file_delete(filename);
  delete[] tab_read;
}


void test_random(int testcases, long max_length) {
  fprintf(stderr,"TEST, testcases = %d, max_n = %ld\n", testcases, max_length);

  int *tab = new int[max_length];

  for (int tc = 0; tc < testcases; ++tc) {
    // Print progress information.
    if (tc % 10 == 0)
      fprintf(stderr,"%d (%.2Lf%%)\r", tc, (tc * 100.L) / testcases);

    // Generate input.
    long length = utils::random_long(1L, max_length);

    for (long i = 0; i < length; ++i)
      tab[i] = utils::random_int(1, 1000000000);

    // Run the test on generated input.
    test(tab, length);
  }

  // Clean up.
  delete[] tab;
}

int main() {
  std::srand(std::time(0) + getpid());

  test_random(10000,  10);
  test_random(10000,  100);
  test_random(1000,   1000);
  test_random(1000,   10000);
  test_random(100,    100000);
  test_random(10,     1000000);

  fprintf(stderr,"All tests passed.\n");
}

