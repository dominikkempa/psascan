#include <cstdio>
#include <algorithm>
#include <ctime>

#include <stdint.h>

#include <unistd.h>

#include "utils.h"
#include "nocache_async_stream_reader.h"

#define NOCACHE_ALIGNMENT 4096

struct uint40 {
  uint32_t    low;
  uint8_t     high;

  uint40() {}
  uint40(const long& a) : low(a & 0xFFFFFFFFL), high((a >> 32) & 0xFF) {}

  inline bool operator == (const uint40 &x) const {
    return low == x.low && high == x.high;
  }
} __attribute__((packed));

struct uint48 {
  uint32_t    low;
  uint16_t    high;

  uint48() {}
  uint48(const long& a) : low(a & 0xFFFFFFFFL), high((a >> 32) & 0xFFFF) {}

  inline bool operator == (const uint48 &x) const {
    return low == x.low && high == x.high;
  }
} __attribute__((packed));

struct uint56 {
  uint32_t    low;
  uint16_t    mid;
  uint8_t     high;

  uint56() {}
  uint56(const long& a) : low(a & 0xFFFFFFFFL), mid((a >> 32) & 0xFFFF), high((a >> 48) & 0xFF) {}

  inline bool operator == (const uint56 &x) const {
    return low == x.low && high == x.high && mid == x.mid;
  }
} __attribute__((packed));



template<typename value_type>
void test(value_type *tab, long length) {
  // Choose the filename as random.
  std::string filename = "testfile." + utils::random_string_hash();

  // Write tab into file.
  utils::write_objects_to_file(tab, length, filename);

  // Create the async stream reader.
  long bufsize = utils::random_long(1, 100000);  // in bytes
  typedef nocache_async_stream_reader<value_type, NOCACHE_ALIGNMENT> reader_type;
  reader_type *reader = new reader_type(filename, bufsize);

  // Read tab into filename using writer.
  value_type *tab_read = new value_type[length];
  for (long i = 0; i < length; ++i)
    tab_read[i] = reader->read();

  // Delete the reader.
  delete reader;

  if (!std::equal(tab, tab + length, tab_read)) {
    fprintf(stderr, "Error!\n");
    if (length < 100) {
      fprintf(stderr, "\tlength = %ld\n", length);
      /*fprintf(stderr, "\tOriginal tab: ");
      for (long i = 0; i < length; ++i)
        fprintf(stderr, "%ld ", (long)tab[i]);
      fprintf(stderr, "\n");

      fprintf(stderr, "Read tab: ");
      for (long i = 0; i < length; ++i)
        fprintf(stderr, "%ld ", (long)tab_read[i]);
      fprintf(stderr, "\n");*/
    }

    if (utils::file_exists(filename))
      utils::file_delete(filename);
    std::exit(EXIT_FAILURE);
  }

  if (utils::file_exists(filename))
    utils::file_delete(filename);
  delete[] tab_read;
}

template<typename value_type>
void test_random_type(int testcases, long max_length) {}

template<>
void test_random_type<int>(int testcases, long max_length) {
  fprintf(stderr,"TEST, testcases = %d, max_n = %ld, type = int\n", testcases, max_length);

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

template<>
void test_random_type<long>(int testcases, long max_length) {
  fprintf(stderr,"TEST, testcases = %d, max_n = %ld,  type = long\n", testcases, max_length);

  long *tab = new long[max_length];

  for (int tc = 0; tc < testcases; ++tc) {
    // Print progress information.
    if (tc % 10 == 0)
      fprintf(stderr,"%d (%.2Lf%%)\r", tc, (tc * 100.L) / testcases);

    // Generate input.
    long length = utils::random_long(1L, max_length);

    for (long i = 0; i < length; ++i)
      tab[i] = utils::random_long(1, 1000000000000000000L);

    // Run the test on generated input.
    test(tab, length);
  }

  // Clean up.
  delete[] tab;
}

template<>
void test_random_type<uint40>(int testcases, long max_length) {
  fprintf(stderr,"test, testcases = %d, max_n = %ld,  type = uint40\n", testcases, max_length);

  uint40 *tab = new uint40[max_length];

  for (int tc = 0; tc < testcases; ++tc) {
    // print progress information.
    if (tc % 10 == 0)
      fprintf(stderr,"%d (%.2Lf%%)\r", tc, (tc * 100.l) / testcases);

    // generate input.
    long length = utils::random_long(1l, max_length);

    for (long i = 0; i < length; ++i)
      tab[i] = utils::random_long(1, 1000000000000L);

    // run the test on generated input.
    test(tab, length);
  }

  // clean up.
  delete[] tab;
}

template<>
void test_random_type<uint48>(int testcases, long max_length) {
  fprintf(stderr,"test, testcases = %d, max_n = %ld,  type = uint48\n", testcases, max_length);

  uint48 *tab = new uint48[max_length];

  for (int tc = 0; tc < testcases; ++tc) {
    // print progress information.
    if (tc % 10 == 0)
      fprintf(stderr,"%d (%.2Lf%%)\r", tc, (tc * 100.l) / testcases);

    // generate input.
    long length = utils::random_long(1l, max_length);

    for (long i = 0; i < length; ++i)
      tab[i] = utils::random_long(1, 200000000000000L);

    // run the test on generated input.
    test(tab, length);
  }

  // clean up.
  delete[] tab;
}

template<>
void test_random_type<uint56>(int testcases, long max_length) {
  fprintf(stderr,"test, testcases = %d, max_n = %ld,  type = uint56\n", testcases, max_length);

  uint56 *tab = new uint56[max_length];

  for (int tc = 0; tc < testcases; ++tc) {
    // print progress information.
    if (tc % 10 == 0)
      fprintf(stderr,"%d (%.2Lf%%)\r", tc, (tc * 100.l) / testcases);

    // generate input.
    long length = utils::random_long(1l, max_length);

    for (long i = 0; i < length; ++i)
      tab[i] = utils::random_long(1, 60000000000000000L);

    // run the test on generated input.
    test(tab, length);
  }

  // clean up.
  delete[] tab;
}


template<typename value_type>
void test_random() {
  test_random_type<value_type>(100000, 100);
  test_random_type<value_type>(100000, 1000);
  test_random_type<value_type>(10000,  10000);
  test_random_type<value_type>(10000,  100000);
  test_random_type<value_type>(1000,   1000000);
  test_random_type<value_type>(100,    10000000);
}

int main() {
  std::srand(std::time(0) + getpid());

  test_random<int>();
  test_random<long>();
  test_random<uint40>();
  test_random<uint48>();
  test_random<uint56>();

  fprintf(stderr,"All tests passed.\n");
}

