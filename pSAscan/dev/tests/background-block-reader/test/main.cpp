#include <cstdio>
#include <algorithm>
#include <ctime>

#include <unistd.h>

#include "utils.h"
#include "background_block_reader.h"


// Each thread randomly selects its own segment size and
// sums up the symbols on the given prefix in segments.
//
// After each segment is processed, it waits until the
// backgroud reader reads the next segment.
void thread_body(background_block_reader &reader, long prefix_length, long &result_sum) {
  long segment_size = utils::random_long(1L, prefix_length);
  result_sum = 0L;

  long pos = 0L;
  while (pos < prefix_length) {
    long toread = std::min(prefix_length - pos, segment_size);
    long target_fetched = pos + toread;
    reader.wait(target_fetched);

    for (long j = pos; j < pos + toread; ++j)
      result_sum += reader.m_data[j];
 
    pos += toread;
  }
}


// Write the text to file. Then create multiple threads, that read
// the file (and compute, say, the sum of symbols on a string prefix)
// and compare the result to brute force. Threads doing the summing
// are using background_block_reader to access text symbols.
void test(unsigned char *text, long length) {
  std::string filename = "testfile." + utils::random_string_hash();
  utils::write_objects_to_file(text, length, filename);

  long block_length = utils::random_long(1L, length);
  long block_beg = utils::random_long(0L, length - block_length);
  unsigned char *block = text + block_beg;
  
  long chunk_size = utils::random_long(1L, block_length);
  background_block_reader *bbr = new background_block_reader(filename, block_beg, block_length, chunk_size);

  // XXX insert random sleep calls into the code of background_block_reader,
  // so that it does not immediatelly read the whole file.

  long n_threads = utils::random_long(1L, 50L);
  long *prefix_length = new long[n_threads];
  for (long j = 0; j < n_threads; ++j)
    prefix_length[j] = utils::random_long(1L, block_length);
  long *computed = new long[n_threads];
  std::thread **threads = new std::thread*[n_threads];

  for (long j = 0; j < n_threads; ++j)
    threads[j] = new std::thread(thread_body, std::ref(*bbr), prefix_length[j], std::ref(computed[j]));

  for (long j = 0; j < n_threads; ++j) threads[j]->join();
  bbr->stop();

  for (long j = 0; j < n_threads; ++j) delete threads[j];
  delete[] threads;

  // Compare.
  long *correct = new long[n_threads];
  for (long j = 0; j < n_threads; ++j) {
    correct[j] = 0;
    for (long jj = 0; jj < prefix_length[j]; ++jj)
      correct[j] += block[jj];
  }

  if (!std::equal(correct, correct + n_threads, computed)) {
    fprintf(stderr, "\nError!\n");
    std::exit(EXIT_FAILURE);
  }

  delete bbr;
  delete[] prefix_length;
  delete[] computed;
  delete[] correct;
  utils::file_delete(filename);
}

void test_random(int testcases, long max_length) {
  fprintf(stderr,"TEST, testcases = %d, max_length = %ld\n", testcases, max_length);

  unsigned char *text = new unsigned char[max_length];

  for (int tc = 0; tc < testcases; ++tc) {
    // Print progress information.
    if (tc % 10 == 0)
      fprintf(stderr,"%d (%.2Lf%%)\r", tc, (tc * 100.L) / testcases);

    // Generate input.
    long length = utils::random_long(1L, max_length);
    long sigma = utils::random_long(1L, 26L);

    utils::fill_random_letters(text, length, sigma);

    // Run the test on generated input.
    test(text, length);
  }

  // Clean up.
  delete[] text;
}

int main() {
  std::srand(std::time(0) + getpid());

  test_random(100000,  10);
  test_random(100000,  100);
  test_random(10000,   1000);
  test_random(10000,   10000);
  test_random(1000,    100000);
  test_random(100,     1000000);
  test_random(10,      10000000);

  fprintf(stderr,"All tests passed.\n");
}

