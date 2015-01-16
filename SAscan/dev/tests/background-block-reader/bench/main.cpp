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

    for (long j = pos; j < pos + toread; ++j/*j += (1L << 20)*/)
      result_sum += reader.m_data[j];
 
    pos += toread;
  }
}


// Write the text to file. Then create multiple threads, that read
// the file (and compute, say, the sum of symbols on a string prefix)
// and compare the result to brute force. Threads doing the summing
// are using background_block_reader to access text symbols.
void test(std::string filename, long block_beg, long block_length) {
  fprintf(stderr, "Running benchmark on %s (block_beg = %ld, block_length = %ld)\n", filename.c_str(), block_beg, block_length);
  background_block_reader *bbr = new background_block_reader(filename, block_beg, block_length);

  // Generate the input.
  long n_threads = 16L;
  long *prefix_length = new long[n_threads];
  for (long j = 0; j < n_threads; ++j)
    prefix_length[j] = utils::random_long(1L, std::max(1L, block_length / 4L));

  // Run the threads.
  long *computed = new long[n_threads];
  std::thread **threads = new std::thread*[n_threads];
  for (long j = 0; j < n_threads; ++j)
    threads[j] = new std::thread(thread_body, std::ref(*bbr), prefix_length[j], std::ref(computed[j]));
  for (long j = 0; j < n_threads; ++j) threads[j]->join();
  bbr->stop();
  delete bbr;

  // Clean up.
  for (long j = 0; j < n_threads; ++j)
    delete threads[j];
  delete[] threads;

#if 0
  unsigned char *block = new unsigned char[block_length];
  utils::read_block(filename, block_beg, block_length, block);

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
  } else fprintf(stderr, "OK.\n");

  delete[] block;
  delete[] correct;
#endif

  delete[] prefix_length;
  delete[] computed;
}

int main(int argc, char **argv) {
  if (argc != 4) {
    fprintf(stderr, "Usage: %s FILE BLOCK_BEG BLOCK_SIZE.\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }


  std::string filename = std::string(argv[1]);
  long block_beg = atol(argv[2]);
  long block_size = atol(argv[3]);

  test(filename, block_beg, block_size);
}

