// Simply testing (but most importantly, already implementing for the final
// FGM) the buffered sparse gap array class. It optimizes random accesses to
// the array.
//
// Note: this optimization is not necessarily resulting in a speedup. If
// the random accesses can be somehow perform out-of-order, perhaps it would be
// better to interleave them than to do all at once.
//
// Note 2: would this trick be effective in other algorithms? Perhaps in some
// of these algorithms, the results of write operations are not necessary right
// away and could be done later 'in one batch', KKP2?
//
// Worth to note that this trick goes into the same group as the 'buffering'
// BWT in lzisa9 or buffering in Juha's implementation of induced sorting.  

#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <string>

#include "utils.h"
#include "buffered_gap_array.h"

void test(int length, int n_query) {
  int *gap_correct = new int[length];
  std::fill(gap_correct, gap_correct + length, 0);
  buffered_gap_array *gap = new buffered_gap_array(length);

  // Generate a sequence of integers to increment
  for (int i = 0; i < n_query; ++i) {
    int q = utils::random_int(0, length - 1);
    ++gap_correct[q];
    gap->increment(q);
  }
  
  // Stream and check the values from the buffered gap array.
  gap->flush();
  std::sort(gap->excess.begin(), gap->excess.end());
  for (int i = 0, pos = 0; i < length; ++i) {
    // Compute gap[i] from buffered gap array: we need to count
    // how many elements with value i is inside gap->excess.
    int c = 0;
    while (pos < (int)gap->excess.size() && gap->excess[pos] == i) ++pos, ++c;
    int gap_i = gap->count[i] + (c << 8);
    if (gap_correct[i] != gap_i) {
      fprintf(stderr, "Error!\n");
      fprintf(stderr, "  correct_gap[%d]=%d, buffered_gap[%d]=%d\n", i,
        gap_correct[i], i, gap_i);
      std::exit(EXIT_FAILURE);
    }
  }

  delete gap;
  delete[] gap_correct;
}

void test_random(int testcases, int max_length, int max_queries) {
  fprintf(stderr,"TEST, testcases = %d, max_length = %d, max_queries = %d\n",
      testcases, max_length, max_queries);
  for (int tc = 0, dbg = 0; tc < testcases; ++tc, ++dbg) {
    // Print progress information.
    if (dbg == 10) {
      fprintf(stderr,"%d (%.2Lf%%)\r", tc, (tc * 100.L) / testcases);
      dbg = 0;
    } 
    int length = utils::random_int(1, max_length);
    int queries = utils::random_int(1, max_queries);
    test(length, queries);
  }
}

int main(int, char **) {
  srand(time(0) + getpid());

  // Run tests.
  fprintf(stderr, "Testing the buffered gap array.\n");
  test_random(10000, 100, 10000);
  test_random(1000, 1000, 100000);
  test_random(100, 10000, 1000000);
  fprintf(stderr,"All tests passed.\n");

  return 0;
}

