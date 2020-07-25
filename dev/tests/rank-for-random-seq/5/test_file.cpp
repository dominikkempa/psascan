#include <ctime>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <unistd.h>
#include <string>

#include "utils.h"
#include "rank.h"
#include "oldrank.h"

// Test rank queries on a given strings.
void test(unsigned char *text, long length, long queries) {
  unsigned char *text_copy = new unsigned char[length + 1];
  std::copy(text, text + length, text_copy);
  text_copy[length] = 0;
  long max_threads = utils::random_long(1, 24);
  context_rank_4n *correct_rank = new context_rank_4n(text, length);
  rank4n<> *rank = new rank4n<>(text, length, max_threads);

  for (long q = 0; q < queries; ++q) {
    long i = utils::random_long(-2 * length, 2 * length);
    unsigned char c = utils::random_int(0, 255);

    // First, compute the correct answer.
    long ans = correct_rank->rank(i, c);
      
    // Now, query the tested data structured.
    long rank_ans = rank->rank(i, c);
    if (rank_ans != ans) {
      fprintf(stderr, "\n\033[22;31mFAILED\033[0m\n");
      fprintf(stderr, "  length = %ld\n", length);
      if (length <= 1000) {
        text[length] = 0;
        fprintf(stderr, "  text = ");
        for (long j = 0; j < length; ++j)
          fprintf(stderr, "%d ", text_copy[j]);
        fprintf(stderr, "\n");
      } else {
        fprintf(stderr, "  text = ");
        for (long j = 0; j < 100; ++j)
          fprintf(stderr, "%d ", text_copy[j]);
        fprintf(stderr, "...\n");

        FILE *fout = std::fopen("seq.txt", "w");
        fwrite(&length, sizeof(long), 1, fout);
        fwrite(text, sizeof(unsigned char), length, fout);
        std::fclose(fout);

      }
      fprintf(stderr, "  max_threads = %ld\n", max_threads);
      fprintf(stderr, "  i = %ld, c = %ld\n", i, (long)c);
      fprintf(stderr, "  correct = %ld, computed = %ld\n",
          ans, rank_ans);
      std::exit(EXIT_FAILURE);
    }
  }

  delete[] text_copy;
  delete rank;
  delete correct_rank;
}

void test_file(std::string filename) {
  std::FILE *f = std::fopen(filename.c_str(), "r");
  long length;
  long rrr = fread(&length, sizeof(long), 1, f);
  if (rrr != 1) {
    fprintf(stderr, "fread error\n");
    std::exit(EXIT_FAILURE);
  }
  unsigned char *text = new unsigned char[length];
  rrr = fread(text, sizeof(unsigned char), length, f);
  if (rrr != length) {
    fprintf(stderr, "fread error\n");
    std::exit(EXIT_FAILURE);
  }
  std::fclose(f);

  fprintf(stderr, "length = %ld\n", length);
  fprintf(stderr, "text = ");
  for (long i = 0; i < 100; ++i)
    fprintf(stderr, "%ld ", (long)text[i]);
  fprintf(stderr, "\n");

  test(text, length, 100000);
  delete[] text;
}

int main(int, char **argv) {
  test_file(argv[1]);
}
