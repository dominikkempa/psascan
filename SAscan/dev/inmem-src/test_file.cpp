#include <cstdio>
#include <cstring>
#include <ctime>
#include <unistd.h>

#include "utils.h"
#include "inmem_sascan.h"
#include "divsufsort_template.h"



template<typename T>
void test(unsigned char *text, T text_length,
    long max_blocks, long max_threads) {
  long double start;

  fprintf(stderr, "Running inmem sascan\n");
  T *computed_sa = new T[text_length];
  start = utils::wclock();
  inmem_sascan(text, text_length, computed_sa, max_blocks, max_threads);
  fprintf(stderr, "Total time: %.2Lf\n", utils::wclock() - start);

  /*fprintf(stderr, "\nRunning divsufsort\n");
  T *correct_sa = new T[text_length];
  start = utils::wclock();
  run_divsufsort(text, correct_sa, text_length);
  fprintf(stderr, "Total time: %.2Lf\n", utils::wclock() - start);

  if (!std::equal(correct_sa, correct_sa + text_length, computed_sa))
    fprintf(stderr, "FAIL\n");
  else fprintf(stderr, "OK\n");

  delete[] correct_sa;*/
  delete[] computed_sa;
}


void test_file(const char *filename) {
  fprintf(stderr, "Reading text: ");
  long length;
  unsigned char *text;
  utils::read_objects_from_file(text, length, filename);
  fprintf(stderr, "DONE\n");

  test<long>(text, length, 24, 16);
//  test<long>(text, length, 24, 24);
//  test<long>(text, length, 24, 32);

  delete[] text;
}

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s <file>\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }

  test_file(argv[1]);
}
