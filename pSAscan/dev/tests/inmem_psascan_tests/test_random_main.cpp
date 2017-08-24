#include <cstdio>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "gt_test_random.hpp"
#include "pgt_test_random.hpp"
#include "sa_test_random.hpp"
#include "psa_test_random.hpp"
#include "bwt_test_random.hpp"
#include "pbwt_test_random.hpp"


int main() {
  std::srand(std::time(0) + getpid());

  // Redirect stdout to /dev/null
  int redir = open("/dev/null", O_WRONLY);
  dup2(redir, 2);
  close(redir);

  fprintf(stdout, "gt test random\n");

#ifndef NDEBUG
  gt_test_random_private::test_random<7> (100,  10, 255);
  gt_test_random_private::test_random<12>(100,  10, 255);
  gt_test_random_private::test_random<7> (100, 100, 255);
  gt_test_random_private::test_random<12>(100, 100, 255);
  gt_test_random_private::test_random<7> (10, 1000, 255);
  gt_test_random_private::test_random<12>(10, 1000, 255);
  gt_test_random_private::test_random<7> (5, 10000, 255);
  gt_test_random_private::test_random<12>(5, 10000, 255);

#else
  gt_test_random_private::test_random<7> (50, 1000, 5);
  gt_test_random_private::test_random<12>(50, 1000, 5);
  gt_test_random_private::test_random<7> (50, 1000, 20);
  gt_test_random_private::test_random<12>(50, 1000, 20);
  gt_test_random_private::test_random<7> (50, 1000, 128);
  gt_test_random_private::test_random<12>(50, 1000, 128);
  gt_test_random_private::test_random<7> (50, 1000, 255);
  gt_test_random_private::test_random<12>(50, 1000, 255);
  
  gt_test_random_private::test_random<7> (10, 10000, 5);
  gt_test_random_private::test_random<12>(10, 10000, 5);
  gt_test_random_private::test_random<7> (10, 10000, 20);
  gt_test_random_private::test_random<12>(10, 10000, 20);
  gt_test_random_private::test_random<7> (10, 10000, 128);
  gt_test_random_private::test_random<12>(10, 10000, 128);
  gt_test_random_private::test_random<7> (10, 10000, 255);
  gt_test_random_private::test_random<12>(10, 10000, 255);

  gt_test_random_private::test_random<7> (5, 1000000, 5);
  gt_test_random_private::test_random<12>(5, 1000000, 5);
  gt_test_random_private::test_random<7> (5, 1000000, 20);
  gt_test_random_private::test_random<12>(5, 1000000, 20);
  gt_test_random_private::test_random<7> (5, 1000000, 128);
  gt_test_random_private::test_random<12>(5, 1000000, 128);
  gt_test_random_private::test_random<7> (5, 1000000, 255);
  gt_test_random_private::test_random<12>(5, 1000000, 255);
#endif

  fprintf(stdout, "All tests passed.\n\n");
  fflush(stdout);

  fprintf(stdout, "pgt test random\n");

#ifndef NDEBUG
  pgt_test_random_private::test_random<uint40, 2>(10,   10,    255);
  pgt_test_random_private::test_random<uint40, 5>(10,   10,    255);
  pgt_test_random_private::test_random<uint40, 8>(10,   10,    255);
  pgt_test_random_private::test_random<int,    2>(10,   10,    255);
  pgt_test_random_private::test_random<int,    5>(10,   10,    255);
  pgt_test_random_private::test_random<int,    8>(10,   10,    255);

  pgt_test_random_private::test_random<uint40, 2>(10,   100,    255);
  pgt_test_random_private::test_random<uint40, 5>(10,   100,    255);
  pgt_test_random_private::test_random<uint40, 8>(10,   100,    255);
  pgt_test_random_private::test_random<int,    2>(10,   100,    255);
  pgt_test_random_private::test_random<int,    5>(10,   100,    255);
  pgt_test_random_private::test_random<int,    8>(10,   100,    255);

  pgt_test_random_private::test_random<uint40, 2>(5,   1000,    255);
  pgt_test_random_private::test_random<uint40, 5>(5,   1000,    255);
  pgt_test_random_private::test_random<uint40, 8>(5,   1000,    255);
  pgt_test_random_private::test_random<int,    2>(5,   1000,    255);
  pgt_test_random_private::test_random<int,    5>(5,   1000,    255);
  pgt_test_random_private::test_random<int,    8>(5,   1000,    255);

  pgt_test_random_private::test_random<uint40, 2>(5,   10000,    255);
  pgt_test_random_private::test_random<uint40, 5>(5,   10000,    255);
  pgt_test_random_private::test_random<uint40, 8>(5,   10000,    255);
  pgt_test_random_private::test_random<int,    2>(5,   10000,    255);
  pgt_test_random_private::test_random<int,    5>(5,   10000,    255);
  pgt_test_random_private::test_random<int,    8>(5,   10000,    255);

#else
  pgt_test_random_private::test_random<uint40, 2>(1000,   10,      5);
  pgt_test_random_private::test_random<uint40, 5>(1000,   10,      5);
  pgt_test_random_private::test_random<uint40, 8>(1000,   10,      5);
  pgt_test_random_private::test_random<uint40, 2>(1000,   10,    255);
  pgt_test_random_private::test_random<uint40, 5>(1000,   10,    255);
  pgt_test_random_private::test_random<uint40, 8>(1000,   10,    255);
  pgt_test_random_private::test_random<int,    2>(1000,   10,      5);
  pgt_test_random_private::test_random<int,    5>(1000,   10,      5);
  pgt_test_random_private::test_random<int,    8>(1000,   10,      5);
  pgt_test_random_private::test_random<int,    2>(1000,   10,    255);
  pgt_test_random_private::test_random<int,    5>(1000,   10,    255);
  pgt_test_random_private::test_random<int,    8>(1000,   10,    255);

  pgt_test_random_private::test_random<uint40, 2>(100,   100,      5);
  pgt_test_random_private::test_random<uint40, 5>(100,   100,      5);
  pgt_test_random_private::test_random<uint40, 8>(100,   100,      5);
  pgt_test_random_private::test_random<uint40, 2>(100,   100,    255);
  pgt_test_random_private::test_random<uint40, 5>(100,   100,    255);
  pgt_test_random_private::test_random<uint40, 8>(100,   100,    255);
  pgt_test_random_private::test_random<int,    2>(100,   100,      5);
  pgt_test_random_private::test_random<int,    5>(100,   100,      5);
  pgt_test_random_private::test_random<int,    8>(100,   100,      5);
  pgt_test_random_private::test_random<int,    2>(100,   100,    255);
  pgt_test_random_private::test_random<int,    5>(100,   100,    255);
  pgt_test_random_private::test_random<int,    8>(100,   100,    255);

  pgt_test_random_private::test_random<uint40, 2>(30,   1000,      5);
  pgt_test_random_private::test_random<uint40, 5>(30,   1000,      5);
  pgt_test_random_private::test_random<uint40, 8>(30,   1000,      5);
  pgt_test_random_private::test_random<uint40, 2>(30,   1000,    255);
  pgt_test_random_private::test_random<uint40, 5>(30,   1000,    255);
  pgt_test_random_private::test_random<uint40, 8>(30,   1000,    255);
  pgt_test_random_private::test_random<int,    2>(30,   1000,      5);
  pgt_test_random_private::test_random<int,    5>(30,   1000,      5);
  pgt_test_random_private::test_random<int,    8>(30,   1000,      5);
  pgt_test_random_private::test_random<int,    2>(30,   1000,    255);
  pgt_test_random_private::test_random<int,    5>(30,   1000,    255);
  pgt_test_random_private::test_random<int,    8>(30,   1000,    255);

  pgt_test_random_private::test_random<uint40, 2>(4,   1000000,      5);
  pgt_test_random_private::test_random<uint40, 5>(4,   1000000,      5);
  pgt_test_random_private::test_random<uint40, 8>(4,   1000000,      5);
  pgt_test_random_private::test_random<uint40, 2>(4,   1000000,    255);
  pgt_test_random_private::test_random<uint40, 5>(4,   1000000,    255);
  pgt_test_random_private::test_random<uint40, 8>(4,   1000000,    255);
  pgt_test_random_private::test_random<int,    2>(4,   1000000,      5);
  pgt_test_random_private::test_random<int,    5>(4,   1000000,      5);
  pgt_test_random_private::test_random<int,    8>(4,   1000000,      5);
  pgt_test_random_private::test_random<int,    2>(4,   1000000,    255);
  pgt_test_random_private::test_random<int,    5>(4,   1000000,    255);
  pgt_test_random_private::test_random<int,    8>(4,   1000000,    255);
#endif

  fprintf(stdout,"All tests passed.\n\n");
  fflush(stdout);

  fprintf(stdout, "sa test random\n");
  
#ifndef NDEBUG
  sa_test_random_private::test_random<2> (20, 100, 255);
  sa_test_random_private::test_random<5> (20, 100, 255);
  sa_test_random_private::test_random<7> (20, 100, 255);
  sa_test_random_private::test_random<2> (10, 1000, 255);
  sa_test_random_private::test_random<5> (10, 1000, 255);
  sa_test_random_private::test_random<7> (10, 1000, 255);
  sa_test_random_private::test_random<2> (5, 10000, 255);
  sa_test_random_private::test_random<5> (5, 10000, 255);
  sa_test_random_private::test_random<7> (5, 10000, 255);

#else
  sa_test_random_private::test_random<1> (30, 300, 5);
  sa_test_random_private::test_random<1> (30, 300, 128);
  sa_test_random_private::test_random<3> (30, 300, 5);
  sa_test_random_private::test_random<3> (30, 300, 128);
  sa_test_random_private::test_random<5> (30, 300, 5);
  sa_test_random_private::test_random<5> (30, 300, 128);
  sa_test_random_private::test_random<12>(30, 300, 5);
  sa_test_random_private::test_random<12>(30, 300, 128);
  sa_test_random_private::test_random<1> (30, 300, 255);
  sa_test_random_private::test_random<3> (30, 300, 255);
  sa_test_random_private::test_random<5> (30, 300, 255);
  sa_test_random_private::test_random<12>(30, 300, 255);

  sa_test_random_private::test_random<1> (10, 1000, 5);
  sa_test_random_private::test_random<1> (10, 1000, 128);
  sa_test_random_private::test_random<3> (10, 1000, 5);
  sa_test_random_private::test_random<3> (10, 1000, 128);
  sa_test_random_private::test_random<5> (10, 1000, 5);
  sa_test_random_private::test_random<5> (10, 1000, 128);
  sa_test_random_private::test_random<12>(10, 1000, 5);
  sa_test_random_private::test_random<12>(10, 1000, 128);
  sa_test_random_private::test_random<1> (10, 1000, 255);
  sa_test_random_private::test_random<3> (10, 1000, 255);
  sa_test_random_private::test_random<5> (10, 1000, 255);
  sa_test_random_private::test_random<12>(10, 1000, 255);

  sa_test_random_private::test_random<7> (20,  10000,   5);
  sa_test_random_private::test_random<12>(20,  10000,   5);
  sa_test_random_private::test_random<7> (20,  10000,   20);
  sa_test_random_private::test_random<12>(20,  10000,   20);
  sa_test_random_private::test_random<7> (20,  10000,   128);
  sa_test_random_private::test_random<12>(20,  10000,   128);
  sa_test_random_private::test_random<7> (20,  10000,   255);
  sa_test_random_private::test_random<12>(20,  10000,   255);

  sa_test_random_private::test_random<7> (4,  1000000, 5);
  sa_test_random_private::test_random<12>(4,  1000000, 5);
  sa_test_random_private::test_random<7> (4,  1000000, 20);
  sa_test_random_private::test_random<12>(4,  1000000, 20);
  sa_test_random_private::test_random<7> (4,  1000000, 128);
  sa_test_random_private::test_random<12>(4,  1000000, 128);
  sa_test_random_private::test_random<7> (4,  1000000, 255);
  sa_test_random_private::test_random<12>(4,  1000000, 255);
#endif

  fprintf(stdout, "All tests passed.\n\n");
  fflush(stdout);

  fprintf(stdout, "psa test random\n");

#ifndef NDEBUG
  psa_test_random_private::test_random<uint40, 2>(10,   10,    255);
  psa_test_random_private::test_random<uint40, 5>(10,   10,    255);
  psa_test_random_private::test_random<uint40, 8>(10,   10,    255);
  psa_test_random_private::test_random<int,    2>(10,   10,    255);
  psa_test_random_private::test_random<int,    5>(10,   10,    255);
  psa_test_random_private::test_random<int,    8>(10,   10,    255);

  psa_test_random_private::test_random<uint40, 2>(10,   50,    255);
  psa_test_random_private::test_random<uint40, 5>(10,   50,    255);
  psa_test_random_private::test_random<uint40, 8>(10,   50,    255);
  psa_test_random_private::test_random<int,    2>(10,   50,    255);
  psa_test_random_private::test_random<int,    5>(10,   50,    255);
  psa_test_random_private::test_random<int,    8>(10,   50,    255);

  psa_test_random_private::test_random<uint40, 2>(10,   100,    255);
  psa_test_random_private::test_random<uint40, 5>(10,   100,    255);
  psa_test_random_private::test_random<uint40, 8>(10,   100,    255);
  psa_test_random_private::test_random<int,    2>(10,   100,    255);
  psa_test_random_private::test_random<int,    5>(10,   100,    255);
  psa_test_random_private::test_random<int,    8>(10,   100,    255);

  psa_test_random_private::test_random<uint40, 2>(5,   1000,    255);
  psa_test_random_private::test_random<uint40, 5>(5,   1000,    255);
  psa_test_random_private::test_random<uint40, 8>(5,   1000,    255);
  psa_test_random_private::test_random<int,    2>(5,   1000,    255);
  psa_test_random_private::test_random<int,    5>(5,   1000,    255);
  psa_test_random_private::test_random<int,    8>(5,   1000,    255);

  psa_test_random_private::test_random<int, 1>(2, 1000, 255);
  psa_test_random_private::test_random<int, 3>(2, 1000, 255);
  psa_test_random_private::test_random<int, 5>(2, 1000, 255);
  psa_test_random_private::test_random<int, 12>(2, 1000, 255);

  psa_test_random_private::test_random<int,    2>(2,   10000,    255);
  psa_test_random_private::test_random<int,    5>(2,   10000,    255);
  psa_test_random_private::test_random<int,    8>(2,   10000,    255);

#else
  psa_test_random_private::test_random<uint40, 2>(1000,   10,      5);
  psa_test_random_private::test_random<uint40, 5>(1000,   10,      5);
  psa_test_random_private::test_random<uint40, 8>(1000,   10,      5);
  psa_test_random_private::test_random<uint40, 2>(1000,   10,    255);
  psa_test_random_private::test_random<uint40, 5>(1000,   10,    255);
  psa_test_random_private::test_random<uint40, 8>(1000,   10,    255);
  psa_test_random_private::test_random<int,    2>(1000,   10,      5);
  psa_test_random_private::test_random<int,    5>(1000,   10,      5);
  psa_test_random_private::test_random<int,    8>(1000,   10,      5);
  psa_test_random_private::test_random<int,    2>(1000,   10,    255);
  psa_test_random_private::test_random<int,    5>(1000,   10,    255);
  psa_test_random_private::test_random<int,    8>(1000,   10,    255);

  psa_test_random_private::test_random<uint40, 2>(200,   50,      5);
  psa_test_random_private::test_random<uint40, 5>(200,   50,      5);
  psa_test_random_private::test_random<uint40, 8>(200,   50,      5);
  psa_test_random_private::test_random<uint40, 2>(200,   50,    255);
  psa_test_random_private::test_random<uint40, 5>(200,   50,    255);
  psa_test_random_private::test_random<uint40, 8>(200,   50,    255);
  psa_test_random_private::test_random<int,    2>(200,   50,      5);
  psa_test_random_private::test_random<int,    5>(200,   50,      5);
  psa_test_random_private::test_random<int,    8>(200,   50,      5);
  psa_test_random_private::test_random<int,    2>(200,   50,    255);
  psa_test_random_private::test_random<int,    5>(200,   50,    255);
  psa_test_random_private::test_random<int,    8>(200,   50,    255);

  psa_test_random_private::test_random<uint40, 2>(100,   100,      5);
  psa_test_random_private::test_random<uint40, 5>(100,   100,      5);
  psa_test_random_private::test_random<uint40, 8>(100,   100,      5);
  psa_test_random_private::test_random<uint40, 2>(100,   100,    255);
  psa_test_random_private::test_random<uint40, 5>(100,   100,    255);
  psa_test_random_private::test_random<uint40, 8>(100,   100,    255);
  psa_test_random_private::test_random<int,    2>(100,   100,      5);
  psa_test_random_private::test_random<int,    5>(100,   100,      5);
  psa_test_random_private::test_random<int,    8>(100,   100,      5);
  psa_test_random_private::test_random<int,    2>(100,   100,    255);
  psa_test_random_private::test_random<int,    5>(100,   100,    255);
  psa_test_random_private::test_random<int,    8>(100,   100,    255);

  psa_test_random_private::test_random<uint40, 2>(20,   1000,      5);
  psa_test_random_private::test_random<uint40, 5>(20,   1000,      5);
  psa_test_random_private::test_random<uint40, 8>(20,   1000,      5);
  psa_test_random_private::test_random<uint40, 2>(20,   1000,    255);
  psa_test_random_private::test_random<uint40, 5>(20,   1000,    255);
  psa_test_random_private::test_random<uint40, 8>(20,   1000,    255);
  psa_test_random_private::test_random<int,    2>(20,   1000,      5);
  psa_test_random_private::test_random<int,    5>(20,   1000,      5);
  psa_test_random_private::test_random<int,    8>(20,   1000,      5);
  psa_test_random_private::test_random<int,    2>(20,   1000,    255);
  psa_test_random_private::test_random<int,    5>(20,   1000,    255);
  psa_test_random_private::test_random<int,    8>(20,   1000,    255);

  psa_test_random_private::test_random<int, 1>(100, 300, 5);
  psa_test_random_private::test_random<int, 1>(100, 300, 128);
  psa_test_random_private::test_random<int, 3>(100, 300, 5);
  psa_test_random_private::test_random<int, 3>(100, 300, 128);
  psa_test_random_private::test_random<int, 5>(100, 300, 5);
  psa_test_random_private::test_random<int, 5>(100, 300, 128);
  psa_test_random_private::test_random<int, 12>(100, 300, 5);
  psa_test_random_private::test_random<int, 12>(100, 300, 128);
  psa_test_random_private::test_random<int, 1>(100, 300, 255);
  psa_test_random_private::test_random<int, 3>(100, 300, 255);
  psa_test_random_private::test_random<int, 5>(100, 300, 255);
  psa_test_random_private::test_random<int, 12>(100, 300, 255);

  psa_test_random_private::test_random<int, 1>(30, 1000, 5);
  psa_test_random_private::test_random<int, 1>(30, 1000, 128);
  psa_test_random_private::test_random<int, 3>(30, 1000, 5);
  psa_test_random_private::test_random<int, 3>(30, 1000, 128);
  psa_test_random_private::test_random<int, 5>(30, 1000, 5);
  psa_test_random_private::test_random<int, 5>(30, 1000, 128);
  psa_test_random_private::test_random<int, 12>(30, 1000, 5);
  psa_test_random_private::test_random<int, 12>(30, 1000, 128);
  psa_test_random_private::test_random<int, 1>(30, 1000, 255);
  psa_test_random_private::test_random<int, 3>(30, 1000, 255);
  psa_test_random_private::test_random<int, 5>(30, 1000, 255);
  psa_test_random_private::test_random<int, 12>(30, 1000, 255);

  psa_test_random_private::test_random<uint40, 2>(4,   1000000,      5);
  psa_test_random_private::test_random<uint40, 5>(4,   1000000,      5);
  psa_test_random_private::test_random<uint40, 8>(4,   1000000,      5);
  psa_test_random_private::test_random<uint40, 2>(4,   1000000,    255);
  psa_test_random_private::test_random<uint40, 5>(4,   1000000,    255);
  psa_test_random_private::test_random<uint40, 8>(4,   1000000,    255);
  psa_test_random_private::test_random<int,    2>(4,   1000000,      5);
  psa_test_random_private::test_random<int,    5>(4,   1000000,      5);
  psa_test_random_private::test_random<int,    8>(4,   1000000,      5);
  psa_test_random_private::test_random<int,    2>(4,   1000000,    255);
  psa_test_random_private::test_random<int,    5>(4,   1000000,    255);
  psa_test_random_private::test_random<int,    8>(4,   1000000,    255);
#endif

  fprintf(stdout,"All tests passed.\n\n");
  std::fflush(stdout);

  fprintf(stdout, "bwt test random\n");

#ifndef NDEBUG
  bwt_test_random_private::test_random<7> (30, 100, 255);
  bwt_test_random_private::test_random<12>(30, 100, 255);
  bwt_test_random_private::test_random<7> (10, 1000, 255);
  bwt_test_random_private::test_random<12>(10, 1000, 255);  
  bwt_test_random_private::test_random<7> (5, 10000, 255);
  bwt_test_random_private::test_random<12>(5, 10000, 255);

#else
  bwt_test_random_private::test_random<7> (1000, 10, 5);
  bwt_test_random_private::test_random<12>(1000, 10, 5);
  bwt_test_random_private::test_random<7> (1000, 10, 20);
  bwt_test_random_private::test_random<12>(1000, 10, 20);
  bwt_test_random_private::test_random<7> (1000, 10, 128);
  bwt_test_random_private::test_random<12>(1000, 10, 128);
  bwt_test_random_private::test_random<7> (1000, 10, 255);
  bwt_test_random_private::test_random<12>(1000, 10, 255);

  bwt_test_random_private::test_random<7> (500, 100, 5);
  bwt_test_random_private::test_random<12>(500, 100, 5);
  bwt_test_random_private::test_random<7> (500, 100, 20);
  bwt_test_random_private::test_random<12>(500, 100, 20);
  bwt_test_random_private::test_random<7> (500, 100, 128);
  bwt_test_random_private::test_random<12>(500, 100, 128);
  bwt_test_random_private::test_random<7> (500, 100, 255);
  bwt_test_random_private::test_random<12>(500, 100, 255);

  bwt_test_random_private::test_random<7> (30, 1000, 5);
  bwt_test_random_private::test_random<12>(30, 1000, 5);
  bwt_test_random_private::test_random<7> (30, 1000, 20);
  bwt_test_random_private::test_random<12>(30, 1000, 20);
  bwt_test_random_private::test_random<7> (30, 1000, 128);
  bwt_test_random_private::test_random<12>(30, 1000, 128);
  bwt_test_random_private::test_random<7> (30, 1000, 255);
  bwt_test_random_private::test_random<12>(30, 1000, 255);
  
  bwt_test_random_private::test_random<7> (5, 10000, 5);
  bwt_test_random_private::test_random<12>(5, 10000, 5);
  bwt_test_random_private::test_random<7> (5, 10000, 20);
  bwt_test_random_private::test_random<12>(5, 10000, 20);
  bwt_test_random_private::test_random<7> (5, 10000, 128);
  bwt_test_random_private::test_random<12>(5, 10000, 128);
  bwt_test_random_private::test_random<7> (5, 10000, 255);
  bwt_test_random_private::test_random<12>(5, 10000, 255);

  bwt_test_random_private::test_random<7> (4, 1000000, 5);
  bwt_test_random_private::test_random<12>(4, 1000000, 5);
  bwt_test_random_private::test_random<7> (4, 1000000, 20);
  bwt_test_random_private::test_random<12>(4, 1000000, 20);
  bwt_test_random_private::test_random<7> (4, 1000000, 128);
  bwt_test_random_private::test_random<12>(4, 1000000, 128);
  bwt_test_random_private::test_random<7> (4, 1000000, 255);
  bwt_test_random_private::test_random<12>(4, 1000000, 255);
#endif

  fprintf(stdout, "All tests passed.\n\n");
  fflush(stdout);

  fprintf(stdout, "pbwt test random\n");

#ifndef NDEBUG
  pbwt_test_random_private::test_random<uint40, 2>(20,   10,    255);
  pbwt_test_random_private::test_random<uint40, 5>(20,   10,    255);
  pbwt_test_random_private::test_random<uint40, 8>(20,   10,    255);
  pbwt_test_random_private::test_random<int,    2>(20,   10,    255);
  pbwt_test_random_private::test_random<int,    5>(20,   10,    255);
  pbwt_test_random_private::test_random<int,    8>(20,   10,    255);

  pbwt_test_random_private::test_random<uint40, 2>(10,   100,    255);
  pbwt_test_random_private::test_random<uint40, 5>(10,   100,    255);
  pbwt_test_random_private::test_random<uint40, 8>(10,   100,    255);
  pbwt_test_random_private::test_random<int,    2>(10,   100,    255);
  pbwt_test_random_private::test_random<int,    5>(10,   100,    255);
  pbwt_test_random_private::test_random<int,    8>(10,   100,    255);

  pbwt_test_random_private::test_random<uint40, 2>(4,   1000,    255);
  pbwt_test_random_private::test_random<uint40, 5>(4,   1000,    255);
  pbwt_test_random_private::test_random<uint40, 8>(4,   1000,    255);
  pbwt_test_random_private::test_random<int,    2>(4,   1000,    255);
  pbwt_test_random_private::test_random<int,    5>(4,   1000,    255);
  pbwt_test_random_private::test_random<int,    8>(4,   1000,    255);

  pbwt_test_random_private::test_random<uint40, 2>(2,   10000,    255);
  pbwt_test_random_private::test_random<uint40, 5>(2,   10000,    255);
  pbwt_test_random_private::test_random<uint40, 8>(2,   10000,    255);
  pbwt_test_random_private::test_random<int,    2>(2,   10000,    255);
  pbwt_test_random_private::test_random<int,    5>(2,   10000,    255);
  pbwt_test_random_private::test_random<int,    8>(2,   10000,    255);

#else
  pbwt_test_random_private::test_random<uint40, 2>(1000,   10,      5);
  pbwt_test_random_private::test_random<uint40, 5>(1000,   10,      5);
  pbwt_test_random_private::test_random<uint40, 8>(1000,   10,      5);
  pbwt_test_random_private::test_random<uint40, 2>(1000,   10,    255);
  pbwt_test_random_private::test_random<uint40, 5>(1000,   10,    255);
  pbwt_test_random_private::test_random<uint40, 8>(1000,   10,    255);
  pbwt_test_random_private::test_random<int,    2>(1000,   10,      5);
  pbwt_test_random_private::test_random<int,    5>(1000,   10,      5);
  pbwt_test_random_private::test_random<int,    8>(1000,   10,      5);
  pbwt_test_random_private::test_random<int,    2>(1000,   10,    255);
  pbwt_test_random_private::test_random<int,    5>(1000,   10,    255);
  pbwt_test_random_private::test_random<int,    8>(1000,   10,    255);

  pbwt_test_random_private::test_random<uint40, 2>(100,   100,      5);
  pbwt_test_random_private::test_random<uint40, 5>(100,   100,      5);
  pbwt_test_random_private::test_random<uint40, 8>(100,   100,      5);
  pbwt_test_random_private::test_random<uint40, 2>(100,   100,    255);
  pbwt_test_random_private::test_random<uint40, 5>(100,   100,    255);
  pbwt_test_random_private::test_random<uint40, 8>(100,   100,    255);
  pbwt_test_random_private::test_random<int,    2>(100,   100,      5);
  pbwt_test_random_private::test_random<int,    5>(100,   100,      5);
  pbwt_test_random_private::test_random<int,    8>(100,   100,      5);
  pbwt_test_random_private::test_random<int,    2>(100,   100,    255);
  pbwt_test_random_private::test_random<int,    5>(100,   100,    255);
  pbwt_test_random_private::test_random<int,    8>(100,   100,    255);

  pbwt_test_random_private::test_random<uint40, 2>(50,   1000,      5);
  pbwt_test_random_private::test_random<uint40, 5>(50,   1000,      5);
  pbwt_test_random_private::test_random<uint40, 8>(50,   1000,      5);
  pbwt_test_random_private::test_random<uint40, 2>(50,   1000,    255);
  pbwt_test_random_private::test_random<uint40, 5>(50,   1000,    255);
  pbwt_test_random_private::test_random<uint40, 8>(50,   1000,    255);
  pbwt_test_random_private::test_random<int,    2>(50,   1000,      5);
  pbwt_test_random_private::test_random<int,    5>(50,   1000,      5);
  pbwt_test_random_private::test_random<int,    8>(50,   1000,      5);
  pbwt_test_random_private::test_random<int,    2>(50,   1000,    255);
  pbwt_test_random_private::test_random<int,    5>(50,   1000,    255);
  pbwt_test_random_private::test_random<int,    8>(50,   1000,    255);

  pbwt_test_random_private::test_random<uint40, 2>(4,   1000000,      5);
  pbwt_test_random_private::test_random<uint40, 5>(4,   1000000,      5);
  pbwt_test_random_private::test_random<uint40, 8>(4,   1000000,      5);
  pbwt_test_random_private::test_random<uint40, 2>(4,   1000000,    255);
  pbwt_test_random_private::test_random<uint40, 5>(4,   1000000,    255);
  pbwt_test_random_private::test_random<uint40, 8>(4,   1000000,    255);
  pbwt_test_random_private::test_random<int,    2>(4,   1000000,      5);
  pbwt_test_random_private::test_random<int,    5>(4,   1000000,      5);
  pbwt_test_random_private::test_random<int,    8>(4,   1000000,      5);
  pbwt_test_random_private::test_random<int,    2>(4,   1000000,    255);
  pbwt_test_random_private::test_random<int,    5>(4,   1000000,    255);
  pbwt_test_random_private::test_random<int,    8>(4,   1000000,    255);
#endif

  fprintf(stdout,"All tests passed.\n");
  fflush(stdout);

  utils::file_delete("supertext.txt");
}

