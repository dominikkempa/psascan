/**
 * @file    src/main.cpp
 * @author  Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * @section LICENCE
 *
 * This file is part of pSAscan v0.1.0
 * See: http://www.cs.helsinki.fi/group/pads/
 *
 * Copyright (C) 2014-2015
 *   Juha Karkkainen <juha.karkkainen (at) cs.helsinki.fi>
 *   Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cstdint>
#include <string>
#include <sstream>
#include <getopt.h>
#include <unistd.h>
#include <omp.h>

#include "psascan_src/psascan.h"
#include "psascan_src/uint40.hpp"
#include "psascan_src/uint48.hpp"


char *program_name;

void usage(int status) {
  printf(
"Usage: %s [OPTION]... FILE\n"
"Construct the suffix array for text stored in FILE.\n"
"\n"
"Mandatory arguments to long options are mandatory for short options too.\n"
"  -g, --gap=GAPFILE       specify the file holding the gap array (default:\n"
"                          OUTFILE.gap, see the -o flag)\n"
"  -h, --help              display this help and exit\n"
"  -i, --intsize=SIZE      use integers of SIZE bytes to encode output.\n"
"                          Supported values 5, 6, and 8 (default: 5)\n"
"  -m, --mem=LIMIT         limit RAM usage to LIMIT MiB (default: 3072)\n"
"  -o, --output=OUTFILE    specify output filename (default: FILE.saX,\n"
"                          where X = integer size, see the -i flag)\n"
"  -v, --verbose           print detailed information during internal sufsort\n",
    program_name);

  std::exit(status);
}

bool file_exists(std::string fname) {
  std::FILE *f = std::fopen(fname.c_str(), "r");
  bool ret = (f != NULL);
  if (f != NULL) std::fclose(f);

  return ret;
}

template<typename int_type>
std::string intToStr(int_type x) {
  std::stringstream ss;
  ss << x;
  return ss.str();
}

int main(int argc, char **argv) {
  srand(time(0) + getpid());
  program_name = argv[0];
  bool verbose = false;

  static struct option long_options[] = {
    {"help",    no_argument,       NULL, 'h'},
    {"verbose", no_argument,       NULL, 'v'},
    {"intsize", required_argument, NULL, 'i'},
    {"mem",     required_argument, NULL, 'm'},
    {"output",  required_argument, NULL, 'o'},
    {"gap",     required_argument, NULL, 'g'},
    {NULL, 0, NULL, 0}
  };

  long ram_use = 3072L << 20;
  std::uint64_t int_size = 5;
  std::string out_filename("");
  std::string gap_filename("");

  // Parse command-line options.
  int c;
  while ((c = getopt_long(argc, argv, "g:hi:m:o:v", long_options, NULL)) != -1) {
    switch(c) {
      case 'i':
        int_size = std::atol(optarg);
        if (!(int_size == 5 || int_size == 6 || int_size == 8)) {
          fprintf(stderr, "\nError: invalid int size: %lu (see the -i flag).\n\n", int_size);
          usage(EXIT_FAILURE);
        }
        break;
      case 'm':
        ram_use = std::atol(optarg) << 20;
        if (ram_use <= 0L) {
          fprintf(stderr, "Error: invalid RAM limit (%ld)\n\n", ram_use);
          usage(EXIT_FAILURE);
        }
        break;
      case 'o':
        out_filename = std::string(optarg);
        break;
      case 'g':
        gap_filename = std::string(optarg);
        break;
      case 'v':
        verbose = true;
        break;
      case 'h':
        usage(EXIT_FAILURE);
      default:
        usage(EXIT_FAILURE);
    }
  }

  if (optind >= argc) {
    fprintf(stderr, "Error: FILE not provided\n\n");
    usage(EXIT_FAILURE);
  }

  // Parse the text filename.
  std::string text_filename = std::string(argv[optind++]);
  if (optind < argc) {
    fprintf(stderr, "Warning: multiple input files provided. "
    "Only the first will be processed.\n");
  }

  // Set default output filename (if not provided).
  if (out_filename.empty())
    out_filename = text_filename + ".sa" + intToStr(int_size);

  // Set default gap filename (if not provided).
  if (gap_filename.empty())
    gap_filename = out_filename;

  // Check if input exists.
  if (!file_exists(text_filename)) {
    fprintf(stderr, "Error: input file (%s) does not exist\n\n",
        text_filename.c_str());
    usage(EXIT_FAILURE);
  }

  if (file_exists(out_filename)) {
    // Output file exists, should we proceed?
    char *line = NULL;
    size_t buflen = 0;
    long len = 0L;

    do {
      printf("Output file (%s) exists. Overwrite? [y/n]: ",
          out_filename.c_str());
      if ((len = getline(&line, &buflen, stdin)) == -1) {
        fprintf(stderr, "\nError: failed to read answer\n\n");
        usage(EXIT_FAILURE);
      }
    } while (len != 2 || (line[0] != 'y' && line[0] != 'n'));

    if (line[0] == 'n') {
      free(line);
      std::exit(EXIT_FAILURE);
    }
    free(line);
  }

  // Find the number of (logical) cores on the machine.
  long max_threads = (long)omp_get_max_threads();

  // Run the pSAscan algorithm.
  if (int_size == 5) pSAscan<uint40>(text_filename, out_filename, gap_filename, ram_use, max_threads, verbose);
  else if (int_size == 6) pSAscan<uint48>(text_filename, out_filename, gap_filename, ram_use, max_threads, verbose);
  else pSAscan<std::uint64_t>(text_filename, out_filename, gap_filename, ram_use, max_threads, verbose);
}
