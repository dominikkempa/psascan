/**
 * @file    src/psascan_src/utils.cpp
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
#include <cstdint>
#include <cstring>
#include <errno.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/time.h>
#include <string>
#include <fstream>
#include <algorithm>

#include "utils.h"


namespace psascan_private {
namespace utils {

long double wclock() {
  timeval tim;
  gettimeofday(&tim, NULL);

  return tim.tv_sec + (tim.tv_usec / 1000000.0L);
}

std::FILE *open_file(std::string fname, std::string mode) {
  std::FILE *f = std::fopen(fname.c_str(), mode.c_str());
  if (f == NULL) {
    std::perror(fname.c_str());
    std::exit(EXIT_FAILURE);
  }

  return f;
}

std::uint64_t file_size(std::string fname) {
  std::FILE *f = open_file(fname, "rt");
  std::fseek(f, 0L, SEEK_END);
  long size = std::ftell(f);
  if (size < 0) {
    fprintf(stderr, "\nError: cannot find file size of %s.\n", fname.c_str());
    std::exit(EXIT_FAILURE);
  }
  std::fclose(f);

  return (std::uint64_t)size;
}

bool file_exists(std::string fname) {
  std::FILE *f = std::fopen(fname.c_str(), "r");
  bool ret = (f != NULL);
  if (f != NULL)
    std::fclose(f);

  return ret;
}

void file_delete(std::string fname) {
  int res = std::remove(fname.c_str());
  if (res) {
    fprintf(stderr, "Failed to delete %s: %s\n",
        fname.c_str(), strerror(errno));
    std::exit(EXIT_FAILURE);
  }
}

std::string absolute_path(std::string fname) {
  char path[1 << 12];
  bool created = false;

  if (!file_exists(fname)) {
    // We need to create the file, since realpath fails on non-existing files.
    std::fclose(open_file(fname, "w"));
    created = true;
  }
  if (!realpath(fname.c_str(), path)) {
    fprintf(stderr, "\nError: realpath failed for %s\n", fname.c_str());
    std::exit(EXIT_FAILURE);
  }

  if (created)
    file_delete(fname);

  return std::string(path);
}

void read_block(std::FILE *f, std::uint64_t beg, std::uint64_t length, unsigned char *b) {
  std::fseek(f, beg, SEEK_SET);
  read_n_objects_from_file<unsigned char>(b, length, f);
}

void read_block(std::string fname, std::uint64_t beg, std::uint64_t length, unsigned char *b) {
  std::FILE *f = open_file(fname.c_str(), "r");
  read_block(f, beg, length, b);
  std::fclose(f);
}

std::int32_t random_int32(std::int32_t p, std::int32_t r) {
  return p + rand() % (r - p + 1);
}

std::int64_t random_int64(std::int64_t p, std::int64_t r) {
  std::int64_t x = random_int32(0, 1000000000);
  std::int64_t y = random_int32(0, 1000000000);
  std::int64_t z = x * 1000000000L + y;
  return p + z % (r - p + 1);
}

void fill_random_string(unsigned char* &s, std::uint64_t length, std::uint64_t sigma) {
  for (std::uint64_t i = 0; i < length; ++i)
    s[i] = random_int32(0, sigma - 1);
}

void fill_random_letters(unsigned char* &s, std::uint64_t length, std::uint64_t sigma) {
  fill_random_string(s, length, sigma);
  for (std::uint64_t i = 0; i < length; ++i)
    s[i] += 'a';
}

std::string random_string_hash() {
  uint64_t hash = (uint64_t)rand() * RAND_MAX + rand();
  std::stringstream ss;
  ss << hash;
  return ss.str();
}

std::uint64_t log2ceil(std::uint64_t x) {
  std::uint64_t pow2 = 1, w = 0;
  while (pow2 < x) { pow2 <<= 1; ++w; }
  return w;
}

std::uint64_t log2floor(std::uint64_t x) {
  std::uint64_t pow2 = 1, w = 0;
  while ((pow2 << 1) <= x) { pow2 <<= 1; ++w; }
  return w;
}

}  // namespace utils
}  // namespace psascan_private
