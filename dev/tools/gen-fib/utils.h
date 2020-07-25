/**
 * @file    utils.h
 * @section LICENCE
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

#ifndef __UTILS_H_INCLUDED
#define __UTILS_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <string>
#include <sstream>


namespace utils {

#define STRX(x) #x
#define STR(x) STRX(x)

// Time
long double wclock();

// Basic file handling
std::FILE *open_file(std::string fname, std::string mode);
std::size_t file_size(std::string fname);
bool file_exists(std::string fname);
void file_delete(std::string fname);
std::string absolute_path(std::string fname);

// File I/O
void read_block(std::string fname, std::size_t beg, std::size_t length, unsigned char *b);
void read_block(std::FILE *f, std::size_t beg, std::size_t length, unsigned char *b);

template<typename value_type>
void write_objects_to_file(const value_type *tab, std::size_t length, std::string fname) {
  std::FILE *f = open_file(fname, "w");
  std::size_t fwrite_ret = std::fwrite(tab, sizeof(value_type), length, f);
  if (fwrite_ret != length) {
    fprintf(stderr, "\nError: fwrite in line %s of %s returned %lu\n",
        STR(__LINE__), STR(__FILE__), fwrite_ret);
    std::exit(EXIT_FAILURE);
  }

  std::fclose(f);
}

template<typename value_type>
void add_objects_to_file(const value_type *tab, std::size_t length, std::FILE *f) {
  std::size_t fwrite_ret = std::fwrite(tab, sizeof(value_type), length, f);
  if (fwrite_ret != length) {
    fprintf(stderr, "\nError: fwrite in line %s of %s returned %lu\n",
        STR(__LINE__), STR(__FILE__), fwrite_ret);
    std::exit(EXIT_FAILURE);
  }
}

template<typename value_type>
void add_objects_to_file(const value_type *tab, std::size_t length, std::string fname) {
  std::FILE *f = utils::open_file(fname.c_str(), "a");
  add_objects_to_file<value_type>(tab, length, f);
  std::fclose(f);
}

template<typename value_type>
void read_n_objects_from_file(value_type* tab, std::size_t length, std::FILE *f) {
  std::size_t fread_ret = std::fread(tab, sizeof(value_type), length, f);
  if (fread_ret != length) {
    fprintf(stderr, "\nError: fread in line %s of %s returned %lu\n",
        STR(__LINE__), STR(__FILE__), fread_ret);
    std::exit(EXIT_FAILURE);
  }
}

template<typename value_type>
void read_n_objects_from_file(value_type* tab, std::size_t length, std::string fname) {
  std::FILE *f = open_file(fname, "r");
  read_n_objects_from_file<value_type>(tab, length, f);
  std::fclose(f);
}

template<typename value_type>
void read_objects_from_file(value_type* &tab, std::size_t &length, std::string fname) {
  std::FILE *f = open_file(fname, "r");
  std::fseek(f, 0L, SEEK_END);
  length = std::ftell(f) / sizeof(value_type);
  std::rewind(f);
  tab = (value_type *)malloc(length * sizeof(value_type));
  read_n_objects_from_file<value_type>(tab, length, f);
  std::fclose(f);
}

// Randomness
int random_int(int p, int r);
long random_long(long p, long r);
void fill_random_string(unsigned char* &s, std::size_t length, int sigma);
void fill_random_letters(unsigned char* &s, std::size_t length, int sigma);
std::string random_string_hash();

// Math
long log2ceil(long x);
long log2floor(long x);

// Misc
template<typename int_type>
std::string intToStr(int_type x) {
  std::stringstream ss;
  ss << x;
  return ss.str();
}

}  // namespace utils

#endif  // __UTILS_H_INCLUDED
