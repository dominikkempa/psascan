// Common routines. (c) Dominik Kempa 2013.
#ifndef __UTILS_H
#define __UTILS_H

#include <cstdio>
#include <ctime>
#include <string>
#include <sstream>

namespace utils {

#define STRX(x) #x
#define STR(x) STRX(x)

/******************************* SYSTEM CALLS *********************************/

void execute(std::string cmd);
std::string get_absolute_dir(std::string fname);
std::string get_absolute_path(std::string fname);

/****************************** MEASURING TIME ********************************/

long double wclock();
long double elapsed(clock_t &ts);

/**************************** FILE MANIPULATION *******************************/

// Basic routines.
FILE *open_file(std::string fname, std::string mode);
long file_size(std::string fname);
bool file_exists(std::string fname);
void file_delete(std::string fname);

// Writing sequences.
void write_text_to_file(unsigned char *text, long length, std::string fname);
void write_ints_to_file(int *tab, long length, std::string fname);
void add_ints_to_file(int *tab, long length, FILE *f);

template<typename value_type>
void write_objects_to_file(value_type *tab, long length, std::string fname) {
  FILE *f = open_file(fname, "w");
  long fwrite_ret = fwrite(tab, sizeof(value_type), length, f);
  if (fwrite_ret != length) {
    fprintf(stderr, "Error: fwrite in line %s of %s returned %ld\n",
        STR(__LINE__), STR(__FILE__), fwrite_ret);
    std::exit(EXIT_FAILURE);
  }
  fclose(f);
}

template<typename value_type>
void add_objects_to_file(value_type *tab, long length, FILE *f) {
  long fwrite_ret = fwrite(tab, sizeof(value_type), length, f);
  if (fwrite_ret != length) {
    fprintf(stderr, "Error: fwrite in line %s of %s returned %ld\n",
        STR(__LINE__), STR(__FILE__), fwrite_ret);
    std::exit(EXIT_FAILURE);
  }
}

// Reading sequences.
void read_text_from_file(unsigned char* &text, long length, std::string fname);
void read_ints_from_file(int* &tab, long length, std::string fname);
void read_file(unsigned char* &text, long &length, std::string fname);
void read_block(std::string fname, long beg, long length, unsigned char *b);

template<typename value_type>
void read_objects_from_file(value_type* &tab, long &length, std::string fname) {
  FILE *f = open_file(fname, "r");
  fseek(f, 0, SEEK_END);
  length = (long)(ftell(f) / sizeof(value_type)); // No length was given.
  rewind(f);
  tab = new value_type[length];
  if (!tab) {
    fprintf(stderr, "Error: read_objects_from_file: alloc failed\n");
    std::exit(EXIT_FAILURE);
  }
  long fread_ret = (long)fread(tab, sizeof(value_type), length, f);
  if (fread_ret != length) {
    fprintf(stderr, "Error: fread in line %s of %s returned %ld\n",
        STR(__LINE__), STR(__FILE__), fread_ret);
    std::exit(EXIT_FAILURE);
  }
  fclose(f);
}

template<typename value_type>
void read_objects_from_file(value_type* &tab, long length, std::FILE *f) {
  size_t fread_ret = fread(tab, sizeof(value_type), length, f);
  if (fread_ret != (size_t)length) {
    fprintf(stderr, "Error: fread in line %s of %s returned %ld\n",
        STR(__LINE__), STR(__FILE__), fread_ret);
    std::exit(EXIT_FAILURE);
  }
}

template<typename value_type>
void read_n_objects_from_file(value_type* &tab, long length, std::string fname) {
  tab = new value_type[length + 5];
  if (!tab) {
    fprintf(stderr, "Error: alloc faild in line %s of %s\n",
        STR(__LINE__), STR(__FILE__));
    std::exit(EXIT_FAILURE);
  }
  FILE *f = open_file(fname, "r");
  long fread_ret = fread(tab, sizeof(value_type), length, f);
  if (fread_ret != length) {
    fprintf(stderr, "Error: fread in line %s of %s returned %ld\n",
        STR(__LINE__), STR(__FILE__), fread_ret);
    std::exit(EXIT_FAILURE);
  }
  fclose(f);
}

// Reading single objects from fie.
long double read_ld_from_file(std::string fname);
int read_int_from_file(std::string fname);

/******************************* RANDOMNESS ***********************************/

int random_int(int p, int r);
long random_long(long p, long r);
void fill_random_string(unsigned char* &s, long length, int sigma);
void fill_random_letters(unsigned char* &s, long n, int sigma);
std::string random_string_hash();

/********************************* MATH ***************************************/

long log2ceil(long x);

/********************************* MISC ***************************************/

template<typename int_type>
std::string intToStr(int_type x) {
  std::stringstream ss;
  ss << x;
  return ss.str();
}

void find_stxxl_config();

} // namespace utils

#endif // __UTILS_H
