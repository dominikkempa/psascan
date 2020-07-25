#ifndef __UTILS_H_INCLUDED
#define __UTILS_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <string>
#include <sstream>


namespace utils {

long double wclock();

std::FILE *file_open(std::string fname, std::string mode);
std::uint64_t file_size(std::string fname);
bool file_exists(std::string fname);
void file_delete(std::string fname);
std::string absolute_path(std::string fname);
void drop_disk_pages(std::string filename);

template<typename value_type>
void write_to_file(const value_type *src, std::uint64_t length, std::FILE *f) {
  std::uint64_t fwrite_ret = std::fwrite(src, sizeof(value_type), length, f);
  if (fwrite_ret != length) {
    fprintf(stderr, "Error: fwrite failed.\n");
    std::exit(EXIT_FAILURE);
  }
}

template<typename value_type>
void write_to_file(const value_type *src, std::uint64_t length, std::string fname) {
  std::FILE *f = file_open(fname, "w");
  write_to_file(src, length, f);
  std::fclose(f);
}

template<typename value_type>
void read_from_file(value_type* dest, std::uint64_t length, std::FILE *f) {
  std::uint64_t fread_ret = std::fread(dest, sizeof(value_type), length, f);
  if (fread_ret != length) {
    fprintf(stderr, "\nError: fread failed.\n");
    std::exit(EXIT_FAILURE);
  }
}

template<typename value_type>
void read_from_file(value_type* dest, std::uint64_t length, std::string fname) {
  std::FILE *f = file_open(fname, "r");
  read_from_file<value_type>(dest, length, f);
  std::fclose(f);
}

template<typename value_type>
void read_at_offset(value_type *dest, std::uint64_t offset,
    std::uint64_t length, std::FILE *f) {
  std::fseek(f, sizeof(value_type) * offset, SEEK_SET);
  read_from_file(dest, length, f);
}

template<typename value_type>
void read_at_offset(value_type *dest, std::uint64_t offset,
    std::uint64_t length, std::string filename) {
  std::FILE *f = file_open(filename, "r");
  read_at_offset(dest, offset, length, f);
  std::fclose(f);
}

std::int32_t random_int32(std::int32_t p, std::int32_t r);
std::int64_t random_int64(std::int64_t p, std::int64_t r);
void fill_random_string(std::uint8_t* &s, std::uint64_t length, std::uint64_t sigma);
void fill_random_letters(std::uint8_t* &s, std::uint64_t length, std::uint64_t sigma);
std::string random_string_hash();

std::uint64_t log2ceil(std::uint64_t x);
std::uint64_t log2floor(std::uint64_t x);

template<typename int_type>
std::string intToStr(int_type x) {
  std::stringstream ss;
  ss << x;
  return ss.str();
}

}  // namespace utils

#endif  // __UTILS_H_INCLUDED
