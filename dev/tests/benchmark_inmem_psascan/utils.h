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

}  // namespace utils

#endif  // __UTILS_H_INCLUDED
