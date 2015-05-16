#include <sys/stat.h>
#include <fcntl.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cerrno>

#include <string>

namespace utils {

int open_file_direct(std::string filename, int mode) {
  int fd = open(filename.c_str(), O_DIRECT, mode);
  if (fd == -1) {
    fprintf(stderr, "\n open() failed with error [%s]\n", strerror(errno));
    std::exit(EXIT_FAILURE);
  }

  return fd;
}

}  // namespace utils

