#ifndef __LINUX_IO_UTILS_H_INCLUDED
#define __LINUX_IO_UTILS_H_INCLUDED

#include <string>

namespace utils {

int open_file_direct(std::string filename, int mode);

}  // namespace utils

#endif  // __LINUX_IO_UTILS_H_INCLUDED
