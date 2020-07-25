#ifndef __SMALLER_SUFFIXES_H_INCLUDED
#define __SMALLER_SUFFIXES_H_INCLUDED

#include <string>

long parallel_smaller_suffixes(unsigned char *block, long block_size,
    std::string text_filename, long suffix_start_pos);

#endif // __SMALLER_SUFFIXES_H_INCLUDED
