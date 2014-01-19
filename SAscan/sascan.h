#ifndef __SASCAN_H
#define __SASCAN_H

#include <string>

// Run SAscan with given RAM limit.
void SAscan(std::string filename, long ram_use);

// Run SAscan with given block size.
void SAscan_block_size(std::string filename, long max_block_size);

#endif // __SASCAN_H
