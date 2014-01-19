#ifndef __FGM_H
#define __FGM_H

#include <string>

// Run FGM with given RAM limit.
void FGM(std::string filename, long ram_use);

// Run FGM with given block size.
void FGM_block_size(std::string filename, long max_block_size);

// Main procedure computing partial suffix arrays.
long partial_sorting(std::string filename, long max_block_size);

#endif // __FGM_H
