#ifndef __SASCAN_H
#define __SASCAN_H

#include <string>

// Run SAscan with given RAM limit.
void SAscan(std::string     filename,
            long            ram_use,
            unsigned char** BWT = NULL,
            bool            compute_bwt = false
);

#endif // __SASCAN_H
