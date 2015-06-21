#ifndef __PSASCAN_SRC_HALF_BLOCK_INFO_H_INCLUDED
#define __PSASCAN_SRC_HALF_BLOCK_INFO_H_INCLUDED

#include <string>

#include "distributed_file.h"


namespace psascan_private {

// Stores the information about half-blocks.
template<typename block_offset_type>
struct half_block_info {
  long beg;
  long end;

  std::string gap_filename;
  distributed_file<block_offset_type> *psa;

  bool operator < (const half_block_info &i) const {
    return beg < i.beg;
  }
};

}  // namespace psascan_private

#endif  // __PSASCAN_SRC_HALF_BLOCK_INFO_H_INCLUDED
