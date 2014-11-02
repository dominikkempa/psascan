#ifndef __HALF_BLOCK_INFO_H_INCLUDED
#define __HALF_BLOCK_INFO_H_INCLUDED

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

#endif  // __HALF_BLOCK_INFO_H_INCLUDED
