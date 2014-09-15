#ifndef __INMEM_SASCAN_BWTSA_H_INCLUDED
#define __INMEM_SASCAN_BWTSA_H_INCLUDED

#include "../uint40.h"

namespace inmem_sascan_private {


template<typename sa_type>
struct bwtsa_t {
  sa_type sa;
  unsigned char bwt;

  inline operator sa_type() const {
    return sa;
  }

  bwtsa_t() {
  }

  bwtsa_t(long x) {
    sa = (sa_type)x;
  }

  bwtsa_t(int x) {
    sa = (sa_type)x;
  }

  bwtsa_t(uint40 x) {
    sa = (sa_type)x;
  }

} __attribute__((packed));


}  // namespace inmem_sascan

#endif  // __BWTSA_H_INCLUDED
