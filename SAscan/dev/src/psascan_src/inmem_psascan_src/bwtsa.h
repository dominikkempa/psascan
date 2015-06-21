#ifndef __PSASCAN_SRC_INMEM_PSASCAN_SRC_BWTSA_H_INCLUDED
#define __PSASCAN_SRC_INMEM_PSASCAN_SRC_BWTSA_H_INCLUDED

#include "../uint40.h"


namespace psascan_private {
namespace inmem_psascan_private {

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

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __PSASCAN_SRC_INMEM_PSASCAN_SRC_BWTSA_H_INCLUDED
