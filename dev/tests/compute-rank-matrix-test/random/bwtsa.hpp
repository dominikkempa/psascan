#ifndef __INMEM_SASCAN_BWTSA_H_INCLUDED
#define __INMEM_SASCAN_BWTSA_H_INCLUDED

#include "uint40.h"

namespace inmem_psascan_private {


template<typename sa_type>
struct bwtsa_t {
  sa_type m_sa;
  unsigned char m_bwt;

  inline operator sa_type() const {
    return m_sa;
  }

  bwtsa_t() {
  }

  bwtsa_t(long x) {
    m_sa = (sa_type)x;
  }

  bwtsa_t(int x) {
    m_sa = (sa_type)x;
  }

  bwtsa_t(uint40 x) {
    m_sa = (sa_type)x;
  }

} __attribute__((packed));


}  // namespace inmem_psascan_private

#endif  // __BWTSA_H_INCLUDED
