#ifndef __BWTSA_H_INCLUDED
#define __BWTSA_H_INCLUDED


template<typename sa_type>
struct bwtsa {
  sa_type sa;
  unsigned char bwt;
} __attribute__((packed));

#endif  // __BWTSA_H_INCLUDED
