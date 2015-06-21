#ifndef __PSASCAN_SRC_INMEM_PSASCAN_SRC_SRANK_AUX_H_INCLUDED
#define __PSASCAN_SRC_INMEM_PSASCAN_SRC_SRANK_AUX_H_INCLUDED


namespace psascan_private {
namespace inmem_psascan_private {

//==============================================================================
// Compute ms-decomposition of text[0..length) from ms-decomposition of
// text[0..length - 1). The result is returned via updated values s, p, r.
//==============================================================================
template<typename T>
inline void update_ms(const unsigned char *text, T length, T &s, T &p) {
  if (length == 1) { s = 0; p = 1; return; }

  T i = length - 1;
  while (i < length) {
    unsigned char a = text[i - p];
    unsigned char b = text[i];

    if (a > b) p = i - s + 1;
    else if (a < b) {
      long r = (i - s);
      while (r >= p) r -= p;
      i -= r;
      s = i;
      p = 1;
    }

    ++i;
  }
}

}  // namespace inmem_psascan_private
}  // namespace psascan_private

#endif  // __PSASCAN_SRC_INMEM_PSASCAN_SRC_SRANK_AUX_H_INCLUDED
