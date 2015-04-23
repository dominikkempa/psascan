#ifndef __INMEM_SASCAN_PRIVATE_SRANK_AUX_H_INCLUDED
#define __INMEM_SASCAN_PRIVATE_SRANK_AUX_H_INCLUDED


namespace inmem_sascan_private {

//==============================================================================
// Compute ms-decomposition of text[0..length) from ms-decomposition of
// text[0..length - 1). The result is returned via updated values s, p, r.
//==============================================================================
template<typename T>
inline void update_ms(unsigned char *text, T length, T &s, T &p) {
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

}  // namespace inmem_sascan


#endif  // __INMEM_SASCAN_PRIVATE_SRANK_AUX_H_INCLUDED
