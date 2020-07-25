#ifndef __SRANK_AUX_H
#define __SRANK_AUX_H


//==============================================================================
// Compute ms-decomposition of text[0..length) from ms-decomposition of
// text[0..length - 1). The result is returned via updated values s, p, r.
//==============================================================================
template<typename T>
void next(unsigned char *text, T length, T &s, T &p, T &r) {
  if (length == 1) { s = 0; p = 1; r = 0; return; }
  T i = length - 1;
  while (i < length) {
    unsigned char a = text[s + r], b = text[i];
    if (a > b) { p = i - s + 1; r = 0; }
    else if (a < b) { i -= r; s = i; p = 1; r = 0; }
    else { ++r; if (r == p) r = 0; } ++i;
  }
}


#endif  // __SRANK_AUX_H
