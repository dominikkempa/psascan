/**
 * @file    src/psascan_src/inmem_psascan_src/srank_aux.hpp
 * @section LICENCE
 *
 * This file is part of pSAscan v0.2.0
 * See: http://www.cs.helsinki.fi/group/pads/
 *
 * Copyright (C) 2014-2017
 *   Juha Karkkainen <juha.karkkainen (at) cs.helsinki.fi>
 *   Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/

#ifndef __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_SRANK_AUX_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_SRANK_AUX_HPP_INCLUDED

#include <cstdint>


namespace psascan_private {
namespace inmem_psascan_private {

//==============================================================================
// Compute ms-decomposition of text[0..length) from ms-decomposition of
// text[0..length - 1). The result is returned via updated values s, p, r.
//==============================================================================
template<typename uint_type>
inline void update_ms(const std::uint8_t *text, uint_type length,
    uint_type &s, uint_type &p) {
  if (length == 1) {
    s = 0;
    p = 1;
    return;
  }

  uint_type i = length - 1;
  while (i < length) {
    std::uint8_t a = text[i - p];
    std::uint8_t b = text[i];

    if (a > b) p = i - s + 1;
    else if (a < b) {
      uint_type r = (i - s);
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

#endif  // __SRC_PSASCAN_SRC_INMEM_PSASCAN_SRC_SRANK_AUX_HPP_INCLUDED
