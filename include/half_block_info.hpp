/**
 * @file    src/psascan_src/half_block_info.hpp
 * @section LICENCE
 *
 * This file is part of pSAscan v0.1.1
 * See: https://github.com/dominikkempa/psascan
 *
 * Copyright (C) 2014-2020
 *   Dominik Kempa <dominik.kempa (at) gmail.com>
 *   Juha Karkkainen <juha.karkkainen (at) cs.helsinki.fi>
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

#ifndef __SRC_PSASCAN_SRC_HALF_BLOCK_INFO_HPP_INCLUDED
#define __SRC_PSASCAN_SRC_HALF_BLOCK_INFO_HPP_INCLUDED

#include <string>

#include "io/distributed_file.hpp"


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

#endif  // __SRC_PSASCAN_SRC_HALF_BLOCK_INFO_HPP_INCLUDED
