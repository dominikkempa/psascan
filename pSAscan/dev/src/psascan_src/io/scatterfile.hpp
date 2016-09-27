/**
 * @file    psascan_src/io/scatterfile.hpp
 * @section LICENCE
 *
 * This file is part of pSAscan v0.2.0
 * See: http://www.cs.helsinki.fi/group/pads/
 *
 * Copyright (C) 2014-2016
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

#ifndef __PSASCAN_SRC_IO_SCATTERFILE_HPP_INCLUDED
#define __PSASCAN_SRC_IO_SCATTERFILE_HPP_INCLUDED

#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>


namespace psascan_private {

template<typename value_type>
struct scatterfile {
  scatterfile() {}
  scatterfile(std::uint64_t max_file_size_in_bytes) {
    m_max_items_per_file = std::max(1UL, max_file_size_in_bytes / sizeof(value_type));
    m_items_written = 0;
  }

  std::uint64_t m_items_written;
  std::uint64_t m_max_items_per_file;
  std::vector<std::string> m_filenames;
};

}  // psascan_private

#endif // __PSASCAN_SRC_IO_SCATTERFILE_HPP_INCLUDED
