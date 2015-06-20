#ifndef __MULTIFILE_H_INCLUDED
#define __MULTIFILE_H_INCLUDED

#include <vector>
#include <string>

#include "utils.h"


struct single_file_info {
  long m_beg;
  long m_end;
  std::string m_filename;

  single_file_info(long beg, long end, std::string filename) {
    m_beg = beg;
    m_end = end;
    m_filename = filename;
  }
};

struct multifile {
  std::vector<single_file_info> files_info;

  void add_file(long beg, long end, std::string filename) {
    files_info.push_back(single_file_info(beg, end, filename));
  }

  ~multifile() {
    for (size_t i = 0; i < files_info.size(); ++i)
      utils::file_delete(files_info[i].m_filename);
  }
};

#endif  // __MULTIFILE_H_INCLUDED
