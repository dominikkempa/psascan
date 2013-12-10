#ifndef __FILE_STREAMER_H
#define __FILE_STREAMER_H

#include <string>

#include "utils.h"

/********************************* usage ***************************************
file_streamer<int> f("in.sa", 1 << 21);
while (!f.empty()) {
  int next = f.read();
  printf("%d ", next);
}
*******************************************************************************/

template<typename data_type>
struct file_streamer {
  file_streamer(std::string fname, int bufsize)
      : m_bufelems(bufsize / sizeof(data_type)) {
    f = utils::open_file(fname, "r");
    buffer = new data_type[m_bufelems];
    if (!buffer) {
      fprintf(stderr, "Error: not enough memory for file_streamer buffer!\n");
      std::exit(EXIT_FAILURE);
    }
    refill();
  }

  ~file_streamer() {
    delete[] buffer;
    std::fclose(f);
  }

  inline data_type read() {
    data_type ret = buffer[pos++];
    if (pos == filled) refill();
    
    return ret;
  }

  inline bool empty() {
    return (!filled && !refill());
  }
  
private:
  int m_bufelems, filled, pos;
  data_type *buffer;
  std::FILE *f;

  int refill() {
    filled = std::fread(buffer, sizeof(data_type), m_bufelems, f);
    pos = 0;
    return filled;
  }
};

#endif // __FILE_STREAMER_H

