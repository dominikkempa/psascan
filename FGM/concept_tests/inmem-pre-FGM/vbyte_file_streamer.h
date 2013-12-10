#ifndef __VBYTE_FILE_STREAMER_H
#define __VBYTE_FILE_STREAMER_H

#include <string>

#include "utils.h"

/********************************* usage ***************************************
vbyte_file_streamer<int> f("in.sa", 1 << 21);
while (!f.empty()) {
  int next = f.read();
  printf("%d ", next);
}
*******************************************************************************/

struct vbyte_file_streamer {
  vbyte_file_streamer(std::string fname, int bufsize)
      : m_bufsize(bufsize) {
    f = utils::open_file(fname, "r");
    buf = new unsigned char[m_bufsize];
    refill();
  }

  ~vbyte_file_streamer() {
    delete[] buf;
    std::fclose(f);
  }

  inline int read() {
    int ret = 0, offset = 0;
    while (buf[pos] & 0x80) {
      ret |= ((buf[pos++] & 0x7f) << offset);
      if (pos == filled) refill();
      offset += 7;
    }
    ret |= (buf[pos++] << offset);
    if (pos == filled) refill();

    return ret;
  }
  
private:
  int m_bufsize, filled, pos;
  unsigned char *buf;
  std::FILE *f;

  inline void refill() {
    filled = fread(buf, 1, m_bufsize, f);
    pos = 0;
  }
};

#endif // __VBYTE_FILE_STREAMER_H

