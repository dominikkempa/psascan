// Simplified version of functions from STXXL library (v.1.4.1).
// See: http://stxxl.sourceforge.net/
#ifndef __ALIGNED_ALLOC_H_INCLUDED
#define __ALIGNED_ALLOC_H_INCLUDED

#include <stdint.h>
#include <cstdlib>
#include <new>

#define DIRECT_ALIGNMENT 4096 

template <typename MustBeInt> struct aligned_alloc_settings { static bool may_use_realloc; };
template <typename MustBeInt> bool aligned_alloc_settings<MustBeInt>::may_use_realloc = true;

template <size_t Alignment>
inline void* aligned_alloc(size_t size) {
  size_t alloc_size = Alignment + sizeof(char *) + size;
  char* buffer = (char *)std::malloc(alloc_size);
  if (buffer == NULL) throw std::bad_alloc();

  char* reserve_buffer = buffer + sizeof(char *);
  char* result = reserve_buffer + Alignment - (((uint64_t)reserve_buffer) % (Alignment));

  // free unused memory behind the data area so access behind the requested size can be recognized
  size_t realloc_size = (result - buffer) + size;
  if (realloc_size < alloc_size && aligned_alloc_settings<int>::may_use_realloc) {
    char* realloced = (char *)std::realloc(buffer, realloc_size);
    if (buffer != realloced) {
      std::free(realloced);
      aligned_alloc_settings<int>::may_use_realloc = false;
      return aligned_alloc<Alignment>(size);
    }
  }

  *(((char **)result) - 1) = buffer;
  return result;
}

inline void aligned_dealloc(void* ptr) {
  if (!ptr) return;
  char* buffer = *(((char **)ptr) - 1);
  std::free(buffer);
}

#endif  // __ALIGNED_ALLOC_H_INCLUDED
