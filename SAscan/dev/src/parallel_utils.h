#ifndef __PARALLEL_UTILS_H_INCLUDED
#define __PARALLEL_UTILS_H_INCLUDED

#include <thread>
#include <algorithm>


namespace parallel_utils {

//==============================================================================
// Encode tab[0..length) using vbyte encoding and write to dest sequentially.
//==============================================================================
void encode_vbyte_slab(long *tab, long length, unsigned char *dest) {
  long ptr = 0L;
  for (long j = 0; j < length; ++j) {
    long x = tab[j];
    while (x > 127) {
      dest[ptr++] = ((x & 0x7f) | 0x80);
      x >>= 7;
    }
    dest[ptr++] = x;
  }
}


//==============================================================================
// Compute the size of vbyte encoding of tab[0..length).
//==============================================================================
void compute_size_of_vbyte_slab(long *tab, long length, long &result) {
  result = 0L;
  for (long j = 0; j < length; ++j) {
    long x = tab[j];
    while (x > 127) {
      ++result;
      x >>= 7;
    }
    ++result;
  }
}


//==============================================================================
// Encode tab[0..length) using v-byte encoding and write to dest in parallel.
// We assume that dest is sufficiently large to hold the output.
// The function returns the length of the slab.
//==============================================================================
long convert_array_to_vbyte_slab(long *tab, long length, unsigned char *dest, long max_threads) {
  long max_block_size = (length + max_threads - 1) / max_threads;
  long n_blocks = (length + max_block_size - 1) / max_block_size;

  // 1
  //
  // Compute the length of slab for each block.
  long *block_slab_length = new long[n_blocks];

  std::thread **threads = new std::thread*[n_blocks];
  for (long t = 0; t < n_blocks; ++t) {
    long block_beg = t * max_block_size;
    long block_end = std::min(block_beg + max_block_size, length);
    long block_size = block_end - block_beg;

    threads[t] = new std::thread(compute_size_of_vbyte_slab,
        tab + block_beg, block_size, std::ref(block_slab_length[t]));
  }

  for (long t = 0; t < n_blocks; ++t) threads[t]->join();
  for (long t = 0; t < n_blocks; ++t) delete threads[t];

  // 2
  //
  // Compute cummulative sum for block slab lengths.
  long total_slab_length = 0L;
  for (long j = 0; j < n_blocks; ++j) {
    long temp = block_slab_length[j];
    block_slab_length[j] = total_slab_length;
    total_slab_length += temp;
  }

  // 3
  //
  // Compute the slabs. Now we know where each slab begins.
  for (long t = 0; t < n_blocks; ++t) {
    long block_beg = t * max_block_size;
    long block_end = std::min(block_beg + max_block_size, length);
    long block_size = block_end - block_beg;

    threads[t] = new std::thread(encode_vbyte_slab,
        tab + block_beg, block_size, dest + block_slab_length[t]);
  }

  for (long t = 0; t < n_blocks; ++t) threads[t]->join();
  for (long t = 0; t < n_blocks; ++t) delete threads[t];
  delete[] threads;
  delete[] block_slab_length;

  return total_slab_length;
}


//==============================================================================
// Copy src[0..length) to dest[0..length).
//==============================================================================
template<typename T, typename S>
void parallel_copy_aux(T *src, S *dest, long length) {
  for (long i = 0; i < length; ++i)
    dest[i] = (S)src[i];
}


//==============================================================================
// Parallel version of std::copy (with slightly different interface).
// Conversion from T to S has to make sense.
//==============================================================================
template<typename T, typename S>
void parallel_copy(T *src, S *dest, long length, long max_threads) {
  long max_block_size = (length + max_threads - 1) / max_threads;
  long n_blocks = (length + max_block_size - 1) / max_block_size;

  std::thread **threads = new std::thread*[n_blocks];
  for (long i = 0; i < n_blocks; ++i) {
    long block_beg = i * max_block_size;
    long block_end = std::min(block_beg + max_block_size, length);
    long block_size = block_end - block_beg;

    threads[i] = new std::thread(parallel_copy_aux<T, S>,
        src + block_beg, dest + block_beg, block_size);
  }

  for (long i = 0; i < n_blocks; ++i) threads[i]->join();
  for (long i = 0; i < n_blocks; ++i) delete threads[i];
  delete[] threads;
}


//==============================================================================
// Set all values in tab[0..length) to x.
//==============================================================================
template<typename T>
void parallel_fill_aux(T *tab, long length, T x) {
  for (long i = 0; i < length; ++i)
    tab[i] = x;
}


//==============================================================================
// Parallel version of std::fill (with slightly different interface).
//==============================================================================
template<typename T>
void parallel_fill(T *tab, long length, T x, long max_threads) {
  long max_block_size = (length + max_threads - 1) / max_threads;
  long n_blocks = (length + max_block_size - 1) / max_block_size;

  std::thread **threads = new std::thread*[n_blocks];
  for (long i = 0; i < n_blocks; ++i) {
    long block_beg = i * max_block_size;
    long block_end = std::min(block_beg + max_block_size, length);
    long block_size = block_end - block_beg;

    threads[i] = new std::thread(parallel_fill_aux<T>,
        tab + block_beg, block_size, x);
  }

  for (long i = 0; i < n_blocks; ++i) threads[i]->join();
  for (long i = 0; i < n_blocks; ++i) delete threads[i];
  delete[] threads;
}

}  // namespace parallel_utils

#endif  // __PARALLEL_UTILS_H_INCLUDED
