#ifndef __MERGE_H_INCLUDED
#define __MERGE_H_INCLUDED

#include <string>
#include <algorithm>

#include "utils.h"
#include "stream.h"
#include "uint40.h"
#include "distributed_file.h"

// Merge partial suffix arrays into final suffix array (stored in normal file).
// INVARIANT: 5 * length <= ram_use.
template<typename offset_type>
void merge(std::string input_filename, long length, long max_block_size, long ram_use) {
  std::string out_filename = input_filename + ".sa5";

  long n_block = (length + max_block_size - 1) / max_block_size;
  long pieces = (1 + sizeof(offset_type)) * n_block - 1 + sizeof(uint40);
  long buffer_size = (ram_use + pieces - 1) / pieces;

  fprintf(stderr, "Buffer size for merging: %ld\n", buffer_size);
  fprintf(stderr, "sizeof(output_type) = %ld\n", sizeof(uint40));

  stream_writer<uint40> *output = new stream_writer<uint40>(out_filename, sizeof(uint40) * buffer_size);
  stream_reader<offset_type> **sparseSA = new stream_reader<offset_type>*[n_block];
  vbyte_stream_reader **gap = new vbyte_stream_reader*[n_block - 1];
  for (long i = 0; i < n_block; ++i) {
    sparseSA[i] = new stream_reader<offset_type>(input_filename + ".partial_sa." + utils::intToStr(i), sizeof(offset_type) * buffer_size);
    if (i + 1 != n_block) gap[i] = new vbyte_stream_reader(input_filename + ".gap." + utils::intToStr(i), buffer_size);
  }

  long *gap_head = new long[n_block];
  for (long i = 0; i + 1 < n_block; ++i)
    gap_head[i] = gap[i]->read();
  gap_head[n_block - 1] = 0;

  fprintf(stderr, "Merging:\r");
  long double merge_start = utils::wclock();
  for (long i = 0, dbg = 0; i < length; ++i, ++dbg) {
    if (dbg == (1 << 23)) {
      long double elapsed = utils::wclock() - merge_start;
      fprintf(stderr, "Merging: %.1Lf%%. Time: %.2Lfs\r", (100.L * i) / length, elapsed);
      dbg = 0;
    }

    long k = 0;
    while (gap_head[k]) --gap_head[k++];
    if (k != n_block - 1) gap_head[k] = gap[k]->read();
    long SA_i = sparseSA[k]->read() + k * max_block_size;
    output->write(SA_i);
  }
  long double merge_time = utils::wclock() - merge_start;
  fprintf(stderr, "Merging: 100.0%%. Time: %.2Lfs\n", merge_time);

  // Clean up.
  delete output;
  for (long i = 0; i < n_block; ++i) {
    delete sparseSA[i];
    if (i + 1 != n_block)
      delete gap[i];
  }

  delete[] sparseSA;
  delete[] gap;
  delete[] gap_head;
  
  for (int i = 0; i < n_block; ++i) {
    utils::file_delete(input_filename + ".partial_sa." + utils::intToStr(i));
    if (i + 1 != n_block)
      utils::file_delete(input_filename + ".gap." + utils::intToStr(i));
  }
}

// Merge partial suffix arrays (inside recursive call) into bigger partial
// suffix array and compute BWT associated with this SA. The resulting SA is
// stored using distrbuted file (which it returned).
// INVARIANT: 5 * length <= ram_use.
//
// XXX: currently distributed file is not used.
template<typename offset_type, typename output_type>
//distributed_file<output_type> *partial_merge(
void partial_merge(
    std::string input_filename,
    long length,
    long max_block_size,
    long ram_use,
    unsigned char **BWT,
    std::string text_filename,
    long text_offset) {
  std::string out_filename = input_filename + ".sa5";
  long n_block = (length + max_block_size - 1) / max_block_size;

  unsigned char *text = NULL;
  long pieces = (1 + sizeof(offset_type)) * n_block - 1 + sizeof(output_type);

  // Read the original block of text
  text = new unsigned char[length];
  utils::read_block(text_filename, text_offset, length, text);

  *BWT = new unsigned char[length - 1];
  long merge_ram_use = ram_use - 2 * length; // > 0
  long buffer_size = (merge_ram_use + pieces - 1) / pieces;

  fprintf(stderr, "Buffer size for merging: %ld\n", buffer_size);
  fprintf(stderr, "sizeof(offset_type) = %ld\n", sizeof(offset_type));
  fprintf(stderr, "sizeof(output_type) = %ld\n", sizeof(output_type));
  
  // TODO: ram_use/10 -> std::max(1<<20, ram_use/10) 
  // distributed_file<output_type> *output = new distributed_file<output_type>(out_filename, ram_use / 10);
  // output->initialize_writing(sizeof(output_Type) * buffer_size);
  stream_writer<uint40> *output = new stream_writer<uint40>(out_filename, sizeof(uint40) * buffer_size);

  stream_reader<offset_type> **sparseSA = new stream_reader<offset_type>*[n_block];
  vbyte_stream_reader **gap = new vbyte_stream_reader*[n_block - 1];
  for (long i = 0; i < n_block; ++i) {
    sparseSA[i] = new stream_reader<offset_type>(input_filename + ".partial_sa." + utils::intToStr(i), sizeof(offset_type) * buffer_size);
    if (i + 1 != n_block) gap[i] = new vbyte_stream_reader(input_filename + ".gap." + utils::intToStr(i), buffer_size);
  }

  long *gap_head = new long[n_block];
  for (long i = 0; i + 1 < n_block; ++i)
    gap_head[i] = gap[i]->read();
  gap_head[n_block - 1] = 0;

  fprintf(stderr, "Merging:\r");
  long double merge_start = utils::wclock();
  for (long i = 0, bwt_ptr = 0, dbg = 0; i < length; ++i, ++dbg) {
    if (dbg == (1 << 23)) {
      long double elapsed = utils::wclock() - merge_start;
      fprintf(stderr, "Merging: %.1Lf%%. Time: %.2Lfs\r",
          (100.L * i) / length, elapsed);
      dbg = 0;
    }
    
    long k = 0;
    while (gap_head[k]) --gap_head[k++];
    if (k != n_block - 1) gap_head[k] = gap[k]->read();
    long SA_i = sparseSA[k]->read() + k * max_block_size;

    output->write(SA_i);
    if (SA_i) (*BWT)[bwt_ptr++] = text[SA_i - 1];
  }
  long double merge_time = utils::wclock() - merge_start;
  fprintf(stderr, "Merging: 100.0%%. Time: %.2Lfs\n", merge_time);

  // Clean up.
  delete[] text;

  for (long i = 0; i < n_block; ++i) {
    delete sparseSA[i];
    if (i + 1 != n_block)
      delete gap[i];
  }

  delete[] sparseSA;
  delete[] gap;
  delete[] gap_head;

  for (int i = 0; i < n_block; ++i) {
    utils::file_delete(input_filename + ".partial_sa." + utils::intToStr(i));
    if (i + 1 != n_block)
      utils::file_delete(input_filename + ".gap." + utils::intToStr(i));
  }

  // output->finish_writing();
  // return output;
}

#endif // __MERGE_H_INCLUDED

