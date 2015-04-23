#ifndef __MERGE_H_INCLUDED
#define __MERGE_H_INCLUDED

#include <string>
#include <vector>
#include <algorithm>

#include "utils.h"
#include "io_streamer.h"
#include "uint40.h"
#include "distributed_file.h"
#include "half_block_info.h"

#include "async_stream_reader.h"
#include "async_stream_writer.h"
#include "async_vbyte_stream_reader.h"


// Merge partial suffix arrays into final suffix array.
// INVARIANT: 5.2 * length <= ram_use.
template<typename block_offset_type>
void merge(std::string output_filename, long ram_use, std::vector<half_block_info<block_offset_type> > &hblock_info) {
  long n_block = (long)hblock_info.size();
  long length = 0;

  std::sort(hblock_info.begin(), hblock_info.end());
  for (size_t j = 0; j < hblock_info.size(); ++j)
    length += hblock_info[j].end - hblock_info[j].beg;

  long pieces = (1 + sizeof(block_offset_type)) * n_block - 1 + sizeof(uint40);
  long buffer_size = (ram_use + pieces - 1) / pieces;

  fprintf(stderr, "\nBuffer size for merging: %ld\n", buffer_size);
  fprintf(stderr, "sizeof(output_type) = %ld\n", sizeof(uint40));

  typedef async_vbyte_stream_reader<long> vbyte_reader_type;
  typedef async_stream_writer<uint40> output_writer_type;

  output_writer_type *output = new output_writer_type(output_filename, sizeof(uint40) * buffer_size);
  vbyte_reader_type **gap = new vbyte_reader_type*[n_block - 1];
  for (long i = 0; i < n_block; ++i) {
    hblock_info[i].psa->initialize_reading(sizeof(block_offset_type) * buffer_size);
    if (i + 1 != n_block)
      gap[i] = new vbyte_reader_type(hblock_info[i].gap_filename, buffer_size);
  }

  long *gap_head = new long[n_block];
  for (long i = 0; i + 1 < n_block; ++i)
    gap_head[i] = gap[i]->read();
  gap_head[n_block - 1] = 0;

  fprintf(stderr, "Merge:\r");
  long double merge_start = utils::wclock();
  for (long i = 0, dbg = 0; i < length; ++i, ++dbg) {
    if (dbg == (1 << 23)) {
      long double elapsed = utils::wclock() - merge_start;
      long double scanned_m = i / (1024.L * 1024);
      long inp_vol = (1L + sizeof(uint40)) * i;
      long out_vol = sizeof(uint40) * i;
      long tot_vol = inp_vol + out_vol;
      long double tot_vol_m = tot_vol / (1024.L * 1024);
      long double io_speed = tot_vol_m / elapsed;
      fprintf(stderr, "Merge: %.1Lf%%, time = %.2Lfs (%.3Lfs/MiB), io = %2.LfMiB/s\r",
          (100.L * i) / length, elapsed, elapsed / scanned_m, io_speed);
      dbg = 0;
    }

    // XXX can the method for finding k be too slow with many blocks?
    // What are other, practically faster, options.
    long k = 0;
    while (gap_head[k]) --gap_head[k++];
    if (k != n_block - 1) gap_head[k] = gap[k]->read();
    long SA_i = hblock_info[k].psa->read() + hblock_info[k].beg;
    output->write(SA_i);
  }
  long double merge_time = utils::wclock() - merge_start;
  fprintf(stderr, "Merge: 100.0%%. Time: %.2Lfs\n", merge_time);

  // Clean up.
  delete output;
  for (long i = 0; i < n_block; ++i) {
    hblock_info[i].psa->finish_reading();
    delete hblock_info[i].psa;
    if (i + 1 != n_block)
      delete gap[i];
  }

  delete[] gap;
  delete[] gap_head;
  
  for (int i = 0; i + 1 < n_block; ++i)
    utils::file_delete(hblock_info[i].gap_filename);
}


#endif  // __MERGE_H_INCLUDED

