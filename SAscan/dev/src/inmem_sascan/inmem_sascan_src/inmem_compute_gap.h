//==============================================================================
// This file implements a method:
//
// template<typename T>
// void compute_gap(
//                   unsigned char*    text,
//                   long              text_length,
//                   long              left_block_beg,
//                   long              left_block_size,
//                   long              right_block_size,
//                   T*              partial_sa,
//                   bitvector*        gt_in,
//                   bitvector*&       gt_out,
//                   bool              compute_gt_out,
//                   inmem_gap_array*  gap,
//                   long              max_threads
// );
//
// DESCRIPTION OF INPUT PARAMETERS:
// --------------------------------
//
//  * The algorithm deals with two substrings of text called left and right
//    block. Right block start immediatelly after left block and has size
//    right_block_size. Left block starts at left_block_size and has length
//    left_block_size. Formally:
//    - left block  = text[left_block_beg . .left_block_beg + left_block_size)
//    - right block = text[left_block_end .. left_block_end + right_block_size)
//      where left_block_end = left_block_beg + left_block_size.
//
//  * partial_sa is an array of size left_block_size that contains permutation
//    of integers 0, 1, .., left_block_size - 1. It contains the ordering of
//    the suffixes starting inside the left block and extending until the end
//    of text.
//
//  * gt_in is a bitvector of length right_block_size + 1. We have gt[i] = 1
//    iff text[left_block_end + i..text_length) > text[left_block_end ..
//    text_length). Formally, we never need gt[0], so it could be omitted.
//
//  * gt_out is a bitvector of length left_block_size + right_blocK_size + 1
//    to be computed. It tells whether suffixes at positions starting inside
//    left and right block (and also on one position after the right block)
//    are greater than the suffix starting at left_block_beg. It is
//    indexed from 0, so we know what gt_out[0] will always be 0, since
//    it compares suffix starting at left_block_beg with itself.
//
//  * bool compute_gt_out tells, whether gt_out should be computed or not.
//
//  * gap is the pointer to the in-memory gap array that is the result of
//    this computation. The gap array is of size left_block_size + 1 and is
//    defined as follows:
//    - for i = 1, .., left_block_size - 1  gap[i] = the number of suffixes
//      starting inside the right block (and extending until the end of the
//      string), that are lexicoographically between suffixes starting at
//      positions left_block_beg + partial_sa[i - 1] and left_block_beg +
//      partial_sa[i] of text (and extending until the end of of the text).
//    - gap[0] is the number of suffixes starting inside the right block
//      that are lexicographically smaller than the suffix of text starting
//      at position left_block_beg + partial_sa[0].
//    - gap[left_block_size] is the number of suffixes starting inside
//      the right block that are lexocpgraphically larger than the suffix
//      of text starting as position left_block_size + partial_sa[
//      left_block_size - 1]
//
//  * max_threads is the desired number of threads used for computation. It
//    should be equal to the number of physical cores installed in a machine.
//    For CPUs with Hyper Threading, the number should be 2 * #cores.
//==============================================================================

#ifndef __INMEM_SASCAN_INMEM_COMPUTE_GAP_H_INCLUDED
#define __INMEM_SASCAN_INMEM_COMPUTE_GAP_H_INCLUDED

#include "../../bitvector.h"
#include "rank.h"
#include "buffer.h"
#include "inmem_gap_array.h"
#include "inmem_smaller_suffixes.h"
#include "inmem_stream.h"
#include "inmem_update.h"
#include "inmem_bwt_from_sa.h"
#include "pagearray.h"
#include "bwtsa.h"
#include "../../multifile.h"

namespace inmem_sascan_private {

template<typename saidx_t, unsigned pagesize_log>
void inmem_compute_gap(unsigned char *text, long text_length, long left_block_beg,
    long left_block_size, long right_block_size,
    const pagearray<bwtsa_t<saidx_t>, pagesize_log> &bwtsa,
    bitvector *gt, inmem_gap_array* &gap, long max_threads, bool need_gt, long i0,
    long stream_buffer_size,
    long double &rank_init_time, long double &streaming_time,
    long text_beg,
    long text_end,
    long supertext_length,
    std::string supertext_filename,
    multifile *tail_gt_begin_reversed) {


  //----------------------------------------------------------------------------
  // STEP 1: build rank data structure over BWT.
  //----------------------------------------------------------------------------
  fprintf(stderr, "    Building rank: ");
  long double start = utils::wclock();
  typedef rank4n<saidx_t, pagesize_log> rank_type;
  rank_type *rank = new rank_type(&bwtsa, left_block_size, max_threads);
  rank_init_time = utils::wclock() - start;
  fprintf(stderr, "total: %.2Lf\n", rank_init_time);


  //----------------------------------------------------------------------------
  // STEP 2: compute symbol counts and the last symbol of the left block.
  //----------------------------------------------------------------------------
  long *count = new long[256];
  unsigned char *left_block = text + left_block_beg;
  std::copy(rank->m_count, rank->m_count + 256, count);
  unsigned char last = left_block[left_block_size - 1];
  ++count[last];
  --count[0];
  for (long i = 0, s = 0, t; i < 256; ++i)
    { t = count[i]; count[i] = s; s += t; }



  //----------------------------------------------------------------------------
  // STEP 3: compute starting positions for all streaming threads.
  //----------------------------------------------------------------------------
  long left_block_end = left_block_beg + left_block_size;
  long right_block_beg = left_block_end;
  long right_block_end = left_block_end + right_block_size;

  long max_stream_block_size = (right_block_size + max_threads - 1) / max_threads;
  while (max_stream_block_size & 7) ++max_stream_block_size;
  long n_threads = (right_block_size + max_stream_block_size - 1) / max_stream_block_size;

  fprintf(stderr, "    Computing initial ranks: ");
  start = utils::wclock();
  std::vector<long> initial_ranks(n_threads);
  std::thread **threads = new std::thread*[n_threads];
  for (long i = 0; i < n_threads; ++i) {
    long stream_block_beg = right_block_beg + i * max_stream_block_size;
    long stream_block_end = std::min(stream_block_beg + max_stream_block_size, right_block_end);

    // The i-th thread streams symbols text[beg..end), right-to-left.
    // where beg = stream_block_beg[i], end = stream_block_end[i];
    typedef pagearray<bwtsa_t<saidx_t>, pagesize_log> pagearray_bwtsa_type;
    threads[i] = new std::thread(inmem_smaller_suffixes<pagearray_bwtsa_type>, text,
        text_length, left_block_beg, left_block_end, stream_block_end,
        std::ref(bwtsa), std::ref(initial_ranks[i]), text_beg, text_end, supertext_length,
        supertext_filename, tail_gt_begin_reversed);
  }
  for (long i = 0; i < n_threads; ++i) threads[i]->join();
  for (long i = 0; i < n_threads; ++i) delete threads[i];
  delete[] threads;
  fprintf(stderr, "%.2Lf\n", utils::wclock() - start);





  //----------------------------------------------------------------------------
  // STEP 4: allocate gap array. As explained in the description of the module
  //         on the top of the page, the gap array is indexed from 0 to
  //         left_block_size so the number of elements is left_block_size + 1.
  //----------------------------------------------------------------------------
  start = utils::wclock();
  gap = new inmem_gap_array(left_block_size + 1);


  //----------------------------------------------------------------------------
  // STEP 5: allocate buffers, buffer polls and auxiliary arrays.
  //----------------------------------------------------------------------------

  // Allocate buffers.
  long n_stream_buffers = 2 * n_threads;
  buffer<saidx_t> **buffers = new buffer<saidx_t>*[n_stream_buffers];
  for (long i = 0; i < n_stream_buffers; ++i)
    buffers[i] = new buffer<saidx_t>(stream_buffer_size, max_threads);

  // Create poll of empty and full buffers.
  buffer_poll<saidx_t> *empty_buffers = new buffer_poll<saidx_t>();
  buffer_poll<saidx_t> *full_buffers = new buffer_poll<saidx_t>(n_threads);

  // Add empty buffers to empty poll.
  for (long i = 0; i < n_stream_buffers; ++i)
    empty_buffers->add(buffers[i]);

  // Allocate temp arrays and oracles.
  long max_buffer_elems = stream_buffer_size / sizeof(saidx_t);
  saidx_t *temp = (saidx_t *)malloc(max_buffer_elems * n_threads * sizeof(saidx_t));
  int *oracle = (int *)malloc(max_buffer_elems * n_threads * sizeof(int));
  long double allocations_time = utils::wclock() - start;
  if (allocations_time > 0.05L)
    fprintf(stderr, "    Allocations: %.2Lf\n", allocations_time);


  //----------------------------------------------------------------------------
  // STEP 6: stream.
  //----------------------------------------------------------------------------

  // Start streaming threads.
  fprintf(stderr, "    Streaming: ");
  start = utils::wclock();
  threads = new std::thread*[n_threads];
  for (long t = 0; t < n_threads; ++t) {
    long beg = right_block_beg + t * max_stream_block_size;
    long end = std::min(beg + max_stream_block_size, right_block_end);

    threads[t] = new std::thread(inmem_parallel_stream<rank_type, saidx_t>,
      text, beg, end, last, count, full_buffers, empty_buffers,
      initial_ranks[t], i0, rank, gap->m_length, max_threads, gt,
      temp + t * max_buffer_elems, oracle + t * max_buffer_elems, need_gt);
  }

  // Start updating thread.
  std::thread *updater = new std::thread(inmem_gap_updater<saidx_t>,
      full_buffers, empty_buffers, gap, max_threads);

  // Wait to all threads to finish.
  for (long t = 0; t < n_threads; ++t) threads[t]->join();
  updater->join();
  streaming_time = utils::wclock() - start;
  long double streaming_speed =
    (right_block_size / (1024.L * 1024)) / streaming_time;
  fprintf(stderr, "%.2Lf (%.2LfMiB/s)\n", streaming_time,
      streaming_speed);

  //----------------------------------------------------------------------------
  // Clean up and sort m_excess. Consider using gnu parallel sort.
  //----------------------------------------------------------------------------
  start = utils::wclock();
  free(oracle);
  free(temp);
  for (long i = 0; i < n_threads; ++i) delete threads[i];
  for (long i = 0; i < n_stream_buffers; ++i) delete buffers[i];
  delete updater;
  delete[] threads;
  delete[] buffers;
  delete empty_buffers;
  delete full_buffers;
  delete rank;
  delete[] count;

  std::sort(gap->m_excess.begin(), gap->m_excess.end());

  long double cleaning_time = utils::wclock() - start;
  if (cleaning_time > 0.1L)
    fprintf(stderr, "    Cleaning: %.2Lf\n", cleaning_time);
}

}  // namespace inmem_sascan
                 
#endif  // __INMEM_SASCAN_INMEM_COMPUTE_GAP_H_INCLUDED
