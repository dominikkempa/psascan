//==============================================================================
// This file implements a method:
//
// void compute_gap(
//                   unsigned char*    text,
//                   long              text_length,
//                   long              left_block_beg,
//                   long              left_block_size,
//                   long              right_block_size,
//                   int*              partial_sa,
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

#ifndef __INMEM_COMPUTE_GAP_H_INCLUDED
#define __INMEM_COMPUTE_GAP_H_INCLUDED

#include "bitvector.h"
#include "rank.h"
#include "buffer.h"
#include "inmem_gap_array.h"
#include "inmem_smaller_suffixes.h"
#include "inmem_stream.h"
#include "inmem_update.h"
#include "inmem_bwt_from_sa.h"
#include "inmem_finalize_gt.h"


void inmem_compute_gap(unsigned char *text, long text_length, long left_block_beg,
    long left_block_size, long right_block_size, int *partial_sa, bitvector *gt_in,
    bitvector* &gt_out, bool compute_gt_out, inmem_gap_array* &gap, long max_threads,
    long stream_buffer_size = (1L << 20)) {

  //----------------------------------------------------------------------------
  // STEP 1: compute BWT from partial_sa. Together with BWT we get the index
  //         i0, such that partial_sa[i0] = 0.
  //----------------------------------------------------------------------------
  unsigned char *left_block = text + left_block_beg;
  unsigned char *bwt = new unsigned char[left_block_size - 1];
  long i0 = bwt_from_sa_into_dest(partial_sa, left_block,
      left_block_size, bwt, max_threads);


  //----------------------------------------------------------------------------
  // STEP 2: build rank data structure over BWT and delete BWT.
  //----------------------------------------------------------------------------
  rank4n<> *rank = new rank4n<>(bwt, left_block_size - 1, max_threads);
  delete[] bwt;


  //----------------------------------------------------------------------------
  // STEP 3: compute symbol counts and the last symbol of the left block.
  //----------------------------------------------------------------------------
  long *count = new long[256];
  std::copy(rank->c_rank, rank->c_rank + 256, count);
  unsigned char last = left_block[left_block_size - 1];
  ++count[last];
  for (long i = 0, s = 0, t; i < 256; ++i)
    { t = count[i]; count[i] = s; s += t; }



  //----------------------------------------------------------------------------
  // STEP 4: compute starting positions for all streaming threads.
  //----------------------------------------------------------------------------
  long left_block_end = left_block_beg + left_block_size;
  long right_block_end = left_block_end + right_block_size;

  //----------------------------------------------------------------------------
  // We need to guarantee that no two threads will attempt to write to the same
  // byte. To do this we separare min(right_block_size, 8 - left_block_size % 8)
  // bits from the right block and designate them to be handled by the first
  // thread.
  //
  // In addition, we make sure that the stream_block_size is a multiple of 8
  // this way, we know all other thread (second, third, etc.) will always start
  // filling bits at position that is a multiple of 8.
  //
  // Making sure that left_block_size is a multiple of 8 and making stream block
  // size to be multiple of 8 would also work, but I feel this is better,
  // as we don't have to make anything else elswhere, in particular making
  // sure that left block size is a multiple of 8 would be a huge pain.
  //----------------------------------------------------------------------------

  // This many bits of the right block will in
  // addition be handled by the first thread.
  long initial_chunk = std::min(right_block_size, 8 - (left_block_size & 7));
  long n_threads = 0L;
  long *stream_block_beg, *stream_block_end;
  if (initial_chunk == right_block_size) {
    // Just one block containing initial chunk.
    n_threads = 1;
    stream_block_beg = new long[n_threads];
    stream_block_end = new long[n_threads];
    stream_block_beg[0] = left_block_end;
    stream_block_end[0] = left_block_end + initial_chunk;
  } else {
    // At least one symbol in addition to initial chunk.
    long reduced_right_block_size = right_block_size - initial_chunk;
    long stream_block_size = (reduced_right_block_size + max_threads - 1) / max_threads;
    while (stream_block_size & 7) ++stream_block_size;
    n_threads = (reduced_right_block_size + stream_block_size - 1) / stream_block_size;
    stream_block_beg = new long[n_threads];
    stream_block_end = new long[n_threads];
    stream_block_beg[0] = left_block_end;
    stream_block_end[0] = std::min(right_block_end,
        left_block_end + initial_chunk + stream_block_size);
    for (long j = 1; j < n_threads; ++j) {
      stream_block_beg[j] = stream_block_end[j - 1];
      stream_block_end[j] = std::min(right_block_end,
          stream_block_beg[j] + stream_block_size);
    }
  }

  std::vector<long> initial_ranks(n_threads);
  std::thread **threads = new std::thread*[n_threads];
  for (long i = 0; i < n_threads; ++i) {
    // The i-th thread streams symbols text[beg..end), right-to-left.
    // where beg = stream_block_beg[i], end = stream_block_end[i];
    threads[i] = new std::thread(inmem_smaller_suffixes<int>, text,
        text_length, left_block_beg, left_block_end, stream_block_end[i],
        partial_sa, std::ref(initial_ranks[i]));
  }
  for (long i = 0; i < n_threads; ++i) threads[i]->join();
  for (long i = 0; i < n_threads; ++i) delete threads[i];
  delete[] threads;



  //----------------------------------------------------------------------------
  // STEP 5: allocate gap array. As explained in the description of the module
  //         on the top of the page, the gap array is indexed from 0 to
  //         left_block_size so the number of elements is left_block_size + 1.
  //----------------------------------------------------------------------------
  gap = new inmem_gap_array(left_block_size + 1, max_threads);


  //----------------------------------------------------------------------------
  // STEP 6: allocate buffers and buffer polls.
  //----------------------------------------------------------------------------

  // Allocate buffers.
  long n_stream_buffers = 2 * n_threads;
  buffer<int> **buffers = new buffer<int>*[n_stream_buffers];
  for (long i = 0; i < n_stream_buffers; ++i)
    buffers[i] = new buffer<int>(stream_buffer_size, max_threads);

  // Create poll of empty and full buffers.
  buffer_poll<int> *empty_buffers = new buffer_poll<int>();
  buffer_poll<int> *full_buffers = new buffer_poll<int>(n_threads);

  // Add empty buffers to empty poll.
  for (long i = 0; i < n_stream_buffers; ++i)
    empty_buffers->add(buffers[i]);


  //----------------------------------------------------------------------------
  // STEP 7: allocate gt_out bitvector.
  //----------------------------------------------------------------------------
  if (!compute_gt_out) gt_out = NULL;
  else {
    gt_out = new bitvector(left_block_size + right_block_size + 1, max_threads);

    // We manually set the last bit.
    if (initial_ranks[n_threads - 1] > i0)
      gt_out->set(left_block_size + right_block_size);

    // We compute the left half of gt_out.
    finalize_gt(text, text_length, left_block_beg, left_block_size,
        gt_in, gt_out, max_threads);
  }


  //----------------------------------------------------------------------------
  // STEP 8: stream.
  //----------------------------------------------------------------------------

  // Start streaming threads.
  threads = new std::thread*[n_threads];
  for (long t = 0; t < n_threads; ++t) {
    long beg = stream_block_beg[t];
    long end = stream_block_end[t];

    threads[t] = new std::thread(inmem_parallel_stream<int>,
      text, beg, end, last, count, full_buffers, empty_buffers,
      initial_ranks[t], i0, rank, gap->m_length, stream_buffer_size,
      max_threads, gt_in, gt_out, compute_gt_out, left_block_beg,
      left_block_end);
  }

  // Start updating thread.
  std::thread *updater = new std::thread(inmem_gap_updater<int>,
      full_buffers, empty_buffers, gap, max_threads);

  // Wait to all threads to finish.
  for (long t = 0; t < n_threads; ++t) threads[t]->join();
  updater->join();
  
  // Clean up.
  for (long i = 0; i < n_threads; ++i) delete threads[i];
  for (long i = 0; i < n_stream_buffers; ++i) delete buffers[i];
  delete updater;
  delete[] threads;
  delete[] buffers;
  delete empty_buffers;
  delete full_buffers;
  delete rank;
  delete[] count;
  delete[] stream_block_beg;
  delete[] stream_block_end;


  //----------------------------------------------------------------------------
  // STEP 9: sort excess values. Consider using gnu parallel sort here.
  //----------------------------------------------------------------------------
  std::sort(gap->m_excess.begin(), gap->m_excess.end());
}
                 
#endif  // __INMEM_COMPUTE_GAP_H_INCLUDED
