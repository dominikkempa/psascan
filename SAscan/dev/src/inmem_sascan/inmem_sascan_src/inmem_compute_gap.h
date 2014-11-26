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

#include <map>

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

inline size_t find(const std::vector<long> &v, long x) {
  size_t left = 0;
  size_t right = v.size();

  while (left != right) {
    // Invariant: The answer is in the range [left..right].
    size_t mid = (left + right) / 2;
    if (v[mid] >= x) right = mid;
    else left = mid + 1;
  }

  if (right != v.size() && v[right] != x)
    right = v.size();
  return right;
}

template<typename saidx_t, unsigned pagesize_log, unsigned filter_block_size_bits>
void answer_isa_queries_aux(const pagearray<bwtsa_t<saidx_t>, pagesize_log> &bwtsa,
    long block_beg, long block_end, const std::vector<long> &queries, std::vector<long> &answers,
    unsigned char *filter) {
  for (long j = block_beg; j < block_end; ++j) {
    long sa_j = bwtsa[j].sa;
    long bit = (sa_j >> filter_block_size_bits);
    if (filter[bit >> 3] & (1 << (bit & 7))) {
      size_t pos = find(queries, sa_j);
      if (pos != queries.size())
        answers[pos] = j;
    }
  }
}

template<typename saidx_t, unsigned pagesize_log>
void answer_isa_queries(const pagearray<bwtsa_t<saidx_t>, pagesize_log> &bwtsa,
    long size, std::vector<long> &queries, std::vector<long> &answers, long max_threads) {
  if (queries.empty()) return;

  std::sort(queries.begin(), queries.end());
  queries.erase(std::unique(queries.begin(), queries.end()), queries.end());
  answers.resize(queries.size());

#if 1
  // 1
  //
  // Compute the filter.
  static const unsigned filter_block_size_bits = 14;  // best value, selected empirically
  static const unsigned filter_block_size = (1U << filter_block_size_bits);

  long n_filter_blocks = (size + filter_block_size - 1) / filter_block_size;
  unsigned char *filter = new unsigned char[(n_filter_blocks + 7) / 8];
  std::fill(filter, filter + (n_filter_blocks + 7) / 8, 0);
  for (size_t j = 0; j < queries.size(); ++j) {
    long x = queries[j];
    long bit = (x >> filter_block_size_bits);
    filter[bit >> 3] |= (1 << (bit & 7));
  }

  // 2
  //
  // Split the partial suffix array into blocks.
  long max_block_size = (size + max_threads - 1) / max_threads;
  long n_blocks = (size + max_block_size - 1) / max_block_size;

  // 3
  //
  // Compute answers to queries in parallel for each block.
  std::thread **threads = new std::thread*[n_blocks];
  for (long t = 0; t < n_blocks; ++t) {
    long block_beg = t * max_block_size;
    long block_end = std::min(block_beg + max_block_size, size);

    threads[t] = new std::thread(answer_isa_queries_aux<saidx_t, pagesize_log, filter_block_size_bits>,
        std::ref(bwtsa), block_beg, block_end, std::ref(queries), std::ref(answers), filter);
  }

  for (long t = 0; t < n_blocks; ++t) threads[t]->join();
  for (long t = 0; t < n_blocks; ++t) delete threads[t];

  // 4
  //
  // Clean up.
  delete[] threads;
  delete[] filter;
#else

#if 0
  for (long j = 0; j < size; ++j) {
    long sa_j = bwtsa[j].sa;
    size_t pos = find(queries, sa_j);
    if (pos != queries.size())
      answers[pos] = j;
  }
#else
  // 1
  //
  // compute the filter.
  static const long block_size_bits = 14;  // this is optimal choice
  static const long block_size = (1L << block_size_bits);

  long n_blocks = (size + block_size - 1) / block_size;
  unsigned char *filter = new unsigned char[(n_blocks + 7) / 8];
  std::fill(filter, filter + (n_blocks + 7) / 8, 0);
  for (size_t j = 0; j < queries.size(); ++j) {
    long x = queries[j];
    long bit = (x >> block_size_bits);
    filter[bit >> 3] |= (1 << (bit & 7));
  }

  for (long j = 0; j < size; ++j) {
    long sa_j = bwtsa[j].sa;
    long bit = (sa_j >> block_size_bits);
    if (filter[bit >> 3] & (1 << (bit & 7))) {
      size_t pos = find(queries, sa_j);
      if (pos != queries.size())
        answers[pos] = j;
    }
  }

  // 3
  //
  // Clean up.
  delete[] filter;
#endif

#endif
}

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
  std::vector<std::pair<long, long> > initial_ranges(n_threads);
  std::thread **threads = new std::thread*[n_threads];

  // 1
  //
  // Compute the last starting position. This
  // in the current form may require accessing disk.
  typedef pagearray<bwtsa_t<saidx_t>, pagesize_log> pagearray_bwtsa_type;
  long last_stream_block_beg = right_block_beg + (n_threads - 1) * max_stream_block_size;
  long last_stream_block_end = right_block_end;

  compute_last_starting_position<pagearray_bwtsa_type>(text, text_length, left_block_beg, left_block_end,
      last_stream_block_end, std::ref(bwtsa), std::ref(initial_ranks[n_threads - 1]), text_beg,
      text_end, supertext_length, supertext_filename, tail_gt_begin_reversed);

  fprintf(stderr, "%.2Lf ", utils::wclock() - start);
  start = utils::wclock();

  // 2
  //
  // Compute the starting position for all
  // starting positions other than the last one.
  long prev_stream_block_size = last_stream_block_end - last_stream_block_beg;
  for (long i = n_threads - 2; i >= 0; --i) {
    long stream_block_beg = right_block_beg + i * max_stream_block_size;
    long stream_block_end = std::min(stream_block_beg + max_stream_block_size, right_block_end);
    long stream_block_size = stream_block_end - stream_block_beg;

    // i-th thread streams text[stream_block_beg[i]..stream_block_end[i]), right-to-left.
    threads[i] = new std::thread(compute_other_starting_position<pagearray_bwtsa_type>,
        text, left_block_beg, left_block_end, stream_block_end, prev_stream_block_size,
        std::ref(bwtsa), std::ref(initial_ranges[i]));

    prev_stream_block_size = stream_block_size;
  }

  for (long i = 0; i + 1 < n_threads; ++i) threads[i]->join();
  for (long i = 0; i + 1 < n_threads; ++i) delete threads[i];
  delete[] threads;

  fprintf(stderr, "%.2Lf ", utils::wclock() - start);
  start = utils::wclock();

  // 2.5
  //
  // XXX shrink ranges to size O(n_threads).
  // XXX also, at this stage we try to minimize the intervals even further,
  //     (if we can). More precisely: when the period is very small, we should
  //     be able to determine the status (simpy by symbol comparisons) evem of
  //     the suffixes in the shrinked range.


////////////////////////////////////////////////////////////////// testing ISA queries (speed).
/*  {
    fprintf(stderr, "\nanswering isa queries: ");
    long double starting = utils::wclock();
    std::vector<long> queries;
    std::vector<long> answers;
    for (long j = 0; j < max_threads * max_threads; ++j)
      queries.push_back(utils::random_long(0, left_block_size - 1));
    answer_isa_queries(bwtsa, left_block_size, queries, answers, max_threads);
    fprintf(stderr, "%.2Lf\n", utils::wclock() - starting);
  }*/
//////////////////////////////////////////////////////////////////



  // 3
  //
  // Prepare ISA queries.
  std::vector<long> isa_queries;
  prev_stream_block_size = last_stream_block_end - last_stream_block_beg;
  for (long i = n_threads - 2; i >= 0; --i) {
    long stream_block_beg = right_block_beg + i * max_stream_block_size;
    long stream_block_end = std::min(stream_block_beg + max_stream_block_size, right_block_end);
    long stream_block_size = stream_block_end - stream_block_beg;

    long beg = initial_ranges[i].first;
    long end = initial_ranges[i].second;
    for (long j = beg + 1; j < end; ++j)
      if ((long)bwtsa[j].sa + prev_stream_block_size < left_block_size)
        isa_queries.push_back((long)bwtsa[j].sa + prev_stream_block_size);

    prev_stream_block_size = stream_block_size;
  }

  // 4
  //
  // Answer ISA queries.
  std::vector<long> isa_answers;
  answer_isa_queries(bwtsa, left_block_size, isa_queries, isa_answers, max_threads);
  std::map<long, long> isa;
  for (size_t j = 0; j < isa_queries.size(); ++j)
    isa[isa_queries[j]] = isa_answers[j];

  // 5
  //
  // Use answers to ISA queries to refine ranges to single elements.
  prev_stream_block_size = last_stream_block_end - last_stream_block_beg;
  long prev_rank = initial_ranks[n_threads - 1];
  for (long i = n_threads - 2; i >= 0; --i) {
    long stream_block_beg = right_block_beg + i * max_stream_block_size;
    long stream_block_end = std::min(stream_block_beg + max_stream_block_size, right_block_end);
    long stream_block_size = stream_block_end - stream_block_beg;
    long suf_start = stream_block_end;

    long beg = initial_ranges[i].first;
    long end = initial_ranges[i].second;

    // Invariant: the answer is in the range (beg..end].
    long ret = beg + 1;
    while (ret < end) {
      // Check if suffix starting at position stream_block_end is larger
      // than the one starting at block_beg + bwtsa[ret].sa in the text.
      // We know they have a common prefix of length prev_stream_block_size.
      if ((long)bwtsa[ret].sa + prev_stream_block_size >= left_block_size) {
        if (gt->get(suf_start + left_block_size - bwtsa[ret].sa - 1)) ++ret;
        else break;
      } else {
        long j = bwtsa[ret].sa + prev_stream_block_size;
        if (isa[j] < prev_rank) ++ret; else break;
      }
    }

    initial_ranks[i] = ret;
    prev_rank = ret;
    prev_stream_block_size = stream_block_size;
  }

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
