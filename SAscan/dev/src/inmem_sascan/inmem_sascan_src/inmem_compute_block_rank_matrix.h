#ifndef __INMEM_COMPUTE_BLOCK_RANK_MATRIX
#define __INMEM_COMPUTE_BLOCK_RANK_MATRIX

#include <cstdio>
#include <cstdlib>
#include <string>

#include "bwtsa.h"
#include "background_block_reader.h"
#include "../../multifile.h"

//#define BLOCK_MATRIX_MODULE_DEBUG_MODE


namespace inmem_sascan_private {


inline int lcp_compare_2(unsigned char *text, long text_length, unsigned char *pat, long pat_length,
    long gt_begin_length, long j, multifile_bit_stream_reader &rev_gt_begin_reader, long &lcp) {
  while (lcp < pat_length && j + lcp < text_length && pat[lcp] == text[j + lcp]) ++lcp;
  if (j + lcp >= text_length) {
    if (rev_gt_begin_reader.access(gt_begin_length - (text_length - j))) return 1;
    else return -1;
  } else if (lcp == pat_length) return 0;
  else {
    if (pat[lcp] < text[j + lcp]) return -1;
    else return 1;
  }
}

inline int lcp_compare_2_simple(unsigned char *text, unsigned char *pat, long pat_length, long j, long &lcp) {
  while (lcp < pat_length && pat[lcp] == text[j + lcp]) ++lcp;
  if (lcp == pat_length) return 0;
  else if (pat[lcp] < text[j + lcp]) return -1;
  else return 1;
}


//------------------------------------------------------------------------------
// On the entry to the function:
// - all suffixes in the range [0..left) are smaller than pat[0..old_pat_length),
// - all suffixes in the range [right..text_length) are larger than the pattern,
// - suffixes in the range [left..right) are unknown -- they can either be
//   larger or smaller than the pattern, or equal -- in any case, they have a
//   common prefix of length `lcp' with the pattern.
//------------------------------------------------------------------------------
template<typename saidx_t>
void refine_range(unsigned char *text, long block_beg, bwtsa_t<saidx_t> *block_psa, long left, long right,
    long old_pat_length, long pat_length, unsigned char *pat, long &newleft, long &newright) {
  long low = left - 1;
  long high = right;
  long llcp = old_pat_length;
  long rlcp = old_pat_length;


#ifdef BLOCK_MATRIX_MODULE_DEBUG_MODE
  long min_discrepancy = utils::random_long(0L, 10L);
  long balancing_factor = utils::random_long(1L, 10L);
#else
  static const long min_discrepancy = (1L << 16);
  static const long balancing_factor = 64L;
#endif


  while (low + 1 != high) {
    // Invariant: newleft is in the range (low, high].
    long lcp = std::min(llcp, rlcp);

    // Compute mid.
    // Valid values for mid are: low + 1, .., high - 1.
    long mid = 0L;
    if (llcp + min_discrepancy < rlcp) {
      // Choose the pivot that split the range into two
      // parts of sizes with ratio equal to logd / d.
      long d = rlcp - llcp;
      long logd = utils::log2ceil(d);
      mid = low + 1 + ((high - low - 1) * balancing_factor * logd) / (d + balancing_factor * logd);
    } else if (rlcp + min_discrepancy < llcp) {
      long d = llcp - rlcp;
      long logd = utils::log2ceil(d);
      mid = high - 1 - ((high - low - 1) * balancing_factor * logd) / (d + balancing_factor * logd);
    } else  // Discrepancy is too small, use standard binary search.
      mid = (low + high) / 2;

    if (lcp_compare_2_simple(text, pat, pat_length, block_beg + block_psa[mid].sa, lcp) <= 0) {
      high = mid;
      rlcp = lcp;
    } else {
      low = mid;
      llcp = lcp;
    }
  }

  newleft = high;

  if (rlcp < pat_length) newright = newleft;
  else {
    high = right;
    rlcp = old_pat_length;

    while (low + 1 != high) {
      // Invariant: newright is in the range (low, high].
      long lcp = std::min(llcp, rlcp);
      long mid = 0L;
      if (llcp + min_discrepancy < rlcp) {
        long d = rlcp - llcp;
        long logd = utils::log2ceil(d);
        mid = low + 1 + ((high - low - 1) * balancing_factor * logd) / (d + balancing_factor * logd);
      } else if (rlcp + min_discrepancy < llcp) {
        long d = llcp - rlcp;
        long logd = utils::log2ceil(d);
        mid = high - 1 - ((high - low - 1) * balancing_factor * logd) / (d + balancing_factor * logd);
      } else mid = (low + high) / 2;

      if (lcp_compare_2_simple(text, pat, pat_length, block_beg + block_psa[mid].sa, lcp) < 0) {
        high = mid;
        rlcp = lcp;
      } else {
        low = mid;
        llcp = lcp;
      }
    }

    newright = high;
  }
}


//=============================================================================
// NOTE: range returned by this function is NOT the range of suffixes
// starting in the block that has pat[0..pat_length) as a prefix. This would
// only be the case, if we didn't use gt_begin for the tail.
// Instead what we get is a subrange of the above: the range of suffixes
// starting inside the block that have pat[0..pat_length) as a prefix and
// we could not determine whether they are larger than the whole pattern
// or not. This is because we stop the symbol comparisons (when computing
// lcp value) at the end of the text and then use gt_begin for the tail.
// This may cause the we classify the suffix with pat[0..pat_length) as a
// prefix already at this point as larger than the whole pattern because
// gt_begin encodes that information.
//=============================================================================
template<typename saidx_t>
void refine_range(unsigned char *text, long text_length, long tail_gt_begin_reversed_length,
    long block_beg, bwtsa_t<saidx_t> *block_psa, long left, long right, multifile *tail_gt_begin_reversed,
    long old_pat_length, long pat_length, unsigned char *pat, long &newleft, long &newright) {
  multifile_bit_stream_reader reader(tail_gt_begin_reversed);


  long low = left - 1;
  long high = right;
  long llcp = old_pat_length;
  long rlcp = old_pat_length;


#ifdef BLOCK_MATRIX_MODULE_DEBUG_MODE
  long min_discrepancy = utils::random_long(0L, 10L);
  long balancing_factor = utils::random_long(1L, 10L);
#else
  static const long min_discrepancy = (1L << 16);
  static const long balancing_factor = 64L;
#endif


  while (low + 1 != high) {
    // Invariant: newleft is in the range (low, high].
    long lcp = std::min(llcp, rlcp);
    long mid = 0L;
    if (llcp + min_discrepancy < rlcp) {
      long d = rlcp - llcp;
      long logd = utils::log2ceil(d);
      mid = low + 1 + ((high - low - 1) * balancing_factor * logd) / (d + balancing_factor * logd);
    } else if (rlcp + min_discrepancy < llcp) {
      long d = llcp - rlcp;
      long logd = utils::log2ceil(d);
      mid = high - 1 - ((high - low - 1) * balancing_factor * logd) / (d + balancing_factor * logd);
    } else mid = (low + high) / 2;

    if (lcp_compare_2(text, text_length, pat, pat_length, tail_gt_begin_reversed_length,
          block_beg + block_psa[mid].sa, reader, lcp) <= 0) {
      high = mid;
      rlcp = lcp;
    } else {
      low = mid;
      llcp = lcp;
    }
  }
  newleft = high;

  if (rlcp < pat_length) newright = newleft;
  else {
    high = right;
    rlcp = old_pat_length;

    while (low + 1 != high) {
      // Invariant: newright is in the range (low, high].
      long lcp = std::min(llcp, rlcp);
      long mid = 0L;
      if (llcp + min_discrepancy < rlcp) {
        long d = rlcp - llcp;
        long logd = utils::log2ceil(d);
        mid = low + 1 + ((high - low - 1) * balancing_factor * logd) / (d + balancing_factor * logd);
      } else if (rlcp + min_discrepancy < llcp) {
        long d = llcp - rlcp;
        long logd = utils::log2ceil(d);
        mid = high - 1 - ((high - low - 1) * balancing_factor * logd) / (d + balancing_factor * logd);
      } else mid = (low + high) / 2;

      if (lcp_compare_2(text, text_length, pat, pat_length, tail_gt_begin_reversed_length,
            block_beg + block_psa[mid].sa, reader, lcp) < 0) {
        high = mid;
        rlcp = lcp;
      } else {
        low = mid;
        llcp = lcp;
      }
    }
    newright = high;
  }
}



//=============================================================================
// Variant 1: compute ranges for columns other than the last two.
//=============================================================================
template<typename saidx_t>
void compute_ranges_1(unsigned char *text, long text_length, bwtsa_t<saidx_t> *bwtsa, long max_block_size,
    std::pair<long, long> **primary_range, std::pair<long, long> **secondary_range, long row, long column) {
  long n_blocks = (text_length + max_block_size - 1) / max_block_size;
  long block_end = text_length - (n_blocks - 1 - row) * max_block_size;
  long block_begin = std::max(0L, block_end - max_block_size);
  long block_size = block_end - block_begin;
  long pat_start = text_length - (n_blocks - 1 - column) * max_block_size;
  unsigned char *pat = text + pat_start;
  bwtsa_t<saidx_t> *block_psa = bwtsa + block_begin;


  // Check that 0 <= row < column < n_blocks - 2 and pat_start + 2 * max_block_size <= text_length.
  if (0 > row || row >= column || column >= n_blocks - 2 || pat_start + 2L * max_block_size > text_length) {
    fprintf(stdout, "\nError: invariant in compute_ranges_1 failed.\n");
    std::fflush(stdout);
    std::exit(EXIT_FAILURE);
  }


  long left = 0L;
  long right = block_size;
  long cur_pat_length = 0L;


  // Compute the primary range.
  {
    long new_pat_length = max_block_size;
    if (left != right && cur_pat_length < new_pat_length) {
      long newleft = 0L;
      long newright = 0L;
      refine_range(text, block_begin, block_psa, left, right,
          cur_pat_length, new_pat_length, pat, newleft, newright);
      left = newleft;
      right = newright;
    }
    cur_pat_length = new_pat_length;
  }
  primary_range[row][column] = std::make_pair(left, right);


#ifdef BLOCK_MATRIX_MODULE_DEBUG_MODE
  // Verify the primary range.
  {
    long smaller = 0L;
    long equal = 0L;
    for (long j = block_begin; j < block_end; ++j) {
      long lcp = 0L;
      while (lcp < max_block_size && text[j + lcp] == pat[lcp]) ++lcp;
      if (lcp == max_block_size) ++equal;
      else if (text[j + lcp] < pat[lcp]) ++smaller;
    }
    long check_left = smaller;
    long check_right = smaller + equal;
    if (primary_range[row][column] != std::make_pair(check_left, check_right)) {
      fprintf(stdout, "\nError: incorrect primary range!\n");
      std::fflush(stdout);
      std::exit(EXIT_FAILURE);
    }
  }
#endif


  // Compute secondary range.
  {
    long new_pat_length = cur_pat_length + max_block_size;
    if (left != right && cur_pat_length < new_pat_length) {
      long newleft = 0L;
      long newright = 0L;
      refine_range(text, block_begin, block_psa, left, right,
          cur_pat_length, new_pat_length, pat, newleft, newright);
      left = newleft;
      right = newright;
    }
    cur_pat_length = new_pat_length;
  }
  secondary_range[row][column] = std::make_pair(left, right);


#ifdef BLOCK_MATRIX_MODULE_DEBUG_MODE
  // Verify the secondary range.
  {
    long smaller = 0L;
    long equal = 0L;
    for (long j = block_begin; j < block_end; ++j) {
      long lcp = 0L;
      while (lcp < cur_pat_length && text[j + lcp] == pat[lcp]) ++lcp;
      if (lcp == cur_pat_length) ++equal;
      else if (text[j + lcp] < pat[lcp]) ++smaller;
    }
    long check_left = smaller;
    long check_right = smaller + equal;
    if (secondary_range[row][column] != std::make_pair(check_left, check_right)) {
      fprintf(stdout, "\nError: incorrect secondary range!\n");
      std::fflush(stdout);
      std::exit(EXIT_FAILURE);
    }
  }
#endif
}


//=============================================================================
// Variant 2: compute primary and secondary range for second to last column.
//=============================================================================
template<typename saidx_t>
void compute_ranges_2(unsigned char *text, long text_length, long text_beg, long supertext_length,
    bwtsa_t<saidx_t> *bwtsa, long max_block_size, background_block_reader *reader, unsigned char *next_block,
    std::pair<long, long> **primary_range, std::pair<long, long> **secondary_range, long row, long column) {
  long text_end = text_beg + text_length;
  long tail_length = supertext_length - text_end;
  long n_blocks = (text_length + max_block_size - 1) / max_block_size;
  long block_end = text_length - (n_blocks - 1 - row) * max_block_size;
  long block_begin = std::max(0L, block_end - max_block_size);
  long block_size = block_end - block_begin;
  long pat_start = text_length - (n_blocks - 1 - column) * max_block_size;
  unsigned char *pat = text + pat_start;
  bwtsa_t<saidx_t> *block_psa = bwtsa + block_begin;


  // Check that 0 <= row < column and column == n_blocks - 2 and pat_start + max_block_size == text_length.
  if (0 > row || row >= column || column != n_blocks - 2 || pat_start + max_block_size != text_length) {
    fprintf(stdout, "\nError: invariant in compute_ranges_2 failed.\n");
    std::fflush(stdout);
    std::exit(EXIT_FAILURE);
  }


  long left = 0L;
  long right = block_size;
  long cur_pat_length = 0L;


  // Compute primary range.
  {
    long new_pat_length = max_block_size;
    if (left != right && cur_pat_length < new_pat_length) {
      long newleft = 0L;
      long newright = 0L;
      refine_range(text, block_begin, block_psa, left, right,
          cur_pat_length, new_pat_length, pat, newleft, newright);
      left = newleft;
      right = newright;
    }
    cur_pat_length = new_pat_length;
  }
  primary_range[row][column] = std::make_pair(left, right);


#ifdef BLOCK_MATRIX_MODULE_DEBUG_MODE
  // Verify the primary range.
  {
    long smaller = 0L;
    long equal = 0L;
    for (long j = block_begin; j < block_end; ++j) {
      long lcp = 0L;
      while (lcp < cur_pat_length && text[j + lcp] == pat[lcp]) ++lcp;
      if (lcp == cur_pat_length) ++equal;
      else if (text[j + lcp] < pat[lcp]) ++smaller;
    }
    long check_left = smaller;
    long check_right = smaller + equal;
    if (primary_range[row][column] != std::make_pair(check_left, check_right)) {
      fprintf(stdout, "\nError: incorrect primary range!\n");
      std::fflush(stdout);
      std::exit(EXIT_FAILURE);
    }
  }
#endif


  // Compute secondary range.
  long pat_length = cur_pat_length + std::min(tail_length, max_block_size);
  if (reader) {
    // The reader != NULL, meaning that we have to gradually refine the range.
    while (left != right && cur_pat_length < pat_length) {
      long next_chunk = std::min(reader->get_chunk_size(), pat_length - cur_pat_length);
      long new_pat_length = cur_pat_length + next_chunk;
      reader->wait(new_pat_length - max_block_size);

      long newleft = 0L;
      long newright = 0L;
      refine_range(text, block_begin, block_psa, left, right, cur_pat_length,
          new_pat_length, reader->m_data - max_block_size, newleft, newright);
      left = newleft;
      right = newright;
      cur_pat_length = new_pat_length;
    }
  } else {
#ifdef BLOCK_MATRIX_MODULE_DEBUG_MODE
    // This version extends the range chunk by chunk (using random chunk lengths)
    // even if the whole next block is available. This is for debugging purpose.
    while (left != right && cur_pat_length < pat_length) {
      long next_chunk = utils::random_long(1L, pat_length - cur_pat_length);
      long new_pat_length = cur_pat_length + next_chunk;

      long newleft = 0L;
      long newright = 0L;
      refine_range(text, block_begin, block_psa, left, right, cur_pat_length,
          new_pat_length, next_block - max_block_size, newleft, newright);
      left = newleft;
      right = newright;
      cur_pat_length = new_pat_length;
    }
#else
    // The whole next block is available, we can just do one binary search.
    long new_pat_length = pat_length;
    if (left != right && cur_pat_length < new_pat_length) {
      long newleft = 0L;
      long newright = 0L;
      refine_range(text, block_begin, block_psa, left, right, cur_pat_length,
          new_pat_length, next_block - max_block_size, newleft, newright);
      left = newleft;
      right = newright;
    }
    cur_pat_length = new_pat_length;
#endif
  }
  secondary_range[row][column] = std::make_pair(left, right);
}



//=============================================================================
// Variant 3: compute primary and secondary range for the last column.
//=============================================================================
template<typename saidx_t>
void compute_ranges_3(unsigned char *text, long text_length, long text_beg, long supertext_length, bwtsa_t<saidx_t> *bwtsa,
    long max_block_size, multifile *tail_gt_begin_reversed, background_block_reader *reader, unsigned char *next_block,
    std::pair<long, long> **primary_range, std::pair<long, long> **secondary_range, long row, long column) {
  long text_end = text_beg + text_length;
  long tail_length = supertext_length - text_end;
  long n_blocks = (text_length + max_block_size - 1) / max_block_size;
  long block_end = text_length - (n_blocks - 1 - row) * max_block_size;
  long block_beg = std::max(0L, block_end - max_block_size);
  long block_size = block_end - block_beg;
  bwtsa_t<saidx_t> *block_psa = bwtsa + block_beg;
  long first_range_pat_length = std::min(max_block_size, tail_length);
  long pat_length = std::min(text_length, tail_length);  // length of text stored in next_block (if not NULL)


  // Note: max_block_size <= text_length thus first_range_pat_length <= pat_length
  // Invariant: one of the folloing case holds:
  //   next_block != NULL and reader == NULL and next_block stores std::min(text_length, tail_length) symbols after text
  //   next_block == NULL and reader != NULL and reader will read std::min(text_length, tail_length) symbols after text


  // Check that 0 <= row < colum and column == n_blocks - 1.
  if (0 > row || row >= column || column != n_blocks - 1) {
    fprintf(stdout, "\nError: invariant 1 in compute_ranges_3 failed.\n");
    std::fflush(stdout);
    std::exit(EXIT_FAILURE);
  }


  long left = 0L;
  long right = block_size;
  long cur_pat_length = 0L;


  // Compute the primary range.
  if (reader) {
    // The reader != NULL, meaning that we have to gradually refine the range.
    while (left != right && cur_pat_length < first_range_pat_length) {
      long next_chunk = std::min(reader->get_chunk_size(), first_range_pat_length - cur_pat_length);
      long new_pat_length = cur_pat_length + next_chunk;
      reader->wait(new_pat_length);

      long newleft = 0L;
      long newright = 0L;
      refine_range(text, text_length, tail_length, block_beg, block_psa, left, right,
          tail_gt_begin_reversed, cur_pat_length, new_pat_length, reader->m_data, newleft, newright);
      left = newleft;
      right = newright;
      cur_pat_length = new_pat_length;
    }
  } else {
#ifdef BLOCK_MATRIX_MODULE_DEBUG_MODE
    // This version extends the range chunk by chunk (using random chunk lengths)
    // even if the whole next block is available. This is for debugging purpose.
    while (left != right && cur_pat_length < first_range_pat_length) {
      long next_chunk = utils::random_long(1L, first_range_pat_length - cur_pat_length);
      long new_pat_length = cur_pat_length + next_chunk;

      long newleft = 0L;
      long newright = 0L;
      refine_range(text, text_length, tail_length, block_beg, block_psa, left, right,
          tail_gt_begin_reversed, cur_pat_length, new_pat_length, next_block, newleft, newright);
      left = newleft;
      right = newright;
      cur_pat_length = new_pat_length;
    }
#else
    // The whole next block is available, we can just do one binary search.
    long new_pat_length = first_range_pat_length;
    if (left != right && cur_pat_length < new_pat_length) {
      long newleft = 0L;
      long newright = 0L;
      refine_range(text, text_length, tail_length, block_beg, block_psa, left, right,
          tail_gt_begin_reversed, cur_pat_length, new_pat_length, next_block, newleft, newright);
      left = newleft;
      right = newright;
    }
    cur_pat_length = new_pat_length;
#endif
  }
  primary_range[row][column] = std::make_pair(left, right);


  // Compute the secondary range.
  if (reader) {
    // The reader != NULL, meaning that we have to gradually refine the range.
    while (left != right && cur_pat_length < pat_length) {
      long next_chunk = std::min(reader->get_chunk_size(), pat_length - cur_pat_length);
      long new_pat_length = cur_pat_length + next_chunk;
      reader->wait(new_pat_length);

      long newleft = 0L;
      long newright = 0L;
      refine_range(text, text_length, tail_length, block_beg, block_psa, left, right,
          tail_gt_begin_reversed, cur_pat_length, new_pat_length, reader->m_data, newleft, newright);
      left = newleft;
      right = newright;
      cur_pat_length = new_pat_length;
    }
  } else {
#ifdef BLOCK_MATRIX_MODULE_DEBUG_MODE
    // This version extends the range chunk by chunk (using random chunk lengths)
    // even if the whole next block is available. This is for debugging purpose.
    while (left != right && cur_pat_length < pat_length) {
      long next_chunk = utils::random_long(1L, pat_length - cur_pat_length);
      long new_pat_length = cur_pat_length + next_chunk;

      long newleft = 0L;
      long newright = 0L;
      refine_range(text, text_length, tail_length, block_beg, block_psa, left, right,
          tail_gt_begin_reversed, cur_pat_length, new_pat_length, next_block, newleft, newright);
      left = newleft;
      right = newright;
      cur_pat_length = new_pat_length;
    }
#else
    // The whole next block is available, we can just do one binary search.
    long new_pat_length = pat_length;
    if (left != right && cur_pat_length < new_pat_length) {
      long newleft = 0L;
      long newright = 0L;
      refine_range(text, text_length, tail_length, block_beg, block_psa, left, right,
          tail_gt_begin_reversed, cur_pat_length, new_pat_length, next_block, newleft, newright);
      left = newleft;
      right = newright;
    }
    cur_pat_length = new_pat_length;
#endif
  }
  secondary_range[row][column] = std::make_pair(left, right);

  if (left != right && text_length <= tail_length) {
    fprintf(stdout, "\nError: left != right && text_length <= tail_length.\n");
    std::fflush(stdout);
    std::exit(EXIT_FAILURE);
  }
}


template<typename saidx_t>
void task_solver_code(unsigned char *text, long text_length, bwtsa_t<saidx_t> *bwtsa, long max_block_size,
    std::pair<long, long> **primary_range, std::pair<long, long> **secondary_range,
    std::vector<std::pair<long, long> > &tasks, std::mutex &tasks_mutex) {
  while (true) {
    // Get a task from the task collection.
    std::pair<long, long> task;
    bool task_avail = true;
    std::unique_lock<std::mutex> lk(tasks_mutex);
    if (tasks.empty()) task_avail = false;
    else {
      task = tasks.back();
      tasks.pop_back();
    }
    lk.unlock();

    if (!task_avail) break;

    // Solve the task and save the answer.
    compute_ranges_1(text, text_length, bwtsa, max_block_size,
        primary_range, secondary_range, task.first, task.second);
  }
}


template<typename saidx_t>
void compute_block_rank_matrix(unsigned char *text, long text_length, 
    bwtsa_t<saidx_t> *bwtsa, long max_block_size, long text_beg,
    long supertext_length, std::string,
    multifile *tail_gt_begin_reversed,  background_block_reader *reader,
    unsigned char *next_block, long **block_rank_matrix) {
  long n_blocks = (text_length + max_block_size - 1) / max_block_size;
  long text_end = text_beg + text_length;
  long tail_length = supertext_length - text_end;


  // Allocate primary and secondary ranges.
  std::pair<long, long> **primary_range = new std::pair<long, long>*[n_blocks];
  std::pair<long, long> **secondary_range = new std::pair<long, long>*[n_blocks];
  for (long row = 0; row < n_blocks; ++row) {
    primary_range[row] = new std::pair<long, long>[n_blocks];
    secondary_range[row] = new std::pair<long, long>[n_blocks];
  }


  //---------------------------------------------------------------------------
  // STEP 1: Start the threads computing ranges for the last column
  //---------------------------------------------------------------------------
  std::thread **threads_last_col = NULL;
  if (n_blocks > 1) {
    threads_last_col = new std::thread*[n_blocks - 1];
    for (long row = 0; row + 1 < n_blocks; ++row) {
      long column = n_blocks - 1;
      threads_last_col[row] = new std::thread(compute_ranges_3<saidx_t>, text, text_length, text_beg, supertext_length,
          bwtsa, max_block_size, tail_gt_begin_reversed, reader, next_block, primary_range, secondary_range, row, column);
    }
  }


  //---------------------------------------------------------------------------
  // STEP 2: Start the threads computing ranges for the second-to-last column.
  //---------------------------------------------------------------------------
  std::thread **threads_second_last_col = NULL;
  if (n_blocks > 2) {
    threads_second_last_col = new std::thread*[n_blocks - 2];
    for (long row = 0; row + 2 < n_blocks; ++row) {
      long column = n_blocks - 2;
      threads_second_last_col[row] = new std::thread(compute_ranges_2<saidx_t>, text, text_length, text_beg,
          supertext_length, bwtsa, max_block_size, reader, next_block, primary_range, secondary_range, row, column);
    }
  }


  //---------------------------------------------------------------------------
  // STEP 3: Start threads computing columns other than the last two.
  //---------------------------------------------------------------------------
  std::vector<std::pair<long, long> > tasks;
  std::mutex tasks_mutex;
  for (long row = 0; row < n_blocks; ++row)
    for (long col = row + 1; col + 2 < n_blocks; ++col)
      tasks.push_back(std::make_pair(row, col));
  std::random_shuffle(tasks.begin(), tasks.end());  // solve in no particular order
  std::thread **threads_other = new std::thread*[n_blocks];
  for (long t = 0; t < n_blocks; ++t)
    threads_other[t] = new std::thread(task_solver_code<saidx_t>, text, text_length, bwtsa,
        max_block_size, primary_range, secondary_range, std::ref(tasks), std::ref(tasks_mutex));


  //---------------------------------------------------------------------------
  // STEP 4: Wait for all threads to finish.
  //---------------------------------------------------------------------------

  // 4.1
  //
  // Wait for the threads computing columns other than last two.
  for (long t = 0; t < n_blocks; ++t) threads_other[t]->join();
  for (long t = 0; t < n_blocks; ++t) delete threads_other[t];
  delete[] threads_other;

  // 4.2
  //
  // Wait for the threads computing second-to-last column to finish.
  if (n_blocks > 2) {
    for (long row = 0; row + 2 < n_blocks; ++row) threads_second_last_col[row]->join();
    for (long row = 0; row + 2 < n_blocks; ++row) delete threads_second_last_col[row];
    delete[] threads_second_last_col;
  }

  // 4.3
  //
  // Wait for the threads computing the last column to finish.
  if (n_blocks > 1) {
    for (long row = 0; row + 1 < n_blocks; ++row) threads_last_col[row]->join();
    for (long row = 0; row + 1 < n_blocks; ++row) delete threads_last_col[row];
    delete[] threads_last_col;
  }


  //----------------------------------------------------------------------------
  // STEP 5: Compute the rank values from primary and secondary ranges.
  //----------------------------------------------------------------------------
  for (long row = n_blocks - 1; row >= 0; --row) {
    for (long col = n_blocks - 1; col > row; --col) {
      long left = secondary_range[row][col].first;
      long right = secondary_range[row][col].second;

      if (col != n_blocks - 1 && (col != n_blocks - 2 || tail_length >= max_block_size)) {
        long cur_block_end = text_length - (n_blocks - 1 - row) * max_block_size;
        long cur_block_beg = std::max(0L, cur_block_end - max_block_size);
        long cur_block_size = cur_block_end - cur_block_beg;
        long shift = max_block_size - cur_block_size;
        long next_block_end = text_length - (n_blocks - 1 - (row + 1)) * max_block_size;
        long next_block_beg = std::max(0L, next_block_end - max_block_size);
        bwtsa_t<saidx_t> *cur_block_psa = bwtsa + cur_block_beg;
        bwtsa_t<saidx_t> *next_block_psa = bwtsa + next_block_beg;

        // Compute the ranges.
        long next_primary_range_beg = primary_range[row + 1][col + 1].first;
        long next_primary_range_end = primary_range[row + 1][col + 1].second;
        long next_primary_range_size = next_primary_range_end - next_primary_range_beg;

        // Compute the difference of the arithmetic progression.
        long delta = 0L;
        long next_psa_first = 0L;
        long next_psa_second = 0L;
        if (next_primary_range_size > 1) {
          next_psa_first = next_block_psa[next_primary_range_beg].sa;
          next_psa_second = next_block_psa[next_primary_range_beg + 1].sa;
          delta = next_psa_second - next_psa_first;
        }

        // Invariant:
        // 1. the primary range of next block contains (possibly
        //    zero) values forming an arithmetic progression,
        // 2. elements in the range [left..right) of the psa of the
        //    current block incremented by `shift' appear in the primary
        //    range of the next block.

#ifdef BLOCK_MATRIX_MODULE_DEBUG_MODE
        // Check that both invariants hold.
        for (long j = next_primary_range_end; j + 1 < next_primary_range_end; ++j)
          if ((long)next_block_psa[j + 1].sa - (long)next_block_psa[j].sa != delta) {
            fprintf(stdout, "Invariant 1 failed.\n"); std::exit(EXIT_FAILURE); }
        for (long j = left; j < right; ++j) {
          long suf = cur_block_psa[j].sa + shift;
          bool found = false;
          for (long jj = next_primary_range_beg; jj < next_primary_range_end; ++jj)
            if ((long)next_block_psa[jj].sa == suf) { found = true; break; }
          if (!found) {
            fprintf(stdout, "Invariant 2 failed.\n");
            std::fflush(stdout);
            std::exit(EXIT_FAILURE);
          }
        }
#endif


        // Keep refining the range [left..right) until it's empty.
        while (left != right) {
          // Valid values for mid are in [left..right).
          long mid = (left + right) / 2;
          long suf = (long)cur_block_psa[mid].sa + shift;

          // Locate suf in next_block_psa using invariants 1. and 2.
          long pos = next_primary_range_beg;
          if (next_primary_range_size > 1)
            pos += (suf - next_psa_first) / delta;

          // Refine the range.
          if (pos < block_rank_matrix[row + 1][col + 1]) left = mid + 1;
          else right = mid;
        }
      }

      block_rank_matrix[row][col] = left;
    }
  }


  // Clean up.
  for (long row = 0; row < n_blocks; ++row) {
    delete[] primary_range[row];
    delete[] secondary_range[row];
  }
  delete[] primary_range;
  delete[] secondary_range;
}

}  // inmem_sascan_private

#endif  // __INMEM_COMPUTE_BLOCK_RANK_MATRIX
