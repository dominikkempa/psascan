#ifndef __INMEM_SASCAN_INMEM_SMALLER_SUFFIXES_H_INCLUDED
#define __INMEM_SASCAN_INMEM_SMALLER_SUFFIXES_H_INCLUDED

#include "../bitvector.h"
#include "pagearray.h"
#include "../multifile_bitvector.h"

namespace inmem_sascan_private {

//==============================================================================
// Return true iff text[i..length) < text[j..length). To speed up the
// computation we are given a lower bound for the lcp values, that is
// lcp <= lcp(i, j).
//==============================================================================
bool lcp_compare(unsigned char *text, long length, long i, long j, long &lcp) {
  while (i + lcp < length && j + lcp < length && text[i + lcp] == text[j + lcp])
    ++lcp;

  return (j + lcp != length &&
      (i + lcp == length || text[i + lcp] < text[j + lcp]));
}

// Return true iff text[i..) (but we always stop the comparison at text_length)
// is smaller than pat[0..pat_length).
bool lcp_compare2(unsigned char *text, long text_length, pattern &pat, long pat_length, long pat_absolute_beg,
    long supertext_length, long j, multifile_bitvector_reader &reader) {
  long lcp = 0;
  while (lcp < pat_length && j + lcp < text_length && pat[lcp] == text[j + lcp]) ++lcp;

  return (
      (lcp == pat_length) ||
      (j + lcp < text_length && pat[lcp] < text[j + lcp]) ||
      (j + lcp == text_length && !(reader.access(supertext_length - pat_absolute_beg - lcp)))
  );
}

template<typename pagearray_type>
void inmem_smaller_suffixes(unsigned char *text, long text_length,
    long block_beg, long block_end, long suf_start, const pagearray_type &partial_bwtsa, long &ret,
    long text_beg, long text_end, long supertext_length,
    std::string supertext_filename, multifile *tail_gt_begin_reversed) {

  bool has_tail = (text_end != supertext_length);

  if (!has_tail) {

    //----------------------------------------------------------------------------
    // Find the smallest j such that text[partial_sa[j]..text_length) >
    // text[suf_start..text_length) using binary search or block_size if there
    // is not such j.
    //----------------------------------------------------------------------------
    long left = -1L, right = block_end - block_beg;
    long llcp = 0L, rlcp = 0L;
    while (left + 1 != right) {
      // Invariant: llcp = lcp(left, suf_start), rlcp = lcp(right, suf_start).
      long mid = (left + right) / 2;
      long lcp = std::min(llcp, rlcp);

      if (lcp_compare(text, text_length, suf_start, block_beg + partial_bwtsa[mid].sa, lcp)) {
        right = mid;
        rlcp = lcp;
      } else {
        left = mid;
        llcp = lcp;
      }
    }

    ret = right;
  } else {

    suf_start += text_beg; // suf_start is now absolute wrt. to supertext.
    pattern pat(supertext_filename, suf_start);
    long pat_length = supertext_length - suf_start;

    multifile_bitvector_reader reader(tail_gt_begin_reversed);

    // XXX we ignore the lcp information, should not slow down too much.
    // XXX in the future we have to fix it.

    long left = -1L, right = block_end - block_beg;
    while (left + 1 != right) {
      long mid = (left + right) / 2;

      if (lcp_compare2(text, text_length, pat, pat_length, suf_start, supertext_length, block_beg + partial_bwtsa[mid].sa, reader))
        right = mid;
      else left = mid;
    }

    ret = right;
  }
}

}  // namespace inmem_sascan

#endif // __INMEM_SMALLER_SUFFIXES_H_INCLUDED
