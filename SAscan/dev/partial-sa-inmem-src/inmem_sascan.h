#ifndef __INMEM_SASCAN_MAIN_HEADER
#define __INMEM_SASCAN_MAIN_HEADER

#include "inmem_sascan/inmem_sascan.h"
#include "bitvector.h"
#include "multifile_bitvector.h"

template<typename saidx_t, unsigned pagesize_log = 12>
void inmem_sascan(unsigned char *text, long text_length, unsigned char *sa_bwt,
    long max_threads = 1, bool compute_bwt = false, bool compute_gt_begin = false,
    bitvector *gt_begin = NULL, long max_blocks = -1,
    long text_beg = 0,
    long text_end = 0,
    long supertext_length = 0,
    std::string supertext_filename = "",
    multifile *tail_gt_begin_reversed = NULL,
    long *i0 = NULL) {
    
  inmem_sascan_private::inmem_sascan<saidx_t, pagesize_log>(text, text_length, sa_bwt,
    max_threads, compute_bwt, compute_gt_begin, gt_begin, max_blocks, text_beg,
    text_end, supertext_length, supertext_filename, tail_gt_begin_reversed, i0);
}

#endif  // __INMEM_SASCAN_MAIN_HEADER
