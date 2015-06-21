#ifndef __INMEM_PSASCAN_H_INCLUDED
#define __INMEM_PSASCAN_H_INCLUDED

#include "../../src/psascan_src/inmem_psascan_src/inmem_psascan.h"
#include "../../src/psascan_src/bitvector.h"
#include "../../src/psascan_src/multifile.h"


template<typename saidx_t, unsigned pagesize_log = 12>
void inmem_psascan(
    unsigned char *text,
    long text_length,
    unsigned char *sa_bwt,
    long max_threads = 1,
    bool compute_bwt = false,
    bool compute_gt_begin = false,
    psascan_private::bitvector *gt_begin = NULL,
    long max_blocks = -1,
    long text_beg = 0,
    long text_end = 0,
    long supertext_length = 0,
    std::string supertext_filename = "",
    const psascan_private::multifile *tail_gt_begin_reversed = NULL,
    long *i0 = NULL,
    unsigned char *next_block = NULL) {
  psascan_private::inmem_psascan_private::inmem_psascan<saidx_t, pagesize_log>(text,
      text_length, sa_bwt, max_threads, compute_bwt, compute_gt_begin, gt_begin,
      max_blocks, text_beg, text_end, supertext_length, supertext_filename,
      tail_gt_begin_reversed, i0, next_block);
}

#endif  // __INMEM_PSASCAN_MAIN_HEADER
