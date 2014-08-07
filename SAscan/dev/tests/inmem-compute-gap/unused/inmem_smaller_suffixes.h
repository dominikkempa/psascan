#ifndef __INMEM_SMALLER_SUFFIXES_H_INCLUDED
#define __INMEM_SMALLER_SUFFIXES_H_INCLUDED

#include "bitvector.h"

void inmem_smaller_suffixes(unsigned char *text, long text_length,
    long block_beg, long block_end, long suf_start,  bitvector *gt,
    long *ret);

#endif // __INMEM_SMALLER_SUFFIXES_H_INCLUDED
