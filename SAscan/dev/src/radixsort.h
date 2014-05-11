#ifndef __RADIXSORT_H_INCLUDED
#define __RADIXSORT_H_INCLUDED

#include <cassert>
#include <cstring>
#include <stdint.h>
#include <algorithm>

template <typename T>
static inline void radixsort8msb_copy2(T *begin, T *end, size_t depth) {
    if (end - begin < 128) {
        std::sort(begin, end);
        return;
    }

    const size_t K = 256;
    size_t bucketsize[K];
    std::memset(bucketsize, 0, K * sizeof(bucketsize[0]));

    // fill oracle and count character occurances
    uint8_t* oracle = (uint8_t*)malloc((end - begin) * sizeof(uint8_t));
    size_t ic = 0;
    for (T *it = begin; it != end; ++it, ++ic) {
        uint8_t v = ((uint8_t*)&(*it))[depth];
        ++bucketsize[ oracle[ic] = v ];
    }

    // prefix sum
    ssize_t bucketindex[K];
    bucketindex[0] = 0;
    for (size_t i = 1; i < K; ++i) {
        bucketindex[i] = bucketindex[i-1] + bucketsize[i-1];
    }

    // out-of-place permutation
    T* sorted = new T[end - begin];
    ic = 0;
    for (T *it = begin; it != end; ++it, ++ic)
        sorted[ bucketindex[ oracle[ic] ]++ ] = *it;

    free(oracle);
    std::copy(sorted, sorted + (end - begin), begin);
    delete [] sorted;

    if (depth == 0) return;

    // recursion into bucket
    T *it = begin;
    for (size_t j = 0; j < K; it += bucketsize[j++]) {
        if (bucketsize[j] <= 1) continue;
        radixsort8msb_copy2(it, it + bucketsize[j], depth - 1);
    }
}

template <typename Iterator, typename Compare>
static inline void radixsort8msb_nocopy2(Iterator begin, Iterator end, size_t depth, Compare& cmp)
{
    typedef typename Iterator::value_type value_type;

    if (end - begin < 128) {
        std::sort(begin, end, cmp);
        return;
    }

    const size_t K = 256;

    size_t bucketsize[K];
    memset(bucketsize, 0, K * sizeof(bucketsize[0]));

    // fill oracle and count character occurances
    uint8_t* oracle = (uint8_t*)malloc((end - begin) * sizeof(uint8_t));
    size_t ic = 0, jc;
    for (Iterator i = begin; i != end; ++i, ++ic) {
        uint8_t v = ((uint8_t*)&(*i))[depth];
        ++bucketsize[ oracle[ic] = v ];
    }

    // prefix sum
    ssize_t bucketindex[K];
    bucketindex[0] = bucketsize[0];
    size_t last_bucket_size = bucketsize[0];
    for (size_t i = 1; i < K; ++i) {
        bucketindex[i] = bucketindex[i-1] + bucketsize[i];
        if (bucketsize[i]) last_bucket_size = bucketsize[i];
    }

    // in-place permutation
    ic = 0;
    for (Iterator i = begin, j; i < end - last_bucket_size; )
    {
        while ( (jc = --bucketindex[ oracle[ic] ]) > ic )
        {
            j = begin + jc;
            assert( j > i );

            std::swap(*i, *j);
            std::swap(oracle[ic], oracle[jc]);
        }
        i  += bucketsize[ oracle[ic] ];
        ic += bucketsize[ oracle[ic] ];
    }
    free(oracle);

    if (depth == 0) return;

    // recursion into bucket
    Iterator i = begin;
    for (size_t j = 0; j < K; i += bucketsize[j++]) {
        if (bucketsize[j] <= 1) continue;
        radixsort8msb_nocopy2(i, i + bucketsize[j], depth-1, cmp);
    }
}

#endif // __RADIXSORT_H_INCLUDED
