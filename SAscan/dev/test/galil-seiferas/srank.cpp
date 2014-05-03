#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>

#include "srank.h"

namespace string_range_matching {

triple contains(std::vector<triple> &S_p, long j) {
  for (long i = 0L; i < (long)S_p.size(); ++i)
    if (S_p[i].b <= j && j < S_p[i].e) return S_p[i];
    else if (S_p[i].b > j) break;
  return triple(0L, 0L, 0L);
}

pair pred(std::vector<pair> &S_n, long j) {
  long i = 0;
  while (i + 1 < (long)S_n.size() && S_n[i + 1].b <= j) ++i;
  return S_n[i];
}


//==============================================================================
// We modify the "Precompute" from the paper. We replace
// 
//   if (i + L == m or Y[i + L] < Y[L]) count = count + 1;
//
// with
//
//   if (i + L == m) return (S_p, S_n); // exit the algorithm
//   else if (Y[i + L] < Y[L]) count = count + 1;
//
// This modification causes the algorithm to potentially compute smaller S_n
// but uncomputed elements would never be used during string range matching.
// This is because there are two cases when i + L == m may occur:
//  * Y[0..i) was a newly added k-hrp. In this case the last scope is [2i, m].
//    This means that the largest match length (during the proper string range
//    matching algorithm) that would cause a 'pred' query has to be smaller
//    than 2i (larger matches fall inside the scope) is at most 2(i-1)/k+1 < i
//    so we never actually need entries larger than i from S_n thus it's safe
//    to exit.
//  * If Y[0..i) was not a newly added k-hrp it means that either k*i > i+L = m
//    or b != 0 was holding after line 5. The second case cannot occur (because
//    it means that the scope of primitive root of Y[0..i) reaches m, i.e. the
//    algorithm should have finished before reaching this point). This leaves
//    the case k*i > m or simply i > m/k. Such i is never an argument for pred
//    thus we can just exit.
//------------------------------------------------------------------------------
// The modification above has an additional feature that it allows continuing
// Precompute procedure after appending few letters and still yields a correct
// result.
//
// Note that the case i + L == m eventually always occurs, because at some
// point the algorithm research i = m, and L = 0.
//------------------------------------------------------------------------------
// After we're done with the algorithm, we can append few letters and continue
// it by simply leaving i as it was and making simply an attempt at extending
// the match.
//------------------------------------------------------------------------------

std::pair<std::vector<triple>, std::vector<pair> > precompute(unsigned char *Y, long m) { // Y[0..m)
  std::vector<triple> S_p;
  std::vector<pair> S_n;
  S_n.push_back(pair(1,0));

  long i = 1L, last = 1L, previ = 1L;
  long L = 0L, count = 0L;
  while (i < m) {
    while (i + L < m && Y[i + L] == Y[L]) ++L;
    triple bec = contains(S_p, L);
    bool new_hrp = false;
    if (3L * i <= i + L && bec.b == 0L) {
      bec = triple(2L * i, i + L + 1L, count);
      S_p.push_back(bec);
      new_hrp = true;
    }
    if (2L * last <= i) {
      S_n.push_back(pair(i, count));
      last = i;
    }

    // if (i + L == m || Y[i + L] < Y[L]) ++count;
    if (i + L == m) {
      //---------------------------
      // Claim:
      // Either new_hrp == true or
      //        3 * i > m
      //
      // Below we verify this claim.
      //---------------------------
      if (new_hrp == false && 3L * i <= m) {
        fprintf(stderr, "Error: assumption [new_hrp == true or 3 * i > m] failed\n");
        std::exit(EXIT_FAILURE);
      }
      return std::make_pair(S_p, S_n);
    }
    else if (Y[i + L] < Y[L]) ++count;

    // Inv: count = rank(Y, Y_{[0..i]}).
    if (bec.b != 0L) {
      // Inv: per(Y[0..L)) = b/2 and
      // c = rank(Y, Y_{[1..b/2)) = rank(Y, Y_{[i+1..i+b/2)}).
      count += bec.c;
      i += bec.b / 2L;
      L -= bec.b / 2L;
    } else {
      // Inv: 3 * per(Y[0..L)) > L.
      pair bc = pred(S_n, (L / 3L) + 1L);
      // Inv: c = rank(Y, Y_{[1..s)}) = rank(Y, Y_{[i+1..i+s)}).
      // and ((j/k)+1)/4 <= s <= (j/k)+1.
      count += bc.c;
      i += bc.b;
      L = 0L;
    }

    if (previ * 2L < i) {
      fprintf(stderr, "Error: i(%ld) is bigger that 2 * prev_i(%ld)\n", i, previ);
      std::exit(EXIT_FAILURE);
    }
    previ = i;
  }

  // check if the list S_n actually satisfies the contraints
  // stated in the paper:
  for (long ii = 0L; ii + 1L < (long)S_n.size(); ++ii) {
    if (!(2L * S_n[ii].b <= S_n[ii + 1L].b && S_n[ii + 1].b < 4L * S_n[ii].b)) {
      fprintf(stderr, "Error: S_n does not satisfye the contraints stated in the paper\n");
      std::exit(EXIT_FAILURE);
    }
  }
  if (4L * S_n.back().b <= m) {
    fprintf(stderr, "Error: S_n does not satisfy the contrains stated in the paper\n");
    std::exit(EXIT_FAILURE);
  }

  return std::make_pair(S_p, S_n);
}

// identical to pseudo-code, k = 3
long strrank(unsigned char *X, long n, unsigned char *Y, long m) { // X[0..n), Y[0..m)
  std::pair<std::vector<triple>, std::vector<pair> > S = precompute(Y, m);
  std::vector<triple> S_p = S.first;
  std::vector<pair> S_n = S.second;

  long count = 0L, i = 0L, L = 0L;
  while (i < n) {
    // Inv: count = rank(Y, X_{[0..i)}).
    while (i + L < n && L < m && X[i + L] == Y[L]) ++L;
    // Inv: L = lcp(X_i, Y).
    triple bec = contains(S_p, L);
    if (L < m && (i + L == n || X[i + L] < Y[L])) ++count;
    // Inv: count = rank(Y, X_{[0..i]}).
    if (bec.b != 0L) {
      // Inv: per(Y[0..L)) = b / 2.
      // c = rank(Y, Y_{[1..b/2)}) = rank(Y, X_{[i+1..i+b/2)}).
      count += bec.c;
      i += bec.b / 2L;
      L -= bec.b / 2L;
    } else {
      // Inv: 3 * per(Y[0..L)) > L.
      pair bc = pred(S_n, (L / 3L) + 1L);
      // Inv: c = rank(Y, Y_{[1..s)}) = rank(Y, Y_{[i+1..i+s)}).
      // and ((j/k)+1)/4 <= s <= (j/k)+1.
      count += bc.c;
      i += bc.b;
      L = 0L;
    }
  }
  return count;
}

// A version of Galil-Seiferas doing string range counting with the
// S_p and S_n computed on demand.
long strrank_ondemand(unsigned char *X, long n, unsigned char *Y, long m) {
  GS_sets GS;
  long count = 0L, i = 0L, L = 0L, maxL = 0L;

  while (i < n) {
    while (i + L < n && L < m && X[i + L] == Y[L]) ++L;
    if (L > maxL) { maxL = L; GS.extend_pattern(X + i, L); }
    triple bec = contains(GS.S_p, L);
    if (L < m && (i + L == n || X[i + L] < Y[L])) ++count;
    if (bec.b != 0L) { count += bec.c; i += bec.b / 2L; L -= bec.b / 2L; }
    else { pair bc = pred(GS.S_n, (L / 3L) + 1L); count += bc.c; i += bc.b; L = 0L; }
  }

  return count;
}

long compute_smaller_suffixes(unsigned char *X, long n, unsigned char *Y, long m, unsigned char *gt) {
  GS_sets GS;
  long count = 0L, i = 0L, L = 0L, maxL = 0L;

  while (i < n) {
    while (i + L < n && L < m && X[i + L] == Y[L]) ++L;
    if (L > maxL) { maxL = L; GS.extend_pattern(X + i, L); }
    triple bec = contains(GS.S_p, L);

    if (L < m && (
                    (i + L < n && X[i + L] < Y[L]) ||
                    (i + L == n && gt[m - L - 1])
                 )
       ) ++count;

    if (bec.b != 0L) { count += bec.c; i += bec.b / 2L; L -= bec.b / 2L; }
    else { pair bc = pred(GS.S_n, (L / 3L) + 1L); count += bc.c; i += bc.b; L = 0L; }
  }

  return count;
}

// number of sufs smaller than X[s..n).
long sufrank(unsigned char *X, long n, long s) {
  return strrank(X, n, X + s, n - s);
}

} // namespace string_range_matching
