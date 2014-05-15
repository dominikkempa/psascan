#ifndef __SUFRANK
#define __SUFRANK

#include <vector>
#include <algorithm>

namespace string_range_matching {

long strrank(unsigned char *t, long n, unsigned char *p, long m);
long strrank_ondemand(unsigned char *X, long n, unsigned char *Y, long m);

long compute_smaller_suffixes(unsigned char *X, long n, unsigned char *Y, long m, unsigned char *gt);

struct triple {
  long b, e, c;
  triple(long b_ = 0L, long e_ = 0L, long c_ = 0L)
    : b(b_), e(e_), c(c_) {}

  bool operator == (const triple &t) const {
    return b == t.b && e == t.e && c == t.c;
  }
};

struct pair {
  long b, c;
  pair(long b_ = 0L, long c_ = 0L)
    : b(b_), c(c_) {}

  bool operator == (const pair &p) const {
    return b == p.b && c == p.c;
  }
};

triple contains(std::vector<triple> &S_p, long j);
pair pred(std::vector<pair> &S_n, long j);

std::pair<std::vector<triple>, std::vector<pair> > precompute(unsigned char *Y, long m);

struct GS_sets {
  GS_sets()
      : i(1L), L(0L), last(1L), count(0L) {
    S_n.push_back(pair(1, 0));
  }

  void extend_pattern(unsigned char *new_pat, long new_length) {
    while (i < new_length) {
      while (i + L < new_length && new_pat[i + L] == new_pat[L]) ++L;

      triple bec = contains(S_p, L);

      if (3L * i <= i + L && (bec.b == 0L || S_p.back().b == 2L * i)) {
        bec = triple(2L * i, i + L + 1L, count);

        // Modification: we check if we are trying to update the most
        // recently added scope. If yes, we drop the old and add the
        // new (extended) one.
        if (!S_p.empty() && S_p.back().b == 2L * i)
          S_p.pop_back();
        S_p.push_back(bec);
      }

      if (2L * last <= i) {
        S_n.push_back(pair(i, count));
        last = i;
      }

      // Next line is crucial to ensure the correctness of this nethod.
      // Normally, it is not necessary in the Galil-Seiferas, but the
      // elements computed without exiting at this point are not used
      // in the GS matching procedure anyway.
      if (i + L == new_length) return;
      else if (new_pat[i + L] < new_pat[L]) ++count;

      if (bec.b != 0L) { count += bec.c; i += bec.b / 2L; L -= bec.b / 2L; }
      else { pair bc = pred(S_n, (L / 3L) + 1L); count += bc.c; i += bc.b; L = 0L; }
    }
  }

  long i, L, last, count;
  std::vector<triple> S_p;
  std::vector<pair> S_n;
};

} // namespace string_range_matching

#endif // __SUFRANK
