/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_RGFMC_H
#define LINE_API_PFQN_PFQN_RGFMC_H

/**
 * Multiclass Recursion by Generating Functions (RGF), with think times.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_rgfmc.m: the normalizing constant
 * of a closed MULTICLASS product-form network by eliminating one class at a
 * time by residues, finishing in the single-class convolution of pfqn_rgf.
 *
 * P. G. Harrison, S. Coury, "On the asymptotic behaviour of closed multiclass
 * queueing networks", Perf. Eval. 47:131-138, 2002, Thm 1, expresses the
 * generating function of a q-class network through those of (q-1)-class ones;
 * P. G. Harrison, T. T. Lee, "A new recursive algorithm for computing
 * generating functions in closed multi-class queueing networks", IEEE MASCOTS
 * 2004, eqs. (4)-(5), turns it into the RGF algorithm, bottoming out in
 * single-class constants memoised by load vector (Sec. 3.4).
 *
 * THINK TIMES ARE NOT IN EITHER PAPER. Both write the generating function as
 * the RATIONAL prod_i (1 - rho_i z)^-m_i, every node a load-independent single
 * server. An infinite server multiplies it by the ENTIRE exp(sum_r Z_r z_r),
 * which breaks the step Thm 1 rests on: G_n(z') = -sum_i r_i holds only because
 * the residues of n(z)/d(z) sum to zero when deg d >= deg n + 2 (Bertozzi and
 * McKenna, SIAM Review 35(2):239-268, 1993, fact (IV), p. 246), and an
 * exponential numerator does not decay at infinity. The delay is carried by
 * their own repair, eqs. (3.19)-(3.21): only the first k_r+1 Taylor
 * coefficients of exp(Z_r z_r) can reach the coefficient of z_r^k_r, so
 * replacing the exponential by that polynomial is EXACT and leaves a rational
 * integrand. The price is that the eliminated class's population re-enters the
 * term count, which is precisely the population-insensitivity Harrison-Lee
 * Sec. 4 advertises; the class kept for the base case pays nothing.
 *
 * DEGENERACY. Thm 1 assumes rho_iq != rho_lq and its Conclusion leaves the tied
 * case open. Two affine forms name the SAME pole only when PROPORTIONAL, so
 * fusing proportional forms into one factor of summed multiplicity disposes of
 * it; a tie in the eliminated class alone leaves a form with no constant term,
 * which the recursion carries unchanged.
 *
 * ARITHMETIC. The elimination is exact in exact arithmetic but is an
 * ALTERNATING sum over residues, so near-coincident loads over an eliminated
 * class destroy significance; the worst cancellation ratio is tracked and the
 * routine REFUSES past maxcancel rather than returning a confidently wrong lG.
 * Coefficients are carried in the log domain with a separate sign, so no
 * Poisson weight or binomial is ever formed as a naive ratio; that is why the
 * routine is gated on num_traits<T>::has_transcendental, as pfqn_rgf is.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/api/pfqn/pfqn_rgf.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_rgfmc, mirroring [G, lG]. */
template <class T>
struct RgfmcResult {
    T G;   ///< normalizing constant
    T lG;  ///< its logarithm
};

namespace detail {

/** coef * prod_t ( F[t][0] + sum_{p>=1} F[t][p] z_p ) ^ (-m[t]). */
template <class T>
struct RgfmcTerm {
    T lc;                          ///< log|coefficient|
    int sc;                        ///< sign of the coefficient
    std::vector<std::vector<T> > F;  ///< affine forms, one row per factor
    std::vector<int> m;            ///< pole orders
};

/** Signed log-domain sum; writes the worst cancellation ratio into cond. */
template <class T>
void rgfmc_slogsum(const std::vector<T>& lv, const std::vector<int>& sv, std::size_t n,
                   T& ls, int& sg, T& cond) {
    using std::exp;
    using std::log;
    const T ninf = num_traits<T>::from_double(-std::numeric_limits<double>::infinity());
    bool any = false;
    T mx = ninf;
    for (std::size_t i = 0; i < n; ++i) {
        if (sv[i] != 0 && lv[i] > ninf && (!any || lv[i] > mx)) {
            mx = lv[i];
            any = true;
        }
    }
    if (!any) {
        ls = ninf;
        sg = 0;
        return;
    }
    T tot = num_traits<T>::from_int(0);
    T abs = num_traits<T>::from_int(0);
    std::size_t cnt = 0;
    for (std::size_t i = 0; i < n; ++i) {
        if (sv[i] == 0 || !(lv[i] > ninf)) continue;
        const T e = exp(T(lv[i] - mx));
        tot = T(tot + num_traits<T>::from_int(sv[i]) * e);
        abs = T(abs + e);
        ++cnt;
    }
    if (tot == num_traits<T>::from_int(0)) {
        ls = ninf;
        sg = 0;
        return;
    }
    const T at = (tot < num_traits<T>::from_int(0)) ? T(-tot) : tot;
    ls = T(mx + log(at));
    sg = (tot > num_traits<T>::from_int(0)) ? 1 : -1;
    if (cnt > 1) {
        const T c = T(mx + log(abs) - ls);
        if (c > cond) cond = c;
    }
}

/** Signed log-domain linear convolution truncated at the common length. */
template <class T>
void rgfmc_slogconv(std::vector<T>& lu, std::vector<int>& su, const std::vector<T>& lv,
                    const std::vector<int>& sv, T& cond) {
    const std::size_t n = lu.size();
    std::vector<T> lo(n), tl(n);
    std::vector<int> so(n), ts(n);
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j <= k; ++j) {
            tl[j] = T(lu[j] + lv[k - j]);
            ts[j] = su[j] * sv[k - j];
        }
        rgfmc_slogsum(tl, ts, k + 1, lo[k], so[k], cond);
    }
    lu = lo;
    su = so;
}

/** log C(n,r); never a factorial quotient. */
template <class T>
T rgfmc_lbinom(const T& n, const T& r) {
    return T(num_factln<T>(n) - num_factln<T>(r) - num_factln<T>(T(n - r)));
}

/** Nonnegative integer rows of length parts summing to total. */
inline std::vector<std::vector<int> > rgfmc_compositions(int total, int parts) {
    std::vector<std::vector<int> > out;
    if (parts == 0) {
        if (total == 0) out.push_back(std::vector<int>());
        return out;
    }
    if (parts == 1) {
        out.push_back(std::vector<int>(1, total));
        return out;
    }
    for (int first = 0; first <= total; ++first) {
        std::vector<std::vector<int> > sub = rgfmc_compositions(total - first, parts - 1);
        for (std::size_t i = 0; i < sub.size(); ++i) {
            std::vector<int> row(1, first);
            row.insert(row.end(), sub[i].begin(), sub[i].end());
            out.push_back(row);
        }
    }
    return out;
}

inline int rgfmc_signpow(bool negative, int k) { return (negative && (k % 2) == 1) ? -1 : 1; }

/** Fuse PROPORTIONAL affine forms into one factor of summed multiplicity. */
template <class T>
void rgfmc_merge(std::vector<std::vector<T> >& F, std::vector<int>& m, const T& tol,
                 T& lc, int& sc) {
    using std::log;
    const T zero = num_traits<T>::from_int(0);
    const std::size_t Tn = F.size();
    std::vector<bool> used(Tn, false);
    std::vector<std::vector<T> > oF;
    std::vector<int> om;
    for (std::size_t i = 0; i < Tn; ++i) {
        if (used[i]) continue;
        used[i] = true;
        std::vector<T> Fi = F[i];
        int mi = m[i];
        std::size_t pi = 0;
        for (std::size_t c = 1; c < Fi.size(); ++c) {
            if (num_abs(Fi[c]) > num_abs(Fi[pi])) pi = c;
        }
        if (Fi[pi] == zero) throw InputError("pfqn_rgfmc: met an identically zero factor");
        for (std::size_t j = i + 1; j < Tn; ++j) {
            if (used[j]) continue;
            T nj = zero;
            for (std::size_t c = 0; c < F[j].size(); ++c) nj = std::max(nj, num_abs(F[j][c]));
            if (!(nj > zero)) continue;
            const T r = T(F[j][pi] / Fi[pi]);
            if (r == zero) continue;
            bool prop = true;
            for (std::size_t c = 0; c < F[j].size(); ++c) {
                if (num_abs(T(F[j][c] - r * Fi[c])) > tol * nj) {
                    prop = false;
                    break;
                }
            }
            if (prop) {
                lc = T(lc - num_traits<T>::from_int(m[j]) * log(num_abs(r)));
                sc *= rgfmc_signpow(r < zero, m[j]);
                mi += m[j];
                used[j] = true;
            }
        }
        oF.push_back(Fi);
        om.push_back(mi);
    }
    F = oF;
    m = om;
}

/** Eliminate the class in column col: Harrison-Coury Thm 1 as a partial fraction. */
template <class T>
std::vector<RgfmcTerm<T> > rgfmc_step(const std::vector<RgfmcTerm<T> >& terms, std::size_t col,
                                      int kr, const T& Zr, const T& tol, std::size_t maxterms,
                                      T& cond) {
    using std::log;
    const T zero = num_traits<T>::from_int(0);
    std::vector<RgfmcTerm<T> > out;
    for (std::size_t it = 0; it < terms.size(); ++it) {
        std::vector<std::vector<T> > F = terms[it].F;
        std::vector<int> m = terms[it].m;
        T lc = terms[it].lc;
        int sc = terms[it].sc;
        rgfmc_merge(F, m, tol, lc, sc);
        const std::size_t Tn = F.size();
        std::vector<std::vector<T> > A(Tn, std::vector<T>(col, zero));
        std::vector<T> B(Tn, zero);
        std::vector<bool> isMono(Tn, false);
        int shift = 0;
        for (std::size_t i = 0; i < Tn; ++i) {
            T scale = zero;
            for (std::size_t c = 0; c < F[i].size(); ++c) scale = std::max(scale, num_abs(F[i][c]));
            if (scale == zero) scale = num_traits<T>::from_int(1);
            bool mono = true;
            for (std::size_t c = 0; c < col; ++c) {
                A[i][c] = F[i][c];
                if (num_abs(F[i][c]) > tol * scale) mono = false;
            }
            B[i] = T(-F[i][col]);
            isMono[i] = mono;
            if (mono) {
                lc = T(lc - num_traits<T>::from_int(m[i]) * log(num_abs(B[i])));
                sc *= rgfmc_signpow(T(-B[i]) < zero, m[i]);
                shift += m[i];
            }
        }
        std::vector<std::size_t> S, P;
        for (std::size_t i = 0; i < Tn; ++i) {
            if (isMono[i]) continue;
            if (B[i] != zero) S.push_back(i);
            else P.push_back(i);
        }
        const int Ntot = kr + shift;
        const int pmax = (Zr > zero) ? Ntot : 0;
        for (int p = 0; p <= pmax; ++p) {
            // Bertozzi-McKenna truncation, in logs: naive Z^p/p! overflows.
            T lcz = lc;
            if (p > 0) {
                const T pT = num_traits<T>::from_int(p);
                lcz = T(lc + pT * log(Zr) - num_factln<T>(pT));
            }
            const int n = Ntot - p;
            if (S.empty()) {
                if (n == 0) {
                    RgfmcTerm<T> nt;
                    nt.lc = lcz;
                    nt.sc = sc;
                    for (std::size_t a = 0; a < P.size(); ++a) {
                        nt.F.push_back(A[P[a]]);
                        nt.m.push_back(m[P[a]]);
                    }
                    out.push_back(nt);
                }
                continue;
            }
            for (std::size_t jj = 0; jj < S.size(); ++jj) {
                const std::size_t j = S[jj];
                std::vector<std::size_t> oth;
                for (std::size_t a = 0; a < S.size(); ++a)
                    if (S[a] != j) oth.push_back(S[a]);
                const std::size_t no = oth.size();
                const T Bj = B[j];
                std::vector<std::vector<T> > Cjl(no, std::vector<T>(col, zero));
                for (std::size_t a = 0; a < no; ++a)
                    for (std::size_t c = 0; c < col; ++c)
                        Cjl[a][c] = T((A[oth[a]][c] * Bj - B[oth[a]] * A[j][c]) / Bj);
                for (int k = 0; k < m[j]; ++k) {
                    const T nT = num_traits<T>::from_int(n);
                    const T lbase = T(lcz - num_traits<T>::from_int(k) * log(num_abs(Bj)) +
                                      rgfmc_lbinom<T>(T(nT + num_traits<T>::from_int(m[j] - k - 1)), nT) +
                                      nT * log(num_abs(Bj)));
                    const int sbase = sc * rgfmc_signpow(T(-Bj) < zero, k) * rgfmc_signpow(Bj < zero, n);
                    std::vector<std::vector<int> > comps = rgfmc_compositions(k, static_cast<int>(no));
                    for (std::size_t cc = 0; cc < comps.size(); ++cc) {
                        T lt = lbase;
                        int st = sbase;
                        for (std::size_t a = 0; a < no; ++a) {
                            if (comps[cc][a] > 0) {
                                const T jlT = num_traits<T>::from_int(comps[cc][a]);
                                lt = T(lt + rgfmc_lbinom<T>(T(num_traits<T>::from_int(m[oth[a]]) + jlT -
                                                              num_traits<T>::from_int(1)), jlT) +
                                       jlT * log(num_abs(B[oth[a]])));
                                st *= rgfmc_signpow(B[oth[a]] < zero, comps[cc][a]);
                            }
                        }
                        RgfmcTerm<T> nt;
                        nt.lc = lt;
                        nt.sc = st;
                        nt.F.push_back(A[j]);
                        nt.m.push_back(n + m[j] - k);
                        for (std::size_t a = 0; a < no; ++a) {
                            nt.F.push_back(Cjl[a]);
                            nt.m.push_back(m[oth[a]] + comps[cc][a]);
                        }
                        for (std::size_t a = 0; a < P.size(); ++a) {
                            nt.F.push_back(A[P[a]]);
                            nt.m.push_back(m[P[a]]);
                        }
                        out.push_back(nt);
                    }
                }
            }
        }
        if (out.size() > maxterms)
            throw InputError("pfqn_rgfmc: exceeded maxterms; the residue term count grows as "
                             "C(S+M-1,M-1) per further elimination, use method 'ca'");
    }
    return out;
}

/** [z^N] exp(Z z) prod_t (1 - p_t z)^-m_t, signed, Coury-Harrison Property 1. */
template <class T>
void rgfmc_base_kernel(const std::vector<T>& loads, const std::vector<int>& mults, int N,
                       const T& Z, T& lg, int& sg, T& cond) {
    using std::log;
    const T zero = num_traits<T>::from_int(0);
    const T ninf = num_traits<T>::from_double(-std::numeric_limits<double>::infinity());
    const std::size_t len = static_cast<std::size_t>(N) + 1;
    std::vector<T> lv(len, ninf), lr(len, zero);
    std::vector<int> sv(len, 0), sr(len, 1);
    lv[0] = zero;
    sv[0] = 1;
    if (Z > zero) {
        for (std::size_t k = 0; k < len; ++k) {
            const T kT = num_traits<T>::from_int(static_cast<long>(k));
            lr[k] = T(kT * log(Z) - num_factln<T>(kT));
            sr[k] = 1;
        }
        rgfmc_slogconv(lv, sv, lr, sr, cond);
    }
    for (std::size_t a = 0; a < loads.size(); ++a) {
        if (loads[a] == zero) continue;
        const T ap = num_abs(loads[a]);
        const int mm = mults[a];
        for (std::size_t k = 0; k < len; ++k) {
            const T kT = num_traits<T>::from_int(static_cast<long>(k));
            if (mm == 1) {
                lr[k] = T(kT * log(ap));
            } else {
                const T mT = num_traits<T>::from_int(mm);
                lr[k] = T(num_lgamma<T>(T(kT + mT)) - num_factln<T>(kT) - num_lgamma<T>(mT) +
                          kT * log(ap));
            }
            sr[k] = (loads[a] > zero || (k % 2) == 0) ? 1 : -1;
        }
        rgfmc_slogconv(lv, sv, lr, sr, cond);
    }
    lg = lv[static_cast<std::size_t>(N)];
    sg = sv[static_cast<std::size_t>(N)];
}

/** Single-class base case, memoised on (loads, multiplicities) per Sec. 3.4. */
template <class T>
void rgfmc_base(const std::vector<RgfmcTerm<T> >& terms, int k1, const T& Z1, const T& tol,
                T& lg, int& sg, T& cond) {
    using std::log;
    const T zero = num_traits<T>::from_int(0);
    std::map<std::string, std::pair<T, int> > cache;
    std::vector<T> lv;
    std::vector<int> sv;
    for (std::size_t it = 0; it < terms.size(); ++it) {
        std::vector<std::vector<T> > F = terms[it].F;
        std::vector<int> m = terms[it].m;
        T lc = terms[it].lc;
        int sc = terms[it].sc;
        rgfmc_merge(F, m, tol, lc, sc);
        int shift = 0;
        std::vector<T> loads;
        std::vector<int> mults;
        for (std::size_t a = 0; a < F.size(); ++a) {
            const T a0 = F[a][0];
            const T a1 = F[a][1];
            const int ma = m[a];
            const T sca = std::max(num_abs(a0), num_abs(a1));
            if (sca == zero) throw InputError("pfqn_rgfmc: met an identically zero factor");
            if (num_abs(a0) <= tol * sca) {
                lc = T(lc - num_traits<T>::from_int(ma) * log(num_abs(a1)));
                sc *= rgfmc_signpow(a1 < zero, ma);
                shift += ma;
            } else {
                lc = T(lc - num_traits<T>::from_int(ma) * log(num_abs(a0)));
                sc *= rgfmc_signpow(a0 < zero, ma);
                if (num_abs(a1) > tol * sca) {
                    loads.push_back(T(-a1 / a0));
                    mults.push_back(ma);
                }
            }
        }
        const int Ntot = k1 + shift;
        std::ostringstream key;
        key << Ntot << '|' << num_traits<T>::to_double(Z1);
        for (std::size_t a = 0; a < loads.size(); ++a)
            key << '|' << num_traits<T>::to_double(loads[a]) << ':' << mults[a];
        typename std::map<std::string, std::pair<T, int> >::iterator hit = cache.find(key.str());
        T klg = zero;
        int ksg = 0;
        if (hit == cache.end()) {
            rgfmc_base_kernel(loads, mults, Ntot, Z1, klg, ksg, cond);
            cache[key.str()] = std::make_pair(klg, ksg);
        } else {
            klg = hit->second.first;
            ksg = hit->second.second;
        }
        if (ksg != 0) {
            lv.push_back(T(lc + klg));
            sv.push_back(sc * ksg);
        }
    }
    rgfmc_slogsum(lv, sv, lv.size(), lg, sg, cond);
}

}  // namespace detail

/**
 * @param L (M x R) service demands
 * @param N (R) populations, nonnegative integers
 * @param Z (R) think times
 * @param tol relative tolerance for calling two affine forms proportional
 * @param maxterms cap on residue terms carried between eliminations
 * @param maxcancel nats of cancellation tolerated before refusing
 */
template <class T>
RgfmcResult<T> pfqn_rgfmc(const Matrix<T>& L, const std::vector<int>& N,
                          const std::vector<T>& Z, const T& tol, std::size_t maxterms,
                          const T& maxcancel) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_rgfmc requires transcendental arithmetic: the residue coefficients are "
                  "carried in the log domain so that no Poisson weight or binomial is ever formed "
                  "as a naive ratio. Use pfqn_ca for the same constant in exact arithmetic");
    using std::exp;
    using std::log;
    const T zero = num_traits<T>::from_int(0);
    const std::size_t Rall = L.cols();
    if (N.size() != Rall || Z.size() != Rall)
        throw InputError("pfqn_rgfmc: requires N and Z to match the number of columns of L");
    std::vector<std::size_t> kc;
    for (std::size_t r = 0; r < Rall; ++r) {
        if (N[r] < 0 || Z[r] < zero) throw InputError("pfqn_rgfmc: requires nonnegative N and Z");
        if (N[r] > 0) kc.push_back(r);
    }
    const std::size_t R = kc.size();
    RgfmcResult<T> res;
    if (R == 0) {
        res.G = num_traits<T>::from_int(1);
        res.lG = zero;
        return res;
    }
    std::vector<std::size_t> kr;
    for (std::size_t i = 0; i < L.rows(); ++i) {
        bool any = false;
        for (std::size_t c = 0; c < R; ++c) {
            const T d = L(i, kc[c]);
            if (d < zero) throw InputError("pfqn_rgfmc: requires nonnegative demands");
            if (d > zero) any = true;
        }
        if (any) kr.push_back(i);
    }
    const std::size_t M = kr.size();
    if (M == 0) {
        T lG = zero;
        for (std::size_t c = 0; c < R; ++c) {
            const T nT = num_traits<T>::from_int(N[kc[c]]);
            lG = T(lG + nT * log(Z[kc[c]]) - detail::num_factln<T>(nT));
        }
        res.lG = lG;
        res.G = exp(lG);
        return res;
    }
    if (R == 1) {
        std::vector<T> col(M);
        for (std::size_t i = 0; i < M; ++i)
            col[i] = L(kr[i], kc[0]);
        RgfResult<T> r1 = pfqn_rgf<T>(col, N[kc[0]], Z[kc[0]]);
        res.G = r1.G;
        res.lG = r1.lG;
        return res;
    }
    // Base class = smallest population: that population is the degree the
    // sign-indefinite base series is carried to, so it drives the cancellation.
    std::vector<std::size_t> ord(R);
    for (std::size_t c = 0; c < R; ++c) ord[c] = c;
    std::stable_sort(ord.begin(), ord.end(),
                     [&](std::size_t a, std::size_t b) { return N[kc[a]] < N[kc[b]]; });
    std::vector<std::vector<T> > Ls(M, std::vector<T>(R, zero));
    std::vector<int> Ns(R);
    std::vector<T> Zs(R);
    for (std::size_t c = 0; c < R; ++c) {
        Ns[c] = N[kc[ord[c]]];
        Zs[c] = Z[kc[ord[c]]];
        for (std::size_t i = 0; i < M; ++i)
            Ls[i][c] = L(kr[i], kc[ord[c]]);
    }
    T lGscale = zero;
    for (std::size_t c = 0; c < R; ++c) {
        T cs = Zs[c];
        for (std::size_t i = 0; i < M; ++i) cs = std::max(cs, Ls[i][c]);
        if (!(cs > zero)) cs = num_traits<T>::from_int(1);
        for (std::size_t i = 0; i < M; ++i) Ls[i][c] = T(Ls[i][c] / cs);
        Zs[c] = T(Zs[c] / cs);
        lGscale = T(lGscale + num_traits<T>::from_int(Ns[c]) * log(cs));
    }
    detail::RgfmcTerm<T> t0;
    t0.lc = zero;
    t0.sc = 1;
    for (std::size_t i = 0; i < M; ++i) {
        std::vector<T> row(R + 1, zero);
        row[0] = num_traits<T>::from_int(1);
        for (std::size_t c = 0; c < R; ++c) row[c + 1] = T(-Ls[i][c]);
        t0.F.push_back(row);
        t0.m.push_back(1);
    }
    std::vector<detail::RgfmcTerm<T> > terms(1, t0);
    T cond = zero;
    // F column p+1 carries class p, so eliminating class p means column p+1.
    for (std::size_t col = R; col >= 2; --col) {
        terms = detail::rgfmc_step(terms, col, Ns[col - 1], Zs[col - 1], tol, maxterms, cond);
        if (terms.empty()) {
            res.G = zero;
            res.lG = num_traits<T>::from_double(-std::numeric_limits<double>::infinity());
            return res;
        }
    }
    T lg = zero;
    int sg = 0;
    detail::rgfmc_base(terms, Ns[0], Zs[0], tol, lg, sg, cond);
    if (sg == 0) {
        res.G = zero;
        res.lG = num_traits<T>::from_double(-std::numeric_limits<double>::infinity());
        return res;
    }
    if (sg < 0 || cond > maxcancel)
        throw InputError("pfqn_rgfmc: the residue sum cancelled past the tolerated nats, so lG "
                         "carries no significant digits. The eliminated classes have "
                         "near-coincident loads over the stations; use method 'ca'");
    res.lG = T(lg + lGscale);
    res.G = exp(res.lG);
    return res;
}

/** Overload with the reference defaults (tol 1e-12, 1e6 terms, 15 nats). */
template <class T>
RgfmcResult<T> pfqn_rgfmc(const Matrix<T>& L, const std::vector<int>& N,
                          const std::vector<T>& Z) {
    return pfqn_rgfmc<T>(L, N, Z, num_traits<T>::from_double(1e-12),
                         static_cast<std::size_t>(1000000), num_traits<T>::from_double(15.0));
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_RGFMC_H
