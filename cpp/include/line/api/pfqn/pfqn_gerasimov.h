/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_GERASIMOV_H
#define LINE_API_PFQN_GERASIMOV_H

/**
 * Gerasimov's residue (closed-form) normalizing constant, generalized to R classes.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_gerasimov.m and
 * jar/src/main/java/jline/api/pfqn/nc/Pfqn_gerasimov.java.
 *
 * A. I. Gerasimov, "On Normalizing Constants in Multiclass Queueing Networks",
 * Operations Research 43(4):704-711, 1995, evaluates
 *
 *   G(N_1,...,N_R) = (2 pi i)^-R int_G1 ... int_GR
 *                      prod_s z_s^{N_s-1} prod_i (1 - sum_s x_is/z_s)^-1
 *
 * by residues, and gives the resulting CLOSED FORM only for R = 1 (Thm 1-2) and
 * R = 2 (Thm 3 for simple poles, Thm 4 for multiple ones), stating that "for
 * three or more classes of customers, the normalizing constants can be found by
 * numerical methods". This header implements the residue elimination itself, so
 * the closed form is produced for ANY R; at R = 2 it reproduces Thm 3/4 term by
 * term.
 *
 * Written as a coefficient of the u_s = 1/z_s series,
 *
 *   G(N) = [prod_s u_s^{N_s}] exp(sum_s Z_s u_s) prod_i (1 - sum_s x_is u_s)^-1,
 *
 * every factor is AFFINE in u, so singling out u_r gives f = A - B u_r with A
 * affine in the surviving variables. Partial fractions in u_r,
 *
 *   [u_r^n] prod_j (A_j - B_j u_r)^-m_j
 *     = sum_j sum_{k=0}^{m_j-1} (-B_j)^-k C(n+m_j-k-1,n) B_j^n A_j^-(n+m_j-k)
 *       [t^k] prod_{l!=j} (C_jl - B_l t)^-m_l,   C_jl = (A_l B_j - B_l A_j)/B_j,
 *
 * map a sum of products of affine powers into another one with one variable
 * fewer, and R-1 such steps leave a univariate coefficient extraction. At R = 2
 * the single step returns one term per station i, with outer factor
 * x_i2^{N_2+M-1}/prod_{k!=i}(x_i2-x_k2), a pole of order N_2+1 at x_i1 and simple
 * poles at the paper's z_1ik = (x_k1 x_i2 - x_i1 x_k2)/(x_i2 - x_k2): exactly
 * Thm 3, with the multiple poles of Thm 4 (his xi_i < M) handled by the same
 * step. Tied x_i2, vanishing x_i2 and identical station rows, all outside the
 * paper's hypotheses, are ordinary cases here.
 *
 * Cost. Let M be the number of stations and order the populations
 * N_(1) <= ... <= N_(R). The first elimination turns the single input term into
 * M, and every later one multiplies the count by C(S+M-1,M-1) + M-1, where S is
 * the total population already eliminated: a pole of order S+1 has to be
 * differentiated against the M-1 remaining ones. The innermost extraction then
 * convolves M series of length N_(1). Hence R = 1 costs O(M N), Buzen's own cost;
 * R = 2 costs O(M^2 N_(1)^2), INDEPENDENT OF N_(2); and R >= 3 costs the same
 * times prod_{r=3}^{R} C(N_(r)+M-1, M-1). The R = 2 line is the reason to reach
 * for this method: a population removed by residues enters only as a pole ORDER,
 * i.e. through binomial coefficients, so it costs nothing at all. For R >= 3 the
 * term count is polynomial in the populations of degree (M-1)(R-2) and
 * exponential in R, which is why the paper stops at two classes and why maxterms
 * exists.
 *
 * Conditioning. The class left for the innermost extraction is the one with the
 * SMALLEST population, because that population is the degree the final,
 * sign-indefinite series is carried to, whereas the eliminated ones cancel
 * nothing; basing on N = 100 instead of N = 6 in one four-station model cost 39
 * nats of lG in double.
 *
 * Arithmetic: RATIONAL. Every operation is a product, a quotient or a binomial
 * coefficient, so the exact backend runs the residue expansion with NO
 * cancellation at all. That is the reason to have this method: in double the
 * alternating sum loses accuracy (~1e-15 typically, ~1e-9 on ill-conditioned
 * demand matrices), and the same expansion in exact arithmetic returns G to the
 * last digit. Two affine forms count as one pole when they are proportional,
 * which in the exact backend means proportional exactly, and in a floating one
 * within the relative tolerance.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

namespace detail {

/** Poisson weight Z^p/p!, in the log domain wherever the type carries logs. */
template <class T>
T geras_poisw(const T& Z, unsigned p) {
    if (p == 0) return num_traits<T>::from_int(1);
    if constexpr (num_traits<T>::has_transcendental) {
        using std::exp;
        using std::log;
        return T(exp(num_traits<T>::from_int(static_cast<long>(p)) * log(Z) -
                     num_lgamma<T>(num_traits<T>::from_int(static_cast<long>(p) + 1))));
    } else {
        return T(num_pow_int(Z, p) / num_factorial<T>(p));
    }
}


/** A product of powers of affine forms, carrying a scalar coefficient. */
template <class T>
struct GerasTerm {
    T c;
    std::vector<std::vector<T>> F;  ///< form j is F[j][0] + sum_s F[j][s] u_s
    std::vector<int> m;             ///< multiplicity of form j
};

/** base^e for a signed exponent. */
template <class T>
T geras_pow(const T& base, int e) {
    if (e >= 0) return num_pow_int(base, static_cast<unsigned>(e));
    return num_traits<T>::from_int(1) / num_pow_int(base, static_cast<unsigned>(-e));
}

/** Binomial coefficient. Exact in the exact backend, and in double corrected to
 *  the integer it is while that integer is representable. */
template <class T>
T geras_binom(int n, int k) {
    if (k < 0 || n < 0 || k > n) return num_traits<T>::from_int(0);
    const int kk = std::min(k, n - k);
    T b = num_traits<T>::from_int(1);
    for (int i = 1; i <= kk; ++i) {
        b *= num_traits<T>::from_int(n - kk + i);
        b /= num_traits<T>::from_int(i);
    }
    if constexpr (std::is_same<T, double>::value) {
        if (b < 9007199254740992.0) b = std::rint(b);
    }
    return b;
}

/** All k-tuples of nonnegative integers summing to n. */
inline std::vector<std::vector<int>> geras_compositions(int n, int k) {
    std::vector<std::vector<int>> out;
    if (k == 0) {
        if (n == 0) out.push_back(std::vector<int>());
        return out;
    }
    if (k == 1) {
        out.push_back(std::vector<int>(1, n));
        return out;
    }
    for (int a = 0; a <= n; ++a) {
        std::vector<std::vector<int>> sub = geras_compositions(n - a, k - 1);
        for (std::size_t i = 0; i < sub.size(); ++i) {
            std::vector<int> row;
            row.reserve(k);
            row.push_back(a);
            row.insert(row.end(), sub[i].begin(), sub[i].end());
            out.push_back(row);
        }
    }
    return out;
}

/**
 * Merge proportional affine forms: f_k = lambda f_j is one pole of order
 * m_j+m_k, not two nearby simple ones, and lambda^-m_k moves into the scalar.
 */
template <class T>
void geras_merge(GerasTerm<T>& t, const T& tol) {
    const std::size_t nf = t.F.size();
    if (nf <= 1) return;
    const std::size_t w = t.F[0].size();
    std::vector<bool> keep(nf, true);
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t j = 0; j < nf; ++j) {
        if (!keep[j]) continue;
        std::size_t pj = 0;
        for (std::size_t s = 1; s < w; ++s)
            if (num_abs(t.F[j][s]) > num_abs(t.F[j][pj])) pj = s;
        if (t.F[j][pj] == zero) continue;
        for (std::size_t k = j + 1; k < nf; ++k) {
            if (!keep[k]) continue;
            const T lam = t.F[k][pj] / t.F[j][pj];
            if (lam == zero) continue;
            T dev = zero, scale = zero;
            for (std::size_t s = 0; s < w; ++s) {
                dev = std::max(dev, num_abs(T(t.F[k][s] - lam * t.F[j][s])));
                scale = std::max(scale, std::max(num_abs(t.F[k][s]), num_abs(t.F[j][s])));
            }
            if (dev <= tol * scale) {
                t.c *= geras_pow(lam, -t.m[k]);
                t.m[j] += t.m[k];
                keep[k] = false;
            }
        }
    }
    std::vector<std::vector<T>> Fn;
    std::vector<int> mn;
    for (std::size_t j = 0; j < nf; ++j) {
        if (keep[j]) {
            Fn.push_back(t.F[j]);
            mn.push_back(t.m[j]);
        }
    }
    t.F.swap(Fn);
    t.m.swap(mn);
}

/**
 * One residue elimination: integrate out u_r and return the surviving sum of
 * products of affine powers, each form narrowed from r+1 to r columns.
 */
template <class T>
std::vector<GerasTerm<T>> geras_step(const std::vector<GerasTerm<T>>& terms, int r, int Nr,
                                     const T& Zr, const T& tol, std::size_t maxterms) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    std::vector<GerasTerm<T>> out;
    for (std::size_t it = 0; it < terms.size(); ++it) {
        GerasTerm<T> t = terms[it];
        geras_merge(t, tol);
        const std::size_t nf = t.F.size();
        std::vector<std::vector<T>> A(nf, std::vector<T>(static_cast<std::size_t>(r), zero));
        std::vector<T> B(nf, zero), scale(nf, zero);
        for (std::size_t j = 0; j < nf; ++j) {
            for (int s = 0; s < r; ++s) A[j][static_cast<std::size_t>(s)] = t.F[j][static_cast<std::size_t>(s)];
            B[j] = -t.F[j][static_cast<std::size_t>(r)];
            for (int s = 0; s <= r; ++s) scale[j] = std::max(scale[j], num_abs(t.F[j][static_cast<std::size_t>(s)]));
        }
        T c = t.c;
        int shift = 0;
        std::vector<bool> drop(nf, false);
        for (std::size_t j = 0; j < nf; ++j) {
            bool mono = true;
            for (int s = 0; s < r && mono; ++s)
                if (num_abs(A[j][static_cast<std::size_t>(s)]) > tol * scale[j]) mono = false;
            if (!mono) continue;
            if (num_abs(B[j]) <= tol * scale[j])
                throw InputError("pfqn_gerasimov: identically zero factor, impossible after merging proportional ones");
            // A factor -B u_r carries no finite pole: it only shifts the exponent.
            c *= geras_pow(T(-B[j]), -t.m[j]);
            shift += t.m[j];
            drop[j] = true;
        }
        std::vector<std::size_t> S, P;  // S carries a pole in u_r, P is free of u_r
        for (std::size_t j = 0; j < nf; ++j) {
            if (drop[j]) continue;
            if (B[j] != zero) S.push_back(j); else P.push_back(j);
        }
        const int Ntot = Nr + shift;
        const int pmax = (Zr > zero) ? Ntot : 0;
        for (int p = 0; p <= pmax; ++p) {
            T cz = c;
            // Poisson weight Z_r^p/p!. In double the naive ratio overflows for
            // p >~ 171, reachable when the eliminated class carries think time;
            // an exact T cannot overflow and must keep the exact quotient.
            if (p > 0) cz = c * geras_poisw<T>(Zr, static_cast<unsigned>(p));
            const int Neff = Ntot - p;
            if (S.empty()) {
                if (Neff == 0) {
                    GerasTerm<T> nt;
                    nt.c = cz;
                    for (std::size_t q = 0; q < P.size(); ++q) {
                        nt.F.push_back(A[P[q]]);
                        nt.m.push_back(t.m[P[q]]);
                    }
                    out.push_back(nt);
                }
                continue;
            }
            for (std::size_t jj = 0; jj < S.size(); ++jj) {
                const std::size_t j = S[jj];
                std::vector<std::size_t> oth;
                for (std::size_t q = 0; q < S.size(); ++q)
                    if (q != jj) oth.push_back(S[q]);
                const std::size_t no = oth.size();
                std::vector<std::vector<T>> Cjl(no, std::vector<T>(static_cast<std::size_t>(r), zero));
                for (std::size_t l = 0; l < no; ++l) {
                    const std::size_t k = oth[l];
                    for (int s = 0; s < r; ++s) {
                        const std::size_t ss = static_cast<std::size_t>(s);
                        Cjl[l][ss] = (A[k][ss] * B[j] - B[k] * A[j][ss]) / B[j];
                    }
                }
                for (int k = 0; k < t.m[j]; ++k) {
                    const std::vector<std::vector<int>> comps = geras_compositions(k, static_cast<int>(no));
                    for (std::size_t ci = 0; ci < comps.size(); ++ci) {
                        const std::vector<int>& nk = comps[ci];
                        T coef = cz * geras_pow(T(-B[j]), -k) * geras_binom<T>(Neff + t.m[j] - k - 1, Neff) *
                                 geras_pow(B[j], Neff);
                        for (std::size_t l = 0; l < no; ++l)
                            coef *= geras_binom<T>(t.m[oth[l]] + nk[l] - 1, nk[l]) * geras_pow(B[oth[l]], nk[l]);
                        if (coef == zero) continue;
                        GerasTerm<T> nt;
                        nt.c = coef;
                        nt.F.push_back(A[j]);
                        nt.m.push_back(Neff + t.m[j] - k);
                        for (std::size_t l = 0; l < no; ++l) {
                            nt.F.push_back(Cjl[l]);
                            nt.m.push_back(t.m[oth[l]] + nk[l]);
                        }
                        for (std::size_t q = 0; q < P.size(); ++q) {
                            nt.F.push_back(A[P[q]]);
                            nt.m.push_back(t.m[P[q]]);
                        }
                        out.push_back(nt);
                    }
                }
            }
        }
        if (out.size() > maxterms)
            throw InputError("pfqn_gerasimov: residue expansion exceeded maxterms; use pfqn_ca or pfqn_nc");
    }
    (void)one;
    return out;
}

/** Truncated convolution of two coefficient sequences, kept to degree n. */
template <class T>
std::vector<T> geras_conv(const std::vector<T>& a, const std::vector<T>& b, int n) {
    const T zero = num_traits<T>::from_int(0);
    std::vector<T> y(static_cast<std::size_t>(n) + 1, zero);
    for (std::size_t i = 0; i < a.size() && static_cast<int>(i) <= n; ++i) {
        if (a[i] == zero) continue;
        for (std::size_t j = 0; j < b.size() && static_cast<int>(i + j) <= n; ++j)
            y[i + j] += a[i] * b[j];
    }
    return y;
}

/**
 * Last class: a univariate coefficient extraction. Summing the residues here too
 * would repeat the step above, but convolving the series of each factor returns
 * the same number without expanding the multiple poles.
 */
template <class T>
T geras_base(const std::vector<GerasTerm<T>>& terms, int N1, const T& Z1, const T& tol) {
    const T zero = num_traits<T>::from_int(0);
    T G = zero;
    for (std::size_t it = 0; it < terms.size(); ++it) {
        GerasTerm<T> t = terms[it];
        geras_merge(t, tol);
        const std::size_t nf = t.F.size();
        T c = t.c;
        int shift = 0;
        std::vector<bool> drop(nf, false);
        for (std::size_t j = 0; j < nf; ++j) {
            const T Aj = t.F[j][0];
            const T Bj = -t.F[j][1];
            const T scale = std::max(num_abs(t.F[j][0]), num_abs(t.F[j][1]));
            if (num_abs(Aj) > tol * scale) continue;
            if (num_abs(Bj) <= tol * scale)
                throw InputError("pfqn_gerasimov: identically zero factor at the innermost coefficient extraction");
            c *= geras_pow(T(-Bj), -t.m[j]);
            shift += t.m[j];
            drop[j] = true;
        }
        const int Ntot = N1 + shift;
        const std::size_t len = static_cast<std::size_t>(Ntot) + 1;
        std::vector<T> s(len, zero);
        s[0] = num_traits<T>::from_int(1);
        if (Z1 > zero) {
            std::vector<T> pois(len, zero);
            for (std::size_t n = 0; n < len; ++n)
                pois[n] = geras_poisw<T>(Z1, static_cast<unsigned>(n));
            s = geras_conv(s, pois, Ntot);
        }
        for (std::size_t j = 0; j < nf; ++j) {
            if (drop[j]) continue;
            const T Aj = t.F[j][0];
            const T Bj = -t.F[j][1];
            c *= geras_pow(Aj, -t.m[j]);
            if (Bj == zero) continue;
            const T ratio = Bj / Aj;
            std::vector<T> seq(len, zero);
            for (std::size_t n = 0; n < len; ++n)
                seq[n] = geras_binom<T>(t.m[j] + static_cast<int>(n) - 1, static_cast<int>(n)) *
                         num_pow_int(ratio, static_cast<unsigned>(n));
            s = geras_conv(s, seq, Ntot);
        }
        G += c * s[static_cast<std::size_t>(Ntot)];
    }
    return G;
}

}  // namespace detail

/**
 * Exact normalizing constant of a closed multiclass product-form network by
 * ITERATED RESIDUES of its rational generating function, one class at a time.
 *
 * @param L  (M x R) service demands, M queueing stations, R classes
 * @param N  (R) population per class, nonnegative
 * @param Z  (K x R) think times, summed over rows; may be empty. A delay
 *           contributes the entire factor exp(sum_s Z_s u_s), handled exactly by
 *           convolving its Poisson coefficients into each elimination.
 * @param tol relative tolerance for declaring two affine forms proportional,
 *           hence one pole rather than two. Ignored (taken as exactly zero) in
 *           the exact backend.
 * @param maxterms cap on the residue terms carried between eliminations.
 *           Exceeding it is an error, not a truncation: a truncated residue sum
 *           is not a bound or an approximation of G, it is a wrong number.
 */
template <class T>
NcResult<T> pfqn_gerasimov(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z,
                           double tol = 1e-12, std::size_t maxterms = 200000) {
    const std::size_t R0 = N.size();
    if (!L.empty() && L.cols() != R0)
        throw InputError("pfqn_gerasimov: L and N disagree on the class count");
    const T zero = num_traits<T>::from_int(0);
    std::vector<T> Zsum(R0, zero);
    if (!Z.empty()) {
        if (Z.cols() != R0) throw InputError("pfqn_gerasimov: Z and N disagree on the class count");
        for (std::size_t k = 0; k < Z.rows(); ++k)
            for (std::size_t r = 0; r < R0; ++r) Zsum[r] += Z(k, r);
    }
    for (std::size_t r = 0; r < R0; ++r) {
        if (N[r] < 0 || Zsum[r] < zero) throw InputError("pfqn_gerasimov: L, N and Z must be nonnegative");
    }
    for (std::size_t i = 0; i < L.rows(); ++i)
        for (std::size_t r = 0; r < R0; ++r)
            if (L(i, r) < zero) throw InputError("pfqn_gerasimov: L, N and Z must be nonnegative");

    // A class with no jobs is eliminated by evaluating the generating function at
    // u_r = 0, i.e. by deleting its column outright.
    std::vector<std::size_t> cls;
    for (std::size_t r = 0; r < R0; ++r)
        if (N[r] > 0) cls.push_back(r);
    if (cls.empty()) return {num_traits<T>::from_int(1), 0.0};
    // Class order, which decides both the cost and the accuracy.
    //  - The class left for the innermost extraction sets the CONDITIONING. Its
    //    population is the degree the final series is carried to, and the poles of
    //    the reduced problem have arbitrary sign, so that series cancels; the
    //    populations eliminated by residues enter only as pole ORDERS, through
    //    binomial coefficients, and cancel nothing. Basing on N = 100 rather than
    //    on N = 6 in one 4-station model cost 39 nats of lG in double. The SMALLEST
    //    population therefore goes to the base.
    //  - Eliminating class r leaves a pole of order N_r+1 that every LATER
    //    elimination has to differentiate, so the remaining classes are eliminated
    //    smallest-first to keep the multiplicities low for as long as possible.
    // Eliminations run from index R down to 2, so indices 2..R hold the remaining
    // populations in DECREASING order and index 1 holds the smallest.
    std::stable_sort(cls.begin(), cls.end(),
                     [&N](std::size_t a, std::size_t b) { return N[a] < N[b]; });
    if (cls.size() > 2) std::reverse(cls.begin() + 1, cls.end());
    const std::size_t R = cls.size();

    // Per-class scaling. The residue coefficients carry x_ir^{N_r+M-1}, which in
    // double overflows well before G itself does: at x = 4 and N_r = 400 the factor
    // alone is 1e240 while G is finite. Dividing column r by c_r divides G by
    // exactly c_r^N_r (substitute u_r -> u_r/c_r in the generating function), so
    // the scaling is exact in every backend and is undone at the end.
    std::vector<T> cs(R, num_traits<T>::from_int(1));
    for (std::size_t r = 0; r < R; ++r) {
        T c = Zsum[cls[r]];
        for (std::size_t i = 0; i < L.rows(); ++i) c = std::max(c, L(i, cls[r]));
        if (c > zero) cs[r] = c;
    }

    // A station with no demand at all contributes the factor 1.
    std::vector<std::vector<T>> rows;
    for (std::size_t i = 0; i < L.rows(); ++i) {
        std::vector<T> row(R + 1, zero);
        row[0] = num_traits<T>::from_int(1);
        bool any = false;
        for (std::size_t r = 0; r < R; ++r) {
            const T v = L(i, cls[r]);
            row[r + 1] = -(v / cs[r]);
            if (v > zero) any = true;
        }
        if (any) rows.push_back(row);
    }

    detail::GerasTerm<T> t0;
    t0.c = num_traits<T>::from_int(1);
    t0.F = rows;
    t0.m.assign(rows.size(), 1);
    std::vector<detail::GerasTerm<T>> terms(1, t0);

    const T tolT = num_traits<T>::is_exact ? zero : num_traits<T>::from_double(tol);
    // class s (1-based) sits in column s of F = [1, -L]; eliminate R, R-1, ..., 2
    for (std::size_t r = R; r >= 2; --r) {
        terms = detail::geras_step(terms, static_cast<int>(r), N[cls[r - 1]],
                                   T(Zsum[cls[r - 1]] / cs[r - 1]), tolT, maxterms);
        if (terms.empty()) return {zero, -std::numeric_limits<double>::infinity()};
    }
    const T Gs = detail::geras_base(terms, N[cls[0]], T(Zsum[cls[0]] / cs[0]), tolT);
    if (Gs == zero) return {zero, -std::numeric_limits<double>::infinity()};
    T unscale = num_traits<T>::from_int(1);
    double lGscale = 0.0;
    for (std::size_t r = 0; r < R; ++r) {
        unscale *= num_pow_int(cs[r], static_cast<unsigned>(N[cls[r]]));
        lGscale += N[cls[r]] * std::log(num_traits<T>::to_double(cs[r]));
    }
    const double lGout = num_traits<T>::log_as_double(Gs) + lGscale;
    if constexpr (num_traits<T>::has_transcendental) {
        // exp(lGout), not Gs*unscale: in double the product Prod cs_r^N_r
        // overflows on its own whenever the scaling is large, even when G is
        // inside the range. An exact T cannot overflow and keeps the product.
        using std::exp;
        return {T(exp(num_traits<T>::from_double(lGout))), lGout};
    } else {
        return {T(Gs * unscale), lGout};
    }
}

/** Delay-free overload. */
template <class T>
NcResult<T> pfqn_gerasimov(const Matrix<T>& L, const std::vector<int>& N) {
    return pfqn_gerasimov(L, N, Matrix<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_GERASIMOV_H
