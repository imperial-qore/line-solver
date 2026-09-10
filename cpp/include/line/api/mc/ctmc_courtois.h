/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_COURTOIS_H
#define LINE_API_MC_CTMC_COURTOIS_H

/**
 * Courtois decomposition of a nearly completely decomposable (NCD) CTMC.
 *
 * Templated port of matlab/src/api/mc/ctmc_courtois.m and
 * jar/src/main/java/jline/api/mc/Ctmc_courtois.java. The states are permuted
 * into macro-state order; the coupling between macro-states is deleted to give
 * a block-diagonal generator Qdec whose blocks are solved independently for
 * the conditional (micro) distributions; the macro-state chain G is assembled
 * from the uniformized matrix weighted by those micro distributions; and the
 * approximation is the product of the macro and micro probabilities, mapped
 * back to the original state ordering.
 *
 * The approximation error is governed by the degree of coupling eps, which is
 * meaningful only when it is small against epsMAX = (1 - max_i mu_i)/2, mu_i
 * being the subdominant eigenvalue modulus of the i-th diagonal block made
 * stochastic. Both diagnostics are always computed here; MATLAB computes them
 * only when more than three outputs are requested and returns eps = 0
 * otherwise, which is a trap for a caller who asks for two outputs and reads
 * the second as a coupling measure.
 *
 * WHICH SUM eps IS, and the divergence that used to live here. The NCD index of
 * Courtois is the largest ROW sum of the coupling matrix, ||B||_inf, the
 * largest probability of leaving a macro-state in one uniformized step. MATLAB
 * once wrote it `max(sum(B))`, and `sum(B)` on a matrix is the vector of COLUMN
 * sums, so its eps was the largest column sum instead; `Ctmc_courtois.java`
 * computed the row sum (`B.sumRows().elementMax()`) and the two separated
 * whenever any state carried coupling to more than one other macro-state.
 * MATLAB is `max(sum(B,2))` as of 2026-08-15 and the three codebases agree, so
 * `eps` here is the ROW sum. `epsColMax` keeps the old column-sum value,
 * because a golden or a comparison recorded before that date holds it: on the
 * six-state fixture of test_mc_aggregation the row sum is 0.0076190 and the
 * column sum 0.013333, a factor of 1.75, and the column sum is the OPTIMISTIC
 * one -- it understates the coupling, which is the wrong direction for a
 * diagnostic answering "is this partition decomposable enough to trust".
 *
 * GATED ON TRANSCENDENTAL ARITHMETIC. epsMAX is a subdominant eigenvalue
 * modulus, obtained by the Francis QR iteration, which is iterative and stops
 * on a tolerance; it has no closed form in the field of the rates (the
 * eigenvalues of a rational matrix are algebraic, not rational). The rest of
 * the method -- permutation, block solves, aggregation, disaggregation -- is
 * finite and exact, so a caller who needs the approximate stationary vector at
 * Rational and no diagnostics can assemble it from ctmc_solve_reducible and
 * dtmc_solve directly.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/mc/ctmc_randomization.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/api/mc/ctmc_solve_reducible.h"
#include "line/api/mc/dtmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

namespace detail {

/**
 * Moduli of the eigenvalues of a real square matrix, ascending.
 *
 * Householder-free elimination to upper Hessenberg form followed by the
 * Francis implicit double-shift QR iteration, in real arithmetic, so a complex
 * conjugate pair is deflated as a 2 x 2 block and never needs complex storage.
 * This is the classical hqr of Wilkinson and Reinsch. The convergence tests are
 * relative to the machine epsilon of T, so raising the precision of T tightens
 * them rather than leaving them at double resolution.
 *
 * Indices are one-based throughout, as in the original algorithm; the working
 * array is therefore (n+1) x (n+1) with row and column zero unused. Departing
 * from that indexing is the classic way to introduce an off-by-one into this
 * particular routine.
 */
template <class T>
std::vector<T> eig_moduli(const Matrix<T>& Ain) {
    const std::size_t n = Ain.rows();
    if (Ain.cols() != n) throw InputError("eig_moduli: matrix is not square");
    const T zero = num_traits<T>::from_int(0);
    if (n == 0) return std::vector<T>();
    if (n == 1) return std::vector<T>(1, num_abs(T(Ain(0, 0))));

    Matrix<T> a(n + 1, n + 1, zero);
    for (std::size_t i = 1; i <= n; ++i)
        for (std::size_t j = 1; j <= n; ++j) a(i, j) = Ain(i - 1, j - 1);

    // Reduction to upper Hessenberg form by elimination with pivoting.
    for (std::size_t m = 2; m < n; ++m) {
        T x = zero;
        std::size_t piv = m;
        for (std::size_t j = m; j <= n; ++j)
            if (num_abs(T(a(j, m - 1))) > num_abs(T(x))) {
                x = a(j, m - 1);
                piv = j;
            }
        if (piv != m) {
            for (std::size_t j = m - 1; j <= n; ++j) std::swap(a(piv, j), a(m, j));
            for (std::size_t j = 1; j <= n; ++j) std::swap(a(j, piv), a(j, m));
        }
        if (x != zero) {
            for (std::size_t i = m + 1; i <= n; ++i) {
                T y = a(i, m - 1);
                if (y == zero) continue;
                y /= x;
                a(i, m - 1) = y;
                for (std::size_t j = m; j <= n; ++j) a(i, j) -= y * a(m, j);
                for (std::size_t j = 1; j <= n; ++j) a(j, m) += y * a(j, i);
            }
        }
    }
    // The elimination multipliers left below the subdiagonal are not part of
    // the Hessenberg matrix and must go before the QR iteration reads them.
    for (std::size_t i = 3; i <= n; ++i)
        for (std::size_t j = 1; j + 1 < i; ++j) a(i, j) = zero;

    const T epsT = std::numeric_limits<T>::epsilon();
    T anorm = zero;
    for (std::size_t i = 1; i <= n; ++i)
        for (std::size_t j = (i > 1 ? i - 1 : 1); j <= n; ++j) anorm += num_abs(T(a(i, j)));

    std::vector<T> wr(n + 1, zero), wi(n + 1, zero);
    std::size_t nn = n;
    T t = zero;
    while (nn >= 1) {
        int its = 0;
        std::size_t l = 0;
        do {
            for (l = nn; l >= 2; --l) {
                T s = num_abs(T(a(l - 1, l - 1))) + num_abs(T(a(l, l)));
                if (s == zero) s = anorm;
                if (num_abs(T(a(l, l - 1))) <= epsT * s) {
                    a(l, l - 1) = zero;
                    break;
                }
            }
            if (l < 2) l = 1;
            T x = a(nn, nn);
            if (l == nn) {
                wr[nn] = x + t;
                wi[nn] = zero;
                --nn;
            } else {
                T y = a(nn - 1, nn - 1);
                T w = a(nn, nn - 1) * a(nn - 1, nn);
                if (l == nn - 1) {
                    const T p = (y - x) / num_traits<T>::from_int(2);
                    const T q = p * p + w;
                    using std::sqrt;
                    T z = sqrt(T(num_abs(T(q))));
                    x += t;
                    if (!(q < zero)) {
                        z = p + (p < zero ? T(-z) : z);
                        wr[nn - 1] = x + z;
                        wr[nn] = (z != zero) ? T(x - w / z) : T(x + z);
                        wi[nn - 1] = zero;
                        wi[nn] = zero;
                    } else {
                        wr[nn - 1] = x + p;
                        wr[nn] = x + p;
                        wi[nn] = z;
                        wi[nn - 1] = -z;
                    }
                    nn -= 2;
                } else {
                    if (its == 60) throw NumericError("eig_moduli: QR iteration did not converge");
                    if (its == 10 || its == 20 || its == 30 || its == 40 || its == 50) {
                        // Exceptional shift, to break a cycle the Wilkinson
                        // shift cannot resolve.
                        t += x;
                        for (std::size_t i = 1; i <= nn; ++i) a(i, i) -= x;
                        const T s = num_abs(T(a(nn, nn - 1))) + num_abs(T(a(nn - 1, nn - 2)));
                        x = num_traits<T>::from_rational(3, 4) * s;
                        y = x;
                        w = -num_traits<T>::from_rational(7, 16) * s * s;
                    }
                    ++its;
                    std::size_t m = nn - 2;
                    T p = zero, q = zero, r = zero;
                    for (; m >= l; --m) {
                        const T z = a(m, m);
                        const T rr = x - z;
                        const T ss = y - z;
                        p = (rr * ss - w) / a(m + 1, m) + a(m, m + 1);
                        q = a(m + 1, m + 1) - z - rr - ss;
                        r = a(m + 2, m + 1);
                        const T s = num_abs(T(p)) + num_abs(T(q)) + num_abs(T(r));
                        p /= s;
                        q /= s;
                        r /= s;
                        if (m == l) break;
                        const T u = num_abs(T(a(m, m - 1))) * (num_abs(T(q)) + num_abs(T(r)));
                        const T v = num_abs(T(p)) *
                                    (num_abs(T(a(m - 1, m - 1))) + num_abs(T(z)) + num_abs(T(a(m + 1, m + 1))));
                        if (u <= epsT * v) break;
                    }
                    for (std::size_t i = m + 2; i <= nn; ++i) {
                        a(i, i - 2) = zero;
                        if (i != m + 2) a(i, i - 3) = zero;
                    }
                    for (std::size_t k = m; k <= nn - 1; ++k) {
                        if (k != m) {
                            p = a(k, k - 1);
                            q = a(k + 1, k - 1);
                            r = zero;
                            if (k != nn - 1) r = a(k + 2, k - 1);
                            x = num_abs(T(p)) + num_abs(T(q)) + num_abs(T(r));
                            if (x != zero) {
                                p /= x;
                                q /= x;
                                r /= x;
                            }
                        }
                        using std::sqrt;
                        T s = sqrt(T(p * p + q * q + r * r));
                        if (p < zero) s = -s;
                        if (s == zero) continue;
                        if (k == m) {
                            if (l != m) a(k, k - 1) = -a(k, k - 1);
                        } else {
                            a(k, k - 1) = -s * x;
                        }
                        p += s;
                        x = p / s;
                        y = q / s;
                        const T z = r / s;
                        q /= p;
                        r /= p;
                        for (std::size_t j = k; j <= nn; ++j) {
                            T pp = a(k, j) + q * a(k + 1, j);
                            if (k != nn - 1) {
                                pp += r * a(k + 2, j);
                                a(k + 2, j) -= pp * z;
                            }
                            a(k + 1, j) -= pp * y;
                            a(k, j) -= pp * x;
                        }
                        const std::size_t mmin = nn < k + 3 ? nn : k + 3;
                        for (std::size_t i = l; i <= mmin; ++i) {
                            T pp = x * a(i, k) + y * a(i, k + 1);
                            if (k != nn - 1) {
                                pp += z * a(i, k + 2);
                                a(i, k + 2) -= pp * r;
                            }
                            a(i, k + 1) -= pp * q;
                            a(i, k) -= pp;
                        }
                    }
                }
            }
        } while (nn >= 2 && l < nn - 1);
    }

    std::vector<T> mod(n);
    using std::sqrt;
    for (std::size_t i = 1; i <= n; ++i) mod[i - 1] = sqrt(T(wr[i] * wr[i] + wi[i] * wi[i]));
    std::sort(mod.begin(), mod.end());
    return mod;
}

/** Concatenation of the macro-state index sets, the permutation v of MATLAB. */
inline std::vector<std::size_t> macrostate_permutation(const std::vector<std::vector<std::size_t>>& MS,
                                                       std::size_t n) {
    std::vector<std::size_t> v;
    std::vector<char> seen(n, 0);
    for (const std::vector<std::size_t>& b : MS)
        for (std::size_t k : b) {
            if (k >= n) throw InputError("ctmc_courtois: macro-state index out of range");
            if (seen[k]) throw InputError("ctmc_courtois: state listed in more than one macro-state");
            seen[k] = 1;
            v.push_back(k);
        }
    if (v.size() != n) throw InputError("ctmc_courtois: the macro-states do not cover every state");
    return v;
}

/** Zeroes every entry outside the diagonal blocks of the given sizes. */
template <class T>
void zero_offblock(Matrix<T>& A, const std::vector<std::vector<std::size_t>>& MS) {
    const std::size_t n = A.rows();
    const T zero = num_traits<T>::from_int(0);
    std::size_t proc = 0;
    for (const std::vector<std::size_t>& b : MS) {
        const std::size_t sz = b.size();
        for (std::size_t row = proc; row < proc + sz; ++row) {
            for (std::size_t col = 0; col < proc; ++col) A(row, col) = zero;
            for (std::size_t col = proc + sz; col < n; ++col) A(row, col) = zero;
        }
        proc += sz;
    }
}

/**
 * Everything the Courtois construction produces before the macro-state chain is
 * solved. ctmc_multi differs from ctmc_courtois only in HOW that chain is
 * solved -- by a second Courtois decomposition rather than directly -- so the
 * shared part lives here and neither routine reimplements it.
 */
template <class T>
struct CourtoisCore {
    std::vector<std::size_t> v;  ///< macro-state-major permutation
    Matrix<T> Qperm, Qdec, P, B;
    std::vector<T> pmicro;  ///< conditional distributions, permuted order
    Matrix<T> G;            ///< macro-state transition matrix
    T eps, epsRowMax, epsColMax, epsMAX, q;
};

template <class T>
CourtoisCore<T> courtois_core(const Matrix<T>& Q, const std::vector<std::vector<std::size_t>>& MS,
                              const T& q) {
    const std::size_t n = Q.rows();
    if (Q.cols() != n) throw InputError("ctmc_courtois: generator is not square");
    if (MS.empty()) throw InputError("ctmc_courtois: no macro-states given");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t nMacro = MS.size();

    CourtoisCore<T> r;
    r.q = q;
    r.v = macrostate_permutation(MS, n);
    r.Qperm = submatrix(Q, r.v);

    r.Qdec = r.Qperm;
    zero_offblock(r.Qdec, MS);
    r.Qdec = ctmc_makeinfgen(r.Qdec);

    r.P = ctmc_randomization(r.Qperm, q).P;

    Matrix<T> A = r.P;
    zero_offblock(A, MS);
    r.B = Matrix<T>(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) r.B(i, j) = r.P(i, j) - A(i, j);

    r.epsColMax = zero;
    for (std::size_t j = 0; j < n; ++j) {
        T cs = zero;
        for (std::size_t i = 0; i < n; ++i) cs += r.B(i, j);
        if (j == 0 || cs > r.epsColMax) r.epsColMax = cs;
    }
    r.eps = zero;
    for (std::size_t i = 0; i < n; ++i) {
        T rs = zero;
        for (std::size_t j = 0; j < n; ++j) rs += r.B(i, j);
        if (i == 0 || rs > r.eps) r.eps = rs;
    }
    r.epsRowMax = r.eps;

    // epsMAX: make each diagonal block stochastic by absorbing the row deficit
    // into its diagonal, then take its subdominant eigenvalue modulus.
    {
        Matrix<T> As = A;
        std::size_t proc = 0;
        for (const std::vector<std::size_t>& b : MS) {
            const std::size_t sz = b.size();
            for (std::size_t i = 0; i < sz; ++i) {
                T off = zero;
                for (std::size_t j = 0; j < sz; ++j)
                    if (j != i) off += As(proc + i, proc + j);
                As(proc + i, proc + i) = one - off;
            }
            proc += sz;
        }
        T maxSub = zero;
        proc = 0;
        for (const std::vector<std::size_t>& b : MS) {
            const std::size_t sz = b.size();
            if (sz > 1) {
                Matrix<T> blk(sz, sz);
                for (std::size_t i = 0; i < sz; ++i)
                    for (std::size_t j = 0; j < sz; ++j) blk(i, j) = As(proc + i, proc + j);
                const std::vector<T> mod = eig_moduli(blk);
                const T sub = mod[mod.size() - 2];
                if (sub > maxSub) maxSub = sub;
            }
            proc += sz;
        }
        r.epsMAX = (one - maxSub) / num_traits<T>::from_int(2);
    }

    // Microprobabilities, one decoupled block at a time.
    r.pmicro.assign(n, zero);
    {
        std::size_t proc = 0;
        for (const std::vector<std::size_t>& b : MS) {
            const std::size_t sz = b.size();
            Matrix<T> blk(sz, sz);
            for (std::size_t i = 0; i < sz; ++i)
                for (std::size_t j = 0; j < sz; ++j) blk(i, j) = r.Qdec(proc + i, proc + j);
            const std::vector<T> pb = ctmc_solve_reducible(blk).pi;
            for (std::size_t i = 0; i < sz; ++i) r.pmicro[proc + i] = pb[i];
            proc += sz;
        }
    }

    // Macro-state chain, weighted by the micro distributions.
    r.G = Matrix<T>(nMacro, nMacro, zero);
    std::size_t procRows = 0;
    for (std::size_t i = 0; i < nMacro; ++i) {
        std::size_t procCols = 0;
        for (std::size_t j = 0; j < nMacro; ++j) {
            if (i != j) {
                T acc = zero;
                for (std::size_t a = 0; a < MS[i].size(); ++a) {
                    T s = zero;
                    for (std::size_t b = 0; b < MS[j].size(); ++b) s += r.P(procRows + a, procCols + b);
                    acc += r.pmicro[procRows + a] * s;
                }
                r.G(i, j) = acc;
            }
            procCols += MS[j].size();
        }
        procRows += MS[i].size();
    }
    for (std::size_t i = 0; i < nMacro; ++i) {
        T rs = zero;
        for (std::size_t j = 0; j < nMacro; ++j)
            if (j != i) rs += r.G(i, j);
        r.G(i, i) = one - rs;
    }
    return r;
}

/** The rate MATLAB derives when none is supplied, (21/20) max|Qperm|. */
template <class T>
T courtois_default_rate(const Matrix<T>& Q, const std::vector<std::vector<std::size_t>>& MS) {
    const std::vector<std::size_t> v = macrostate_permutation(MS, Q.rows());
    const T m = ctmc_maxabs(submatrix(Q, v));
    if (m == num_traits<T>::from_int(0)) throw InputError("ctmc_courtois: the generator has no transitions");
    return T(m * num_traits<T>::from_rational(21, 20));
}

/** Scatters a permuted vector back to the original state ordering. */
template <class T>
std::vector<T> unpermute_states(const std::vector<T>& pperm, const std::vector<std::size_t>& v) {
    std::vector<T> p(pperm.size(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < v.size(); ++i) p[v[i]] = pperm[i];
    return p;
}

}  // namespace detail

template <class T>
struct CourtoisResult {
    std::vector<T> p;     ///< approximate stationary vector, ORIGINAL state ordering
    std::vector<std::size_t> v;  ///< macro-state-major permutation used
    Matrix<T> Qperm;      ///< Q reordered by macro-state
    Matrix<T> Qdec;       ///< block-diagonal generator of the decoupled chain
    Matrix<T> P;          ///< uniformized Qperm
    Matrix<T> B;          ///< the coupling part of P, P minus its block diagonal
    T C;                  ///< the reference's degenerate output, identically zero
    T eps;                ///< NCD index: largest ROW sum of B, ||B||_inf (MATLAB and the JAR)
    T epsRowMax;          ///< the same quantity under the name it had when only the JAR computed it
    T epsColMax;          ///< largest COLUMN sum of B: what MATLAB reported before 2026-08-15
    T epsMAX;             ///< (1 - max subdominant block eigenvalue modulus) / 2
    T q;                  ///< uniformization rate used
};

/**
 * @param Q  generator
 * @param MS macro-states, MS[i] listing the states of macro-state i; the sets
 *           must partition 0..n-1
 * @param q  uniformization rate
 */
template <class T>
CourtoisResult<T> ctmc_courtois(const Matrix<T>& Q, const std::vector<std::vector<std::size_t>>& MS,
                                const T& q) {
    static_assert(num_traits<T>::has_transcendental,
                  "ctmc_courtois requires transcendental arithmetic: epsMAX is a subdominant "
                  "eigenvalue modulus, computed by an iterative QR that stops on a tolerance and "
                  "has no rational closed form");
    const detail::CourtoisCore<T> c = detail::courtois_core(Q, MS, q);
    const std::vector<T> pMacro = dtmc_solve(c.G);

    std::vector<T> pperm(Q.rows(), num_traits<T>::from_int(0));
    std::size_t proc = 0;
    for (std::size_t i = 0; i < MS.size(); ++i) {
        for (std::size_t a = 0; a < MS[i].size(); ++a) pperm[proc + a] = pMacro[i] * c.pmicro[proc + a];
        proc += MS[i].size();
    }

    CourtoisResult<T> r;
    r.p = detail::unpermute_states(pperm, c.v);
    r.v = c.v;
    r.Qperm = c.Qperm;
    r.Qdec = c.Qdec;
    r.P = c.P;
    r.B = c.B;
    r.C = num_traits<T>::from_int(0);
    r.eps = c.eps;
    r.epsRowMax = c.epsRowMax;
    r.epsColMax = c.epsColMax;
    r.epsMAX = c.epsMAX;
    r.q = c.q;
    return r;
}

/** Overload deriving the rate as MATLAB does, q = (21/20) max|Qperm|. */
template <class T>
CourtoisResult<T> ctmc_courtois(const Matrix<T>& Q, const std::vector<std::vector<std::size_t>>& MS) {
    return ctmc_courtois(Q, MS, detail::courtois_default_rate(Q, MS));
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_COURTOIS_H
