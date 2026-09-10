/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_PASSAGE_H
#define LINE_API_MC_CTMC_PASSAGE_H

/**
 * First passage times into a target STATE SET, for Markov and semi-Markov
 * chains.
 *
 * Port of matlab/src/api/mc/ctmc_passage_*.m and smp_passage_*.m, which are the
 * reference. Source: P. G. Harrison and W. J. Knottenbelt, "Passage Time
 * Distributions in Large Markov Chains", 2002 -- Eqs. 1-3 for the Markov case,
 * Eqs. 4-8 for the semi-Markov one.
 *
 * THE IDENTITY THE WHOLE FAMILY RESTS ON. The first passage time from an
 * initial law pi0 into a target set B is PHASE-TYPE. With A the complement,
 *
 *     S     = Q(A,A)        sub-generator: the passage has not completed
 *     s0    = -S*1          exit vector, equal to Q(A,B)*1
 *     alpha = pi0(A)        UNNORMALIZED, see below
 *     atom  = sum pi0(B)
 *
 * so L(s) = alpha (sI-S)^{-1} s0 + atom and F(t) = 1 - alpha exp(St) 1. The
 * paper writes the same system as n scalar equations with L_i = 1 on B.
 *
 * ALPHA IS DELIBERATELY NOT NORMALIZED. Its mass is 1 - atom; the missing mass
 * is the ATOM AT ZERO carried by initial states already inside the target. A
 * caller that normalizes alpha and forgets the atom reports F(0) = 0 for a
 * passage that has already completed with probability atom.
 *
 * THIS IS NOT THE SPLIT `solver_ctmc_cdf.h` USES. That one is by EVENT (the
 * tagged job arriving at or departing from a station, through the filtration),
 * this one is by STATE SET. The two are complementary and must not be merged.
 *
 * ARITHMETIC: needs `expm` and transcendentals, so exact arithmetic is refused
 * by name rather than silently producing a wrong type.
 */

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <functional>
#include <limits>
#include <vector>

#include "line/api/lti/laplace_invert.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

/** The phase-type form of a first passage time. */
template <class T>
struct PassagePh {
    std::vector<T> alpha;        ///< pi0 restricted to the non-target block, UNNORMALIZED
    Matrix<T> S;                 ///< sub-generator Q(A,A)
    std::vector<T> s0;           ///< exit vector -S*1
    std::vector<std::size_t> keep;  ///< row of S -> state index of Q
    T atom;                      ///< mass of pi0 already inside the target: F(0)
};

namespace passage_detail {

inline std::vector<std::size_t> unique_target(const std::vector<std::size_t>& target,
                                              std::size_t n, const char* fn) {
    std::vector<std::size_t> t = target;
    std::sort(t.begin(), t.end());
    t.erase(std::unique(t.begin(), t.end()), t.end());
    if (t.empty())
        throw InputError(std::string(fn) +
                         ": the target state set is empty: a first passage time into no state "
                         "is undefined");
    if (t.back() >= n)
        throw InputError(std::string(fn) + ": a target state index is outside the state space");
    return t;
}

/**
 * Gaussian elimination with partial pivoting for a COMPLEX system.
 *
 * `line::solve` cannot serve here: its pivot search compares magnitudes with
 * `operator>`, which std::complex does not have, so instantiating it on
 * complex is a compile error rather than a silent wrong answer. The transform
 * routes need complex coefficients, so they carry their own solve.
 */
inline std::vector<std::complex<double>> solve_cplx(Matrix<std::complex<double>> A,
                                                    std::vector<std::complex<double>> b,
                                                    const char* fn) {
    using C = std::complex<double>;
    const std::size_t n = A.rows();
    if (A.cols() != n || b.size() != n)
        throw InputError(std::string(fn) + ": complex system is not square");
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t piv = k;
        double best = std::abs(A(k, k));
        for (std::size_t i = k + 1; i < n; ++i) {
            const double m = std::abs(A(i, k));
            if (m > best) {
                best = m;
                piv = i;
            }
        }
        if (!(best > 0.0))
            throw NumericError(std::string(fn) + ": singular transform matrix");
        if (piv != k) {
            for (std::size_t j = 0; j < n; ++j) std::swap(A(k, j), A(piv, j));
            std::swap(b[k], b[piv]);
        }
        for (std::size_t i = k + 1; i < n; ++i) {
            const C f = A(i, k) / A(k, k);
            if (f == C(0.0, 0.0)) continue;
            for (std::size_t j = k; j < n; ++j) A(i, j) -= f * A(k, j);
            b[i] -= f * b[k];
        }
    }
    std::vector<C> x(n, C(0.0, 0.0));
    for (std::size_t i = n; i-- > 0;) {
        C s = b[i];
        for (std::size_t j = i + 1; j < n; ++j) s -= A(i, j) * x[j];
        x[i] = s / A(i, i);
    }
    return x;
}

/** Backward reachability closure over the transition graph. */
template <class T>
std::vector<bool> reaches_target(const Matrix<T>& Q, const std::vector<std::size_t>& keep,
                                 const std::vector<std::size_t>& target) {
    const std::size_t n = Q.rows();
    std::vector<bool> seen(n, false);
    std::vector<std::size_t> frontier;
    for (std::size_t k : target) {
        seen[k] = true;
        frontier.push_back(k);
    }
    while (!frontier.empty()) {
        std::vector<std::size_t> next;
        for (std::size_t j : frontier)
            for (std::size_t i = 0; i < n; ++i)
                if (i != j && !seen[i] && Q(i, j) != T(0)) {
                    seen[i] = true;
                    next.push_back(i);
                }
        frontier.swap(next);
    }
    std::vector<bool> out(keep.size(), false);
    for (std::size_t a = 0; a < keep.size(); ++a) out[a] = seen[keep[a]];
    return out;
}

}  // namespace passage_detail

/**
 * Phase-type representation of the first passage time from `pi0` into `target`.
 *
 * @param pi0 empty selects the conditional stationary law on the complement
 * @param target 0-based state indices (1-based in the MATLAB reference)
 */
template <class T>
PassagePh<T> ctmc_passage_ph(const Matrix<T>& Q, const std::vector<T>& pi0,
                             const std::vector<std::size_t>& target) {
    const std::size_t n = Q.rows();
    if (Q.cols() != n) throw InputError("ctmc_passage_ph: the generator must be square");
    const std::vector<std::size_t> tgt =
        passage_detail::unique_target(target, n, "ctmc_passage_ph");

    double scale = 1.0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            scale = std::max(scale, std::abs(num_traits<T>::to_double(Q(i, j))));
    for (std::size_t i = 0; i < n; ++i) {
        T rs = T(0);
        for (std::size_t j = 0; j < n; ++j) rs += Q(i, j);
        if (std::abs(num_traits<T>::to_double(rs)) > 1e-8 * scale)
            throw InputError(
                "ctmc_passage_ph: Q is not an infinitesimal generator: its rows do not sum to "
                "zero. Pass it through ctmc_makeinfgen first");
    }

    std::vector<bool> is_target(n, false);
    for (std::size_t k : tgt) is_target[k] = true;
    PassagePh<T> out;
    for (std::size_t i = 0; i < n; ++i)
        if (!is_target[i]) out.keep.push_back(i);

    const std::size_t nA = out.keep.size();
    out.S = Matrix<T>(nA, nA);
    for (std::size_t a = 0; a < nA; ++a)
        for (std::size_t c = 0; c < nA; ++c) out.S(a, c) = Q(out.keep[a], out.keep[c]);
    out.s0.assign(nA, T(0));
    for (std::size_t a = 0; a < nA; ++a) {
        T r = T(0);
        for (std::size_t c = 0; c < nA; ++c) r += out.S(a, c);
        out.s0[a] = -r;
    }

    out.alpha.assign(nA, T(0));
    out.atom = T(0);
    if (pi0.empty()) {
        // An empty initial law selects the conditional stationary one on the
        // complement of the target set, the contract of the MATLAB reference
        // and of the python and JAR api twins.
        const std::vector<T> p = ctmc_solve(Q);
        T mass = T(0);
        for (std::size_t a = 0; a < nA; ++a) mass += p[out.keep[a]];
        if (!(num_traits<T>::to_double(mass) > 0.0))
            throw InputError(
                "ctmc_passage_ph: the stationary law puts no mass outside the target set, so "
                "there is no passage to time");
        for (std::size_t a = 0; a < nA; ++a) out.alpha[a] = p[out.keep[a]] / mass;
        return out;
    }
    if (pi0.size() != n)
        throw InputError(
            "ctmc_passage_ph: pi0 must be a distribution over the state space, one entry per "
            "state");
    for (std::size_t a = 0; a < nA; ++a) out.alpha[a] = pi0[out.keep[a]];
    for (std::size_t k : tgt) out.atom += pi0[k];
    return out;
}

/**
 * L(s) = alpha (sI-S)^{-1} s0 + atom at the (complex) points `s`. Eqs. 1-2: one
 * linear system per value of s.
 *
 * ONE SOLVE PER s, NOT PER (s,t) PAIR. The saving over a dense matrix
 * exponential is that the solves are sparse, so this route reaches chains a
 * dense expm cannot hold. It is NOT a saving in the number of time points:
 * every Abate-Whitt inverter places its nodes at s = beta/t, so a grid of T
 * points costs T*|beta| solves.
 */
template <class T>
std::vector<std::complex<double>> ctmc_passage_lst(const Matrix<T>& Q, const std::vector<T>& pi0,
                                                   const std::vector<std::size_t>& target,
                                                   const std::vector<std::complex<double>>& s) {
    const PassagePh<T> ph = ctmc_passage_ph(Q, pi0, target);
    const std::size_t nA = ph.S.rows();
    using C = std::complex<double>;
    std::vector<C> out(s.size(), C(0.0, 0.0));
    for (std::size_t is = 0; is < s.size(); ++is) {
        Matrix<C> A(nA, nA);
        for (std::size_t i = 0; i < nA; ++i)
            for (std::size_t j = 0; j < nA; ++j)
                A(i, j) = (i == j ? s[is] : C(0.0, 0.0)) -
                          C(num_traits<T>::to_double(ph.S(i, j)), 0.0);
        std::vector<C> b(nA);
        for (std::size_t i = 0; i < nA; ++i) b[i] = C(num_traits<T>::to_double(ph.s0[i]), 0.0);
        const std::vector<C> x = passage_detail::solve_cplx(A, b, "ctmc_passage_lst");
        C acc(num_traits<T>::to_double(ph.atom), 0.0);
        for (std::size_t i = 0; i < nA; ++i)
            acc += C(num_traits<T>::to_double(ph.alpha[i]), 0.0) * x[i];
        out[is] = acc;
    }
    return out;
}

/** Per-source and pi0-weighted passage moments. */
template <class T>
struct PassageMoments {
    Matrix<T> mall;      ///< (nstates x nmax), zero on the target, inf where unreachable
    std::vector<T> m;    ///< the pi0-weighted moment vector
};

/**
 * Moments of order 1..nmax of the first passage time into `target`.
 *
 * This is Eq. 3, -q_ii M_i(n) = sum_{k not in B} q_ik M_k(n) + n M_i(n-1),
 * i.e. (-S) M(n) = n M(n-1) with M(0) = 1: nmax linear solves and no transform
 * inversion at all. The equivalent closed form n! alpha (-S)^{-n} 1 is NOT how
 * it is evaluated -- forming the inverse of the sub-generator destroys the
 * sparsity the recursion preserves.
 */
template <class T>
PassageMoments<T> ctmc_passage_moments(const Matrix<T>& Q, const std::vector<T>& pi0,
                                       const std::vector<std::size_t>& target,
                                       std::size_t nmax = 1) {
    static_assert(num_traits<T>::has_transcendental,
                  "ctmc_passage_moments marks unreachable states with an infinity and therefore "
                  "requires an arithmetic that has one");
    if (nmax == 0) throw InputError("ctmc_passage_moments: nmax must be positive");
    const PassagePh<T> ph = ctmc_passage_ph(Q, pi0, target);
    const std::size_t n = Q.rows();
    const std::vector<std::size_t> tgt =
        passage_detail::unique_target(target, n, "ctmc_passage_moments");
    const std::size_t nA = ph.keep.size();

    PassageMoments<T> out;
    out.mall = Matrix<T>(n, nmax);
    out.m.assign(nmax, T(0));
    if (nA == 0) return out;

    // A state that cannot reach the target has an infinite passage time; the
    // sub-generator is singular on that block, and a solve that ignored this
    // would return a finite number instead of saying so.
    const std::vector<bool> reach = passage_detail::reaches_target(Q, ph.keep, tgt);
    bool all_reach = true;
    for (bool r : reach) all_reach = all_reach && r;

    Matrix<T> A(nA, nA);
    for (std::size_t i = 0; i < nA; ++i)
        for (std::size_t j = 0; j < nA; ++j) A(i, j) = -ph.S(i, j);

    const T inf = num_traits<T>::from_double(std::numeric_limits<double>::infinity());
    std::vector<T> x(nA, T(1));
    for (std::size_t k = 1; k <= nmax; ++k) {
        std::vector<T> rhs(nA);
        for (std::size_t i = 0; i < nA; ++i) rhs[i] = num_traits<T>::from_double(double(k)) * x[i];
        if (all_reach) {
            x = solve(A, rhs);
        } else {
            // Restrict to the reachable block: the unreachable rows are exactly
            // the singular ones, and they are reported as infinite rather than
            // regularized away.
            std::vector<std::size_t> idx;
            for (std::size_t i = 0; i < nA; ++i)
                if (reach[i]) idx.push_back(i);
            Matrix<T> Ar(idx.size(), idx.size());
            std::vector<T> br(idx.size());
            for (std::size_t a = 0; a < idx.size(); ++a) {
                for (std::size_t c = 0; c < idx.size(); ++c) Ar(a, c) = A(idx[a], idx[c]);
                br[a] = rhs[idx[a]];
            }
            const std::vector<T> xr = solve(Ar, br);
            x.assign(nA, inf);
            for (std::size_t a = 0; a < idx.size(); ++a) x[idx[a]] = xr[a];
        }
        for (std::size_t i = 0; i < nA; ++i)
            out.mall(ph.keep[i], k - 1) = reach[i] ? x[i] : inf;
        if (!all_reach)
            for (std::size_t i = 0; i < nA; ++i)
                if (!reach[i]) x[i] = T(1);  // keeps the recursion finite on the reachable block
    }

    bool unreachable_start = false;
    for (std::size_t i = 0; i < nA; ++i)
        if (!reach[i] && num_traits<T>::to_double(ph.alpha[i]) > 0.0) unreachable_start = true;
    for (std::size_t k = 0; k < nmax; ++k) {
        if (unreachable_start) {
            out.m[k] = inf;
            continue;
        }
        T acc = T(0);
        for (std::size_t i = 0; i < nA; ++i)
            if (reach[i]) acc += ph.alpha[i] * out.mall(ph.keep[i], k);
        out.m[k] = acc;
    }
    return out;
}

/**
 * Mean time to reach any state in `target` from each state of a CTMC.
 *
 * Continuous-time twin of `dtmc_hitting_time` and the first-moment special case
 * of `ctmc_passage_moments`: (-S) h = 1 on the non-target block, where
 * `dtmc_hitting_time` solves (I - P_NT) h = 1. Unreachable states give infinity.
 */
template <class T>
std::vector<T> ctmc_hitting_time(const Matrix<T>& Q, const std::vector<std::size_t>& target) {
    const std::size_t n = Q.rows();
    // mall does not depend on the initial law, so a uniform one is passed
    // rather than requiring the caller to invent one.
    std::vector<T> pi0(n, num_traits<T>::from_double(1.0 / double(n)));
    const PassageMoments<T> pm = ctmc_passage_moments(Q, pi0, target, 1);
    std::vector<T> h(n);
    for (std::size_t i = 0; i < n; ++i) h[i] = pm.mall(i, 0);
    return h;
}

/** A passage-time law on a grid. */
template <class T>
struct PassageCurve {
    std::vector<double> t;
    std::vector<double> F;
    std::vector<double> f;
    double atom = 0.0;
};

/**
 * CDF and density of the first passage time on the grid `tset`:
 * F(t) = 1 - alpha exp(St) 1 and f(t) = alpha exp(St) s0.
 *
 * `method` is "expm" (default) or "lt". The transform route exists for chains
 * whose non-target block is too large for a dense exp(St), not because it needs
 * fewer time points; on a small chain "expm" is both faster and more accurate,
 * which is why it is the default.
 */
template <class T>
PassageCurve<T> ctmc_passage_time(const Matrix<T>& Q, const std::vector<T>& pi0,
                                  const std::vector<std::size_t>& target,
                                  const std::vector<double>& tset,
                                  const std::string& method = "expm",
                                  const std::string& lti_method = "euler") {
    static_assert(num_traits<T>::has_transcendental,
                  "ctmc_passage_time takes a matrix exponential and therefore requires an "
                  "arithmetic with transcendental functions");
    const PassagePh<T> ph = ctmc_passage_ph(Q, pi0, target);
    const std::size_t nA = ph.S.rows();
    PassageCurve<T> out;
    out.t = tset;
    out.F.assign(tset.size(), 0.0);
    out.f.assign(tset.size(), 0.0);
    out.atom = num_traits<T>::to_double(ph.atom);

    if (method == "expm") {
        bool uniform = tset.size() > 2;
        const double dt = tset.size() > 1 ? tset[1] - tset[0] : 0.0;
        if (uniform && !(dt > 0.0)) uniform = false;
        for (std::size_t i = 1; uniform && i + 1 < tset.size(); ++i)
            if (std::abs((tset[i + 1] - tset[i]) - dt) > 1e-12 * std::max(1.0, std::abs(dt)))
                uniform = false;

        std::vector<double> v(nA, 0.0);
        Matrix<double> Sd(nA, nA);
        for (std::size_t i = 0; i < nA; ++i)
            for (std::size_t j = 0; j < nA; ++j) Sd(i, j) = num_traits<T>::to_double(ph.S(i, j));
        std::vector<double> s0d(nA), ad(nA);
        for (std::size_t i = 0; i < nA; ++i) {
            s0d[i] = num_traits<T>::to_double(ph.s0[i]);
            ad[i] = num_traits<T>::to_double(ph.alpha[i]);
        }

        auto step = [&](const Matrix<double>& E, std::vector<double>& w) {
            std::vector<double> z(nA, 0.0);
            for (std::size_t j = 0; j < nA; ++j) {
                double acc = 0.0;
                for (std::size_t i = 0; i < nA; ++i) acc += w[i] * E(i, j);
                z[j] = acc;
            }
            w.swap(z);
        };

        if (uniform) {
            // One exponential, then propagate: recomputing expm(S*t) at every
            // grid point is the same answer at a cost linear in the grid.
            Matrix<double> Sdt = Sd;
            for (std::size_t i = 0; i < nA; ++i)
                for (std::size_t j = 0; j < nA; ++j) Sdt(i, j) = Sd(i, j) * dt;
            const Matrix<double> E = expm(Sdt);
            Matrix<double> S0 = Sd;
            for (std::size_t i = 0; i < nA; ++i)
                for (std::size_t j = 0; j < nA; ++j) S0(i, j) = Sd(i, j) * tset[0];
            const Matrix<double> E0 = expm(S0);
            v.assign(nA, 0.0);
            for (std::size_t j = 0; j < nA; ++j) {
                double acc = 0.0;
                for (std::size_t i = 0; i < nA; ++i) acc += ad[i] * E0(i, j);
                v[j] = acc;
            }
            for (std::size_t i = 0; i < tset.size(); ++i) {
                if (i > 0) step(E, v);
                double sF = 0.0, sf = 0.0;
                for (std::size_t j = 0; j < nA; ++j) {
                    sF += v[j];
                    sf += v[j] * s0d[j];
                }
                out.F[i] = 1.0 - sF;
                out.f[i] = sf;
            }
        } else {
            for (std::size_t i = 0; i < tset.size(); ++i) {
                if (tset[i] < 0.0) continue;
                Matrix<double> St = Sd;
                for (std::size_t a = 0; a < nA; ++a)
                    for (std::size_t b = 0; b < nA; ++b) St(a, b) = Sd(a, b) * tset[i];
                const Matrix<double> E = expm(St);
                double sF = 0.0, sf = 0.0;
                for (std::size_t j = 0; j < nA; ++j) {
                    double acc = 0.0;
                    for (std::size_t a = 0; a < nA; ++a) acc += ad[a] * E(a, j);
                    sF += acc;
                    sf += acc * s0d[j];
                }
                out.F[i] = 1.0 - sF;
                out.f[i] = sf;
            }
        }
    } else if (method == "lt") {
        const double atom = out.atom;
        const lti::LaplaceFn L = [&](std::complex<double> s) {
            std::vector<std::complex<double>> sv(1, s);
            return ctmc_passage_lst(Q, pi0, target, sv)[0];
        };
        const lti::LaplaceMethod lm = lti::laplace_method(lti_method);
        out.F = lti::laplace_invert_cdf(L, tset, lm);
        const lti::LaplaceFn Ld = [&](std::complex<double> s) { return L(s) - atom; };
        out.f = lti::laplace_invert_pdf(Ld, tset, lm);
    } else {
        throw InputError("ctmc_passage_time: unknown method '" + method +
                         "', expected expm or lt");
    }

    for (std::size_t i = 0; i < out.F.size(); ++i) {
        out.F[i] = std::min(1.0, std::max(0.0, out.F[i]));
        out.f[i] = std::max(0.0, out.f[i]);
    }
    return out;
}

// ---------------------------------------------------------------------------
// Semi-Markov chains (Sec. 3)
// ---------------------------------------------------------------------------

/**
 * Moments of the semi-Markov first passage time from the per-state holding
 * moments m_i(r), Eq. 7 with the u_i(r) recurrence of Eq. 8:
 *
 *     u_i(r) = -sum_{j=1..r} C(r,j) m_i(j) u_i(r-j),   u_i(0) = 1,
 *
 * which are the derivatives at the origin of 1/h*_i(s). Cheaper than Eq. 6
 * because it needs no per-pair moments.
 *
 * @param hmom (nstates x nmax): hmom(i,r-1) is the r-th moment of the sojourn
 *             in state i
 */
template <class T>
PassageMoments<T> smp_passage_moments(const Matrix<T>& P, const Matrix<T>& hmom,
                                      const std::vector<T>& pi0,
                                      const std::vector<std::size_t>& target,
                                      std::size_t nmax = 1) {
    const std::size_t n = P.rows();
    if (P.cols() != n)
        throw InputError("smp_passage_moments: the embedded transition matrix must be square");
    if (nmax == 0) throw InputError("smp_passage_moments: nmax must be positive");
    for (std::size_t i = 0; i < n; ++i) {
        T rs = T(0);
        for (std::size_t j = 0; j < n; ++j) rs += P(i, j);
        if (std::abs(num_traits<T>::to_double(rs) - 1.0) > 1e-8)
            throw InputError(
                "smp_passage_moments: the embedded transition matrix rows must sum to one");
    }
    if (hmom.rows() != n)
        throw InputError("smp_passage_moments: hmom must carry one row per state");
    if (hmom.cols() < nmax)
        throw InputError(
            "smp_passage_moments: hmom must carry at least nmax holding-time moments per state");
    const std::vector<std::size_t> tgt =
        passage_detail::unique_target(target, n, "smp_passage_moments");

    std::vector<bool> is_target(n, false);
    for (std::size_t k : tgt) is_target[k] = true;
    std::vector<std::size_t> A;
    for (std::size_t i = 0; i < n; ++i)
        if (!is_target[i]) A.push_back(i);
    const std::size_t nA = A.size();

    PassageMoments<T> out;
    out.mall = Matrix<T>(n, nmax);
    out.m.assign(nmax, T(0));
    if (nA == 0) return out;

    // Eq. 8, per state.
    Matrix<T> u(nA, nmax);
    for (std::size_t r = 1; r <= nmax; ++r) {
        for (std::size_t a = 0; a < nA; ++a) {
            T acc = T(0);
            for (std::size_t j = 1; j <= r; ++j) {
                const T base = (r - j == 0) ? T(1) : u(a, r - j - 1);
                double c = 1.0;
                for (std::size_t q = 0; q < j; ++q)
                    c = c * double(r - q) / double(q + 1);
                acc += num_traits<T>::from_double(c) * hmom(A[a], j - 1) * base;
            }
            u(a, r - 1) = -acc;
        }
    }

    Matrix<T> IPAA(nA, nA);
    for (std::size_t a = 0; a < nA; ++a)
        for (std::size_t c = 0; c < nA; ++c)
            IPAA(a, c) = (a == c ? T(1) : T(0)) - P(A[a], A[c]);

    Matrix<T> M(nA, nmax);
    for (std::size_t q = 1; q <= nmax; ++q) {
        std::vector<T> b(nA, T(0));
        for (std::size_t r = 1; r <= q; ++r) {
            double c = 1.0;
            for (std::size_t k = 0; k < r; ++k) c = c * double(q - k) / double(k + 1);
            for (std::size_t a = 0; a < nA; ++a) {
                const T base = (r < q) ? M(a, q - r - 1) : T(1);
                b[a] -= num_traits<T>::from_double(c) * u(a, r - 1) * base;
            }
        }
        const std::vector<T> x = solve(IPAA, b);
        for (std::size_t a = 0; a < nA; ++a) M(a, q - 1) = x[a];
    }

    for (std::size_t a = 0; a < nA; ++a)
        for (std::size_t q = 0; q < nmax; ++q) out.mall(A[a], q) = M(a, q);
    if (pi0.size() == n)
        for (std::size_t q = 0; q < nmax; ++q) {
            T acc = T(0);
            for (std::size_t i = 0; i < n; ++i) acc += pi0[i] * out.mall(i, q);
            out.m[q] = acc;
        }
    return out;
}

/**
 * L(s) of the semi-Markov first passage time, Eqs. 4-5:
 *
 *     L_i(s) = sum_{k not in B} r*_ik(s) L_k(s) + sum_{k in B} r*_ik(s),
 *
 * so (I - R*_AA(s)) L_A(s) = R*_AB(s) 1, one linear system per value of s.
 *
 * @param hlst per-state sojourn transforms h*_i(s), so r*_ik(s) = P(i,k) h*_i(s)
 *             and the complex numbers stay on the DIAGONAL of the system
 */
template <class T>
std::vector<std::complex<double>> smp_passage_lst(
    const Matrix<T>& P, const std::vector<std::function<std::complex<double>(std::complex<double>)>>& hlst,
    const std::vector<T>& pi0, const std::vector<std::size_t>& target,
    const std::vector<std::complex<double>>& s) {
    const std::size_t n = P.rows();
    const std::vector<std::size_t> tgt = passage_detail::unique_target(target, n, "smp_passage_lst");
    if (hlst.size() != n)
        throw InputError("smp_passage_lst: hlst must carry one transform per state");
    std::vector<bool> is_target(n, false);
    for (std::size_t k : tgt) is_target[k] = true;
    std::vector<std::size_t> A;
    for (std::size_t i = 0; i < n; ++i)
        if (!is_target[i]) A.push_back(i);
    const std::size_t nA = A.size();

    using C = std::complex<double>;
    double atom = 0.0;
    if (pi0.size() == n)
        for (std::size_t k : tgt) atom += num_traits<T>::to_double(pi0[k]);

    std::vector<C> out(s.size(), C(0.0, 0.0));
    for (std::size_t is = 0; is < s.size(); ++is) {
        std::vector<C> h(nA);
        for (std::size_t a = 0; a < nA; ++a) h[a] = hlst[A[a]](s[is]);
        Matrix<C> M(nA, nA);
        std::vector<C> b(nA, C(0.0, 0.0));
        for (std::size_t a = 0; a < nA; ++a) {
            for (std::size_t c = 0; c < nA; ++c)
                M(a, c) = (a == c ? C(1.0, 0.0) : C(0.0, 0.0)) -
                          h[a] * C(num_traits<T>::to_double(P(A[a], A[c])), 0.0);
            double pb = 0.0;
            for (std::size_t k : tgt) pb += num_traits<T>::to_double(P(A[a], k));
            b[a] = h[a] * C(pb, 0.0);
        }
        const std::vector<C> x = passage_detail::solve_cplx(M, b, "smp_passage_lst");
        C acc(atom, 0.0);
        for (std::size_t a = 0; a < nA; ++a)
            acc += C(pi0.size() == n ? num_traits<T>::to_double(pi0[A[a]]) : 1.0 / double(n), 0.0) *
                   x[a];
        out[is] = acc;
    }
    return out;
}

/**
 * CDF and density of the semi-Markov first passage time, by inverting
 * `smp_passage_lst` through api/lti.
 *
 * There is no matrix-exponential route here: a semi-Markov chain has no
 * generator to exponentiate, which is exactly the case uniformization does not
 * reach and the transform does.
 *
 * `lti_method` defaults to "euler" RATHER THAN "weeks". Semi-Markov passage
 * densities are the case Sec. 4.2 singles out as slow-converging for a Laguerre
 * series, and `laplace_weeks_scaling` then refuses by name rather than
 * returning noise.
 */
template <class T>
PassageCurve<T> smp_passage_time(
    const Matrix<T>& P, const std::vector<std::function<std::complex<double>(std::complex<double>)>>& hlst,
    const std::vector<T>& pi0, const std::vector<std::size_t>& target,
    const std::vector<double>& tset, const std::string& lti_method = "euler") {
    const std::size_t n = P.rows();
    const std::vector<std::size_t> tgt =
        passage_detail::unique_target(target, n, "smp_passage_time");
    double atom = 0.0;
    if (pi0.size() == n)
        for (std::size_t k : tgt) atom += num_traits<T>::to_double(pi0[k]);

    const lti::LaplaceFn L = [&](std::complex<double> s) {
        std::vector<std::complex<double>> sv(1, s);
        return smp_passage_lst(P, hlst, pi0, target, sv)[0];
    };
    const lti::LaplaceMethod lm = lti::laplace_method(lti_method);
    PassageCurve<T> out;
    out.t = tset;
    out.atom = atom;
    out.F = lti::laplace_invert_cdf(L, tset, lm);
    const lti::LaplaceFn Ld = [&](std::complex<double> s) { return L(s) - atom; };
    out.f = lti::laplace_invert_pdf(Ld, tset, lm);
    return out;
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_PASSAGE_H
