/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_MARIE_H
#define LINE_API_PFQN_MARIE_H

/**
 * Marie's iterative aggregation-decomposition for closed networks with FCFS
 * general (Coxian) service.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_marie.m. There is no JAR
 * counterpart, so MATLAB is the only reference.
 *
 * SINGLE CLASS (R = 1). Each station's service is given a Coxian phase
 * representation matching its mean L(i) and its squared coefficient of
 * variation scv(i). The aggregate model is the exact load-dependent product
 * form pfqn_mvald driven by a multiplier lattice mu(i,n), from which the
 * marginal queue-length distribution P_i(n) at the full population is read off.
 * Flow balance across the n <-> n+1 cut of that birth-death marginal gives the
 * complementary arrival rate seen by station i,
 *
 *   lambda_i(n) = (mu(i,n+1)/L(i)) P_i(n+1) / P_i(n),   n = 0, ..., N-1,
 *
 * and the lambda(n)/Cox/1(-m) isolation chain -- states (n,k) with n present
 * and the head job in phase k, phase rates scaled by min(n,m) for m servers --
 * is solved for its stationary distribution. Its conditional departure rate
 *
 *   mu_i(n) = sum_k p(n,k) rate_k phi_k min(n,m) / sum_k p(n,k)
 *
 * is converted back to a multiplier (multiplied by L(i)) and fed to the next
 * aggregate solve. The iteration stops when max|mu_new - mu| < tol. For
 * exponential service (scv == 1) the isolation chain is the M/M/1(-m) queue,
 * mu_i(n) = min(n,m)/L(i), so the multiplier lattice is the initial one, the
 * first iteration already reproduces exact product form and the loop exits on
 * the second pass.
 *
 * MULTIPLE CLASSES (R > 1). Exponential and class-independent demands at every
 * station is genuine BCMP FCFS and is dispatched to exact pfqn_mva. Otherwise
 * the aggregate is a Schweitzer-style multiclass AMVA in which the class-r
 * demand at station i is divided by a class-dependent scaling
 * beta_{i,r}(nvec) = muCox_{i,r}(nvec) / muExp_{i,r}(nvec), the ratio of the
 * conditional class-r throughput of a multiclass Cox/1 FCFS isolation chain to
 * that of the same chain with exponential service of the same means. beta == 1
 * therefore recovers plain FCFS AMVA and carries only the non-exponential
 * correction. The isolation chain is fed the aggregate per-class throughput
 * (Baynat-Dallery isolation) and the outer loop iterates to a fixed point on X.
 * beta is tabulated on the integer population box and read at the real-valued
 * arrival-instant populations by multilinear interpolation; the reference
 * guards a non-positive or non-finite ratio by falling back to 1 and clamps the
 * result to [1e-3, 1e3], and both are reproduced.
 *
 * COXIAN FIT. Coxian.fitMeanAndSCV is closed form and is reproduced here as
 * marie_cox_fit, with tol = GlobalConstants.CoarseTol = 1e-3 (matlab/lineStart.m):
 *
 *   |scv - 1| <= tol         exponential, n = 1, mu = [1/mean], phi = [1]
 *   0.5 + tol < scv < 1 - tol  hypoexponential, n = 2,
 *                              mu = 2/mean/(1 +- sqrt(2 scv - 1)), phi = [0,1]
 *   scv <= 0.5 + tol         Erlang, n = ceil(1/scv), mu = (n/mean) 1, phi = e_n
 *   scv > 1 + tol            Coxian-2, mu = [2/mean, 1/(scv mean)],
 *                              phi = [1 - 1/(2 scv), 1]
 *
 * The exponential, hypoexponential and Coxian-2 branches match the requested
 * mean and scv exactly. The Erlang branch does NOT: ceil(1/scv) is an integer,
 * so the fitted scv is 1/ceil(1/scv) <= scv, with equality only when 1/scv is
 * an integer. That is the reference's behaviour and is reproduced, not
 * improved. The phase count also goes through a double-valued ceiling, which is
 * the one place the fit is not a pure field computation.
 *
 * DIVERGENCE from the reference, deliberate and documented. MATLAB solves the
 * isolation chains as the overdetermined least-squares system
 * [Q'; ones] p = [0; 1] via backslash. This port calls mc::ctmc_solve, which
 * replaces one balance equation by the normalization and solves the resulting
 * square system, and which additionally splits a reducible generator into its
 * weakly connected components. The two are the same computation whenever the
 * isolation chain is irreducible, which is the case whenever every
 * complementary arrival rate lambda_i(n) is positive -- the only regime in
 * which the reference's own least-squares answer is a probability vector at
 * all.
 *
 * Arithmetic: INEXACT BY CONSTRUCTION, and additionally gated on
 * has_transcendental. The method is a fixed-point decomposition stopped on a
 * tolerance, the Coxian fit of an scv in (0.5, 1) needs a square root that has
 * no exact rational counterpart, and the multiclass path reports throughputs
 * from an approximate MVA. It is therefore unavailable at T = Rational.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/mc/ctmc_solve.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/api/pfqn/pfqn_mvams.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Coxian phase representation: phase rates and per-phase completion probabilities. */
template <class T>
struct MarieCoxFit {
    std::vector<T> mu;   ///< phase rates
    std::vector<T> phi;  ///< completion probability out of each phase, phi.back() == 1
};

/**
 * Closed-form Coxian fit of a mean and an SCV (matlab/src/lang/processes/Coxian.m,
 * fitMeanAndSCV), with the branch thresholds at CoarseTol = 1e-3.
 *
 * @param mean strictly positive mean
 * @param scv  strictly positive squared coefficient of variation
 */
template <class T>
MarieCoxFit<T> marie_cox_fit(const T& mean, const T& scv) {
    static_assert(num_traits<T>::has_transcendental,
                  "marie_cox_fit requires transcendental arithmetic: the hypoexponential branch "
                  "matches the second moment through a square root");
    using std::sqrt;

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T half = num_traits<T>::from_rational(1, 2);
    const T coarse = num_traits<T>::from_double(1e-3);  // GlobalConstants.CoarseTol

    if (!(mean > zero)) throw InputError("marie_cox_fit: the mean must be strictly positive");
    if (!(scv > zero)) throw InputError("marie_cox_fit: the SCV must be strictly positive");

    MarieCoxFit<T> f;
    if (scv >= T(one - coarse) && scv <= T(one + coarse)) {
        // Exponential.
        f.mu.push_back(T(one / mean));
        f.phi.push_back(one);
    } else if (scv > T(half + coarse) && scv < T(one - coarse)) {
        // Hypoexponential: two phases in series, neither completing early.
        const T d = T(sqrt(T(two * scv - one)));
        f.mu.push_back(T(two / mean / T(one + d)));
        f.mu.push_back(T(two / mean / T(one - d)));
        f.phi.push_back(zero);
        f.phi.push_back(one);
    } else if (scv <= T(half + coarse)) {
        // Erlang-n with n = ceil(1/scv); the fitted SCV is 1/n, not scv.
        const double inv = num_traits<T>::to_double(T(one / scv));
        const long n = static_cast<long>(std::ceil(inv));
        if (n < 1 || n > 100000)
            throw InputError("marie_cox_fit: the Erlang branch needs an unreasonable phase count");
        const T rate = num_traits<T>::from_int(n) / mean;
        f.mu.assign(static_cast<std::size_t>(n), rate);
        f.phi.assign(static_cast<std::size_t>(n), zero);
        f.phi.back() = one;
    } else {
        // Coxian-2, the hyperexponential rewritten in Coxian form.
        const T mu1 = T(two / mean);
        const T mu2 = T(mu1 / T(two * scv));
        f.mu.push_back(mu1);
        f.mu.push_back(mu2);
        f.phi.push_back(T(one - mu2 / mu1));
        f.phi.push_back(one);
    }
    return f;
}

namespace detail {

/**
 * Stationary analysis of a lambda(n)/Cox/1(-m) queue in isolation, returning the
 * conditional departure rate mu(n) given n present, n = 1, ..., N.
 *
 * @param lam  (N) arrival rate with n present, lam[n] for n = 0, ..., N-1
 * @param rate (P) Coxian phase rates
 * @param phi  (P) per-phase completion probabilities
 * @param m    server count; the phase rate at population n is scaled by min(n,m)
 * @param N population of the isolated station
 */
template <class T>
std::vector<T> marie_isol_condtput(const std::vector<T>& lam, const std::vector<T>& rate,
                                   const std::vector<T>& phi, int N, int m) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t P = rate.size();
    const std::size_t S = 1 + static_cast<std::size_t>(N) * P;

    // State layout: 0 = empty; (n,k) -> 1 + (n-1)P + (k-1), n = 1..N, k = 1..P.
    const auto idx = [P](int n, std::size_t k) {
        return 1 + static_cast<std::size_t>(n - 1) * P + k;
    };

    Matrix<T> Gq(S, S, zero);
    Gq(0, idx(1, 0)) += lam[0];
    for (int n = 1; n <= N; ++n) {
        const T sc = num_traits<T>::from_int(n < m ? n : m);
        for (std::size_t k = 0; k < P; ++k) {
            const std::size_t r = idx(n, k);
            if (n < N) Gq(r, idx(n + 1, k)) += lam[static_cast<std::size_t>(n)];
            const T compl_ = rate[k] * phi[k] * sc;
            const T adv = rate[k] * T(one - phi[k]) * sc;
            if (adv > zero && k + 1 < P) Gq(r, idx(n, k + 1)) += adv;
            if (compl_ > zero) {
                if (n > 1)
                    Gq(r, idx(n - 1, 0)) += compl_;
                else
                    Gq(r, 0) += compl_;
            }
        }
    }

    const std::vector<T> p = mc::ctmc_solve(Gq);

    std::vector<T> muvec(static_cast<std::size_t>(N), zero);
    for (int n = 1; n <= N; ++n) {
        const T sc = num_traits<T>::from_int(n < m ? n : m);
        T Pn = zero, dep = zero;
        for (std::size_t k = 0; k < P; ++k) {
            const T pk = p[idx(n, k)];
            Pn += pk;
            dep += pk * rate[k] * phi[k] * sc;
        }
        if (Pn > zero) {
            muvec[static_cast<std::size_t>(n - 1)] = dep / Pn;
        } else {
            // Fallback: the exponential-equivalent rate of the phase chain.
            T mean = zero;
            for (std::size_t k = 0; k < P; ++k) mean += one / rate[k];
            muvec[static_cast<std::size_t>(n - 1)] = sc / mean;
        }
    }
    return muvec;
}

}  // namespace detail

/**
 * Class-dependent scaling of one station, tabulated on the integer population
 * box. Both the Coxian and the exponential conditional throughputs are kept, so
 * that the ratio is formed AFTER interpolation exactly as the reference does.
 * An empty table means beta == 1, the initial value of the outer iteration.
 */
template <class T>
struct MarieCdScaling {
    std::vector<int> N;                 ///< box bounds; the lattice is [0..N]
    std::vector<std::vector<T>> muCox;  ///< muCox[r][lin], lin the row-major box index
    std::vector<std::vector<T>> muExp;  ///< muExp[r][lin]

    bool identity() const { return muCox.empty(); }
};

namespace detail {

/** Row-major linear index over the box [0..N]. */
inline std::size_t marie_box_index(const std::vector<int>& n, const std::vector<int>& N) {
    std::size_t lin = 0;
    for (std::size_t d = 0; d < N.size(); ++d)
        lin = lin * static_cast<std::size_t>(N[d] + 1) + static_cast<std::size_t>(n[d]);
    return lin;
}

inline std::size_t marie_box_size(const std::vector<int>& N) {
    std::size_t sz = 1;
    for (std::size_t d = 0; d < N.size(); ++d) sz *= static_cast<std::size_t>(N[d] + 1);
    return sz;
}

/**
 * Multilinear interpolation of a table over the box [0..N] at a real point,
 * clamped to the box (ndlininterp in the reference).
 */
template <class T>
T marie_ndlininterp(const std::vector<T>& A, const std::vector<T>& x, const std::vector<int>& N) {
    const std::size_t R = N.size();
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    std::vector<int> lo(R), hi(R);
    std::vector<T> fr(R);
    for (std::size_t d = 0; d < R; ++d) {
        T xd = x[d];
        if (xd < zero) xd = zero;
        const T ub = num_traits<T>::from_int(N[d]);
        if (xd > ub) xd = ub;
        const double f = std::floor(num_traits<T>::to_double(xd));
        lo[d] = static_cast<int>(f);
        if (lo[d] > N[d]) lo[d] = N[d];
        hi[d] = lo[d] + 1 < N[d] ? lo[d] + 1 : N[d];
        fr[d] = T(xd - num_traits<T>::from_int(lo[d]));
    }
    T v = zero;
    std::vector<int> sub(R);
    const unsigned long corners = 1ul << R;
    for (unsigned long mask = 0; mask < corners; ++mask) {
        T w = one;
        for (std::size_t d = 0; d < R; ++d) {
            if ((mask >> d) & 1ul) {
                sub[d] = hi[d];
                w *= fr[d];
            } else {
                sub[d] = lo[d];
                w *= T(one - fr[d]);
            }
        }
        if (w == zero) continue;
        v += w * A[marie_box_index(sub, N)];
    }
    return v;
}

}  // namespace detail

/**
 * Evaluate a class-dependent scaling at a real-valued population vector
 * (cdscale_eval in the reference): the interpolated ratio, guarded against a
 * non-positive or non-finite value and clamped to [1e-3, 1e3].
 */
template <class T>
std::vector<T> marie_cd_eval(const MarieCdScaling<T>& cd, const std::vector<T>& nv) {
    const std::size_t R = nv.size();
    const T one = num_traits<T>::from_int(1);
    std::vector<T> be(R, one);
    if (cd.identity()) return be;
    const T zero = num_traits<T>::from_int(0);
    const T lo = num_traits<T>::from_double(1e-3);
    const T hi = num_traits<T>::from_double(1e3);
    const T inf = num_traits<T>::from_double(std::numeric_limits<double>::infinity());
    for (std::size_t r = 0; r < R; ++r) {
        const T num = detail::marie_ndlininterp(cd.muCox[r], nv, cd.N);
        const T den = detail::marie_ndlininterp(cd.muExp[r], nv, cd.N);
        const bool okNum = num == num && num > zero && num < inf;
        const bool okDen = den == den && den > zero && den < inf;
        be[r] = (okNum && okDen) ? T(num / den) : one;
        if (be[r] < lo) be[r] = lo;
        if (be[r] > hi) be[r] = hi;
    }
    return be;
}

/** Result of pfqn_marie, mirroring the six MATLAB outputs. */
template <class T>
struct MarieResult {
    std::vector<T> X;  ///< (R) per-class throughput
    Matrix<T> Q;       ///< (M x R) mean queue length
    Matrix<T> U;       ///< (M x R) utilization
    /**
     * Single class: (1 x 1), the CYCLE time returned by pfqn_mvald. Multiple
     * classes: (M x R) per-station residence times. The shapes differ because
     * the reference's two paths return different quantities under the same
     * name; that is reproduced rather than harmonized.
     */
    Matrix<T> C;
    int it;      ///< iterations performed
    Matrix<T> mu;  ///< single class: (M x N) converged multiplier lattice; empty for R > 1
    std::vector<MarieCdScaling<T>> cds;  ///< R > 1: converged per-station cd scalings
};

namespace detail {

/**
 * Multiclass Schweitzer AMVA with a class-dependent service-rate scaling
 * (amva_qd in the reference). Gauss-Seidel in the class index, so the queue
 * lengths written for class r are visible to class r+1 within the same sweep;
 * that update order is part of the fixed point and is reproduced.
 */
template <class T>
void marie_amva_qd(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z,
                   const std::vector<MarieCdScaling<T>>& cds, std::vector<T>& X, Matrix<T>& Q,
                   Matrix<T>& U, Matrix<T>& C) {
    const std::size_t M = L.rows();
    const std::size_t R = L.cols();
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T tol = num_traits<T>::from_double(1e-9);

    Q = Matrix<T>(M, R, zero);
    const T Mt = num_traits<T>::from_int(static_cast<long>(M > 1 ? M : 1));
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) Q(i, r) = num_traits<T>::from_int(N[r]) / Mt;

    Matrix<T> W(M, R, zero);
    X.assign(R, zero);
    U = Matrix<T>(M, R, zero);
    Matrix<T> Qprev(M, R, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) Qprev(i, r) = T(Q(i, r) + one);

    std::vector<T> nv(R), Leff(R);
    int it = 0;
    for (;;) {
        T delta = zero;
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) {
                const T d = num_abs(T(Q(i, r) - Qprev(i, r)));
                if (d > delta) delta = d;
            }
        if (!(delta > tol) || it >= 5000) break;
        ++it;
        Qprev = Q;
        for (std::size_t r = 0; r < R; ++r) {
            for (std::size_t i = 0; i < M; ++i) {
                for (std::size_t s = 0; s < R; ++s) nv[s] = Q(i, s);
                if (N[r] > 0)
                    nv[r] = Q(i, r) * num_traits<T>::from_int(N[r] - 1) /
                            num_traits<T>::from_int(N[r]);
                const std::vector<T> be = marie_cd_eval(cds[i], nv);
                T w = zero;
                for (std::size_t s = 0; s < R; ++s) {
                    Leff[s] = L(i, s) / be[s];
                    w += Leff[s] * nv[s];
                }
                W(i, r) = Leff[r] + w;
            }
            T denom = Z[r];
            for (std::size_t i = 0; i < M; ++i) denom += W(i, r);
            X[r] = denom > zero ? T(num_traits<T>::from_int(N[r]) / denom) : zero;
            for (std::size_t i = 0; i < M; ++i) Q(i, r) = X[r] * W(i, r);
        }
    }
    C = W;
    for (std::size_t r = 0; r < R; ++r)
        for (std::size_t i = 0; i < M; ++i) U(i, r) = X[r] * L(i, r);
}

/**
 * Stationary analysis of a multiclass lambda_r/Cox/1 FCFS queue in isolation
 * over the population box [0..N] (isol_mc in the reference). The head-of-line
 * job is tracked as (class, phase); on a departure the next head class is drawn
 * in random order, with probability n_c / sum(n). Returns, per class, the
 * conditional class-r throughput tabulated on the box.
 */
template <class T>
std::vector<std::vector<T>> marie_isol_mc(const std::vector<T>& lam,
                                          const std::vector<MarieCoxFit<T>>& fits,
                                          const std::vector<int>& N) {
    const std::size_t R = N.size();
    const T zero = num_traits<T>::from_int(0);
    const std::size_t npops = marie_box_size(N);

    // Enumerate the states: 0 is the empty station, then (pop, head class, phase).
    struct StateId {
        std::size_t pop;
        std::size_t cls;
        std::size_t phase;
    };
    std::vector<StateId> ids;
    ids.push_back(StateId{marie_box_index(std::vector<int>(R, 0), N), 0, 0});
    // id[pop][c][k] with a flat index; -1 where the state does not exist.
    std::vector<std::vector<long>> byPop(npops);
    std::vector<std::vector<int>> popVec(npops, std::vector<int>(R, 0));
    {
        std::vector<int> n(R, 0);
        for (std::size_t p = 0; p < npops; ++p) {
            // Row-major enumeration matching marie_box_index.
            std::size_t rem = p;
            for (std::size_t d = R; d-- > 0;) {
                const std::size_t w = static_cast<std::size_t>(N[d] + 1);
                n[d] = static_cast<int>(rem % w);
                rem /= w;
            }
            popVec[p] = n;
        }
    }
    for (std::size_t p = 0; p < npops; ++p) {
        int tot = 0;
        for (std::size_t r = 0; r < R; ++r) tot += popVec[p][r];
        if (tot == 0) continue;
        std::vector<long> slot;
        for (std::size_t c = 0; c < R; ++c) {
            const std::size_t P = fits[c].mu.size();
            for (std::size_t k = 0; k < P; ++k) {
                if (popVec[p][c] > 0) {
                    slot.push_back(static_cast<long>(ids.size()));
                    ids.push_back(StateId{p, c, k});
                } else {
                    slot.push_back(-1);
                }
            }
        }
        byPop[p] = slot;
    }
    std::vector<std::size_t> phaseOff(R + 1, 0);
    for (std::size_t c = 0; c < R; ++c) phaseOff[c + 1] = phaseOff[c] + fits[c].mu.size();

    const auto getid = [&](std::size_t p, std::size_t c, std::size_t k) -> std::size_t {
        const long v = byPop[p][phaseOff[c] + k];
        if (v < 0) throw NumericError("pfqn_marie: isolation state does not exist");
        return static_cast<std::size_t>(v);
    };

    const std::size_t S = ids.size();
    Matrix<T> Gq(S, S, zero);
    std::vector<int> nn(R);
    for (std::size_t s = 0; s < S; ++s) {
        const StateId& st = ids[s];
        const bool empty = (s == 0);
        const std::vector<int>& nvec = empty ? popVec[ids[0].pop] : popVec[st.pop];
        for (std::size_t r = 0; r < R; ++r) {
            if (nvec[r] < N[r] && lam[r] > zero) {
                nn = nvec;
                nn[r] += 1;
                const std::size_t pn = marie_box_index(nn, N);
                if (empty)
                    Gq(s, getid(pn, r, 0)) += lam[r];
                else
                    Gq(s, getid(pn, st.cls, st.phase)) += lam[r];
            }
        }
        if (empty) continue;
        const std::size_t c = st.cls, k = st.phase;
        const T rate = fits[c].mu[k];
        const T compl_ = rate * fits[c].phi[k];
        const T adv = rate * T(num_traits<T>::from_int(1) - fits[c].phi[k]);
        if (adv > zero && k + 1 < fits[c].mu.size()) Gq(s, getid(st.pop, c, k + 1)) += adv;
        if (compl_ > zero) {
            nn = nvec;
            nn[c] -= 1;
            int tot = 0;
            for (std::size_t r = 0; r < R; ++r) tot += nn[r];
            if (tot == 0) {
                Gq(s, 0) += compl_;
            } else {
                const std::size_t pn = marie_box_index(nn, N);
                for (std::size_t cp = 0; cp < R; ++cp)
                    if (nn[cp] > 0)
                        Gq(s, getid(pn, cp, 0)) +=
                            compl_ * num_traits<T>::from_int(nn[cp]) / num_traits<T>::from_int(tot);
            }
        }
    }

    const std::vector<T> p = mc::ctmc_solve(Gq);

    std::vector<T> Ppop(npops, zero);
    std::vector<std::vector<T>> dep(R, std::vector<T>(npops, zero));
    for (std::size_t s = 1; s < S; ++s) {
        const StateId& st = ids[s];
        Ppop[st.pop] += p[s];
        dep[st.cls][st.pop] += p[s] * fits[st.cls].mu[st.phase] * fits[st.cls].phi[st.phase];
    }
    std::vector<std::vector<T>> mumat(R, std::vector<T>(npops, zero));
    for (std::size_t r = 0; r < R; ++r)
        for (std::size_t q = 0; q < npops; ++q)
            if (Ppop[q] > zero) mumat[r][q] = dep[r][q] / Ppop[q];
    return mumat;
}

/** Multiclass path of pfqn_marie (marie_multi in the reference). */
template <class T>
MarieResult<T> marie_multi(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z,
                           const Matrix<T>& scv, double tol, int maxiter) {
    const std::size_t M = L.rows();
    const std::size_t R = L.cols();
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    MarieResult<T> res;
    res.it = 0;

    // Exact product-form dispatch: exponential AND class-independent demands.
    bool isPF = true;
    for (std::size_t i = 0; i < M && isPF; ++i)
        for (std::size_t r = 0; r < R; ++r)
            if (scv(i, r) != one) {
                isPF = false;
                break;
            }
    if (isPF) {
        const T eps = num_traits<T>::from_double(1e-12);
        for (std::size_t i = 0; i < M && isPF; ++i) {
            T lo = L(i, 0), hi = L(i, 0);
            for (std::size_t r = 1; r < R; ++r) {
                if (L(i, r) < lo) lo = L(i, r);
                if (L(i, r) > hi) hi = L(i, r);
            }
            if (T(hi - lo) > eps) isPF = false;
        }
    }
    if (isPF) {
        Matrix<T> Zmat(1, R);
        for (std::size_t r = 0; r < R; ++r) Zmat(0, r) = Z[r];
        const MvaResult<T> m = pfqn_mva(L, N, Zmat);
        res.X = m.XN;
        res.Q = m.QN;
        res.U = m.UN;
        res.C = m.CN;
        return res;
    }

    // Coxian representation, plus the exponential reference with the same means.
    std::vector<std::vector<MarieCoxFit<T>>> phCox(M, std::vector<MarieCoxFit<T>>(R));
    std::vector<std::vector<MarieCoxFit<T>>> phExp(M, std::vector<MarieCoxFit<T>>(R));
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) {
            phCox[i][r] = marie_cox_fit(L(i, r), scv(i, r));
            MarieCoxFit<T> e;
            e.mu.push_back(T(one / L(i, r)));
            e.phi.push_back(one);
            phExp[i][r] = e;
        }

    res.cds.assign(M, MarieCdScaling<T>());
    std::vector<T> Xprev(R, num_traits<T>::from_double(std::numeric_limits<double>::infinity()));
    const T tolT = num_traits<T>::from_double(tol);
    res.X.assign(R, zero);
    res.Q = Matrix<T>(M, R, zero);
    res.U = Matrix<T>(M, R, zero);
    res.C = Matrix<T>(M, R, zero);
    while (res.it < maxiter) {
        ++res.it;
        marie_amva_qd(L, N, Z, res.cds, res.X, res.Q, res.U, res.C);
        for (std::size_t i = 0; i < M; ++i) {
            MarieCdScaling<T> cd;
            cd.N = N;
            cd.muCox = marie_isol_mc(res.X, phCox[i], N);
            cd.muExp = marie_isol_mc(res.X, phExp[i], N);
            res.cds[i] = cd;
        }
        T delta = zero;
        for (std::size_t r = 0; r < R; ++r) {
            const T d = num_abs(T(res.X[r] - Xprev[r]));
            if (d > delta) delta = d;
        }
        if (delta < tolT) break;
        Xprev = res.X;
    }
    return res;
}

}  // namespace detail

/**
 * Marie's method for a closed network with FCFS Coxian service.
 *
 * @param L        (M x R) service demands, every entry strictly positive
 * @param N        (R) closed population vector
 * @param Z        (R) aggregated think times, one per class; may be empty
 * @param scv      (M x R) per-station per-class squared coefficients of
 *                 variation; empty for all-exponential service
 * @param tol      convergence tolerance (reference default 1e-8)
 * @param maxiter  iteration cap (reference default 1000)
 * @param nservers (M) server counts, single class only; empty for all single
 *                 server. The reference does not support a multiserver
 *                 multiclass isolation chain and neither does this port.
 */
template <class T>
MarieResult<T> pfqn_marie(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z,
                          const Matrix<T>& scv, double tol, int maxiter,
                          const std::vector<int>& nservers) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_marie requires transcendental arithmetic: it is a fixed-point "
                  "decomposition stopped on a tolerance, and its Coxian fit needs a square root");

    const std::size_t M = L.rows();
    const std::size_t R = L.cols();
    if (R != N.size()) throw InputError("pfqn_marie: L and N disagree on the class count");
    if (M == 0) throw InputError("pfqn_marie: no stations");
    if (!Z.empty() && Z.size() != R) throw InputError("pfqn_marie: Z has the wrong length");
    if (!scv.empty() && (scv.rows() != M || scv.cols() != R))
        throw InputError("pfqn_marie: scv has the wrong shape");
    if (maxiter < 1) throw InputError("pfqn_marie: maxiter must be at least one");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r)
            if (!(L(i, r) > zero))
                throw InputError("pfqn_marie: every service demand must be strictly positive");

    Matrix<T> SCV = scv;
    if (SCV.empty()) SCV = Matrix<T>(M, R, one);
    std::vector<T> Zv = Z;
    if (Zv.empty()) Zv.assign(R, zero);

    if (R > 1) {
        if (!nservers.empty())
            for (std::size_t i = 0; i < nservers.size(); ++i)
                if (nservers[i] != 1)
                    throw InputError(
                        "pfqn_marie: the multiclass path has no multiserver isolation chain, as in "
                        "the reference");
        return detail::marie_multi(L, N, Zv, SCV, tol, maxiter);
    }

    // ------------------------- single class -------------------------------
    const int Nt = N[0];
    if (Nt < 1) throw InputError("pfqn_marie: the population must be at least one");
    std::vector<int> ns = nservers;
    if (ns.empty()) ns.assign(M, 1);
    if (ns.size() == 1 && M > 1) ns.assign(M, ns[0]);
    if (ns.size() != M) throw InputError("pfqn_marie: nservers has the wrong length");
    for (std::size_t i = 0; i < M; ++i)
        if (ns[i] < 1) throw InputError("pfqn_marie: the server count must be at least one");

    T Ztot = zero;
    for (std::size_t r = 0; r < Zv.size(); ++r) Ztot += Zv[r];
    Matrix<T> Zmat(1, 1);
    Zmat(0, 0) = Ztot;

    std::vector<MarieCoxFit<T>> fit(M);
    for (std::size_t i = 0; i < M; ++i) fit[i] = marie_cox_fit(L(i, 0), SCV(i, 0));

    // Initial multipliers relative to the base rate 1/L(i): min(n, m).
    Matrix<T> mu(M, static_cast<std::size_t>(Nt), one);
    for (std::size_t i = 0; i < M; ++i)
        for (int n = 1; n <= Nt; ++n)
            mu(i, static_cast<std::size_t>(n - 1)) = num_traits<T>::from_int(n < ns[i] ? n : ns[i]);

    MarieResult<T> res;
    res.it = 0;
    res.X.assign(1, zero);
    res.Q = Matrix<T>(M, 1, zero);
    res.U = Matrix<T>(M, 1, zero);
    res.C = Matrix<T>(1, 1, zero);
    const T tolT = num_traits<T>::from_double(tol);
    std::vector<T> lam(static_cast<std::size_t>(Nt), zero);

    while (res.it < maxiter) {
        ++res.it;
        const MvaLdResult<T> ld = pfqn_mvald(L, N, Zmat, mu);
        Matrix<T> muNew = mu;
        for (std::size_t i = 0; i < M; ++i) {
            for (int n = 0; n < Nt; ++n) {
                const T pn = ld.PI(i, static_cast<std::size_t>(n));
                lam[static_cast<std::size_t>(n)] =
                    pn > zero ? T(mu(i, static_cast<std::size_t>(n)) / L(i, 0) *
                                  ld.PI(i, static_cast<std::size_t>(n) + 1) / pn)
                              : zero;
            }
            const std::vector<T> muabs =
                detail::marie_isol_condtput(lam, fit[i].mu, fit[i].phi, Nt, ns[i]);
            for (int n = 0; n < Nt; ++n)
                muNew(i, static_cast<std::size_t>(n)) = muabs[static_cast<std::size_t>(n)] * L(i, 0);
        }
        T delta = zero;
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t n = 0; n < static_cast<std::size_t>(Nt); ++n) {
                const T d = num_abs(T(muNew(i, n) - mu(i, n)));
                if (d > delta) delta = d;
            }
        mu = muNew;
        res.X[0] = ld.XN[0];
        for (std::size_t i = 0; i < M; ++i) {
            res.Q(i, 0) = ld.QN(i, 0);
            res.U(i, 0) = ld.UN[i];
        }
        res.C(0, 0) = ld.CN[0];
        if (delta < tolT) break;
    }
    res.mu = mu;
    return res;
}

/** Reference defaults: tol 1e-8, maxiter 1000, single server everywhere. */
template <class T>
MarieResult<T> pfqn_marie(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z,
                          const Matrix<T>& scv) {
    return pfqn_marie(L, N, Z, scv, 1e-8, 1000, std::vector<int>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_MARIE_H
