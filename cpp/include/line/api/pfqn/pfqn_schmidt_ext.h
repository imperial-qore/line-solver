/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_SCHMIDT_EXT_H
#define LINE_API_PFQN_SCHMIDT_EXT_H

/**
 * Extended Schmidt MVA with queue-aware alpha corrections.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_schmidt_ext.m, cross-checked
 * against the second half of jar/src/main/java/jline/api/pfqn/mva/
 * Pfqn_schmidt_amva.java (that file carries both the plain and the extended
 * recursion). The Java entry point takes SERVICE RATES and inverts them to
 * demands internally, where MATLAB and this port take demands; that is an API
 * divergence, not an algorithmic one. Reference:
 * R. Schmidt, "An approximate MVA algorithm for exponential, class-dependent
 * multiple server stations", Performance Evaluation 29(4), 1997.
 *
 * The population recursion is the one of pfqn_schmidt. What the extension adds
 * is the mean service time B_c(n) used at a class-dependent multiserver FCFS
 * station. Plain Schmidt weights the composition n by the demands themselves,
 *
 *   B_c(n) = D(i,c) + max(0, |n| - s) / (s (|n| - 1)) (sum_t n_t D(i,t) - D(i,c));
 *
 * the extension replaces the second term by max(0, |n| - s) times a mean
 * INTERDEPARTURE time read off an auxiliary solve. For each class r a tagged
 * single-job class R+1 is appended at station i with that station's class-r
 * demand, the population is dropped to N - e_r, plain pfqn_schmidt is run on
 * the (R+1)-class model, and its utilizations u give
 *
 *   alpha(i)     = sum_{s<=R} u(i,s) - u(i,R+1),
 *   1/interdep   = s sum_{s: D(i,s)>0} (u(i,s)/alpha(i)) / D(i,s),
 *
 * i.e. a demand-weighted harmonic mean over the classes actually competing for
 * the servers, with the tagged job's own utilization removed. That is the
 * "queue-aware" correction: the departure rate seen by a waiting job reflects
 * the class mix at the station rather than its own demand.
 *
 * Reference behaviour preserved verbatim: the alphas are computed only for
 * FCFS stations whose demands are class dependent, and used only where
 * N_c > 1; the class-independent multiserver marginal is the binomial form of
 * the reference, not the recursive one of plain Schmidt; the idle-state
 * probability is floored at 1e-12 (2.2e-16 for the PS branch); and the servers
 * used in the queue-length update block are read with the class index left
 * over from the preceding loop, which for an (M x R) server matrix is the LAST
 * class. That last item is a reference quirk, reproduced because it selects
 * the numbers the reference produces; it is inert whenever the server counts
 * do not vary by class.
 *
 * The extra SchedStrategy.FCFS that the reference appends to `sched` before
 * the auxiliary solve is dropped here: it lengthens the vector to M + 1 for an
 * M-station model and pfqn_schmidt never reads past M, so it is inert.
 *
 * The auxiliary solve widens the model to R + 1 classes, so a per-class server
 * matrix cannot be carried into it; the reference indexes S(ist,c) with
 * c = R + 1 and would raise an out-of-bounds error. This port therefore
 * requires a per-station (M x 1) server vector whenever any alpha is needed,
 * and says so rather than inventing a widening rule.
 *
 * Arithmetic: INEXACT BY CONSTRUCTION. The correction is an approximation
 * whose alpha factors come from an auxiliary approximate solve, and the
 * marginal probabilities use non-integer binomial powers; it is gated on
 * has_transcendental accordingly. Plain pfqn_schmidt, being a finite rational
 * recursion, stays instantiable at Rational and is not gated.
 */

#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_amva_common.h"
#include "line/api/pfqn/pfqn_schmidt.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_schmidt_ext, mirroring [XN,QN,UN,CN]. */
template <class T>
struct SchmidtExtResult {
    std::vector<T> XN;  ///< (R) per-class throughput
    Matrix<T> QN;       ///< (M x R) mean queue length
    Matrix<T> UN;       ///< (M x R) utilization, D X / s
    Matrix<T> CN;       ///< (M x R) residence time
};

namespace detail {

/** matlab getBcn: the plain Schmidt composition-weighted mean service time. */
template <class T>
T schmidt_ext_bcn(const Matrix<T>& D, std::size_t i, std::size_t c, const std::vector<int>& nvec,
                  int ns) {
    T bcn = D(i, c);
    long nsum = 0;
    for (int t : nvec) nsum += t;
    if (nsum > 1) {
        const T eps = num_traits<T>::from_double(1e-12);
        T sumVal = num_traits<T>::from_int(0);
        for (std::size_t t = 0; t < nvec.size(); ++t)
            sumVal += num_traits<T>::from_int(nvec[t]) * D(i, t);
        const T num = num_traits<T>::from_int(nsum - ns > 0 ? nsum - ns : 0);
        const T den0 = num_traits<T>::from_int(ns * (nsum - 1));
        const T den = den0 > eps ? den0 : eps;
        bcn += num / den * (sumVal - D(i, c));
    }
    return bcn;
}

/** matlab getBcnExt: the queue-aware mean interdeparture correction. */
template <class T>
T schmidt_ext_bcn_ext(const Matrix<T>& u, const Matrix<T>& D, std::size_t i, std::size_t c,
                      const std::vector<int>& nvec, std::size_t C, int ns) {
    const T zero = num_traits<T>::from_int(0);
    T nonPinned = zero;
    for (std::size_t s = 0; s < C; ++s) nonPinned += u(i, s);
    nonPinned -= u(i, C);  // the tagged class sits in column C (0-based)

    T weighted = zero;
    if (nonPinned > zero)
        for (std::size_t s = 0; s < C; ++s)
            if (D(i, s) > zero) weighted += (u(i, s) / nonPinned) / D(i, s);

    const T interdep =
        weighted > zero ? T(num_traits<T>::from_int(1) / (num_traits<T>::from_int(ns) * weighted))
                        : zero;
    T bcn = D(i, c);
    long nsum = 0;
    for (int t : nvec) nsum += t;
    if (nsum > 1) bcn += num_traits<T>::from_int(nsum - ns > 0 ? nsum - ns : 0) * interdep;
    return bcn;
}

}  // namespace detail

/**
 * @param D     (M x R) service demands
 * @param N     (R) population per class
 * @param S     (M x 1) or (M x R) server counts; (M x 1) is required when any
 *              class-dependent FCFS multiserver station needs an alpha
 * @param sched (M) scheduling discipline per station
 */
template <class T>
SchmidtExtResult<T> pfqn_schmidt_ext(const Matrix<T>& D, const std::vector<int>& N,
                                     const Matrix<int>& S,
                                     const std::vector<SchedStrategy>& sched) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_schmidt_ext requires transcendental arithmetic: its alpha correction is "
                  "an approximation drawn from an auxiliary approximate solve, and its marginal "
                  "probabilities use non-integer binomial powers");

    const std::size_t M = D.rows();
    const std::size_t R = N.size();
    if (!D.empty() && D.cols() != R)
        throw InputError("pfqn_schmidt_ext: D and N disagree on the class count");
    if (sched.size() != M) throw InputError("pfqn_schmidt_ext: sched has the wrong station count");
    if (S.rows() != M || (S.cols() != R && S.cols() != 1))
        throw InputError("pfqn_schmidt_ext: server-count matrix has the wrong shape");
    for (int n : N)
        if (n < 0) throw InputError("pfqn_schmidt_ext: negative population");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T tiny = num_traits<T>::from_double(1e-12);
    const T epsT = num_traits<T>::from_double(2.220446049250313e-16);

    SchmidtExtResult<T> res;
    res.XN.assign(R, zero);
    res.QN = Matrix<T>(M, R, zero);
    res.UN = Matrix<T>(M, R, zero);
    res.CN = Matrix<T>(M, R, zero);
    if (M == 0 || R == 0) return res;

    const auto nserv = [&](std::size_t i, std::size_t c) {
        return S.cols() == 1 ? S(i, 0) : S(i, c);
    };
    std::vector<bool> classIndep(M, true);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 1; c < R; ++c)
            if (D(i, c) != D(i, 0)) classIndep[i] = false;

    // ---- alphas: one auxiliary Schmidt solve per class-dependent FCFS pair --
    std::vector<Matrix<T>> alphas(M * R);
    for (std::size_t i = 0; i < M; ++i) {
        if (sched[i] != SchedStrategy::FCFS || classIndep[i]) continue;
        if (S.cols() != 1)
            throw InputError(
                "pfqn_schmidt_ext: a class-dependent FCFS station needs the alpha correction, "
                "whose auxiliary solve widens the model to R+1 classes; pass a per-station "
                "(M x 1) server vector");
        for (std::size_t r = 0; r < R; ++r) {
            Matrix<T> Dmod(M, R + 1, zero);
            std::vector<int> Nmod(R + 1, 0);
            for (std::size_t k = 0; k < R; ++k) {
                Nmod[k] = k == r ? N[k] - 1 : N[k];
                for (std::size_t j = 0; j < M; ++j) {
                    Dmod(j, k) = D(j, k);
                    Dmod(j, R) = j == i ? D(j, r) : zero;
                }
            }
            Nmod[R] = 1;
            if (Nmod[r] < 0)
                throw InputError(
                    "pfqn_schmidt_ext: a class-dependent FCFS station needs the alpha correction "
                    "at N - e_r, which is negative for an empty class");
            alphas[i * R + r] = pfqn_schmidt(Dmod, Nmod, S, sched).UN;
        }
    }

    // ---- population recursion ---------------------------------------------
    const std::vector<std::size_t> prods = plane_sizes(N);
    const std::size_t total = population_count(N);
    long Ntot = 0;
    for (int n : N) Ntot += n;

    std::vector<detail::SchmidtPc> kind(M, detail::SchmidtPc::None);
    for (std::size_t i = 0; i < M; ++i) {
        bool single = true;
        for (std::size_t c = 0; c < (S.cols() == 1 ? std::size_t(1) : R); ++c)
            if (nserv(i, c) != 1) single = false;
        switch (sched[i]) {
            case SchedStrategy::INF:
                break;
            case SchedStrategy::PS:
                if (!single) kind[i] = detail::SchmidtPc::Scalar;
                break;
            case SchedStrategy::FCFS:
                if (classIndep[i]) {
                    if (!single) kind[i] = detail::SchmidtPc::Scalar;
                } else {
                    kind[i] = detail::SchmidtPc::Vector;
                }
                break;
        }
    }

    std::vector<Matrix<T>> Lq(M, Matrix<T>(R, total, zero));
    std::vector<Matrix<T>> Pc(M);
    for (std::size_t i = 0; i < M; ++i) {
        if (kind[i] == detail::SchmidtPc::Scalar)
            Pc[i] = Matrix<T>(static_cast<std::size_t>(1 + Ntot), total, zero);
        else if (kind[i] == detail::SchmidtPc::Vector)
            Pc[i] = Matrix<T>(total, total, zero);
        if (kind[i] != detail::SchmidtPc::None) Pc[i](0, 0) = one;
    }

    Matrix<T> xtab(R, total, zero);
    std::vector<T> w(M * R * total, zero);

    std::vector<int> kvec(R, 0);
    std::size_t hlast = 0;
    bool more = true;
    while (more) {
        const std::size_t hk = pop_index(kvec, prods);
        hlast = hk;
        long kpop = 0;
        for (int t : kvec) kpop += t;

        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t c = 0; c < R; ++c) {
                if (kvec[c] <= 0) continue;
                const std::size_t hkc = hk - prods[c];
                const int ns = nserv(i, c);
                T& wi = w[(i * R + c) * total + hk];
                if (sched[i] == SchedStrategy::INF) {
                    wi = D(i, c);
                    continue;
                }
                const bool vectorPc = kind[i] == detail::SchmidtPc::Vector;
                if (!vectorPc || ns == 1) {
                    T qtot = zero;
                    for (std::size_t r = 0; r < R; ++r) qtot += Lq[i](r, hkc);
                    if (ns == 1) {
                        wi = D(i, c) * (one + qtot);
                    } else {
                        const T nsT = num_traits<T>::from_int(ns);
                        wi = D(i, c) / nsT * (one + qtot);
                        for (int j = 1; j <= ns - 1; ++j)
                            wi += num_traits<T>::from_int(ns - j) *
                                  Pc[i](static_cast<std::size_t>(j - 1), hkc) * (D(i, c) / nsT);
                    }
                } else {
                    std::vector<int> nvec(R, 0);
                    bool more_n = true;
                    while (more_n) {
                        if (nvec[c] > 0) {
                            const std::size_t hnc = pop_index(nvec, prods) - prods[c];
                            const Matrix<T>& al = alphas[i * R + c];
                            const T Bcn = (N[c] > 1 && !al.empty())
                                              ? detail::schmidt_ext_bcn_ext(al, D, i, c, nvec, R, ns)
                                              : detail::schmidt_ext_bcn(D, i, c, nvec, ns);
                            wi += Bcn * Pc[i](hnc, hkc);
                        }
                        more_n = next_pop(nvec, kvec);
                    }
                }
            }

        for (std::size_t c = 0; c < R; ++c) {
            T denom = zero;
            for (std::size_t i = 0; i < M; ++i) denom += w[(i * R + c) * total + hk];
            xtab(c, hk) = denom > zero ? T(num_traits<T>::from_int(kvec[c]) / denom) : zero;
        }

        for (std::size_t i = 0; i < M; ++i) {
            for (std::size_t c = 0; c < R; ++c)
                Lq[i](c, hk) = xtab(c, hk) * w[(i * R + c) * total + hk];

            // The reference reads the server count with the class index left
            // over from the loop above, which is the last class.
            const int nsLast = nserv(i, R - 1);
            const int s0 = nserv(i, 0);

            if (sched[i] == SchedStrategy::PS) {
                if (nsLast > 1 && kind[i] == detail::SchmidtPc::Scalar) {
                    const long top = s0 < kpop ? s0 : kpop;
                    for (long n = 1; n <= top; ++n)
                        for (std::size_t c = 0; c < R; ++c) {
                            if (kvec[c] <= 0) continue;
                            const std::size_t hkc = hk - prods[c];
                            Pc[i](static_cast<std::size_t>(n), hk) +=
                                D(i, c) / num_traits<T>::from_int(n) * xtab(c, hk) *
                                Pc[i](static_cast<std::size_t>(n - 1), hkc);
                        }
                    if (top >= 1) {
                        T acc = zero;
                        for (long n = 1; n <= top; ++n)
                            acc += Pc[i](static_cast<std::size_t>(n), hk);
                        const T p0 = one - acc;
                        Pc[i](0, hk) = p0 > epsT ? p0 : epsT;
                    }
                }
            } else if (sched[i] == SchedStrategy::FCFS) {
                if (kind[i] == detail::SchmidtPc::Vector) {
                    T sumAll = zero;
                    std::vector<int> nvec(R, 0);
                    bool more_n = next_pop(nvec, kvec);  // skip the zero vector
                    while (more_n) {
                        const std::size_t hn = pop_index(nvec, prods);
                        long nsum = 0;
                        for (int t : nvec) nsum += t;
                        T prob = zero;
                        for (std::size_t r = 0; r < R; ++r) {
                            if (nvec[r] <= 0) continue;
                            const std::size_t hnc = hn - prods[r];
                            const std::size_t hkc = hk - prods[r];
                            const int nsr = nserv(i, r);
                            const Matrix<T>& al = alphas[i * R + r];
                            const T Bcn =
                                (N[r] > 1 && !al.empty())
                                    ? detail::schmidt_ext_bcn_ext(al, D, i, r, nvec, R, nsr)
                                    : detail::schmidt_ext_bcn(D, i, r, nvec, nsr);
                            prob += Bcn / num_traits<T>::from_int(nsum) * xtab(r, hk) *
                                    Pc[i](hnc, hkc);
                        }
                        Pc[i](hn, hk) = prob;
                        sumAll += prob;
                        more_n = next_pop(nvec, kvec);
                    }
                    const T p0 = one - sumAll;
                    Pc[i](0, hk) = p0 > tiny ? p0 : tiny;
                } else if (nsLast > 1 && kind[i] == detail::SchmidtPc::Scalar) {
                    // Class-independent multiserver: binomial marginal.
                    long Kj = 0;
                    for (std::size_t r = 0; r < R; ++r)
                        if (D(i, r) > zero) Kj += N[r];
                    T meanQ = zero;
                    for (std::size_t r = 0; r < R; ++r) meanQ += Lq[i](r, hk);
                    const long top = (nsLast < kpop ? nsLast : kpop) - 1;
                    for (long n = 1; n <= top; ++n) {
                        if (Kj <= 0 || n > Kj) continue;
                        const T frac = meanQ / num_traits<T>::from_int(Kj);
                        Pc[i](static_cast<std::size_t>(n), hk) =
                            num_nck<T>(static_cast<int>(Kj), static_cast<int>(n)) *
                            num_pow_int(frac, static_cast<unsigned>(n)) *
                            num_pow_int(T(one - frac), static_cast<unsigned>(Kj - n));
                    }
                    T sum1 = zero, sum2 = zero;
                    for (std::size_t r = 0; r < R; ++r) sum1 += D(i, r) * xtab(r, hk);
                    for (long n = 0; n <= nsLast - 1; ++n)
                        sum2 += num_traits<T>::from_int(nsLast - n) *
                                Pc[i](static_cast<std::size_t>(n), hk);
                    const T p0 = one - (sum1 + sum2) / num_traits<T>::from_int(nsLast);
                    Pc[i](0, hk) = p0 > tiny ? p0 : tiny;
                }
            }
        }

        more = next_pop(kvec, N);
    }

    for (std::size_t c = 0; c < R; ++c) res.XN[c] = xtab(c, hlast);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < R; ++c) {
            res.UN(i, c) = D(i, c) * res.XN[c] / num_traits<T>::from_int(nserv(i, 0));
            res.CN(i, c) = w[(i * R + c) * total + hlast];
            res.QN(i, c) = Lq[i](c, hlast);
        }
    return res;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_SCHMIDT_EXT_H
