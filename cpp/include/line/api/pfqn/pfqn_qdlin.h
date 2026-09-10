/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_QDLIN_H
#define LINE_API_PFQN_QDLIN_H

/**
 * QD-LIN: the Linearizer arm of AMVA-LD, on a plain demand matrix.
 *
 * Array-level twin of what SolverMVA computes for method='qdlin': the
 * Linearizer of Chandy and Neuse, Commun. ACM 25(2), 1982, run inside the
 * queue-dependent AMVA framework of Casale, Perez and Wang (IFIP PERFORMANCE
 * 2015), so the load-dependent term g_k is evaluated at the CORRECTED
 * arrival-instant queue rather than at the plain one.
 *
 * Port of python/line_solver/api/pfqn/qdlin.py, which is itself a transcription
 * of solver_amvald restricted to the domain a demand matrix describes: closed
 * classes only, one chain per class, unit visits, PS queueing stations and one
 * optional delay carrying Z. Within that domain it reproduces the native-Python
 * SolverMVA(model,'qdlin') to machine precision on random instances, which is
 * what this kernel is for; it is NOT an independent re-derivation.
 *
 * TWO PROPERTIES OF THE REFERENCE ARE REPRODUCED DELIBERATELY, not inherited by
 * accident, and a caller comparing against a textbook Linearizer sees both:
 *
 *  1. THE GAMMA CORRECTION IS CLASS-AGGREGATE, STORED IN SLICE 0. solver_amvald
 *     allocates the (K, M, K) per-class Linearizer array for qdlin but writes
 *     the class-aggregate correction into it with a two-subscript assignment,
 *     gamma(s,k) = sum_r Q_s(k,r)/(Nt-1) - sum_r Q(k,r)/Nt, which MATLAB
 *     linear-indexes to (s,k,1). Slices 1..K-1 stay zero while every reader
 *     indexes gamma per class, so the correction reaching the residence time is
 *     N_0*gamma(r,k,0) - [r==0]*gamma(r,k,0): the aggregate correction scaled by
 *     the population of CHAIN 0 alone, with the self term removed only for
 *     chain 0. It coincides with the queue-dependent AMVA form (Nt-1)*gamma_agg
 *     iff K = 1, so single-chain models are unaffected and multichain ones are
 *     not. method='lin' takes the per-class form instead.
 *  2. A SINGLE-SERVER STATION STILL CARRIES A SOFTMIN TERM. The multiserver
 *     factor is pfqn_lldfun(1 + arrival-instant total, {}, nservers), whose
 *     softmin at c = 1 is not exactly 1, so qdlin does not reduce to a textbook
 *     single-server AMVA even when every station has one server.
 *
 * MU AND NSERVERS ARE DIFFERENT MECHANISMS, unlike in pfqn_qdamva, which folds
 * the multiserver curve into mu. Here mu is sn.lldscaling, an interpolated rate
 * multiplier per station, and nservers is the server count feeding the softmin
 * term. A c-server station is nservers[k] = c, NOT a mu row of min(1..smax, c);
 * passing the latter reproduces Queue.setLoadDependence, a different station.
 *
 * THE WAIT-FACTOR FLOOR IS A SEPARATE KNOB from the convergence tolerance, and
 * is the native-Python solver's: MATLAB and the C++ solver_amvald do not clamp
 * the wait factor at all. It is load-bearing for qdlin, whose class-aggregate
 * correction drives the factor negative at a lightly loaded station.
 *
 * Arithmetic: TRANSCENDENTAL-GATED, inherited whole from pfqn_lldfun.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/pfqn/pfqn_lldfun.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** What `pfqn_qdlin` returns: the fixed point and how it was reached. */
template <class T>
struct QdLinResult {
    Matrix<T> Q;  ///< (M x R) mean queue lengths at the queueing stations
    Matrix<T> U;  ///< (M x R) per-class utilizations, the analytic T*S/c
    Matrix<T> R;  ///< (M x R) per-class residence times
    Matrix<T> X;  ///< (1 x R) per-class throughputs
    Matrix<T> C;  ///< (1 x R) per-class cycle times, think time included
    std::size_t iter = 0;  ///< number of forward evaluations performed
};

namespace detail {

/** One forward evaluation, solver_amvald_forward restricted to PS and INF. */
template <class T>
void qdlin_forward(const Matrix<T>& ST, const std::vector<double>& srv,
                   const std::vector<char>& isdelay, const Matrix<T>& mu,
                   const std::vector<Matrix<T>>& gamma, const Matrix<T>& Qin,
                   const std::vector<T>& Nin, const std::vector<std::size_t>& nnz, double wtol,
                   Matrix<T>& W, Matrix<T>& STeff) {
    const std::size_t Ms = static_cast<std::size_t>(ST.rows());
    const std::size_t K = static_cast<std::size_t>(ST.cols());
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    T Ntin = zero;
    for (std::size_t r = 0; r < K; ++r) Ntin += Nin[r];
    const T delta = (num_traits<T>::to_double(Ntin) > 0.0) ? T((Ntin - one) / Ntin) : one;

    std::vector<T> dcl(K, one);
    for (std::size_t r : nnz) dcl[r] = T((Nin[r] - one) / Nin[r]);

    // Arrival-instant queue lengths, class-aggregate and per class.
    std::vector<T> interp(Ms, zero);
    Matrix<T> totArvl(Ms, K, zero);
    for (std::size_t k = 0; k < Ms; ++k) {
        T sumQk = zero;
        for (std::size_t r : nnz) sumQk += Qin(k, r);
        interp[k] = T(delta * sumQk);
        for (std::size_t r : nnz) totArvl(k, r) = T(dcl[r] * Qin(k, r) + sumQk - Qin(k, r));
    }

    // Load-dependent term, evaluated at the gamma-corrected arrival-instant queue.
    Matrix<T> lldterm(Ms, K, one);
    const std::vector<double> noservers;
    if (!nnz.empty()) {
        std::vector<T> arg(Ms, zero);
        for (std::size_t r : nnz) {
            for (std::size_t k = 0; k < Ms; ++k) {
                T corr = zero;
                for (std::size_t c : nnz) corr += T(Nin[c] * gamma[r](k, c));
                corr -= gamma[r](k, r);
                arg[k] = T(one + interp[k] + corr);
            }
            const std::vector<T> g = pfqn_lldfun(arg, mu, noservers);
            for (std::size_t k = 0; k < Ms; ++k) lldterm(k, r) = g[k];
        }
    } else {
        std::vector<T> arg(Ms, zero);
        for (std::size_t k = 0; k < Ms; ++k) arg[k] = T(one + interp[k]);
        const std::vector<T> g = pfqn_lldfun(arg, mu, noservers);
        for (std::size_t k = 0; k < Ms; ++k)
            for (std::size_t r = 0; r < K; ++r) lldterm(k, r) = g[k];
    }

    // Multiserver term; the 'default' rule leaves PS on the softmin arm.
    std::vector<T> msarg(Ms, zero);
    if (!nnz.empty() && num_traits<T>::to_double(Ntin) > 0.0) {
        const T shrink = T((Ntin - one) / Ntin);
        for (std::size_t k = 0; k < Ms; ++k) {
            T acc = zero;
            for (std::size_t s : nnz) {
                T gs = zero;
                for (std::size_t r : nnz) gs += T(shrink * Nin[r] * gamma[s](k, r));
                acc += gs;
            }
            const T mean = T(acc / num_traits<T>::from_int(static_cast<long>(nnz.size())));
            msarg[k] = T(one + interp[k] + mean);
        }
    } else {
        for (std::size_t k = 0; k < Ms; ++k) msarg[k] = T(one + interp[k]);
    }
    const Matrix<T> nolld;
    const std::vector<T> msterm = pfqn_lldfun(msarg, nolld, srv);

    STeff = Matrix<T>(Ms, K, zero);
    for (std::size_t r : nnz)
        for (std::size_t k = 0; k < Ms; ++k)
            STeff(k, r) = T(ST(k, r) * lldterm(k, r) * msterm[k]);

    W = Matrix<T>(Ms, K, zero);
    const T floor = num_traits<T>::from_double(wtol);
    for (std::size_t r : nnz) {
        for (std::size_t k = 0; k < Ms; ++k) {
            if (isdelay[k]) {
                W(k, r) = STeff(k, r);
                continue;
            }
            T corr = zero;
            for (std::size_t c : nnz) corr += T(Nin[c] * gamma[r](k, c));
            corr -= gamma[r](k, r);
            const T factor = T(one + totArvl(k, r) + corr);
            W(k, r) = T(STeff(k, r) * (factor < floor ? floor : factor));
        }
    }
}

}  // namespace detail

/**
 * @param L        (M x R) service demand matrix, queueing stations only.
 * @param N        (R) population vector, finite.
 * @param Z        (R) think time vector; a delay station carrying it is
 *                 appended to the station list when any entry is positive,
 *                 exactly as the equivalent Network would hold one. Empty
 *                 means no think time.
 * @param mu       (M x smax) load-dependent rate multipliers, sn.lldscaling;
 *                 empty means none.
 * @param nservers (M) server counts; empty means one server everywhere.
 * @param tol      convergence tolerance on the queue lengths, LINE's iter_tol.
 * @param maxiter  iteration budget, LINE's iter_max. The outer sweep and each
 *                 inner sweep are capped at sqrt(maxiter) and the total number
 *                 of forward evaluations at min(maxiter, 10000).
 * @param wtol     floor on the AMVA wait factor, LINE's options.tol. A
 *                 DIFFERENT knob from tol, with its own default: SolverMVA
 *                 passes iter_tol to the fixed point but never sets
 *                 options.tol, so the floor stays at the lineDefaults 1e-4
 *                 while the fixed point converges to 1e-6.
 */
template <class T>
QdLinResult<T> pfqn_qdlin(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z,
                          const Matrix<T>& mu, const std::vector<double>& nservers,
                          double tol = 1e-6, std::size_t maxiter = 1000, double wtol = 1e-4) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_qdlin requires transcendental arithmetic");
    const std::size_t M = static_cast<std::size_t>(L.rows());
    const std::size_t K = static_cast<std::size_t>(L.cols());
    if (N.size() != K)
        throw InputError("pfqn_qdlin: the population vector must have one entry per class");
    if (!Z.empty() && Z.size() != K)
        throw InputError("pfqn_qdlin: the think-time vector must have one entry per class");
    if (!nservers.empty() && nservers.size() != M)
        throw InputError("pfqn_qdlin: the server-count vector must have one entry per station");
    for (std::size_t r = 0; r < K; ++r)
        if (!std::isfinite(num_traits<T>::to_double(N[r])))
            throw InputError(
                "pfqn_qdlin: an infinite population is not supported, closed classes only");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    QdLinResult<T> out;
    out.Q = Matrix<T>(M, K, zero);
    out.U = Matrix<T>(M, K, zero);
    out.R = Matrix<T>(M, K, zero);
    out.X = Matrix<T>(1, K, zero);
    out.C = Matrix<T>(1, K, zero);

    T Ntot = zero;
    for (std::size_t r = 0; r < K; ++r) Ntot += N[r];
    if (!(num_traits<T>::to_double(Ntot) > 0.0)) return out;

    // Station list: the delay, when there is one, then the queueing stations.
    bool hasDelay = false;
    for (std::size_t r = 0; r < K && !Z.empty(); ++r)
        if (num_traits<T>::to_double(Z[r]) > 0.0) hasDelay = true;
    const std::size_t off = hasDelay ? 1 : 0;
    const std::size_t Ms = M + off;

    Matrix<T> ST(Ms, K, zero);
    std::vector<double> srv(Ms, 1.0);
    std::vector<char> isdelay(Ms, 0);
    if (hasDelay) {
        for (std::size_t r = 0; r < K; ++r) ST(0, r) = Z[r];
        srv[0] = std::numeric_limits<double>::infinity();
        isdelay[0] = 1;
    }
    for (std::size_t k = 0; k < M; ++k) {
        for (std::size_t r = 0; r < K; ++r) ST(k + off, r) = L(k, r);
        srv[k + off] = nservers.empty() ? 1.0 : nservers[k];
    }
    Matrix<T> muFull;
    if (!mu.empty()) {
        const std::size_t smax = static_cast<std::size_t>(mu.cols());
        muFull = Matrix<T>(Ms, smax, one);
        for (std::size_t k = 0; k < M; ++k)
            for (std::size_t j = 0; j < smax; ++j) muFull(k + off, j) = mu(k, j);
    }

    std::vector<std::size_t> nnz;
    for (std::size_t r = 0; r < K; ++r)
        if (num_traits<T>::to_double(N[r]) > 0.0) nnz.push_back(r);

    // Balanced initialization, as in solver_amvald.
    Matrix<T> Q(Ms, K, zero);
    const T Msf = num_traits<T>::from_int(static_cast<long>(Ms));
    const T share = T(one / Msf);
    for (std::size_t r : nnz)
        for (std::size_t k = 0; k < Ms; ++k) Q(k, r) = T(share * N[r]);
    std::vector<T> X(K, zero);
    for (std::size_t r = 0; r < K; ++r) {
        T col = zero;
        for (std::size_t k = 0; k < Ms; ++k) col += ST(k, r);
        if (num_traits<T>::to_double(col) > 0.0) X[r] = T(one / col);
    }

    Matrix<T> Umat(Ms, K, zero);
    for (std::size_t r : nnz)
        for (std::size_t k = 0; k < Ms; ++k)
            Umat(k, r) = std::isfinite(srv[k])
                             ? T(ST(k, r) * X[r] / num_traits<T>::from_double(srv[k]))
                             : T(ST(k, r) * X[r]);

    std::vector<Matrix<T>> gamma(K, Matrix<T>(Ms, K, zero));
    Matrix<T> Tput(Ms, K, zero);

    const double omicron = 0.5;
    const T om = num_traits<T>::from_double(omicron);
    const T omc = num_traits<T>::from_double(1.0 - omicron);
    const double maxSweep = std::sqrt(static_cast<double>(maxiter));
    const std::size_t maxTotiter = std::min<std::size_t>(maxiter, 10000);

    Matrix<T> W, STeff(Ms, K, zero);
    Matrix<T> Qouter = Q;
    std::size_t outerIter = 0;
    while (static_cast<double>(outerIter) < maxSweep && out.iter <= maxTotiter) {
        if (outerIter >= 2) {
            double gap = 0.0;
            for (std::size_t k = 0; k < Ms; ++k)
                for (std::size_t r = 0; r < K; ++r)
                    gap = std::max(gap, std::fabs(num_traits<T>::to_double(Q(k, r)) -
                                                  num_traits<T>::to_double(Qouter(k, r))));
            if (gap <= tol) break;
        }
        ++outerIter;
        Qouter = Q;

        // Linearizer recursion: one sweep at each reduced population N - 1_s.
        bool exhausted = false;
        for (std::size_t s = 0; s < K && !exhausted; ++s) {
            if (!(num_traits<T>::to_double(N[s]) > 0.0)) continue;
            std::vector<T> Ns = N;
            Ns[s] = T(Ns[s] - one);
            const T shrink = T((Ntot - one) / Ntot);
            Matrix<T> Qs = Q;
            for (std::size_t k = 0; k < Ms; ++k)
                for (std::size_t r = 0; r < K; ++r) Qs(k, r) = T(Qs(k, r) * shrink);
            std::vector<T> Xs = X;
            for (std::size_t r = 0; r < K; ++r) Xs[r] = T(Xs[r] * shrink);

            std::size_t iterS = 0;
            Matrix<T> QsPrev = Qs;
            while (static_cast<double>(iterS) <= maxSweep) {
                if (iterS >= 2) {
                    double gap = 0.0;
                    for (std::size_t k = 0; k < Ms; ++k)
                        for (std::size_t r = 0; r < K; ++r)
                            gap = std::max(gap, std::fabs(num_traits<T>::to_double(Qs(k, r)) -
                                                          num_traits<T>::to_double(QsPrev(k, r))));
                    if (gap <= tol) break;
                }
                ++iterS;
                QsPrev = Qs;
                const std::vector<T> XsPrev = Xs;

                detail::qdlin_forward(ST, srv, isdelay, muFull, gamma, QsPrev, Ns, nnz, wtol, W,
                                      STeff);
                ++out.iter;
                if (out.iter >= maxTotiter) {
                    exhausted = true;
                    break;
                }

                for (std::size_t r : nnz) {
                    T wsum = zero;
                    for (std::size_t k = 0; k < Ms; ++k) wsum += W(k, r);
                    if (num_traits<T>::to_double(wsum) == 0.0 ||
                        !(num_traits<T>::to_double(Ns[r]) > 0.0)) {
                        Xs[r] = zero;
                    } else if (num_traits<T>::to_double(wsum) > 1e-14) {
                        Xs[r] = T(om * Ns[r] / wsum + omc * XsPrev[r]);
                    } else {
                        Xs[r] = XsPrev[r];
                    }
                    for (std::size_t k = 0; k < Ms; ++k)
                        Qs(k, r) = T(om * Xs[r] * W(k, r) + omc * QsPrev(k, r));
                }
            }

            // Class-aggregate correction into slice 0, see the header.
            if (num_traits<T>::to_double(Ntot) > 1.0) {
                for (std::size_t k = 0; k < Ms; ++k) {
                    T a = zero, b = zero;
                    for (std::size_t r = 0; r < K; ++r) {
                        a += QsPrev(k, r);
                        b += Qouter(k, r);
                    }
                    gamma[s](k, 0) = T(a / (Ntot - one) - b / Ntot);
                }
            } else {
                for (std::size_t k = 0; k < Ms; ++k) gamma[s](k, 0) = zero;
            }
        }
        if (exhausted) break;

        // Sweep at the full population N.
        std::size_t innerIter = 0;
        Matrix<T> Qprev = Q;
        while (static_cast<double>(innerIter) <= maxSweep) {
            if (innerIter >= 2) {
                double gap = 0.0;
                for (std::size_t k = 0; k < Ms; ++k)
                    for (std::size_t r = 0; r < K; ++r)
                        gap = std::max(gap, std::fabs(num_traits<T>::to_double(Q(k, r)) -
                                                      num_traits<T>::to_double(Qprev(k, r))));
                if (gap <= tol) break;
            }
            ++innerIter;
            Qprev = Q;
            const std::vector<T> Xprev = X;
            const Matrix<T> Uprev = Umat;

            detail::qdlin_forward(ST, srv, isdelay, muFull, gamma, Qprev, N, nnz, wtol, W, STeff);
            ++out.iter;
            if (out.iter >= maxTotiter) {
                exhausted = true;
                break;
            }

            for (std::size_t r : nnz) {
                T wsum = zero;
                for (std::size_t k = 0; k < Ms; ++k) wsum += W(k, r);
                if (num_traits<T>::to_double(wsum) == 0.0) {
                    X[r] = zero;
                } else {
                    out.C(0, r) = wsum;
                    if (num_traits<T>::to_double(wsum) > 1e-14)
                        X[r] = T(om * N[r] / wsum + omc * Xprev[r]);
                    else
                        X[r] = Xprev[r];
                }
                for (std::size_t k = 0; k < Ms; ++k) {
                    Q(k, r) = T(om * X[r] * W(k, r) + omc * Qprev(k, r));
                    Tput(k, r) = X[r];
                    Umat(k, r) = T(om * STeff(k, r) * X[r] + omc * Uprev(k, r));
                }
            }
        }
        if (exhausted) break;
    }
    // Utilization capping, as in solver_amvald: a queueing station whose class
    // utilizations sum above one has them renormalized in proportion to STeff.
    // Delay stations are exempt.
    for (std::size_t k = 0; k < Ms; ++k) {
        if (isdelay[k]) continue;
        T usum = zero;
        for (std::size_t r = 0; r < K; ++r) usum += Umat(k, r);
        if (!(num_traits<T>::to_double(usum) > 1.0)) continue;
        T denom = zero;
        for (std::size_t r = 0; r < K; ++r) denom += T(STeff(k, r) * X[r]);
        if (!(num_traits<T>::to_double(denom) > 0.0)) continue;
        for (std::size_t r = 0; r < K; ++r)
            if (num_traits<T>::to_double(STeff(k, r)) > 0.0)
                Umat(k, r) = T(STeff(k, r) * X[r] / denom);
    }

    // WHICH UTILIZATION SolverMVA REPORTS DEPENDS ON THE MODEL. Its analyzer
    // forwards the iterated Uchain to sn_deaggregate_chain_results ONLY under
    // lld, cd or jd scaling; with none of those the deaggregation recomputes
    // T*S/c from the NOMINAL demand instead, and the two differ by the
    // iteration residual. Reproduced on the same test, mu being the only one of
    // the three a demand matrix can carry.
    const bool iteratedU = !mu.empty();
    for (std::size_t r : nnz) {
        out.X(0, r) = X[r];
        for (std::size_t k = 0; k < M; ++k) {
            const std::size_t ks = k + off;
            out.Q(k, r) = Q(ks, r);
            out.U(k, r) = iteratedU ? Umat(ks, r)
                                    : (std::isfinite(srv[ks])
                                           ? T(ST(ks, r) * X[r] /
                                               num_traits<T>::from_double(srv[ks]))
                                           : T(ST(ks, r) * X[r]));
            out.R(k, r) = (num_traits<T>::to_double(Tput(ks, r)) > 0.0)
                              ? T(Q(ks, r) / Tput(ks, r))
                              : zero;
        }
    }
    return out;
}

/** Overload without a load-dependent lattice or explicit server counts. */
template <class T>
QdLinResult<T> pfqn_qdlin(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z,
                          double tol = 1e-6, std::size_t maxiter = 1000, double wtol = 1e-4) {
    return pfqn_qdlin(L, N, Z, Matrix<T>(), std::vector<double>(), tol, maxiter, wtol);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_QDLIN_H
