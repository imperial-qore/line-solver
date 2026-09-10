/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_NC_SOLVER_NC_H
#define LINE_SOLVERS_NC_SOLVER_NC_H

/**
 * Port of `solver_nc.m`: the load-INDEPENDENT normalizing-constant analyzer.
 *
 * WHAT IT COMPUTES. One evaluation of the normalizing constant G(N) settles the
 * whole product-form network: the chain throughput is X_c = G(N - 1_c) / G(N),
 * and the chain queue length at station i follows from a constant evaluated on
 * the model with station i REPLICATED into a private class,
 *
 *     Q_ic = Zms_ic X_c + Lms_ic exp(lG_ar(i,c) - lG),
 *
 * which is the arrival-theorem identity written in constants rather than in a
 * recursion. That is why this file calls `pfqn_nc` once per chain and once per
 * (station, chain) pair rather than solving anything itself: the algorithms
 * live in `api/pfqn`, and choosing among them is `pfqn_nc`'s job.
 *
 * MULTISERVER IS SEIDMANN'S APPROXIMATION, NOT AN EXACT SOLVE. A station with
 * c > 1 servers is split into a demand L/c at a single server plus a pure delay
 * L (c-1)/c. It is exact only in the limit; the reference warns and continues,
 * and `@@SolverNC/runAnalyzer` converts genuine multiservers to load-dependent
 * rates BEFORE reaching here whenever the model is product-form, so this branch
 * serves the models that have no exact load-dependent form.
 *
 * WHY THERE IS AN OUTER LOOP AT ALL. `SolverOptions('NC')` sets
 * `config.highvar = 'interp'`, so a non-product-form FCFS station has its
 * service time rescaled by `npfqn_nonexp_approx` after each pass and the
 * analyzer re-solves until the decay rates eta stop moving. With no FCFS
 * station the reference pins iter_max to 1 and the loop runs once.
 *
 * ARITHMETIC. Every measure here is formed from differences of LOGARITHMS of
 * normalizing constants, which is the only way the ratios stay representable;
 * there is no exact-field formulation of exp(lG' - lG), so a non-transcendental
 * backend is refused by name rather than narrowed.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "line/api/npfqn/npfqn_nonexp_approx.h"
#include "line/api/pfqn/pfqn_nc.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/nc/solver_nc_lcfsqn.h"
#include "line/solvers/mva/mva_types.h"
#include "line/solvers/mva/sn_chain.h"
#include "line/solvers/nc/nc_types.h"
#include "line/solvers/nc/solver_nc_sdr.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace nc {

using lang::GlobalConstants;
using qn::SchedStrategy;

namespace detail {

/** `sn.njobs` per chain, infinite when any class of the chain is open. */
template <class T>
std::vector<double> chain_population(const qn::NetworkStruct<T>& sn) {
    std::vector<double> Nchain(sn.nchains, 0.0);
    for (std::size_t c = 0; c < sn.nchains; ++c) {
        double n = 0.0;
        bool open = false;
        for (std::size_t k : sn.inchain[c]) {
            const double p = sn.classes[k - 1].population;
            if (std::isinf(p)) open = true;
            n += p;
        }
        Nchain[c] = open ? std::numeric_limits<double>::infinity() : n;
    }
    return Nchain;
}

/**
 * The population vector `pfqn_nc` wants: integral counts, with an OPEN chain
 * marked by a negative entry.
 *
 * MATLAB marks an open class by N_r = Inf and zeroes it inside pfqn_nc; the
 * exact-capable port has no infinity and uses the sign instead, which is also
 * the .qn interchange convention.
 */
inline std::vector<int> nc_population(const std::vector<double>& Nchain) {
    std::vector<int> N(Nchain.size(), 0);
    for (std::size_t c = 0; c < Nchain.size(); ++c)
        N[c] = std::isinf(Nchain[c]) ? -1 : static_cast<int>(std::llround(Nchain[c]));
    return N;
}

/**
 * `cellsum(sn.visits)` at STATION level: the visit ratios summed over chains.
 *
 * The reference passes this to npfqn_nonexp_approx, which never reads it; it is
 * built anyway so the argument list matches the reference one for one.
 */
template <class T>
Matrix<T> station_visits(const qn::NetworkStruct<T>& sn) {
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> V(sn.nstations, sn.nclasses, zero);
    for (std::size_t c = 0; c < sn.nchains; ++c)
        for (std::size_t i = 0; i < sn.nstations; ++i) {
            const std::size_t sf = sn.stateful_of_station(i + 1);
            for (std::size_t k = 0; k < sn.nclasses; ++k)
                V(i, k) = T(V(i, k) + sn.visits[c](sf - 1, k));
        }
    return V;
}

/** `oner(N, r)`: one job of chain r removed. */
inline std::vector<int> oner(std::vector<int> N, std::size_t r) {
    if (N[r] > 0) --N[r];
    return N;
}

/** Map `options.method` onto the pfqn dispatcher, refusing an unknown name. */
inline pfqn::NcMethod nc_pfqn_method(const std::string& method) {
    return pfqn::nc_method_of(method);
}

}  // namespace detail

/**
 * Port of `solver_nc.m`.
 *
 * @param sn  the refreshed struct
 * @param opt solver controls
 * @return the class-level measures, the log normalizing constant and the
 *         concrete algorithm the constants were computed with
 */
template <class T>
NcSolution<T> solver_nc(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt) {
    NcSolution<T> out;
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)sn;
        (void)opt;
        throw UnsupportedError(
            "solver_nc: the normalizing-constant analyzer forms X = exp(lG(N-1_c) - lG(N)) and "
            "needs transcendental arithmetic; this backend has none");
    } else {
        // Krzesinski state-dependent routing: the model has its own product
        // form (eq. 16), so it is intercepted before the convolution and MVA
        // analyzers, which assume state-independent routing; see
        // _kb/16-state-dependent-routing.md
        if (!sn.sdr.empty()) return solver_nc_sdr(sn, opt);

        const T zero = num_traits<T>::from_int(0);
        const T one = num_traits<T>::from_int(1);
        const std::size_t M = sn.nstations, K = sn.nclasses, C = sn.nchains;

        std::vector<double> nservers(M, 1.0);
        std::vector<bool> isFCFS(M, false);
        bool anyFCFS = false;
        for (std::size_t i = 0; i < M; ++i) {
            nservers[i] = sn.stations[i].nservers;
            isFCFS[i] = sn.stations[i].sched == SchedStrategy::FCFS;
            if (isFCFS[i]) anyFCFS = true;
        }

        // LCFS paired with LCFS-PR is a two-station special case with its own
        // closed form; nothing else in this file can represent it, since a
        // non-preemptive LCFS station is not a BCMP station at all.
        std::vector<std::size_t> lcfsStat, lcfsprStat;
        for (std::size_t i = 0; i < M; ++i) {
            if (sn.stations[i].sched == SchedStrategy::LCFS) lcfsStat.push_back(i + 1);
            if (sn.stations[i].sched == SchedStrategy::LCFSPR) lcfsprStat.push_back(i + 1);
        }
        if (!lcfsStat.empty() && !lcfsprStat.empty()) {
            if (lcfsStat.size() != 1 || lcfsprStat.size() != 1)
                throw UnsupportedError(
                    "solver_nc: LCFS NC requires exactly one LCFS and one LCFS-PR station");
            for (std::size_t c = 0; c < C; ++c) {
                double nc_pop = 0.0;
                for (std::size_t r : sn.inchain[c]) nc_pop += sn.classes[r - 1].population;
                if (std::isinf(nc_pop))
                    throw UnsupportedError(
                        "solver_nc: LCFS NC requires a closed queueing network");
            }
            // The closed form places every job on one boundary-crossing sequence,
            // which a self-loop would let a job re-enter without crossing.
            for (std::size_t i : {lcfsStat[0], lcfsprStat[0]}) {
                const std::size_t sf = sn.stateful_of_station(i) - 1;
                for (std::size_t r = 0; r < K; ++r)
                    if (sn.rt(sf * K + r, sf * K + r) > zero)
                        throw UnsupportedError(
                            "solver_nc: LCFS NC does not support self-loops at stations");
            }
            out = solver_nc_lcfsqn(sn, opt, lcfsStat[0], lcfsprStat[0]);
            // STeff is the reference's `ST = 1./sn.rates` with NaN zeroed: an
            // INFINITE service time (a zero rate) is kept, only a NaN is not.
            out.STeff = Matrix<T>(M, K, zero);
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t r = 0; r < K; ++r) {
                    const double mu = num_traits<T>::to_double(sn.rates(i, r));
                    if (std::isnan(mu)) continue;
                    out.STeff(i, r) = mu == 0.0
                                          ? num_traits<T>::from_double(
                                                std::numeric_limits<double>::infinity())
                                          : T(one / sn.rates(i, r));
                }
            return out;
        }
        if (!lcfsStat.empty())
            throw UnsupportedError(
                "solver_nc: LCFS scheduling requires a paired LCFS-PR station");

        const std::vector<double> Nchain = detail::chain_population(sn);
        const std::vector<int> Nnc = detail::nc_population(Nchain);
        std::vector<std::size_t> openChains, closedChains;
        for (std::size_t c = 0; c < C; ++c)
            (std::isinf(Nchain[c]) ? openChains : closedChains).push_back(c);

        const pfqn::NcMethod pmethod = detail::nc_pfqn_method(opt.method);
        pfqn::NcOptions popt;
        popt.samples = opt.samples;
        popt.seed = opt.seed;
        popt.tol = opt.tol;
        popt.aghq_nodes = opt.aghq_nodes;
        popt.mcmc_batches = opt.mcmc_batches;
        popt.mcmc_burnin = opt.mcmc_burnin;
        const T atol = num_traits<T>::from_double(opt.tol);

        // Chain aggregation, and the class-level service times the outer loop
        // rescales. ST0 is the untouched original: npfqn_nonexp_approx is always
        // handed the original times, never its own previous output.
        mva::ChainDemands<T> d = mva::sn_get_demands_chain(sn);
        Matrix<T> ST = d.ST, ST0 = d.ST;
        const Matrix<T> V = detail::station_visits(sn);
        // sn.scv with the disabled pairs marked NOT FINITE, which is what
        // isfinite(SCV) selects on in the reference: sn.scv is NaN there, and
        // this port carries the marker in `disabled` instead.
        Matrix<T> SCVnan = sn.scv;
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t k = 0; k < K; ++k)
                if (sn.disabled[i][k])
                    SCVnan(i, k) = num_traits<T>::from_double(
                        std::numeric_limits<double>::quiet_NaN());

        std::vector<T> gamma(M, zero);
        std::vector<T> nserv_t(M, one);
        for (std::size_t i = 0; i < M; ++i) nserv_t[i] = num_traits<T>::from_double(nservers[i]);

        std::vector<T> lambda(C, zero);
        std::vector<T> eta(M, one), eta_1(M, zero);
        const int iter_max = anyFCFS ? opt.iter_max : 1;
        int it = 0;
        double lG = 0.0;
        std::string actualmethod = opt.method;
        mva::ClassResults<T> cls;

        while (it < iter_max) {
            {  // the reference's max(abs(1 - eta./eta_1)) > iter_tol test
                double dev = 0.0;
                for (std::size_t i = 0; i < M; ++i) {
                    const double e1 = num_traits<T>::to_double(eta_1[i]);
                    const double e = num_traits<T>::to_double(eta[i]);
                    const double v = std::fabs(1.0 - e / e1);
                    if (!(v <= dev)) dev = v;  // NaN and Inf both count as "not converged"
                }
                if (!(dev > opt.iter_tol)) break;
            }
            ++it;
            eta_1 = eta;

            if (it == 1) {
                // An open chain's reference station is its Source, whose chain
                // service time is 1 / total arrival rate.
                for (std::size_t c = 0; c < C; ++c) {
                    if (!std::isinf(Nchain[c])) continue;
                    const std::size_t rst = sn.classes[sn.inchain[c][0] - 1].refstat;
                    if (d.STchain(rst - 1, c) != zero)
                        lambda[c] = T(one / d.STchain(rst - 1, c));
                }
            } else {
                for (std::size_t c = 0; c < C; ++c)
                    for (std::size_t i = 0; i < M; ++i) {
                        T st = zero;
                        for (std::size_t k : sn.inchain[c]) st += ST(i, k - 1) * d.alpha(i, k - 1);
                        d.STchain(i, c) = st;
                        d.Lchain(i, c) = T(d.Vchain(i, c) * st);
                    }
            }
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t c = 0; c < C; ++c) {
                    if (!std::isfinite(num_traits<T>::to_double(d.STchain(i, c))))
                        d.STchain(i, c) = zero;
                    if (!std::isfinite(num_traits<T>::to_double(d.Lchain(i, c))))
                        d.Lchain(i, c) = zero;
                }

            // Seidmann's approximation: a c-server station becomes a demand L/c
            // at one server plus a pure delay L (c-1)/c.
            Matrix<T> Lms(M, C, zero), Z(M, C, zero), Zms(M, C, zero);
            std::vector<std::size_t> infServers;
            for (std::size_t i = 0; i < M; ++i) {
                if (std::isinf(nservers[i])) {
                    infServers.push_back(i);
                    for (std::size_t c = 0; c < C; ++c) Z(i, c) = d.Lchain(i, c);
                } else {
                    const T cs = num_traits<T>::from_double(nservers[i]);
                    for (std::size_t c = 0; c < C; ++c) {
                        Lms(i, c) = T(d.Lchain(i, c) / cs);
                        Zms(i, c) = T(d.Lchain(i, c) *
                                      num_traits<T>::from_double(nservers[i] - 1.0) / cs);
                    }
                }
            }
            Matrix<T> Ztot(1, C, zero);
            for (std::size_t c = 0; c < C; ++c)
                for (std::size_t i = 0; i < M; ++i) Ztot(0, c) += Z(i, c) + Zms(i, c);

            // step 1: the base constant
            const pfqn::NcDispatchResult<T> base =
                pfqn::pfqn_nc(lambda, Lms, Nnc, Ztot, pmethod, atol, popt);
            if (!base.valid) {
                // The method declined this model (the reference's empty lG).
                // MATLAB returns every output empty here and `getAvg` renders a
                // table of ZEROS while reporting a completed analysis; that
                // convention is reproduced deliberately, ruled by the user on
                // 2026-07-25 (register row N1). The consequence the ruling
                // accepts: a caller cannot tell "no jobs" from "declined".
                out.sol.Q = Matrix<T>(M, K, zero);
                out.sol.U = Matrix<T>(M, K, zero);
                out.sol.R = Matrix<T>(M, K, zero);
                out.sol.Tp = Matrix<T>(M, K, zero);
                out.sol.X.assign(K, zero);
                out.sol.C.assign(K, zero);
                out.sol.lG = 0.0;
                out.sol.iter = 1;
                out.sol.method = opt.method;
                out.actualmethod = opt.method;
                out.STeff = ST;
                return out;
            }
            lG = base.lG;
            actualmethod = base.method;
            std::vector<T> Xchain = base.X;
            Matrix<T> Qchain = base.Q;

            // Seidmann's delay term makes the by-product measures inconsistent
            // with the approximated model, so the reference discards them and
            // takes the constant-ratio route below. 'mcmc' is EXEMPT: it obtains
            // X and Q from ONE simulation of the regularized network, so
            // discarding them would cost R + M*R further simulations to recover
            // the same means by differencing lG. The surrogate delay is added
            // back to Qchain in the fill branch below instead.
            bool allZms = true;
            for (std::size_t c = 0; c < C; ++c) {
                T s = zero;
                for (std::size_t i = 0; i < M; ++i) s += Zms(i, c);
                if (!(num_traits<T>::to_double(s) > GlobalConstants::FineTol)) allZms = false;
            }
            if (allZms && pmethod != pfqn::NcMethod::Mcmc) {
                Xchain.clear();
                Qchain = Matrix<T>();
            }

            if (Xchain.empty()) {
                Xchain = lambda;
                Qchain = Matrix<T>(M, C, zero);
                for (std::size_t r : closedChains) {
                    const std::vector<int> Nr = detail::oner(Nnc, r);
                    const pfqn::NcDispatchResult<T> sub =
                        pfqn::pfqn_nc(lambda, Lms, Nr, Ztot, pmethod, atol, popt);
                    if (!sub.valid) {
                        // The reference returns empty from the SUBPROBLEMS too,
                        // rather than assigning [] into a scalar slot.
                        out.sol.Q = Matrix<T>(M, K, zero);
                        out.sol.U = Matrix<T>(M, K, zero);
                        out.sol.R = Matrix<T>(M, K, zero);
                        out.sol.Tp = Matrix<T>(M, K, zero);
                        out.sol.X.assign(K, zero);
                        out.sol.C.assign(K, zero);
                        out.sol.lG = 0.0;
                        out.sol.iter = 1;
                        out.sol.method = opt.method;
                        out.actualmethod = opt.method;
                        out.STeff = ST;
                        return out;
                    }
                    const double lGr = sub.lG;
                    Xchain[r] = num_traits<T>::from_double(std::exp(lGr - lG));
                    for (std::size_t i = 0; i < M; ++i) {
                        if (!(d.Lchain(i, r) > zero)) continue;
                        if (std::isinf(nservers[i])) {
                            Qchain(i, r) = T(d.Lchain(i, r) * Xchain[r]);
                            continue;
                        }
                        // Station i replicated into a private class holding one
                        // job: the arrival theorem in normalizing-constant form.
                        Matrix<T> Lar(M, C + 1, zero);
                        std::size_t row = 0;
                        for (std::size_t i2 = 0; i2 < M; ++i2) {
                            if (i2 == i) continue;
                            for (std::size_t c = 0; c < C; ++c) Lar(row, c) = Lms(i2, c);
                            ++row;
                        }
                        for (std::size_t c = 0; c < C; ++c) Lar(row, c) = Lms(i, c);
                        Lar(row, C) = one;
                        std::vector<T> lam_ar = lambda;
                        lam_ar.push_back(zero);
                        std::vector<int> N_ar = Nr;
                        N_ar.push_back(1);
                        Matrix<T> Z_ar(1, C + 1, zero);
                        for (std::size_t c = 0; c < C; ++c) Z_ar(0, c) = Ztot(0, c);
                        const pfqn::NcDispatchResult<T> ar =
                            pfqn::pfqn_nc(lam_ar, Lar, N_ar, Z_ar, pmethod, atol, popt);
                        // The most costly call of the whole solve, so it is the
                        // one whose algorithm is reported.
                        actualmethod = ar.method;
                        Qchain(i, r) = T(Zms(i, r) * Xchain[r] +
                                         Lms(i, r) * num_traits<T>::from_double(
                                                         std::exp(ar.lG - lG)));
                    }
                }
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t c = 0; c < C; ++c)
                        if (std::isnan(num_traits<T>::to_double(Qchain(i, c)))) Qchain(i, c) = zero;

                // An open chain sees the closed population as extra load: the
                // mixed-network queue length of Bruell and Balbo.
                for (std::size_t r : openChains)
                    for (std::size_t i = 0; i < M; ++i) {
                        T load = zero;
                        for (std::size_t o : openChains) load += lambda[o] * d.Lchain(i, o);
                        const T den =
                            std::isinf(nservers[i])
                                ? one
                                : T(one - load / num_traits<T>::from_double(nservers[i]));
                        if (den == zero) continue;
                        T qc = zero;
                        for (std::size_t cc : closedChains) qc += Qchain(i, cc);
                        Qchain(i, r) = T(lambda[r] * d.Lchain(i, r) / den * (one + qc));
                    }
            } else {
                // The method returned the measures, so they are KEPT: only the delay
                // stations, which it does not model, are filled in, and the population
                // Seidmann's surrogate delay holds outside a multiserver queueing station
                // is restored. Zms is zero at a single-server station, so the second arm
                // is a no-op on a single-server model. Rebuilding Qchain from zero here
                // instead would throw away the queueing-station queue lengths the method
                // just produced, which is what the reference keeps.
                if (Qchain.rows() != M || Qchain.cols() != C) Qchain = Matrix<T>(M, C, zero);
                for (std::size_t c = 0; c < C; ++c)
                    for (std::size_t i = 0; i < M; ++i) {
                        if (!(d.Lchain(i, c) > zero)) continue;
                        if (std::isinf(nservers[i]))
                            Qchain(i, c) = T(d.Lchain(i, c) * Xchain[c]);
                        else if (nservers[i] > 1.0)
                            Qchain(i, c) = T(Qchain(i, c) + Zms(i, c) * Xchain[c]);
                    }
            }

            Matrix<T> Rchain(M, C, zero), Tchain(M, C, zero);
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t c = 0; c < C; ++c) {
                    if (Xchain[c] != zero && d.Vchain(i, c) != zero)
                        Rchain(i, c) = T(Qchain(i, c) / Xchain[c] / d.Vchain(i, c));
                    Tchain(i, c) = T(Xchain[c] * d.Vchain(i, c));
                }
            for (std::size_t i : infServers)
                for (std::size_t c = 0; c < C; ++c)
                    Rchain(i, c) = d.Vchain(i, c) == zero
                                       ? zero
                                       : T(d.Lchain(i, c) / d.Vchain(i, c));

            d.ST = ST;
            cls = mva::sn_deaggregate_chain_results(sn, d, Matrix<T>(), Matrix<T>(), Rchain,
                                                    Tchain, Xchain);
            out.STeff = ST;  // the effective service times of THIS pass

            const npfqn::NonexpApproxResult<T> na = npfqn::npfqn_nonexp_approx(
                opt.highvar, isFCFS, sn.rates, ST0, V, SCVnan, cls.Tp, cls.U, gamma, nserv_t);
            ST = na.ST;
            gamma = na.gamma;
            eta = na.eta;
        }

        if (it == 0)
            throw UnsupportedError(
                "solver_nc: the analyzer made no pass; iter_max must be at least one");

        Matrix<T> Q = cls.Q, U = cls.U, R = cls.R, Tp = cls.Tp;
        std::vector<T> X = cls.X;
        auto sanitize = [&](Matrix<T>& A) {
            for (std::size_t i = 0; i < A.rows(); ++i)
                for (std::size_t j = 0; j < A.cols(); ++j) {
                    if (A(i, j) < zero) A(i, j) = T(-A(i, j));
                    if (!std::isfinite(num_traits<T>::to_double(A(i, j)))) A(i, j) = zero;
                }
        };
        sanitize(Q);
        sanitize(U);
        sanitize(R);
        for (T& x : X) {
            if (x < zero) x = T(-x);
            if (!std::isfinite(num_traits<T>::to_double(x))) x = zero;
        }

        // Renormalize: the approximations above do not conserve the closed
        // population exactly, and a chain whose queue lengths do not add up to
        // N would report a throughput off by the same factor.
        for (std::size_t c = 0; c < C; ++c) {
            if (std::isinf(Nchain[c])) continue;
            T qden = zero;
            for (std::size_t k : sn.inchain[c])
                for (std::size_t i = 0; i < M; ++i) qden += Q(i, k - 1);
            const T ratio = qden > zero
                                ? T(num_traits<T>::from_double(Nchain[c]) / qden)
                                : zero;
            for (std::size_t k : sn.inchain[c]) {
                X[k - 1] = T(ratio * X[k - 1]);
                for (std::size_t i = 0; i < M; ++i) {
                    Q(i, k - 1) = T(ratio * Q(i, k - 1));
                    Tp(i, k - 1) = T(ratio * Tp(i, k - 1));
                    U(i, k - 1) = T(ratio * U(i, k - 1));
                    R(i, k - 1) = Tp(i, k - 1) == zero ? zero : T(Q(i, k - 1) / Tp(i, k - 1));
                }
            }
        }

        out.sol.Q = Q;
        out.sol.U = U;
        out.sol.R = R;
        out.sol.Tp = Tp;
        out.sol.X = X;
        out.sol.C.assign(K, zero);
        for (std::size_t k = 0; k < K; ++k) {
            const double njobs = sn.classes[k].population;
            if (std::isfinite(njobs) && X[k] != zero)
                out.sol.C[k] = T(num_traits<T>::from_double(njobs) / X[k]);
            else
                out.sol.C[k] = cls.C[k];
        }
        out.sol.lG = lG;
        out.sol.iter = it;
        out.sol.method = actualmethod;
        out.actualmethod = actualmethod;
        return out;
    }
}

}  // namespace nc
}  // namespace line

#endif  // LINE_SOLVERS_NC_SOLVER_NC_H
