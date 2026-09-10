/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MVA_SOLVER_QNA_H
#define LINE_SOLVERS_MVA_SOLVER_QNA_H

/**
 * QNA, the two-moment open-network decomposition analyzer.
 *
 * Templated port of `matlab/src/solvers/MVA/solver_qna.m`, which implements
 * Section 7.2.3 of N. Gautam, "Analysis of Queues", CRC Press 2012, with
 * corrections the reference records as agreed with the author.
 *
 * The structure is decomposition-aggregation. Each sweep
 *
 *   1. SUPERPOSES the flows into every station: a1(i,r) is the arrival rate of
 *      class r at station i and a2(i,r) the squared coefficient of variation of
 *      that superposed stream, both accumulated over the routing;
 *   2. SOLVES each station in isolation -- an infinite server passes its input
 *      SCV straight through, a PS station takes a geometric-bound queue length,
 *      and an FCFS station takes the multiserver GI/G/c waiting-time
 *      approximation with the Whitt alpha(m) correction;
 *   3. SPLITS the departure stream along each outgoing arc, `f2 = 1 + p (d2-k)`
 *      with k the round-robin split degree (`npfqn_traffic_split_rr`): at k = 1
 *      this is the exact SCV of a Bernoulli-thinned renewal stream, and at k > 1
 *      that of a one-in-k deterministic dispatch, which is less variable;
 *
 * and `da_fpi` drives the sweeps to a fixed point on the queue lengths.
 *
 * SCOPE. OPEN CHAINS ONLY. The reference reaches `State.toMarginal` to seed a
 * closed chain's queue length (`solver_qna.m:94`), which needs a state-encoding
 * layer this port does not have, and `SolverMVA.listValidMethods` advertises
 * `qna` only for a fully open model in the first place. A closed chain is
 * refused by name rather than seeded with a guess.
 *
 * REFERENCE INDEXING, REPRODUCED AND CHECKED. `solver_qna.m` indexes `sn.rt`,
 * which is STATEFUL-indexed, with STATION indices, and reads `sn.njobs`, which
 * is CLASS-indexed, with a CHAIN index. Both are silently correct exactly when
 * the stations coincide with the stateful nodes and the chains with the
 * classes, which is the shape of the models the method is advertised for. This
 * port checks that instead of assuming it, and refuses by name otherwise -- the
 * alternative is reading an unrelated row and returning a number.
 *
 * SELF-LOOPING CLASSES. The reference special-cases `sn.isslc` in three places.
 * This port has no SelfLoopingClass (`JobClassType` is OPEN or CLOSED only), so
 * those branches are unreachable and are not transcribed; a class type this
 * port cannot build cannot reach here.
 *
 * Arithmetic: TRANSCENDENTAL. `da_fpi` stops on a tolerance, and the alpha(m)
 * correction takes a real power of the utilization.
 */

#include <cmath>
#include <string>
#include <vector>

#include "line/api/da/da_fpi.h"
#include "line/api/da/da_traffic_superpos.h"
#include "line/api/npfqn/npfqn_traffic_split_rr.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mva/mva_types.h"

namespace line {
namespace mva {

/**
 * Port of `solver_qna.m`.
 *
 * @param L    the model; open chains only
 * @param opt  tol drives the saturation test, iter_max / iter_tol the sweeps
 */
template <class T>
MvaSolution<T> solver_qna(const qn::NetworkStruct<T>& L, const MvaOptions& opt) {
    using qn::NodeType;
    using qn::SchedStrategy;
    if constexpr (!num_traits<T>::has_transcendental) {
        (void)L;
        (void)opt;
        throw UnsupportedError(
            "solver_qna: the decomposition sweeps stop on a tolerance and the multiserver "
            "correction takes a real power, so exact rational arithmetic cannot run it");
    } else {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations, K = L.nclasses, C = L.nchains;

    for (const auto& c : L.classes)
        if (std::isfinite(c.population))
            throw UnsupportedError(
                "solver_qna: QNA is an open-network decomposition; a closed chain needs the "
                "state-encoding layer the reference seeds it from, which is not ported");
    if (L.nof_stateful() != M)
        throw UnsupportedError(
            "solver_qna: the reference indexes the stateful-indexed sn.rt with station indices, "
            "which is only correct when every stateful node is a station; this model has " +
            std::to_string(L.nof_stateful()) + " stateful nodes and " + std::to_string(M) +
            " stations");
    if (C != K)
        throw UnsupportedError(
            "solver_qna: the reference reads the class-indexed sn.njobs with a chain index, which "
            "is only correct when each chain holds exactly one class");
    // One predicate for the gate and the run: qn::mva_feature_set withholds 'qna'
    // for a discipline the station update has no arm for, and the update refuses
    // it by name rather than leaving that station's row of Q, U, R and T at zero
    // and returning the table as a solution.
    {
        const std::string qna_reason = mva_qna_scheduling_reason(L);
        if (!qna_reason.empty()) throw UnsupportedError(qna_reason);
    }

    // S = 1/rates, with a NaN scv read as zero, as the reference does
    Matrix<T> S(M, K, zero), scv(M, K, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r) {
            const T mu = L.rates(i, r);
            S(i, r) = (mu > zero) ? T(one / mu) : zero;
            const double v = num_traits<T>::to_double(L.scv(i, r));
            scv(i, r) = std::isnan(v) ? zero : L.scv(i, r);
        }

    // V(i,k): the visits summed over chains, at station level
    Matrix<T> V(M, K, zero);
    for (std::size_t c = 0; c < C; ++c)
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t k = 0; k < K; ++k)
                V(i, k) = T(V(i, k) + L.visits[c](L.stateful_of_station(i + 1) - 1, k));

    const Matrix<T>& rt = L.rt;
    auto RT = [&](std::size_t i, std::size_t r, std::size_t j, std::size_t s) -> T {
        return rt(i * K + r, j * K + s);
    };

    // The deterministic (round-robin) split degrees, 1 where the split is
    // Markovian; a flow thinned out of a one-in-k dispatch is the k-fold
    // convolution, whose SCV at d2 = 1 is 1 + p (1 - k), below the renewal 1.
    const Matrix<T> kRR = npfqn::npfqn_traffic_split_rr(L);

    // f2(i,r -> j,s): the SCV of each flow, initialised on every arc that does
    // not end at a Source
    Matrix<T> f2(M * K, M * K, zero);
    std::vector<bool> is_source(M, false);
    for (std::size_t i = 0; i < M; ++i)
        is_source[i] = (L.stations[i].nodetype == NodeType::Source);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < M; ++j) {
            if (is_source[j]) continue;
            for (std::size_t r = 0; r < K; ++r)
                for (std::size_t s = 0; s < K; ++s)
                    if (RT(i, r, j, s) > zero)
                        f2(i * K + r, j * K + s) =
                            T(one + RT(i, r, j, s) * T(one - kRR(i, r)));
        }

    // per-chain arrival rate and the SCV of the superposed source stream
    std::vector<T> lambda(C, zero), d2c(C, zero);
    Matrix<T> Tp(M, K, zero), Q(M, K, zero), U(M, K, zero), Rt(M, K, zero);
    for (std::size_t c = 0; c < C; ++c) {
        const std::size_t ref = L.classes[L.inchain[c][0] - 1].refstat;
        std::vector<T> lam_in, scv_in;
        for (std::size_t k : L.inchain[c]) {
            lam_in.push_back(L.rates(ref - 1, k - 1));
            scv_in.push_back(scv(ref - 1, k - 1));
        }
        T s = zero;
        for (const T& x : lam_in)
            if (std::isfinite(num_traits<T>::to_double(x))) s += x;
        lambda[c] = s;
        d2c[c] = da::da_traffic_superpos(lam_in, scv_in);
        for (std::size_t a = 0; a < L.inchain[c].size(); ++a)
            Tp(ref - 1, L.inchain[c][a] - 1) = lam_in[a];
    }
    // d2(i): the SCV of the departure stream of station i. Every reference
    // station is seeded with the SAME rate-weighted mean over the chains, which
    // is what the reference's `d2(refStatIdx) = d2c*lambda'/sum(lambda)` writes
    // once per chain.
    std::vector<T> d2(M, zero);
    {
        T num = zero, den = zero;
        for (std::size_t c = 0; c < C; ++c) {
            num += T(d2c[c] * lambda[c]);
            den += lambda[c];
        }
        const T seed = (den > zero) ? T(num / den) : zero;
        for (std::size_t c = 0; c < C; ++c)
            d2[L.classes[L.inchain[c][0] - 1].refstat - 1] = seed;
    }

    Matrix<T> a1(M, K, zero), a2(M, K, zero);
    const T tol = num_traits<T>::from_double(opt.tol);

    // One decomposition sweep. `da_fpi` iterates it on the flattened queue
    // lengths; every other quantity is carried in the enclosing scope, as the
    // reference carries it in the enclosing function's workspace.
    auto sweep = [&](const std::vector<T>& qin,
                     std::size_t itnum) -> std::pair<std::vector<T>, std::vector<T>> {
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t k = 0; k < K; ++k) Q(i, k) = qin[i * K + k];

        // THE RENORMALISATION IS LOAD BEARING, AND IT IS WHAT ENDS THE LOOP.
        // The reference rescales each chain's queue lengths to its population,
        // `Q(:,c) = njobs(c) * Q(:,c) / sum(Q(:,c))`. On an OPEN chain njobs is
        // infinite and the initial iterate is zero, so this is Inf * 0 / 0 =
        // NaN; the reference point handed back to da_fpi is therefore NaN, the
        // convergence measure is NaN, and `da_nanstop` -- which the reference
        // sets precisely to reproduce its legacy while-loop -- stops the
        // iteration after ONE sweep. QNA as the reference ships it does not
        // iterate its fixed point on an open model, and the numbers it reports
        // are the first sweep's, computed with every flow SCV still at its
        // initial 1. Reproduced rather than corrected: skipping the NaN lets
        // the sweeps run to their actual fixed point and changes every result
        // (Q1 = 0.75 instead of the reference's 0.875 on the Erlang tandem).
        std::vector<T> qref(M * K, zero);
        for (std::size_t c = 0; c < C; ++c) {
            const double nc = L.classes[c].population;
            T colsum = zero;
            for (std::size_t i = 0; i < M; ++i) colsum += Q(i, c);
            for (std::size_t i = 0; i < M; ++i)
                Q(i, c) = num_traits<T>::from_double(nc * num_traits<T>::to_double(Q(i, c)) /
                                                     num_traits<T>::to_double(colsum));
        }
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t k = 0; k < K; ++k) qref[i * K + k] = Q(i, k);

        if (itnum == 1)
            for (std::size_t c = 0; c < C; ++c)
                for (std::size_t m = 0; m < M; ++m)
                    for (std::size_t k : L.inchain[c])
                        Tp(m, k - 1) = T(V(m, k - 1) * lambda[c]);

        // ---- superposition ------------------------------------------------
        for (std::size_t i = 0; i < M; ++i) {
            T lambda_i = zero;
            for (std::size_t k = 0; k < K; ++k) lambda_i += Tp(i, k);
            for (std::size_t r = 0; r < K; ++r) {
                a1(i, r) = zero;
                a2(i, r) = zero;
            }
            for (std::size_t j = 0; j < M; ++j)
                for (std::size_t r = 0; r < K; ++r)
                    for (std::size_t s = 0; s < K; ++s) {
                        const T p = RT(j, s, i, r);
                        if (!(p > zero)) continue;
                        a1(i, r) = T(a1(i, r) + Tp(j, s) * p);
                        if (lambda_i > zero)
                            a2(i, r) = T(a2(i, r) + T(one / lambda_i) * f2(j * K + s, i * K + r) *
                                                        Tp(j, s) * p);
                    }
        }

        // ---- solve each station in isolation ------------------------------
        for (std::size_t i = 0; i < M; ++i) {
            if (L.stations[i].nodetype == NodeType::Join) continue;  // no-op, as in the reference
            switch (L.stations[i].sched) {
                case SchedStrategy::INF: {
                    // departure SCV of a delay: MATLAB writes d2(ist,s)=a2(ist,s)
                    // but every downstream read is the scalar d2(ist), i.e. the
                    // FIRST class column a2(ist,1); JAR reads the same. Use a2(i,0).
                    d2[i] = a2(i, 0);
                    for (std::size_t k = 0; k < K; ++k) {
                        Tp(i, k) = a1(i, k);
                        Q(i, k) = T(Tp(i, k) * S(i, k) * V(i, k));
                        U(i, k) = Q(i, k);
                        Rt(i, k) = Tp(i, k) > zero ? T(Q(i, k) / Tp(i, k)) : zero;
                    }
                    break;
                }
                case SchedStrategy::PS: {
                    for (std::size_t c = 0; c < C; ++c) {
                        for (std::size_t k : L.inchain[c]) {
                            Tp(i, k - 1) = T(lambda[c] * V(i, k - 1));
                            U(i, k - 1) = T(S(i, k - 1) * Tp(i, k - 1));
                        }
                        T usum = zero;
                        for (std::size_t k = 0; k < K; ++k) usum += U(i, k);
                        const T uden = usum < T(one - tol) ? usum : T(one - tol);
                        // Nc is infinite on an open chain, so U^(Nc+1) vanishes
                        // for a stable station and the geometric bound is
                        // U/(1-Uden); the reference writes the finite-population
                        // form, whose limit this is.
                        for (std::size_t k : L.inchain[c]) {
                            Q(i, k - 1) = T(one - uden) > zero
                                              ? T(U(i, k - 1) / T(one - uden))
                                              : zero;
                            Rt(i, k - 1) =
                                Tp(i, k - 1) > zero ? T(Q(i, k - 1) / Tp(i, k - 1)) : zero;
                        }
                    }
                    break;
                }
                case SchedStrategy::FCFS: {
                    using std::pow;
                    using std::sqrt;
                    const T ftol = num_traits<T>::from_double(GlobalConstants::FineTol);
                    std::vector<T> mu(K, zero), rho_cls(K, zero);
                    T lambda_i = zero;
                    for (std::size_t r = 0; r < K; ++r) {
                        const double m = num_traits<T>::to_double(L.rates(i, r));
                        mu[r] = std::isnan(m) ? zero : L.rates(i, r);
                        const T den = T(ftol + mu[r]);
                        rho_cls[r] = den > zero ? T(a1(i, r) / den) : zero;
                        if (std::isnan(num_traits<T>::to_double(rho_cls[r]))) rho_cls[r] = zero;
                        lambda_i += a1(i, r);
                    }
                    const T mi = num_traits<T>::from_double(L.stations[i].nservers);
                    T rho = zero;
                    for (std::size_t r = 0; r < K; ++r) rho += rho_cls[r];
                    rho = T(rho / mi);
                    if (rho < T(one - tol)) {
                        // Whitt's alpha(m): the heavier form above rho = 0.7
                        const T alpha = num_traits<T>::to_double(rho) > 0.7
                                            ? T(T(pow(rho, mi) + rho) / num_traits<T>::from_int(2))
                                            : T(pow(rho, T(T(mi + one) / num_traits<T>::from_int(2))));
                        const T mubar = rho > zero ? T(lambda_i / rho) : zero;
                        T c2 = T(-one);
                        for (std::size_t r = 0; r < K; ++r) {
                            if (!(mu[r] > zero) || !(lambda_i > zero)) continue;
                            const T q = T(mubar / mi / mu[r]);
                            c2 += T(a1(i, r) / lambda_i * q * q * T(scv(i, r) + one));
                        }
                        T a2sum = zero;
                        for (std::size_t r = 0; r < K; ++r) a2sum += a2(i, r);
                        const T Wiq = mubar > zero
                                          ? T(T(alpha / mubar) * T(one / T(one - rho)) *
                                              T(T(a2sum + c2) / num_traits<T>::from_int(2)))
                                          : zero;
                        for (std::size_t k = 0; k < K; ++k)
                            Q(i, k) = mu[k] > zero ? T(a1(i, k) / mu[k] + a1(i, k) * Wiq) : zero;
                        d2[i] = T(one + T(rho * rho * T(c2 - one) / sqrt(mi)) +
                                  T(T(one - rho * rho) * T(a2sum - one)));
                    } else {
                        // saturated: the reference parks the queue length at the
                        // class population, which is infinite on an open class
                        for (std::size_t k = 0; k < K; ++k)
                            Q(i, k) = num_traits<T>::from_double(L.classes[k].population);
                        d2[i] = one;
                    }
                    for (std::size_t k = 0; k < K; ++k) {
                        Tp(i, k) = a1(i, k);
                        U(i, k) = T(Tp(i, k) * S(i, k) / mi);
                        Rt(i, k) = Tp(i, k) > zero ? T(Q(i, k) / Tp(i, k)) : zero;
                    }
                    break;
                }
                default:
                    // EXT (a Source or a Sink) carries no queue; every other
                    // discipline is one the reference's switch does not handle,
                    // and falling through would leave the station unsolved.
                    if (L.stations[i].sched != SchedStrategy::EXT)
                        throw UnsupportedError(
                            std::string("solver_qna: no isolated-station solution for ") +
                            lang::sched_to_text(L.stations[i].sched) + " scheduling");
                    break;
            }
        }

        // ---- splitting ----------------------------------------------------
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t j = 0; j < M; ++j) {
                if (is_source[j]) continue;
                for (std::size_t r = 0; r < K; ++r)
                    for (std::size_t s = 0; s < K; ++s)
                        if (RT(i, r, j, s) > zero)
                            // k-fold convolution then Bernoulli thinning at
                            // q = k p: C^2 = 1 + p (d2 - k), the Markovian
                            // 1 + p (d2 - 1) at k = 1
                            f2(i * K + r, j * K + s) =
                                T(one + RT(i, r, j, s) * T(d2[i] - kRR(i, r)));
            }

        std::vector<T> qnew(M * K, zero);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t k = 0; k < K; ++k) qnew[i * K + k] = Q(i, k);
        return std::make_pair(qnew, qref);
    };

    da::FpiOptions fo;
    // the legacy while-loop ran one extra sweep at the cap, and exited on a
    // non-finite convergence measure
    fo.iter_max = static_cast<std::size_t>(opt.iter_max) + 1;
    fo.iter_tol = opt.iter_tol;
    fo.nanstop = true;
    const da::FpiResult<T> fr = da::da_fpi<T>(sweep, std::vector<T>(M * K, zero), fo);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k) Q(i, k) = fr.x[i * K + k];

    // an infinite server's utilization IS its queue length
    for (std::size_t i = 0; i < M; ++i)
        if (L.stations[i].sched == SchedStrategy::INF)
            for (std::size_t k = 0; k < K; ++k) U(i, k) = Q(i, k);

    MvaSolution<T> out;
    out.Q = Q;
    out.U = U;
    out.R = Rt;
    out.Tp = Tp;
    out.C.assign(K, zero);
    out.X.assign(K, zero);
    for (std::size_t k = 0; k < K; ++k)
        for (std::size_t i = 0; i < M; ++i) out.C[k] += Rt(i, k);
    // the reference takes |Q| and maps every NaN to zero on the way out
    auto clean = [](Matrix<T>& A) {
        for (std::size_t i = 0; i < A.rows(); ++i)
            for (std::size_t j = 0; j < A.cols(); ++j)
                if (std::isnan(num_traits<T>::to_double(A(i, j))))
                    A(i, j) = num_traits<T>::from_int(0);
    };
    for (std::size_t i = 0; i < out.Q.rows(); ++i)
        for (std::size_t j = 0; j < out.Q.cols(); ++j)
            if (out.Q(i, j) < zero) out.Q(i, j) = T(-out.Q(i, j));
    clean(out.Q);
    clean(out.U);
    clean(out.R);
    for (std::size_t k = 0; k < K; ++k)
        if (std::isnan(num_traits<T>::to_double(out.C[k]))) out.C[k] = zero;
    out.method = "qna";
    out.iter = static_cast<int>(fr.iterations);
    out.lG = 0.0;
    return out;
    }
}

}  // namespace mva
}  // namespace line

#endif  // LINE_SOLVERS_MVA_SOLVER_QNA_H
