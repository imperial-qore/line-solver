/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MVA_SN_CHAIN_H
#define LINE_SOLVERS_MVA_SN_CHAIN_H

/**
 * Chain aggregation and de-aggregation.
 *
 * Ports of matlab/src/api/sn/sn_get_demands_chain.m,
 * sn_get_product_form_chain_params.m and sn_deaggregate_chain_results.m. These
 * are the bridge between the class-level layer struct and the chain-level
 * product-form solvers: a chain is a set of classes a job moves between while
 * circulating, so the solvers work on chains and the results are pushed back to
 * classes with the visit-share weights `alpha`.
 *
 * NON-FINITE VALUES. MATLAB reaches this code with NaN and Inf in `sn.rates`
 * (NaN for a station-class pair the class never visits, Inf for an immediate
 * service) and clears them afterwards with `X(~isfinite(X)) = 0`. The exact
 * backend has neither value, so the two are carried as explicit flags -- the
 * layer's `disabled` mask, and a local `st_inf` mask for a zero rate -- and the
 * same entries end up zero. Every place a MATLAB `isfinite` guard fires, the
 * corresponding flag test fires here.
 */

#include <algorithm>
#include <cmath>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mva {

using lang::GlobalConstants;
using lang::SchedStrategy;

/** The chain-level view of a layer, as sn_get_demands_chain returns it. */
template <class T>
struct ChainDemands {
    Matrix<T> Lchain;              ///< (M x C) demand
    Matrix<T> STchain;             ///< (M x C) mean service time
    Matrix<T> Vchain;              ///< (M x C) visits
    Matrix<T> alpha;               ///< (M x K) class share of its chain's visits at a station
    std::vector<double> Nchain;    ///< (C) population, infinite for an open chain
    Matrix<T> SCVchain;            ///< (M x C)
    std::vector<std::size_t> refstatchain;  ///< (C) 1-based reference station
    Matrix<T> ST;                  ///< (M x K) class-level mean service time, 0 where disabled
};

/**
 * Port of sn_get_demands_chain.
 *
 * @param L the refreshed layer
 */
template <class T>
ChainDemands<T> sn_get_demands_chain(const qn::NetworkStruct<T>& L) {
    const std::size_t M = L.nstations, K = L.nclasses, C = L.nchains;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    ChainDemands<T> d;
    d.ST = Matrix<T>(M, K, zero);
    std::vector<std::vector<bool>> st_inf(M, std::vector<bool>(K, false));
    Matrix<T> SCV(M, K, one);
    std::vector<std::vector<bool>> scv_known(M, std::vector<bool>(K, true));
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r) {
            if (L.disabled[i][r]) {
                // MATLAB: ST = 1./rates gives NaN, then ST(isnan(ST)) = 0;
                // SCV(isnan(SCV)) = 1
                d.ST(i, r) = zero;
                SCV(i, r) = one;
                continue;
            }
            if (L.rates(i, r) == zero) {
                // 1/0 = Inf in MATLAB; every use is guarded by isfinite and
                // ends up zeroed, so the flag records it and the value stays 0
                st_inf[i][r] = true;
                d.ST(i, r) = zero;
            } else {
                d.ST(i, r) = T(one / L.rates(i, r));
            }
            SCV(i, r) = L.scv(i, r);
        }

    d.alpha = Matrix<T>(M, K, zero);
    d.Vchain = Matrix<T>(M, C, zero);
    std::vector<std::vector<bool>> v_bad(M, std::vector<bool>(C, false));

    for (std::size_t c = 0; c < C; ++c) {
        const std::vector<std::size_t>& ic = L.inchain[c];
        const std::size_t rstat = L.classes[ic[0] - 1].refstat;
        // denominator: the reference class's visits when the chain has one,
        // otherwise the whole chain's visits at the reference station
        T den = zero;
        const std::size_t rsf = L.stateful_of_station(rstat) - 1;
        if (L.refclass[c] > 0) {
            den = L.visits[c](rsf, L.refclass[c] - 1);
        } else {
            for (std::size_t k : ic) den += L.visits[c](rsf, k - 1);
        }
        for (std::size_t i = 0; i < M; ++i) {
            const std::size_t sf = L.stateful_of_station(i + 1) - 1;
            T num = zero;
            for (std::size_t k : ic) num += L.visits[c](sf, k - 1);
            if (den == zero)
                v_bad[i][c] = true;  // MATLAB: 0/0 or x/0 -> NaN or Inf, then zeroed
            else
                d.Vchain(i, c) = T(num / den);
            if (num != zero)
                for (std::size_t k : ic)
                    d.alpha(i, k - 1) = T(d.alpha(i, k - 1) + L.visits[c](sf, k - 1) / num);
        }
    }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t c = 0; c < C; ++c)
            if (v_bad[i][c]) d.Vchain(i, c) = zero;

    // renormalise each chain's visits by the value at its reference station
    for (std::size_t c = 0; c < C; ++c) {
        const std::size_t rstat = L.classes[L.inchain[c][0] - 1].refstat;
        const T vref = d.Vchain(rstat - 1, c);
        if (vref == zero) continue;  // MATLAB divides by 0 and clears below
        for (std::size_t i = 0; i < M; ++i) d.Vchain(i, c) = T(d.Vchain(i, c) / vref);
    }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k)
            if (num_traits<T>::to_double(d.alpha(i, k)) < GlobalConstants::Zero)
                d.alpha(i, k) = zero;

    d.Lchain = Matrix<T>(M, C, zero);
    d.STchain = Matrix<T>(M, C, zero);
    d.SCVchain = Matrix<T>(M, C, zero);
    d.Nchain.assign(C, 0.0);
    d.refstatchain.assign(C, 1);

    for (std::size_t c = 0; c < C; ++c) {
        const std::vector<std::size_t>& ic = L.inchain[c];
        bool open = false;
        double n = 0.0;
        for (std::size_t k : ic) {
            const double p = L.classes[k - 1].population;
            if (std::isinf(p)) open = true;
            n += p;
        }
        d.Nchain[c] = open ? std::numeric_limits<double>::infinity() : n;
        const std::size_t rstat = L.classes[ic[0] - 1].refstat;
        d.refstatchain[c] = rstat;

        for (std::size_t i = 0; i < M; ++i) {
            bool bad = false;
            T st = zero;
            for (std::size_t k : ic) {
                if (st_inf[i][k - 1] && d.alpha(i, k - 1) != zero) bad = true;
                st += d.ST(i, k - 1) * d.alpha(i, k - 1);
            }
            if (open && i + 1 == rstat) {
                // a source station: the chain service time is 1 / total arrival rate
                T lam = zero;
                for (std::size_t k : ic)
                    if (!L.disabled[i][k - 1]) lam += L.rates(i, k - 1);
                if (lam == zero)
                    bad = true;
                else
                    st = T(num_traits<T>::from_int(1) / lam);
            }
            d.STchain(i, c) = bad ? zero : st;
            d.Lchain(i, c) = T(d.Vchain(i, c) * d.STchain(i, c));

            T ach = zero;
            for (std::size_t k : ic) ach += d.alpha(i, k - 1);
            if (ach > zero) {
                T s = zero;
                for (std::size_t k : ic) s += SCV(i, k - 1) * d.alpha(i, k - 1);
                d.SCVchain(i, c) = T(s / ach);
            }
        }
    }
    return d;
}

/** Class-level results, as sn_deaggregate_chain_results returns them. */
template <class T>
struct ClassResults {
    Matrix<T> Q, U, R, Tp;
    std::vector<T> C, X;
};

/**
 * Port of sn_deaggregate_chain_results.
 *
 * OPEN CHAINS. A closed class's share of its chain's throughput is its share of
 * the visits at the reference station; an open class's is its share of the
 * DEPARTURES, i.e. the chain's visit count at the Sink, because an open job
 * completes by leaving rather than by returning. Only X and C are computed this
 * way; Q, U, R and T follow the same path for both.
 *
 * The Sink visits come from this layer's own `nodevisits`, normalised by the
 * chain's reference NODE. MATLAB normalises the node-level visits by
 * `statefulToNode(refstat(...))`, which feeds a STATION index to a STATEFUL ->
 * node map and lands on an unrelated node, so its Sink visits carry an
 * arbitrary scale factor (3 rather than 1 on the fork layer of lqn_workflows).
 * The correct normalisation is used here. It reaches X and C only, which
 * SolverLN does not read back, so no layered metric depends on the choice.
 */
template <class T>
ClassResults<T> sn_deaggregate_chain_results(const qn::NetworkStruct<T>& L, const ChainDemands<T>& d,
                                             const Matrix<T>& Qchain, const Matrix<T>& Uchain,
                                             const Matrix<T>& Rchain, const Matrix<T>& Tchain,
                                             const std::vector<T>& Xchain) {
    const std::size_t M = L.nstations, K = L.nclasses;
    const T zero = num_traits<T>::from_int(0);
    ClassResults<T> r;
    r.Q = Matrix<T>(M, K, zero);
    r.U = Matrix<T>(M, K, zero);
    r.R = Matrix<T>(M, K, zero);
    r.Tp = Matrix<T>(M, K, zero);
    r.C.assign(K, zero);
    r.X.assign(K, zero);

    std::vector<T> Vsink(K, zero);
    if (L.sinkNode > 0)
        for (std::size_t c = 0; c < L.nchains; ++c)
            for (std::size_t k = 0; k < K; ++k) Vsink[k] += L.nodevisits[c](L.sinkNode - 1, k);

    for (std::size_t c = 0; c < L.nchains; ++c) {
        const std::vector<std::size_t>& ic = L.inchain[c];
        const bool open = std::isinf(d.Nchain[c]);
        if (open && L.sinkNode == 0)
            throw UnsupportedError(
                "sn_deaggregate_chain_results: an open chain needs the Sink node visits, and this "
                "layer has no Sink node to read them from");
        for (std::size_t kk : ic) {
            const std::size_t k = kk - 1;
            const std::size_t rstat = L.classes[k].refstat;
            r.X[k] = open ? T(Xchain[c] * Vsink[k]) : T(Xchain[c] * d.alpha(rstat - 1, k));
            const T vref = d.Vchain(rstat - 1, c);
            for (std::size_t i = 0; i < M; ++i) {
                const double S = L.stations[i].nservers;
                if (vref != zero) {
                    const T base = T(d.ST(i, k) * (Xchain[c] * d.Vchain(i, c) / vref) *
                                     d.alpha(i, k));
                    if (std::isinf(S)) {
                        r.U(i, k) = base;
                    } else if (Uchain.rows() == 0) {
                        r.U(i, k) = T(base / num_traits<T>::from_double(S));
                    } else {
                        r.U(i, k) = T(Uchain(i, c) * d.alpha(i, k));
                    }
                }
                if (d.Lchain(i, c) > zero) {
                    if (Qchain.rows() > 0) {
                        r.Q(i, k) = T(Qchain(i, c) * d.alpha(i, k));
                    } else if (d.STchain(i, c) != zero && vref != zero) {
                        r.Q(i, k) = T(Rchain(i, c) * d.ST(i, k) / d.STchain(i, c) * Xchain[c] *
                                      d.Vchain(i, c) / vref * d.alpha(i, k));
                    }
                    r.Tp(i, k) = T(Tchain(i, c) * d.alpha(i, k));
                    r.R(i, k) = r.Tp(i, k) == zero ? zero : T(r.Q(i, k) / r.Tp(i, k));
                }
            }
            const double njobs = L.classes[k].population;
            r.C[k] = r.X[k] == zero ? zero
                                    : T(num_traits<T>::from_double(njobs) / r.X[k]);
        }
    }
    auto absify = [&](Matrix<T>& A) {
        for (std::size_t i = 0; i < A.rows(); ++i)
            for (std::size_t j = 0; j < A.cols(); ++j)
                if (A(i, j) < zero) A(i, j) = T(-A(i, j));
    };
    absify(r.Q);
    absify(r.U);
    absify(r.R);
    absify(r.Tp);
    for (T& x : r.C)
        if (x < zero) x = T(-x);
    for (T& x : r.X)
        if (x < zero) x = T(-x);
    return r;
}

/** Product-form chain parameters, the split into queueing and delay stations. */
template <class T>
struct PfChainParams {
    Matrix<T> D;                 ///< (Mq x C) demands at queueing stations
    Matrix<T> Z;                 ///< (Md x C) demands at delay stations
    std::vector<T> lambda;       ///< (C) chain arrival rates, zero on a closed chain
    std::vector<double> N;       ///< (C) populations
    std::vector<double> S;       ///< (Mq) server counts
    std::vector<std::size_t> queue_stations;  ///< 1-based
    std::vector<std::size_t> delay_stations;  ///< 1-based
};

/**
 * Port of sn_get_product_form_chain_params.
 *
 * `lambda` is the chain-level arrival rate the reference builds by summing the class
 * rates of `sn_get_product_form_params` over `sn.inchain{c}`: the Source's rate for an
 * open class, nothing for a closed one. Only the mixed `lin` branch of solver_amva
 * reads it, so a closed model leaves the vector at zero.
 */
template <class T>
PfChainParams<T> sn_get_product_form_chain_params(const qn::NetworkStruct<T>& L, const ChainDemands<T>& d) {
    PfChainParams<T> p;
    for (std::size_t i = 0; i < L.nstations; ++i) {
        if (L.stations[i].nodetype == qn::NodeType::Queue) p.queue_stations.push_back(i + 1);
        else if (L.stations[i].nodetype == qn::NodeType::Delay) p.delay_stations.push_back(i + 1);
    }
    const T zero = num_traits<T>::from_int(0);
    p.D = Matrix<T>(p.queue_stations.size(), L.nchains, zero);
    p.Z = Matrix<T>(p.delay_stations.size(), L.nchains, zero);
    for (std::size_t c = 0; c < L.nchains; ++c) {
        for (std::size_t a = 0; a < p.queue_stations.size(); ++a)
            p.D(a, c) = d.Lchain(p.queue_stations[a] - 1, c);
        for (std::size_t a = 0; a < p.delay_stations.size(); ++a)
            p.Z(a, c) = d.Lchain(p.delay_stations[a] - 1, c);
    }
    p.N = d.Nchain;

    // An open class contributes the Source's rate for it, a closed one nothing. A
    // disabled pair contributes nothing either, which is where MATLAB's 'omitnan' sum
    // lands: sn.rates is NaN there.
    std::size_t source_station = 0;
    for (std::size_t i = 0; i < L.nstations; ++i)
        if (L.stations[i].nodetype == qn::NodeType::Source) source_station = i + 1;
    p.lambda.assign(L.nchains, zero);
    if (source_station > 0)
        for (std::size_t c = 0; c < L.nchains; ++c)
            for (std::size_t r = 0; r < L.nclasses; ++r) {
                if (!L.chains[c][r]) continue;
                if (!std::isinf(L.classes[r].population)) continue;
                if (L.disabled[source_station - 1][r]) continue;
                p.lambda[c] += L.rates(source_station - 1, r);
            }

    for (std::size_t a = 0; a < p.queue_stations.size(); ++a)
        p.S.push_back(L.stations[p.queue_stations[a] - 1].nservers);
    return p;
}

/**
 * Aggregate a class-indexed interlock matrix to the chain basis the MVA analyzers
 * work in. IL[r][s] is the share of the class-s queue that a class-r arrival must
 * not see, the interlocked flow of Franks (1999), Eq. (4.7). Two classes of the same
 * chain belong to the same client, so the diagonal blocks carry no information and
 * the chain diagonal stays zero: an arrival always sees its own chain in full.
 * Returns an empty matrix when nothing is interlocked.
 */
template <class T>
Matrix<T> sn_interlock_chain(const qn::NetworkStruct<T>& L,
                             const std::vector<std::vector<double>>& ILclass) {
    if (ILclass.empty()) return Matrix<T>();
    const std::size_t R = L.nclasses, K = L.nchains;
    if (ILclass.size() != R)
        throw InputError("sn_interlock_chain: the interlock matrix disagrees with the class count");
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> ILchain(K, K, zero);
    bool any = false;
    for (std::size_t cr = 0; cr < K; ++cr)
        for (std::size_t cs = 0; cs < K; ++cs) {
            if (cr == cs) continue;
            double best = 0.0;
            for (std::size_t r = 0; r < R; ++r) {
                if (!L.chains[cr][r]) continue;
                if (ILclass[r].size() != R)
                    throw InputError("sn_interlock_chain: the interlock matrix is not square");
                for (std::size_t sIl = 0; sIl < R; ++sIl) {
                    if (!L.chains[cs][sIl]) continue;
                    if (ILclass[r][sIl] > best) best = ILclass[r][sIl];
                }
            }
            if (best > 0.0) {
                ILchain(cr, cs) = num_traits<T>::from_double(best);
                any = true;
            }
        }
    return any ? ILchain : Matrix<T>();
}

}  // namespace mva
}  // namespace line

#endif  // LINE_SOLVERS_MVA_SN_CHAIN_H
