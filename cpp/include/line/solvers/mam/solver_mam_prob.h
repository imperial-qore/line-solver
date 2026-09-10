/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MAM_SOLVER_MAM_PROB_H
#define LINE_SOLVERS_MAM_SOLVER_MAM_PROB_H

/**
 * Port of `@@SolverMAM/getProb.m` and `@@SolverMAM/getProbMarg.m`: the joint
 * (level, phase) and marginal queue-length distributions at the model's single
 * queue, plus `@@SolverMAM/getMAMResult.m`, which exposes the matrix-analytic
 * internals themselves.
 *
 * THE SINGLE-QUEUE RESTRICTION IS THE REFERENCE'S OWN, and it is stated there
 * in prose: "the MAM solver uses QBD analysis, which is fundamentally a
 * single-queue method". A model with two or more Queue nodes is refused by
 * both, naming SolverCTMC and SolverSSA as the alternatives. The port carries
 * the same gate and the same advice.
 *
 * WHAT `getProb` ACTUALLY RETURNS, and why the doc comment matters more than
 * usual here. The reference computes the queue-length MARGINAL and multiplies
 * it by an arrival-weighted average of the service-phase initial vectors:
 *
 *     P(level, phase) ~= P(level) * avgPie(phase)
 *
 * That is a PRODUCT-FORM APPROXIMATION of the joint law, not the joint law of
 * the QBD, and the reference says so at its line 155 ("For now, approximate by
 * assuming phases are independent of level"). Reproduced exactly, including the
 * fallback to a uniform phase distribution when the weighted average degenerates.
 * The port does not silently improve it: a caller comparing against a CTMC
 * joint distribution needs to know which object this is.
 *
 * THE CLOSED BRANCH is likewise an approximation the reference is explicit
 * about: the arrival process is a single-phase MMAP built from the CONVERGED
 * PER-CLASS THROUGHPUTS, i.e. a Poisson surrogate at the fixed point, and
 * `getProbMarg` then reads class r's marginal off the aggregate distribution
 * truncated at N(r). Both are reproduced; neither is repaired.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mam/mmap_assemble.h"
#include "line/api/mam/mmapph1fcfs.h"
#include "line/api/qsys/qsys_bmapm1.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mam/mam_types.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

namespace prob_detail {

/** The reference's queue-station count, and the gate both accessors share. */
template <class T>
void require_single_queue(const qn::NetworkStruct<T>& L, const char* fn) {
    std::size_t nq = 0;
    for (std::size_t i = 0; i < L.nstations; ++i)
        if (L.stations[i].nodetype == qn::NodeType::Queue) ++nq;
    if (nq > 1)
        throw UnsupportedError(
            std::string(fn) +
            " is not supported for networks with multiple queues in SolverMAM. The MAM solver "
            "uses QBD (quasi-birth-death) analysis, which is fundamentally a single-queue "
            "method. Use SolverCTMC or SolverSSA for networks with multiple queues");
    if (nq == 0)
        throw UnsupportedError(std::string(fn) +
                               ": the model does not contain any queue stations");
}

/**
 * The per-class service phase-type pairs at station `ist`, scaled by 1/(rate c)
 * exactly as both accessors do.
 */
template <class T>
std::vector<PhService<T>> station_services(const qn::NetworkStruct<T>& L, std::size_t ist) {
    using lang::GlobalConstants;
    const T one = num_traits<T>::from_int(1);
    const std::size_t K = L.nclasses;
    const double ns = L.stations[ist - 1].nservers;
    std::vector<PhService<T>> svc(K);
    for (std::size_t r = 0; r < K; ++r) {
        if (L.disabled[ist - 1][r] || !(L.rates(ist - 1, r) > num_traits<T>::from_int(0))) {
            // The reference's isnan(D0) guard: an unserved class becomes Immediate.
            const T imm = num_traits<T>::from_double(GlobalConstants::Immediate);
            svc[r].sigma.assign(1, one);
            svc[r].S = Matrix<T>(1, 1, T(-imm));
            continue;
        }
        const T target =
            std::isfinite(ns)
                ? T(T(one / L.rates(ist - 1, r)) / num_traits<T>::from_double(ns))
                : T(one / L.rates(ist - 1, r));
        const Map<T> sc = map_scale(lang::dist_to_map(L.service[ist - 1][r]), target);
        svc[r].sigma = map_pie(sc);
        svc[r].S = sc.D0;
    }
    return svc;
}

/** The aggregate queue-length marginal, normalized, as both accessors take it. */
template <class T>
std::vector<T> level_marginal(const qn::NetworkStruct<T>& L, const MamOptions& opt,
                              std::size_t ist, const mva::AvgResult<T>& avg,
                              std::vector<PhService<T>>& svc, std::vector<T>& classLambda,
                              bool& closed) {
    using lang::GlobalConstants;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t K = L.nclasses;
    svc = station_services(L, ist);
    classLambda.assign(K, zero);

    closed = true;
    for (const qn::JobClass& c : L.classes)
        if (std::isinf(c.population)) closed = false;

    Mmap<T> arr;
    std::size_t maxLevel;
    if (closed) {
        // Poisson surrogate at the converged per-class throughputs.
        std::size_t Ntot = 0;
        for (const qn::JobClass& c : L.classes)
            Ntot += static_cast<std::size_t>(std::llround(c.population));
        maxLevel = Ntot + 1;
        T tot = zero;
        for (std::size_t r = 0; r < K; ++r) {
            classLambda[r] = avg.TN(ist - 1, r);
            tot += classLambda[r];
        }
        if (!(num_traits<T>::to_double(tot) >= GlobalConstants::FineTol))
            return std::vector<T>();  // no traffic: the caller emits the point mass at 0
        arr.D0 = Matrix<T>(1, 1, T(-tot));
        arr.D1 = Matrix<T>(1, 1, tot);
        for (std::size_t r = 0; r < K; ++r)
            arr.Dc.push_back(Matrix<T>(1, 1, classLambda[r]));
    } else {
        const std::size_t src = L.sourceIdx;
        if (src == 0) throw InputError("getProb: the open model has no Source");
        maxLevel = opt.cutoff > 0 ? opt.cutoff : static_cast<std::size_t>(100);
        std::vector<Map<T>> arrMaps(K);
        T tot = zero;
        for (std::size_t r = 0; r < K; ++r) {
            if (L.disabled[src - 1][r]) {
                arrMaps[r] = Map<T>{Matrix<T>(1, 1, zero), Matrix<T>(1, 1, zero)};
                continue;
            }
            arrMaps[r] = lang::dist_to_map(L.service[src - 1][r]);
            classLambda[r] = map_lambda(arrMaps[r]);
            tot += classLambda[r];
        }
        if (!(num_traits<T>::to_double(tot) >= GlobalConstants::FineTol))
            return std::vector<T>();
        if (K == 1) {
            arr.D0 = arrMaps[0].D0;
            arr.D1 = arrMaps[0].D1;
            arr.Dc.assign(1, arrMaps[0].D1);
        } else {
            // The reference superposes the per-class MAPs and then splits the
            // aggregate D1 by arrival-rate share. NOTE: MATLAB writes
            // `map_super({superMAP, arrMaps{k}})`, passing ONE cell argument to
            // a two-argument function, so this branch RAISES there; see the
            // divergence note in _kb. The port does what the code plainly
            // intends, which is the two-argument superposition.
            Map<T> super = arrMaps[0];
            for (std::size_t r = 1; r < K; ++r) {
                Map<T> m;
                m.D0 = krons(super.D0, arrMaps[r].D0);
                m.D1 = krons(super.D1, arrMaps[r].D1);
                super = map_normalize(m);
            }
            arr.D0 = super.D0;
            arr.D1 = super.D1;
            for (std::size_t r = 0; r < K; ++r) {
                Matrix<T> Dr = super.D1;
                const T w = T(classLambda[r] / tot);
                for (std::size_t i = 0; i < Dr.rows(); ++i)
                    for (std::size_t j = 0; j < Dr.cols(); ++j) Dr(i, j) *= w;
                arr.Dc.push_back(Dr);
            }
        }
    }

    // Both accessors capture ONE output of MMAPPH1FCFS('ncDistr'), i.e. class
    // 1's marginal, and renormalize it. Reproduced.
    const std::vector<std::vector<T>> d = mmapph1fcfs_ncdistr(arr, svc, maxLevel);
    std::vector<T> p(d[0].size(), zero);
    T mass = zero;
    for (std::size_t n = 0; n < d[0].size(); ++n) {
        p[n] = num_abs(d[0][n]);
        mass += p[n];
    }
    if (mass > zero)
        for (T& v : p) v /= mass;
    (void)one;
    return p;
}

}  // namespace prob_detail

/** The joint (level, phase) table `getProb` returns: rows levels, cols phases. */
template <class T>
struct ProbTable {
    Matrix<T> P;
};

/**
 * Port of `@@SolverMAM/getProb.m`.
 *
 * @param node 1-based NODE index; must be a station
 * @param avg  the converged `getAvg` result, which the closed branch reads
 *             throughputs from
 * @param L the refreshed struct
 * @param opt SolverMAM's options
 */
template <class T>
ProbTable<T> solver_mam_getprob(const qn::NetworkStruct<T>& L, const MamOptions& opt,
                                std::size_t node, const mva::AvgResult<T>& avg) {
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "getProb: the queue-length distribution comes from the MMAP[K]/PH[K]/1 age process, "
            "whose first-return matrix is a tolerance-terminated Riccati doubling; rerun with "
            "--arith double or --arith real");
    } else {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (node == 0 || node > L.nof_nodes())
        throw InputError("getProb: node number exceeds the number of nodes in the model");
    const std::size_t ist = L.nodes[node - 1].station;
    if (ist == 0) throw InputError("getProb: the specified node is not a station");
    prob_detail::require_single_queue(L, "getProb");

    const std::size_t K = L.nclasses;
    std::vector<PhService<T>> svc;
    std::vector<T> classLambda;
    bool closed = false;
    const std::vector<T> p =
        prob_detail::level_marginal(L, opt, ist, avg, svc, classLambda, closed);

    ProbTable<T> out;
    if (p.empty()) {
        // No traffic at this station: all the probability sits at (0, 1).
        std::size_t levels = 1;
        if (closed) {
            std::size_t Ntot = 0;
            for (const qn::JobClass& c : L.classes)
                Ntot += static_cast<std::size_t>(std::llround(c.population));
            levels = Ntot + 1;
        } else {
            levels = opt.cutoff > 0 ? opt.cutoff : static_cast<std::size_t>(100);
        }
        out.P = Matrix<T>(levels, 1, zero);
        out.P(0, 0) = one;
        return out;
    }

    std::size_t nPhases = 1;
    for (std::size_t r = 0; r < K; ++r) nPhases = std::max(nPhases, svc[r].S.rows());

    // The arrival-weighted average of the TIME-STATIONARY phase distributions,
    // map_prob, not the initial vectors sigma. sigma is the embedded
    // equilibrium at departure instants -- the phase a service STARTS in -- so
    // for an Erlang-2 it is [1 0] and the joint gave P(phase 2) = 0 at every
    // level for a server that spends half its busy time in phase 2. A PH law
    // (sigma, S) is the MAP {S, (-S e) sigma}, so the stationary vector needs no
    // extra state: PhService stays as MMAPPH1FCFS requires it.
    T tot = zero;
    for (std::size_t r = 0; r < K; ++r) tot += classLambda[r];
    std::vector<T> avgPie(nPhases, zero);
    if (tot > zero)
        for (std::size_t r = 0; r < K; ++r) {
            const std::size_t m = svc[r].S.rows();
            Map<T> ph;
            ph.D0 = svc[r].S;
            ph.D1 = Matrix<T>(m, m, zero);
            for (std::size_t i = 0; i < m; ++i) {
                T exit = zero;
                for (std::size_t j = 0; j < m; ++j) exit -= svc[r].S(i, j);
                for (std::size_t j = 0; j < m && j < svc[r].sigma.size(); ++j)
                    ph.D1(i, j) = T(exit * svc[r].sigma[j]);
            }
            const std::vector<T> piq = map_prob(ph);
            for (std::size_t i = 0; i < piq.size() && i < nPhases; ++i)
                avgPie[i] += T(piq[i] * classLambda[r] / tot);
        }
    T s = zero;
    for (const T& v : avgPie) s += v;
    if (!(num_traits<T>::to_double(s) > 0.0) || !std::isfinite(num_traits<T>::to_double(s))) {
        // The reference's fallback: a uniform phase distribution.
        for (T& v : avgPie) v = T(one / num_traits<T>::from_int(static_cast<int>(nPhases)));
    } else {
        for (T& v : avgPie) v /= s;
    }

    out.P = Matrix<T>(p.size(), nPhases, zero);
    for (std::size_t n = 0; n < p.size(); ++n)
        for (std::size_t i = 0; i < nPhases; ++i) out.P(n, i) = T(p[n] * avgPie[i]);
    return out;
    }  // if constexpr has_transcendental
}

/**
 * Port of `@@SolverMAM/getProbMarg.m`: P(n jobs of class `jobclass`) at station
 * `ist`, for n = 0..N(jobclass) closed, or over the whole cutoff range open.
 */
template <class T>
std::vector<T> solver_mam_getprobmarg(const qn::NetworkStruct<T>& L, const MamOptions& opt,
                                      std::size_t ist, std::size_t jobclass,
                                      const mva::AvgResult<T>& avg) {
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "getProbMarg: the queue-length distribution comes from the MMAP[K]/PH[K]/1 age "
            "process, whose first-return matrix is a tolerance-terminated Riccati doubling; "
            "rerun with --arith double or --arith real");
    } else {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (ist == 0 || ist > L.nstations)
        throw InputError("getProbMarg: station number exceeds the number of stations");
    if (jobclass == 0 || jobclass > L.nclasses)
        throw InputError("getProbMarg: job class index exceeds the number of classes");
    prob_detail::require_single_queue(L, "getProbMarg");

    std::vector<PhService<T>> svc;
    std::vector<T> classLambda;
    bool closed = false;
    const std::vector<T> p =
        prob_detail::level_marginal(L, opt, ist, avg, svc, classLambda, closed);

    if (p.empty()) {
        const std::size_t levels =
            closed ? static_cast<std::size_t>(std::llround(L.classes[jobclass - 1].population)) + 1
                   : static_cast<std::size_t>(100);
        std::vector<T> out(levels, zero);
        out[0] = one;
        return out;
    }
    if (!closed) return p;

    const double Nr = L.classes[jobclass - 1].population;
    if (!(Nr > 0.0))
        throw UnsupportedError(
            "getProbMarg: class " + std::to_string(jobclass) +
            " has zero population, so it has no marginal queue-length distribution");
    const std::size_t Nk = static_cast<std::size_t>(std::llround(Nr));
    std::vector<T> out(Nk + 1, zero);
    T mass = zero;
    for (std::size_t n = 0; n <= Nk && n < p.size(); ++n) {
        out[n] = p[n];
        mass += out[n];
    }
    if (mass > zero)
        for (T& v : out) v /= mass;
    return out;
    }  // if constexpr has_transcendental
}

/**
 * Port of `@@SolverMAM/getMAMResult.m`: the matrix-analytic internals of a
 * single-queue model, for a BMAP (or MAP) arrival stream into an exponential
 * single server.
 *
 * The reference's other arm handles a retrial station through
 * `qsys_bmapphnn_retrial`. That routine IS ported, but the C++ `NetworkStruct`
 * has no retrial fields at all, so no model this port can build reaches it; the
 * arm is recorded rather than written.
 */
template <class T>
qsys::BmapM1Result<T> solver_mam_getmamresult(const qn::NetworkStruct<T>& L) {
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "getMAMResult: qsys_bmapm1 runs a functional iteration for G and an adaptive level "
            "truncation, both tolerance-terminated; rerun with --arith double or --arith real");
    } else {
    std::size_t src = 0, q = 0;
    for (std::size_t i = 1; i <= L.nstations; ++i) {
        const qn::NodeType nt = L.nodes[L.node_of_station(i) - 1].nodetype;
        if (nt == qn::NodeType::Source) src = i;
        else if (nt == qn::NodeType::Queue) {
            if (q != 0)
                throw UnsupportedError(
                    "getMAMResult exposes the matrix-analytic internals of a single-queue model "
                    "only");
            q = i;
        }
    }
    if (src == 0 || q == 0)
        throw UnsupportedError(
            "getMAMResult requires an open model with one Source and one Queue");
    if (L.nclasses > 1)
        throw UnsupportedError(
            "getMAMResult exposes the matrix-analytic internals of a single-class model only");
    if (L.stations[q - 1].nservers != 1.0)
        throw UnsupportedError("getMAMResult requires a single-server queue");

    const Map<T> arv = lang::dist_to_map(L.service[src - 1][0]);
    const Map<T> svc = lang::dist_to_map(L.service[q - 1][0]);
    if (svc.D0.rows() != 1)
        throw UnsupportedError(
            "getMAMResult exposes the M/G/1-type internals for exponential service only; the "
            "queue has a multi-phase service process");
    const T mu = T(-svc.D0(0, 0));

    // A MAP is the batch-size-one BMAP; a genuine BMAP would carry Dmark.
    std::vector<Matrix<T>> D;
    D.push_back(arv.D0);
    const std::vector<Matrix<T>>& batches = L.service[src - 1][0].Dmark;
    if (batches.empty()) {
        D.push_back(arv.D1);
    } else {
        for (const Matrix<T>& Dk : batches) D.push_back(Dk);
    }
    return qsys::qsys_bmapm1(D, mu);
    }  // if constexpr has_transcendental
}

}  // namespace mam
}  // namespace line

#endif  // LINE_SOLVERS_MAM_SOLVER_MAM_PROB_H
