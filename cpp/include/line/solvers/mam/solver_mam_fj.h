/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MAM_SOLVER_MAM_FJ_H
#define LINE_SOLVERS_MAM_SOLVER_MAM_FJ_H

/**
 * Port of `solver_mam_fj.m`, the fork-join route of SolverMAM.
 *
 * The reference analyzer is a THIN WRAPPER around one third-party engine:
 * `mainFJ` of `matlab/lib/thirdparty/FJ_codes`, the response-time-tail
 * approximation of Z. Qiu, J. F. Perez and P. Harrison, "Beyond the Mean in
 * Fork-Join Queues: Efficient Approximation for Response-Time Tails" (IFIP
 * Performance 2015). Everything `solver_mam_fj.m` computes on its own is either
 * a topology check or an M/M/1 closed form; the one quantity that makes it a
 * FORK-JOIN analyzer -- the synchronisation delay charged at the Join -- is
 * `mainFJ`'s mean response time minus the branch response time.
 *
 * WHAT IS PORTED HERE: the gate, the parameter extraction and the wrapper,
 * i.e. `fj_is_homogeneous.m`, `fj_extract_params.m` (on top of the already-
 * ported `api/fj/fj_dist2fj.h`) and `solver_mam_fj.m` itself.
 * `mam_fj_is_homogeneous` is the 2a predicate of `solver_mam_analyzer.m` and is
 * what separates the fork-join models this analyzer claims from the ones that
 * go to `solver_mam_basic_mmap`.
 *
 * THE ENGINE IS `api/fj/fj_codes.h` and `api/fj/fj_codes_matrices.h`: `mainFJ`
 * and the chain it drives -- `returnRT1`, `returnRT2`, `computeT`,
 * `computeT_NARE`, `constructSRK`, `build_SA`, `computePi`, `returnWait`,
 * `returnPer`, `generateService`. It is double only, because the T matrix is
 * the stable invariant subspace of a Riccati pencil (ordered real Schur) and
 * the two Sylvester equations are of order (C + 1) m^2 ma; this file refuses by
 * name at any other arithmetic.
 *
 * WHY THE MMT TRANSFORM IS NOT SUBSTITUTED FOR IT. `mva/fj_mmt.h` plus
 * `mva/fj_driver.h` do solve fork-join models, and driving them with
 * `solver_mam_basic` as the inner solve would produce a full metric tuple. It
 * would not be this analyzer's tuple, and on this analyzer's models it would
 * not be sound either: `fj_mmt` REUSES the model's own Source and Sink when it
 * has them and adds a pair only for a closed layer, so on an open model the
 * auxiliary stream is injected at the SAME Source that carries the real open
 * chain -- the transform then reports one merged arrival stream where this
 * analyzer's decomposition keeps the branch streams apart. The transform is
 * exercised by SolverMVA and
 * SolverNC on CLOSED layers, which is the case it was written for. Reporting
 * its output as `solver_mam_fj` would be a different algorithm on a corrupted
 * struct, under the reference's method name.
 *
 * WHAT THE `api/fj/` FAMILY IS NOT. This codebase does carry ported fork-join
 * approximations -- `fj_synch_delay`, `fj_respt_2way`, `fj_respt_nt`,
 * `fj_respt_varki`, `fj_respt_vm`, `fj_rmax`, `fj_order_stat`. They are
 * DIFFERENT published approximations. Substituting one of them for `mainFJ`
 * would answer the model under a method name that promises the Qiu-Perez-
 * Harrison tail estimate, which is precisely the substitution
 * `solver_mam_runner.h` refuses to make for non-Markovian service.
 *
 * REFERENCE DEFECTS in solver_mam_fj.m. 1, 2 and 4 are reproduced -- they are
 * what the reference reports and changing them would answer differently under
 * its method name; 3 is not, see `solver_mam_fj` for why that is not a change
 * of result:
 *
 *  1. THE BRANCH METRICS IGNORE THE SERVICE REPRESENTATION ENTIRELY. Every
 *     parallel queue is reported at U = lambda/mu, Q = rho/(1 - rho) and
 *     R = 1/(mu - lambda), the M/M/1 closed forms, even though the gate has
 *     just admitted Erlang(2), HyperExp(2) and MAP(2) service and the arrival
 *     may be a MAP(2). For anything but Exp/Exp those are the numbers of a
 *     different queue, and they are not flagged as approximate.
 *  2. NO STABILITY GUARD ON THE REPORTED MEANS. When lambda >= mu the same
 *     formulas return a negative queue length and a negative response time.
 *     `mainFJ` itself errors on load >= 1 and `fj_extract_params` warns, so the
 *     model is stopped -- but only by the engine, one call later, and the two
 *     disagree about whether it is an error.
 *  3. `mainFJ` IS CALLED TWICE PER CLASS, once on a 21-point percentile grid
 *     for the mean and once on the four stored percentiles, and the second call
 *     recomputes the entire T matrix to read four points off a curve the first
 *     call already produced.
 *  4. THE MEAN IS A TRAPEZOID OVER THE INVERSE CDF with an invented tail: the
 *     grid stops at the 99.9th percentile and the integration closes at p = 1
 *     with `RTp(end) * 1.1`. The 1.1 has no derivation, and for a heavy tail it
 *     is the dominant term of the mean.
 *
 * ARITHMETIC. `fj_dist2fj` needs one linear solve per process and stays in the
 * field, so the GATE and the EXTRACTION hold at every instantiation and an
 * unsupported topology or distribution is still named exactly at Rational. The
 * ENGINE is double only and `solver_mam_fj` refuses past that point, which is
 * why the two checks run before the arithmetic gate rather than after it.
 */

#include <cmath>
#include <cstddef>
#include <string>
#include <type_traits>
#include <vector>

#include "line/api/fj/fj_codes.h"
#include "line/api/fj/fj_dist2fj.h"
#include "line/api/sn/sn_join_droprate.h"
#include "line/api/mam/map_moment.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mam/mam_types.h"
#include "line/util/error.h"

namespace line {
namespace mam {

/** What `fj_is_homogeneous.m` returns: the fork-join pair, or why there is none. */
struct MamFjInfo {
    bool ok = false;
    std::size_t forkNode = 0;               ///< 1-based node index
    std::size_t joinNode = 0;               ///< 1-based node index
    std::vector<std::size_t> queueNodes;    ///< the K parallel branches, 1-based node indices
    std::size_t K = 0;
    std::string why;                        ///< the reference's `errorMsg`, empty when ok
};

/** What `fj_extract_params.m` returns: one arrival and one service per class. */
template <class T>
struct MamFjParams {
    std::vector<fj::FjDist<T> > arrival;
    std::vector<fj::FjDist<T> > service;
    std::size_t K = 0;
};

namespace fj_detail {

/**
 * The class-aggregated node routing the reference reads as `sn.rtnodes(i, j)`:
 * is there ANY class pair routed from i to j? A branch is identified by the
 * topology, not by which class traverses it.
 */
template <class T>
bool routes_between(const qn::NetworkStruct<T>& L, std::size_t i, std::size_t j) {
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t r = 1; r <= L.nclasses; ++r)
        for (std::size_t s = 1; s <= L.nclasses; ++s)
            if (L.get_route(r, s, i, j) > zero) return true;
    return false;
}

/** The reference's entrywise (D0, D1) comparison at GlobalConstants::FineTol. */
template <class T>
bool same_map(const Map<T>& a, const Map<T>& b) {
    if (a.order() != b.order()) return false;
    const double tol = lang::GlobalConstants::FineTol;
    for (std::size_t i = 0; i < a.order(); ++i)
        for (std::size_t j = 0; j < a.order(); ++j) {
            if (std::fabs(num_traits<T>::to_double(a.D0(i, j)) -
                          num_traits<T>::to_double(b.D0(i, j))) > tol)
                return false;
            if (std::fabs(num_traits<T>::to_double(a.D1(i, j)) -
                          num_traits<T>::to_double(b.D1(i, j))) > tol)
                return false;
        }
    return true;
}

/**
 * `sn.procid` as the distribution code FJ_codes switches on. The two enums
 * share MATLAB's numbering, so this is a check that the type is in the accepted
 * set rather than a translation.
 */
inline bool to_fj_proc_type(lang::ProcessType p, fj::FjProcType& out) {
    switch (p) {
        case lang::ProcessType::EXP: out = fj::FjProcType::Exp; return true;
        case lang::ProcessType::ERLANG: out = fj::FjProcType::Erlang; return true;
        case lang::ProcessType::HYPEREXP: out = fj::FjProcType::HyperExp; return true;
        case lang::ProcessType::MAP: out = fj::FjProcType::Map; return true;
        default: return false;
    }
}

}  // namespace fj_detail

/**
 * Port of `fj_is_homogeneous.m`.
 *
 * This is NOT a test for the presence of a Fork: most fork-join models fail it.
 * It tests membership in the homogeneous class the FJ_codes approximation is
 * defined on -- one fork-join pair, K parallel queues between them carrying
 * identical service, open classes, FCFS or PS -- and it is the 2a/2b branch
 * predicate of the analyzer dispatch. It returns rather than throws for that
 * reason, with the rejection reason in `why`.
 */
template <class T>
MamFjInfo mam_fj_is_homogeneous(const qn::NetworkStruct<T>& L) {
    using lang::SchedStrategy;
    MamFjInfo info;

    if (!L.is_open_model()) {
        info.why = "FJ_codes supports open queueing models only";
        return info;
    }
    if (L.fj.empty()) {
        info.why = "the network contains no Fork-Join structure";
        return info;
    }

    std::vector<std::size_t> forks, joins;
    for (std::size_t i = 1; i <= L.nodes.size(); ++i) {
        if (L.nodes[i - 1].nodetype == qn::NodeType::Fork) forks.push_back(i);
        if (L.nodes[i - 1].nodetype == qn::NodeType::Join) joins.push_back(i);
    }
    if (forks.empty() || joins.empty()) {
        info.why = "the network must contain both a Fork and a Join node";
        return info;
    }
    if (forks.size() > 1 || joins.size() > 1) {
        info.why = "FJ_codes supports a single Fork-Join pair; this model has " +
                   std::to_string(forks.size()) + " forks and " + std::to_string(joins.size()) +
                   " joins";
        return info;
    }
    // FJ_codes computes the response-time tail of an AND-join: `mainFJ`
    // synchronises on the LAST branch and has no parameter for a quorum, so a
    // k-of-n model routed here would come back with the all-join tail under a
    // quorum's name -- the same number for every k. Refuse it by name;
    // `fj_tail_ordstat` covers the quorum.
    {
        const std::size_t nsib = L.join_siblings(joins[0]);
        if (nsib > 0 && sn::sn_join_quorum(L, joins[0], nsib) < nsib) {
            info.why = "the Join fires on a quorum; FJ_codes synchronises on every branch and "
                       "has no quorum, so use the order-statistic tail fj_tail_ordstat instead";
            return info;
        }
    }
    info.forkNode = forks[0];
    info.joinNode = joins[0];

    bool paired = false;
    for (std::size_t p = 0; p < L.fj.size(); ++p)
        if (L.fj[p].first == info.forkNode && L.fj[p].second == info.joinNode) paired = true;
    if (!paired) {
        info.why = "the Fork and the Join of this model are not paired with each other";
        return info;
    }

    // A branch is a Queue the fork routes INTO and that routes into the join.
    for (std::size_t i = 1; i <= L.nodes.size(); ++i) {
        if (L.nodes[i - 1].nodetype != qn::NodeType::Queue) continue;
        if (fj_detail::routes_between(L, info.forkNode, i) &&
            fj_detail::routes_between(L, i, info.joinNode))
            info.queueNodes.push_back(i);
    }
    if (info.queueNodes.empty()) {
        info.why = "no Queue node lies between the Fork and the Join";
        return info;
    }
    info.K = info.queueNodes.size();

    // Homogeneous branches: the approximation solves ONE branch queue and
    // raises its response time to the K-branch maximum, so branches that differ
    // have no single queue to solve.
    const std::size_t st0 = L.nodes[info.queueNodes[0] - 1].station;
    for (std::size_t r = 1; r <= L.nclasses; ++r) {
        if (st0 == 0 || L.disabled[st0 - 1][r - 1]) {
            info.why = "branch queue '" + L.nodes[info.queueNodes[0] - 1].name +
                       "' has no service distribution for class '" + L.classes[r - 1].name + "'";
            return info;
        }
        const Map<T> first = lang::dist_to_map(L.service[st0 - 1][r - 1]);
        for (std::size_t k = 1; k < info.K; ++k) {
            const std::size_t st = L.nodes[info.queueNodes[k] - 1].station;
            if (st == 0 || L.disabled[st - 1][r - 1]) {
                info.why = "branch queue '" + L.nodes[info.queueNodes[k] - 1].name +
                           "' has no service distribution for class '" + L.classes[r - 1].name +
                           "'";
                return info;
            }
            if (!fj_detail::same_map(first, lang::dist_to_map(L.service[st - 1][r - 1]))) {
                info.why = "the parallel queues have heterogeneous service distributions for "
                           "class '" + L.classes[r - 1].name + "'; FJ_codes requires homogeneous "
                           "servers";
                return info;
            }
        }
    }

    for (std::size_t k = 0; k < info.K; ++k) {
        const std::size_t st = L.nodes[info.queueNodes[k] - 1].station;
        const SchedStrategy sc = L.stations[st - 1].sched;
        if (sc != SchedStrategy::FCFS && sc != SchedStrategy::PS) {
            info.why = "branch queue '" + L.nodes[info.queueNodes[k] - 1].name +
                       "' uses a scheduling strategy FJ_codes does not support; it supports FCFS "
                       "and PS only";
            return info;
        }
    }

    info.ok = true;
    return info;
}

/**
 * Port of `fj_extract_params.m`: the arrival descriptor from the Source and the
 * service descriptor from the first branch, per class.
 *
 * The reference WARNS on an unstable class and calls `mainFJ` anyway, which
 * then errors with "System not stable"; the two disagree about the severity of
 * the same condition. This refuses at the earlier point, which is where the
 * model can still be named.
 */
template <class T>
MamFjParams<T> mam_fj_extract_params(const qn::NetworkStruct<T>& L, const MamFjInfo& info) {
    if (!info.ok) throw InputError("mam_fj_extract_params: the topology gate has not passed");

    std::size_t src = 0;
    for (std::size_t i = 1; i <= L.nstations; ++i)
        if (L.stations[i - 1].nodetype == qn::NodeType::Source) {
            src = i;
            break;
        }
    if (src == 0) throw InputError("SolverMAM: the fork-join model has no Source node");
    const std::size_t q0 = L.nodes[info.queueNodes[0] - 1].station;

    MamFjParams<T> par;
    par.K = info.K;
    par.arrival.reserve(L.nclasses);
    par.service.reserve(L.nclasses);
    for (std::size_t r = 1; r <= L.nclasses; ++r) {
        fj::FjProcType at, st;
        if (!fj_detail::to_fj_proc_type(L.service[src - 1][r - 1].type, at))
            throw UnsupportedError(
                "SolverMAM: the arrival process of class '" + L.classes[r - 1].name +
                "' is " + lang::process_to_text(L.service[src - 1][r - 1].type) +
                ", and the FJ_codes approximation is defined for Exp, HyperExp(2), Erlang(2) and "
                "MAP(2) arrivals only");
        if (!fj_detail::to_fj_proc_type(L.service[q0 - 1][r - 1].type, st))
            throw UnsupportedError(
                "SolverMAM: the branch service process of class '" + L.classes[r - 1].name +
                "' is " + lang::process_to_text(L.service[q0 - 1][r - 1].type) +
                ", and the FJ_codes approximation is defined for Exp, HyperExp(2) and Erlang(2) "
                "service only");
        const fj::FjDist<T> a =
            fj::fj_dist2fj(lang::dist_to_map(L.service[src - 1][r - 1]), fj::FjDistKind::Arrival,
                           at);
        const fj::FjDist<T> s = fj::fj_dist2fj(lang::dist_to_map(L.service[q0 - 1][r - 1]),
                                               fj::FjDistKind::Service, st);
        // mainFJ's own first act: load = lambda / mu, and load >= 1 is an error
        // there. Every branch sees the full arrival stream, so this is the
        // per-branch load and not an aggregate.
        if (!(num_traits<T>::to_double(a.lambda) < num_traits<T>::to_double(s.mu)))
            throw InputError(
                "SolverMAM: class '" + L.classes[r - 1].name +
                "' offers each fork-join branch a load of at least one (arrival rate " +
                std::to_string(num_traits<T>::to_double(a.lambda)) + " against service rate " +
                std::to_string(num_traits<T>::to_double(s.mu)) +
                "); the response-time tail of an unstable branch does not exist");
        par.arrival.push_back(a);
        par.service.push_back(s);
    }
    return par;
}

/** The four percentiles `solver_mam_fj.m` stores for `getPerctRespT`. */
inline const std::vector<double>& mam_fj_stored_percentiles() {
    static const std::vector<double> p = {0.50, 0.90, 0.95, 0.99};
    return p;
}

/**
 * The 21-point grid `solver_mam_fj.m` inverts for the MEAN, `[0.01:0.05:0.95,
 * 0.99, 0.999]`. It is a percentile grid and not a quadrature rule; see defect
 * 4 in the header for what closing it at p = 1 costs.
 */
inline const std::vector<double>& mam_fj_dense_percentiles() {
    static const std::vector<double> p = [] {
        std::vector<double> v;
        for (int k = 0; k < 19; ++k) v.push_back(0.01 + 0.05 * static_cast<double>(k));
        v.push_back(0.99);
        v.push_back(0.999);
        return v;
    }();
    return p;
}

namespace fj_detail {

/** `options.config.fj_accuracy`, the reference's `Cs`, default 100. */
inline std::size_t fj_accuracy(const MamOptions& opt) {
    if (opt.fj_accuracy < 1)
        throw InputError("SolverMAM: config.fj_accuracy is the FJ_codes truncation C of the "
                         "queue-length difference and must be at least 1");
    return opt.fj_accuracy;
}

/**
 * The trapezoid over the inverse CDF that `solver_mam_fj.m` calls the mean.
 *
 * Reproduced including the invented tail point at p = 1: the grid stops at the
 * 99.9th percentile and the integration closes with `RTp(end) * 1.1`, a factor
 * with no derivation behind it. For a heavy tail that last panel is the
 * dominant term, so the number is the reference's and not a mean.
 */
inline double fj_mean_from_percentiles(const std::vector<double>& pers,
                                       const std::vector<double>& rtp) {
    std::vector<double> p, v;
    p.push_back(0.0);
    v.push_back(0.0);
    for (std::size_t i = 0; i < pers.size(); ++i) {
        p.push_back(pers[i]);
        v.push_back(rtp[i]);
    }
    p.push_back(1.0);
    v.push_back(rtp.empty() ? 0.0 : rtp.back() * 1.1);
    double acc = 0.0;
    for (std::size_t i = 0; i + 1 < p.size(); ++i)
        acc += 0.5 * (p[i + 1] - p[i]) * (v[i + 1] + v[i]);
    return acc;
}

}  // namespace fj_detail

/**
 * Port of `solver_mam_fj.m`.
 *
 * The analyzer is a thin wrapper: `mainFJ` returns the response-time
 * percentiles of the whole fork-join subnetwork, the mean is read off that
 * curve by trapezoid, and every station metric the reference reports is an
 * M/M/1 closed form in (lambda, mu) EXCEPT the Join, which carries the
 * synchronisation delay, i.e. the fork-join mean minus one branch response
 * time. Defects 1, 2 and 4 of the header are reproduced deliberately: they are
 * what the reference reports, and the alternative would be a different set of
 * numbers under its method name.
 *
 * DEFECT 3 IS NOT REPRODUCED. The reference calls `mainFJ` TWICE per class, on
 * the dense grid for the mean and again on the four stored percentiles, and the
 * second call rebuilds the entire T matrix to read four points off a curve the
 * first call already produced. Here the two grids are CONCATENATED into one
 * call, which is the same inversion of the same phase-type law on a union of
 * grids: `returnPer` treats every requested percentile independently, so a
 * point's value does not depend on which other points were asked for.
 */
template <class T>
mva::MvaSolution<T> solver_mam_fj(const qn::NetworkStruct<T>& L, const MamOptions& opt,
                                  std::vector<std::vector<T> >* percentiles_out) {
    const MamFjInfo info = mam_fj_is_homogeneous(L);
    if (!info.ok)
        throw UnsupportedError(
            "SolverMAM: model '" + L.name +
            "' is not a homogeneous fork-join network (" + info.why +
            "). The reference routes such a model to solver_mam_basic_mmap, the MMAP fork-join "
            "decomposition; select it with method 'dec.source.mmap'");
    if constexpr (!std::is_same<T, double>::value) {
        (void)opt;
        (void)percentiles_out;
        throw UnsupportedError(
            "SolverMAM: the FJ_codes fork-join engine solves a non-symmetric algebraic Riccati "
            "equation by an ORDERED REAL SCHUR factorization and two Sylvester equations of order "
            "(C + 1) m^2 ma by Bartels-Stewart, all of which are LAPACK and therefore double only. "
            "Re-run this model at double; there is no exact or extended-precision route to the "
            "Qiu-Perez-Harrison approximation");
    } else {
        const MamFjParams<T> par = mam_fj_extract_params(L, info);
        const std::size_t C = fj_detail::fj_accuracy(opt);
        const fj::FjTMode mode = fj::fj_parse_tmode(opt.fj_tmode);

        const std::vector<double>& dense = mam_fj_dense_percentiles();
        const std::vector<double>& stored = mam_fj_stored_percentiles();
        std::vector<double> grid = dense;
        grid.insert(grid.end(), stored.begin(), stored.end());

        const std::size_t M = L.nstations, K = L.nclasses;
        mva::MvaSolution<T> sol;
        sol.Q = Matrix<T>(M, K, num_traits<T>::from_int(0));
        sol.U = Matrix<T>(M, K, num_traits<T>::from_int(0));
        sol.R = Matrix<T>(M, K, num_traits<T>::from_int(0));
        sol.Tp = Matrix<T>(M, K, num_traits<T>::from_int(0));
        // CN and XN are zero in the reference too: solver_mam_fj.m allocates
        // them and fills neither.
        sol.C.assign(K, num_traits<T>::from_int(0));
        sol.X.assign(K, num_traits<T>::from_int(0));
        sol.iter = 0;  // FJ_codes is direct, as the reference's totiter = 0 records
        if (percentiles_out != nullptr)
            percentiles_out->assign(K, std::vector<T>(stored.size(),
                                                      num_traits<T>::from_int(0)));

        std::size_t src = 0, snk = 0;
        for (std::size_t i = 1; i <= M; ++i)
            if (L.stations[i - 1].nodetype == qn::NodeType::Source) src = i;
        for (std::size_t i = 1; i <= L.nodes.size(); ++i)
            if (L.nodes[i - 1].nodetype == qn::NodeType::Sink) snk = i;

        for (std::size_t r = 1; r <= K; ++r) {
            const std::vector<std::size_t> Kv(1, info.K);
            const std::vector<std::size_t> Cv(1, C);
            const std::vector<fj::FjCodesPercentiles> res =
                fj::fj_main(par.arrival[r - 1], par.service[r - 1], grid, Kv, Cv, mode);

            std::vector<double> dense_rt(dense.size()), stored_rt(stored.size());
            for (std::size_t i = 0; i < dense.size(); ++i) dense_rt[i] = res[0].RTp[i];
            for (std::size_t i = 0; i < stored.size(); ++i)
                stored_rt[i] = res[0].RTp[dense.size() + i];
            if (percentiles_out != nullptr)
                for (std::size_t i = 0; i < stored.size(); ++i)
                    (*percentiles_out)[r - 1][i] = stored_rt[i];

            const double mean_fj_rt = fj_detail::fj_mean_from_percentiles(dense, dense_rt);
            const double lambda = par.arrival[r - 1].lambda;
            const double mu = par.service[r - 1].mu;
            const double rho = lambda / mu;
            const double branch_rt = 1.0 / (mu - lambda);

            for (std::size_t k = 0; k < info.K; ++k) {
                const std::size_t st = L.nodes[info.queueNodes[k] - 1].station;
                sol.U(st - 1, r - 1) = rho;
                sol.Tp(st - 1, r - 1) = lambda;
                sol.Q(st - 1, r - 1) = rho / (1.0 - rho);
                sol.R(st - 1, r - 1) = branch_rt;
            }
            const std::size_t fst = L.nodes[info.forkNode - 1].station;
            if (fst != 0) sol.Tp(fst - 1, r - 1) = lambda;
            const std::size_t jst = L.nodes[info.joinNode - 1].station;
            if (jst != 0) {
                sol.Tp(jst - 1, r - 1) = lambda;
                const double sync = mean_fj_rt - branch_rt;
                sol.R(jst - 1, r - 1) = sync;
                sol.Q(jst - 1, r - 1) = lambda * sync;
            }
            if (src != 0) sol.Tp(src - 1, r - 1) = lambda;
            if (snk != 0) {
                const std::size_t sst = L.nodes[snk - 1].station;
                if (sst != 0) sol.Tp(sst - 1, r - 1) = lambda;
            }
        }
        return sol;
    }
}

/** `solver_mam_fj.m` without the percentile side channel. */
template <class T>
mva::MvaSolution<T> solver_mam_fj(const qn::NetworkStruct<T>& L, const MamOptions& opt) {
    return solver_mam_fj<T>(L, opt, nullptr);
}

/**
 * `@@SolverMAM/getPerctRespT.m`'s fork-join path: the percentiles
 * `solver_mam_fj.m` stores in `percResults.RT`, one row per class, at the four
 * levels `mam_fj_stored_percentiles()` names.
 *
 * A REQUESTED LEVEL IS INTERPOLATED, NOT RE-SOLVED, and outside [0.50, 0.99] it
 * is EXTRAPOLATED off the end segment: the reference reads its stored table
 * with `interp1(..., 'linear', 'extrap')`. That is reproduced, defect and all.
 * Asking for the 0.999 quantile of a heavy tail therefore continues the
 * 0.95-to-0.99 chord rather than inverting the law again, and the further out
 * the level the worse the estimate -- `returnPer` would answer it exactly, for
 * the price of another solve. The reference's choice is kept because a caller
 * comparing MATLAB against this port must see the same number.
 */
template <class T>
std::vector<std::vector<T> > solver_mam_fj_percentiles(const qn::NetworkStruct<T>& L,
                                                       const MamOptions& opt,
                                                       const std::vector<double>& percentiles) {
    const std::vector<double>& stored = mam_fj_stored_percentiles();
    std::vector<std::vector<T> > all;
    solver_mam_fj<T>(L, opt, &all);

    std::vector<std::vector<T> > out(all.size(), std::vector<T>(percentiles.size(),
                                                                num_traits<T>::from_int(0)));
    for (std::size_t r = 0; r < all.size(); ++r)
        for (std::size_t i = 0; i < percentiles.size(); ++i) {
            const double p = percentiles[i];
            // interp1 linear with extrapolation: the bracketing segment, or the
            // first / last one continued.
            std::size_t s = 0;
            while (s + 2 < stored.size() && p > stored[s + 1]) ++s;
            const double p0 = stored[s], p1 = stored[s + 1];
            const T v0 = all[r][s], v1 = all[r][s + 1];
            out[r][i] = T(v0 + (v1 - v0) * num_traits<T>::from_double((p - p0) / (p1 - p0)));
        }
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_SOLVERS_MAM_SOLVER_MAM_FJ_H
