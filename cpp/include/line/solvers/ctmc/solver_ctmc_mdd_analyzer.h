/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_CTMC_SOLVER_CTMC_MDD_ANALYZER_H
#define LINE_SOLVERS_CTMC_SOLVER_CTMC_MDD_ANALYZER_H

/**
 * The `mdd` method of SolverCTMC: stationary analysis of a closed single-class
 * network whose state space is held in a decision diagram and solved by level
 * aggregation.
 *
 * Port of matlab/src/solvers/CTMC/solver_ctmc_mdd_analyzer.m, after
 * A.S. Miner, G. Ciardo, S. Donatelli, "Using the exact state space of a Markov
 * model to compute approximate stationary measures", SIGMETRICS 2000.
 *
 * IT NEVER FORMS THE |S|-STATE GENERATOR, which is the whole point and also why
 * it sits beside `solver_ctmc_analyzer` rather than inside it: the reachable set
 * lives in an MDD and K coupled level-CTMCs are iterated to a fixed point, so
 * the memory cost is O(sum_k |M_k|) rather than O(|S|). The saving grows with
 * the number of stations and is NEGATIVE at K = 3, where the diagram compresses
 * nothing.
 *
 * EXACTNESS. The single approximation is Pr{i_k | alpha} = Pr{i_k | p}. It is
 * EXACT on product-form networks (paper Sec. 5), which covers exponential
 * service under any work-conserving discipline and general service at PS or IS
 * stations (BCMP types 2 and 3). It is an approximation otherwise, notably
 * phase-type service at FCFS or LCFS, where errors of a fraction of a percent on
 * the mean queue lengths have been observed. `no_aggregation` additionally
 * certifies exactness STRUCTURALLY: when no diagram node is shared, conditioning
 * on the node equals conditioning on the whole path and the approximation is an
 * identity. False means "not certified", never "approximate" -- a product-form
 * model is exact however much its diagram shares.
 *
 * WHICH ENCODING. Exponential service is discipline-insensitive for the
 * queue-length law, so the compact count encoding of `mdd_descriptor` serves any
 * work-conserving station. Phase-type service is not: `mdd_descriptor`'s
 * count-plus-one-phase local state is NON-preemptive, while a shared server has
 * every job present in service holding its own phase and needs the per-phase
 * counts of `mdd_ps`. A model that mixes the two cases is refused rather than
 * modelled under whichever encoding happens to be picked.
 */

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mdd/mdd.h"
#include "line/api/mdd/mdd_descriptor.h"
#include "line/api/mdd/mdd_mcd.h"
#include "line/api/mdd/mdd_ps.h"
#include "line/api/spn/spn_mdd.h"
#include "line/api/mdd/mdd_reachset.h"
#include "line/api/mdd/mdd_types.h"
#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/solvers/ctmc/solver_ctmc.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace ctmc {

// NetworkStruct, NodeType and SchedStrategy are already in scope from
// solver_ctmc.h; re-declaring them here would be a redundant using.

/** What one `mdd` solve produces beside the means, i.e. the reference's INFO. */
template <class T>
struct CtmcMddSolution {
    CtmcAvg<T> avg;
    /** |S|, counted in the diagram without ever listing a state. */
    long long num_states = 0;
    /** |M_k| per paper level; their sum is what the diagram actually holds. */
    std::vector<std::size_t> level_sizes;
    /** Coupled fixed-point sweeps performed. */
    int iters = 0;
    /** max |A(p)| per paper level; 1 everywhere means no node is shared. */
    std::vector<double> paths_per_level;
    /** True certifies the answer is exact structurally; see the file header. */
    bool no_aggregation = false;
    /** Which local encoding was picked, "np" or "ps". */
    std::string encoding;
    std::string actualmethod = "mdd";
};

namespace mdd_analyzer_detail {

/**
 * Station-to-station routing of the single class, read off `sn.rt`.
 *
 * `rt` is indexed by STATEFUL node and class-major, so the row of station i is
 * `(stateful_of_station(i) - 1) * R + 0` with R = 1 here. A row that does not
 * sum to one means a completion can leave the network or stay put, neither of
 * which the closed descriptor represents.
 */
template <class T>
Matrix<T> mdd_station_routing(const NetworkStruct<T>& sn) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = sn.nstations, R = sn.nclasses;
    Matrix<T> P(M, M, zero);
    for (std::size_t i = 0; i < M; ++i) {
        const std::size_t si = sn.stateful_of_station(i + 1) - 1;
        for (std::size_t j = 0; j < M; ++j) {
            const std::size_t sj = sn.stateful_of_station(j + 1) - 1;
            P(i, j) = sn.rt(si * R, sj * R);
        }
    }
    for (std::size_t i = 0; i < M; ++i) {
        T s = zero;
        for (std::size_t j = 0; j < M; ++j) s += P(i, j);
        if (std::fabs(num_traits<T>::to_double(s) - 1.0) > 1e-8)
            throw InputError(
                "SolverCTMC(mdd): the station-to-station routing chain is not stochastic at "
                "station " +
                std::to_string(i + 1) +
                "; the mdd method needs every completion to move the job to another station");
    }
    return P;
}

/** Disciplines under which every job present is in service, each with its own phase. */
inline bool mdd_is_shared(SchedStrategy s) {
    return s == SchedStrategy::PS || s == SchedStrategy::DPS || s == SchedStrategy::GPS ||
           s == SchedStrategy::INF;
}

/** Disciplines the count-plus-one-phase local state represents, single server only. */
inline bool mdd_is_nonpreemptive(SchedStrategy s, double servers) {
    return (s == SchedStrategy::FCFS || s == SchedStrategy::LCFS || s == SchedStrategy::SIRO ||
            s == SchedStrategy::HOL) &&
           servers == 1.0;
}


/**
 * The Petri-net route of the `mdd` method: `spn_mdd` supplies the reachable set
 * and the Kronecker descriptor, `mdd_mcd` aggregates, and the measures are read
 * back per place.
 *
 * The approximation is the same Eq. 5 as for a queueing network, and it is exact
 * on a product-form net, which `solver_nc_spn_analyzer` solves exactly and far
 * more cheaply -- the aggregation earns its place on the nets that have NO
 * product form.
 *
 * A net carries no per-station service rate, so `mdd_mcd` returns only the level
 * marginals. The token throughput is assembled here from the mode rates and
 * those marginals, under the SAME independence across levels that the
 * aggregation already assumes: it is the method's own approximation applied once
 * more, not a second one layered on top. TN counts FIRING EVENTS, which is what
 * the explicit CTMC path reports and what `sn_pn_avg_rates` converts to a token
 * rate afterwards.
 */
template <class T>
CtmcMddSolution<T> spn_route(const NetworkStruct<T>& sn, const CtmcOptions& opt,
                             const mdd::MddMcdOptions& mcdopt) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = sn.nstations, R = sn.nclasses;

    spn::SpnOptions mddopt;
    const spn::SpnResult<T> net = spn::spn_mdd<T>(sn, mddopt);
    const mdd::MddMcdResult<T> out = mdd::mdd_mcd<T>(net.mdds, net.desc, mcdopt);
    const std::size_t L = net.info.nplacelevels;

    CtmcMddSolution<T> sol;
    sol.avg.QN = Matrix<T>(M, R, zero);
    sol.avg.UN = Matrix<T>(M, R, zero);
    sol.avg.RN = Matrix<T>(M, R, zero);
    sol.avg.TN = Matrix<T>(M, R, zero);
    sol.avg.XN.assign(R, zero);
    sol.avg.CN.assign(R, zero);
    sol.avg.StartN = Matrix<T>(M, R, zero);
    sol.avg.PreemptN = Matrix<T>(M, R, zero);

    for (std::size_t pp = 0; pp < net.info.places.size(); ++pp) {
        const std::size_t ist = sn.nodes[net.info.places[pp] - 1].station;
        if (ist < 1) continue;
        sol.avg.QN(ist - 1, 0) = out.QLen[pp];
        sol.avg.UN(ist - 1, 0) = out.QLen[pp];   // a Place is an INF station: U = Q
    }

    // P(level l = v) from the converged level chains; mdd_mcd works in the
    // paper's orientation, paper level k = level K+1-l.
    std::vector<std::vector<T>> pl(L);
    for (std::size_t l = 0; l < L; ++l) {
        const std::size_t k = net.mdds.K - 1 - l;
        pl[l].assign(net.mdds.domain[l], zero);
        for (std::size_t r = 0; r < out.Mrows[k].size(); ++r)
            pl[l][out.Mrows[k][r].second] += out.pik[k][r];
        T tot = zero;
        for (std::size_t v = 0; v < pl[l].size(); ++v) tot += pl[l][v];
        if (tot > zero)
            for (std::size_t v = 0; v < pl[l].size(); ++v) pl[l][v] = T(pl[l][v] / tot);
    }

    for (std::size_t e = 0; e < net.info.modes.size(); ++e) {
        const spn::SpnMode<T>& mde = net.info.modes[e];
        if (mde.nph > 1) continue;               // no single rate; read the phase level
        // E[min(enabling degree, servers)] under independence across the inputs
        std::vector<std::size_t> lv;
        for (std::size_t l = 0; l < L; ++l)
            if (mde.enab[l] > 0) lv.push_back(l);
        T nsrv = num_traits<T>::from_int(1);
        if (!lv.empty()) {
            double kmax = std::numeric_limits<double>::infinity();
            for (std::size_t i = 0; i < lv.size(); ++i)
                kmax = std::min(kmax, std::floor((pl[lv[i]].size() - 1) / mde.enab[lv[i]]));
            if (std::isfinite(mde.srv)) kmax = std::min(kmax, mde.srv);
            nsrv = zero;
            for (int k = 1; k <= static_cast<int>(kmax); ++k) {
                T ge = num_traits<T>::from_int(1);
                for (std::size_t i = 0; i < lv.size(); ++i) {
                    const std::size_t l = lv[i];
                    const std::size_t thr = static_cast<std::size_t>(k * mde.enab[l]);
                    if (thr >= pl[l].size()) {
                        ge = zero;
                        break;
                    }
                    T s = zero;
                    for (std::size_t v = thr; v < pl[l].size(); ++v) s += pl[l][v];
                    ge = T(ge * s);
                }
                nsrv += ge;                      // E[min] = sum_k P(min >= k)
            }
        }
        const T x = T(mde.D1(0, 0) * nsrv);
        for (std::size_t l = 0; l < L; ++l) {
            if (mde.enab[l] <= 0) continue;
            const std::size_t ist = sn.nodes[net.info.places[l] - 1].station;
            if (ist >= 1) sol.avg.TN(ist - 1, 0) = T(sol.avg.TN(ist - 1, 0) + x);
        }
    }
    for (std::size_t i = 0; i < M; ++i)
        if (sol.avg.TN(i, 0) > zero)
            sol.avg.RN(i, 0) = T(sol.avg.QN(i, 0) / sol.avg.TN(i, 0));

    const std::size_t ref = sn.classes.empty() ? 0 : sn.classes[0].refstat;
    if (ref >= 1 && ref <= M) sol.avg.XN[0] = sol.avg.TN(ref - 1, 0);
    T Nk = zero;
    for (std::size_t i = 0; i < M; ++i) Nk += sol.avg.QN(i, 0);
    if (sol.avg.XN[0] > zero && Nk > zero) sol.avg.CN[0] = T(Nk / sol.avg.XN[0]);

    sol.num_states = static_cast<long long>(net.info.diagram.cardinality());
    sol.level_sizes = out.level_sizes;
    sol.iters = out.iters;
    sol.paths_per_level = out.paths_per_level;
    sol.no_aggregation = out.no_aggregation;
    sol.encoding = "spn";
    return sol;
}

}  // namespace mdd_analyzer_detail

/**
 * Can the `mdd` decision-diagram method be asked for this model?
 *
 * The model-shape gate asked as a predicate rather than thrown.
 * `solver_ctmc_mdd_analyzer` refuses with it before it builds anything, and a
 * REPORT reaches the same call, so a pair the report offers is a pair the
 * method runs. One predicate with two callers is what stops the two from
 * disagreeing about which models the method serves.
 *
 * A STOCHASTIC PETRI NET IS EXEMPT: a Place model is read through
 * `spn_route`, which builds the reachable set and the Kronecker descriptor from
 * the marking rather than from the (station,class) encoding, so neither the
 * single-class rule nor the closed-population rule applies to it.
 *
 * The deeper refusals the analyzer still raises -- a station-to-station chain
 * that is not stochastic, and a phase-type law at a discipline neither local
 * encoding represents -- are not restated here: they are decided from
 * quantities the analyzer computes on its way through, not from the model
 * shape, so a caller cannot be told about them without doing the work.
 *
 * @param sn the refreshed struct of the model
 * @return an empty string when the method may run, else the refusal
 */
template <class T>
std::string solver_ctmc_mdd_supports(const NetworkStruct<T>& sn) {
    for (std::size_t ind = 0; ind < sn.nodes.size(); ++ind)
        if (sn.nodes[ind].nodetype == NodeType::Place) return std::string();
    // A FORK-JOIN MODEL IS NEITHER of the two shapes this method serves. The
    // tag augmentation a fork needs adds one auxiliary class per branch, so the
    // struct that reaches the analyzer is never single-class however the model
    // was written, and the level decomposition has no meaning for a firing that
    // does not conserve the per-chain population. Also stated in
    // `qn::ctmc_feature_set("mdd")`, which drops the Fork/Join names.
    for (std::size_t ind = 0; ind < sn.nodes.size(); ++ind)
        if (sn.nodes[ind].nodetype == NodeType::Fork || sn.nodes[ind].nodetype == NodeType::Join)
            return "SolverCTMC(mdd): the mdd method analyses closed single-class networks and "
                   "stochastic Petri nets; a fork-join model is neither, and its tag "
                   "augmentation adds one auxiliary class per branch";
    if (sn.nclasses != 1)
        return "SolverCTMC(mdd): the mdd method analyses single-class networks; this model has " +
               std::to_string(sn.nclasses) +
               " classes. The Kronecker descriptor would need one level per (station,class)";
    for (std::size_t ind = 0; ind < sn.nodes.size(); ++ind)
        if (sn.nodes[ind].nodetype == NodeType::Source || sn.nodes[ind].nodetype == NodeType::Sink)
            return "SolverCTMC(mdd): the mdd method analyses CLOSED networks; an open stream "
                   "makes the marking unbounded, so the reachable set has no finite decision "
                   "diagram";
    const std::vector<double> njobs_gate = sn.njobs();
    if (!std::isfinite(njobs_gate[0]) || njobs_gate[0] <= 0.0)
        return "SolverCTMC(mdd): the mdd method needs a finite positive closed population";
    return std::string();
}

/**
 * Solve with the `mdd` method.
 *
 * @param sn the refreshed struct of a CLOSED single-class network
 * @param opt the SolverCTMC knobs; only `method` is read, the state-space ones
 *        having no meaning for a solve that enumerates nothing
 * @param mcdopt the level-iteration knobs, an INNER numerical solve whose
 *        tolerance must stay far tighter than any solver-level `iter_tol`
 */
template <class T>
CtmcMddSolution<T> solver_ctmc_mdd_analyzer(const NetworkStruct<T>& sn, const CtmcOptions& opt,
                                            const mdd::MddMcdOptions& mcdopt =
                                                mdd::MddMcdOptions()) {
    // EVERY STEP HERE IS A FIELD OPERATION, including the level solve.
    // `mcd_solve_stat` picks its backend by the arithmetic: Householder QR in
    // floating point, where the reflector norms are square roots and the point
    // is to avoid squaring the condition number, and `line::lstsq` under exact
    // arithmetic, where there is no condition number to square. So `--method
    // mdd` runs under `--arith exact` like `--method default`, and returns the
    // exact rational fixed point of the level system. The one thing exactness
    // does NOT buy is agreement with the enumerated chain: the level
    // aggregation is an approximation away from product form whatever the
    // arithmetic, and `no_aggregation` is what certifies otherwise.
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = sn.nstations, R = sn.nclasses;

    // The model-shape rules live in `solver_ctmc_mdd_supports`, which a report
    // asks too: the analyzer must refuse exactly what the report refuses, and
    // one predicate with two callers is what keeps the two from drifting apart.
    // All three come back as one exception type now; the population case used to
    // be an InputError, and no catch site anywhere distinguishes the two.
    const std::string shape_reason = solver_ctmc_mdd_supports(sn);
    if (!shape_reason.empty()) throw UnsupportedError(shape_reason);

    // A net holding Places takes the SPN translation instead of mdd_descriptor:
    // the levels are places, the measures come back per place, and the same Eq. 5
    // approximation applies. See mdd_analyzer_detail::spn_route.
    for (std::size_t ind = 0; ind < sn.nodes.size(); ++ind)
        if (sn.nodes[ind].nodetype == NodeType::Place)
            return mdd_analyzer_detail::spn_route<T>(sn, opt, mcdopt);

    const std::vector<double> njobs = sn.njobs();
    const int N = static_cast<int>(njobs[0]);

    // ---- service laws: the rate by LINE convention, the (D0,D1) pair when the
    // station is phase-type. A one-phase entry stays exponential so that the
    // descriptor takes its compact index = population encoding.
    std::vector<T> mu(M, zero);
    std::vector<double> servers(M, 1.0);
    std::vector<mdd::MddServiceLaw<T>> proc(M);
    std::vector<SchedStrategy> sched(M, SchedStrategy::FCFS);
    bool any_ph = false;
    for (std::size_t i = 0; i < M; ++i) {
        mu[i] = sn.rates(i, 0);
        servers[i] = sn.stations[i].nservers;
        sched[i] = sn.stations[i].sched;
        const lang::Distrib<T>& d = sn.service[i][0];
        if (d.phases() > 1 && d.D0.rows() > 0) {
            proc[i].D0 = d.D0;
            proc[i].D1 = d.D1;
            proc[i].present = true;
            any_ph = true;
        }
    }

    // ---- which local encoding represents these disciplines exactly
    std::string encoding = "np";
    if (any_ph) {
        bool all_shared = true;
        for (std::size_t i = 0; i < M; ++i)
            if (!mdd_analyzer_detail::mdd_is_shared(sched[i])) all_shared = false;
        if (all_shared) {
            encoding = "ps";
        } else {
            // Every PHASE-TYPE station must be non-preemptive for the compact
            // encoding; a shared-server one among them is the mixed case, which
            // no single encoding covers.
            bool all_np_at_ph = true;
            std::size_t first_ph = M, bad = M;
            for (std::size_t i = 0; i < M; ++i) {
                if (!proc[i].present) continue;
                if (first_ph == M) first_ph = i;
                if (mdd_analyzer_detail::mdd_is_nonpreemptive(sched[i], servers[i])) continue;
                all_np_at_ph = false;
                if (bad == M && !mdd_analyzer_detail::mdd_is_shared(sched[i])) bad = i;
            }
            if (!all_np_at_ph) {
                if (bad == M) bad = first_ph;
                throw UnsupportedError(
                    "SolverCTMC(mdd): station " + std::to_string(bad + 1) +
                    " combines a phase-type service law with a discipline that neither local "
                    "encoding represents: the count-plus-phase encoding is non-preemptive, and "
                    "the per-phase-count encoding covers only shared servers (PS/DPS/GPS/INF). "
                    "Mixing a shared and a non-preemptive phase-type station in one model is "
                    "likewise unsupported");
            }
        }
    }

    const Matrix<T> P = mdd_analyzer_detail::mdd_station_routing(sn);

    // ---- descriptor, reachable set, level aggregation
    mdd::MddDescriptor<T> desc;
    if (encoding == "ps") {
        desc = mdd::mdd_ps(mu, P, servers, N, proc);
    } else {
        std::vector<std::string> schedname(M);
        for (std::size_t i = 0; i < M; ++i) schedname[i] = lang::sched_to_text(sched[i]);
        desc = mdd::mdd_descriptor(mu, P, servers, N, proc, schedname);
    }
    const mdd::MDD diagram = mdd::mdd_reachset(desc.domain, desc.init, desc.nextfun);
    const mdd::MddMcdResult<T> out = mdd::mdd_mcd(diagram.to_struct(), desc, mcdopt);

    // ---- pack the analyzer contract
    CtmcMddSolution<T> s;
    s.avg.QN = Matrix<T>(M, R, zero);
    s.avg.UN = Matrix<T>(M, R, zero);
    s.avg.RN = Matrix<T>(M, R, zero);
    s.avg.TN = Matrix<T>(M, R, zero);
    s.avg.XN.assign(R, zero);
    s.avg.CN.assign(R, zero);
    for (std::size_t i = 0; i < M; ++i) {
        s.avg.QN(i, 0) = out.QLen[i];
        s.avg.UN(i, 0) = out.U[i];
        s.avg.TN(i, 0) = out.X[i];
        if (out.X[i] > zero) s.avg.RN(i, 0) = T(out.QLen[i] / out.X[i]);  // Little's law
    }

    // System throughput at the reference station, per unit visit. `visits` is
    // indexed by STATEFUL node, so the station index has to be converted; a
    // model with a Router has more nodes than stations and reading the station
    // index straight into that array would charge the wrong row.
    const std::size_t ref = sn.classes[0].refstat;
    T vis = num_traits<T>::from_int(1);
    if (!sn.visits.empty() && sn.visits[0].rows() > 0) {
        const std::size_t rsf = sn.stateful_of_station(ref) - 1;
        if (sn.visits[0](rsf, 0) > zero) vis = sn.visits[0](rsf, 0);
    }
    s.avg.XN[0] = T(s.avg.TN(ref - 1, 0) / vis);
    if (s.avg.XN[0] > zero) s.avg.CN[0] = T(num_traits<T>::from_int(N) / s.avg.XN[0]);

    s.num_states = diagram.cardinality();
    s.level_sizes = out.level_sizes;
    s.iters = out.iters;
    s.paths_per_level = out.paths_per_level;
    s.no_aggregation = out.no_aggregation;
    s.encoding = encoding;
    s.actualmethod = opt.method.empty() ? std::string("mdd") : opt.method;
    return s;
}

/** Solve with `mdd` and format, so a caller with no use for the diagram has one call. */
template <class T>
mva::AvgResult<T> solver_ctmc_mdd_run_analyzer(const NetworkStruct<T>& sn, const CtmcOptions& opt,
                                      const mdd::MddMcdOptions& mcdopt = mdd::MddMcdOptions()) {
    const CtmcMddSolution<T> s = solver_ctmc_mdd_analyzer(sn, opt, mcdopt);
    CtmcSolution<T> d;
    d.avg = s.avg;
    d.actualmethod = s.actualmethod;
    return solver_ctmc_avg_table(sn, d, opt.method);
}

}  // namespace ctmc
}  // namespace line

#endif  // LINE_SOLVERS_CTMC_SOLVER_CTMC_MDD_ANALYZER_H
