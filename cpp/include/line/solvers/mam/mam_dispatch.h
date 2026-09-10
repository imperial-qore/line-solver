/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MAM_MAM_DISPATCH_H
#define LINE_SOLVERS_MAM_MAM_DISPATCH_H

/**
 * Port of `solver_mam_analyzer.m`: one inner solve, choosing the analyzer that
 * fits the model and the requested method.
 *
 * THE ORDER IS THE CONTRACT, exactly as in `mva_dispatch.h`. The reference's
 * sequence, top to bottom:
 *
 *  -1  discrete-time (slotted) models -> solver_mam_dt              <- ported
 *   0  exact MAP/MAP/1 fast path, tried before any method dispatch  <- ported
 *   1  method 'dec.mmap'                       -> solver_mam        <- ported
 *   2  method 'default' / 'dec.source':
 *        2a homogeneous Fork-Join              -> solver_mam_fj  ("qiu") <- ported
 *        2b any other Fork-Join, open          -> solver_mam_basic_mmap <- ported
 *        2c BMAP/PH/N/N bufferless retrial     -> solver_mam_retrial <- ported
 *        2d reneging (MAP/M/s+G)               -> solver_mam_retrial
 *        2e single-class closed Delay+Queue    -> solver_mam_ldqbd   <- ported
 *        2f otherwise                          -> solver_mam_basic   <- ported
 *   3  method 'dec.poisson'  -> solver_mam_basic with space_max = 1  <- ported
 *   4  method 'mna'          -> solver_mna_open / solver_mna_closed        <- ported
 *   5  method 'ldqbd'        -> solver_mam_ldqbd                     <- ported
 *   5b method 'bgchain'      -> solver_mam_bgchain                   <- ported
 *   6  methods 'inap' / 'inapplus' / 'inapinf' -> moved to SolverAG (ag_dispatch.h)
 *   6b method 'exact'        -> removed from the reference (SolverMAM.m:39-40:
 *        "'exact' method removed - autocat moved to line-legacy.git"); not in
 *        list_valid_methods either, so unreachable through solver_mam_solve
 *   7  method 'dec.source.mmap' -> solver_mam_basic_mmap             <- ported
 *
 * 2a IS DOUBLE ONLY. The FJ_codes engine behind it needs an ordered real Schur
 * factorization and two Bartels-Stewart Sylvester solves, i.e. LAPACK, so
 * `solver_mam_fj` refuses at Rational and Real AFTER running the reference's two
 * validation steps -- a model outside the homogeneous class is still named as
 * such at every arithmetic.
 *
 * WHAT IS NOT PORTED IS REFUSED BY NAME and never allowed to fall through to
 * `solver_mam_basic`. A fork-join model or a level-dependent closed model solved
 * as an ordinary open decomposition returns numbers that are simply not the
 * model's. The rule extends ONE step past the reference at 2b: a CLOSED
 * fork-join model, which the reference itself lets fall through to
 * `solver_mam_basic`, is refused here for exactly that reason; see the branch.
 *
 * AFTER the analyzer, the reference overwrites the throughput of every EXT
 * (Source) station with `sn.rates`, and zeroes every NaN across its six metrics
 * (the four matrices Q, U, R, Tp and the two vectors C, X). Both are reproduced
 * in `finish_dispatch`, which every branch below step 0 returns through. Step 0
 * bypasses it, as the reference's own early return does; see `finish_dispatch`.
 * The C++ sweep tests `!isfinite` rather than NaN alone, so it also zeroes an
 * infinity the reference would keep -- deliberate, because the `disabled` flags
 * this port carries mean a non-finite entry here is a division artefact and
 * never a modelled infinity.
 */

#include <limits>
#include <cmath>
#include <string>
#include <vector>

#include "line/api/sn/sn_is_discrete_time.h"
#include "line/api/sn/sn_predicates.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mam/mam_types.h"
#include "line/solvers/mam/solver_mam_basic.h"
#include "line/solvers/mam/solver_mam_basic_mmap.h"
#include "line/solvers/mam/solver_mam_decmmap.h"
#include "line/solvers/mam/solver_mam_fj.h"
#include "line/solvers/mam/solver_mam_bgchain.h"
#include "line/solvers/mam/solver_mam_ldqbd.h"
#include "line/solvers/mam/solver_mam_dt.h"
#include "line/solvers/mam/solver_mam_mapmap1_exact.h"
#include "line/solvers/mam/solver_mam_retrial.h"
#include "line/solvers/mam/solver_mna.h"
#include "line/util/error.h"

namespace line {
namespace mam {

namespace dispatch_detail {

/** Is there any Fork node? A Join alone cannot occur without one. */
template <class T>
bool has_fork(const qn::NetworkStruct<T>& L) {
    return L.has_fork();
}

/**
 * The reference's `sn_has_fork_join`, which is `any(sn.fj(:) > 0)` and NOT a
 * scan for Fork nodes: it is the branch-2b gate, and a Fork whose Join the
 * model never declared leaves `fj` empty and fails it.
 */
template <class T>
bool has_fork_join(const qn::NetworkStruct<T>& L) {
    return !L.fj.empty();
}

/**
 * `sn_is_closed_model`, which is `all(isfinite(sn.njobs))`. NetworkStruct
 * carries `is_open_model` but no counterpart, and the mna branch needs both:
 * the reference's own else is the mixed-model refusal, so "not open" is not the
 * same test as "closed". The empty-class guard mirrors `is_open_model`, where
 * MATLAB's `all([])` would answer true for both.
 */
template <class T>
bool is_closed_model(const qn::NetworkStruct<T>& L) {
    for (const qn::JobClass& c : L.classes)
        if (!std::isfinite(c.population)) return false;
    return !L.classes.empty();
}

/** The reference's `isClosedDelayQueue`: one class, two stations, one INF one FCFS. */
template <class T>
bool is_closed_delay_queue(const qn::NetworkStruct<T>& L) {
    if (L.nclasses != 1 || L.nstations != 2) return false;
    if (!std::isfinite(L.classes[0].population)) return false;
    std::size_t ndelay = 0, nqueue = 0;
    for (const qn::Station<T>& st : L.stations) {
        if (st.sched == lang::SchedStrategy::INF) ++ndelay;
        else if (st.sched == lang::SchedStrategy::FCFS) ++nqueue;
    }
    return ndelay == 1 && nqueue == 1;
}

/**
 * Whether a setup/delay-off station, if there is one, is inside the exact
 * regime of `solver_mam_ldqbd`: one server, exponential service, no load
 * dependence. True when the model declares no setup at all.
 */
template <class T>
bool ldqbd_setup_ok(const qn::NetworkStruct<T>& L) {
    if (L.setupparam.empty()) return true;
    for (typename std::map<std::size_t, qn::SetupDelayOffParam<T>>::const_iterator it =
             L.setupparam.begin();
         it != L.setupparam.end(); ++it) {
        const std::size_t ist = it->first;
        if (ist == 0 || ist > L.stations.size()) return false;
        lang::Distrib<T> su, doff;
        if (!it->second.last(su, doff) || doff.disabled) continue;  // not actually declared
        const qn::Station<T>& st = L.stations[ist - 1];
        if (st.nservers > 1) return false;
        // NETWORK-WIDE, as the MATLAB reference's sn_has_load_dependence: the
        // exact chain assumes a plain delay away from this station, so a
        // load-dependent station ANYWHERE takes the model out of the regime, not
        // only a load-dependent setup station. Reading st.lldscaling alone
        // admitted a load-dependent Delay that MATLAB and the JAR refuse.
        if (api::sn_has_load_dependence(L)) return false;
        if (L.service.size() < ist || L.service[ist - 1].empty()) return false;
        if (lang::dist_to_map(L.service[ist - 1][0]).D0.rows() != 1) return false;
    }
    return true;
}

/**
 * The reference's `mnaApplies`: whether the closed MNA analyzer covers this
 * model. The round-robin split is carried by the open traffic equations only, a
 * self-looping class has no inter-station flow to decompose, and a station whose
 * discipline the flow sweep does not update would keep a zero queue length. Any
 * of those routes `default` to dec.source instead.
 */
template <class T>
bool mna_applies(const qn::NetworkStruct<T>& L) {
    // solver_mna_closed drives its bisection over CLASSES but stores the throughput in
    // the CHAIN-indexed lambda, and renormalizes chain c's queue lengths with the
    // class-indexed njobs(c). Both are only correct when each chain holds exactly one
    // class, and the analyzer refuses otherwise by name -- so the DEFAULT has to stop
    // here rather than propagate that refusal (cqn_twoclass_hyperl, 1 chain over 2
    // classes, died at MAM under lang='cpp' while MATLAB answered dec.source).
    if (L.nchains != L.nclasses) return false;

    for (const qn::NodeDef& nd : L.nodes)
        for (std::size_t r = 0; r < nd.routing.size(); ++r)
            if (nd.routing[r] == lang::RoutingStrategy::RROBIN) return false;

    for (const qn::Station<T>& st : L.stations) {
        if (st.sched != lang::SchedStrategy::INF && st.sched != lang::SchedStrategy::PS &&
            st.sched != lang::SchedStrategy::FCFS && st.sched != lang::SchedStrategy::EXT)
            return false;
        // The PS branch of solver_mna_closed forms U = S*T and the geometric bound
        // from it WITHOUT dividing by the number of servers, so a multiserver PS
        // station is misrepresented; dec.source is exact on the non-queueing regime
        // (c >> N) that shape usually stands for.
        if (st.sched == lang::SchedStrategy::PS && st.nservers > 1) return false;
        // A multiclass FCFS station makes the flow sweep superpose one MMAP per class
        // and then solve MMAPPH1FCFS at level sum(N)+1: measured against an exact CTMC
        // that costs 12-61s where dec.source costs 0.03s and is not more accurate
        // (mean relative error 0.11-0.38 against 0.04-0.19). The single-class case is
        // both cheap and better, so keep only that one.
        if (st.sched == lang::SchedStrategy::FCFS && L.nclasses > 1) return false;
    }

    const Matrix<T> V = basic_detail::station_visits(L);
    for (std::size_t k = 0; k < L.nclasses; ++k) {
        if (!std::isfinite(L.classes[k].population)) continue;
        std::size_t seen = 0, at = 0;
        for (std::size_t i = 0; i < L.nstations; ++i)
            if (num_traits<T>::to_double(V(i, k)) > 1e-8) {
                ++seen;
                at = i;
            }
        if (seen == 1 && L.stations[at].sched != lang::SchedStrategy::INF &&
            L.stations[at].sched != lang::SchedStrategy::EXT)
            return false;
    }
    return true;
}

/**
 * The reference's `mnaConserves`: whether the closed MNA outer bisection closed on
 * N. `solver_mna_closed` rescales each chain onto its population as a last step, so
 * a diverged bisection still returns queue lengths that sum to N and the failure is
 * invisible in Q. R and Tp are NOT rescaled, so Little's law over the whole network,
 * sum_i Tp(i,k)*R(i,k) = N_k, still reads the raw iterate: a converged run lands
 * within 5e-4 of N and a diverged one is orders of magnitude out, or negative.
 */
template <class T, class Sol>
bool mna_conserves(const qn::NetworkStruct<T>& L, const Sol& sol) {
    for (std::size_t i = 0; i < L.nstations; ++i)
        for (std::size_t k = 0; k < L.nclasses; ++k)
            if (!std::isfinite(num_traits<T>::to_double(sol.Q(i, k))) ||
                !std::isfinite(num_traits<T>::to_double(sol.R(i, k))) ||
                !std::isfinite(num_traits<T>::to_double(sol.Tp(i, k))))
                return false;

    for (std::size_t k = 0; k < L.nclasses; ++k) {
        const double nk = L.classes[k].population;
        if (!std::isfinite(nk) || nk <= 0) continue;
        double npred = 0;
        for (std::size_t i = 0; i < L.nstations; ++i)
            npred += num_traits<T>::to_double(sol.Tp(i, k)) * num_traits<T>::to_double(sol.R(i, k));
        if (std::fabs(npred - nk) > 0.01 * nk) return false;
    }
    return true;
}

/**
 * The reference's `bgchainApplies`: no class priorities, no fork-join, and at
 * least one station visited by the closed classes, which is what the background
 * chain is built over.
 */
template <class T>
bool bgchain_applies(const qn::NetworkStruct<T>& L, const MamOptions& opt) {
    bool prio_sched = false;
    for (const qn::Station<T>& st : L.stations)
        if (st.sched == lang::SchedStrategy::HOL || st.sched == lang::SchedStrategy::FCFSPRPRIO)
            prio_sched = true;
    if (prio_sched)
        for (std::size_t k = 1; k < L.nclasses; ++k)
            if (L.classes[k].prio != L.classes[0].prio) return false;

    if (has_fork_join(L)) return false;

    const mva::ChainDemands<T> dem = mva::sn_get_demands_chain(L);
    bool visited = false;
    for (std::size_t c = 0; c < L.nchains && !visited; ++c) {
        if (!std::isfinite(num_traits<T>::to_double(dem.Nchain[c])) ||
            num_traits<T>::to_double(dem.Nchain[c]) <= 0.0)
            continue;
        for (std::size_t i = 0; i < L.nstations; ++i)
            if (num_traits<T>::to_double(dem.Vchain(i, c)) > 1e-14) {
                visited = true;
                break;
            }
    }
    if (!visited) return false;

    // The background chain enumerates the closed population vector over the
    // stations the closed classes visit, so a large closed population makes
    // mam_bgchain_ctmc refuse the model outright. That refusal is right when the
    // user asked for bgchain by name and wrong as a default, which must land on
    // a method that answers: size the chain first and leave those models to
    // dec.source. The limit is the one mam_bgchain_ctmc enforces.
    const double states_max = (opt.bgstates_max > 0) ? static_cast<double>(opt.bgstates_max) : 20000.0;
    return bgchain_states(L, opt) <= states_max;
}

/**
 * The reference's `bgchainClosedExact`: whether the background chain represents
 * this closed model's service laws exactly.
 *
 * A station that is not PS or INF must satisfy BOTH conditions below, because the
 * background chain makes two separate first-moment substitutions there.
 *
 *   1  mam_bgchain_ctmc builds its generator from the MEAN service time alone.
 *      That is exact at a PS or INF station, which is insensitive to the service
 *      law beyond its first moment, and exact under any discipline when the law
 *      IS exponential. Measured on a closed Delay+FCFS cycle with Erlang-3
 *      service, the mean-only chain reads 2.8% off SolverCTMC.
 *   2  The capacity a station's closed jobs hold is split over the background
 *      classes in proportion to their COUNTS, which is service in random order.
 *      That is exact under PS, and exact under FCFS only when the classes are
 *      served at the SAME rate -- an FCFS station with class-dependent rates
 *      reads 25.2% off SolverCTMC on a two-chain closed cycle, against 3.5e-16
 *      when the two rates are made equal.
 *
 * solver_mna_closed carries the phase-type representation instead, so neither
 * surrogate may be chosen as the DEFAULT. Asking for bgchain by name still gets
 * it, with both approximations documented in solver_mam_bgchain.
 */
template <class T>
bool bgchain_closed_exact(const qn::NetworkStruct<T>& L) {
    for (std::size_t i = 0; i < L.nstations; ++i) {
        const lang::SchedStrategy sc = L.stations[i].sched;
        if (sc == lang::SchedStrategy::INF || sc == lang::SchedStrategy::PS ||
            sc == lang::SchedStrategy::EXT)
            continue;
        double rate_here = std::numeric_limits<double>::quiet_NaN();
        for (std::size_t r = 0; r < L.nclasses; ++r) {
            const double rate = num_traits<T>::to_double(L.rates(i, r));
            if (!std::isfinite(rate) || rate <= 0.0) continue;
            if (L.service[i][r].type != lang::ProcessType::EXP) return false;
            if (std::isnan(rate_here)) {
                rate_here = rate;
            } else if (std::fabs(rate - rate_here) > 1e-3 * rate_here) {
                return false;
            }
        }
    }
    return true;
}

/** Is any station a Cache, a Place or a Transition, i.e. outside the MAM envelope? */
template <class T>
void reject_unsupported_nodes(const qn::NetworkStruct<T>& L) {
    for (const qn::NodeDef& nd : L.nodes) {
        switch (nd.nodetype) {
            case qn::NodeType::Cache:
                throw UnsupportedError(
                    "SolverMAM: Cache nodes are outside the MAM feature set; use SolverMVA, "
                    "SolverNC, SolverCTMC or SolverLDES");
            case qn::NodeType::Place:
            case qn::NodeType::Transition:
                throw UnsupportedError(
                    "SolverMAM: stochastic Petri net models are outside the MAM feature set; "
                    "use SolverCTMC, SolverSSA or SolverLDES");
            default:
                break;
        }
    }
}

}  // namespace dispatch_detail


/**
 * What `solver_mam_analyzer.m` does AFTER whichever analyzer ran: pin the
 * throughput of every EXT (Source) station to its declared rate, and zero every
 * non-finite entry.
 *
 * Factored out because THREE branches return early and must still take this
 * tail: 2a (`qiu`), 2c (`retrial`) and the 2e `ldqbd` preference on `default`.
 * The reference has no early return at any of them -- its `case` bodies fall
 * through to the tail -- so routing them here is what keeps the C++ ladder
 * equivalent rather than an optimization.
 *
 * THE ONE EARLY RETURN THAT DOES NOT COME HERE is the exact MAP/MAP/1 fast
 * path, and it matches the reference: `solver_mam_analyzer.m` lines 10-17
 * return before the tail as well. It is not an omission -- that analyzer sets
 * the Source throughput itself (`solver_mam_mapmap1_exact.h`, `Tp(src,0) =
 * lambda`), so there is nothing for the tail to pin.
 */
template <class T>
MamSolution<T>& finish_dispatch(const qn::NetworkStruct<T>& L, MamSolution<T>& out) {
    using lang::SchedStrategy;
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < L.nstations; ++i)
        if (L.stations[i].sched == SchedStrategy::EXT)
            for (std::size_t r = 0; r < L.nclasses; ++r)
                out.sol.Tp(i, r) = L.disabled[i][r] ? zero : L.rates(i, r);
    // The reference's terminal NaN sweep. The C++ layer carries `disabled`
    // flags rather than NaN sentinels, so nothing here should be non-finite;
    // the sweep stays because a division by a zero iterate can still produce
    // one at double, and reporting NaN would be worse than reporting zero.
    Matrix<T>* mats[4] = {&out.sol.Q, &out.sol.U, &out.sol.R, &out.sol.Tp};
    for (Matrix<T>* m : mats)
        for (std::size_t i = 0; i < m->rows(); ++i)
            for (std::size_t j = 0; j < m->cols(); ++j)
                if (!std::isfinite(num_traits<T>::to_double((*m)(i, j)))) (*m)(i, j) = zero;
    for (T& v : out.sol.C)
        if (!std::isfinite(num_traits<T>::to_double(v))) v = zero;
    for (T& v : out.sol.X)
        if (!std::isfinite(num_traits<T>::to_double(v))) v = zero;
    return out;
}

/** The ladder. */
template <class T>
MamSolution<T> mam_dispatch(const qn::NetworkStruct<T>& L, const MamOptions& opt_in) {
    using lang::SchedStrategy;
    MamOptions opt = opt_in;
    MamSolution<T> out;

    dispatch_detail::reject_unsupported_nodes(L);

    // -1. the discrete-time (slotted) path, ahead of EVERYTHING. It has to run
    // before any phase-type conversion: by the time a law reaches `proc` a
    // Geometric has already been fitted to a CONTINUOUS MAP and the lattice is
    // no longer visible, so the test reads procid/rates/scv instead.
    if constexpr (num_traits<T>::has_transcendental) {
        api::DiscreteTimeOptions dtopt;
        dtopt.timescale = opt.timescale;
        dtopt.slotlength = opt.slotlength;
        double slot = 1.0;
        api::DiscreteTimeInfo dtinfo;
        if (api::sn_is_discrete_time(L, dtopt, &slot, &dtinfo)) {
            out = solver_mam_dt(L, opt, slot);
            out.sol.method = opt.method;
            return out;
        }
    }

    // 0. the exact MAP/MAP/1 fast path, ahead of any method dispatch. A station
    // that powers down is NOT an M/M/1 whatever its shape looks like: the setup
    // is extra work the QBD below does not carry, so it must not be taken here.
    if (L.setupparam.empty()) {
        const MapMap1Exact<T> ex = solver_mam_mapmap1_exact(L);
        if (ex.ok) {
            out.sol = ex.sol;
            out.sol.method = opt.method;
            out.actualmethod = "exact.mapmap1";
            return out;
        }
    }

    const std::string& method = opt.method;

    if (method == "dec.mmap") {
        out.sol = solver_mam_decmmap(L, opt);
        out.actualmethod = "dec.mmap";
        return finish_dispatch(L, out);
    }

    if (method == "default" || method == "dec.source" || method == "dec.poisson") {
        if (method == "dec.poisson") opt.space_max = 1;
        if (method != "dec.poisson") {
          if constexpr (!num_traits<T>::has_transcendental) {
            if (dispatch_detail::has_fork(L) || dispatch_detail::has_fork_join(L))
                throw UnsupportedError(
                    "SolverMAM: model '" + L.name +
                    "' forks, and every fork-join route here first classifies the branches by "
                    "fitting a phase-type to each, which is transcendental and has no exact "
                    "counterpart; use --arith double or real");
          } else {
            // 2a: `fj_is_homogeneous` -- one fork-join pair, K identical open
            // FCFS/PS branches -- is the ONLY class FJ_codes is defined on, and
            // most fork-join models fail it, in which case 2b takes them. See
            // solver_mam_fj.h for why no other FJ approximation in this tree may
            // be substituted under this method name.
            // THE PREDICATE IS NOT ARITHMETIC-NEUTRAL. mam_fj_is_homogeneous
            // compares branches by fitting a PH to each (dist_to_map ->
            // aph_fit), which static_asserts on transcendental arithmetic, so
            // calling it unconditionally made the WHOLE MAM dispatch
            // uninstantiable at Rational -- a model with no Fork at all stopped
            // compiling. It is therefore reached only from the exact-capable
            // branch below, and exact arithmetic refuses a forking model BY
            // NAME rather than skipping the gate and answering it as ordinary.
            const MamFjInfo fjinfo = mam_fj_is_homogeneous(L);
            if (fjinfo.ok) {
                out.sol = solver_mam_fj(L, opt);
                out.actualmethod = "qiu";
                return finish_dispatch(L, out);
            }
            // 2b: every other OPEN fork-join model, which the reference sends to
            // the MMAP decomposition.
            if (dispatch_detail::has_fork_join(L) && L.is_open_model()) {
                out.sol = solver_mam_basic_mmap(L, opt);
                out.actualmethod = "dec.source.mmap";
                return finish_dispatch(L, out);
            }
            // DELIBERATE DEVIATION, and the only one in this ladder. A CLOSED
            // fork-join model passes both reference gates and falls through to
            // solver_mam_basic, whose dec.source decomposition has no Join
            // synchronization at all: it zeroes the Join's own metrics and
            // charges no waiting for the slowest sibling branch, so the response
            // times it reports are the branches' and not the model's, with
            // nothing in the output to say so. Refused by name instead.
            //
            // SANCTIONED (2026-07-28), not provisional: the refusal is preferred
            // over bit-parity with MATLAB because the alternative here is a
            // SILENTLY WRONG answer, not a less accurate one. Do not "restore
            // parity" by deleting this without reopening that decision.
            if (dispatch_detail::has_fork(L))
                throw UnsupportedError(
                    "SolverMAM: model '" + L.name +
                    "' carries a Fork but is not an open fork-join network, and "
                    "solver_mam_analyzer.m has no branch for it: it falls through to "
                    "solver_mam_basic, whose dec.source decomposition never charges the "
                    "synchronization delay at the Join, so the response times would be those of "
                    "the branches alone. Call method 'dec.source.mmap', which routes the model to "
                    "the MMAP decomposition's closed wrapper and does charge the join, or use "
                    "SolverMVA, SolverNC, SolverJMT or SolverLDES");
          }

            // 2c: the BMAP/PH/N/N bufferless retrial queue of Dudin et al. A
            // station is claimed only when it declares an orbit or is
            // bufferless AND carries a RETRIAL drop rule, exactly as
            // qsys_is_retrial.m requires; nothing is inferred from capacity
            // alone, which would claim finite buffers that are not orbits.
            // The DETECTION is arithmetic-neutral (it reads retrialparam and
            // the drop rule), but the ANALYZER fits a phase-type to the arrival
            // and service processes at solver_mam_retrial.h:360, so the orbit
            // route is transcendental even though qsys_bmapphnn_retrial itself
            // is exact. Refuse the orbit BY NAME under an exact type rather
            // than fall through to 2e/2f, which would answer it as an ordinary
            // waiting line -- the very substitution wiring this branch removed.
            const MamRetrialInfo retinfo = mam_retrial_detect(L);
            if (retinfo.ok) {
                if constexpr (!num_traits<T>::has_transcendental) {
                    throw UnsupportedError(
                        "SolverMAM: model '" + L.name +
                        "' declares a retrial orbit, whose analyzer fits a phase-type to the "
                        "arrival and service processes and has no exact counterpart; use --arith "
                        "double or real");
                } else {
                    out.sol = solver_mam_retrial(L, opt);
                    out.actualmethod = "retrial";
                    return finish_dispatch(L, out);
                }
            }
            // 2d: reneging (`hasRenegingPatience`) reaches the same analyzer's
            // OTHER route, the MAP/M/s+G fluid queue. Its gate reads
            // sn.impatienceClass, sn.patienceProc and ImpatienceType.RENEGING,
            // none of which NetworkStruct carries, so no model this port can
            // build enters it. Recorded rather than approximated: a waiting line
            // whose jobs abandon is not an orbit, and 2c must not claim it.
            //
            // 2e: the level-dependent QBD is EXACT for this shape, and the
            // reference prefers it over dec.source on 'default'. The OPEN
            // Source+Queue regime of solver_mam_ldqbd is deliberately NOT
            // claimed here: it truncates at options.cutoff, so it is not
            // unconditionally better than dec.source and stays opt-in through
            // method='ldqbd', exactly as the reference's isClosedDelayQueue
            // comment states.
            // A closed setup task IS a Delay+Queue tandem by shape, and ldqbd
            // now models it exactly -- but only at a single server with
            // exponential service and no load dependence, which is what
            // qbd_setupdelayoff_closed covers. Outside that it refuses by name,
            // so the shape test has to exclude those or the default would reach
            // the refusal instead of falling through to the decomposition.
            if (method == "default" && dispatch_detail::ldqbd_setup_ok(L) &&
                dispatch_detail::is_closed_delay_queue(L)) {
                out.sol = solver_mam_ldqbd(L, opt).sol;
                out.actualmethod = "ldqbd";
                return finish_dispatch(L, out);
            }

            // 2e-bis: a closed model is the degenerate case of the background chain.
            // With no open work to take a share of the servers the chain is the EXACT
            // closed CTMC at chain granularity, so it dominates the mna fixed point
            // wherever its state space fits. bgchain_applies sizes that state space
            // and bgchain_closed_exact checks the chain is built from a service law it
            // represents exactly.
            // bgchain DROPS setup/delay-off: it would answer with the
            // always-warm chain and nothing would say so. Only solver_mam_basic
            // and solver_mam_ldqbd read it, so a setup model must not reach here.
            if (method == "default" && dispatch_detail::is_closed_model(L) &&
                L.setupparam.empty() &&
                dispatch_detail::bgchain_applies(L, opt) &&
                dispatch_detail::bgchain_closed_exact(L)) {
                out.sol = solver_mam_bgchain(L, opt);
                out.actualmethod = "bgchain";
                return finish_dispatch(L, out);
            }

            // 2f: a closed model has no arrival stream for dec.source to build its
            // Poisson surrogate from -- it replaces each closed chain by a source at
            // the current throughput iterate and never enforces the population, so
            // the answer neither conserves N nor separates the classes. MNA closes
            // the same traffic equations by bisecting the per-class throughput
            // against N.
            // mna DROPS setup/delay-off for the same reason as bgchain.
            if (method == "default" && dispatch_detail::is_closed_model(L) &&
                L.setupparam.empty() &&
                dispatch_detail::mna_applies(L)) {
                out.sol = solver_mna_closed(L, opt);
                out.actualmethod = "mna";
                if (!dispatch_detail::mna_conserves(L, out.sol)) {
                    // The outer bisection did not close on the population: the last
                    // step rescales each chain onto N regardless, so the failure is
                    // invisible in Q alone and only Little's law on the unrescaled R
                    // and Tp still shows it.
                    out.sol = solver_mam_basic(L, opt);
                    out.actualmethod = "dec.source";
                }
                return finish_dispatch(L, out);
            }

            // 2g: mixed. The closed classes are solved exactly as a background chain
            // and the open ones as QBDs driven by it, which is 4-5 significant
            // digits against CTMC where dec.source is 10-24% out.
            // The mixed branch drops setup/delay-off exactly as the closed one
            // does, so it takes the same guard.
            if (method == "default" && !L.is_open_model() && L.setupparam.empty() &&
                !dispatch_detail::is_closed_model(L) && dispatch_detail::bgchain_applies(L, opt)) {
                out.sol = solver_mam_bgchain(L, opt);
                out.actualmethod = "bgchain";
                return finish_dispatch(L, out);
            }
        }
        out.sol = solver_mam_basic(L, opt);
        out.actualmethod = (method == "default") ? "dec.source" : method;
    } else if (method == "mna") {
        // The two analyzers are disjoint by population regime, so unlike the
        // fork-join and retrial branches this pair carries no order hazard; the
        // reference's open-then-closed sequence is kept anyway.
        if (L.is_open_model()) {
            out.sol = solver_mna_open(L, opt);
        } else if (dispatch_detail::is_closed_model(L)) {
            out.sol = solver_mna_closed(L, opt);
        } else {
            // The reference's own line_error. `SolverMAM.supportsModelMethod`
            // already rejects a mixed model one level up, in
            // runner_detail::check_model_method, so this is reachable only by
            // calling mam_dispatch directly; it is transcribed rather than
            // dropped because that caller exists.
            throw UnsupportedError(
                "SolverMAM: the mna method in SolverMAM does not support mixed models");
        }
        out.actualmethod = "mna";
    } else if (method == "bgchain") {
        // Mixed networks: the closed classes are a background modulating chain,
        // the open classes are QBDs driven by it. With several closed chains an
        // outer iteration tags one chain at a time and aggregates the rest, so
        // the background chain always carries two classes.
        out.sol = solver_mam_bgchain(L, opt);
        out.actualmethod = "bgchain";
    } else if (method == "ldqbd") {
        out.sol = solver_mam_ldqbd(L, opt).sol;
        out.actualmethod = "ldqbd";
    } else if (method == "inap" || method == "inapplus" || method == "inapinf"
               || method == "exact") {
        // The RCAT methods moved to SolverAG. Name them rather than reporting an
        // unknown method, so a caller carrying an old options.method is told
        // where they went.
        throw UnsupportedError(
            "SolverMAM: the '" + method + "' method moved to SolverAG: RCAT decomposes the "
            "model into cooperating agents rather than decomposing traffic, and no MAM "
            "algorithm shares its machinery. Solve it with -s ag, or call "
            "line::ag::solver_ag directly");
    } else if (method == "retrial") {
        // The BMAP/PH/N/N retrial analyzer BY NAME, which branch 2c also
        // resolves to from 'default' on a retrial topology. The reference
        // advertises the name (SolverMAM.m) and dispatches it
        // (solver_mam_analyzer.m 'case retrial'), so it must be reachable here
        // and must refuse off that topology rather than answer a waiting line.
        {
            // The predicate `check_model_method` and `auto_family_refusal` ask,
            // so the method the report offers and the method that runs are the
            // same set, with the same sentence when they are not.
            const std::string retrialWhy = mam_retrial_refusal(L);
            if (!retrialWhy.empty()) throw UnsupportedError(retrialWhy);
            if constexpr (!num_traits<T>::has_transcendental) {
                throw UnsupportedError(
                    "SolverMAM: the retrial analyzer fits a phase-type to the arrival and "
                    "service processes and has no exact counterpart; use --arith double or real");
            } else {
                out.sol = solver_mam_retrial(L, opt);
                out.actualmethod = "retrial";
            }
        }
    } else if (method == "dec.source.mmap") {
        // solver_mam_basic_mmap.m: the inner MMAP decomposition for an open
        // model, the throughput bisection around it for a closed one. Unlike
        // branch 2b this is reachable BY NAME on a model with no fork at all,
        // which is what the reference's own dispatcher allows.
        out.sol = solver_mam_basic_mmap(L, opt);
        out.actualmethod = "dec.source.mmap";
    } else {
        throw UnsupportedError("SolverMAM: unknown method '" + method + "'");
    }

    return finish_dispatch(L, out);
}

}  // namespace mam
}  // namespace line

#endif  // LINE_SOLVERS_MAM_MAM_DISPATCH_H
