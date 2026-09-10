/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MAM_SOLVER_MAM_RETRIAL_H
#define LINE_SOLVERS_MAM_SOLVER_MAM_RETRIAL_H

/**
 * Port of `solver_mam_retrial.m`, the customer-impatience analyzer of SolverMAM.
 *
 * The reference file carries TWO analyzers behind one name, and tries them in
 * this order:
 *
 *   1. RENEGING, a MAP/M/s+G queue whose waiting jobs abandon after a generally
 *      distributed patience (Gursoy, Mehr and Akar). Solved as a multi-regime
 *      Markov fluid queue through `MRMFQSolver`.
 *   2. RETRIAL, the BMAP/PH/N/N bufferless orbit queue with flexible retrials
 *      admission control (Dudin et al., Mathematics 2025, 13(9), 1434). Solved
 *      by `qsys_bmapphnn_retrial`, which IS ported (api/qsys).
 *
 * ONLY ROUTE 2 IS REACHABLE FROM C++, AND ROUTE 1 IS REFUSED BY NAME. The
 * reneging gate reads `sn.impatienceClass`, `sn.patienceProc` and
 * `ImpatienceType.RENEGING`, none of which the C++ `NetworkStruct` carries, so
 * no model this port can BUILD is a reneging model. `mam_retrial_detect` says so
 * explicitly instead of letting such a model fall through to the retrial
 * generator, which would answer a queue with a waiting line as if it had an
 * orbit.
 *
 * WHAT THE ORBIT IS. A retrial station has no waiting room: the buffer equals
 * the server count. A job that finds every server busy joins an ORBIT and
 * re-attempts at rate alpha per orbiting job; a completion does NOT promote
 * from the orbit. That is why the queue length reported here is
 * `L_orbit + N_server` and not the marginal of an ordered buffer.
 *
 * THE TRUNCATION IS THE ANSWER'S ACCURACY, AND IT IS ITERATED HERE. The orbit
 * is unbounded, so `qsys_bmapphnn_retrial` truncates it; the reference does not
 * accept the first truncation but doubles the level until the relative tail
 * contribution `truncLevel * mass(top level) / L_orbit` falls under
 * `TailTolerance`. The ported engine implements ONE solve at a given level and
 * returns `truncLevel`, `topLevelMass` and `L_orbit`, which is exactly the
 * triple that refinement test needs, so the loop lives here rather than being
 * dropped. Without it a model is answered at the reference's FIXED default
 * level, whose formula depends on neither alpha nor gamma -- the two parameters
 * that actually govern the tail decay -- and the error is silent.
 *
 * WHERE THIS IS STRICTER THAN THE REFERENCE, and why each is the honest outcome:
 *
 *  - A NON-CONVERGED TRUNCATION IS AN ERROR, NOT A WARNING. When doubling hits
 *    the dimension cap before reaching the tolerance, `solver_mam_retrial.m`
 *    emits `line_warning` and returns the underestimate. The C++ port has no
 *    warning channel, so the same situation would return a number with nothing
 *    attached to it; it refuses instead, and names `orbit_maxlevel` as the way
 *    to ask for that number deliberately.
 *  - A NON-PHASE-TYPE SERVICE IS REFUSED, NOT APPROXIMATED. The reference's
 *    `warnIfApproximated` lets a Det, Replayer, Uniform, Gamma, Pareto, Weibull
 *    or Lognormal service through under an Erlang fit of matching mean and SCV,
 *    and says so in a warning. Everything but Det is already gated out by
 *    `runner_detail::check_processes`; Det reaches here, and `dist_to_map` would
 *    silently give it a 20-phase Erlang -- which is not only a different model
 *    but blows the per-level block up by C(n + 19, 19) states per busy server.
 *  - A BOUNDED RETRY COUNT IS REFUSED. `RetrialParam::max_attempts` and
 *    `DropStrategy::RETRIAL_WITH_LIMIT` say a job gives up after k failed
 *    retries. The generator has no attempt counter in its state descriptor, so
 *    there is no level at which that job could be told apart from a persistent
 *    one.
 *
 * WHAT THE C++ MODEL LAYER CANNOT SAY, so what is pinned at its reference
 * default rather than read: the orbit abandonment rate gamma (`sn.orbitImpatience`,
 * default 0), the batch rejection probability p (`sn.batchRejectProb`, default
 * 0), the finite orbit `sn.orbitMaxJobs`, the CONSTANT retrial policy
 * (`sn.retrialPolicy`; the ported engine is LINEAR-only, which is the
 * reference's own default), and the admission threshold R, which the reference
 * lowers below N-1 only from a Finite Capacity Region and `NetworkStruct` has no
 * region field. Each is stated here so that a later model-layer addition has a
 * list of what to wire, and none of them is invented from a capacity heuristic.
 *
 * REFERENCE DEFECTS in solver_mam_retrial.m and qsys_bmapphnn_retrial.m:
 *
 *  1. `Tolerance` IS PARSED AND NEVER USED. `solver_mam_retrial.m` reads
 *     `options.tol`, defaults it to 1e-10 and forwards it as 'Tolerance';
 *     `qsys_bmapphnn_retrial.m` assigns it on line 83 and no later line reads
 *     it. The single dense solve behind the generator is direct, so there is no
 *     iteration for it to terminate. Not propagated: `MamOptions::tol` is
 *     likewise not forwarded here, and the comment records why rather than
 *     leaving a caller to wonder why lowering the tolerance changes nothing.
 *  2. THE DEFAULT RETRIAL RATE IS A MAGIC 0.1. When no retrial process is
 *     attached, `alpha = 0.1` -- a rate with no relation to the model. It is
 *     reproduced, because a station carrying only a RETRIAL drop rule is a
 *     model the reference accepts and would otherwise be answered differently
 *     by the two codebases, but a caller reaching it has almost certainly
 *     forgotten `setRetrial`.
 *  3. `extractPHParams` NORMALIZES A BAD beta INSTEAD OF REPORTING IT. When the
 *     entry vector recovered from D1 does not sum to one it is rescaled, and
 *     when it sums to zero it is replaced by the uniform vector. Both hide a
 *     malformed service representation. Not reproduced: `map_pie` is the
 *     codebase's own embedded entry vector and is exact for any Markovian pair,
 *     and a pair that is not Markovian is refused above by name.
 *
 * ARITHMETIC. The generator is a rational expression in its inputs and the
 * stationary law is one linear solve, so nothing here needs a transcendental
 * function and this analyzer carries no `has_transcendental` gate -- unlike
 * `solver_mam_basic`, whose station solves run tolerance-terminated Riccati
 * iterations. `qsys_bmapphnn_retrial` documents the same property, and the
 * exact instantiation returns the stationary law of the truncated chain
 * exactly.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/qsys/qsys_bmapphnn_retrial.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mam/mam_types.h"
#include "line/solvers/mam/solver_mam_basic.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * The knobs `solver_mam_retrial.m` reads from `options.config` and from model
 * fields the C++ `NetworkStruct` does not carry.
 *
 * Kept as a separate argument rather than added to `MamOptions`, so that the
 * options every MAM analyzer shares do not grow a field only one of them reads.
 * The defaults are the reference's.
 */
struct MamRetrialConfig {
    /** `options.config.orbit_maxlevel`; 0 asks for the adaptive refinement. */
    std::size_t orbit_maxlevel = 0;
    /** `options.config.orbit_tailtol`, the relative orbit-truncation target. */
    double orbit_tailtol = 1e-6;
    /** `'MaxDim'`: the largest generator dimension the refinement may reach. */
    double max_dim = 2e5;
    /** `'MaxBlockSize'`: the largest per-level block V*d that may be assembled. */
    double max_block = 5e3;
};

/** What `qsys_is_retrial.m` returns: the station it found, or why it found none. */
struct MamRetrialInfo {
    bool ok = false;
    std::size_t station = 0;  ///< 1-based station index of the retrial queue
    std::size_t source = 0;   ///< 1-based station index of the Source
    std::size_t cls = 1;      ///< 1-based class index; the reference is single-class
    int nservers = 0;         ///< N
    long R = 0;               ///< admission threshold, N-1 unless an FCR lowers it
    std::string why;          ///< the reference's `errorMsg`, empty when ok
};

namespace retrial_detail {

/** Does any class at this station declare a retrial process? */
template <class T>
bool has_retrial_proc(const qn::NetworkStruct<T>& L, std::size_t ist) {
    const typename std::map<std::size_t, qn::RetrialParam<T> >::const_iterator it =
        L.retrialparam.find(ist);
    if (it == L.retrialparam.end()) return false;
    for (std::size_t r = 0; r < it->second.retrial_proc.size(); ++r)
        if (!it->second.retrial_proc[r].disabled) return true;
    return false;
}

/**
 * The reference's `assertMarkovian`, as a refusal naming the station and the
 * class. A process that is not a (D0, D1) pair as declared has no place in a
 * matrix-analytic generator, and the C++ layer would hand it over as an Erlang
 * of matching moments instead of reporting it.
 */
template <class T>
void assert_markovian(const qn::NetworkStruct<T>& L, std::size_t ist, std::size_t r,
                      const char* role) {
    const lang::Distrib<T>& d = L.service[ist - 1][r - 1];
    if (L.disabled[ist - 1][r - 1] || d.disabled)
        throw UnsupportedError(std::string("SolverMAM: the ") + role + " process of class '" +
                               L.classes[r - 1].name + "' at station '" + L.stations[ist - 1].name +
                               "' is disabled; the matrix-analytic retrial solver has no state for "
                               "a class that does not visit its own station");
    if (!basic_detail::is_markovian_type(d.type))
        throw UnsupportedError(
            std::string("SolverMAM: the ") + role + " process of class '" + L.classes[r - 1].name +
            "' at station '" + L.stations[ist - 1].name + "' is " + lang::process_to_text(d.type) +
            ", which has no Markovian (D0,D1) representation. The reference replaces it by an "
            "Erlang of matching mean and squared coefficient of variation and warns that the "
            "answer is approximate; this port refuses instead, because the substitute is a "
            "different model under the same method name. Use SolverCTMC, SolverSSA or SolverLDES");
    const Map<T> m = lang::dist_to_map(d);
    if (!basic_detail::is_markovian_map(m))
        throw UnsupportedError(std::string("SolverMAM: the ") + role + " process of class '" +
                               L.classes[r - 1].name + "' at station '" +
                               L.stations[ist - 1].name +
                               "' has a (D0,D1) pair that is not a Markovian generator (a negative "
                               "off-diagonal, a negative arrival entry, or rows that do not sum to "
                               "zero)");
}

/**
 * The per-level block width V*d of the generator, where d counts the multisets
 * of service phases over the busy servers.
 *
 * Computed BEFORE the first solve, because the engine assembles one dense
 * (truncLevel + 1) * V*d square matrix and an oversized block has to be
 * reported rather than attempted. In double throughout: the binomials overflow
 * an integer long before they reach a block size anyone would wait for.
 */
inline double retrial_block_width(std::size_t V, std::size_t M, int N) {
    double d = 0.0, term = 1.0;  // term = C(n + M - 1, M - 1), starting at n = 0
    for (int n = 0; n <= N; ++n) {
        if (n > 0) term *= static_cast<double>(n + static_cast<int>(M) - 1) / static_cast<double>(n);
        d += term;
        if (d > 1e12) return 1e12;  // already past any usable cap; stop growing
    }
    return static_cast<double>(V) * d;
}

/**
 * `orbitTruncationError`: the relative contribution the truncated tail would
 * add to the mean orbit length, which the mass sitting at the top level bounds.
 */
template <class T>
double orbit_trunc_error(const qsys::BmapPhNnRetrialResult<T>& r) {
    const double lorbit = num_traits<T>::to_double(r.L_orbit);
    const double top = num_traits<T>::to_double(r.topLevelMass);
    const double denom = lorbit > std::numeric_limits<double>::min()
                             ? lorbit
                             : std::numeric_limits<double>::min();
    return static_cast<double>(r.truncLevel) * top / denom;
}

}  // namespace retrial_detail

/**
 * Port of `qsys_is_retrial.m`, plus the reneging gate of `solver_mam_retrial.m`.
 *
 * Returns rather than throws, because `solver_mam_analyzer.m` uses it as the
 * 2c/2d branch predicate and only errors once every branch has declined. The
 * rejection reason travels in `why` so the dispatch can quote it.
 */
template <class T>
MamRetrialInfo mam_retrial_detect(const qn::NetworkStruct<T>& L) {
    using lang::DropStrategy;
    MamRetrialInfo info;

    // A mixed model has a closed chain whose population the orbit generator has
    // no level for, so the predicate must be "every class open", not "any".
    if (!L.is_open_model()) {
        info.why = "the BMAP/PH/N/N retrial solver requires an open queueing model";
        return info;
    }
    if (L.nclasses != 1) {
        info.why = "the BMAP/PH/N/N retrial solver supports a single class only";
        return info;
    }
    info.cls = 1;

    // A retrial station is one that DECLARES an orbit, or one that is
    // bufferless: capacity equal to the server count leaves no waiting line, so
    // a blocked job has nowhere to go but an orbit.
    std::vector<std::size_t> bufferless;
    for (std::size_t i = 1; i <= L.nstations; ++i) {
        if (L.stations[i - 1].nodetype != qn::NodeType::Queue) continue;
        const double cap = L.cap[i - 1];
        const bool isbufferless = std::isfinite(cap) && cap == L.stations[i - 1].nservers;
        if (retrial_detail::has_retrial_proc(L, i) || isbufferless) bufferless.push_back(i);
    }
    if (bufferless.empty()) {
        info.why = "no retrial queue found (declare one with setRetrial, or give a station a "
                   "capacity equal to its server count)";
        return info;
    }

    for (std::size_t i : bufferless) {
        for (std::size_t r = 0; r < L.nclasses; ++r)
            if (L.droprule[i - 1][r] == DropStrategy::RETRIAL ||
                L.droprule[i - 1][r] == DropStrategy::RETRIAL_WITH_LIMIT) {
                info.station = i;
                break;
            }
        if (info.station != 0) break;
    }
    if (info.station == 0) {
        info.why = "no RETRIAL drop strategy is configured on a bufferless queue";
        return info;
    }

    for (std::size_t i = 1; i <= L.nstations; ++i)
        if (L.stations[i - 1].nodetype == qn::NodeType::Source) {
            info.source = i;
            break;
        }
    if (info.source == 0) {
        info.why = "no Source node found";
        return info;
    }

    const double ns = L.stations[info.station - 1].nservers;
    if (!std::isfinite(ns) || ns < 1.0) {
        info.why = "the retrial station must have a finite, positive server count";
        return info;
    }
    info.nservers = static_cast<int>(ns);
    // The reference lowers R below N-1 only from a Finite Capacity Region, and
    // NetworkStruct has no region field, so N-1 is the only reachable value:
    // a retrial succeeds exactly when a server is free.
    info.R = static_cast<long>(info.nservers) - 1;
    info.ok = true;
    return info;
}

/**
 * The 'retrial' method's applicability as one sentence; empty when applicable.
 *
 * A "MUST BE PRESENT" RULE, which is why it cannot live in a feature set: a
 * FeatureSet says "I accept this construct", so it can refuse a model for HAVING
 * something and never for LACKING it. This analyzer needs the BMAP/PH/N/N
 * bufferless retrial topology to work on, and a model without one is not a
 * smaller retrial model, it is a different one.
 *
 * ONE PREDICATE, THREE CALLERS: `mam_dispatch`'s 'retrial' arm raises it,
 * `runner_detail::check_model_method` raises it ahead of the dispatch, and
 * `autosolver::auto_family_refusal` returns it through
 * `mam_model_method_refusal`, so the method the report offers and the method
 * that runs are the same set.
 */
template <class T>
std::string mam_retrial_refusal(const qn::NetworkStruct<T>& L) {
    const MamRetrialInfo info = mam_retrial_detect(L);
    if (info.ok) return std::string();
    return "SolverMAM: the 'retrial' method needs a BMAP/PH/N/N bufferless retrial topology; "
           "this model declares no orbit (" + info.why + "). Use method 'default'";
}

/**
 * Port of `solver_mam_retrial.m`.
 *
 * @param L   the refreshed struct, with non-Markovian processes already gated
 *            out by the runner
 * @param opt the MAM options; none of them reaches the engine, see defect 1
 * @param cfg the orbit truncation controls
 */
template <class T>
mva::MvaSolution<T> solver_mam_retrial(const qn::NetworkStruct<T>& L, const MamOptions& opt,
                                       const MamRetrialConfig& cfg = MamRetrialConfig()) {
    using lang::DropStrategy;
    (void)opt;  // see defect 1: the reference forwards options.tol and nothing reads it
    const T zero = num_traits<T>::from_int(0);

    const MamRetrialInfo info = mam_retrial_detect(L);
    if (!info.ok)
        throw UnsupportedError(
            "SolverMAM: no valid impatience configuration was detected (" + info.why +
            "). The reneging route of solver_mam_retrial.m, the MAP/M/s+G call-centre model of "
            "Gursoy, Mehr and Akar solved through MRMFQSolver, is not reachable from C++ at all: "
            "it is gated on sn.impatienceClass, sn.patienceProc and ImpatienceType.RENEGING, and "
            "the C++ NetworkStruct carries none of the three");

    const std::size_t q = info.station, src = info.source, r = info.cls;

    retrial_detail::assert_markovian(L, src, r, "arrival");
    retrial_detail::assert_markovian(L, q, r, "service");

    // A bounded retry count is a per-job attribute; the orbit level counts jobs
    // and nothing else, so a job on its last attempt is indistinguishable there.
    if (L.droprule[q - 1][r - 1] == DropStrategy::RETRIAL_WITH_LIMIT)
        throw UnsupportedError(
            "SolverMAM: station '" + L.stations[q - 1].name +
            "' declares DropStrategy.RETRIAL_WITH_LIMIT. The BMAP/PH/N/N generator indexes its "
            "levels by the orbit population alone and has no per-job attempt counter, so a job "
            "about to give up cannot be told from a persistent one. Use SolverCTMC or SolverSSA");

    // ---- the arrival BMAP ---------------------------------------------------
    // A genuine BATCH arrival would give D = {D0, D1, ..., DK}; lang::Distrib
    // holds a single (D0, D1) pair, so every model the C++ layer can build is
    // the K = 1 case. The engine takes the general list, so nothing is lost
    // here beyond the batch models the model layer cannot express.
    const Map<T> arv = lang::dist_to_map(L.service[src - 1][r - 1]);
    std::vector<Matrix<T> > D;
    D.push_back(arv.D0);
    D.push_back(arv.D1);

    // ---- the PH service -----------------------------------------------------
    const Map<T> svc = lang::dist_to_map(L.service[q - 1][r - 1]);
    const Matrix<T> S = svc.D0;
    const std::vector<T> beta = map_pie(svc);

    // ---- the orbit parameters ----------------------------------------------
    // The retrial process is an inter-retry TIME, and the generator carries a
    // single scalar rate per orbiting job, so only a one-phase process can be
    // read. The reference takes -D0(1,1) whatever the order is, which for a
    // multi-phase process is the exit rate of its FIRST phase and not the mean
    // retry rate of anything.
    T alpha = num_traits<T>::from_double(0.1);  // see defect 2
    const typename std::map<std::size_t, qn::RetrialParam<T> >::const_iterator rit =
        L.retrialparam.find(q);
    if (rit != L.retrialparam.end()) {
        const qn::RetrialParam<T>& rp = rit->second;
        if (rp.max_attempts.size() >= r && rp.max_attempts[r - 1] > 0)
            throw UnsupportedError(
                "SolverMAM: station '" + L.stations[q - 1].name + "' bounds the retrial count at " +
                std::to_string(rp.max_attempts[r - 1]) +
                " attempts. The BMAP/PH/N/N generator has no per-job attempt counter in its state "
                "descriptor. Use SolverCTMC or SolverSSA");
        if (rp.retrial_proc.size() >= r && !rp.retrial_proc[r - 1].disabled) {
            const Map<T> ret = lang::dist_to_map(rp.retrial_proc[r - 1]);
            if (ret.order() != 1)
                throw UnsupportedError(
                    "SolverMAM: the retrial process of class '" + L.classes[r - 1].name +
                    "' at station '" + L.stations[q - 1].name + "' has " +
                    std::to_string(ret.order()) +
                    " phases. The generator retries each orbiting job at a single scalar rate, so "
                    "only an exponential inter-retry time can be represented; the reference reads "
                    "-D0(1,1) regardless, which is the exit rate of the first phase and not the "
                    "retry rate of the process");
            alpha = T(zero - ret.D0(0, 0));
        } else if (rp.retrial_rate.size() >= r && rp.retrial_rate[r - 1] > zero) {
            alpha = rp.retrial_rate[r - 1];
        }
    }
    if (!(alpha > zero))
        throw InputError("SolverMAM: the retrial rate at station '" + L.stations[q - 1].name +
                         "' must be positive; a zero rate means an orbiting job never retries and "
                         "the orbit is not positive recurrent at any load");

    const T gamma = zero;  // sn.orbitImpatience, absent from NetworkStruct
    const T p = zero;      // sn.batchRejectProb, absent from NetworkStruct

    // ---- the truncation, and the refinement the engine leaves to its caller --
    const double blockw =
        retrial_detail::retrial_block_width(arv.order(), svc.order(), info.nservers);
    if (blockw > cfg.max_block)
        throw UnsupportedError(
            "SolverMAM: the per-level block of the retrial generator at station '" +
            L.stations[q - 1].name + "' is " + std::to_string(static_cast<long long>(blockw)) +
            " states wide, past the " + std::to_string(static_cast<long long>(cfg.max_block)) +
            "-state cap: the service has " + std::to_string(svc.order()) +
            " phases and the station " + std::to_string(info.nservers) +
            " servers, and the busy servers' phase multisets grow as C(n + phases - 1, phases - 1). "
            "Fit the service with fewer phases, reduce the server count, or raise "
            "MamRetrialConfig::max_block if the memory cost is acceptable");

    const std::vector<long> Rv(1, info.R);
    qsys::BmapPhNnRetrialOptions qopt;
    qopt.maxLevel = cfg.orbit_maxlevel;
    qsys::BmapPhNnRetrialResult<T> res =
        qsys::qsys_bmapphnn_retrial(D, beta, S, info.nservers, alpha, gamma, p, Rv, qopt);

    if (cfg.orbit_maxlevel == 0) {
        // The reference doubles the level until the tail contribution is under
        // TailTolerance. `blockw` is the per-level width, so (level + 1) * blockw
        // is the dimension the next solve would assemble.
        while (retrial_detail::orbit_trunc_error(res) > cfg.orbit_tailtol) {
            const std::size_t next = 2 * res.truncLevel;
            if (static_cast<double>(next + 1) * blockw > cfg.max_dim)
                throw NumericError(
                    "SolverMAM: the orbit truncation of station '" + L.stations[q - 1].name +
                    "' did not reach the requested accuracy: the tail still contributes " +
                    std::to_string(retrial_detail::orbit_trunc_error(res)) + " of the mean orbit "
                    "length at level " + std::to_string(res.truncLevel) +
                    ", against a target of " + std::to_string(cfg.orbit_tailtol) +
                    ", and doubling would pass the dimension cap. The orbit measures would be "
                    "underestimated, which the reference reports as a warning and this port has no "
                    "channel for. Raise MamRetrialConfig::max_dim, or set "
                    "MamRetrialConfig::orbit_maxlevel to accept a fixed truncation");
            qopt.maxLevel = next;
            res = qsys::qsys_bmapphnn_retrial(D, beta, S, info.nservers, alpha, gamma, p, Rv, qopt);
        }
    }

    // ---- the LINE metric tuple ---------------------------------------------
    mva::MvaSolution<T> out;
    const std::size_t M = L.nstations, K = L.nclasses;
    out.Q = Matrix<T>(M, K, zero);
    out.U = Matrix<T>(M, K, zero);
    out.R = Matrix<T>(M, K, zero);
    out.Tp = Matrix<T>(M, K, zero);
    out.C.assign(K, zero);
    out.X.assign(K, zero);

    // The orbit is part of the station: a job waiting to retry has not left it.
    out.Q(q - 1, r - 1) = T(res.L_orbit + res.N_server);
    out.U(q - 1, r - 1) = res.utilization;
    out.Tp(q - 1, r - 1) = res.throughput;
    if (res.throughput > zero)
        out.R(q - 1, r - 1) = T(out.Q(q - 1, r - 1) / res.throughput);
    // else the reference writes Inf, and the NaN sweep at the end of
    // mam_dispatch zeroes every non-finite entry before any caller sees it; a
    // zero written here is the same value by a route that also holds at
    // T = Rational, which has no infinity to write.

    out.X[r - 1] = res.throughput;
    out.C[r - 1] = out.R(q - 1, r - 1);
    // The reference's `totiter` for this analyzer is the truncation level that
    // produced the answer, not an iteration count.
    out.iter = static_cast<int>(res.truncLevel);
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_SOLVERS_MAM_SOLVER_MAM_RETRIAL_H
