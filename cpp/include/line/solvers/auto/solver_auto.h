/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_AUTO_SOLVER_AUTO_H
#define LINE_SOLVERS_AUTO_SOLVER_AUTO_H

/**
 * The SolverAUTO chooser: which solver a model is handed to.
 *
 * WHAT THIS PORTS. `matlab/src/solvers/AUTO/@@SolverAUTO/` selects in two
 * layers, and both are here:
 *
 *   solverTraits.m computes the structural traits ONCE per call, so that the
 *     choosers stay tables of rankings rather than a second place where model
 *     inspection is written. `auto_traits` is that function.
 *
 *   chooseSolverRanked.m walks a list of candidate slots and returns the FIRST
 *     that exists and whose feature set accepts the model, or nothing, so the
 *     caller can fall back rather than hand an infeasible solver to the
 *     delegate. `auto_ranked` is that function, and it is what makes every
 *     ranking below a preference rather than a claim.
 *
 *   chooseSolverHeur.m / chooseAvgSolverHeur.m / chooseSolverExact.m /
 *     chooseSolverSim.m and the `fast` and `accurate` arms of chooseSolver.m
 *     are the rankings themselves, ported list for list, in the reference's own
 *     order. The Network, LayeredNetwork and Environment arms all have a
 *     counterpart here.
 *
 * THE FEATURE-SET GATE IS THE WHOLE POINT OF THE REWRITE. The previous port
 * transcribed the SUPERSEDED feature cascade (population thresholds 30/10/5
 * over a first-match tree) and consulted no feature set at all, so it named a
 * solver that could then refuse the model. Here `auto_supports` asks the same
 * question `Solver.supports` asks -- is every feature the model uses declared
 * by that solver? -- from `used_lang_features` and the per-solver sets in
 * `solver_feature_sets.h`, and the answer decides.
 *
 * TWO GATES ARE FINER THAN A FLAT FEATURE SET, and the reference states both
 * outside `getFeatureSet` for that reason:
 *
 *   'exact' needs a product-form solution (SolverMVA.supportsExactness,
 *     SolverNC.supportsModelMethod), with the order-independent and
 *     pass-and-swap stations exempt for MVA because solver_mva_oi is exact for
 *     them regardless. Reproduced in `auto_supports`.
 *
 *   CTMC needs its chain to FIT. The reference screens the slot with
 *     SolverCTMC.isStateSpaceTractable, which prices the worst-case state space
 *     against host memory through a profiled power law. This port has no such
 *     calibration and its generator has a hard cap instead, so the screen here
 *     compares the SAME estimator -- `ctmc_state_space_logsize`, ported
 *     factor for factor -- against that cap. The gate is therefore
 *     host-independent where the reference's is host-dependent: the estimate is
 *     identical, the budget it is compared to is this port's own.
 *
 * AN ABSENT ENGINE IS SKIPPED, AND SAID SO. In the reference an unavailable
 * candidate is an empty slot and chooseSolverRanked skips it silently; that is
 * how SolverLQNS behaves when the binary is not installed. The MAM layer engine
 * is absent from this port in exactly that sense, and LDES and LQNS are absent
 * only where their engines are not installed, so all three are skipped the same
 * way -- but skipping changes which engine answers, so every choice carries
 * `skipped`, the slots that outranked the winner and had no engine behind them.
 * The CLI prints it. Silence there would report the second choice as if it had
 * been the first.
 *
 * JMT IS NOT A CANDIDATE AT ALL, in the reference either (SolverAUTO.m:45-47):
 * LDES subsumes its feature set, so automatic selection never dispatches to the
 * external simulator and SolverJMT stays reachable only through an explicit
 * token. The enum has no JMT slot for that reason.
 *
 * THE HOMOGENEOUS-SCHEDULING PREDICATE COLLAPSES. `has_homogeneous_scheduling`
 * reproduces the reference's findstring defect and degenerates to
 * `nstations == 1` for every discipline (see NetworkStruct). Only the
 * response-time-CDF branch and one avgOrder arm consult it now, and both are
 * written as the reference writes them.
 *
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/solver_feature_sets.h"
#include "line/solvers/wrappers/ldes/ldes_probe.h"
#include "line/solvers/wrappers/lqns/lqns_probe.h"
#include "line/util/error.h"

namespace line {
namespace autosolver {

/**
 * The Network candidate slots, in the reference's slot order (SolverAUTO.m:41-50),
 * which is also the order the delegate retries in. There is no JMT slot: the
 * reference removed it as a candidate.
 */
enum class AutoSolver { MVA = 0, NC, MAM, FLUID, SSA, CTMC, LDES };

/** The LayeredNetwork candidate slots (SolverAUTO.m:52-56). */
enum class AutoLayered { LQNS = 0, LN_NC, LN_MVA, LN_MAM, LN_FLUID };

/** The Environment candidate slots (SolverAUTO.m:58-60). */
enum class AutoEnv { ENV_MVA = 0, ENV_NC, ENV_FLUID };

/** The selection intents of `SolverAUTO.selectionIntents`, less 'bound'. */
enum class AutoMode { HEUR, EXACT, SIM, FAST, ACCURATE };

/** Population at or below which an exact solver is preferred (EXACT_POPULATION_MAX). */
const double kAutoExactPopulationMax = 5.0;

inline const char* auto_solver_name(AutoSolver s) {
    switch (s) {
        case AutoSolver::MVA: return "mva";
        case AutoSolver::NC: return "nc";
        case AutoSolver::MAM: return "mam";
        case AutoSolver::FLUID: return "fld";
        case AutoSolver::SSA: return "ssa";
        case AutoSolver::CTMC: return "ctmc";
        case AutoSolver::LDES: return "ldes";
    }
    return "";
}

/**
 * The layered names are the CLI's own tokens, because that is what the choice
 * is spent on: `ln.comom` runs the layers under NC and `ln.mva` under MVA.
 */
inline const char* auto_layered_name(AutoLayered s) {
    switch (s) {
        case AutoLayered::LQNS: return "lqns";
        case AutoLayered::LN_NC: return "ln.comom";
        case AutoLayered::LN_MVA: return "ln.mva";
        case AutoLayered::LN_MAM: return "ln.mam";
        case AutoLayered::LN_FLUID: return "ln.fluid";
    }
    return "";
}

inline const char* auto_env_name(AutoEnv s) {
    switch (s) {
        case AutoEnv::ENV_MVA: return "env.mva";
        case AutoEnv::ENV_NC: return "env.nc";
        case AutoEnv::ENV_FLUID: return "env.fluid";
    }
    return "";
}

/**
 * True when this port has an engine behind the slot at all.
 *
 * THE LDES SLOT IS MACHINE-DEPENDENT, like the LQNS one below and for the same
 * reason: the simulator is a client of an engine LINE ships beside the binary
 * (`common/ldes`, `common/ldes.jar`), and `ldes_is_available` answers whether
 * this machine has one. A port that answered "never" would silently disagree
 * with the reference wherever the engine IS present -- LDES leads the
 * loss-metric, aggregate-sampling and `sim` rankings, so the disagreement would
 * be a different engine answering, not a missing option.
 */
inline bool auto_solver_is_available(AutoSolver s) {
    if (s == AutoSolver::LDES) return ldes::ldes_is_available();
    return true;
}

/**
 * The LN layer engines are the four `--layer-solver` takes (mva, nc, fluid,
 * ssa), so there is no MAM-layer engine to select.
 *
 * LQNS IS AVAILABLE ONLY WHERE ITS BINARY IS. That makes this choice
 * machine-dependent, and deliberately so: `chooseAvgSolverHeur.m` gates its
 * LQNS candidate on `SolverLQNS.isAvailable()` for the same reason, because
 * LINE ships no LQNS binary. A port that answered "never" here would silently
 * disagree with the reference on every machine that has one installed.
 */
inline bool auto_layered_is_available(AutoLayered s) {
    if (s == AutoLayered::LQNS) return lqns::lqns_is_available();
    return s != AutoLayered::LN_MAM;
}

/**
 * SolverENV solves a stage with the fluid analyzer under every coupling but the
 * state-vector one, which uniformizes a CTMC; neither an MVA nor an NC stage
 * solver exists here, so those two slots have no engine.
 */
inline bool auto_env_is_available(AutoEnv s) { return s == AutoEnv::ENV_FLUID; }

// ---------------------------------------------------------------------------
// Method method names: a selection intent, or a method family
// ---------------------------------------------------------------------------

/**
 * `SolverAUTO.selectionIntents`, less 'bound'.
 *
 * 'bound' IS accepted, but it is not a ranking mode: the reference sets
 * selectionMode='bound' and options.method='auto', i.e. it names SolverBA
 * outright. `auto_resolve_token` therefore rewrites it to the 'ba' family
 * before this test runs, which is why it is absent here.
 */
inline bool auto_is_selection_intent(const std::string& token) {
    return token.empty() || token == "default" || token == "auto" || token == "heur" ||
           token == "sim" || token == "exact" || token == "fast" ||
           token == "accurate";
}

inline AutoMode auto_mode_of_token(const std::string& token) {
    if (token.empty() || token == "default" || token == "auto" || token == "heur")
        return AutoMode::HEUR;
    if (token == "exact") return AutoMode::EXACT;
    if (token == "sim") return AutoMode::SIM;
    if (token == "fast") return AutoMode::FAST;
    if (token == "accurate") return AutoMode::ACCURATE;
    throw InputError("SolverAUTO: '" + token + "' is not a selection intent");
}

/**
 * `SolverAUTO.familyAlias`: the canonical family of a method name, or "" when the
 * method name names none. A family method name is not a selection intent -- it asks for a
 * named engine and bypasses the ranking, as `-s auto --method nc.comom` does.
 */
inline std::string auto_family_alias(const std::string& name) {
    if (name == "mam") return "mam";
    if (name == "ag") return "ag";
    if (name == "mva") return "mva";
    if (name == "nc") return "nc";
    if (name == "fluid" || name == "fld") return "fluid";
    if (name == "jmt") return "jmt";
    if (name == "ssa") return "ssa";
    if (name == "ctmc") return "ctmc";
    if (name == "ldes" || name == "des") return "ldes";
    if (name == "ba") return "ba";
    if (name == "env") return "env";
    if (name == "ln") return "ln";
    if (name == "lqns" || name == "lqsim") return "lqns";
    if (name == "qns") return "qns";
    if (name == "uq") return "uq";
    return "";
}

/**
 * `resolveMethodToken`, minus the unqualified-algorithm-name arm.
 *
 * A method name is an intent, or `family[.submethod]`. The reference has a third form
 * -- a bare algorithm name such as `comom`, resolved by asking every family for
 * its `listValidMethods` -- which needs a per-solver method-name registry this
 * port does not have; it is refused by name here, with the qualified spelling
 * in the message, rather than guessed at.
 */
struct AutoToken {
    bool is_intent = true;
    AutoMode mode = AutoMode::HEUR;
    std::string family;     ///< empty when is_intent
    std::string submethod;  ///< the method handed to the family, "default" when bare
};

inline AutoToken auto_resolve_token(const std::string& raw) {
    AutoToken t;
    const std::string token = raw.empty() ? std::string("default") : raw;
    // 'bound' is a selection intent in the reference's own list, and it selects
    // SolverBA with method 'auto' rather than picking a ranking
    // (SolverAUTO.m:116-121). Rewriting it to the family method name here reuses the
    // family path below and keeps the intent spelling usable, which is what the
    // other three codebases accept.
    if (token == "bound") {
        t.is_intent = false;
        t.family = "ba";
        t.submethod = "auto";
        return t;
    }
    if (auto_is_selection_intent(token)) {
        t.is_intent = true;
        t.mode = auto_mode_of_token(token);
        return t;
    }
    const std::size_t dot = token.find('.');
    const std::string head = dot == std::string::npos ? token : token.substr(0, dot);
    const std::string rest = dot == std::string::npos ? std::string() : token.substr(dot + 1);
    const std::string fam = auto_family_alias(head);
    if (fam.empty())
        throw InputError(
            "SolverAUTO: '" + token +
            "' is neither a selection intent (default, heur, exact, sim, fast, accurate, "
            "bound) "
            "nor a method family (mva, nc, ctmc, fluid, mam, ag, ba, ssa, ldes, jmt, qns, ln, "
            "env, lqns, uq). A bare algorithm name must be qualified by its family here, as in "
            "'nc.comom': resolving it needs the per-family method registry "
            "(listValidMethods) that this port does not carry");
    t.is_intent = false;
    t.family = fam;
    t.submethod = rest.empty() ? std::string("default") : rest;
    return t;
}

// ---------------------------------------------------------------------------
// solverTraits.m
// ---------------------------------------------------------------------------

/** The structural traits the rankings are keyed on. Port of `solverTraits.m`. */
struct AutoTraits {
    bool has_cache = false;
    bool has_fcr = false;
    bool has_fork = false;
    bool has_map = false;
    /** "none", "preempt", "ps" or "hol", in the reference's own precedence. */
    std::string prio = "none";
    bool is_closed = false;
    bool is_open = false;
    bool is_mixed = false;
    bool has_multi_server = false;
    bool is_product_form = false;
    double pop_per_chain = 0.0;
    double total_jobs = 0.0;
    bool single_chain = false;
};

template <class T>
AutoTraits auto_traits(const qn::NetworkStruct<T>& sn) {
    using qn::NodeType;
    using lang::ProcessType;
    using lang::SchedStrategy;

    AutoTraits t;
    for (const qn::NodeDef& nd : sn.nodes)
        if (nd.nodetype == NodeType::Cache) t.has_cache = true;
    t.has_fcr = !sn.regions.empty();
    t.has_fork = sn.has_fork();
    t.has_multi_server = sn.has_multi_server();
    t.is_product_form = sn.has_product_form();
    t.single_chain = sn.nchains == 1;

    bool has_open = false, has_closed = false;
    for (const qn::JobClass& c : sn.classes) {
        if (std::isinf(c.population)) has_open = true;
        else has_closed = true;
    }
    t.is_open = has_open && !has_closed;
    t.is_closed = has_closed && !has_open;
    t.is_mixed = has_open && has_closed;

    t.total_jobs = sn.total_jobs();
    if (sn.nchains > 0) t.pop_per_chain = t.total_jobs / static_cast<double>(sn.nchains);

    // Autocorrelated arrival or service: only MAM keeps the correlation, every
    // other analytical solver sees the marginal only.
    for (std::size_t i = 1; i <= sn.nstations && !t.has_map; ++i)
        for (std::size_t r = 1; r <= sn.nclasses; ++r) {
            const ProcessType p = sn.procid(i, r);
            if (p == ProcessType::MAP || p == ProcessType::MMPP2) {
                t.has_map = true;
                break;
            }
        }

    // Preemptive priority is a strictly narrower capability than HOL, so the
    // two rank differently and are distinguished here rather than downstream.
    bool preempt = false, psprio = false, hol = false;
    for (const qn::Station<T>& s : sn.stations) {
        if (s.sched == SchedStrategy::FCFSPRPRIO || s.sched == SchedStrategy::FCFSPIPRIO ||
            s.sched == SchedStrategy::LCFSPRPRIO || s.sched == SchedStrategy::LCFSPIPRIO)
            preempt = true;
        if (s.sched == SchedStrategy::PSPRIO || s.sched == SchedStrategy::DPSPRIO ||
            s.sched == SchedStrategy::GPSPRIO)
            psprio = true;
        if (s.sched == SchedStrategy::HOL || s.sched == SchedStrategy::LCFSPRIO ||
            s.sched == SchedStrategy::SRPTPRIO)
            hol = true;
    }
    if (preempt) t.prio = "preempt";
    else if (psprio) t.prio = "ps";
    else if (hol) t.prio = "hol";
    return t;
}

// ---------------------------------------------------------------------------
// The CTMC screen: ctmc_state_space_logsize.m
// ---------------------------------------------------------------------------

/**
 * Worst-case log-size of the CTMC state space induced by `sn`, summed in log
 * space over the reference's four factors: job placements (stars and bars, per
 * class, open classes truncated at the cutoff, over the stations that keep no
 * ordered buffer), the class-sequence multiplicity of every order-preserving
 * buffer, service phases, and one routing pointer per round-robin (node,
 * class).
 *
 * `cutoff` negative selects the analyzer's own default for open and mixed
 * models, `ceil(6000^(1/(M*K)))`.
 *
 * MATLAB reads `sn.phasessz`; this port reads `phases_of`, which is the same
 * quantity floored at one -- a single-phase representation contributes no
 * factor either way.
 */
template <class T>
double auto_ctmc_state_space_logsize(const qn::NetworkStruct<T>& sn, double cutoff = -1.0) {
    using lang::RoutingStrategy;
    using lang::SchedStrategy;

    const std::size_t M = sn.nstations, K = sn.nclasses;
    if (M == 0 || K == 0) return 0.0;
    if (!(cutoff > 0.0))
        cutoff = std::ceil(std::pow(6000.0, 1.0 / static_cast<double>(M * K)));

    const auto is_share = [](SchedStrategy s) {
        return s == SchedStrategy::INF || s == SchedStrategy::PS || s == SchedStrategy::DPS ||
               s == SchedStrategy::GPS || s == SchedStrategy::PSPRIO ||
               s == SchedStrategy::DPSPRIO || s == SchedStrategy::GPSPRIO ||
               s == SchedStrategy::LPS;
    };
    std::vector<bool> is_buffered(K, true);
    std::size_t Kb = 0;
    for (std::size_t r = 0; r < K; ++r) {
        if (r < sn.issignal.size() && sn.issignal[r]) is_buffered[r] = false;
        if (is_buffered[r]) ++Kb;
    }
    std::size_t n_ord = 0;
    if (Kb > 1) {
        for (std::size_t i = 0; i < M; ++i) {
            const SchedStrategy s = sn.stations[i].sched;
            if (s == SchedStrategy::EXT || is_share(s)) continue;
            ++n_ord;
        }
    }
    const std::size_t m_place = (M > n_ord) ? (M - n_ord) : 0;

    double log_n = 0.0;
    std::vector<double> nk_eff(K, 0.0);
    for (std::size_t r = 0; r < K; ++r) {
        const double nk = std::isinf(sn.classes[r].population) ? cutoff : sn.classes[r].population;
        nk_eff[r] = nk;
        if (m_place >= 1) {
            const double m = static_cast<double>(m_place);
            log_n += std::lgamma(1.0 + nk + m - 1.0) - std::lgamma(m) - std::lgamma(1.0 + nk);
        }
    }

    if (n_ord > 0) {
        double tot_jobs = 0.0;
        for (std::size_t r = 0; r < K; ++r)
            if (is_buffered[r]) tot_jobs += nk_eff[r];
        const double log_k = std::log(static_cast<double>(Kb));
        const double log_seq = (tot_jobs + 1.0) * log_k - std::log(static_cast<double>(Kb) - 1.0) +
                               std::log1p(-std::exp(-(tot_jobs + 1.0) * log_k));
        log_n += static_cast<double>(n_ord) * log_seq;
    }

    for (std::size_t i = 1; i <= M; ++i) {
        const SchedStrategy sched = sn.stations[i - 1].sched;
        const bool shares = sched == SchedStrategy::INF || sched == SchedStrategy::PS ||
                            sched == SchedStrategy::DPS || sched == SchedStrategy::GPS ||
                            sched == SchedStrategy::PSPRIO || sched == SchedStrategy::DPSPRIO ||
                            sched == SchedStrategy::GPSPRIO || sched == SchedStrategy::LPS;
        for (std::size_t r = 1; r <= K; ++r) {
            const double p = static_cast<double>(std::max<std::size_t>(sn.phases_of(i, r), 1));
            if (p <= 1.0) continue;
            double m;
            if (sched == SchedStrategy::EXT) m = 1.0;
            else if (shares) m = nk_eff[r - 1];
            else m = std::min(nk_eff[r - 1], sn.stations[i - 1].nservers);
            if (!std::isfinite(m)) m = nk_eff[r - 1];
            log_n += std::lgamma(1.0 + m + p - 1.0) - std::lgamma(p) - std::lgamma(1.0 + m);
        }
    }

    // The routing pointers. MATLAB counts the out-degree from `sn.connmatrix`,
    // which this port does not carry; `rtnodes` is the same graph after the
    // refresh has resolved the strategies, and is what every other consumer
    // here reads for the same purpose (see NetworkStruct::downstream_stations).
    const std::size_t N = sn.nodes.size(), R = sn.nclasses;
    if (sn.rtnodes.rows() >= N * R) {
        for (std::size_t i = 1; i <= N; ++i) {
            std::size_t nout = 0;
            for (std::size_t j = 1; j <= N; ++j) {
                bool linked = false;
                for (std::size_t r = 0; r < R && !linked; ++r)
                    for (std::size_t s = 0; s < R && !linked; ++s)
                        if (num_traits<T>::to_double(
                                sn.rtnodes((i - 1) * R + r, (j - 1) * R + s)) > 0.0)
                            linked = true;
                if (linked) ++nout;
            }
            if (nout <= 1) continue;
            std::size_t nrr = 0;
            const std::vector<RoutingStrategy>& rt = sn.nodes[i - 1].routing;
            for (std::size_t r = 0; r < rt.size(); ++r)
                if (rt[r] == RoutingStrategy::RROBIN || rt[r] == RoutingStrategy::WRROBIN) ++nrr;
            if (nrr > 0) log_n += static_cast<double>(nrr) * std::log(static_cast<double>(nout));
        }
    }
    return log_n;
}

/**
 * The cap `reachable_space_generator` enforces (solver_ctmc.h, `maxst`). It is
 * this port's budget, standing in for the reference's memory model.
 */
const double kAutoCtmcStateCap = 3000000.0;

template <class T>
bool auto_ctmc_is_tractable(const qn::NetworkStruct<T>& sn, double cutoff = -1.0) {
    return auto_ctmc_state_space_logsize(sn, cutoff) <= std::log(kAutoCtmcStateCap);
}

// ---------------------------------------------------------------------------
// chooseSolverRanked.m
// ---------------------------------------------------------------------------

/**
 * `Solver.supports(model)` for a candidate slot, tightened by the two
 * method-level rules a flat feature set cannot express.
 *
 * `method name` is the method the slot would run: "" or "default" for the solver's
 * own default, "exact" for the exactness-gated request.
 */
template <class T>
bool auto_supports(AutoSolver s, const qn::NetworkStruct<T>& sn, const std::string& token) {
    using lang::SchedStrategy;
    if (!auto_solver_is_available(s)) return false;
    const std::string method = token.empty() ? std::string("default") : token;

    qn::FeatureSet declared;
    switch (s) {
        case AutoSolver::MVA: declared = qn::mva_feature_set(method); break;
        case AutoSolver::NC: declared = qn::nc_feature_set(method); break;
        case AutoSolver::MAM: declared = qn::mam_feature_set(method); break;
        case AutoSolver::FLUID: declared = qn::fluid_feature_set(method); break;
        case AutoSolver::SSA: declared = qn::ssa_feature_set(method); break;
        case AutoSolver::CTMC: declared = qn::ctmc_feature_set(method); break;
        // The LDES slot is gated like any other now that an engine can stand
        // behind it: `auto_solver_is_available` has already answered whether
        // this machine has one, and what remains is the reference's own
        // declaration, which the client honours by forwarding the model to the
        // engine that implements it.
        case AutoSolver::LDES: declared = qn::ldes_feature_set(method); break;
    }
    if (!qn::feature_set_supports(auto_solver_name(s), declared, qn::used_lang_features(sn)).ok)
        return false;

    if (method == "exact" && (s == AutoSolver::MVA || s == AutoSolver::NC)) {
        if (!sn.has_product_form()) {
            // solver_mva_oi is exact for the order-independent and
            // pass-and-swap stations whatever the product-form test says; NC
            // has no such exemption.
            bool oi = false;
            if (s == AutoSolver::MVA)
                for (const qn::Station<T>& st : sn.stations)
                    if (st.sched == SchedStrategy::OI || st.sched == SchedStrategy::PAS) oi = true;
            if (!oi) return false;
        }
    }
    if (s == AutoSolver::CTMC && !auto_ctmc_is_tractable(sn)) return false;
    return true;
}

/** What a ranking resolved to, and what it had to skip to get there. */
struct AutoChoice {
    AutoSolver solver = AutoSolver::MVA;
    /** Slots that outranked `solver` and have no engine in this port. */
    std::vector<AutoSolver> skipped;
    /** The method the choice was gated on: "" for the default, "exact". */
    std::string method;
};

/**
 * `chooseSolverRanked`: the first slot in ORDER that exists and accepts the
 * model. Returns false when none qualifies, so the caller can fall back.
 */
template <class T>
bool auto_ranked(const std::vector<AutoSolver>& order, const qn::NetworkStruct<T>& sn,
                 const std::string& token, AutoChoice& out) {
    std::vector<AutoSolver> skipped;
    for (std::size_t k = 0; k < order.size(); ++k) {
        if (!auto_solver_is_available(order[k])) {
            skipped.push_back(order[k]);
            continue;
        }
        if (!auto_supports(order[k], sn, token)) continue;
        out.solver = order[k];
        out.skipped = skipped;
        out.method = (token == "default") ? std::string() : token;
        return true;
    }
    return false;
}

/**
 * The candidate pool in slot order, filtered by `supports`: what the delegate
 * retries through after the chosen solver fails (SolverAUTO.m:139-147).
 */
template <class T>
std::vector<AutoSolver> auto_candidates(const qn::NetworkStruct<T>& sn) {
    static const AutoSolver kSlots[] = {AutoSolver::MVA,  AutoSolver::NC,   AutoSolver::MAM,
                                        AutoSolver::FLUID, AutoSolver::SSA, AutoSolver::CTMC,
                                        AutoSolver::LDES};
    std::vector<AutoSolver> out;
    for (std::size_t i = 0; i < sizeof(kSlots) / sizeof(*kSlots); ++i)
        if (auto_solver_is_available(kSlots[i]) && auto_supports(kSlots[i], sn, "default"))
            out.push_back(kSlots[i]);
    return out;
}

// ---------------------------------------------------------------------------
// chooseAvgSolverHeur.m, the Network arm
// ---------------------------------------------------------------------------

namespace detail {

/**
 * `avgOrder` in chooseAvgSolverHeur.m, arm for arm and in its order. The arms
 * OVERLAP, so the first match wins and reordering silently rehomes models.
 *
 * `homogeneous_inf` is the penultimate arm's predicate, passed in because it
 * needs the struct and this table does not otherwise.
 */
inline std::vector<AutoSolver> avg_order(const AutoTraits& t, bool homogeneous_inf) {
    if (t.has_cache)
        return {AutoSolver::NC, AutoSolver::MVA, AutoSolver::FLUID, AutoSolver::CTMC,
                AutoSolver::LDES};
    if (t.has_fcr) {
        if (t.total_jobs <= 10.0) return {AutoSolver::NC, AutoSolver::CTMC, AutoSolver::LDES};
        return {AutoSolver::NC, AutoSolver::LDES, AutoSolver::CTMC};
    }
    if (t.prio == "preempt")
        return {AutoSolver::MAM, AutoSolver::LDES, AutoSolver::CTMC, AutoSolver::SSA};
    if (t.prio == "ps") return {AutoSolver::CTMC, AutoSolver::LDES, AutoSolver::SSA};
    if (t.has_map)
        return {AutoSolver::MAM, AutoSolver::MVA, AutoSolver::FLUID, AutoSolver::LDES};
    if (t.prio == "hol")
        return {AutoSolver::MVA, AutoSolver::MAM, AutoSolver::FLUID, AutoSolver::CTMC,
                AutoSolver::LDES};
    if (t.pop_per_chain > 30.0) return {AutoSolver::FLUID, AutoSolver::MVA, AutoSolver::NC};
    // No exact solver was available at this population (the exact-first pass
    // above tried), so keep the exact-leaning approximate order.
    if (t.total_jobs > 0.0 && t.total_jobs <= kAutoExactPopulationMax)
        return {AutoSolver::NC, AutoSolver::MVA, AutoSolver::MAM};
    if (homogeneous_inf) return {AutoSolver::MVA, AutoSolver::NC, AutoSolver::FLUID};
    return {AutoSolver::MVA, AutoSolver::NC, AutoSolver::MAM, AutoSolver::FLUID, AutoSolver::LDES};
}

inline UnsupportedError no_solver(const std::string& what,
                                  const std::vector<AutoSolver>& skipped) {
    std::string msg = "SolverAUTO: no solver supports this model" + what;
    if (!skipped.empty()) {
        msg += " (the ranking preferred ";
        for (std::size_t i = 0; i < skipped.size(); ++i) {
            if (i) msg += ", ";
            msg += auto_solver_name(skipped[i]);
        }
        msg += ", which this port does not build)";
    }
    return UnsupportedError(msg);
}

}  // namespace detail

/**
 * `chooseAvgSolverHeur`, the Network arm: exact first at small populations,
 * then the trait-keyed ranking, then the whole pool in the global order.
 *
 * The INF branch of avgOrder consults `has_homogeneous_scheduling`, which is
 * `nstations == 1` here and in MATLAB alike; see the header note.
 */
template <class T>
AutoChoice auto_choose_avg_solver_ex(const qn::NetworkStruct<T>& sn) {
    using lang::SchedStrategy;
    const AutoTraits t = auto_traits(sn);
    AutoChoice out;

    // Small populations: an approximation buys nothing there, so take an exact
    // solver whenever one is available. The 'exact' method name is what makes this a
    // claim rather than a preference -- MVA and NC reject it without a
    // product-form solution, and CTMC is screened for a chain that fits.
    if (t.total_jobs > 0.0 && t.total_jobs <= kAutoExactPopulationMax) {
        const std::vector<AutoSolver> exact_order =
            t.has_cache
                ? std::vector<AutoSolver>{AutoSolver::NC, AutoSolver::MVA, AutoSolver::CTMC}
                : std::vector<AutoSolver>{AutoSolver::MVA, AutoSolver::NC, AutoSolver::CTMC};
        if (auto_ranked(exact_order, sn, "exact", out)) return out;
    }

    const std::vector<AutoSolver> order =
        detail::avg_order(t, sn.has_homogeneous_scheduling(SchedStrategy::INF));
    if (auto_ranked(order, sn, "default", out)) return out;

    // Nothing in the ranked list is feasible: fall back to the whole pool in
    // the global order rather than returning nothing.
    const std::vector<AutoSolver> pool = {AutoSolver::MVA,  AutoSolver::NC,   AutoSolver::MAM,
                                          AutoSolver::FLUID, AutoSolver::LDES, AutoSolver::CTMC,
                                          AutoSolver::SSA};
    if (auto_ranked(pool, sn, "default", out)) return out;
    throw detail::no_solver("", pool);
}

template <class T>
AutoSolver auto_choose_avg_solver(const qn::NetworkStruct<T>& sn) {
    return auto_choose_avg_solver_ex(sn).solver;
}

// ---------------------------------------------------------------------------
// chooseSolverHeur.m, the Network arm
// ---------------------------------------------------------------------------

namespace detail {

inline bool in_list(const std::string& m, const char* const* tab, std::size_t n) {
    for (std::size_t i = 0; i < n; ++i)
        if (m == tab[i]) return true;
    return false;
}

/** The average-metric getters, which defer to the feature tree. */
inline bool is_avg_method(const std::string& m) {
    static const char* kAvg[] = {
        "getAvgChainTable", "getAvgTputTable", "getAvgRespTTable", "getAvgUtilTable",
        "getAvgSysTable", "getAvgNodeTable", "getAvgTable", "getAvgTableLayered", "getAvg",
        "getAvgChain", "getAvgSys", "getAvgNode", "getAvgNodeChain", "getAvgArvRChain",
        "getAvgQLenChain", "getAvgUtilChain", "getAvgRespTChain", "getAvgTputChain",
        "getAvgSysRespT", "getAvgSysTput", "getAvgQLen", "getAvgUtil", "getAvgRespT",
        "getAvgResidT", "getAvgWaitT", "getAvgTput", "getAvgArvR", "getAvgQLenTable",
        "getAvgResidTChain", "getAvgNodeQLenChain", "getAvgNodeUtilChain",
        "getAvgNodeRespTChain", "getAvgNodeResidTChain", "getAvgNodeTputChain",
        "getAvgNodeArvRChain", "getAvgNodeChainTable", "getResults", "hasResults",
        "getAvgHandles", "getTranHandles", "getAvgQLenHandles", "getAvgUtilHandles",
        "getAvgRespTHandles", "getAvgTputHandles", "getAvgArvRHandles", "getAvgResidTHandles",
        "getMethodFeatureSet", "supportsModelMethod", "isStochasticMethod", "libraries",
        "showLibraryAttribution", "citations"};
    return in_list(m, kAvg, sizeof(kAvg) / sizeof(*kAvg));
}

/** The ensemble getters, which exist for LayeredNetwork models only. */
inline bool is_ensemble_method(const std::string& m) {
    static const char* kEns[] = {"getEnsembleAvg", "getEnsembleAvgTables", "getSolver",
                                 "setSolver", "getNumberOfModels", "getIteration",
                                 "get_state", "set_state", "update_solver"};
    return in_list(m, kEns, sizeof(kEns) / sizeof(*kEns));
}

inline bool is_cdf_method(const std::string& m) {
    return m == "getCdfRespT" || m == "getCdfPassT" || m == "getPerctRespT";
}

inline bool is_tran_prob_method(const std::string& m) {
    return m == "getTranProb" || m == "getTranProbSys" || m == "getTranProbAggr" ||
           m == "getTranProbSysAggr";
}

inline bool is_prob_method(const std::string& m) {
    return m == "getProb" || m == "getProbAggr" || m == "getProbSys" || m == "getProbSysAggr" ||
           m == "getProbMarg" || m == "getProbNormConstAggr";
}

inline bool is_sample_method(const std::string& m) { return m == "sample" || m == "sampleSys"; }

inline bool is_sample_aggr_method(const std::string& m) {
    return m == "sampleAggr" || m == "sampleSysAggr";
}

inline bool is_cache_metric_method(const std::string& m) {
    static const char* kTab[] = {"getAvgCacheTable", "getAvgCacheT", "getAvgItemTable",
                                 "getAvgItemT", "cacheAvgT", "itemAvgT", "aCaT", "aIT"};
    return in_list(m, kTab, sizeof(kTab) / sizeof(*kTab));
}

inline bool is_loss_metric_method(const std::string& m) {
    static const char* kTab[] = {"getAvgLossTable",       "getAvgLossT", "getAvgRegionLossTable",
                                 "getAvgRegionLossT",     "lossAvgT",    "regionLossAvgT",
                                 "aLT",                   "aRLT"};
    return in_list(m, kTab, sizeof(kTab) / sizeof(*kTab));
}

inline bool is_orbit_metric_method(const std::string& m) {
    static const char* kTab[] = {"getAvgOrbitTable", "getAvgOrbitT", "getAvgOrbit", "orbitAvgT",
                                 "aOT"};
    return in_list(m, kTab, sizeof(kTab) / sizeof(*kTab));
}

inline bool is_moment_method(const std::string& m) {
    static const char* kTab[] = {"getMomentTable",   "getMomentChainTable", "getMomentStationTable",
                                 "getMomentT",       "getMomentChainT",     "getMomentStationT",
                                 "momentT",          "momentChainT",        "momentStationT",
                                 "mT",               "mCT",                 "mST"};
    return in_list(m, kTab, sizeof(kTab) / sizeof(*kTab));
}

inline bool is_sens_method(const std::string& m) {
    static const char* kTab[] = {"getSensitivityTable", "getSensitivityT", "sensitivityT", "sT",
                                 "supportsExactSensitivity"};
    return in_list(m, kTab, sizeof(kTab) / sizeof(*kTab));
}

}  // namespace detail

/**
 * `chooseNetworkSolver` in chooseSolverHeur.m: the getter names the metric
 * family, the family names a ranking, and the average family alone consults the
 * feature tree. A ranking that yields nothing falls back to the average
 * heuristic, as the reference does, rather than refusing.
 */
template <class T>
AutoChoice auto_choose_solver_heur(const qn::NetworkStruct<T>& sn, const std::string& method) {
    using lang::SchedStrategy;
    if (detail::is_avg_method(method)) return auto_choose_avg_solver_ex(sn);
    if (detail::is_ensemble_method(method))
        throw InputError("SolverAUTO: method '" + method +
                         "' is only available for LayeredNetwork models");

    std::vector<AutoSolver> order;
    if (method == "getTranAvg") {
        order = {AutoSolver::FLUID, AutoSolver::LDES};
    } else if (detail::is_cdf_method(method)) {
        // NC gives the exact passage-time distribution on FCFS product form;
        // otherwise Fluid is the smooth approximation, then the simulators.
        if (sn.has_homogeneous_scheduling(SchedStrategy::FCFS) && sn.has_product_form())
            order = {AutoSolver::NC, AutoSolver::FLUID, AutoSolver::LDES};
        else
            order = {AutoSolver::FLUID, AutoSolver::LDES};
    } else if (method == "getTranCdfPassT" || method == "getTranCdfRespT") {
        order = {AutoSolver::FLUID, AutoSolver::LDES};
    } else if (detail::is_tran_prob_method(method)) {
        order = {AutoSolver::CTMC};
    } else if (detail::is_sample_method(method)) {
        order = {AutoSolver::SSA, AutoSolver::LDES};
    } else if (detail::is_sample_aggr_method(method)) {
        order = {AutoSolver::LDES, AutoSolver::SSA};
    } else if (detail::is_prob_method(method)) {
        if (sn.has_product_form())
            order = {AutoSolver::NC, AutoSolver::CTMC, AutoSolver::LDES};
        else
            order = {AutoSolver::CTMC, AutoSolver::LDES};
    } else if (detail::is_cache_metric_method(method)) {
        order = {AutoSolver::NC, AutoSolver::MVA, AutoSolver::FLUID, AutoSolver::CTMC,
                 AutoSolver::LDES};
    } else if (detail::is_loss_metric_method(method)) {
        order = {AutoSolver::LDES, AutoSolver::CTMC, AutoSolver::SSA};
    } else if (detail::is_orbit_metric_method(method)) {
        order = {AutoSolver::MVA, AutoSolver::CTMC, AutoSolver::LDES};
    } else if (detail::is_moment_method(method)) {
        order = {AutoSolver::MVA, AutoSolver::NC, AutoSolver::CTMC, AutoSolver::LDES};
    } else if (detail::is_sens_method(method)) {
        order = {AutoSolver::FLUID, AutoSolver::MVA, AutoSolver::NC};
    } else {
        return auto_choose_avg_solver_ex(sn);
    }

    AutoChoice out;
    if (auto_ranked(order, sn, "default", out)) return out;
    // No solver in the metric's ranking supports the model: the average
    // heuristic is the floor.
    return auto_choose_avg_solver_ex(sn);
}

/** `chooseSolverExact`: the ranking restricted to solvers that answer exactly. */
template <class T>
AutoChoice auto_choose_solver_exact(const qn::NetworkStruct<T>& sn, const std::string& method) {
    const AutoTraits t = auto_traits(sn);
    std::vector<AutoSolver> order;
    if (detail::is_tran_prob_method(method) || detail::is_cdf_method(method) ||
        method == "getTranCdfPassT" || method == "getTranCdfRespT" || method == "getTranAvg") {
        order = {AutoSolver::CTMC};
    } else if (detail::is_sample_method(method) || detail::is_sample_aggr_method(method)) {
        // A sample path is exact in distribution, not in the mean.
        order = {AutoSolver::SSA, AutoSolver::LDES};
    } else if (detail::is_prob_method(method)) {
        if (t.is_product_form) order = {AutoSolver::NC, AutoSolver::CTMC};
        else order = {AutoSolver::CTMC};
    } else if (t.is_product_form && !t.has_multi_server) {
        order = {AutoSolver::NC, AutoSolver::CTMC};
    } else {
        order = {AutoSolver::CTMC, AutoSolver::NC};
    }
    AutoChoice out;
    if (auto_ranked(order, sn, "exact", out)) return out;
    throw detail::no_solver(" exactly for method '" + method +
                                "'; use the default heuristic for the approximation",
                            order);
}

/** `chooseSolverSim`: the ranking restricted to simulators. */
template <class T>
AutoChoice auto_choose_solver_sim(const qn::NetworkStruct<T>& sn, const std::string& method) {
    const std::vector<AutoSolver> order =
        detail::is_sample_method(method)
            ? std::vector<AutoSolver>{AutoSolver::SSA, AutoSolver::LDES}
            : std::vector<AutoSolver>{AutoSolver::LDES, AutoSolver::SSA};
    AutoChoice out;
    if (auto_ranked(order, sn, "default", out)) return out;
    throw detail::no_solver(" by simulation", order);
}

/**
 * `chooseSolver`: the selection mode picks the ranking, and every mode but the
 * two learned ones keeps the heuristic as its floor.
 */
template <class T>
AutoChoice auto_choose_solver_mode(const qn::NetworkStruct<T>& sn, const std::string& method,
                                   AutoMode mode) {
    AutoChoice out;
    switch (mode) {
        case AutoMode::EXACT: return auto_choose_solver_exact(sn, method);
        case AutoMode::SIM: return auto_choose_solver_sim(sn, method);
        case AutoMode::FAST:
            // Cheapest analytical answer; the heuristic is the floor for
            // metrics no mean-value solver can serve.
            if (auto_ranked({AutoSolver::MVA, AutoSolver::NC, AutoSolver::FLUID, AutoSolver::MAM},
                            sn, "default", out))
                return out;
            return auto_choose_solver_heur(sn, method);
        case AutoMode::ACCURATE:
            // A smooth or matrix-analytic answer preferred over the fastest.
            if (auto_ranked({AutoSolver::FLUID, AutoSolver::MAM, AutoSolver::CTMC,
                             AutoSolver::LDES},
                            sn, "default", out))
                return out;
            return auto_choose_solver_heur(sn, method);
        case AutoMode::HEUR: break;
    }
    return auto_choose_solver_heur(sn, method);
}

/** The heuristic Network arm, by getter name. */
template <class T>
AutoSolver auto_choose_solver(const qn::NetworkStruct<T>& sn, const std::string& method) {
    return auto_choose_solver_heur(sn, method).solver;
}

/**
 * `delegate`'s proposed order: the chosen solver, then every feasible candidate
 * in slot order. Duplicates are dropped, since the reference retries the chosen
 * solver only once in practice.
 */
template <class T>
std::vector<AutoSolver> auto_proposed_solvers(const qn::NetworkStruct<T>& sn,
                                              const std::string& method, AutoMode mode) {
    const AutoChoice chosen = auto_choose_solver_mode(sn, method, mode);
    std::vector<AutoSolver> out(1, chosen.solver);
    const std::vector<AutoSolver> cand = auto_candidates(sn);
    for (std::size_t i = 0; i < cand.size(); ++i)
        if (std::find(out.begin(), out.end(), cand[i]) == out.end()) out.push_back(cand[i]);
    return out;
}

// ---------------------------------------------------------------------------
// The LayeredNetwork and Environment arms
// ---------------------------------------------------------------------------

struct AutoLayeredChoice {
    AutoLayered solver = AutoLayered::LN_MVA;
    std::vector<AutoLayered> skipped;
};

namespace detail {

inline bool layered_ranked(const std::vector<AutoLayered>& order, AutoLayeredChoice& out) {
    std::vector<AutoLayered> skipped;
    for (std::size_t k = 0; k < order.size(); ++k) {
        if (!auto_layered_is_available(order[k])) {
            skipped.push_back(order[k]);
            continue;
        }
        out.solver = order[k];
        out.skipped = skipped;
        return true;
    }
    return false;
}

}  // namespace detail

/**
 * `chooseLayeredSolver` plus the LayeredNetwork arm of `chooseAvgSolverHeur`.
 *
 * There is no feature-set gate on this path: the reference gates a layered
 * candidate with `SolverLN.supports(model)`, a LayeredNetwork-level check this
 * port does not carry, so availability is the only screen. A cache task is what
 * inverts the analytical order -- the cache layer is where NC beats MVA.
 */
inline AutoLayeredChoice auto_choose_layered_solver(const std::string& method,
                                                    bool has_cache_task,
                                                    AutoMode mode = AutoMode::HEUR) {
    AutoLayeredChoice out;
    std::vector<AutoLayered> order;
    if (mode == AutoMode::EXACT) {
        // No layered solver is exact; NC layers are the closest available.
        order = {AutoLayered::LN_NC, AutoLayered::LN_MVA};
    } else if (mode == AutoMode::SIM) {
        // lqsim is the layered simulator; the LN solvers are the fallback.
        order = {AutoLayered::LQNS, AutoLayered::LN_MVA, AutoLayered::LN_NC};
    } else if (method == "getTranAvg" || detail::is_cdf_method(method) ||
               method == "getTranCdfPassT" || method == "getTranCdfRespT") {
        order = {AutoLayered::LN_FLUID, AutoLayered::LN_MVA, AutoLayered::LN_NC};
    } else if (detail::is_sample_method(method) || detail::is_sample_aggr_method(method)) {
        order = {AutoLayered::LQNS, AutoLayered::LN_MVA, AutoLayered::LN_NC};
    } else if (detail::is_prob_method(method) || detail::is_tran_prob_method(method)) {
        order = {AutoLayered::LN_NC, AutoLayered::LN_MVA};
    } else if (detail::is_ensemble_method(method)) {
        order = {AutoLayered::LN_MVA, AutoLayered::LN_NC, AutoLayered::LN_MAM,
                 AutoLayered::LN_FLUID, AutoLayered::LQNS};
    } else if (has_cache_task) {
        // The LayeredNetwork arm of chooseAvgSolverHeur, which returns the NC
        // layer solver outright on a cache task.
        out.solver = AutoLayered::LN_NC;
        return out;
    } else {
        order = {AutoLayered::LQNS, AutoLayered::LN_MVA, AutoLayered::LN_NC, AutoLayered::LN_MAM,
                 AutoLayered::LN_FLUID};
    }
    if (has_cache_task && mode != AutoMode::SIM) {
        // A layered cache model needs the NC layer solver whichever metric was
        // asked for: it is the only one that solves the cache layer.
        std::vector<AutoLayered> promoted(1, AutoLayered::LN_NC);
        for (std::size_t i = 0; i < order.size(); ++i)
            if (order[i] != AutoLayered::LN_NC) promoted.push_back(order[i]);
        order = promoted;
    }
    if (detail::layered_ranked(order, out)) return out;
    throw UnsupportedError(
        "SolverAUTO: no LayeredNetwork solver in this port serves method '" + method +
        "'; the reference's ranking is served by SolverLQNS, which shells out to an external "
        "binary this port does not wrap");
}

struct AutoEnvChoice {
    AutoEnv solver = AutoEnv::ENV_FLUID;
    std::vector<AutoEnv> skipped;
};

/**
 * The Environment arm of `chooseSolverHeur` / `chooseSolverExact` /
 * `chooseSolverSim`.
 *
 * Fluid leads the heuristic ranking because the blending method is transient:
 * it restarts each stage from the mean state the previous one left, and an
 * inner solver without a transient analysis returns zeros for every blended
 * metric. That is also the only stage engine this port builds.
 */
inline AutoEnvChoice auto_choose_env_solver(const std::string& method,
                                            AutoMode mode = AutoMode::HEUR) {
    std::vector<AutoEnv> order;
    if (mode == AutoMode::EXACT) order = {AutoEnv::ENV_NC, AutoEnv::ENV_MVA};
    else if (mode == AutoMode::SIM) order = {AutoEnv::ENV_MVA, AutoEnv::ENV_NC, AutoEnv::ENV_FLUID};
    else order = {AutoEnv::ENV_FLUID, AutoEnv::ENV_MVA, AutoEnv::ENV_NC};

    AutoEnvChoice out;
    std::vector<AutoEnv> skipped;
    for (std::size_t k = 0; k < order.size(); ++k) {
        if (!auto_env_is_available(order[k])) {
            skipped.push_back(order[k]);
            continue;
        }
        out.solver = order[k];
        out.skipped = skipped;
        return out;
    }
    throw UnsupportedError(
        "SolverAUTO: the Environment ranking for method '" + method +
        "' selects an MVA or NC stage solver, and SolverENV in this port solves a stage with the "
        "fluid analyzer (or, under the state-vector coupling, an explicit chain); rerun with "
        "-s env, whose default coupling is the mean-field one");
}

}  // namespace autosolver
}  // namespace line

#endif  // LINE_SOLVERS_AUTO_SOLVER_AUTO_H
