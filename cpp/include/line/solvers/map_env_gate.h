/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MAP_ENV_GATE_H
#define LINE_SOLVERS_MAP_ENV_GATE_H

/**
 * @file
 * @ingroup line_solvers
 * The map_env decision, as a predicate over a feature set.
 *
 * WHAT IT DECIDES. `@NetworkSolver/needsMapEnv.m` asks one question: are the
 * ONLY features this solver cannot consume non-renewal arrival or service
 * processes? If so the model becomes supported the moment each modulated
 * process is frozen into an exponential environment stage, and `map2renv` plus
 * SolverENV can answer it. If anything else is unsupported the environment
 * image would not make the model solvable, so the original refusal stands and
 * is raised by the runner in its own words.
 *
 * WHY IT IS ITS OWN HEADER, and a small one. Two callers need it and they sit
 * at opposite ends of the include graph: the driver in `map_env.h`, which pulls
 * in `env/` and therefore everything below it, and the AUTO / `find_solver`
 * reporting layer, which only wants to know whether a solver would answer. This
 * header depends on `lang/qn/feature_set.h` alone, so the reporting layer can
 * ask the question without linking a single solver.
 *
 * THE FOUR NAMES ARE THE REFERENCE'S OWN `mapTokens` and are a closed list:
 * MAP, MMPP2, MMAP and MPH. DMAP, BMAP, MAPt, PHt, MMAPt, MPHt and BMMAPt are
 * deliberately NOT here. `map2renv` freezes a stationary, continuous-time,
 * single-arrival modulating chain into one exponential stage per phase; a
 * discrete-time chain has no such stage, a batch process releases more than one
 * job per epoch and so is not described by a rate, and a time-inhomogeneous one
 * has no stationary phase to freeze at all. Adding a name here without adding
 * the transform would send a model into an image that does not represent it.
 */

#include <string>
#include <vector>

#include "line/lang/qn/feature_set.h"

namespace line {
namespace solvers {

/**
 * The caller-facing map_env knobs, `options.config.map_env` and friends.
 *
 * `mode` is `"default"` (fall back when the gate says so) or `"off"` (keep the
 * plain rejection). `method` is the SolverENV coupling the driver asks for,
 * with `"auto"` choosing between the mean-field one and the closed-form limits
 * the way `selectEnvLimit` does. `max_stages` caps how many stages the image
 * may have, since one stage per modulating phase of every non-renewal process
 * is a product and can be large.
 */
struct MapEnvConfig {
    std::string mode = "default";
    std::string method = "auto";
    std::size_t max_stages = 0;  ///< 0 = the transform's own default cap
    bool off() const { return mode == "off"; }
};

/** The four process names an environment image can represent; the order is the registry's. */
inline const std::vector<qn::Feature>& map_env_tokens() {
    static const std::vector<qn::Feature> t = {qn::Feature::MAP, qn::Feature::MMAP,
                                               qn::Feature::MMPP2, qn::Feature::MPH};
    return t;
}

/** True when `f` is one of the four. */
inline bool is_map_env_token(qn::Feature f) {
    const std::vector<qn::Feature>& t = map_env_tokens();
    for (std::size_t i = 0; i < t.size(); ++i)
        if (t[i] == f) return true;
    return false;
}

/** What the gate decided, and what it decided it about. */
struct MapEnvDecision {
    bool needed = false;
    /** The unsupported features, in registry order; all four-token when `needed`. */
    std::vector<qn::Feature> tokens;
    /** Those names, comma-separated, for a message. */
    std::string token_names() const {
        std::string s;
        for (std::size_t i = 0; i < tokens.size(); ++i) {
            if (i) s += ", ";
            s += qn::feature_name(tokens[i]);
        }
        return s;
    }
};

/**
 * `needsMapEnv`: does this model need the environment image, and would the image
 * make it solvable?
 *
 * `declared` is the feature set of the RESOLVED method, not of the solver: MVA's
 * `default` upgrades to `rqna` on a bursty open model and rqna declares MAP, so
 * asking the base envelope would send a model into an image the run did not need.
 * The callers resolve the method first for exactly that reason.
 *
 * FALSE FOR A SUPPORTED MODEL, which is not the same as "no MAP in the model":
 * a solver that declares the process consumes it natively and must not be
 * diverted through an approximation of it.
 */
template <class T>
MapEnvDecision needs_map_env(const qn::FeatureSet& declared, const qn::NetworkStruct<T>& sn,
                             const MapEnvConfig& cfg = MapEnvConfig()) {
    MapEnvDecision d;
    if (cfg.off()) return d;
    const qn::SupportResult r =
        qn::feature_set_supports("SolverENV", declared, qn::used_lang_features(sn));
    if (r.ok || r.missing.empty()) return d;
    for (std::size_t i = 0; i < r.missing.size(); ++i)
        if (!is_map_env_token(r.missing[i])) return d;
    d.needed = true;
    d.tokens = r.missing;
    return d;
}

/**
 * Which ENV stage backend a solver's name maps to, or empty when it has none.
 *
 * NOT THE SAME QUESTION AS `supports_transient_analysis`. `SolverEnv::init`
 * admits `stage_solver` in {fluid, ctmc} only, so LDES and JMT answer the
 * transient capability TRUE and still have no backend in this port; asking the
 * capability alone would pick the mean-field coupling for a solver that cannot
 * run a stage, and the coupling would then refuse by name. Keyed on the label
 * string the feature gate already takes ("SolverMVA", "SolverNC", ...), because
 * that is the only solver identity that exists at this level.
 */
inline std::string map_env_stage_backend(const std::string& solver) {
    if (solver == "SolverFLD" || solver == "SolverFluid") return "fluid";
    if (solver == "SolverCTMC") return "ctmc";
    return std::string();
}

/**
 * `supportsTransientAnalysis`: does this solver produce transient averages?
 *
 * A CAPABILITY CLAIM, answered before any run, because the driver reads it to
 * decide whether the stages can be coupled by the mean-field analyzer (which
 * integrates each stage over its sojourn and so needs a transient) or only by
 * the two steady-state limits. The four that populate `result.Tran.Avg` in the
 * reference are FLD, CTMC, LDES and JMT.
 */
inline bool supports_transient_analysis(const std::string& solver) {
    return solver == "SolverFLD" || solver == "SolverFluid" || solver == "SolverCTMC" ||
           solver == "SolverLDES" || solver == "SolverJMT";
}

}  // namespace solvers
}  // namespace line

#endif  // LINE_SOLVERS_MAP_ENV_GATE_H
