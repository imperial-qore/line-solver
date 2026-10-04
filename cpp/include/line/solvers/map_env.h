/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MAP_ENV_H
#define LINE_SOLVERS_MAP_ENV_H

/**
 * @file
 * @ingroup line_solvers
 * Port of `@NetworkSolver/mapEnvApprox.m`: the solver-agnostic
 * random-environment approximation of a network whose arrival or service
 * processes are non-renewal.
 *
 * `map2renv` turns each modulated process into a set of environment stages in
 * which that process is exponential with the phase-conditional intensity, and
 * SolverENV recombines the stages. The stage models are exponential and carry
 * the base model's structure unchanged, so THE STAGE SOLVER IS THE CALLER
 * ITSELF: nothing here reads a solver internal, and any runner that can answer
 * `getAvg` on a struct can be driven through it.
 *
 * WHERE IT SITS. The reference's `runAnalyzer` still REFUSES a MAP model; only
 * `getAvg` falls back (`@NetworkSolver/getAvg.m` lines 150-151). That placement
 * is forced here as well as faithful: `env/solver_env.h` includes
 * `ln/solver_ln.h`, which includes the MVA and NC runners, so the gate cannot
 * live where the feature gate lives and has to stand in front of the runners
 * instead. C++ has no `getAvg`, so the three ladders that play its part -- the
 * facade's `avg_table`, the CLI's `run_avg_engine`, and the CLI's per-solver
 * `-a avg` arms -- each wrap their runner in `run_avg` below. All three, or
 * `line-cli -a avg` and `-a cache` would accept different models for one
 * command line.
 *
 * WHAT IS DOUBLE-ONLY, and why the refusal stays. `SolverEnvLimit::init` and
 * `SolverEnv::init` both refuse `T != double` -- a fluid stage integrates its
 * drift with LSODA and a mean-field sojourn is a `double` quadrature -- while
 * the MVA and NC runners this wraps are multiprecision. So `--arith rational`
 * on a MAP model keeps the PLAIN rejection rather than falling back, and says
 * so by name. See `_kb/14-cpp-multiprecision.md`.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/io/map2renv.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/env/env_dispatch.h"
#include "line/solvers/map_env_gate.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/util/error.h"

namespace line {
namespace solvers {

namespace map_env_detail {

/**
 * `mapImageKind`: an MMPP image preserves the modulating chain exactly; a
 * general MAP image aggregates the event-epoch phase jumps into the phase
 * generator, and the correlation between the jumps and the event stream is lost.
 */
template <class T>
const char* image_kind(const io::Map2RenvInfo<T>& info) {
    return info.is_mmpp ? "exact-modulation" : "intensity-matched";
}

/** Longest mean stage sojourn of the environment, the mean-field horizon's scale. */
template <class T>
double max_hold_time(const env::Environment<T>& e) {
    double t = 0.0;
    for (std::size_t s = 0; s < e.nstages(); ++s)
        t = std::max(t, mam::map_mean(e.hold_time[s].map()));
    if (!std::isfinite(t) || !(t > 0.0))
        throw InputError(
            "map_env: the environment image has no finite stage sojourn, so no integration "
            "horizon can be set");
    return t;
}

/**
 * `selectEnvLimit`: the timescale test between the two closed-form limits.
 *
 * Compare the mean stage holding time of the environment with the relaxation
 * time of the model, taken as the time the slowest station needs to clear the
 * jobs it can hold. A stage that OUTLIVES the relaxation time lets each stage
 * reach its own steady state, which is the quasi-stationary regime (`dec`); a
 * stage that expires first leaves the model responding to the mean rate only,
 * which is the rate-averaged regime (`avg`). The closed population enters the
 * relaxation time because a closed queue drains in N services, so a populous
 * model relaxes far more slowly than one service time.
 */
template <class T>
std::string select_env_limit(const qn::NetworkStruct<T>& sn, const env::Environment<T>& e) {
    const std::size_t E = e.nstages();
    std::vector<double> exit_rate;
    for (std::size_t a = 0; a < E; ++a) {
        double r = 0.0;
        for (std::size_t b = 0; b < E; ++b) {
            if (a == b) continue;
            const env::EnvArc<T>& arc = e.arc(a, b);
            if (!arc.enabled) continue;
            // `1 / env{e,h}.getMean()`, the reference's `getRate()`: the arc's
            // own rate, whatever family carries it.
            const double m = num_traits<T>::to_double(arc.dist.mean);
            if (m > 0.0) r += 1.0 / m;
        }
        if (r > 0.0) exit_rate.push_back(r);
    }
    // An ABSORBING environment never leaves the stage it is in, so every stage
    // IS its own steady state and the quasi-stationary limit is exact.
    if (exit_rate.empty()) return "dec";
    double tau_env = 0.0;
    for (std::size_t i = 0; i < exit_rate.size(); ++i) tau_env += 1.0 / exit_rate[i];
    tau_env /= static_cast<double>(exit_rate.size());

    double minrate = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < sn.rates.rows(); ++i)
        for (std::size_t r = 0; r < sn.rates.cols(); ++r) {
            const double v = num_traits<T>::to_double(sn.rates(i, r));
            if (std::isfinite(v) && v > 0.0) minrate = std::min(minrate, v);
        }
    if (!std::isfinite(minrate)) return "dec";
    const std::vector<double> nj = sn.njobs();
    double totn = 0.0;
    for (std::size_t k = 0; k < nj.size(); ++k)
        if (std::isfinite(nj[k])) totn += nj[k];
    const double tau_sys = (1.0 + totn) / minrate;
    return tau_env >= tau_sys ? "dec" : "avg";
}

/**
 * The base-to-stage index correspondence `map_env_approx` relies on, asserted
 * rather than assumed.
 *
 * `map2renv` copies the struct and rewrites only service entries, so station,
 * class and node indices are identical between the base and every stage. Every
 * column of the blend is read at the base's index, so a later change there
 * would shift each one SILENTLY; this is what makes the correspondence an
 * invariant rather than a coincidence.
 */
template <class T>
void assert_stage_correspondence(const qn::NetworkStruct<T>& sn, const env::Environment<T>& e) {
    for (std::size_t s = 0; s < e.nstages(); ++s) {
        const qn::NetworkStruct<T>& st = e.stage(s).model;
        if (st.nstations != sn.nstations || st.nclasses != sn.nclasses ||
            st.nodes.size() != sn.nodes.size())
            throw InputError("map_env: stage " + std::to_string(s + 1) +
                             " of the environment image has a different station, class or node "
                             "count from the base model; every metric is blended at the base's "
                             "own index");
        for (std::size_t i = 0; i < sn.nodes.size(); ++i)
            if (st.nodes[i].name != sn.nodes[i].name)
                throw InputError("map_env: node " + std::to_string(i + 1) + " of stage " +
                                 std::to_string(s + 1) + " is '" + st.nodes[i].name +
                                 "' where the base model has '" + sn.nodes[i].name +
                                 "'; the image must keep the base's node order");
    }
}

}  // namespace map_env_detail

/**
 * `mapEnvApprox`: solve the model through the random-environment image of its
 * non-renewal processes.
 *
 * @param sn       the base model
 * @param solver   the calling solver's label ("SolverMVA", "SolverNC", ...),
 *                 which decides the stage backend and the transient capability
 * @param cfg      the caller's map_env knobs
 * @param stage_fn runs ONE stage with the calling solver and returns its
 *                 Q/U/T and cache surface; the C++ spelling of the reference's
 *                 `feval(class(self), stageModel, innerOptions)` factory
 * @param requested_method the method the caller asked for, echoed in the result
 *
 * THE MEAN-FIELD COUPLING DOES NOT TAKE `stage_fn`, and that is not an omission:
 * it carries the marginal means and the RMF cache transient across a switch,
 * objects only its own backend produces, so it runs `SolverEnv`'s own fluid or
 * ctmc stage instead. A solver with neither backend asking for `meanfield`
 * EXPLICITLY is refused by name; `auto` never picks it for such a solver.
 */
template <class T, class StageFn>
mva::AvgResult<T> map_env_approx(const qn::NetworkStruct<T>& sn, const std::string& solver,
                                 const MapEnvConfig& cfg, StageFn stage_fn,
                                 const std::string& requested_method = "default") {
    if (!std::is_same<T, double>::value)
        throw UnsupportedError(
            "map_env: the random-environment fallback for a MAP/MMPP/MMAP/MPH process solves "
            "each stage through SolverENV, whose couplings are double-only (a fluid stage "
            "integrates with LSODA and a mean-field sojourn is a double quadrature). Rerun with "
            "--arith double, or set map_env='off' to keep the plain rejection");

    io::Map2RenvInfo<T> info;
    env::Environment<T> e =
        io::map2renv(sn, &info, cfg.max_stages > 0 ? cfg.max_stages : static_cast<std::size_t>(64));
    map_env_detail::assert_stage_correspondence(sn, e);

    const std::string backend = map_env_stage_backend(solver);
    std::string method = cfg.method.empty() ? std::string("auto") : cfg.method;
    if (method == "auto") {
        // The mean-field coupling is the only one that MODELS the phase switch
        // rather than taking a limit of it, so it is preferred wherever it can
        // run. It can run only where BOTH hold: the solver produces transient
        // averages, and this port has a stage backend for it.
        method = (supports_transient_analysis(solver) && !backend.empty())
                     ? std::string("meanfield")
                     : map_env_detail::select_env_limit(sn, e);
    }
    if (method != "dec" && method != "avg" && method != "meanfield")
        throw UnsupportedError("map_env: map_env_method='" + method +
                               "' is not a supported environment recombination; use 'meanfield', "
                               "'dec', 'avg' or 'auto'");
    if (method == "meanfield" && backend.empty())
        throw UnsupportedError(
            "map_env: the mean-field environment coupling integrates each stage over its sojourn, "
            "so it needs a stage backend that produces transient averages, and this port wires "
            "'fluid' and 'ctmc' only -- " +
            solver +
            " has neither. Use map_env_method='dec' or 'avg', which solve each stage in steady "
            "state with the calling solver itself");

    env::EnvOptions o;
    o.method = method == "meanfield" ? "meanfield" : method;
    if (method == "meanfield") {
        o.stage_solver = backend;
        // The stage transients are weighted by the sojourn density over the
        // integration grid, so the horizon must cover the sojourn distribution;
        // beyond it the weights vanish and the extra span is inert.
        o.timespan_end = 20.0 * map_env_detail::max_hold_time(e);
    }
    const env::EnvAnalyzerSolution<T> r =
        method == "meanfield" ? env::solver_env(e, o)
                              : env::solver_env(e, o, env::EnvStageAvgFn<T>(stage_fn));

    mva::AvgResult<T> out;
    out.QN = r.QN;
    out.UN = r.UN;
    out.TN = r.TN;
    out.cache = r.cache;
    out.method = requested_method;
    out.actualmethod = "env." + method;

    const std::size_t M = sn.nstations, K = sn.nclasses;
    const T zero = num_traits<T>::from_int(0);
    out.RN = Matrix<T>(M, K, zero);
    out.AN = Matrix<T>(M, K, num_traits<T>::from_double(std::numeric_limits<double>::quiet_NaN()));
    out.XN.assign(K, zero);
    out.CN.assign(K, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k)
            if (num_traits<T>::to_double(out.TN(i, k)) > 0.0)
                out.RN(i, k) = T(out.QN(i, k) / out.TN(i, k));
    for (std::size_t k = 0; k < K; ++k) {
        const std::size_t ist = sn.classes[k].refstat;
        if (ist >= 1 && ist <= M) out.XN[k] = out.TN(ist - 1, k);
        if (num_traits<T>::to_double(out.XN[k]) > 0.0) {
            T q = zero;
            for (std::size_t i = 0; i < M; ++i) q = T(q + out.QN(i, k));
            out.CN[k] = T(q / out.XN[k]);
        }
    }
    out.WN = mva::sn_get_residt_from_respt<T>(sn, out.RN);

    // THE APPROXIMATION NOTICE IS NOT OPTIONAL. The caller asked for a solver
    // that cannot consume this model's processes and is being answered by an
    // approximation of a different model; `AvgResult::warning` exists for
    // exactly this and is printed by the CLI and the avg table.
    out.warning = "This solver has no native support for the non-renewal (MAP/MMPP) processes of "
                  "this model; the reported averages come from its " +
                  method + " random-environment approximation (" + std::to_string(info.nstages) +
                  " stages, " + map_env_detail::image_kind(info) +
                  " image). Set map_env='off' to reject the model instead.";

    // The system metrics read the reference station, which for an open class is
    // the Source. A stage solver whose transient does not report Source
    // throughput leaves it at zero under the mean-field coupling, and the zero
    // propagates into XN and CN. REPORT THAT rather than substituting the
    // arrival rate, which would hide whose metric is missing.
    const std::vector<double> nj = sn.njobs();
    std::string openzero;
    for (std::size_t k = 0; k < K; ++k) {
        if (k < nj.size() && std::isfinite(nj[k])) continue;
        if (num_traits<T>::to_double(out.XN[k]) != 0.0) continue;
        bool any = false;
        for (std::size_t i = 0; i < M; ++i)
            if (num_traits<T>::to_double(out.QN(i, k)) > 0.0) any = true;
        const std::size_t ist = sn.classes[k].refstat;
        if (!any && ist >= 1 && ist <= sn.rates.rows() && k < sn.rates.cols()) {
            const double lam = num_traits<T>::to_double(sn.rates(ist - 1, k));
            any = std::isfinite(lam) && lam > 0.0;
        }
        if (!any) continue;
        if (!openzero.empty()) openzero += ", ";
        openzero += std::to_string(k + 1);
    }
    if (!openzero.empty())
        out.warning += " The " + method +
                       " environment coupling returned no reference-station throughput for open "
                       "class(es) " +
                       openzero +
                       ", so their system throughput and system response time are reported as "
                       "zero; use map_env_method='dec' or 'avg' for system-level metrics.";
    return out;
}

/**
 * The `getAvg` funnel: run the model, or its environment image when the ONLY
 * thing in the way is a non-renewal process.
 *
 * THREE OUTCOMES, AND THE THIRD IS THE SUBTLE ONE. Supported -> run. Needs the
 * image -> `map_env_approx`. Otherwise -> RUN ANYWAY, and let the runner raise
 * its own refusal in its own words. Raising a second copy of the refusal here
 * would mean two messages for one condition, drifting apart as the feature sets
 * move; `needs_map_env` is deliberately a question about whether the IMAGE
 * helps, not a gate on whether the model is supported.
 *
 * NOT WIRED INTO THE `-s ba` ARM, deliberately. A bound request must be
 * answered with a bound, and the environment image is an approximation of the
 * model, so its bounds do not bracket the original one. The reference refuses
 * `SolverBA` inside `needsMapEnv` itself; here the exclusion is the absence of
 * this wrapper at the three BA call sites, each of which says so.
 */
template <class T, class Run, class StageFn>
mva::AvgResult<T> run_avg(const qn::NetworkStruct<T>& sn, const std::string& solver,
                          const qn::FeatureSet& declared, const MapEnvConfig& cfg, Run run,
                          StageFn stage_fn, const std::string& requested_method = "default") {
    const MapEnvDecision d = needs_map_env(declared, sn, cfg);
    if (!d.needed) return run(sn);
    return map_env_approx<T>(sn, solver, cfg, stage_fn, requested_method);
}

}  // namespace solvers
}  // namespace line

#endif  // LINE_SOLVERS_MAP_ENV_H
