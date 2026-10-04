/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_ENV_ENV_DISPATCH_H
#define LINE_SOLVERS_ENV_ENV_DISPATCH_H

/**
 * @file
 * @ingroup line_solvers
 * The SolverENV entry surface: a port of the analyzer selection that
 * `@@SolverENV/SolverENV.m`'s `init` performs (lines 311-327) and of the
 * compression switch `setCompression` arms (line 160, applied at line 412).
 *
 * WHAT SELECTS WHAT. The reference holds ONE solver object and swaps the
 * function handle it delegates `pre`/`analyze`/`post`/`finish` to:
 * `method` == statevec installs `solver_env_statevec_analyzer`, and every other
 * method installs `solver_env_meanfield_analyzer` -- `blend` among them since
 * 2026-09-13, when that name was repointed from the state-vector coupling to
 * the mean-field one. This port has two solver CLASSES instead of one object
 * with two handles, because the two couplings carry different objects across
 * an environment switch -- the mean-field one carries the marginal mean queue
 * lengths, the state-vector one the whole joint distribution -- and so they
 * need different options, different stage solvers and different result types.
 * This file is the switch that stands in front of them, and it is the only
 * place that knows both exist.
 *
 * WHY `SolverEnv` STILL REFUSES `statevec` BY NAME. It is the mean-field
 * coupling and nothing else; asked for the state-vector one it cannot answer,
 * and answering with the mean-field collapse would be the silent fallback the
 * whole port avoids. The refusal is therefore kept where it is and this entry
 * routes AROUND it rather than through it.
 *
 * THE THIRD THING THIS SWITCH SELECTS is neither coupling: `avg` and `dec`, the
 * closed-form fast/slow limits of `solveEnvLimit`, which bypass the fixed point
 * entirely and are ported as `SolverEnvLimit`. The reference routes them the
 * same way, at the TOP of `runAnalyzer` and before `iterate`, because there is
 * no inter-stage coupling to choose when nothing crosses a switch.
 *
 * WHICH COUPLINGS TAKE A LAYERED STAGE. Only the mean-field one, and that is
 * the reference's split rather than this port's: a `statevec` run needs
 * one generator per stage to propagate a joint law across, which an LQN does not
 * have (it decomposes into a network per layer), and the `avg`/`dec` limits read
 * a stage's station rate table, which an LQN does not carry. Both refuse a
 * LayeredNetwork stage by name, in their own `init`, through
 * `Environment::reject_lqn_stages`; so does the compression, whose macro-state
 * is a network at averaged rates. The one thing still refused everywhere is the
 * cache aggregation of the two analyzers.
 *
 * WHAT COMPRESSION IS, AND WHY IT IS A SEPARATE OVERLOAD. `setCompression` is
 * an explicit opt-in flag in the reference, defaulting to false, and the knobs
 * it reads live in `options.config` rather than beside `method`. Passing the
 * `EnvCompressOptions` is that opt-in here: there is no way to ask for the
 * compressed solve without saying which decomposition kernel and which
 * partition search it should use. The reference applies compression inside
 * `init`, i.e. BEFORE the analyzer runs and independently of which analyzer was
 * installed, so both couplings can run on a compressed environment and both do
 * here.
 */

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

#include "line/lang/qn/environment.h"
#include "line/num/number.h"
#include "line/solvers/env/solver_env.h"
#include "line/solvers/env/solver_env_limit.h"
#include "line/solvers/env/solver_env_meanfield.h"
#include "line/solvers/env/solver_env_statevec.h"
#include "line/util/matrix.h"

namespace line {
namespace env {

/**
 * What the ENV entry reports: the environment-blended metrics, and the whole
 * result of whichever coupling produced them.
 *
 * The sub-result is kept rather than reduced away because the two couplings
 * report different things beside the means -- the mean-field one the per-stage
 * exit metrics and the entry queue lengths its fixed point converged to, the
 * state-vector one the entry and sojourn-end DISTRIBUTIONS and the cache blend
 * -- and a caller who chose a coupling chose it for those.
 */
template <class T>
struct EnvAnalyzerSolution {
    /** Environment-averaged metrics, (nstations x nclasses). */
    Matrix<T> QN, UN, TN;
    /** What ran: `meanfield`, `statevec`, or the limit `avg` / `dec`. */
    std::string method;
    /** True when the environment was aggregated before the coupling ran. */
    bool compressed = false;
    /** The aggregation, when there was one; its `eps` against `epsMAX` says whether it was meaningful. */
    EnvCompression<T> compression{};
    /** Populated on the mean-field path. */
    EnvMeanfieldSolution<T> meanfield{};
    /** Populated on the state-vector path. */
    EnvStatevecSolution<T> statevec{};
    /** Populated on the closed-form limit path (`avg`, `dec`). */
    EnvLimitSolution limit{};
    /**
     * The environment-blended cache surface, whichever coupling produced it.
     *
     * A Cache does not report in (station, class) and so cannot ride in QN/UN/
     * TN; all three couplings now report it in `CacheMetrics`' one shape, keyed
     * by the Cache's NAME. EMPTY on a model with no cache, and within an entry
     * every field is optional: ABSENT MEANS NOT COMPUTED, never zero.
     */
    solvers::CacheMetrics<T> cache{};
    /**
     * ZERO ON THE LIMIT PATH, and that is the answer rather than a gap: a limit
     * reads the environment once and iterates nothing, so reporting an
     * iteration count would describe a fixed point that never ran.
     */
    int iterations = 0;
    /**
     * True on the limit path: a closed form has converged by construction. The
     * flag says whether the numbers can be trusted as the method's own answer,
     * not whether a loop ended, and a limit's answer is always its own.
     */
    bool converged = false;
};

namespace dispatch_detail {

/** The mean-field metrics, which are double by construction, in the caller's arithmetic. */
template <class T>
Matrix<T> env_widen(const Matrix<double>& A) {
    Matrix<T> B(A.rows(), A.cols(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) B(i, j) = num_traits<T>::from_double(A(i, j));
    return B;
}

/**
 * The limit path's cache surface, which is double by construction, in the
 * caller's arithmetic.
 *
 * `EnvLimitSolution` is deliberately not templated -- `SolverEnvLimit::init`
 * refuses `T != double` -- so the widening happens here rather than by
 * templating a struct that can only ever hold one arithmetic.
 */
template <class T>
solvers::CacheMetrics<T> env_widen_cache(const solvers::CacheMetrics<double>& c) {
    solvers::CacheMetrics<T> out;
    for (std::size_t i = 0; i < c.caches.size(); ++i) {
        const solvers::CacheNodeMetrics<double>& m = c.caches[i];
        solvers::CacheNodeMetrics<T> w;
        w.node = m.node;
        w.name = m.name;
        w.itemcap = m.itemcap;
        w.itemsize = m.itemsize;
        w.nitems = m.nitems;
        for (std::size_t k = 0; k < m.hitprob.size(); ++k)
            w.hitprob.push_back(num_traits<T>::from_double(m.hitprob[k]));
        for (std::size_t k = 0; k < m.missprob.size(); ++k)
            w.missprob.push_back(num_traits<T>::from_double(m.missprob[k]));
        for (std::size_t k = 0; k < m.delayedprob.size(); ++k)
            w.delayedprob.push_back(num_traits<T>::from_double(m.delayedprob[k]));
        for (std::size_t k = 0; k < m.latency.size(); ++k)
            w.latency.push_back(num_traits<T>::from_double(m.latency[k]));
        for (std::size_t k = 0; k < m.listcost.size(); ++k)
            w.listcost.push_back(num_traits<T>::from_double(m.listcost[k]));
        w.hitproblist = env_widen<T>(m.hitproblist);
        w.itemprob = env_widen<T>(m.itemprob);
        out.caches.push_back(w);
    }
    return out;
}

/**
 * `EnvOptions` as the state-vector coupling reads it.
 *
 * `stage_solver` is carried through UNTRANSLATED, so a caller who left it at
 * the mean-field default of `fluid` is refused by name inside
 * `SolverEnvStatevec::init` -- which is the reference's own error, "this method
 * requires SolverENV to be instantiated with the CTMC solver"
 * (`@@SolverENV/SolverENV.m` line 763). Silently substituting `ctmc` would run a
 * different inner solver from the one asked for.
 *
 * TWO THINGS `EnvOptions` CANNOT CARRY, and neither is invented here: the
 * per-transition state reset maps (functions on a state vector, which have no
 * counterpart among the mean-field resets, as `ResetStateVec` explains) and the
 * per-stage horizon overrides. A caller who needs either builds an
 * `EnvStatevecOptions` and calls `solver_env_statevec` directly. The horizon
 * that IS carried is `timespan_end`, and it should be a few mean holding times
 * rather than a large safe number, for the ode23 reason `EnvStatevecOptions`
 * gives at length.
 */
template <class T>
EnvStatevecOptions<T> env_statevec_options(const EnvOptions& o) {
    EnvStatevecOptions<T> s;
    s.iter_max = o.iter_max;
    s.iter_tol = o.iter_tol;
    s.stage_solver = o.stage_solver;
    s.sojourn = o.sojourn;
    s.timespan_end = o.timespan_end;
    // `stage_cutoff` bounds an OPEN stage's level count under either backend:
    // it truncates the enumerated chain for the CTMC one and the LD-QBD level
    // recursion for the MAM one. A closed stage takes its population instead
    // and ignores it in both.
    if (o.stage_cutoff > 0) s.mam_cutoff = static_cast<std::size_t>(o.stage_cutoff);
    return s;
}

/** Read the state-vector result into the common shape. */
template <class T>
void env_take_statevec(EnvAnalyzerSolution<T>& out) {
    out.method = "statevec";
    out.QN = out.statevec.QN;
    out.UN = out.statevec.UN;
    out.TN = out.statevec.TN;
    out.iterations = out.statevec.iterations;
    out.converged = out.statevec.converged;
    out.cache = out.statevec.cache;
}

/**
 * Read the mean-field result into the common shape.
 *
 * `asked` is reported rather than the bare analyzer name, because `smp` and
 * `statedep` also run this analyzer and are not the same run: `statedep`
 * rewrote the environment as it went, and a result labelled `meanfield` would
 * describe an environment that was an input when it was an output.
 */
template <class T>
void env_take_meanfield(EnvAnalyzerSolution<T>& out, const std::string& asked = "meanfield") {
    out.method = asked.empty() ? "meanfield" : asked;
    out.QN = env_widen<T>(out.meanfield.avg.QN);
    out.UN = env_widen<T>(out.meanfield.avg.UN);
    out.TN = env_widen<T>(out.meanfield.avg.TN);
    out.iterations = out.meanfield.avg.iterations;
    out.converged = out.meanfield.avg.converged;
    out.cache = out.meanfield.cache;
}

/** Read a closed-form limit result into the common shape. */
template <class T>
void env_take_limit(EnvAnalyzerSolution<T>& out) {
    out.method = out.limit.method;
    out.QN = env_widen<T>(out.limit.QN);
    out.UN = env_widen<T>(out.limit.UN);
    out.TN = env_widen<T>(out.limit.TN);
    out.iterations = 0;
    out.converged = true;
    out.cache = env_widen_cache<T>(out.limit.cache);
}

/** True for the one name the reference's `init` installs the state-vector analyzer on. */
inline bool env_is_statevec(const std::string& m) { return m == "statevec"; }

/** True for the two names `runAnalyzer` sends to `solveEnvLimit`. */
inline bool env_is_limit(const std::string& m) { return m == "avg" || m == "dec"; }

}  // namespace dispatch_detail

/**
 * Port of `SolverENV.listValidMethods`.
 *
 * Every name here selects a COUPLING -- what crosses an environment switch --
 * and each is dispatched by this file or by `SolverEnv`'s own ladder:
 * `default`/`meanfield`/`mean`/`blend` carry the marginal means, `meancov` the
 * same coupling carrying a COVARIANCE beside the mean (the exit second moment
 * over the sojourn law, mixed over the origins of a switch by the law of total
 * variance, handed to the next stage as `FluidOptions::init_qlen`/`init_qcov`),
 * `statevec` the whole joint distribution, `smp` the marginal means over the
 * semi-Markov stage probabilities of the embedded jump chain, `statedep` makes the transition depend on the
 * state it leaves, and `avg`/`dec` are the closed-form fast/slow limits of
 * `solveEnvLimit`.
 *
 * `default` is the reference's spelling and `meanfield` this port's; both name
 * the mean-field coupling and `SolverEnv::init` accepts either, as it does
 * `blend`.
 *
 * `blend` NAMES THE MEAN-FIELD COUPLING SINCE 2026-09-13 and named the
 * state-vector one before that. It is kept as a spelling of `default` so an
 * existing script still runs, but it no longer reaches `SolverEnvStatevec` and
 * the numbers it returns move: ask for `statevec` by name to carry the joint
 * distribution. The word survives inside the analyzers as the
 * environment-averaged BLEND they both compute, which is what it described.
 */
inline std::vector<std::string> env_list_valid_methods() {
    return {"default", "meanfield", "mean",     "meancov",  "blend",
            "blending", "smp",      "statedep", "statevec", "avg",   "dec"};
}

/** Port of `runAnalyzerChecks`' method gate: an unlisted method is refused. */
inline void env_check_method(const std::string& method) {
    const std::vector<std::string> valid = env_list_valid_methods();
    if (std::find(valid.begin(), valid.end(), method) != valid.end()) return;
    throw UnsupportedError("SolverENV: the '" + method + "' method is unsupported by this solver");
}

/**
 * `SolverENV.init`'s analyzer selection: solve the environment with the
 * coupling `o.method` names.
 *
 * An unknown method is refused by `SolverEnv`'s ladder, which is reached
 * because everything that is not the state-vector coupling IS the mean-field
 * one in the reference -- the `else` of its `if`, not a separate case.
 */
template <class T>
EnvAnalyzerSolution<T> solver_env(Environment<T>& e, const EnvOptions& o) {
    EnvAnalyzerSolution<T> out;
    if (dispatch_detail::env_is_limit(o.method)) {
        out.limit = solver_env_limit(e, o);
        dispatch_detail::env_take_limit(out);
        return out;
    }
    if (dispatch_detail::env_is_statevec(o.method)) {
        out.statevec = solver_env_statevec(e, dispatch_detail::env_statevec_options<T>(o));
        dispatch_detail::env_take_statevec(out);
        return out;
    }
    out.meanfield = solver_env_meanfield(e, o);
    dispatch_detail::env_take_meanfield(out, o.method);
    return out;
}

/**
 * The same, with the stage solver supplied by the caller.
 *
 * WHAT THIS IS FOR. The reference's `SolverENV(renv, @(m) SolverX(m, opts))`
 * takes a FACTORY and runs every stage with whatever it returns, and
 * `mapEnvApprox` hands it `feval(class(self), ...)` so that the environment
 * image is solved by the solver the user called. Without an injected callable
 * this port could only ever run fluid stages, and "the stage solver is the
 * calling solver" would simply be untrue here.
 *
 * WHY ONLY THE LIMITS TAKE IT, and why the other two refuse by name rather than
 * ignoring it. The two couplings carry an object across an environment switch
 * that only their own backend can produce -- the mean-field one the RMF cache
 * transient and the marginal means `initFromMarginal` reads, the state-vector
 * one a per-stage generator and the joint law over it -- so an arbitrary
 * `getAvg` callable cannot drive either. Accepting the function and quietly
 * solving with the built-in backend would answer a different question from the
 * one asked, which is the silent substitution this port refuses everywhere.
 *
 * The limits have no such object: `solveEnvLimit` asks each stage for its
 * steady state and blends, so any solver that can answer `getAvg` will do.
 */
template <class T>
EnvAnalyzerSolution<T> solver_env(Environment<T>& e, const EnvOptions& o,
                                  EnvStageAvgFn<T> stage_fn) {
    if (!stage_fn) return solver_env(e, o);
    if (!dispatch_detail::env_is_limit(o.method))
        throw UnsupportedError(
            "SolverENV: the '" + o.method +
            "' coupling cannot run an injected stage solver; it carries a per-stage object across "
            "an environment switch that only its own backend produces (the mean-field coupling "
            "the RMF transient and the marginal means, the state-vector one a generator and the "
            "joint law over it). The closed-form limits 'avg' and 'dec' take one");
    EnvAnalyzerSolution<T> out;
    SolverEnvLimit<T> s(e, o, stage_fn);
    out.limit = s.solve();
    dispatch_detail::env_take_limit(out);
    return out;
}

/**
 * The same, on a COMPRESSED environment: aggregate the stages first, then run
 * the coupling over the macro-states.
 *
 * The mean-field path delegates to `solver_env_meanfield`, which already owns
 * the construct-then-apply-then-solve order that `env_apply_macro_probabilities`
 * requires. The state-vector path repeats that order here rather than reaching
 * for a wrapper of its own: `SolverEnvStatevec`'s constructor calls
 * `Environment::init()`, which recomputes probEnv and probOrig from the macro
 * arcs, so the macro probabilities have to be written back AFTER it and BEFORE
 * `solve()`, exactly as on the mean-field side.
 */
template <class T>
EnvAnalyzerSolution<T> solver_env(Environment<T>& e, const EnvOptions& o,
                                  const EnvCompressOptions& c) {
    EnvAnalyzerSolution<T> out;
    out.compressed = true;
    // A limit reads probEnv and the stage models, both of which the aggregation
    // rewrites, so it runs on the macro-states like the couplings do; the macro
    // probabilities must still be written back after `SolverEnvLimit`'s own
    // `Environment::init()` and before it reads them, which is why the construct
    // and the solve are separated here.
    if (dispatch_detail::env_is_limit(o.method)) {
        out.compression = env_compress(e, c);
        SolverEnvLimit<T> s(*out.compression.env, o);
        env_apply_macro_probabilities(*out.compression.env, out.compression);
        out.limit = s.solve();
        dispatch_detail::env_take_limit(out);
        return out;
    }
    if (dispatch_detail::env_is_statevec(o.method)) {
        out.compression = env_compress(e, c);
        SolverEnvStatevec<T> s(*out.compression.env, dispatch_detail::env_statevec_options<T>(o));
        env_apply_macro_probabilities(*out.compression.env, out.compression);
        out.statevec = s.solve();
        dispatch_detail::env_take_statevec(out);
        return out;
    }
    out.meanfield = solver_env_meanfield(e, o, c);
    out.compression = out.meanfield.compression;
    dispatch_detail::env_take_meanfield(out, o.method);
    return out;
}

}  // namespace env
}  // namespace line

#endif  // LINE_SOLVERS_ENV_ENV_DISPATCH_H
