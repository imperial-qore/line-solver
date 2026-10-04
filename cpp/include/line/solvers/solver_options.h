/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_SOLVER_OPTIONS_H
#define LINE_SOLVERS_SOLVER_OPTIONS_H

/**
 * @file
 * @ingroup line_solvers
 * @ingroup line_public
 * The keyword arguments a solver takes, as one struct.
 *
 * Python builds a solver as `MVA(model, method='lin')` or
 * `CTMC(model, cutoff=3, seed=23000)` -- one shared keyword surface across
 * every solver, narrowed per solver by which keywords it reads. This is that
 * surface: a field per keyword, and a fluent `set_*` per field so a C++ call
 * reads like the Python one it was ported from.
 *
 *     SolverOptions().set_cutoff(3).set_seed(23000)
 *
 * A field left at its default is NOT forwarded, so the engine's own default
 * stands: that is what the negative/empty sentinels below mean. This mirrors
 * MATLAB's `Solver.defaultOptions` merge rather than overwriting it.
 */

#include <cstddef>
#include <string>
#include <vector>

namespace line {

/** The knobs a solver reads; a negative or empty field keeps the engine default. */
struct SolverOptions {
    std::string method = "default";
    double tol = -1.0;
    double iter_tol = -1.0;
    int iter_max = -1;
    std::size_t samples = 0;
    unsigned long seed = 0;
    /** CTMC `options.cutoff`, scalar or per class. */
    double cutoff = -1.0;
    std::vector<std::size_t> cutoff_vec;
    /** CTMC refusal threshold on the state-space size. */
    std::size_t state_max = 0;
    /** Fluid `options.timespan(2)` and warm start. */
    double timespan_end = -1.0;
    std::vector<double> init_sol;
    bool stiff = false;
    /** The `options.config` fields of the MVA / NC / fluid families. */
    std::string multiserver;
    std::string highvar;
    std::string np_priority;
    /** MVA / NC `options.config.fork_join`: 'default'/'mmt'/'fjt' or 'ht'. */
    std::string fork_join;
    /** SSA `options.config.state_space_gen`. */
    std::string state_space_gen;
    /**
     * `options.config.symbolic`: the computer-algebra backend to resolve.
     *
     * Empty keeps the engine default, which is `auto`. Other values are a URL,
     * an image name, or `none` to stay local. Hoisted out of `config` like
     * every other field here: this struct is flat by design.
     */
    std::string symbolic;
    /**
     * `options.config.symbolic_timeout`, seconds per symbolic request.
     *
     * Non-positive keeps the engine default of 300. It is a real knob and not
     * a formality: a symbolic solve grows superpolynomially in the state count,
     * and the fluid equilibria of a smoothed softmin drift are transcendental
     * in exp(-pstar*x) and outrun 300 s on models that are otherwise ordinary.
     */
    int symbolic_timeout = -1;
    /**
     * `options.config.map_env`: `"default"` lets a model whose ONLY unsupported
     * features are non-renewal processes fall back to its random-environment
     * image, `"off"` keeps the plain rejection.
     *
     * CALLER-FACING ONLY, and deliberately not on `MvaOptions`, `NcSolverOptions`
     * or `FluidOptions`: no `runAnalyzer` reads it in the reference either. The
     * fallback happens ABOVE the runner, in `getAvg`'s place, because the runner
     * must go on refusing the model it cannot consume.
     */
    std::string map_env;
    /** `options.config.map_env_method`: `auto`, `meanfield`, `dec` or `avg`. */
    std::string map_env_method;
    /** `options.config.map_env_maxstages`: 0 keeps `map2renv`'s own cap of 64. */
    std::size_t map_env_maxstages = 0;
    /** SolverBA `options.level`. */
    int level = 2;
    /**
     * JMT `options.config.replications`: the independent JSIM runs the transient ensemble of
     * `default` averages over a finite timespan. Non-positive keeps the engine default of 10.
     */
    int replications = -1;
    /** JMT `options.keep`: leave the scratch directory in place after the solve. */
    bool keep = false;
    bool verbose = false;

    SolverOptions& set_method(const std::string& v) { method = v; return *this; }
    SolverOptions& set_tol(double v) { tol = v; return *this; }
    SolverOptions& set_iter_tol(double v) { iter_tol = v; return *this; }
    SolverOptions& set_iter_max(int v) { iter_max = v; return *this; }
    SolverOptions& set_replications(int v) { replications = v; return *this; }
    SolverOptions& set_samples(std::size_t v) { samples = v; return *this; }
    SolverOptions& set_seed(unsigned long v) { seed = v; return *this; }
    SolverOptions& set_cutoff(double v) { cutoff = v; return *this; }
    SolverOptions& set_cutoff_vec(const std::vector<std::size_t>& v) { cutoff_vec = v; return *this; }
    SolverOptions& set_state_max(std::size_t v) { state_max = v; return *this; }
    SolverOptions& set_timespan_end(double v) { timespan_end = v; return *this; }
    SolverOptions& set_init_sol(const std::vector<double>& v) { init_sol = v; return *this; }
    SolverOptions& set_stiff(bool v) { stiff = v; return *this; }
    SolverOptions& set_multiserver(const std::string& v) { multiserver = v; return *this; }
    SolverOptions& set_highvar(const std::string& v) { highvar = v; return *this; }
    SolverOptions& set_np_priority(const std::string& v) { np_priority = v; return *this; }
    SolverOptions& set_fork_join(const std::string& v) { fork_join = v; return *this; }
    SolverOptions& set_state_space_gen(const std::string& v) { state_space_gen = v; return *this; }
    SolverOptions& set_symbolic(const std::string& v) { symbolic = v; return *this; }
    SolverOptions& set_symbolic_timeout(int v) { symbolic_timeout = v; return *this; }
    SolverOptions& set_map_env(const std::string& v) { map_env = v; return *this; }
    SolverOptions& set_map_env_method(const std::string& v) { map_env_method = v; return *this; }
    SolverOptions& set_map_env_maxstages(std::size_t v) { map_env_maxstages = v; return *this; }
    SolverOptions& set_level(int v) { level = v; return *this; }
    SolverOptions& set_keep(bool v) { keep = v; return *this; }
    SolverOptions& set_verbose(bool v) { verbose = v; return *this; }
};

}  // namespace line

#endif  // LINE_SOLVERS_SOLVER_OPTIONS_H
