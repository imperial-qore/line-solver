/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MAP_ENV_STAGES_H
#define LINE_SOLVERS_MAP_ENV_STAGES_H

/**
 * @file
 * @ingroup line_solvers
 * The stage solvers `map_env_approx` injects, one per runner that can be a
 * caller.
 *
 * WHAT THEY ARE. The C++ spelling of `feval(class(self), stageModel,
 * self.options)` in `mapEnvApprox.m`: the environment image is solved with the
 * solver that was asked for, at the settings it was asked for. Binding the
 * caller's own options into the callable is the whole point -- a hard-coded
 * backend silently discards the method, the tolerances and the iteration caps
 * the user set, which is the defect the injection exists to close.
 *
 * WHY A SHARED HEADER RATHER THAN THREE LAMBDAS PER CALLER. Two translation
 * units enter the fallback, `solver_facade.cpp` (`avg_table`) and
 * `line_cli.cpp` (`run_avg_engine` plus the per-solver `-a avg` arms), and they
 * must agree on what a stage solve reports or the same model answers differently
 * through the facade and through the CLI. That is the divergence the Knobs
 * struct exists to prevent one layer up, and it applies here for the same
 * reason.
 *
 * They are `double`-only because the driver is: `map_env_approx` and
 * `SolverEnvLimit::init` both refuse `T != double`, so a multiprecision caller
 * keeps the plain refusal rather than falling back.
 */

#include "line/lang/qn/network_struct.h"
#include "line/solvers/cache_metrics.h"
#include "line/solvers/env/solver_env_limit.h"
#include "line/solvers/fluid/fluid_runner.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/solvers/nc/solver_nc_runner.h"
#include "line/util/matrix.h"

namespace line {
namespace solvers {

/** MVA stages, bound to the caller's own `MvaOptions`. */
inline env::EnvStageAvgFn<double> mva_stage_fn(const mva::MvaOptions& opt) {
    return [opt](const qn::NetworkStruct<double>& stage) {
        Matrix<double> init;
        const mva::AvgResult<double> r = mva::solver_mva_run_analyzer(stage, opt, init);
        env::EnvStageAvg<double> s;
        s.QN = r.QN;
        s.UN = r.UN;
        s.TN = r.TN;
        s.cache = r.cache;
        return s;
    };
}

/** NC stages, bound to the caller's own `NcSolverOptions`. */
inline env::EnvStageAvgFn<double> nc_stage_fn(const nc::NcSolverOptions& opt) {
    return [opt](const qn::NetworkStruct<double>& stage) {
        const mva::AvgResult<double> r = nc::solver_nc_run_analyzer(stage, opt);
        env::EnvStageAvg<double> s;
        s.QN = r.QN;
        s.UN = r.UN;
        s.TN = r.TN;
        s.cache = r.cache;
        return s;
    };
}

/** Fluid stages, bound to the caller's own `FluidOptions`. */
inline env::EnvStageAvgFn<double> fluid_stage_fn(const fluid::FluidOptions& opt) {
    return [opt](const qn::NetworkStruct<double>& stage) {
        // The five-argument overload is the only one that reports cache metrics;
        // the short one drops them, which on a cache-bearing stage is the whole
        // quantity the environment blend is trying to recombine.
        solvers::CacheMetrics<double> cache;
        const fluid::FluidSolution r =
            fluid::solver_fluid_run_analyzer<double>(stage, opt, nullptr, nullptr, &cache);
        env::EnvStageAvg<double> s;
        s.QN = r.QN;
        s.UN = r.UN;
        s.TN = r.TN;
        s.cache = cache;
        return s;
    };
}

}  // namespace solvers
}  // namespace line

#endif  // LINE_SOLVERS_MAP_ENV_STAGES_H
