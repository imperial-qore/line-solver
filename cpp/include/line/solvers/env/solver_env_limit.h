/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_ENV_SOLVER_ENV_LIMIT_H
#define LINE_SOLVERS_ENV_SOLVER_ENV_LIMIT_H

/**
 * The two CLOSED-FORM environment limits: `SolverENV.solveEnvLimit`, reached by
 * `options.method` in {`avg`, `dec`}.
 *
 * Both replace the fixed point with a single reading, and each is exact at one
 * end of the time-scale separation between the environment and the network:
 *
 *   `avg`, the FAST-environment limit. The environment switches so much faster
 *   than the network responds that the network only ever sees the AVERAGE of the
 *   modulated rates. One rate-averaged model is built -- every station-class
 *   rate that varies across stages replaced by its probEnv-weighted mean, as an
 *   exponential -- and solved once, in steady state. Exact as the stage-switch
 *   rate goes to infinity.
 *
 *   `dec`, the SLOW-environment limit, i.e. quasi-stationary decomposition. The
 *   environment stays in a stage long enough for the network to reach that
 *   stage's own steady state, so each stage is solved independently and the
 *   metrics are averaged with weights probEnv. Exact as the stage-switch rate
 *   goes to zero.
 *
 * NEITHER CARRIES ANYTHING ACROSS A SWITCH, which is what makes them closed
 * form and also what they give up: no entry state, no reset policy, no
 * transient. A reset policy declared on an arc is therefore inert here, exactly
 * as it is in the reference, whose `solveEnvLimit` never reads `resetFun`.
 *
 * WHAT IS AVERAGED, AND WHAT IS LEFT ALONE (`buildRateAveragedModel`). Only a
 * Source, a Queue or a Delay carries a rate to average. A station-class pair
 * that is disabled in any stage, or whose rate does not actually vary across
 * them, keeps its ORIGINAL distribution rather than being rewritten as an
 * exponential of its own mean -- so a non-modulated Erlang stays an Erlang, and
 * the base model is preserved exactly outside the modulated rates.
 *
 * WHAT IS REFUSED. A stage holding a Cache: the reference aggregates the hit and
 * miss ratios over the stages (`SolverENV.accumCacheMetric`) and writes them
 * back onto the stage-1 model, and this port has no cache metric to aggregate,
 * so it would report Q, U and T while silently dropping the answer a cache model
 * is asked for.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <type_traits>
#include <vector>

#include "line/lang/qn/environment.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/solvers/env/solver_env.h"
#include "line/solvers/fluid/solver_fluid.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace env {

/** What a limit solve reports. */
struct EnvLimitSolution {
    /** Environment-averaged metrics, (nstations x nclasses). */
    Matrix<double> QN, UN, TN;
    /** The limit that ran: `avg` or `dec`. */
    std::string method;
    /** `dec` only: each stage's own steady state, before the probEnv blend. */
    std::vector<Matrix<double> > QStage, UStage, TStage;
    /** The stage probabilities the blend used. */
    std::vector<double> prob_env;
};

/**
 * The limit solver.
 *
 * `envObj` must already carry a model per stage; `init()` is called here, as
 * `solveEnvLimit` calls `self.init()` before reading probEnv.
 */
template <class T>
class SolverEnvLimit {
public:
    SolverEnvLimit(Environment<T>& e, const EnvOptions& o) : envObj(e), opt(o) { init(); }

    EnvLimitSolution solve() { return dec_ ? solve_dec() : solve_avg(); }

private:
    void init() {
        if (opt.method != "avg" && opt.method != "dec")
            throw UnsupportedError("SolverENV limit: '" + opt.method +
                                   "' is not a closed-form environment limit; the two are 'avg' "
                                   "(fast environment) and 'dec' (slow environment)");
        dec_ = (opt.method == "dec");
        if (opt.stage_solver != "fluid")
            throw UnsupportedError(
                "SolverENV limit: stage solver '" + opt.stage_solver +
                "' is not available; the limits solve each stage in STEADY STATE and the fluid "
                "analyzer is the one this port wires into the environment");
        if (!std::is_same<T, double>::value)
            throw UnsupportedError(
                "SolverENV limit: a fluid stage integrates its drift with LSODA, which is double "
                "only; rerun with --arith double");

        // The two limits read the stage NETWORKS directly -- `dec` solves each
        // in steady state, `avg` builds one network at the probEnv-weighted
        // rates -- and a layered stage has no station rate table to weight. The
        // reference has no branch for it either: `solveEnvLimit` reaches for
        // `self.ensemble{1}.nodes`, which a LayeredNetwork does not carry.
        envObj.reject_lqn_stages(
            "SolverENV limit",
            "the fast/slow limits read a stage's station rates directly -- 'dec' solves each "
            "stage network in steady state and 'avg' builds one network at the probEnv-weighted "
            "rates -- and a layered model has no such rate table, only the layers SolverLN "
            "derives from it");
        envObj.init();
        const std::size_t E = envObj.nstages();
        M = envObj.stage(0).model.nstations;
        K = envObj.stage(0).model.nclasses;
        for (std::size_t e = 1; e < E; ++e)
            if (envObj.stage(e).model.nstations != M || envObj.stage(e).model.nclasses != K)
                throw InputError(
                    "SolverENV limit: every stage must have the same stations and classes; the "
                    "metrics are blended entrywise across them");
        for (std::size_t e = 0; e < E; ++e) {
            const qn::NetworkStruct<T>& sn = envObj.stage(e).model;
            for (std::size_t ind = 0; ind < sn.nodes.size(); ++ind)
                if (sn.nodes[ind].nodetype == lang::NodeType::Cache)
                    throw UnsupportedError(
                        "SolverENV limit: stage " + std::to_string(e + 1) +
                        " holds a Cache, whose environment-blended hit and miss ratios come from "
                        "SolverENV.accumCacheMetric over the per-stage cache results; no cache "
                        "metric is reported by the stage solver here, so the blend would drop "
                        "the answer the model is asked for");
        }
    }

    /** The slow-environment limit: each stage on its own, blended by probEnv. */
    EnvLimitSolution solve_dec() {
        const std::size_t E = envObj.nstages();
        EnvLimitSolution out;
        out.method = "dec";
        out.prob_env = envObj.prob_env;
        out.QN = Matrix<double>(M, K, 0.0);
        out.UN = Matrix<double>(M, K, 0.0);
        out.TN = Matrix<double>(M, K, 0.0);
        out.QStage.assign(E, Matrix<double>(M, K, 0.0));
        out.UStage.assign(E, Matrix<double>(M, K, 0.0));
        out.TStage.assign(E, Matrix<double>(M, K, 0.0));
        for (std::size_t e = 0; e < E; ++e) {
            const fluid::FluidSolution s = stage_steady_state(envObj.stage(e).model);
            const double p = envObj.prob_env[e];
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t r = 0; r < K; ++r) {
                    out.QStage[e](i, r) = s.QN(i, r);
                    out.UStage[e](i, r) = s.UN(i, r);
                    out.TStage[e](i, r) = s.TN(i, r);
                    out.QN(i, r) += p * s.QN(i, r);
                    out.UN(i, r) += p * s.UN(i, r);
                    out.TN(i, r) += p * s.TN(i, r);
                }
        }
        return out;
    }

    /** The fast-environment limit: one rate-averaged model, solved once. */
    EnvLimitSolution solve_avg() {
        EnvLimitSolution out;
        out.method = "avg";
        out.prob_env = envObj.prob_env;
        const qn::NetworkStruct<T> avg = rate_averaged_model();
        const fluid::FluidSolution s = stage_steady_state(avg);
        out.QN = Matrix<double>(M, K, 0.0);
        out.UN = Matrix<double>(M, K, 0.0);
        out.TN = Matrix<double>(M, K, 0.0);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r) {
                out.QN(i, r) = s.QN(i, r);
                out.UN(i, r) = s.UN(i, r);
                out.TN(i, r) = s.TN(i, r);
            }
        return out;
    }

    fluid::FluidSolution stage_steady_state(const qn::NetworkStruct<T>& sn) const {
        fluid::FluidOptions fo = opt.stage;
        // A limit starts from nothing carried over -- there is no entry state to
        // seed with, which is the whole content of the approximation.
        fo.init_sol.clear();
        return fluid::solver_fluid(sn, fo);
    }

    /**
     * `buildRateAveragedModel`: stage 1's network with every MODULATED rate
     * replaced by its probEnv-weighted mean, as an exponential.
     *
     * The refresh chain is rerun because `rates` and everything derived from it
     * are read off the service table; editing the table alone would leave the
     * struct describing stage 1.
     */
    qn::NetworkStruct<T> rate_averaged_model() const {
        const std::size_t E = envObj.nstages();
        qn::NetworkStruct<T> sn = envObj.stage(0).model;
        std::vector<double> r(E, 0.0);
        for (std::size_t i = 0; i < M; ++i) {
            const lang::NodeType nt = sn.stations[i].nodetype;
            if (nt != lang::NodeType::Source && nt != lang::NodeType::Queue &&
                nt != lang::NodeType::Delay)
                continue;
            for (std::size_t k = 0; k < K; ++k) {
                bool ok = true;
                for (std::size_t e = 0; e < E && ok; ++e) {
                    const qn::NetworkStruct<T>& se = envObj.stage(e).model;
                    // The reference's `any(isnan(r)) || any(r<=0)`: a rate this
                    // port reports as disabled is MATLAB's NaN, and either way
                    // the pair is left as configured.
                    if (se.disabled[i][k]) {
                        ok = false;
                        break;
                    }
                    r[e] = num_traits<T>::to_double(se.rates(i, k));
                    if (!(r[e] > 0.0) || !std::isfinite(r[e])) ok = false;
                }
                if (!ok) continue;
                double lo = r[0], hi = r[0], avg = 0.0;
                for (std::size_t e = 0; e < E; ++e) {
                    lo = std::min(lo, r[e]);
                    hi = std::max(hi, r[e]);
                    avg += envObj.prob_env[e] * r[e];
                }
                // Not modulated: keep the original distribution, which may carry
                // a shape the exponential of its mean would throw away.
                if (hi - lo <= 1e-12 * std::max(1.0, hi)) continue;
                sn.set_service(i + 1, k + 1,
                               lang::Distrib<T>::exp_rate(num_traits<T>::from_double(avg)));
            }
        }
        sn.refresh_struct();
        return sn;
    }

    Environment<T>& envObj;
    EnvOptions opt;
    std::size_t M = 0, K = 0;
    bool dec_ = false;
};

/** `solveEnvLimit` on the original stages. */
template <class T>
EnvLimitSolution solver_env_limit(Environment<T>& e, const EnvOptions& o) {
    SolverEnvLimit<T> s(e, o);
    return s.solve();
}

}  // namespace env
}  // namespace line

#endif  // LINE_SOLVERS_ENV_SOLVER_ENV_LIMIT_H
