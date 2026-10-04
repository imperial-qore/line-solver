/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_ENV_SOLVER_ENV_H
#define LINE_SOLVERS_ENV_SOLVER_ENV_H

/**
 * @file
 * @ingroup line_solvers
 * SolverENV: a queueing network in a random environment.
 *
 * Port of `matlab/src/solvers/@@SolverENV/SolverENV.m` driving
 * `solver_env_meanfield_analyzer.m`, the default (and historical) coupling.
 *
 * THE FIXED POINT, in one paragraph. Each stage is solved TRANSIENTLY, not in
 * steady state, because a stage does not last long enough to reach one -- the
 * environment switches first. What the network carries across a switch is its
 * queue lengths, so stage e must be started from the mean queue lengths it is
 * ENTERED with, and those are the exit queue lengths of whichever stage
 * preceded it. That is a fixed point over the entry vectors, and the iteration
 * below is a Picard iteration on it:
 *
 *     Qentry[e] = sum_h probOrig(h, e) * reset_{h->e}( Qexit[h][e] )
 *
 * where Qexit[h][e] is the transient of stage h AVERAGED OVER WHEN the h -> e
 * switch happens -- the increments of that transition's CDF are the weights.
 * The environment-averaged answer at the end weights each stage's transient
 * over its own HOLDING time CDF instead (any exit, not a particular one), and
 * blends the stages by their stationary probabilities probEnv.
 *
 * WHY THE WEIGHTS ARE CDF INCREMENTS AND NOT A DENSITY. The transient comes
 * back on a grid, so the reference integrates the metric against the measure
 * the transition induces on that same grid: w_j = F(t_j) - F(t_{j-1}), with
 * w_1 = 0, normalized by their sum. That is a Riemann-Stieltjes sum, and it
 * needs no density -- which matters, because a deterministic transition has
 * none.
 *
 * THE MEAN-FIELD COLLAPSE is the approximation, and it is worth naming: only
 * the MARGINAL MEAN queue lengths cross a switch, so any correlation between
 * stations at the moment of the switch is discarded.
 * `solver_env_statevec_analyzer.m` carries the whole joint distribution
 * instead; it is ported in `solver_env_statevec.h`, as a separate class with a
 * CTMC stage solver, and `env_dispatch.h` chooses between the two the way
 * `SolverENV.init` does. `method = "statevec"` is refused HERE by name rather
 * than silently served by the mean-field coupling.
 *
 * A LAYERED STAGE IS THE SAME FIXED POINT WITH A DIFFERENT STAGE SOLVER. When
 * a stage holds a `lqn::LqnStruct` rather than a flat network, its transient
 * comes from a `SolverLN` over the stage's own layers instead of from one fluid
 * integration, and the (station, class) view the coupling blends over is the
 * BLOCK-DIAGONAL UNION of those layers -- `LayeredNetwork.layerBlocks` in the
 * reference, `SolverLN::layer_blocks` here. Nothing else changes: the same
 * marginal mean queue lengths cross a switch, split into per-layer blocks on
 * the way in (`SolverLN::init_from_marginal`) and reassembled on the way out.
 * The one thing to know about the handoff is that it must be replayed AFTER the
 * layered fixed point and BEFORE the layer transients, because the fixed point
 * resets its layers as it converges; `run_layer_transient` is where that
 * happens, and a warm start installed anywhere earlier is silently inert.
 *
 * WHAT IS PORTED, and what is refused by name:
 *   ported   `meanfield` (the default), `smp` (the mean-field analyzer over the
 *            semi-Markov stage probabilities -- see `smp_stage_probabilities`), `statedep` (state-dependent environment rates,
 *            `resetEnvRates`), stochastic and deterministic sojourn,
 *            per-transition reset policies including the named `keep`/`clear`
 *            of a node breakdown, fluid stages, CTMC stages, and LayeredNetwork
 *            stages over fluid layers
 *   refused  `statevec`/`blend` (this class is the mean-field coupling; use
 *            env::solver_env for them), `avg`/`dec` (the closed-form fast/slow
 *            limits of `solveEnvLimit`, ported as env::SolverEnvLimit and
 *            likewise reached through env::solver_env), and the cache
 *            aggregation of `aggregateCacheMeanfield_`.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "line/api/mam/map_cdf.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mc/dtmc_solve.h"
#include "line/lang/qn/environment.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/ctmc/solver_ctmc_transient.h"
#include "line/solvers/fluid/fluid_dae.h"
#include "line/solvers/fluid/fluid_kp.h"
#include "line/solvers/fluid/solver_fluid.h"
#include "line/solvers/ln/solver_ln.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace env {

/**
 * `LnOptions` as a LAYERED environment stage is solved with.
 *
 * Only the layer engine differs from the SolverLN default, and it differs
 * because the mean-field coupling has no use for a stage it cannot integrate:
 * see `EnvOptions::lqn`.
 */
inline ln::LnOptions env_default_lqn_options() {
    ln::LnOptions o;
    o.layer_solver = "fluid";
    return o;
}

/** Options of SolverENV. Defaults are `Solver.defaultOptions`. */
struct EnvOptions {
    int iter_max = 100;
    double iter_tol = 1e-4;
    /**
     * Which solver runs each FLAT stage: the fluid transient or the enumerated
     * CTMC. A LAYERED stage is run by SolverLN whatever this says, and its own
     * layer engine is named by `lqn.layer_solver`; asking for `ctmc` alongside a
     * layered stage is refused rather than reinterpreted, since an LQN has no
     * single generator to enumerate.
     */
    std::string stage_solver = "fluid";
    /** The inter-stage coupling: `meanfield` is the reference's default. */
    std::string method = "meanfield";
    /** `options.sojourn`: `stochastic` (default) or `deterministic`. */
    std::string sojourn = "stochastic";
    /** Options handed to each stage solver. */
    fluid::FluidOptions stage;
    /**
     * `options.cutoff` of a CTMC stage, read only when `stage_solver` is ctmc.
     * An open stage needs one; a closed one enumerates its own population.
     */
    double stage_cutoff = -1.0;
    /** `options.timespan(2)` of the inner solver: the transient horizon. */
    double timespan_end = 100.0;
    /**
     * Points on a UNIFORM transient grid, used only where `stage_grid` declines
     * to build one (a stage whose holding time is the disabled 1 x 1 zero pair).
     *
     * It used to be the accuracy knob of the whole method, because the exit
     * metrics are a Riemann-Stieltjes sum over whatever grid the stage was
     * integrated on and a uniform grid resolves the HORIZON rather than the
     * sojourn. Since 2026-08-11 each stage is integrated on the grid the
     * sojourn asks for -- 90% of the points under `5*E[S]` -- so the answer no
     * longer moves with this: on renv_node_breakdown the horizon may be 100 or
     * 1000 and Server QLen is 0.458854 or 0.459191.
     */
    std::size_t tran_points = 1001;
    /**
     * Options of the `SolverLN` that runs a LAYERED stage.
     *
     * The reference names the stage solver by handing SolverENV a FACTORY --
     * `ENV(env, @(m) LN(m, @(mm) FLD(mm), 'timespan', [0 T]))` -- and this field
     * is that factory's argument list, because a C++ template cannot take a
     * MATLAB function handle. `timespan_end`, `tran_points` and `tran_grid` are
     * OVERWRITTEN by the environment for every stage: the horizon and the grid
     * the exit average is summed on belong to the coupling, not to a stage.
     *
     * `layer_solver` defaults to `fluid` here where `LnOptions` alone defaults
     * to `mva`, and the difference is not a preference: the coupling carries
     * queue lengths across a switch and therefore needs each stage's
     * TRANSIENT, which only the fluid layer engine produces. A caller who names
     * another engine is refused by `init`, not quietly served the fluid one.
     */
    ln::LnOptions lqn = env_default_lqn_options();
};

/** What SolverENV reports. */
struct EnvSolution {
    /** Environment-averaged metrics, (nstations x nclasses). */
    Matrix<double> QN, UN, TN;
    /** Per-stage, sojourn-averaged metrics. */
    std::vector<Matrix<double>> QExit, UExit, TExit;
    /** The entry queue lengths the fixed point converged to. */
    std::vector<Matrix<double>> Qentry;
    /**
     * `meancov` only: the environment-wide queue-length COVARIANCE over the
     * (station, class) pairs, (M*K)-by-(M*K) and indexed `ir = r*M + i`, and its
     * diagonal reshaped to (nstations x nclasses). Empty under every other
     * coupling, which carries a first moment alone.
     */
    Matrix<double> QCov, QVar;
    int iterations = 0;
    bool converged = false;
};

namespace detail {

/**
 * `maxpe(approx, exact)`: max |1 - approx/exact| over the entries where exact
 * is nonzero, which is what the reference's convergence test compares.
 */
inline double env_maxpe(const std::vector<double>& approx, const std::vector<double>& exact) {
    double worst = -1.0;
    for (std::size_t i = 0; i < approx.size() && i < exact.size(); ++i) {
        if (exact[i] == 0.0) continue;
        const double e = std::fabs(1.0 - approx[i] / exact[i]);
        if (e > worst) worst = e;
    }
    return worst;  // negative means "no comparable entry", the reference's empty
}

/** Linear interpolation of a transient metric at `d`, clamped to the grid. */
inline double env_det_eval(const std::vector<double>& t, const std::vector<double>& metric,
                           double d) {
    if (t.empty()) return 0.0;
    d = std::max(t.front(), std::min(d, t.back()));
    for (std::size_t j = 1; j < t.size(); ++j) {
        if (d <= t[j]) {
            const double dt = t[j] - t[j - 1];
            if (!(dt > 0.0)) return metric[j];
            const double a = (d - t[j - 1]) / dt;
            return metric[j - 1] * (1.0 - a) + metric[j] * a;
        }
    }
    return metric.back();
}

}  // namespace detail

/**
 * The environment solver.
 *
 * `envObj` must already carry a model per stage; `init()` is called here, as
 * `SolverENV.init` does.
 */
template <class T>
class SolverEnv {
public:
    SolverEnv(Environment<T>& e, const EnvOptions& o) : envObj(e), opt(o) { init(); }

    EnvSolution solve() {
        const std::size_t E = envObj.nstages();
        EnvSolution out;
        out.QN = Matrix<double>(M, K, 0.0);
        out.UN = Matrix<double>(M, K, 0.0);
        out.TN = Matrix<double>(M, K, 0.0);

        pre();
        std::vector<std::vector<double>> qfirst_prev(E), qfirst_curr(E);
        int it = 0;
        for (it = 1; it <= opt.iter_max; ++it) {
            for (std::size_t e = 0; e < E; ++e) analyze(e);
            // The reference's convergence test compares the FIRST point of the
            // transient -- the entry queue length -- across iterations.
            qfirst_prev = qfirst_curr;
            qfirst_curr.assign(E, std::vector<double>());
            for (std::size_t e = 0; e < E; ++e) {
                qfirst_curr[e].reserve(M * K);
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t r = 0; r < K; ++r) qfirst_curr[e].push_back(tranQ[e][i][r][0]);
            }
            bool conv = it > 1;
            if (conv)
                for (std::size_t e = 0; e < E && conv; ++e) {
                    const double d = detail::env_maxpe(qfirst_curr[e], qfirst_prev[e]);
                    if (d < 0.0) continue;  // nothing comparable, treated as converged
                    if (!std::isfinite(d) || d >= opt.iter_tol) conv = false;
                }
            post();
            if (conv) {
                out.converged = true;
                break;
            }
        }
        out.iterations = std::min(it, opt.iter_max);
        finish(out);
        return out;
    }

private:
    void init() {
        if (opt.stage_solver != "fluid" && opt.stage_solver != "ctmc")
            throw UnsupportedError(
                "SolverENV: stage solver '" + opt.stage_solver +
                "' is not available; the environment coupling needs a TRANSIENT stage solve and "
                "only the fluid analyzer and the enumerated CTMC provide one in this port");
        ctmc_stages = (opt.stage_solver == "ctmc");
        if (opt.method == "statevec")
            throw UnsupportedError(
                "SolverENV: the state-vector analyzer (solver_env_statevec_analyzer.m) carries the "
                "full joint distribution across a switch, which THIS class does not: it is the "
                "mean-field coupling and carries the marginal means. The state-vector coupling is "
                "ported as env::SolverEnvStatevec, with its own options and its CTMC stage solver; "
                "reach it through env::solver_env (env_dispatch.h), which is the analyzer "
                "selection of SolverENV.init");
        if (opt.method == "avg" || opt.method == "dec")
            throw UnsupportedError(
                "SolverENV: the closed-form fast/slow environment limits (SolverENV.solveEnvLimit) "
                "replace this fixed point rather than configure it -- neither carries anything "
                "across a switch. They are ported as env::SolverEnvLimit, with their own result "
                "type; reach them through env::solver_env (env_dispatch.h), which is the method "
                "dispatch of SolverENV.runAnalyzer");
        // `smp` selects no analyzer of its own: it runs the mean-field one after
        // `smp_stage_probabilities` recomputes prob_env and prob_orig from the
        // embedded jump chain of the arc CDFs, as the JAR's SolverENV.init does.
        statedep = (opt.method == "statedep");
        // "default" IS the reference's spelling of this coupling: SolverENV.m
        // installs solver_env_meanfield_analyzer for every method that is not
        // statevec, and its listValidMethods leads with 'default'. Only the CLI
        // normalised it away, so a caller constructing SolverEnv directly with
        // the reference's own name was refused. "blend" is a third spelling of
        // it since 2026-09-13, when that name was repointed here from the
        // state-vector coupling, and this whitelist is why it has to be named:
        // falling past the statevec branch reaches the mean-field analyzer, but
        // an unlisted name is refused before it gets there.
        if (opt.method != "meanfield" && opt.method != "default" && opt.method != "mean"
            && opt.method != "meancov" && opt.method != "blend" && opt.method != "blending"
            && opt.method != "smp"
            && !statedep)
            throw UnsupportedError("SolverENV: unknown method '" + opt.method + "'");
        // 'meancov' IS THIS COUPLING, carrying a COVARIANCE beside the mean: the
        // propagation is the mean-field one and only the object crossing a switch
        // is richer, so it is a flag here rather than an analyzer of its own.
        meancov = (opt.method == "meancov");
        // THE STAGE METHOD IS THE CALLER'S, as it is in the reference, where a
        // stage is an ordinary SolverFluid the caller constructed: `kp` integrates
        // the Ko-Pender fluid AND diffusion limits, and is the only route in this
        // port that carries a second moment along a stage trajectory. Every other
        // fluid method integrates the closing drift, which carries a first moment
        // alone; `meancov` over such stages keeps the TIMING variance of the
        // sojourn and the spread of the per-origin means, and nothing of the
        // within-stage spread, exactly as the reference does for a stage solver
        // that reports none.
        {
            std::string sm = opt.stage.method;
            if (sm.size() > 4 && sm.compare(0, 4, "fld.") == 0) sm = sm.substr(4);
            kp_stages = (!ctmc_stages && sm == "kp");
            // `dae` is the OTHER route that carries a second moment along a stage
            // trajectory: it integrates the linear-noise covariance beside the
            // min-normal mean, and unlike `kp` its DRIFT reads that covariance, so
            // a stage restarting at Sigma(0) = 0 gets the mean wrong as well as
            // the variance. It needs a stage route of its own: `solver_fluid_transient`
            // integrates the CLOSING drift whatever the method name says (only
            // `solver_fluid_run_transient` dispatches, and it takes no output grid),
            // so a `dae` stage sent there would be answered by a different method,
            // report no covariance at all, and drop `init_qcov` in silence.
            dae_stages = (!ctmc_stages && !kp_stages && sm == "dae");
        }
        // A CTMC STAGE IS NOT REFUSED under `meancov`: it takes a lattice STATE
        // rather than a distribution, so the carried covariance simply does not
        // seed it, which is what the reference does too. The coupling then keeps
        // the timing variance of the sojourn and the spread of the per-origin
        // means, and reports both.
        if (!(opt.timespan_end > 0.0) || !std::isfinite(opt.timespan_end))
            throw InputError(
                "SolverENV: the stage transient needs a finite positive horizon, "
                "options.timespan(2); the mean-field coupling integrates each stage's "
                "TRAJECTORY against the holding-time CDF and has nothing to integrate over an "
                "infinite one");
        if (opt.tran_points < 2)
            throw InputError("SolverENV: the transient grid needs at least two points");
        if (!std::is_same<T, double>::value)
            throw UnsupportedError(
                "SolverENV: a fluid stage integrates its drift with LSODA, which is double only; "
                "rerun with --arith double");

        envObj.init();
        if (opt.method == "smp") smp_stage_probabilities();
        const std::size_t E = envObj.nstages();
        build_lqn_stages();
        M = stage_M(0);
        K = stage_K(0);
        for (std::size_t e = 1; e < E; ++e)
            if (stage_M(e) != M || stage_K(e) != K)
                throw InputError(
                    "SolverENV: every stage must have the same stations and classes; the metrics "
                    "are blended entrywise across them. A LAYERED stage counts the block-diagonal "
                    "union of the layers SolverLN builds for it, which is what layerBlocks "
                    "reports, so two layered stages agree exactly when they have the same layer "
                    "shape -- not merely the same LQN element count");
        entry.assign(E, Matrix<double>(M, K, 0.0));
        centry.assign(E, Matrix<double>(M * K, M * K, 0.0));
        tranC.assign(E, std::vector<Matrix<double>>());
        tranT.assign(E, std::vector<double>());
        tranQ.assign(E, std::vector<std::vector<std::vector<double>>>());
        tranU.assign(E, std::vector<std::vector<std::vector<double>>>());
        tranTp.assign(E, std::vector<std::vector<std::vector<double>>>());
        steady.assign(E, StageSteady());
        steady_asked.assign(E, 0);
        det_sojourn = (opt.sojourn == "deterministic");
        dvals.assign(E, 0.0);
        refresh_sojourn_means();
        if (statedep) {
            bool any = false;
            for (std::size_t e = 0; e < E; ++e)
                for (std::size_t h = 0; h < E; ++h)
                    if (envObj.arc(e, h).enabled && envObj.arc(e, h).reset_rates) any = true;
            if (!any)
                throw InputError(
                    "SolverENV: method 'statedep' updates each environment transition from the "
                    "state its stage is left in, and no arc carries a rate reset "
                    "(Environment::set_env_rate_reset, resetEnvRatesFun in the reference); the "
                    "run would be the mean-field fixed point under a different name");
        }
    }

    // ---- layered (LayeredNetwork) stages ---------------------------------

    /**
     * Build one persistent `SolverLN` per LAYERED stage.
     *
     * ONE SOLVER FOR THE WHOLE RUN, not one per sweep, and that is the
     * reference's arrangement too (`self.solvers{e}` outlives the iteration).
     * It matters twice over: the layer ensemble is what carries the aggregate
     * (station, class) SHAPE this coupling blends over, and it has to be known
     * before any stage is solved; and the layered fixed point is solved ONCE,
     * because a steady solve ignores the state it is started from -- only the
     * transient reads it -- so re-running it every environment sweep would
     * recompute the same numbers.
     */
    void build_lqn_stages() {
        const std::size_t E = envObj.nstages();
        lnsolv.assign(E, std::shared_ptr<ln::SolverLN<T>>());
        lnblk.assign(E, ln::LnLayerBlocks());
        if (!envObj.has_lqn_stages()) return;
        if (ctmc_stages)
            throw UnsupportedError(
                "SolverENV: a LayeredNetwork stage is solved by SolverLN over its layers, not by "
                "an enumerated CTMC over one stage generator -- an LQN has no single generator to "
                "enumerate. Leave stage_solver at 'fluid', which names the engine each LAYER is "
                "integrated with (EnvOptions::lqn.layer_solver)");
        if (opt.lqn.layer_solver != "fluid")
            throw UnsupportedError(
                "SolverENV: a LayeredNetwork stage is carried across an environment switch by its "
                "queue lengths, so the coupling needs the stage's TRANSIENT; among the layer "
                "engines only the fluid one produces one, and this ensemble asks for '" +
                opt.lqn.layer_solver + "' layers (EnvOptions::lqn.layer_solver)");
        for (std::size_t e = 0; e < E; ++e) {
            if (!envObj.is_lqn(e)) continue;
            ln::LnOptions lo = opt.lqn;
            // THE HORIZON AND THE GRID BELONG TO THE COUPLING. A stage lasts as
            // long as the environment lets it, and the exit average is summed
            // against the holding-time CDF, so neither is a property of the LQN
            // and a caller's value for either would silently change what the
            // Stieltjes sum integrates over. `tran_grid` is installed per sweep
            // by `lqn_stage_transient`, since a state-dependent method may move
            // the holding time between sweeps.
            lo.timespan_end = opt.timespan_end;
            lo.tran_points = opt.tran_points;
            lo.tran_grid.clear();
            lnsolv[e] = std::make_shared<ln::SolverLN<T>>(envObj.stage(e).lqn_model, lo);
            lnblk[e] = lnsolv[e]->layer_blocks();
        }
    }

    /** Stations of stage `e`: its own, or the block-diagonal union of its layers. */
    std::size_t stage_M(std::size_t e) const {
        return envObj.is_lqn(e) ? lnblk[e].M : envObj.stage(e).model.nstations;
    }

    /** Classes of stage `e`: its own, or the block-diagonal union of its layers. */
    std::size_t stage_K(std::size_t e) const {
        return envObj.is_lqn(e) ? lnblk[e].K : envObj.stage(e).model.nclasses;
    }

    /**
     * One LAYERED stage's transient, assembled block-diagonally into the same
     * `tranT`/`tranQ`/`tranU`/`tranTp` the flat path fills.
     *
     * `seed` says whether the fixed point has an entry vector yet, exactly as it
     * does on the CTMC path: `pre()` runs before it has one, and warm-starting
     * from an all-zero matrix would empty every closed layer rather than leave
     * it on its default state.
     *
     * THE OFF-DIAGONAL BLOCKS ARE ZERO AND ARE ALLOCATED ANYWAY. They pair a
     * station of one layer with a class of another and stand for nothing, which
     * is why the reference leaves them as empty cells; here they must still be
     * present and sized, because the convergence test reads the FIRST point of
     * every (station, class) series and an absent series has no first point.
     * Zero series contribute nothing to the blend and are skipped by `maxpe`,
     * which ignores an exact value of zero -- so the two ports agree on the
     * numbers as well as on the shape.
     */
    void lqn_stage_transient(std::size_t e, bool seed) {
        ln::SolverLN<T>& s = *lnsolv[e];
        s.init_from_marginal(seed ? entry[e] : Matrix<double>());
        s.set_tran_grid(stage_grid(e));
        const ln::LnTranSolution tr = s.get_tran_avg();
        const ln::LnLayerBlocks& b = lnblk[e];

        // Every layer was integrated on the same grid -- one horizon, one
        // out_grid -- so layer one's is the stage's time base.
        tranT[e].clear();
        if (!tr.layers.empty()) tranT[e] = tr.layers[0].t;
        const std::size_t P = std::max<std::size_t>(tranT[e].size(), 1);
        tranQ[e].assign(M, std::vector<std::vector<double>>(K, std::vector<double>(P, 0.0)));
        tranU[e].assign(M, std::vector<std::vector<double>>(K, std::vector<double>(P, 0.0)));
        tranTp[e].assign(M, std::vector<std::vector<double>>(K, std::vector<double>(P, 0.0)));
        if (tranT[e].empty()) return;

        // `b` and the blocks were both derived from the SAME layer ensemble --
        // it is built once, in SolverLN's constructor, and never resized -- so
        // `msz`/`ksz` are the exact extents of each block and the offsets tile
        // the aggregate exactly. Nothing here is a bounds guard.
        for (std::size_t l = 0; l < tr.layers.size(); ++l) {
            const ln::LnTranLayer& L = tr.layers[l];
            for (std::size_t i = 0; i < b.msz[l]; ++i)
                for (std::size_t r = 0; r < b.ksz[l]; ++r) {
                    const std::size_t row = b.roff[l] + i, col = b.coff[l] + r;
                    copy_series(L.QN[i][r], tranQ[e][row][col]);
                    copy_series(L.UN[i][r], tranU[e][row][col]);
                    copy_series(L.TN[i][r], tranTp[e][row][col]);
                }
        }
    }

    /**
     * Copy a layer series onto the stage grid.
     *
     * The two lengths agree by construction -- every layer was integrated on the
     * grid this stage's time base came from -- and the loop is written over
     * whichever is shorter so that a layer whose integration was cut short
     * leaves the tail at zero rather than reading past its own trajectory.
     */
    static void copy_series(const std::vector<double>& src, std::vector<double>& dst) {
        for (std::size_t j = 0; j < dst.size() && j < src.size(); ++j) dst[j] = src[j];
    }

    /** F_eh(t) of one arc, `evalCDF` in the reference; zero for a disabled arc. */
    double arc_cdf(std::size_t e, std::size_t h, double t) const {
        const EnvArc<T>& a = envObj.arc(e, h);
        if (!a.enabled) return 0.0;
        return num_traits<T>::to_double(lang::dist_cdf(a.dist, num_traits<T>::from_double(t)));
    }

    /**
     * Semi-Markov stage probabilities for method `smp`, as `SolverENV.init` does in the JAR and
     * MATLAB. P(k,e) = int dF_ke(t) prod_{h!=k,e} (1 - F_kh(t)) is the embedded jump chain,
     * integrated on N = max(1000, 100T) intervals with the survival at the midpoint, T doubling
     * from 1 until F_ke(T) >= 1 - 1e-8. The mean holding time integrates the sojourn survival
     * prod_{h!=k} (1 - F_kh(t)) by composite Simpson (N = 10000) over [0, U], U doubling from 10
     * until the survival is <= 1e-8. prob_env(k) is then proportional to pie_dtmc(k) * hold(k),
     * and prob_orig(k,e) = prob_env(k) E0(k,e) / sum_{h!=e} prob_env(h) E0(h,e), E0 the arc rates.
     */
    void smp_stage_probabilities() {
        const std::size_t E = envObj.nstages();
        const double eps = 1e-8;
        Matrix<double> E0(E, E, 0.0);
        for (std::size_t k = 0; k < E; ++k)
            for (std::size_t h = 0; h < E; ++h) {
                const EnvArc<T>& a = envObj.arc(k, h);
                if (!a.enabled) continue;
                const double m = num_traits<T>::to_double(lang::dist_moment(a.dist, 1));
                E0(k, h) = (m > 0.0) ? 1.0 / m : 0.0;
            }
        Matrix<double> P(E, E, 0.0);
        for (std::size_t k = 0; k < E; ++k)
            for (std::size_t e = 0; e < E; ++e) {
                if (k == e || !envObj.arc(k, e).enabled) continue;
                double Tup = 1.0;
                while (arc_cdf(k, e, Tup) < 1.0 - eps) {
                    Tup *= 2.0;
                    if (Tup > 1e6) break;
                }
                const std::size_t N = std::max<std::size_t>(1000, static_cast<std::size_t>(std::lround(Tup * 100.0)));
                const double dt = Tup / static_cast<double>(N);
                double sum = 0.0, Fprev = arc_cdf(k, e, 0.0);
                for (std::size_t i = 0; i < N; ++i) {
                    const double t1 = static_cast<double>(i + 1) * dt;
                    const double Fnext = arc_cdf(k, e, t1);
                    const double tmid = t1 - 0.5 * dt;
                    double surv = 1.0;
                    for (std::size_t h = 0; h < E; ++h)
                        if (h != k && h != e && envObj.arc(k, h).enabled) surv *= 1.0 - arc_cdf(k, h, tmid);
                    sum += (Fnext - Fprev) * surv;
                    Fprev = Fnext;
                }
                P(k, e) = sum;
            }
        const std::vector<double> pie = mc::dtmc_solve(P);

        std::vector<double> hold(E, 0.0);
        const std::size_t Nh = 10000;
        for (std::size_t k = 0; k < E; ++k) {
            auto surv = [&](double t) {
                double s = 1.0;
                for (std::size_t h = 0; h < E; ++h)
                    if (h != k) s *= 1.0 - arc_cdf(k, h, t);
                return s;
            };
            double U = 10.0;
            while (surv(U) > eps) {
                U *= 2.0;
                if (U > 1e6) break;
            }
            const double dt = U / static_cast<double>(Nh);
            double integral = 0.0;
            for (std::size_t i = 0; i < Nh; ++i) {
                const double t0 = static_cast<double>(i) * dt, t1 = t0 + dt;
                integral += (surv(t0) + 4.0 * surv(0.5 * (t0 + t1)) + surv(t1)) * dt / 6.0;
            }
            hold[k] = integral;
        }

        double denom = 0.0;
        for (std::size_t e = 0; e < E; ++e) denom += pie[e] * hold[e];
        std::vector<double> pi(E, 0.0);
        for (std::size_t k = 0; k < E; ++k) pi[k] = pie[k] * hold[k] / denom;
        envObj.prob_env = pi;
        Matrix<double> emb(E, E, 0.0);
        for (std::size_t e = 0; e < E; ++e) {
            double s = 0.0;
            for (std::size_t h = 0; h < E; ++h)
                if (h != e) s += pi[h] * E0(h, e);
            if (s > 0.0)
                for (std::size_t k = 0; k < E; ++k)
                    if (k != e) emb(k, e) = pi[k] * E0(k, e) / s;
        }
        envObj.prob_orig = emb;
    }

    /** The mean holding times the deterministic sojourn evaluates at. */
    void refresh_sojourn_means() {
        if (!det_sojourn) return;
        for (std::size_t e = 0; e < envObj.nstages(); ++e)
            dvals[e] = mam::map_mean(envObj.hold_time[e].map());
    }

    /**
     * `pre_`: seed every stage from its own solve at the first iteration, so
     * the fixed point starts somewhere the model could be.
     *
     * Which solve is the reference's branch on the STAGE horizon: a finite one
     * takes the LAST POINT of the stage transient, and only an infinite one --
     * which this coupling refuses in `init`, since it has no transient to
     * integrate -- takes the steady state. The distinction is not cosmetic on a
     * stage that is unstable on its own -- a broken server offered more than it
     * can serve -- where the steady state is whatever the integrator's own
     * horizon happened to reach and is far from the transient at
     * `timespan_end`. The seed only moves by the drift per environment cycle,
     * so a seed off by a factor of ten costs iterations proportionally.
     */
    void pre() {
        const std::size_t E = envObj.nstages();
        for (std::size_t e = 0; e < E; ++e) {
            if (envObj.is_lqn(e)) {
                // The layered seed is the reference's own: the LAST POINT of the
                // stage transient run from the layers' default state, which for
                // a finite horizon is what `pre_` takes for every stage solver.
                lqn_stage_transient(e, false);
                if (tranT[e].empty()) continue;
                const std::size_t last = tranT[e].size() - 1;
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t r = 0; r < K; ++r) entry[e](i, r) = tranQ[e][i][r][last];
                continue;
            }
            if (ctmc_stages) {
                // A CTMC stage has no "seed from nowhere" solve to run: it
                // starts from the model's OWN initial state, which is what
                // `initDefault` put there, so the first sweep enters at that
                // state's queue lengths rather than at zero.
                const ctmc::CtmcTransient<T> tr = ctmc_stage_transient(e, false, std::vector<double>());
                if (tr.t.empty()) continue;
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t r = 0; r < K; ++r)
                        entry[e](i, r) = num_traits<T>::to_double(tr.QNt[i][r].back());
                continue;
            }
            if (kp_stages) {
                kp_stage_transient(e, false);
                if (tranT[e].empty()) continue;
                const std::size_t last = tranT[e].size() - 1;
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t r = 0; r < K; ++r) entry[e](i, r) = tranQ[e][i][r][last];
                continue;
            }
            fluid::FluidOptions fo = opt.stage;
            fo.init_sol.clear();
            // THE SAME METHOD AS `analyze`: the entry vector this seeds the fixed
            // point with is read off a stage trajectory, so taking it from the
            // closing drift and then iterating with `dae` would start the loop at
            // another method's answer.
            const std::vector<fluid::FluidTranPoint> tr =
                dae_stages ? fluid::solver_fluid_dae_transient(envObj.stage(e).model, fo,
                                                               opt.timespan_end, opt.tran_points)
                           : fluid::solver_fluid_transient(envObj.stage(e).model, fo,
                                                           opt.timespan_end, opt.tran_points);
            if (tr.empty()) continue;
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t r = 0; r < K; ++r) entry[e](i, r) = tr.back().QN(i, r);
        }
    }

    /** `analyze_`: the transient of one stage from its entry queue lengths. */
    void analyze(std::size_t e) {
        if (envObj.is_lqn(e)) {
            lqn_stage_transient(e, true);
            return;
        }
        tranT[e].clear();
        tranC[e].clear();
        tranQ[e].assign(M, std::vector<std::vector<double>>(K));
        tranU[e].assign(M, std::vector<std::vector<double>>(K));
        tranTp[e].assign(M, std::vector<std::vector<double>>(K));
        if (ctmc_stages) {
            // THE REFINED GRID, and it was MEASURED against the alternative.
            // The reference forms its Stieltjes sum on whatever points ode23
            // stopped at, so reporting on those looks like the faithful choice;
            // on renv_threestages_repairmen it is the worse one -- Queue1 QLen
            // 0.8344 against the reference's 0.83053, where the refined grid
            // gives 0.83092. The refined grid puts 90% of its points under
            // 5*E[S], which is where the holding-time CDF puts its mass, so it
            // resolves the integrand rather than the trajectory. The residual
            // 4.7e-4 is what is left of this quadrature difference.
            const ctmc::CtmcTransient<T> tr = ctmc_stage_transient(e, true, stage_grid(e));
            for (std::size_t j = 0; j < tr.t.size(); ++j) {
                tranT[e].push_back(num_traits<T>::to_double(tr.t[j]));
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t r = 0; r < K; ++r) {
                        tranQ[e][i][r].push_back(num_traits<T>::to_double(tr.QNt[i][r][j]));
                        tranU[e][i][r].push_back(num_traits<T>::to_double(tr.UNt[i][r][j]));
                        tranTp[e][i][r].push_back(num_traits<T>::to_double(tr.TNt[i][r][j]));
                    }
            }
            return;
        }
        const qn::NetworkStruct<T>& sn = envObj.stage(e).model;
        if (kp_stages) {
            kp_stage_transient(e, true);
            return;
        }
        fluid::FluidOptions fo = opt.stage;
        fo.init_sol = initsol_from_marginal(sn, entry[e]);
        // A `dae` stage takes the entry covariance as well as the entry mean, and
        // gives its own back. Only that method reads `init_qcov`: offering it to a
        // closing stage would be dropped in silence, which is the failure this
        // predicate exists to avoid.
        if (meancov && dae_stages) fo.init_qcov = centry[e];
        const std::vector<fluid::FluidTranPoint> tr =
            dae_stages ? fluid::solver_fluid_dae_transient(sn, fo, opt.timespan_end,
                                                           opt.tran_points, stage_grid(e))
                       : fluid::solver_fluid_transient(sn, fo, opt.timespan_end,
                                                       opt.tran_points, stage_grid(e));
        for (const fluid::FluidTranPoint& p : tr) {
            tranT[e].push_back(p.t);
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t r = 0; r < K; ++r) {
                    tranQ[e][i][r].push_back(p.QN(i, r));
                    tranU[e][i][r].push_back(p.UN(i, r));
                    tranTp[e][i][r].push_back(p.TN(i, r));
                }
            if (meancov && dae_stages && p.QCov.rows() == M * K) tranC[e].push_back(p.QCov);
        }
        // ALL OR NOTHING: `exit_cov` reads tranC[e] only when it has one matrix per
        // point, so a run where some points carried a covariance and others did not
        // must contribute none rather than a ragged series.
        if (tranC[e].size() != tranT[e].size()) tranC[e].clear();
    }

    /**
     * One stage's transient under the KO-PENDER limits, mean AND covariance.
     *
     * THE ONLY STAGE ROUTE THAT CARRIES A SECOND MOMENT in this port: the
     * closing drift integrates a first moment alone, so a `meancov` run over
     * closing stages keeps the TIMING variance of the sojourn and nothing of the
     * within-stage spread. `kp` integrates both on one grid, which is why the
     * covariance is read off the same solve as the mean rather than asked for
     * again. Selected by `EnvOptions::stage.method`, so a caller asks for the
     * limit by name exactly as a standalone fluid solve does.
     *
     * `seed` says whether the fixed point has an entry vector yet, as on the
     * layered and CTMC paths: `pre()` runs before it has one, and seeding from
     * an all-zero mean and covariance would start every stage empty rather than
     * at this method's own stationary arrival phase.
     */
    void kp_stage_transient(std::size_t e, bool seed) {
        tranT[e].clear();
        tranC[e].clear();
        const qn::NetworkStruct<T>& sn = envObj.stage(e).model;
        fluid::FluidOptions fo = opt.stage;
        fo.method = "kp";
        fo.init_sol.clear();
        fo.timespan_end = opt.timespan_end;
        if (seed) {
            fo.init_qlen = entry[e];
            if (meancov) fo.init_qcov = centry[e];
        }
        fluid::FluidKpTransient tr;
        fluid::solver_fluid_kp_core(sn, fo, &tr, stage_grid(e));
        tranT[e] = tr.t;
        const std::size_t P = tr.t.size();
        tranQ[e].assign(M, std::vector<std::vector<double>>(K, std::vector<double>(P, 0.0)));
        tranU[e].assign(M, std::vector<std::vector<double>>(K, std::vector<double>(P, 0.0)));
        tranTp[e].assign(M, std::vector<std::vector<double>>(K, std::vector<double>(P, 0.0)));
        // Mean and covariance come off the SAME integration, on the same grid, so
        // the two moments this coupling carries need no interpolation onto one
        // another and the covariance costs nothing beyond the aggregation.
        for (std::size_t n = 0; n < P; ++n)
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t r = 0; r < K; ++r) {
                    tranQ[e][i][r][n] = tr.QN[n](i, r);
                    tranU[e][i][r][n] = tr.UN[n](i, r);
                    tranTp[e][i][r][n] = tr.TN[n](i, r);
                }
        if (meancov) tranC[e] = tr.QCov;
    }

    /**
     * One stage's transient under an ENUMERATED CTMC, from its entry marginal.
     *
     * THE MARGINAL IS ROUNDED HERE and is not for a fluid stage, which is the
     * reference's own split (`roundMarginalForDiscreteSolver` runs for every
     * stage solver except SolverFluid): a chain has no state holding 1.7 jobs,
     * so the fixed point's continuous iterate has to be placed on the lattice
     * before it can be a state at all. The rounding conserves each closed
     * class's population, so the placed state is one the chain actually
     * contains.
     *
     * The state is written onto a COPY of the stage struct as its one-row
     * declared state space, which is the same channel `setState` reaches
     * `solver_ctmc_transient_analyzer` through -- so this seeds the stage
     * exactly as a caller would, rather than through a private door.
     */
    ctmc::CtmcTransient<T> ctmc_stage_transient(std::size_t e, bool seed,
                                                const std::vector<double>& grid) const {
        qn::NetworkStruct<T> sn = envObj.stage(e).model;
        // `seed` is explicit and not inferred from the arguments: `pre()` runs
        // BEFORE the fixed point has an entry vector, and its all-zero matrix
        // would place every job nowhere -- a state the chain does not have.
        if (seed) seed_ctmc_state(sn, entry[e]);
        std::vector<T> g;
        g.reserve(grid.size());
        for (std::size_t j = 0; j < grid.size(); ++j) g.push_back(num_traits<T>::from_double(grid[j]));
        return ctmc::solver_ctmc_transient_analyzer(
            sn, ctmc_stage_options(), num_traits<T>::from_int(0),
            num_traits<T>::from_double(opt.timespan_end), g);
    }

    /** The CTMC options a stage is solved with; the cutoff is the caller's. */
    ctmc::CtmcOptions ctmc_stage_options() const {
        ctmc::CtmcOptions co;
        co.cutoff = opt.stage_cutoff;
        return co;
    }

    /**
     * Write the rounded entry marginal onto `sn` as its declared state.
     *
     * A row that the encoding cannot realize leaves the station alone, so the
     * analyzer falls back to that station's default marking rather than being
     * handed a state the chain does not contain -- an error the coupling could
     * not act on, where the default is at least a state of the same model.
     */
    void seed_ctmc_state(qn::NetworkStruct<T>& sn, const Matrix<double>& Q) const {
        if (Q.empty()) return;
        const std::vector<std::vector<std::size_t>> nir = round_marginal(sn, Q);
        for (std::size_t i = 0; i < sn.nstations && i < nir.size(); ++i) {
            const std::size_t ind = sn.station_to_node[i];  // 1-based node index
            // A SOURCE'S STATE IS NOT A FUNCTION OF ITS QUEUE LENGTH: `State.fromMarginal`'s EXT
            // branch ignores n and seeds every generated class in phase one, which is the default
            // row. Built from the zero marginal it is all-zero, a state the chain does not contain.
            if (sn.stations[i].nodetype == lang::NodeType::Source) continue;
            std::vector<std::size_t> ph(sn.nclasses, 1);
            for (std::size_t r = 0; r < sn.nclasses; ++r) ph[r] = sn.phasessz_of(i + 1, r + 1);
            std::vector<T> row;
            if (!qn::from_marginal_node_first(sn, ind, nir[i], ph, row))
                continue;
            Matrix<T> space(1, row.size());
            for (std::size_t c = 0; c < row.size(); ++c) space(0, c) = row[c];
            sn.statespace[ind] = space;
            sn.stateprior[ind] = std::vector<T>(1, num_traits<T>::from_int(1));
        }
    }

    /**
     * `roundMarginalForDiscreteSolver`: the continuous per-class marginal on the
     * integer lattice, conserving each CLOSED class's population.
     *
     * Rounding entrywise does not conserve it -- three stations holding 0.5 jobs
     * of a one-job class round to 0, 0, 0 and the class disappears from the
     * chain -- so the residual is handed to the stations with the largest
     * fractional parts, which is the reference's own rule.
     */
    std::vector<std::vector<std::size_t>> round_marginal(const qn::NetworkStruct<T>& sn,
                                                         const Matrix<double>& Q) const {
        std::vector<std::vector<std::size_t>> n(sn.nstations, std::vector<std::size_t>(sn.nclasses, 0));
        for (std::size_t r = 0; r < sn.nclasses; ++r) {
            const double pop = sn.classes[r].population;
            std::vector<double> frac(sn.nstations, 0.0);
            long placed = 0;
            for (std::size_t i = 0; i < sn.nstations; ++i) {
                const double q = (i < Q.rows() && r < Q.cols()) ? std::max(0.0, Q(i, r)) : 0.0;
                const double fl = std::floor(q);
                n[i][r] = static_cast<std::size_t>(fl);
                frac[i] = q - fl;
                placed += static_cast<long>(fl);
            }
            if (!std::isfinite(pop)) continue;  // an open class has no population to conserve
            long want = static_cast<long>(std::llround(pop));
            while (placed < want) {
                std::size_t best = 0;
                double bv = -1.0;
                for (std::size_t i = 0; i < sn.nstations; ++i)
                    if (frac[i] > bv) { bv = frac[i]; best = i; }
                ++n[best][r];
                frac[best] = -1.0;
                ++placed;
                if (bv < 0.0) break;  // nowhere left to place; the loop would not terminate
            }
            while (placed > want) {
                std::size_t best = 0;
                double bv = 2.0;
                bool any = false;
                for (std::size_t i = 0; i < sn.nstations; ++i)
                    if (n[i][r] > 0 && frac[i] < bv) { bv = frac[i]; best = i; any = true; }
                if (!any) break;
                --n[best][r];
                frac[best] = 2.0;
                --placed;
            }
        }
        return n;
    }

    /** One stage's steady tables, in the (station, class) space of the blend. */
    struct StageSteady {
        Matrix<double> QN, UN, TN;
    };

    /**
     * The stage's own STEADY tables, which are the exit value of every metric
     * its transient carries no trajectory for.
     *
     * A METRIC WITH NO TRAJECTORY IS NOT A METRIC WORTH ZERO. A stage solver
     * advances the quantities its integration actually carries; what it leaves
     * out it holds CONSTANT over the stage, so the value at the exit instant is
     * that constant whatever the sojourn was. `post()` and `finish()` used to
     * skip such a stage outright and blend the zeros, reporting a stage whose
     * metrics never move as ABSENT rather than as constant. That is the defect
     * `stageSteady_` closes in the reference (`solver_env_meanfield_analyzer.m`)
     * and `stageSteady`/`_stage_steady` close in the JAR and native python.
     *
     * THE HORIZON IS PUT ASIDE FOR THE ASK, and here it is not a refusal to
     * dodge but a DIFFERENT ANSWER to avoid. The reference's `getAvg` rejects a
     * finite `options.timespan` by name; `solver_fluid` instead ACCEPTS one and
     * quietly stops its window iteration there, disarming both the fixed-point
     * and the geometric-tail tests. Asked with the stage horizon still on the
     * options it would hand back the state at `t = timespan_end` -- the very
     * transient this seed exists to complement -- instead of the constant.
     *
     * THE ASK GOES TO THE ENGINE THAT PRODUCED THE TRAJECTORY, branch for branch
     * with `analyze`: the closing drift's own fixed point for a closing stage
     * (`solver_fluid_transient` integrates that drift whatever the method name
     * says), `kp` and `dae` through their own steady entries, and the enumerated
     * chain's stationary law for a CTMC stage. A seed taken from another closure
     * would not be the constant the trajectory holds.
     *
     * A LAYERED STAGE IS NOT ASKED, exactly as the JAR's `stageSteadyAsk`
     * refuses a non-`Network` stage by name: `SolverLN`'s steady table is over
     * the LQN's OWN nodes -- hosts, tasks, entries, activities -- and not over
     * the block-diagonal (station, class) aggregate this coupling blends in, so
     * a leading-block copy would read LQN node k as class k. Its transient is
     * already assembled in the right layout by `lqn_stage_transient`, and every
     * cell it owns carries one.
     *
     * ASKED ONCE PER STAGE, AND ONLY WHERE IT IS READ. The stage networks do not
     * change across sweeps -- `statedep` rewrites the environment arcs, not the
     * models -- so the constants are the same every time; and a stage whose
     * trajectory IS readable has every cell of the seed overwritten from it, so
     * asking then would cost a steady solve per sweep to be discarded. That is
     * this port's form of the reference's "read `result` first".
     */
    const StageSteady& stage_steady(std::size_t e) {
        if (steady_asked[e]) return steady[e];
        steady_asked[e] = 1;
        StageSteady& out = steady[e];
        out.QN = Matrix<double>(M, K, 0.0);
        out.UN = Matrix<double>(M, K, 0.0);
        out.TN = Matrix<double>(M, K, 0.0);
        if (envObj.is_lqn(e)) return out;
        const qn::NetworkStruct<T>& sn = envObj.stage(e).model;
        if (ctmc_stages) {
            const ctmc::CtmcAvg<T> a = ctmc::solver_ctmc_analyzer(sn, ctmc_stage_options()).avg;
            copy_finite(a.QN, out.QN);
            copy_finite(a.UN, out.UN);
            copy_finite(a.TN, out.TN);
            return out;
        }
        fluid::FluidOptions fo = opt.stage;
        // NOTHING CARRIED IN: the constant is what the stage holds on its own,
        // and the entry vector the handoff computes is the transient's business.
        fo.init_sol.clear();
        fo.init_qlen = Matrix<double>();
        fo.init_qcov = Matrix<double>();
        fo.timespan_end = std::numeric_limits<double>::infinity();
        fluid::FluidSolution s;
        if (kp_stages) {
            fo.method = "kp";
            s = fluid::solver_fluid_kp<T>(sn, fo);
        } else if (dae_stages) {
            fo.method = "dae";
            s = fluid::solver_fluid_dae<T>(sn, fo);
        } else {
            fo.method = "closing";
            s = fluid::solver_fluid<T>(sn, fo);
        }
        copy_finite(s.QN, out.QN);
        copy_finite(s.UN, out.UN);
        copy_finite(s.TN, out.TN);
        return out;
    }

    /**
     * Copy `src` onto `dst`, leaving every non-finite entry at zero: a NaN is
     * ABSENT, not a number, and must not enter the probEnv blend.
     *
     * The two shapes agree by construction -- `init` refuses an environment
     * whose stages differ in stations or classes, and a flat stage's own tables
     * are (nstations, nclasses) -- so this is a copy and not a fit.
     */
    template <class S>
    static void copy_finite(const Matrix<S>& src, Matrix<double>& dst) {
        for (std::size_t i = 0; i < dst.rows(); ++i)
            for (std::size_t r = 0; r < dst.cols(); ++r) {
                const double v = num_traits<S>::to_double(src(i, r));
                if (std::isfinite(v)) dst(i, r) = v;
            }
    }

    /**
     * The exit tables of stage `e`: its transient averaged over the sojourn on
     * the cells that carry one, and its steady constants on the cells that do
     * not. `want_all` asks for U and T beside Q, which only a state-dependent
     * rate reads in `post()` and which `finish()` always reports.
     *
     * THE WEIGHT IS THE STAGE SOJOURN, NOT THE e -> h CLOCK, and the reference
     * says so where it builds it: competing exponentials leave the exit TIME
     * independent of which destination won, so the exit average does not depend
     * on h at all -- h enters only through the reset applied on the way in.
     * Weighting by `proc[e][h]` instead read the transient over the mean of ONE
     * risk rather than of their minimum: on renv_twostages_repairmen, whose
     * Stage2 competes a 0.5 self arc with a 0.5 arc back, that is a mean of 2
     * against the sojourn's 1, and it reported Queue1 QLen 0.55882 against
     * MATLAB's 0.55550 with the two stations' throughputs 1.4% apart in a closed
     * cycle that admits one throughput.
     *
     * A TRAJECTORY THE COUPLING CANNOT AVERAGE IS NOT AN ANSWER OF ZERO. Two
     * cases reach `stage_steady` rather than a blend: a stage that reported no
     * transient at all, and one whose sojourn puts NO MASS on the grid, which is
     * this port's form of the reference's single-point transient (`tR = 1`, the
     * table `cdf_weights` can form no increment over). Both mean the same thing
     * -- nothing moved that the exit instant could be read off -- and the stage's
     * own constant is the exit value in both.
     */
    void stage_exit(std::size_t e, bool want_all, Matrix<double>& Qe, Matrix<double>& Ue,
                    Matrix<double>& Te) {
        std::vector<double> w;
        double wsum = 0.0;
        if (!tranT[e].empty() && !det_sojourn) {
            w = cdf_weights(envObj.hold_time[e].map(), tranT[e]);
            for (double v : w) wsum += v;
        }
        if (tranT[e].empty() || (!det_sojourn && !(wsum > 0.0))) {
            const StageSteady& st = stage_steady(e);
            Qe = st.QN;
            if (want_all) {
                Ue = st.UN;
                Te = st.TN;
            }
            return;
        }
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r) {
                if (det_sojourn) {
                    // Deterministic sojourn: the exit metrics are the transient
                    // read at t = d_e.
                    Qe(i, r) = detail::env_det_eval(tranT[e], tranQ[e][i][r], dvals[e]);
                    if (!want_all) continue;
                    Ue(i, r) = detail::env_det_eval(tranT[e], tranU[e][i][r], dvals[e]);
                    Te(i, r) = detail::env_det_eval(tranT[e], tranTp[e][i][r], dvals[e]);
                    continue;
                }
                Qe(i, r) = weighted(tranQ[e][i][r], w, wsum);
                if (!want_all) continue;
                Ue(i, r) = weighted(tranU[e][i][r], w, wsum);
                Te(i, r) = weighted(tranTp[e][i][r], w, wsum);
            }
    }

    /**
     * `post_`: average each stage's transient over WHEN it hands over to each
     * destination, then re-seed every stage from its predecessors.
     */
    void post() {
        const std::size_t E = envObj.nstages();
        std::vector<std::vector<Matrix<double>>> Qexit(
            E, std::vector<Matrix<double>>(E, Matrix<double>(M, K, 0.0)));
        // The other two exit metrics are read only by a state-dependent rate,
        // which is given all three; without one they are never looked at, and
        // the reference computes them in the same loop regardless.
        std::vector<std::vector<Matrix<double>>> Uexit, Texit;
        if (statedep) {
            Uexit.assign(E, std::vector<Matrix<double>>(E, Matrix<double>(M, K, 0.0)));
            Texit.assign(E, std::vector<Matrix<double>>(E, Matrix<double>(M, K, 0.0)));
        }
        for (std::size_t e = 0; e < E; ++e) {
            // Seeded with the stage's own steady state rather than skipped: a
            // metric its transient carries no trajectory for is one it holds
            // CONSTANT, and calling that zero is what reported an open cache
            // model's Q, U and T as identically zero. see `stage_steady`
            Matrix<double> Qd(M, K, 0.0), Ud(M, K, 0.0), Td(M, K, 0.0);
            stage_exit(e, statedep, Qd, Ud, Td);
            for (std::size_t h = 0; h < E; ++h) {
                Qexit[e][h] = Qd;
                if (!statedep) continue;
                Uexit[e][h] = Ud;
                Texit[e][h] = Td;
            }
        }

        // 'meancov' carries a COVARIANCE beside the mean, and it is taken over the
        // same sojourn weights the exit mean uses -- independent of the
        // destination, exactly as Qexit is, since competing exponentials leave the
        // exit time independent of where the switch goes.
        const std::size_t n2 = M * K;
        std::vector<Matrix<double>> Cexit;
        if (meancov) {
            Cexit.assign(E, Matrix<double>(n2, n2, 0.0));
            for (std::size_t e = 0; e < E; ++e)
                if (!tranT[e].empty()) Cexit[e] = stage_exit_cov(e);
        }

        for (std::size_t e = 0; e < E; ++e) {
            if (tranT[e].empty()) continue;
            Matrix<double> Qe(M, K, 0.0);
            std::vector<double> ment(meancov ? n2 : 0, 0.0);
            Matrix<double> sent(meancov ? n2 : 0, meancov ? n2 : 0, 0.0);
            for (std::size_t h = 0; h < E; ++h) {
                const double p = envObj.prob_orig(h, e);
                if (!(p > 0.0)) continue;
                const ResetMarginal& f = envObj.arc(h, e).reset;
                const Matrix<double> reset = f ? f(Qexit[h][e]) : Qexit[h][e];
                if (reset.rows() != M || reset.cols() != K)
                    throw InputError(
                        "SolverENV: a reset policy returned a matrix of the wrong shape");
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t r = 0; r < K; ++r) Qe(i, r) += p * reset(i, r);
                if (!meancov) continue;
                // The reset is an arbitrary map on the means, so the covariance
                // crosses it by the delta method. The MIXTURE over origins is then
                // taken on the SECOND MOMENT, not on the covariances: a convex
                // combination of the C_h alone drops the spread of the per-origin
                // means, which is most of the variance when the stages differ.
                const Matrix<double> R = reset_jacobian(f, Qexit[h][e]);
                const Matrix<double> Ch = congruence(R, Cexit[h]);
                for (std::size_t a = 0; a < n2; ++a) {
                    const double ma = reset(a % M, a / M);
                    ment[a] += p * ma;
                    for (std::size_t b = 0; b < n2; ++b)
                        sent(a, b) += p * (Ch(a, b) + ma * reset(b % M, b / M));
                }
            }
            entry[e] = Qe;
            if (!meancov) continue;
            Matrix<double> Ce(n2, n2, 0.0);
            for (std::size_t a = 0; a < n2; ++a)
                for (std::size_t b = 0; b < n2; ++b)
                    Ce(a, b) = 0.5 * ((sent(a, b) - ment[a] * ment[b]) +
                                      (sent(b, a) - ment[b] * ment[a]));
            centry[e] = Ce;
        }

        // The state-dependent rates come AFTER the entry update, as they do in
        // the reference: the entries just computed used the probOrig of the
        // environment as it was during this iteration, and rewriting the arcs
        // first would blend them with weights from an environment the stage
        // transients were never solved under.
        if (!statedep) return;
        bool touched = false;
        for (std::size_t e = 0; e < E; ++e)
            for (std::size_t h = 0; h < E; ++h) {
                const EnvArc<T>& a = envObj.arc(e, h);
                if (!a.enabled || !a.reset_rates) continue;
                envObj.set_transition_dist(e, h,
                                           a.reset_rates(a.dist, Qexit[e][h], Uexit[e][h],
                                                         Texit[e][h]));
                touched = true;
            }
        if (!touched) return;
        // Everything the analyzer integrates against -- the marked transition
        // processes, the superposed holding times, probEnv and probOrig -- is
        // derived from the arc distributions, so the environment is rebuilt
        // whole rather than patched.
        envObj.init();
        refresh_sojourn_means();
    }

    /**
     * `finish_`: average each stage over its own holding time and blend the
     * stages by their stationary probabilities.
     */
    void finish(EnvSolution& out) {
        const std::size_t E = envObj.nstages();
        out.QExit.assign(E, Matrix<double>(M, K, 0.0));
        out.UExit.assign(E, Matrix<double>(M, K, 0.0));
        out.TExit.assign(E, Matrix<double>(M, K, 0.0));
        for (std::size_t e = 0; e < E; ++e)
            // Seeded with the stage's own steady state rather than skipped, for
            // the reason `post()` is: a metric with no trajectory is a constant
            // and not a zero, and here the zeros would be what the run REPORTS.
            // see `stage_steady`
            stage_exit(e, true, out.QExit[e], out.UExit[e], out.TExit[e]);
        for (std::size_t e = 0; e < E; ++e) {
            const double p = envObj.prob_env[e];
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t r = 0; r < K; ++r) {
                    out.QN(i, r) += p * out.QExit[e](i, r);
                    out.UN(i, r) += p * out.UExit[e](i, r);
                    out.TN(i, r) += p * out.TExit[e](i, r);
                }
        }
        out.Qentry = entry;
        if (!meancov) return;
        // 'meancov' also REPORTS the second moment, mixed over the stages by the
        // same law of total variance the handoff uses. The term
        // sum_e p_e (m_e - m)(m_e - m)' is what the environment itself
        // contributes: two stages with identical within-stage variance but
        // different means still leave the queue length varying, and averaging the
        // per-stage covariances alone would report none of it.
        const std::size_t n = M * K;
        std::vector<double> mval(n, 0.0);
        Matrix<double> sval(n, n, 0.0);
        for (std::size_t e = 0; e < E; ++e) {
            const double p = envObj.prob_env[e];
            if (!(p > 0.0) || tranT[e].empty()) continue;
            const Matrix<double> Ce = stage_exit_cov(e);
            for (std::size_t a = 0; a < n; ++a) {
                const double ma = out.QExit[e](a % M, a / M);
                mval[a] += p * ma;
                for (std::size_t b = 0; b < n; ++b)
                    sval(a, b) += p * (Ce(a, b) + ma * out.QExit[e](b % M, b / M));
            }
        }
        out.QCov = Matrix<double>(n, n, 0.0);
        out.QVar = Matrix<double>(M, K, 0.0);
        for (std::size_t a = 0; a < n; ++a) {
            for (std::size_t b = 0; b < n; ++b)
                out.QCov(a, b) = 0.5 * ((sval(a, b) - mval[a] * mval[b]) +
                                        (sval(b, a) - mval[b] * mval[a]));
            out.QVar(a % M, a / M) = std::max(0.0, out.QCov(a, a));
        }
    }

    /**
     * The grid an exit average is summed on, refined where the SOJOURN WEIGHT
     * CANNOT SEE the integrator's own grid.
     *
     * The exit metric is `sum_k m(t_k) * [F(t_k) - F(t_{k-1})]` over the ODE
     * solver's OUTPUT grid, and that grid is chosen for the horizon rather than
     * for the sojourn: a stage integrated over [0,1e3] and read through an
     * Exp(1) clock puts almost every point where the weight is zero, so the
     * answer becomes an artifact of step placement.
     *
     * The grid is rebuilt UNCONDITIONALLY -- 90% of the points under `5*E[S]`
     * and the rest across the tail -- rather than only when the solver's own
     * grid looks too coarse. A "50 points inside the support is enough" escape
     * (native python's, before this) stops wherever the integrator's steps
     * happened to fall and does not converge: on renv_node_breakdown the sum
     * runs 0.462260, 0.460580, 0.459704, 0.459272, 0.459138, 0.459122 as the
     * point count goes 500 to 5e4, and this engine on a 1e5-point uniform grid
     * answers 0.459171. Rebuilding always is also what makes the four codebases
     * sum the SAME points, which is the property parity needs.
     */
    static constexpr std::size_t kCdfInterp = 5000;

    static std::vector<double> refine_grid(const mam::Map<double>& m,
                                           const std::vector<double>& t) {
        if (t.size() < 2) return t;
        // A DISABLED ARC is the 1 x 1 zero pair, and `map_mean` refuses it by
        // name ("zero arrival rate") rather than returning an infinity. Its
        // weights are identically zero, so the grid it would be summed on
        // cannot matter; leave it alone, exactly as cdf_weights does.
        bool all_zero = true;
        for (std::size_t a = 0; a < m.D1.rows() && all_zero; ++a)
            for (std::size_t b = 0; b < m.D1.cols() && all_zero; ++b)
                if (m.D1(a, b) != 0.0) all_zero = false;
        if (all_zero) return t;
        const double t0 = t.front();
        const double tend = t.back();
        double mean_sojourn = mam::map_mean(m);
        if (!(mean_sojourn > 0.0) || !std::isfinite(mean_sojourn))
            mean_sojourn = (tend - t0) / 10.0;
        double tcdf = std::min(tend, 5.0 * mean_sojourn);
        if (tcdf <= t0) tcdf = tend;
        const std::size_t ndense = static_cast<std::size_t>(0.9 * kCdfInterp);
        const std::size_t ntail = kCdfInterp - ndense;
        const bool with_tail = tcdf < tend && ntail > 1;
        std::vector<double> fine;
        fine.reserve(with_tail ? ndense + ntail : ndense);
        for (std::size_t k = 0; k < ndense; ++k)
            fine.push_back(t0 + (tcdf - t0) * static_cast<double>(k) /
                                    static_cast<double>(ndense - 1));
        if (with_tail)
            for (std::size_t k = 1; k <= ntail; ++k)
                fine.push_back(tcdf + (tend - tcdf) * static_cast<double>(k) /
                                          static_cast<double>(ntail));
        return fine;
    }

    /**
     * The OUTPUT GRID stage `e` is integrated on: the refined one, so the exit
     * average is summed over points the trajectory was actually evaluated at.
     *
     * Interpolating instead cannot recover resolution the trajectory never had
     * -- over [0,1e3] a uniform 1001-point grid carries SIX samples below
     * `5*E[S]` for an `Exp(1)` sojourn, and a piecewise-linear reading of six
     * samples is still six samples' worth of information; it reported a
     * throughput of 0.8099 for a source admitting 0.8. LSODA takes an arbitrary
     * increasing output vector, so the points are simply asked for.
     *
     * The scale is the HOLDING time, which is the sojourn regardless of which
     * destination fires (competing exponentials), so one grid serves post()'s
     * per-destination weights and finish()'s holding-time ones alike.
     */
    std::vector<double> stage_grid(std::size_t e) const {
        std::vector<double> ends(2);
        ends[0] = 0.0;
        ends[1] = opt.timespan_end;
        const std::vector<double> g = refine_grid(envObj.hold_time[e].map(), ends);
        // A disabled holding time leaves `refine_grid` with nothing to say; fall
        // back to the uniform grid the option asks for.
        if (g.size() < 3) return std::vector<double>();
        return g;
    }

    /**
     * The covariance of the station-class queue lengths at the instant stage `e`
     * is left, by the law of total variance over the random sojourn T:
     *
     *   Cov[Q(T)] = E_T[Cov(Q(t)|t)] + Cov_T[E(Q(t)|t)]
     *             = sum_n w_n (C(t_n) + m(t_n) m(t_n)')/sum(w) - m_exit m_exit'
     *
     * with the SAME weights the exit mean uses, so the two are consistent by
     * construction. `C(t)` is the within-stage covariance the stage solver
     * integrated and is absent -- hence zero -- for a stage carrying a first
     * moment only; the timing term survives regardless. A deterministic sojourn
     * contributes no timing variance, so the answer is `C(d)` alone.
     */
    Matrix<double> stage_exit_cov(std::size_t e) const {
        const std::size_t n = M * K;
        Matrix<double> C(n, n, 0.0);
        const std::vector<double>& t = tranT[e];
        if (t.size() < 2) return C;
        const bool have_c = tranC[e].size() == t.size();
        if (det_sojourn) {
            if (!have_c) return C;
            // Clamped to the stage's own grid rather than extrapolated: linear
            // extrapolation of a covariance can leave the positive semidefinite
            // cone, while a convex combination of two members stays inside it.
            double d = std::max(t.front(), std::min(dvals[e], t.back()));
            std::size_t j = t.size() - 1;
            for (std::size_t k = 1; k < t.size(); ++k)
                if (d <= t[k]) {
                    j = k;
                    break;
                }
            const double dt = t[j] - t[j - 1];
            const double a = (dt > 0.0) ? (d - t[j - 1]) / dt : 0.0;
            for (std::size_t p = 0; p < n; ++p)
                for (std::size_t q = 0; q < n; ++q)
                    C(p, q) = (1.0 - a) * tranC[e][j - 1](p, q) + a * tranC[e][j](p, q);
            return C;
        }
        const std::vector<double> w = cdf_weights(envObj.hold_time[e].map(), t);
        double wsum = 0.0;
        for (double v : w) wsum += v;
        if (!(wsum > 0.0)) return C;
        std::vector<double> mexit(n, 0.0);
        for (std::size_t a = 0; a < n; ++a)
            mexit[a] = weighted(tranQ[e][a % M][a / M], w, wsum);
        for (std::size_t a = 0; a < n; ++a) {
            const std::vector<double>& qa = tranQ[e][a % M][a / M];
            for (std::size_t b = 0; b < n; ++b) {
                const std::vector<double>& qb = tranQ[e][b % M][b / M];
                double acc = 0.0;
                for (std::size_t k = 0; k < w.size() && k < qa.size() && k < qb.size(); ++k)
                    acc += w[k] * qa[k] * qb[k];
                if (have_c)
                    for (std::size_t k = 0; k < w.size(); ++k) acc += w[k] * tranC[e][k](a, b);
                C(a, b) = acc / wsum - mexit[a] * mexit[b];
            }
        }
        for (std::size_t a = 0; a < n; ++a)
            for (std::size_t b = a + 1; b < n; ++b) {
                const double v = 0.5 * (C(a, b) + C(b, a));  // drop the rounding asymmetry
                C(a, b) = v;
                C(b, a) = v;
            }
        return C;
    }

    /**
     * The Jacobian of a reset policy at the exit mean, so a covariance can cross
     * the switch as `R C R'`.
     *
     * A reset policy is an arbitrary map on the (station x class) mean queue
     * lengths, and there is no general way to push a second moment through one.
     * The delta method is the first-order image, which is the order the whole
     * mean-field coupling works to. The two NAMED policies are linear and R is
     * then EXACT: the identity for `keep`, zero for `clear`. Linear indexing is
     * column-major, `ir = r*M + i`, the same index space as the covariance.
     */
    Matrix<double> reset_jacobian(const ResetMarginal& f, const Matrix<double>& qexit) const {
        const std::size_t n = M * K;
        Matrix<double> R(n, n, 0.0);
        if (!f) {  // the identity map: no reset on this arc
            for (std::size_t a = 0; a < n; ++a) R(a, a) = 1.0;
            return R;
        }
        const Matrix<double> base = f(qexit);
        if (base.rows() != M || base.cols() != K) return R;
        double scale = 1.0;
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r) scale = std::max(scale, std::fabs(qexit(i, r)));
        const double step = 1e-6 * scale;
        for (std::size_t j = 0; j < n; ++j) {
            Matrix<double> qp = qexit;
            qp(j % M, j / M) += step;
            const Matrix<double> pert = f(qp);
            if (pert.rows() != M || pert.cols() != K) continue;
            for (std::size_t a = 0; a < n; ++a)
                R(a, j) = (pert(a % M, a / M) - base(a % M, a / M)) / step;
        }
        return R;
    }

    /** The congruence `R C R'`, the image of a covariance under a linear map. */
    static Matrix<double> congruence(const Matrix<double>& R, const Matrix<double>& C) {
        const std::size_t n = R.rows();
        Matrix<double> RC(n, n, 0.0), out(n, n, 0.0);
        for (std::size_t a = 0; a < n; ++a)
            for (std::size_t b = 0; b < n; ++b) {
                double acc = 0.0;
                for (std::size_t c = 0; c < n; ++c) acc += R(a, c) * C(c, b);
                RC(a, b) = acc;
            }
        for (std::size_t a = 0; a < n; ++a)
            for (std::size_t b = 0; b < n; ++b) {
                double acc = 0.0;
                for (std::size_t c = 0; c < n; ++c) acc += RC(a, c) * R(b, c);
                out(a, b) = acc;
            }
        return out;
    }

    /** The Stieltjes weights of a transition over the transient grid. */
    std::vector<double> cdf_weights(const mam::Map<double>& m, const std::vector<double>& t) const {
        std::vector<double> w(t.size(), 0.0);
        if (t.size() < 2) return w;
        // A disabled arc is the 1 x 1 zero pair, whose CDF is identically zero.
        bool all_zero = true;
        for (std::size_t a = 0; a < m.D1.rows() && all_zero; ++a)
            for (std::size_t b = 0; b < m.D1.cols() && all_zero; ++b)
                if (m.D1(a, b) != 0.0) all_zero = false;
        if (all_zero) return w;
        const std::vector<double> F = mam::map_cdf(m, t);
        for (std::size_t j = 1; j < t.size(); ++j) w[j] = F[j] - F[j - 1];
        return w;
    }

    static double weighted(const std::vector<double>& v, const std::vector<double>& w,
                           double wsum) {
        double s = 0.0;
        for (std::size_t j = 0; j < v.size() && j < w.size(); ++j) s += v[j] * w[j];
        return s / wsum;
    }

    /**
     * `initFromMarginal` for a fluid stage: the mean queue length of every
     * (station, class) enters in PHASE ONE.
     *
     * The reference does NOT round this for a fluid stage
     * (`roundMarginalForDiscreteSolver` is skipped when the stage solver is
     * SolverFluid), because a fluid state is continuous -- rounding it would
     * quantize the very quantity the fixed point is iterating on.
     */
    std::vector<double> initsol_from_marginal(const qn::NetworkStruct<T>& sn,
                                              const Matrix<double>& Q) const {
        const fluid::FluidLayout L = fluid::fluid_layout(sn);
        std::vector<double> y(L.nstates, 0.0);
        for (std::size_t i = 0; i < sn.nstations && i < Q.rows(); ++i)
            for (std::size_t r = 0; r < sn.nclasses && r < Q.cols(); ++r) {
                if (!L.enabled[i][r]) continue;
                y[L.qidx[i][r]] = std::max(0.0, Q(i, r));
            }
        return y;
    }

    Environment<T>& envObj;
    EnvOptions opt;
    std::size_t M = 0, K = 0;
    bool det_sojourn = false;
    bool statedep = false;  ///< `method = "statedep"`: the arcs are rewritten each sweep
    bool ctmc_stages = false;  ///< `stage_solver = "ctmc"`: enumerate each stage's chain
    bool meancov = false;      ///< `method = "meancov"`: a covariance crosses beside the mean
    bool kp_stages = false;    ///< `stage.method = "kp"`: the Ko-Pender limits run each stage
    bool dae_stages = false;   ///< `stage.method = "dae"`: the linear-noise DAE runs each stage
    std::vector<double> dvals;
    std::vector<Matrix<double>> entry;  ///< per-stage entry queue lengths
    /**
     * `meancov` only: per-stage entry COVARIANCE over the (station, class)
     * pairs, the companion of `entry`, and the WITHIN-STAGE covariance the stage
     * solver integrated, one matrix per point of `tranT`. `tranC[e]` is empty
     * for a stage solver carrying a first moment only, and the coupling then
     * keeps the timing variance of the sojourn alone.
     */
    std::vector<Matrix<double>> centry;
    std::vector<std::vector<Matrix<double>>> tranC;
    std::vector<std::vector<double>> tranT;
    std::vector<std::vector<std::vector<std::vector<double>>>> tranQ, tranU, tranTp;
    /**
     * Per stage, its steady tables and whether they have been asked for: the
     * exit value of every metric its transient carries no trajectory for. Asked
     * at most once per run and only where one is read; see `stage_steady`.
     */
    std::vector<StageSteady> steady;
    std::vector<char> steady_asked;
    /** Per LAYERED stage, the SolverLN that runs it; null for a flat stage. */
    std::vector<std::shared_ptr<ln::SolverLN<T>>> lnsolv;
    /** Per LAYERED stage, where each of its layers sits in the aggregate view. */
    std::vector<ln::LnLayerBlocks> lnblk;
};

}  // namespace env
}  // namespace line

#endif  // LINE_SOLVERS_ENV_SOLVER_ENV_H
