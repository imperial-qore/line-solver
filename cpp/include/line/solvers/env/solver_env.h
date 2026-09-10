/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_ENV_SOLVER_ENV_H
#define LINE_SOLVERS_ENV_SOLVER_ENV_H

/**
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
 *   ported   `meanfield` (the default), `smp` (which selects no analyzer of its
 *            own -- see `init`), `statedep` (state-dependent environment rates,
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
#include "line/lang/qn/environment.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/ctmc/solver_ctmc_transient.h"
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
        if (opt.method == "statevec" || opt.method == "blend")
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
        // `smp` is the reference's escape hatch for a NON-MARKOVIAN environment
        // transition, and it selects no analyzer of its own: it only lifts the
        // constructor check that every arc is Markovian, after which the
        // mean-field analyzer runs unchanged. An arc here is a (D0, D1) pair by
        // construction, so there is no check to lift and nothing else to do.
        statedep = (opt.method == "statedep");
        // "default" IS the reference's spelling of this coupling: SolverENV.m
        // installs solver_env_meanfield_analyzer for every method that is not
        // statevec/blend, and its listValidMethods leads with 'default'. Only
        // the CLI normalised it away, so a caller constructing SolverEnv
        // directly with the reference's own name was refused.
        if (opt.method != "meanfield" && opt.method != "default" && opt.method != "smp"
            && !statedep)
            throw UnsupportedError("SolverENV: unknown method '" + opt.method + "'");
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
        tranT.assign(E, std::vector<double>());
        tranQ.assign(E, std::vector<std::vector<std::vector<double>>>());
        tranU.assign(E, std::vector<std::vector<std::vector<double>>>());
        tranTp.assign(E, std::vector<std::vector<std::vector<double>>>());
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
            fluid::FluidOptions fo = opt.stage;
            fo.init_sol.clear();
            const std::vector<fluid::FluidTranPoint> tr = fluid::solver_fluid_transient(
                envObj.stage(e).model, fo, opt.timespan_end, opt.tran_points);
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
        fluid::FluidOptions fo = opt.stage;
        fo.init_sol = initsol_from_marginal(sn, entry[e]);
        const std::vector<fluid::FluidTranPoint> tr = fluid::solver_fluid_transient(
            sn, fo, opt.timespan_end, opt.tran_points, stage_grid(e));
        for (const fluid::FluidTranPoint& p : tr) {
            tranT[e].push_back(p.t);
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t r = 0; r < K; ++r) {
                    tranQ[e][i][r].push_back(p.QN(i, r));
                    tranU[e][i][r].push_back(p.UN(i, r));
                    tranTp[e][i][r].push_back(p.TN(i, r));
                }
        }
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
            if (tranT[e].empty()) continue;
            if (det_sojourn) {
                Matrix<double> Qd(M, K, 0.0), Ud(M, K, 0.0), Td(M, K, 0.0);
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t r = 0; r < K; ++r) {
                        Qd(i, r) = detail::env_det_eval(tranT[e], tranQ[e][i][r], dvals[e]);
                        if (!statedep) continue;
                        Ud(i, r) = detail::env_det_eval(tranT[e], tranU[e][i][r], dvals[e]);
                        Td(i, r) = detail::env_det_eval(tranT[e], tranTp[e][i][r], dvals[e]);
                    }
                for (std::size_t h = 0; h < E; ++h) {
                    Qexit[e][h] = Qd;
                    if (!statedep) continue;
                    Uexit[e][h] = Ud;
                    Texit[e][h] = Td;
                }
                continue;
            }
            // THE WEIGHT IS THE STAGE SOJOURN, NOT THE e -> h CLOCK, and the
            // reference says so where it builds it: competing exponentials leave
            // the exit TIME independent of which destination won, so the exit
            // average does not depend on h at all -- h enters only through the
            // reset applied on the way in. Weighting by `proc[e][h]` instead read
            // the transient over the mean of ONE risk rather than of their
            // minimum: on renv_twostages_repairmen, whose Stage2 competes a 0.5
            // self arc with a 0.5 arc back, that is a mean of 2 against the
            // sojourn's 1, and it reported Queue1 QLen 0.55882 against MATLAB's
            // 0.55550 with the two stations' throughputs 1.4% apart in a closed
            // cycle that admits one throughput.
            const std::vector<double> w = cdf_weights(envObj.hold_time[e].map(), tranT[e]);
            double wsum = 0.0;
            for (double v : w) wsum += v;
            if (!(wsum > 0.0)) continue;  // the stage never leaves
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t r = 0; r < K; ++r) {
                    const double q = weighted(tranQ[e][i][r], w, wsum);
                    double u = 0.0, t = 0.0;
                    if (statedep) {
                        u = weighted(tranU[e][i][r], w, wsum);
                        t = weighted(tranTp[e][i][r], w, wsum);
                    }
                    for (std::size_t h = 0; h < E; ++h) {
                        Qexit[e][h](i, r) = q;
                        if (!statedep) continue;
                        Uexit[e][h](i, r) = u;
                        Texit[e][h](i, r) = t;
                    }
                }
        }

        for (std::size_t e = 0; e < E; ++e) {
            if (tranT[e].empty()) continue;
            Matrix<double> Qe(M, K, 0.0);
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
            }
            entry[e] = Qe;
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
        for (std::size_t e = 0; e < E; ++e) {
            if (tranT[e].empty()) continue;
            if (det_sojourn) {
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t r = 0; r < K; ++r) {
                        out.QExit[e](i, r) = detail::env_det_eval(tranT[e], tranQ[e][i][r], dvals[e]);
                        out.UExit[e](i, r) = detail::env_det_eval(tranT[e], tranU[e][i][r], dvals[e]);
                        out.TExit[e](i, r) = detail::env_det_eval(tranT[e], tranTp[e][i][r], dvals[e]);
                    }
                continue;
            }
            const std::vector<double> w = cdf_weights(envObj.hold_time[e].map(), tranT[e]);
            double wsum = 0.0;
            for (double v : w) wsum += v;
            if (!(wsum > 0.0)) continue;
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t r = 0; r < K; ++r) {
                    out.QExit[e](i, r) = weighted(tranQ[e][i][r], w, wsum);
                    out.UExit[e](i, r) = weighted(tranU[e][i][r], w, wsum);
                    out.TExit[e](i, r) = weighted(tranTp[e][i][r], w, wsum);
                }
        }
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
    std::vector<double> dvals;
    std::vector<Matrix<double>> entry;  ///< per-stage entry queue lengths
    std::vector<std::vector<double>> tranT;
    std::vector<std::vector<std::vector<std::vector<double>>>> tranQ, tranU, tranTp;
    /** Per LAYERED stage, the SolverLN that runs it; null for a flat stage. */
    std::vector<std::shared_ptr<ln::SolverLN<T>>> lnsolv;
    /** Per LAYERED stage, where each of its layers sits in the aggregate view. */
    std::vector<ln::LnLayerBlocks> lnblk;
};

}  // namespace env
}  // namespace line

#endif  // LINE_SOLVERS_ENV_SOLVER_ENV_H
