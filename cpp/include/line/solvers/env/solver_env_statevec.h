/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_ENV_SOLVER_ENV_STATEVEC_H
#define LINE_SOLVERS_ENV_SOLVER_ENV_STATEVEC_H

/**
 * SolverENV, `method = "statevec"`: a port of
 * `matlab/src/solvers/ENV/solver_env_statevec_analyzer.m`.
 *
 * WHAT IT DOES DIFFERENTLY FROM THE MEAN-FIELD COUPLING. `solver_env.h` carries
 * only the MARGINAL MEAN queue lengths across an environment switch, so any
 * correlation between stations at the moment of the switch is thrown away. This
 * analyzer carries the whole JOINT distribution instead: each stage is an
 * explicitly enumerated CTMC, and what crosses a switch is the full state
 * probability vector. The two agree exactly when the marginal collapse happens
 * to be lossless and diverge when it is not, which is the entire reason the
 * reference keeps both.
 *
 * THE FIXED POINT, and it is the same shape as the mean-field one. For stage e
 * with generator Q_e and entry distribution pi_enter[e], propagate the transient
 * pi(t) = pi_enter[e] exp(Q_e t) over the stage's time span, then read off
 *   pi_exit[e][h]  the distribution AT the e -> h switch, the expectation of
 *                  pi(t) under the e -> h transition time,
 *   pi_timeavg[e]  the distribution at the END of the sojourn, under the
 *                  superposed holding time,
 * and chain the entries as
 *   pi_enter[e] = sum_h prob_orig(h, e) reset_{h->e}( pi_exit[h][e] ),
 * renormalized, iterated to an L1 fixed point. The blend at the end weights each
 * stage's pi_timeavg by prob_env and maps it to means.
 *
 * THREE SOJOURN REGIMES, and only the third one integrates anything. When every
 * outgoing transition of a stage is exponential the exit distribution equals the
 * time-average and both are the RESOLVENT s pi_0 (sI - Q)^{-1}, a single linear
 * solve with no quadrature -- so that regime stays exact even in rational
 * arithmetic. A deterministic sojourn is uniformization at one point. Only a
 * general phase-type sojourn needs the adaptive transient, and that branch is
 * gated below.
 *
 * WHY THIS FILE DOES NOT DECODE MARGINALS ITSELF. There is a recorded defect
 * here: `_kb/06-solver-catalog.md`, "ENV state-vector Util had drifted from the
 * CTMC analyzer". MATLAB, Java and Python each grew a THIRD copy of the
 * state-to-means reduction beside `solver_ctmc_analyzer`, all three drifted
 * identically in the lld/cd branch, and no parity check saw it because they
 * drifted together -- Util was overstated by 44% on a closed lld model. The
 * reduction here is therefore `ctmc::solver_ctmc_avg_from_pi`, the same function
 * the CTMC analyzer calls, invoked on the same `CtmcResult`. Nothing about the
 * discipline, the load dependence or the loss guards is re-derived.
 *
 * THE ORACLE THAT NEEDS NO EXTERNAL REFERENCE, named in that same entry: an
 * environment whose stages are IDENTICAL cannot change anything the network
 * does, so this solver must reproduce the plain SolverCTMC solution of that one
 * model. It does so EXACTLY here, not approximately, because the stationary
 * distribution is a fixed point of all three sojourn regimes: pi Q = 0 makes the
 * resolvent, the uniformized time-average and the transient all return pi
 * unchanged, and `pre()` seeds from pi.
 *
 * WHAT IS PORTED, and what is refused by name:
 *   ported   the CTMC backend, all three sojourn regimes, per-transition state
 *            reset policies, the cache hit/miss blend of `aggregateCacheBlend_`
 *   ported   both stage backends -- the enumerated CTMC and the MAM/LDQBD one
 *            (`solver_mam_ldqbd` + `solver_mam_ldqbd_flatten`, reduced by
 *            `solver_mam_ldqbd_avg`) -- all three sojourn regimes, per-transition
 *            state reset policies, the cache hit/miss blend of
 *            `aggregateCacheBlend_`
 *   refused  an infinite inner-solver timespan
 *
 * THE TWO BACKENDS DIFFER IN WHAT A STATE IS, and everything downstream of that
 * is shared. A CTMC stage's state is an enumerated row of `sn.space`; a MAM
 * stage's is a (level, phase) pair of the flattened QBD, with no marginal to
 * decode and no cache content in it. So the propagation, the fixed point and the
 * reset policies are backend-agnostic -- they act on a probability vector and a
 * generator -- while the two ENDS are not: `pre()` builds the generator
 * differently and `finish()` reduces it differently
 * (`solver_ctmc_avg_from_pi` against `solver_mam_ldqbd_avg`).
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "line/api/mam/map_cdf.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mc/ctmc_transient.h"
#include "line/api/mc/ctmc_uniformization.h"
#include "line/lang/qn/environment.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/solvers/ctmc/solver_ctmc.h"
#include "line/api/mc/ctmc_solve_reducible.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/mam/solver_mam_ldqbd_flatten.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace env {

/**
 * `resetStateFun{h,e}` of the reference: the state distribution of stage h at
 * the h -> e switch, mapped onto the state space of stage e.
 *
 * This is NOT `ResetMarginal`, which acts on the (nstations x nclasses) mean
 * queue lengths the mean-field coupling carries. The two coexist because the
 * couplings carry different objects, and a policy expressed on means has no
 * canonical lift to a joint distribution.
 */
template <class T>
using ResetStateVec = std::function<std::vector<T>(const std::vector<T>&)>;

/** Options of the state-vector coupling. */
template <class T>
struct EnvStatevecOptions {
    int iter_max = 100;
    double iter_tol = 1e-4;
    /**
     * Which solver runs each stage: `"ctmc"`, the enumerated chain, or `"mam"`,
     * which flattens a level-dependent QBD into a generator over its
     * level/phase states. The MAM backend only accepts what `solver_mam_ldqbd`
     * accepts -- one class and either Delay+Queue or Source+Queue -- and it
     * carries NO cache state, so the hit/miss blend is skipped there, as the
     * reference skips it.
     */
    std::string stage_solver = "ctmc";
    /**
     * `cutoff` handed to `solver_mam_ldqbd` for an OPEN stage, which is where
     * the level count comes from; a closed stage takes its population instead
     * and ignores this.
     */
    std::size_t mam_cutoff = 10;
    /** `options.sojourn`: `stochastic` (default) or `deterministic`. */
    std::string sojourn = "stochastic";
    /** Options handed to each stage's SolverCTMC. */
    ctmc::CtmcOptions stage;
    /**
     * `options.timespan` of the inner solver. The end DEFAULTS TO INFINITY so
     * that a caller who never set one is refused by name rather than served a
     * silently invented horizon -- the reference's own check, and the reason
     * `CTMC(model,'timespan',[0,T])` is spelled out in its error message.
     *
     * SET IT TO A FEW MEAN HOLDING TIMES, NOT TO A LARGE SAFE NUMBER. Only the
     * general phase-type branch reads it, and that branch integrates with
     * `ctmc_transient`, i.e. with ode23, whose maximum step is one tenth of the
     * span. An explicit Runge-Kutta pair is unstable once that step exceeds the
     * chain's fastest time scale, and the local error estimate CANNOT see it:
     * near an equilibrium the estimate vanishes with the derivative, so every
     * oversized step is accepted and roundoff is amplified instead. On the
     * three-job closed model of `test_env_statevec.cpp` (max |Q_ii| = 4) a span
     * of 5 reproduces the stationary law to the last bit while a span of 50
     * drifts to an L1 error of 7e-4. The holding-time CDF has decayed long
     * before the horizon anyway, so a longer one buys nothing.
     */
    double timespan_start = 0.0;
    double timespan_end = std::numeric_limits<double>::infinity();
    /** Per-stage override of `timespan_end`, or empty for the global one. */
    std::vector<double> stage_timespan_end;
    /** `resetStateFun[h][e]`; an empty entry, or an empty table, is identity. */
    std::vector<std::vector<ResetStateVec<T>>> reset_state;
};

/** What the state-vector coupling reports. */
template <class T>
struct EnvStatevecSolution {
    /** Environment-averaged metrics, (nstations x nclasses). */
    Matrix<T> QN, UN, TN;
    /** Per-stage, sojourn-averaged metrics. */
    std::vector<Matrix<T>> QStage, UStage, TStage;
    /** The entry distributions the fixed point converged to, per stage. */
    std::vector<std::vector<T>> pi_enter;
    /** The sojourn-end distributions the blend was taken over, per stage. */
    std::vector<std::vector<T>> pi_timeavg;
    /**
     * `aggregateCacheBlend_`: environment-blended hit and miss probabilities,
     * keyed by the 1-based node index of each Cache. The reference writes these
     * onto the stage-one node objects with `setResultHitProb`; a struct-level
     * port has no node object to write to, so they are reported here instead.
     * An entry is NaN for a class the cache never serves.
     */
    std::map<std::size_t, std::vector<T>> hit_prob, miss_prob;
    int iterations = 0;
    bool converged = false;
};

/**
 * The state-vector environment solver.
 *
 * `envObj` must already carry a model per stage; `init()` is called here, as
 * `SolverENV.init` does.
 */
template <class T>
class SolverEnvStatevec {
public:
    SolverEnvStatevec(Environment<T>& e, const EnvStatevecOptions<T>& o) : envObj(e), opt(o) {
        init();
    }

    EnvStatevecSolution<T> solve() {
        const std::size_t E = envObj.nstages();
        EnvStatevecSolution<T> out;

        pre();
        int it = 0;
        for (it = 1; it <= opt.iter_max; ++it) {
            for (std::size_t e = 0; e < E; ++e) analyze(e);
            post();
            if (converged()) {
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
        mam_backend = (opt.stage_solver == "mam");
        if (!mam_backend && opt.stage_solver != "ctmc")
            throw UnsupportedError(
                "SolverENV statevec: stage solver '" + opt.stage_solver +
                "' is not available; the state-vector coupling needs an EXPLICIT generator over a "
                "state space it can propagate a distribution across, which SolverCTMC exposes by "
                "enumeration and SolverMAM by flattening its level-dependent QBD blocks");
        if (opt.sojourn != "stochastic" && opt.sojourn != "deterministic")
            throw InputError("SolverENV statevec: unknown sojourn '" + opt.sojourn + "'");

        // A LAYERED STAGE HAS NO SINGLE GENERATOR to propagate a distribution
        // across -- an LQN decomposes into one network per layer, and the joint
        // law over the union of them is not what any layer solver produces. The
        // reference refuses it in the same place and for the same reason
        // (`@@SolverENV/SolverENV.m`, where `method` selects the state-vector
        // analyzer), so this is its error and not a limit of this port.
        envObj.reject_lqn_stages(
            "SolverENV statevec",
            "the state-vector coupling propagates a joint distribution over ONE stage "
            "generator, which a layered model does not have -- it decomposes into a network "
            "per layer");
        envObj.init();
        const std::size_t E = envObj.nstages();
        M = envObj.stage(0).model.nstations;
        K = envObj.stage(0).model.nclasses;
        for (std::size_t e = 1; e < E; ++e)
            if (envObj.stage(e).model.nstations != M || envObj.stage(e).model.nclasses != K)
                throw InputError(
                    "SolverENV statevec: every stage must have the same stations and classes; the "
                    "metrics are blended entrywise across them");

        // The reference checks the horizon PER STAGE and names the stage that
        // is missing one, because each stage carries its own inner solver.
        tspan_end.assign(E, opt.timespan_end);
        for (std::size_t e = 0; e < E; ++e) {
            if (e < opt.stage_timespan_end.size()) tspan_end[e] = opt.stage_timespan_end[e];
            if (!std::isfinite(tspan_end[e]) || !(tspan_end[e] > opt.timespan_start))
                throw InputError(
                    "SolverENV statevec: the statevec analyzer requires a finite inner-solver "
                    "timespan for stage " +
                    std::to_string(e + 1) + ", e.g. CTMC(model,'timespan',[0,T])");
        }

        det_sojourn = (opt.sojourn == "deterministic");
        stages.assign(E, ctmc::CtmcSolution<T>());
        mam_ld.assign(E, mam::LdqbdBlocks<T>());
        mam_flat.assign(E, mam::LdqbdFlat<T>());
        pi_enter.assign(E, std::vector<T>());
        pi_enter_prev.assign(E, std::vector<T>());
        pi_exit.assign(E, std::vector<std::vector<T>>());
        pi_timeavg.assign(E, std::vector<T>());
        dvals.assign(E, 0.0);
        if (det_sojourn)
            for (std::size_t e = 0; e < E; ++e)
                dvals[e] = std::max(mam::map_mean(envObj.hold_time[e].map()),
                                    std::numeric_limits<double>::epsilon());
    }

    /**
     * `pre_`: build each stage's chain once and warm-start its entry
     * distribution from that stage's own stationary law, which is a valid
     * probability vector over its state space whatever the environment does.
     *
     * THE CHAIN COMES FROM THE ANALYZER, not from a bare `solver_ctmc` call as
     * the reference's `pre_` makes. The analyzer restricts a reducible generator
     * to the weakly connected component of the model's initial state; solving
     * the unrestricted generator instead spreads mass over states the model can
     * never occupy. That restriction is also what makes the identical-stage
     * oracle hold EXACTLY rather than up to the mass on unreachable states.
     */
    void pre() {
        const std::size_t E = envObj.nstages();
        for (std::size_t e = 0; e < E; ++e) {
            if (mam_backend) {
                // The MAM backend has no enumerated space: the stage is a
                // level-dependent QBD, and flattening its blocks is what turns
                // it into something a distribution can be propagated across.
                // `Nlev` is finite by construction here -- the closed population,
                // or the open truncation `cutoff` -- so the flat matrix exists.
                mam::MamOptions mo;
                mo.method = "default";
                mo.cutoff = opt.mam_cutoff;
                mam_ld[e] = mam::solver_mam_ldqbd(envObj.stage(e).model, mo).ld;
                mam_flat[e] = mam::solver_mam_ldqbd_flatten(mam_ld[e]);
            } else {
                stages[e] = ctmc::solver_ctmc_analyzer(envObj.stage(e).model, opt.stage);
            }
        }
        seed_entry_distributions();
        pi_enter_prev = pi_enter;
    }

    /**
     * Warm start. A stage's OWN stationary law is not a usable seed: a stage
     * that is individually unstable or critical (arrival rate >= its own
     * service rate) has no stationary law at all, and the reducible solver then
     * returns the stationary law of the TRUNCATED generator, which piles mass
     * against the truncation wall and whose mean grows linearly with the cutoff
     * -- for a critical M/M/1 truncated at N it is uniform, with mean N/2.
     *
     * The fixed point chained in `post` is EXACT: it is the stationary equation
     * of the joint (queue,stage) chain, phi_e = (sum_h phi_h q_he)(s_e I-Q_e)^-1.
     * It does contract to the right answer from that seed, but the number of
     * sweeps it needs grows with the cutoff, so at a finite `iter_max` the
     * reported result drifts FURTHER from the truth as the cutoff is RAISED --
     * the natural response to a suspect number makes it worse.
     *
     * Seed instead from the environment-averaged generator sum_e prob_env[e]*Q_e,
     * positive recurrent exactly when the model is stable on average, which is
     * the regime in which the answer exists at all; its stationary law is
     * therefore cutoff-independent. Averaging needs one common state space, so
     * when the stages differ in size fall back to the per-stage law, which is
     * the best available and no worse than before.
     */
    void seed_entry_distributions() {
        const std::size_t E = envObj.nstages();

        // The per-stage law, which is what the seed falls back to and also what
        // selects the recurrent class below. THE CHAIN COMES FROM THE ANALYZER,
        // not from a bare solve: the analyzer restricts a reducible generator to
        // the weakly connected component of the model's initial state AND uses
        // that state to pick among the BSCCs that survive, so `stages[e].pi` is
        // not reproducible from `gen(e)` alone. Re-deriving it here is what
        // would spread mass over states the model can never occupy.
        std::vector<std::vector<T>> own(E);
        for (std::size_t e = 0; e < E; ++e) {
            // The MAM backend has no analyzer law, and the REDUCIBLE solver:
            // the flat generator of a truncated open QBD need not be irreducible.
            own[e] = mam_backend ? normalized(mc::ctmc_solve_reducible(mam_flat[e].Q).pi)
                                 : normalized(stages[e].pi);
        }

        bool sameSpace = E > 1;
        for (std::size_t e = 1; e < E && sameSpace; ++e)
            sameSpace = (state_count(e) == state_count(0));

        std::vector<T> shared;
        if (sameSpace) {
            std::vector<double> w(E, 1.0 / static_cast<double>(E));
            double wsum = 0.0;
            bool usable = envObj.prob_env.size() == E;
            for (std::size_t e = 0; e < E && usable; ++e) {
                if (!std::isfinite(envObj.prob_env[e]) || envObj.prob_env[e] < 0.0) usable = false;
                else wsum += envObj.prob_env[e];
            }
            // Stage probabilities unavailable or degenerate: weight stages equally.
            if (usable && wsum > 0.0)
                for (std::size_t e = 0; e < E; ++e) w[e] = envObj.prob_env[e] / wsum;

            const std::size_t n = state_count(0);
            Matrix<T> Qbar(n, n, num_traits<T>::from_double(0.0));
            for (std::size_t e = 0; e < E; ++e) {
                const Matrix<T>& Qe = gen(e);
                const T we = num_traits<T>::from_double(w[e]);
                for (std::size_t i = 0; i < n; ++i)
                    for (std::size_t j = 0; j < n; ++j) Qbar(i, j) += we * Qe(i, j);
            }
            // Seeded with the blended per-stage law so the recurrent class the
            // ANALYZER chose is the one carried forward. This is also what keeps
            // the identical-stage oracle EXACT: with equal stages Qbar == Q_0 and
            // pi0 == pi, pi is already stationary for Q_0, so the solve returns
            // it unchanged and the fixed point is still reached at sweep one.
            std::vector<T> pi0(n, num_traits<T>::from_double(0.0));
            for (std::size_t e = 0; e < E; ++e) {
                const T we = num_traits<T>::from_double(w[e]);
                for (std::size_t i = 0; i < n; ++i) pi0[i] += we * own[e][i];
            }
            shared = normalized(mc::ctmc_solve_reducible(Qbar, pi0).pi);
        }

        for (std::size_t e = 0; e < E; ++e) {
            pi_enter[e] = shared.empty() ? own[e] : shared;
        }
    }

    /** The stage's generator, whichever backend built it. */
    const Matrix<T>& gen(std::size_t e) const {
        return mam_backend ? mam_flat[e].Q : stages[e].chain.Q;
    }

    /** The size of the stage's state space, whichever backend built it. */
    std::size_t state_count(std::size_t e) const {
        return mam_backend ? mam_flat[e].levelOf.size() : stages[e].chain.space.size();
    }

    /**
     * `analyze_`: propagate the entry distribution of stage e through its
     * sojourn, recording where it lands at each destination and at the end.
     */
    void analyze(std::size_t e) {
        const std::size_t E = envObj.nstages();
        const Matrix<T>& Q = gen(e);
        const std::vector<T> pi0 = pi_enter[e];

        if (det_sojourn) {
            // A deterministic sojourn has no CDF to integrate against: the exit
            // is pi0 exp(Q d) and the blend is the exact time-average over
            // [0, d], both from uniformization, and both the same toward every
            // destination since the switch instant does not depend on where it
            // goes.
            std::vector<T> avg, ex;
            time_average(pi0, Q, dvals[e], avg, ex);
            pi_exit[e].assign(E, std::vector<T>());
            for (std::size_t h = 0; h < E; ++h)
                if (envObj.arc(e, h).enabled) pi_exit[e][h] = ex;
            pi_timeavg[e] = avg;
            return;
        }

        double s_e = 0.0;
        if (exp_sojourn(e, s_e)) {
            // Every outgoing transition exponential: the memorylessness makes
            // the distribution at the switch equal to the time-average, and both
            // are the resolvent s pi0 (sI - Q)^{-1}. One linear solve, no
            // quadrature, and no dependence on the horizon at all.
            const std::vector<T> res = resolvent(pi0, Q, s_e);
            pi_exit[e].assign(E, std::vector<T>());
            for (std::size_t h = 0; h < E; ++h)
                if (envObj.arc(e, h).enabled) pi_exit[e][h] = res;
            pi_timeavg[e] = res;
            return;
        }

        // General phase-type sojourn: the transient on the integrator's own
        // adaptive grid, averaged against the transition and holding-time CDFs.
        if constexpr (!num_traits<T>::has_transcendental) {
            throw UnsupportedError(
                "SolverENV statevec: a non-exponential environment sojourn needs the transient "
                "pi0 exp(Qt) from ctmc_transient, which is an adaptive approximation governed by "
                "a tolerance and has no exact value in rational arithmetic; rerun with "
                "--arith double, or use exponential transitions, whose resolvent is exact");
        } else {
            const mc::TransientResult<T> tr =
                mc::ctmc_transient(Q, pi0, num_traits<T>::from_double(opt.timespan_start),
                                   num_traits<T>::from_double(tspan_end[e]));
            std::vector<double> t(tr.t.size());
            for (std::size_t j = 0; j < tr.t.size(); ++j) t[j] = num_traits<T>::to_double(tr.t[j]);

            pi_exit[e].assign(E, std::vector<T>());
            for (std::size_t h = 0; h < E; ++h)
                pi_exit[e][h] = stieltjes(tr.pi, cdf_weights(envObj.proc[e][h].map(), t));

            const std::vector<T> avg =
                stieltjes(tr.pi, cdf_weights(envObj.hold_time[e].map(), t));
            if (!avg.empty()) {
                pi_timeavg[e] = avg;
            } else {
                // Degenerate holding time on this grid: fall back to the
                // terminal distribution, as the reference does.
                const std::size_t last = tr.pi.rows() - 1, n = tr.pi.cols();
                std::vector<T> tail(n);
                for (std::size_t j = 0; j < n; ++j) tail[j] = tr.pi(last, j);
                pi_timeavg[e] = tail;
            }
        }
    }

    /**
     * `post_`: chain the entry distributions, carrying each stage's exit
     * distributions into the stages they feed, weighted by prob_orig.
     */
    void post() {
        const std::size_t E = envObj.nstages();
        pi_enter_prev = pi_enter;
        std::vector<std::vector<T>> next(E);

        for (std::size_t e = 0; e < E; ++e) {
            const std::size_t n = state_count(e);
            std::vector<T> acc(n, num_traits<T>::from_int(0));
            double wsum = 0.0;
            for (std::size_t h = 0; h < E; ++h) {
                const double po = envObj.prob_orig(h, e);
                if (!(po > 0.0)) continue;
                if (pi_exit[h].size() <= e || pi_exit[h][e].empty()) continue;
                const std::vector<T> pex = reset_apply(h, e, pi_exit[h][e]);
                if (pex.size() != n)
                    throw InputError(
                        "SolverENV statevec: reset_state[" + std::to_string(h) + "][" +
                        std::to_string(e) + "] returned a " + std::to_string(pex.size()) +
                        "-element vector but stage " + std::to_string(e + 1) + " has " +
                        std::to_string(n) +
                        " states; supply a reset that maps the state space of stage " +
                        std::to_string(h + 1) + " onto that of stage " + std::to_string(e + 1));
                const T w = num_traits<T>::from_double(po);
                for (std::size_t s = 0; s < n; ++s) acc[s] += T(w * pex[s]);
                wsum += po;
            }
            if (wsum > 0.0) {
                const T w = num_traits<T>::from_double(wsum);
                for (std::size_t s = 0; s < n; ++s) acc[s] = T(acc[s] / w);
            } else {
                acc = pi_enter[e];  // no inflow this cycle: retain the estimate
            }
            next[e] = normalized(acc);
        }
        pi_enter = next;
    }

    /**
     * `converged_`: the max L1 change of any entry distribution over a full
     * cycle. It is the DISTRIBUTIONS that are compared and not the means,
     * because two different joint laws can share every marginal mean.
     */
    bool converged() const {
        const std::size_t E = envObj.nstages();
        double l1 = 0.0;
        for (std::size_t e = 0; e < E; ++e) {
            const std::vector<T>& a = pi_enter[e];
            const std::vector<T>& b = pi_enter_prev[e];
            if (a.empty() || b.empty() || a.size() != b.size()) return false;
            double d = 0.0;
            for (std::size_t s = 0; s < a.size(); ++s)
                d += std::fabs(num_traits<T>::to_double(a[s]) - num_traits<T>::to_double(b[s]));
            l1 = std::max(l1, d);
        }
        if (!std::isfinite(l1)) return false;
        return l1 < opt.iter_tol;
    }

    /**
     * `finish_`: map each stage's sojourn-end distribution to means with the
     * CTMC analyzer's own reduction, and blend by prob_env.
     */
    void finish(EnvStatevecSolution<T>& out) {
        const std::size_t E = envObj.nstages();
        const T zero = num_traits<T>::from_int(0);
        out.QN = Matrix<T>(M, K, zero);
        out.UN = Matrix<T>(M, K, zero);
        out.TN = Matrix<T>(M, K, zero);
        out.QStage.assign(E, Matrix<T>(M, K, zero));
        out.UStage.assign(E, Matrix<T>(M, K, zero));
        out.TStage.assign(E, Matrix<T>(M, K, zero));

        for (std::size_t e = 0; e < E; ++e) {
            if (pi_timeavg[e].empty()) continue;
            Matrix<T> QNe, UNe, TNe;
            if (mam_backend) {
                // The LD-QBD reduction reports one column, the model being
                // single-class by construction; K is 1 here for the same reason.
                const mam::LdqbdAvg<T> a =
                    mam::solver_mam_ldqbd_avg(mam_ld[e], pi_timeavg[e], mam_flat[e].levelOf);
                QNe = a.QN;
                UNe = a.UN;
                TNe = a.TN;
            } else {
                const ctmc::CtmcAvg<T> a = ctmc::solver_ctmc_avg_from_pi(
                    envObj.stage(e).model, stages[e].chain, pi_timeavg[e]);
                QNe = a.QN;
                UNe = a.UN;
                TNe = a.TN;
            }
            out.QStage[e] = QNe;
            out.UStage[e] = UNe;
            out.TStage[e] = TNe;
            const T p = num_traits<T>::from_double(envObj.prob_env[e]);
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t r = 0; r < K && r < QNe.cols(); ++r) {
                    out.QN(i, r) += T(p * QNe(i, r));
                    out.UN(i, r) += T(p * UNe(i, r));
                    out.TN(i, r) += T(p * TNe(i, r));
                }
        }
        out.pi_enter = pi_enter;
        out.pi_timeavg = pi_timeavg;
        cache_blend(out);
    }

    /**
     * `aggregateCacheBlend_`: hit and miss probabilities as the ratio of the
     * environment-blended hit and miss THROUGHPUTS.
     *
     * The ratio has to be taken after the blend and not before it: a per-stage
     * ratio weighted by prob_env would be an average of ratios, which is not
     * the ratio the environment actually exhibits unless every stage carries
     * the same total rate.
     */
    void cache_blend(EnvStatevecSolution<T>& out) const {
        // A MAM stage's state is a (level, phase) pair of the flattened QBD: it
        // carries no cache content and no departure-rate table to read hit and
        // miss throughputs off, so there is nothing to blend. The reference
        // returns here for the same reason.
        if (mam_backend) return;
        const std::size_t E = envObj.nstages();
        const qn::NetworkStruct<T>& sn1 = envObj.stage(0).model;
        const T zero = num_traits<T>::from_int(0);

        for (const auto& kv : sn1.nodeparam) {
            const std::size_t ind = kv.first;
            const std::size_t isf = sn1.stateful_index(ind);
            if (isf == 0) continue;
            std::vector<T> hitT(K, zero), missT(K, zero);
            for (std::size_t e = 0; e < E; ++e) {
                if (pi_timeavg[e].empty()) continue;
                const std::vector<T> pv = normalized(pi_timeavg[e]);
                const qn::NetworkStruct<T>& sne = envObj.stage(e).model;
                const auto it = sne.nodeparam.find(ind);
                if (it == sne.nodeparam.end()) continue;
                const qn::CacheParam<T>& np = it->second;
                const T w = num_traits<T>::from_double(envObj.prob_env[e]);
                const auto& dr = stages[e].chain.dep_rates;
                for (std::size_t k = 0; k < K && k < np.hitclass.size(); ++k) {
                    const std::size_t hc = np.hitclass[k];
                    const std::size_t mc = k < np.missclass.size() ? np.missclass[k] : 0;
                    if (hc == 0 || mc == 0) continue;
                    for (std::size_t s = 0; s < pv.size(); ++s) {
                        hitT[k] += T(w * pv[s] * dr[s][isf - 1][hc - 1]);
                        missT[k] += T(w * pv[s] * dr[s][isf - 1][mc - 1]);
                    }
                }
            }
            const T nan = num_traits<T>::from_double(std::numeric_limits<double>::quiet_NaN());
            std::vector<T> hp(K, nan), mp(K, nan);
            for (std::size_t k = 0; k < K; ++k) {
                const T tot = T(hitT[k] + missT[k]);
                if (num_traits<T>::to_double(tot) > 0) {
                    hp[k] = T(hitT[k] / tot);
                    mp[k] = T(missT[k] / tot);
                }
            }
            out.hit_prob[ind] = hp;
            out.miss_prob[ind] = mp;
        }
    }

    // ---- the three sojourn regimes ---------------------------------------

    /** True when every enabled transition out of e is exponential; `s` is their total rate. */
    bool exp_sojourn(std::size_t e, double& s) const {
        const std::size_t E = envObj.nstages();
        s = 0.0;
        bool any = false;
        for (std::size_t h = 0; h < E; ++h) {
            const EnvArc<T>& a = envObj.arc(e, h);
            if (!a.enabled) continue;
            if (a.dist.type != lang::ProcessType::EXP) return false;
            s += num_traits<T>::to_double(a.dist.D1(0, 0));
            any = true;
        }
        return any && s > 0.0;
    }

    /** `s pi0 (sI - Q)^{-1}`, solved as `(sI - Q)^T x = s pi0` on the columns. */
    static std::vector<T> resolvent(const std::vector<T>& pi0, const Matrix<T>& Q, double s) {
        const std::size_t n = Q.rows();
        const T sT = num_traits<T>::from_double(s);
        Matrix<T> A(n, n, num_traits<T>::from_int(0));
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                A(j, i) = T((i == j ? sT : num_traits<T>::from_int(0)) - Q(i, j));
        std::vector<T> b(n);
        for (std::size_t i = 0; i < n; ++i) b[i] = T(sT * pi0[i]);
        return line::solve(A, b);
    }

    /** `ctmc_timeaverage` at the mean holding time, gated on transcendentals. */
    static void time_average(const std::vector<T>& pi0, const Matrix<T>& Q, double d,
                             std::vector<T>& avg, std::vector<T>& ex) {
        if constexpr (!num_traits<T>::has_transcendental) {
            throw UnsupportedError(
                "SolverENV statevec: a deterministic sojourn needs ctmc_timeaverage, whose "
                "uniformization carries Poisson weights exp(-qt) that are not rational; rerun "
                "with --arith double");
        } else {
            const mc::TimeAverageResult<T> r =
                mc::ctmc_timeaverage(pi0, Q, num_traits<T>::from_double(d));
            avg = r.piTimeAvg;
            ex = r.piExit;
        }
    }

    /**
     * The Stieltjes weights of a transition over the transient grid, w_j =
     * F(t_j) - F(t_{j-1}) with w_1 = 0.
     *
     * NO DENSITY IS INVOLVED, deliberately: the increments of the CDF are the
     * measure the transition induces on this grid, and a deterministic
     * transition has increments but no density.
     */
    static std::vector<double> cdf_weights(const mam::Map<double>& m,
                                           const std::vector<double>& t) {
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

    /** `(w' * pit) / sum(w)`; empty when the weights carry no mass or a NaN. */
    static std::vector<T> stieltjes(const Matrix<T>& pit, const std::vector<double>& w) {
        double sw = 0.0;
        for (double v : w) {
            if (std::isnan(v)) return std::vector<T>();
            sw += v;
        }
        if (!(sw > 0.0)) return std::vector<T>();
        const std::size_t n = pit.cols();
        std::vector<T> out(n, num_traits<T>::from_int(0));
        for (std::size_t j = 0; j < w.size() && j < pit.rows(); ++j) {
            if (w[j] == 0.0) continue;
            const T wj = num_traits<T>::from_double(w[j] / sw);
            for (std::size_t s = 0; s < n; ++s) out[s] += T(wj * pit(j, s));
        }
        return out;
    }

    // ---- small helpers ----------------------------------------------------

    std::vector<T> reset_apply(std::size_t h, std::size_t e, const std::vector<T>& p) const {
        if (h < opt.reset_state.size() && e < opt.reset_state[h].size() && opt.reset_state[h][e])
            return opt.reset_state[h][e](p);
        return p;
    }

    /** Clamp the numerical dust below zero and renormalize to a distribution. */
    static std::vector<T> normalized(const std::vector<T>& p) {
        const T zero = num_traits<T>::from_int(0);
        std::vector<T> q = p;
        T tot = zero;
        for (std::size_t s = 0; s < q.size(); ++s) {
            if (num_traits<T>::to_double(q[s]) < 0) q[s] = zero;
            tot += q[s];
        }
        if (num_traits<T>::to_double(tot) > 0)
            for (std::size_t s = 0; s < q.size(); ++s) q[s] = T(q[s] / tot);
        return q;
    }

    Environment<T>& envObj;
    EnvStatevecOptions<T> opt;
    std::size_t M = 0, K = 0;
    bool det_sojourn = false;
    std::vector<double> dvals, tspan_end;
    std::vector<ctmc::CtmcSolution<T>> stages;
    bool mam_backend = false;
    std::vector<mam::LdqbdBlocks<T>> mam_ld;
    std::vector<mam::LdqbdFlat<T>> mam_flat;
    std::vector<std::vector<T>> pi_enter, pi_enter_prev, pi_timeavg;
    std::vector<std::vector<std::vector<T>>> pi_exit;  ///< pi_exit[e][h]
};

/** Solve in one call, for a caller with no use for the solver object. */
template <class T>
EnvStatevecSolution<T> solver_env_statevec(Environment<T>& e, const EnvStatevecOptions<T>& o) {
    return SolverEnvStatevec<T>(e, o).solve();
}

}  // namespace env
}  // namespace line

#endif  // LINE_SOLVERS_ENV_SOLVER_ENV_STATEVEC_H
