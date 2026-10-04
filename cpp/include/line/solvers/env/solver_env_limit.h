/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_ENV_SOLVER_ENV_LIMIT_H
#define LINE_SOLVERS_ENV_SOLVER_ENV_LIMIT_H

/**
 * @file
 * @ingroup line_solvers
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
 * A STAGE HOLDING A CACHE is solved, not refused. The hit, miss and delayed
 * fractions each stage reports are blended the way the reference's
 * `accumCacheMetric` does after 2026-09-13: as RATES weighted by that stage's
 * per-class arrival INTO THE CACHE, normalised once at the end. Rates and not
 * ratios, because a stage's hit ratio is conditional on arriving in that stage,
 * so a probEnv-weighted mean of ratios disagrees with the Sink-row hit
 * throughput -- `sum_e p_e*lambda_e*h_e` -- in the SAME node table. It is also
 * the convention `aggregate_cache_meanfield` and the statevec `cache_blend`
 * already hold, each with a comment saying so.
 *
 * WHICH SOLVER RUNS A STAGE. The reference hands `SolverENV` a FACTORY and calls
 * it for every stage, so "the stage solver is the calling solver" is literally
 * true there. The second constructor here takes an `EnvStageAvgFn`, the same
 * thing as a callable, and the no-factory path runs `stage_solver`: the fluid
 * analyzer, or with `ctmc` the enumerated chain at `stage_cutoff` (an ensemble built
 * on SolverCTMC stages). Without a factory an MVA or NC caller would silently be
 * answered by a fluid solve, and even an FLD caller would lose its own
 * `FluidOptions` (method, timespan, tolerances) on the way in.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include "line/api/sn/sn_node_metrics.h"
#include "line/lang/qn/environment.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/solvers/cache_metrics.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/env/solver_env.h"
#include "line/solvers/fluid/fluid_runner.h"
#include "line/solvers/fluid/solver_fluid.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace env {

/**
 * What one stage solve reports back to the limits.
 *
 * The literal translation of what the reference's `[Qe,Ue,~,Te] = se.getAvg()`
 * plus `se.getAvgNode()` hand back, minus the node table: a Cache reports on the
 * NODE, so `cache` carries what `getAvgNode`'s resync would have written there.
 */
template <class T>
struct EnvStageAvg {
    Matrix<double> QN, UN, TN;
    solvers::CacheMetrics<T> cache;
};

/**
 * The stage solver, as a callable: the C++ spelling of the MATLAB function
 * handle `SolverENV(renv, @(m) SolverX(m, opts))` passes.
 *
 * Keeping it a `std::function` is what lets `env/` name a stage solver without
 * including one: this header sits ABOVE the MVA and NC runners in the include
 * graph (`env/solver_env.h` -> `ln/solver_ln.h` -> `mva/solver_mva.h`), so a
 * direct include would be a cycle.
 */
template <class T>
using EnvStageAvgFn = std::function<EnvStageAvg<T>(const qn::NetworkStruct<T>&)>;

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
    /**
     * The environment-blended cache surface, one entry per Cache node.
     *
     * EMPTY on a model with no cache. Within an entry every field is optional
     * and ABSENT MEANS NOT COMPUTED, never zero, exactly as `CacheNodeMetrics`
     * states: a class no stage measured stays NaN rather than reading as "never
     * hits". Kept `double`, since `init()` refuses `T != double` anyway.
     */
    solvers::CacheMetrics<double> cache;
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

    /**
     * With an injected stage solver: `stage_fn` runs every stage, exactly as the
     * reference calls its factory. `opt.stage_solver` is then not consulted,
     * since the caller has already named the solver by handing one over.
     */
    SolverEnvLimit(Environment<T>& e, const EnvOptions& o, EnvStageAvgFn<T> stage_fn)
        : envObj(e), opt(o), stage_fn_(stage_fn) {
        init();
    }

    EnvLimitSolution solve() { return dec_ ? solve_dec() : solve_avg(); }

private:
    void init() {
        if (opt.method != "avg" && opt.method != "dec")
            throw UnsupportedError("SolverENV limit: '" + opt.method +
                                   "' is not a closed-form environment limit; the two are 'avg' "
                                   "(fast environment) and 'dec' (slow environment)");
        dec_ = (opt.method == "dec");
        if (!stage_fn_ && opt.stage_solver != "fluid" && opt.stage_solver != "ctmc")
            throw UnsupportedError(
                "SolverENV limit: stage solver '" + opt.stage_solver +
                "' is not available; the limits solve each stage in STEADY STATE, by the fluid "
                "analyzer or the enumerated CTMC. Hand a stage solver to the EnvStageAvgFn "
                "constructor to run another one");
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
    }

public:
    /**
     * The RATE accumulator behind the cache blend, `SolverENV.accumCacheMetric`. Public, since
     * the mean-field coupling's `substituteCacheBlend_` weights by the same rule.
     *
     * Keyed by the Cache node's NAME, never by its position: the struct's node
     * order is not the declaration order, and indexing by position is what once
     * made a host write cache results onto a Sink (see `cache_metrics.h`).
     */
    class CacheAccum {
    public:
        /** One stage's contribution, weighted by its arrival into each cache. */
        void add(const qn::NetworkStruct<T>& sn, const EnvStageAvg<T>& stage, double w) {
            const std::vector<double> arv_all = cache_arrival(sn, stage.TN);
            for (std::size_t c = 0; c < stage.cache.caches.size(); ++c) {
                const solvers::CacheNodeMetrics<T>& cm = stage.cache.caches[c];
                Entry& en = entry(cm);
                const std::size_t K = std::max(cm.hitprob.size(), cm.missprob.size());
                std::vector<double> arv(std::max(K, en.arv.size()), 0.0);
                for (std::size_t r = 0; r < arv.size(); ++r) {
                    const std::size_t ind = cm.node;  // 1-based node index
                    const std::size_t off = (ind >= 1 ? (ind - 1) * sn.nclasses + r : arv_all.size());
                    if (off < arv_all.size()) arv[r] = w * arv_all[off];
                }
                grow(en.arv, arv.size());
                for (std::size_t r = 0; r < arv.size(); ++r) en.arv[r] += arv[r];
                accum(en.hit, cm.hitprob, arv);
                accum(en.miss, cm.missprob, arv);
                accum(en.dhit, cm.delayedprob, arv);
                accum_rows(en.hitl, cm.hitproblist, arv);
            }
        }

        /** Divide the accumulated rates by the accumulated arrivals, ONCE. */
        solvers::CacheMetrics<double> finish() const {
            solvers::CacheMetrics<double> out;
            for (std::size_t c = 0; c < order_.size(); ++c) {
                const Entry& en = order_[c];
                solvers::CacheNodeMetrics<double> cm;
                cm.node = en.node;
                cm.name = en.name;
                cm.itemcap = en.itemcap;
                cm.nitems = en.nitems;
                cm.hitprob = divide(en.hit, en.arv);
                cm.missprob = divide(en.miss, en.arv);
                cm.delayedprob = divide(en.dhit, en.arv);
                if (en.hitl.rows() > 0) {
                    cm.hitproblist = Matrix<double>(en.hitl.rows(), en.hitl.cols(), nan_());
                    for (std::size_t r = 0; r < en.hitl.rows(); ++r) {
                        const double d = r < en.arv.size() ? en.arv[r] : 0.0;
                        for (std::size_t l = 0; l < en.hitl.cols(); ++l)
                            cm.hitproblist(r, l) = d > 0.0 ? en.hitl(r, l) / d : nan_();
                    }
                }
                out.caches.push_back(cm);
            }
            return out;
        }

    private:
        struct Entry {
            std::size_t node = 0;
            std::string name;
            std::vector<double> itemcap;
            std::size_t nitems = 0;
            std::vector<double> arv;
            /** EMPTY stays empty: a quantity no stage measured must not become 0. */
            std::vector<double> hit, miss, dhit;
            Matrix<double> hitl;
        };

        static double nan_() { return std::numeric_limits<double>::quiet_NaN(); }

        static void grow(std::vector<double>& v, std::size_t n) {
            if (v.size() < n) v.resize(n, 0.0);
        }

        Entry& entry(const solvers::CacheNodeMetrics<T>& cm) {
            for (std::size_t i = 0; i < order_.size(); ++i)
                if (order_[i].name == cm.name) return order_[i];
            Entry en;
            en.node = cm.node;
            en.name = cm.name;
            en.itemcap = cm.itemcap;
            en.nitems = cm.nitems;
            order_.push_back(en);
            return order_.back();
        }

        /** One quantity: an ABSENT ratio contributes nothing and leaves it absent. */
        static void accum(std::vector<double>& a, const std::vector<T>& ratio,
                          const std::vector<double>& arv) {
            if (ratio.empty()) return;
            grow(a, ratio.size());
            for (std::size_t r = 0; r < ratio.size() && r < arv.size(); ++r) {
                const double v = num_traits<T>::to_double(ratio[r]);
                if (std::isfinite(v)) a[r] += arv[r] * v;
            }
        }

        static void accum_rows(Matrix<double>& a, const Matrix<T>& hl,
                               const std::vector<double>& arv) {
            if (hl.rows() == 0 || hl.cols() == 0) return;
            if (a.rows() == 0) a = Matrix<double>(hl.rows(), hl.cols(), 0.0);
            for (std::size_t r = 0; r < hl.rows() && r < a.rows(); ++r) {
                const double w = r < arv.size() ? arv[r] : 0.0;
                for (std::size_t l = 0; l < hl.cols() && l < a.cols(); ++l) {
                    const double v = num_traits<T>::to_double(hl(r, l));
                    if (std::isfinite(v)) a(r, l) += w * v;
                }
            }
        }

        /**
         * A class with no arrival anywhere divides by nothing and STAYS NaN: its
         * ratio is undefined, and 0 would read as "never hits".
         */
        static std::vector<double> divide(const std::vector<double>& a,
                                          const std::vector<double>& d) {
            std::vector<double> r;
            if (a.empty()) return r;
            r.assign(a.size(), nan_());
            for (std::size_t i = 0; i < a.size(); ++i)
                if (i < d.size() && d[i] > 0.0) r[i] = a[i] / d[i];
            return r;
        }

        /**
         * The per-class arrival INTO EACH NODE, flattened node-major.
         *
         * The weight is the cache's own arrival and not the Source throughput,
         * because the latter is 0 on a closed cache model. The station table
         * handed in is all zeros on purpose: only the Cache rows are read, and
         * those come from the chain reference throughput and the node visits,
         * never from a station row.
         */
        static std::vector<double> cache_arrival(const qn::NetworkStruct<T>& sn,
                                                 const Matrix<double>& TN) {
            const std::size_t I = sn.nodes.size(), R = sn.nclasses, M = sn.nstations;
            std::vector<double> out(I * R, 0.0);
            if (TN.rows() == 0) return out;
            Matrix<T> TNt(TN.rows(), TN.cols(), num_traits<T>::from_int(0));
            for (std::size_t i = 0; i < TN.rows(); ++i)
                for (std::size_t r = 0; r < TN.cols(); ++r)
                    TNt(i, r) = num_traits<T>::from_double(TN(i, r));
            const Matrix<T> AN(M, R, num_traits<T>::from_int(0));
            const Matrix<T> ANn = api::sn_get_node_arvr_from_tput<T>(sn, TNt, AN);
            if (ANn.rows() != I) return out;
            for (std::size_t ind = 0; ind < I; ++ind)
                for (std::size_t r = 0; r < R; ++r) {
                    const double v = num_traits<T>::to_double(ANn(ind, r));
                    out[ind * R + r] = std::isfinite(v) ? v : 0.0;
                }
            return out;
        }

        std::vector<Entry> order_;
    };

    /**
     * One stage in steady state under the enumerated CTMC at `cutoff`, with its cache split,
     * which is exact (read off the chain). An ensemble built on SolverCTMC stages.
     */
    static EnvStageAvg<T> ctmc_stage_steady(const qn::NetworkStruct<T>& sn, double cutoff) {
        ctmc::CtmcOptions co;
        co.cutoff = cutoff;
        const mva::AvgResult<T> r = ctmc::solver_ctmc_run_analyzer<T>(sn, co);
        EnvStageAvg<T> out;
        out.QN = to_double(r.QN);
        out.UN = to_double(r.UN);
        out.TN = to_double(r.TN);
        out.cache = r.cache;
        return out;
    }

private:
    /** A stage table widened to double, which is what the limits blend in. */
    static Matrix<double> to_double(const Matrix<T>& A) {
        Matrix<double> out(A.rows(), A.cols(), 0.0);
        for (std::size_t i = 0; i < A.rows(); ++i)
            for (std::size_t j = 0; j < A.cols(); ++j) out(i, j) = num_traits<T>::to_double(A(i, j));
        return out;
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
        CacheAccum acc;
        for (std::size_t e = 0; e < E; ++e) {
            const EnvStageAvg<T> s = stage_steady_state(envObj.stage(e).model);
            const double p = envObj.prob_env[e];
            acc.add(envObj.stage(e).model, s, p);
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
        out.cache = acc.finish();
        return out;
    }

    /** The fast-environment limit: one rate-averaged model, solved once. */
    EnvLimitSolution solve_avg() {
        EnvLimitSolution out;
        out.method = "avg";
        out.prob_env = envObj.prob_env;
        const qn::NetworkStruct<T> avg = rate_averaged_model();
        const EnvStageAvg<T> s = stage_steady_state(avg);
        // One model, so the weight is the arrival itself and the normalisation
        // returns exactly what the single solve reported.
        CacheAccum acc;
        acc.add(avg, s, 1.0);
        out.cache = acc.finish();
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

    /**
     * One stage, in steady state: the injected solver when there is one, else
     * the fluid analyzer, in the form that also REPORTS ITS CACHE SPLIT when the
     * stage holds a Cache.
     *
     * WHICH ENTRY POINT IS A MODEL QUESTION AND NOT A PREFERENCE, because the
     * two resolve `default` differently and therefore answer with DIFFERENT
     * CLOSURES. `solver_fluid` is `solver_fluid_analyzer.m`: it hands the method
     * to `fluid_dispatch`, which takes `default` to the matrix method (to
     * `closing` on a DPS model). `solver_fluid_run_analyzer` is the port of
     * `runAnalyzer`'s resolution and re-resolves `default` through the
     * minnormal/dae ladder. On a stage sitting at its saturation knee the two
     * disagree outright -- N/(Z+D) = 1/D is a continuum of first-order
     * equilibria, and which point a closure picks is the whole answer -- so
     * routing every stage through the runner silently moved the limits off the
     * matrix closure they are pinned against.
     *
     * A CACHE STAGE IS THE ONE THAT MUST TAKE THE RUNNER, and for the reference's
     * own reason rather than for a reporting convenience: `runAnalyzer` resolves
     * a Cache model to `rmf`, which lives in `fluid_cacheqn.h` and which
     * `fluid_dispatch` refuses BY NAME (it cannot call it without a cyclic
     * include). That route is also the one that reports the split, which is the
     * silence this file used to refuse a cache model over.
     */
    EnvStageAvg<T> stage_steady_state(const qn::NetworkStruct<T>& sn) const {
        if (stage_fn_) return stage_fn_(sn);
        if (opt.stage_solver == "ctmc") return ctmc_stage_steady(sn, opt.stage_cutoff);
        fluid::FluidOptions fo = opt.stage;
        // A limit starts from nothing carried over -- there is no entry state to
        // seed with, which is the whole content of the approximation.
        fo.init_sol.clear();
        EnvStageAvg<T> out;
        const fluid::FluidSolution s =
            fluid::detail::fluid_has_cache(sn)
                ? fluid::solver_fluid_run_analyzer<T>(sn, fo, nullptr, nullptr, &out.cache)
                : fluid::solver_fluid<T>(sn, fo);
        out.QN = s.QN;
        out.UN = s.UN;
        out.TN = s.TN;
        return out;
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
    EnvStageAvgFn<T> stage_fn_;
    std::size_t M = 0, K = 0;
    bool dec_ = false;
};

/** `solveEnvLimit` on the original stages. */
template <class T>
EnvLimitSolution solver_env_limit(Environment<T>& e, const EnvOptions& o) {
    SolverEnvLimit<T> s(e, o);
    return s.solve();
}

/** `solveEnvLimit` with the CALLING solver running every stage. */
template <class T>
EnvLimitSolution solver_env_limit(Environment<T>& e, const EnvOptions& o,
                                  EnvStageAvgFn<T> stage_fn) {
    SolverEnvLimit<T> s(e, o, stage_fn);
    return s.solve();
}

}  // namespace env
}  // namespace line

#endif  // LINE_SOLVERS_ENV_SOLVER_ENV_LIMIT_H
