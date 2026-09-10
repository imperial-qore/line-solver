/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_ENV_SOLVER_ENV_MEANFIELD_H
#define LINE_SOLVERS_ENV_SOLVER_ENV_MEANFIELD_H

/**
 * SolverENV, the default mean-field path: ENVIRONMENT COMPRESSION.
 *
 * WHAT THIS FILE IS, and what it deliberately is not. The mean-field coupling
 * itself -- `solver_env_meanfield_analyzer.m`, the `pre_`/`analyze_`/`post_`/
 * `finish_`/`converged_` cycle over per-stage transient means -- is ALREADY
 * ported, in `solver_env.h`. What was missing is the other half of the
 * mean-field path, the half that lives in `@@SolverENV/SolverENV.m` rather than
 * in the analyzer file: the environment rate matrix `E0`, its generator
 * `Eutil`, `ctmc_decompose`, `findBestPartition`, `beamSearchPartition`,
 * `computeMacroRate` and `applyCompression`. That is what this header adds, and
 * it composes with `SolverEnv<T>` instead of duplicating it.
 *
 * WHAT COMPRESSION BUYS. The mean-field fixed point costs one transient stage
 * solve per stage per iteration, so an environment with many stages is
 * expensive in the number of stages and not in their size. When the environment
 * is NEARLY COMPLETELY DECOMPOSABLE -- stages fall into groups that switch
 * rapidly among themselves and rarely across groups -- each group behaves, on
 * the slow time scale the network actually feels, like a single stage running
 * at the group's conditionally-averaged rates. Aggregating the groups replaces
 * E stage solves by E' < E of them, and the error is governed by the degree of
 * coupling eps against the admissible epsMAX that the kernels report.
 *
 * THE KERNELS ARE NOT RE-DERIVED HERE. `ctmc_courtois`, `ctmc_kms`,
 * `ctmc_takahashi` and `ctmc_multi` are already ported under `line/api/mc/`,
 * with their reference-defect histories recorded in their own headers; this
 * file is the dispatch `SolverENV.ctmc_decompose` puts in front of them and
 * nothing more. The same discipline applies downstream: the metrics still come
 * out of `SolverEnv<T>`, which is the single mean-field implementation, so
 * there is no second copy of the Util computation to drift -- the failure mode
 * `_kb/06-solver-catalog.md` records under "ENV state-vector Util had drifted
 * from the CTMC analyzer".
 *
 * ONE DELIBERATE DIVERGENCE FROM THE REFERENCE, and it is a correctness fix
 * rather than a preference. `applyCompression` replaces `self.ensemble`,
 * `self.solvers` and `self.sn` with their E' macro versions but leaves
 * `self.envObj` at its original E stages, so `envObj.proc{e}{h}` and
 * `envObj.holdTime{e}` still describe MICRO transitions while every consumer
 * now indexes them as macro ones. The analyzer's `post_` and `finish_` weight
 * the macro transients by those stale micro CDFs, silently, and only when
 * E' == E do the two agree. Here `env_compress` builds a genuinely compressed
 * `Environment<T>`: E' stages, and arcs carrying the aggregated macro rates, so
 * the holding times the analyzer integrates against are the ones its stages
 * actually have.
 *
 * WHAT IS PORTED, and what is refused by name:
 *   ported   `E0`/`Eutil`, `ctmc_decompose` over all four kernels,
 *            `findBestPartition`, `beamSearchPartition`, `computeMacroRate`,
 *            `applyCompression`, and the mean-field solve on top of the result
 *   ported   `aggregateCacheMeanfield_` too, as of 2026-08-15: it is a SECOND
 *            mean-field fixed point, nested beside the queue-length one and
 *            carrying each Cache's own occupancy across a switch by prob_orig,
 *            because a cache's state is a per-item occupancy that no marginal
 *            queue length encodes. Its per-sweep transient is
 *            `solver_fld_cacheqn_tran`, whose grid runs in PER-REQUEST time --
 *            it is divided by the cache's total arrival rate before the
 *            holding-time CDF is evaluated on it, or every stage is weighted as
 *            though its cache saw the same request rate
 *   refused  a Cache under a NON-FLUID stage solver (only the fluid stage
 *            exposes that RMF transient), compression of a non-exponential
 *            environment, and everything `solver_env.h` already refuses by name
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "line/api/mc/ctmc_courtois.h"
#include "line/api/mc/ctmc_kms.h"
#include "line/api/mc/ctmc_multi.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/api/mc/ctmc_takahashi.h"
#include "line/lang/qn/environment.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/api/mam/map_cdf.h"
#include "line/api/mam/map_moment.h"
#include "line/solvers/env/solver_env.h"
#include "line/solvers/fluid/fluid_cacheqn.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace env {

/** A partition of the stage indices 0..E-1 into macro-states, MATLAB's `MS`. */
using MacroPartition = std::vector<std::vector<std::size_t>>;

/** The knobs `applyCompression` reads out of `options.config`. */
struct EnvCompressOptions {
    /** `options.config.da`: courtois (the reference's default), kms, takahashi, multi. */
    std::string da = "courtois";
    /** `options.config.da_iter`: sweeps for the two iterative kernels. */
    std::size_t da_iter = 10;
    /** `options.config.env_alpha`: the beam search's per-depth merge penalty. */
    double env_alpha = 0.01;
    /** Beam width, the reference's hard-coded B = 3. */
    std::size_t beam_width = 3;
    /**
     * Stage count above which the reference switches from the pairwise search
     * to the beam search. The two are genuinely different searches, not one
     * search with a budget, so the threshold changes the answer and is exposed.
     */
    std::size_t beam_above_stages = 10;
    /**
     * A partition supplied by the caller, which SKIPS the search entirely.
     * Worth having because both searches evaluate a decomposition per candidate
     * pair, and a caller who already knows the group structure -- a repair model
     * whose stages are (working, degraded) times (peak, offpeak), say -- should
     * not pay for rediscovering it.
     */
    MacroPartition partition;
};

/** What `SolverENV.ctmc_decompose` returns: `[p, eps, epsMax, q]`. */
template <class T>
struct EnvDecomp {
    std::vector<T> p;  ///< approximate stationary vector of the environment
    T eps, epsMAX, q;
};

namespace meanfield_detail {

/** Every kernel is seeded by Courtois, whose epsMAX is an eigenvalue modulus. */
inline void refuse_inexact(const std::string& da) {
    throw UnsupportedError(
        "SolverENV compression: the '" + da +
        "' decomposition needs transcendental arithmetic, because every one of ctmc_courtois, "
        "ctmc_kms, ctmc_takahashi and ctmc_multi reports epsMAX, a subdominant eigenvalue "
        "modulus with no rational closed form; rerun with --arith double or real");
}

/** MS must be a partition of 0..n-1, which every kernel assumes and none checks. */
inline void check_partition(const MacroPartition& MS, std::size_t n) {
    std::vector<bool> seen(n, false);
    std::size_t total = 0;
    for (const std::vector<std::size_t>& blk : MS) {
        if (blk.empty()) throw InputError("SolverENV compression: a macro-state is empty");
        for (std::size_t s : blk) {
            if (s >= n)
                throw InputError("SolverENV compression: a macro-state names stage " +
                                 std::to_string(s + 1) + ", which does not exist");
            if (seen[s])
                throw InputError("SolverENV compression: stage " + std::to_string(s + 1) +
                                 " appears in more than one macro-state");
            seen[s] = true;
            ++total;
        }
    }
    if (total != n)
        throw InputError(
            "SolverENV compression: the macro-states must cover every stage; the partition "
            "covers " +
            std::to_string(total) + " of " + std::to_string(n));
}

/** Singletons, the starting point of both searches and the no-compression fallback. */
inline MacroPartition singletons(std::size_t E) {
    MacroPartition MS(E);
    for (std::size_t i = 0; i < E; ++i) MS[i] = std::vector<std::size_t>{i};
    return MS;
}

/** `trial`: merge blocks i and j of `ms`, the merged block taking i's position. */
inline MacroPartition merge_blocks(const MacroPartition& ms, std::size_t i, std::size_t j) {
    MacroPartition out;
    out.reserve(ms.size() - 1);
    for (std::size_t k = 0; k < ms.size(); ++k) {
        if (k == i) {
            std::vector<std::size_t> m = ms[i];
            m.insert(m.end(), ms[j].begin(), ms[j].end());
            out.push_back(m);
        } else if (k != j) {
            out.push_back(ms[k]);
        }
    }
    return out;
}

}  // namespace meanfield_detail

/**
 * Port of `SolverENV.ctmc_decompose`: one NCD decomposition, by whichever
 * kernel `options.config.da` names.
 *
 * The uniformization rate the three iterative kernels report is the reference's
 * `1.05 * max(max(abs(Q)))` rather than anything the kernel itself derived,
 * which is why it is recomputed here instead of read back off the result.
 */
template <class T>
EnvDecomp<T> env_ctmc_decompose(const Matrix<T>& Q, const MacroPartition& MS,
                                const EnvCompressOptions& opt) {
    meanfield_detail::check_partition(MS, Q.rows());
    EnvDecomp<T> out;
    if constexpr (!num_traits<T>::has_transcendental) {
        meanfield_detail::refuse_inexact(opt.da);
    } else {
        // 1.05 max|Q|, the rate the reference hands back for every kernel that
        // does not report one of its own.
        T qmax = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < Q.rows(); ++i)
            for (std::size_t j = 0; j < Q.cols(); ++j) {
                const T a = num_abs(T(Q(i, j)));
                if (a > qmax) qmax = a;
            }
        const T qdefault = T(num_traits<T>::from_rational(21, 20) * qmax);

        if (opt.da == "courtois") {
            const mc::CourtoisResult<T> r = mc::ctmc_courtois(Q, MS);
            out.p = r.p;
            out.eps = r.eps;
            out.epsMAX = r.epsMAX;
            out.q = r.q;
        } else if (opt.da == "kms") {
            const mc::KmsResult<T> r = mc::ctmc_kms(Q, MS, opt.da_iter);
            out.p = r.p;
            out.eps = r.eps;
            out.epsMAX = r.epsMAX;
            out.q = qdefault;
        } else if (opt.da == "takahashi") {
            const mc::TakahashiResult<T> r = mc::ctmc_takahashi(Q, MS, opt.da_iter);
            out.p = r.p;
            out.eps = r.eps;
            out.epsMAX = r.epsMAX;
            out.q = qdefault;
        } else if (opt.da == "multi") {
            // The coarse partition defaults to singletons over the macro-states,
            // so the second level decouples nothing: the reference exposes no
            // way to supply a real macro-macro partition, and a two-level method
            // whose coarse level is singletons is Courtois plus one extra solve.
            MacroPartition MSS(MS.size());
            for (std::size_t i = 0; i < MS.size(); ++i) MSS[i] = std::vector<std::size_t>{i};
            const mc::MultiResult<T> r = mc::ctmc_multi(Q, MS, MSS);
            out.p = r.p;
            out.eps = r.eps;
            out.epsMAX = r.epsMAX;
            out.q = qdefault;
        } else {
            throw UnsupportedError(
                "SolverENV compression: unknown decomposition '" + opt.da +
                "'; options.config.da is one of courtois, kms, takahashi, multi");
        }
    }
    return out;
}

/**
 * `E0`, the environment's rate matrix: `E0(e,h) = env{e,h}.getRate()`.
 *
 * `getRate()` is the RECIPROCAL MEAN of the transition, so a general Markovian
 * arc collapses to a single rate here and everything downstream treats the
 * environment as a CTMC. That is the reference's own reading and it is why
 * `env_compress` refuses a non-exponential environment by name: the collapse is
 * harmless for the NCD diagnostics, which only ever look at Eutil, but it is
 * not harmless once the macro arcs are rebuilt from it.
 */
template <class T>
Matrix<T> env_rate_matrix(const Environment<T>& e) {
    const std::size_t E = e.nstages();
    Matrix<T> E0(E, E, num_traits<T>::from_int(0));
    for (std::size_t a = 0; a < E; ++a)
        for (std::size_t b = 0; b < E; ++b)
            if (e.arc(a, b).enabled) E0(a, b) = e.arc(a, b).dist.rate();
    return E0;
}

/**
 * Port of `findBestPartition`, the small-environment search.
 *
 * ITS COMMENT CLAIMS AN EXHAUSTIVE SEARCH OVER ALL PARTITIONS AND THE CODE DOES
 * NOT DO THAT. It evaluates the singletons and then every single pairwise merge
 * of them, so it explores E(E-1)/2 + 1 partitions out of the Bell number of
 * them and can never return a macro-state of more than two stages. The port is
 * literal, because the alternative is a different method wearing the reference's
 * name; a caller who wants deeper merging has `beam_above_stages` and
 * `EnvCompressOptions::partition`.
 */
template <class T>
MacroPartition env_find_best_partition(const Matrix<T>& Eutil, const EnvCompressOptions& opt) {
    const std::size_t E = Eutil.rows();
    MacroPartition best = meanfield_detail::singletons(E);
    EnvDecomp<T> b = env_ctmc_decompose(Eutil, best, opt);
    double best_eps = num_traits<T>::to_double(b.eps);
    if (std::isnan(best_eps)) return best;

    for (std::size_t i = 0; i < E; ++i)
        for (std::size_t j = i + 1; j < E; ++j) {
            const MacroPartition trial =
                meanfield_detail::merge_blocks(meanfield_detail::singletons(E), i, j);
            const EnvDecomp<T> t = env_ctmc_decompose(Eutil, trial, opt);
            const double te = num_traits<T>::to_double(t.eps);
            if (!std::isnan(te) && te < best_eps) {
                best_eps = te;
                best = trial;
            }
        }
    return best;
}

/**
 * Port of `beamSearchPartition`, the large-environment search: repeatedly merge
 * two blocks, keeping the `beam_width` cheapest partitions at each depth.
 *
 * THE COST AND THE INCUMBENT ARE NOT THE SAME QUANTITY, in the reference. The
 * incumbent `bestEps` is seeded with the raw eps of the singleton partition,
 * and thereafter compared against `childEps - childEpsMax + alpha * depth`,
 * which is a penalized score and not an eps at all. A merge is therefore
 * adopted partly on the strength of its epsMAX and of how deep it sits, against
 * a threshold that measured neither. This is ported literally rather than
 * repaired: the search is a heuristic whose output is checked afterwards
 * against eps <= epsMAX, so the comparison decides which candidate is tried and
 * not whether the result is admissible.
 */
template <class T>
MacroPartition env_beam_search_partition(const Matrix<T>& Eutil, const EnvCompressOptions& opt) {
    const std::size_t E = Eutil.rows();
    std::vector<MacroPartition> beam{meanfield_detail::singletons(E)};
    MacroPartition best = beam[0];
    double best_cost = num_traits<T>::to_double(env_ctmc_decompose(Eutil, best, opt).eps);

    for (std::size_t depth = 1; depth < E; ++depth) {
        std::vector<std::pair<double, MacroPartition>> cand;
        for (const MacroPartition& ms : beam) {
            for (std::size_t i = 0; i < ms.size(); ++i)
                for (std::size_t j = i + 1; j < ms.size(); ++j) {
                    const MacroPartition trial = meanfield_detail::merge_blocks(ms, i, j);
                    const EnvDecomp<T> t = env_ctmc_decompose(Eutil, trial, opt);
                    const double te = num_traits<T>::to_double(t.eps);
                    if (std::isnan(te) || !(te > 0.0)) continue;
                    const double cost = te - num_traits<T>::to_double(t.epsMAX) +
                                        opt.env_alpha * static_cast<double>(depth);
                    cand.push_back(std::make_pair(cost, trial));
                    if (cost < best_cost) {
                        best_cost = cost;
                        best = trial;
                    }
                }
        }
        if (cand.empty()) break;
        // Stable, so that ties keep the order the merges were generated in and
        // the search is reproducible across runs.
        std::stable_sort(cand.begin(), cand.end(),
                         [](const std::pair<double, MacroPartition>& a,
                            const std::pair<double, MacroPartition>& b) { return a.first < b.first; });
        beam.clear();
        for (std::size_t i = 0; i < opt.beam_width && i < cand.size(); ++i)
            beam.push_back(cand[i].second);
    }
    return best;
}

/** Everything `applyCompression` computes, plus the compressed environment. */
template <class T>
struct EnvCompression {
    MacroPartition MS;
    /**
     * The compressed environment. HELD BY SHARED POINTER because `SolverEnv<T>`
     * stores a reference to the environment it solves, so the compressed one has
     * to outlive the solver; returning it by value would make that the caller's
     * problem to get right, and getting it wrong is a dangling reference rather
     * than a wrong number.
     */
    std::shared_ptr<Environment<T>> env;
    Matrix<T> E0, Eutil;
    Matrix<T> macro_rate;    ///< `computeMacroRate(i,j)`
    std::vector<T> p;        ///< micro stationary vector from the decomposition
    std::vector<T> pmicro;   ///< within-macro-state conditional probabilities
    std::vector<T> pmacro;   ///< `pMacro`, the macro-state probabilities
    Matrix<double> prob_orig;  ///< the macro embedding weights `newEmbweight`
    T eps, epsMAX, q;
    /** eps <= epsMAX: below this the aggregation is meaningful, above it is not. */
    bool compressible = false;
};

/**
 * Port of `applyCompression`: pick a partition, decompose, and build the
 * macro-state environment.
 *
 * WHY THE MACRO SERVICE RATES ARE A pmicro-WEIGHTED AVERAGE. Within a
 * macro-state the environment switches fast compared with the network, so the
 * network sees the group's rates averaged over the CONDITIONAL distribution of
 * being in each micro-stage given the group -- which is exactly pmicro. That
 * average is over rates and not over distributions, so a phase-type service
 * collapses to an exponential of the same mean: the compression keeps the first
 * moment and discards the SCV, as the reference's `Exp(rateSum)` does.
 */
template <class T>
EnvCompression<T> env_compress(const Environment<T>& e0, const EnvCompressOptions& opt) {
    const std::size_t E = e0.nstages();
    const T zero = num_traits<T>::from_int(0);

    // A MACRO-STATE IS A NETWORK AT AVERAGED RATES, so every stage merged into
    // one has to have a station rate table to average. A layered stage does not:
    // its stations are the layers SolverLN derives from it, and averaging those
    // would aggregate an artifact of the layering rather than the model.
    e0.reject_lqn_stages(
        "SolverENV compression",
        "a macro-state is built as one stage network carrying the pmicro-weighted average of "
        "its members' station rates, and a layered model has no such rate table -- only the "
        "layers SolverLN derives from it");

    // The whole construction reads the environment as a CTMC, so a transition
    // that is not exponential cannot survive it: the macro arc would be built as
    // an Exp of the aggregated rate, and the analyzer would then integrate the
    // stage transient against a holding-time CDF the model never had.
    for (std::size_t a = 0; a < E; ++a)
        for (std::size_t b = 0; b < E; ++b) {
            if (!e0.arc(a, b).enabled) continue;
            if (e0.arc(a, b).dist.type != lang::ProcessType::EXP)
                throw UnsupportedError(
                    "SolverENV compression: the transition from stage " + std::to_string(a + 1) +
                    " to " + std::to_string(b + 1) +
                    " is not exponential, and the NCD decomposition reads the environment as a "
                    "CTMC through E0 = getRate(); aggregating it would silently replace the "
                    "transition by an exponential of the same mean, so it is refused instead");
        }

    EnvCompression<T> c;
    c.E0 = env_rate_matrix(e0);
    c.Eutil = mc::ctmc_makeinfgen(c.E0);

    if (!opt.partition.empty()) {
        meanfield_detail::check_partition(opt.partition, E);
        c.MS = opt.partition;
    } else if (E <= opt.beam_above_stages) {
        c.MS = env_find_best_partition(c.Eutil, opt);
    } else {
        c.MS = env_beam_search_partition(c.Eutil, opt);
    }
    const std::size_t Ec = c.MS.size();

    const EnvDecomp<T> d = env_ctmc_decompose(c.Eutil, c.MS, opt);
    c.p = d.p;
    c.eps = d.eps;
    c.epsMAX = d.epsMAX;
    c.q = d.q;
    // The reference warns and continues. The flag is reported rather than
    // thrown for the same reason: an environment that does not decompose still
    // has an answer, it is simply the answer to a model the aggregation moved.
    c.compressible = num_traits<T>::to_double(d.eps) <= num_traits<T>::to_double(d.epsMAX);

    c.pmacro.assign(Ec, zero);
    for (std::size_t i = 0; i < Ec; ++i)
        for (std::size_t s : c.MS[i]) c.pmacro[i] += c.p[s];
    c.pmicro.assign(E, zero);
    for (std::size_t i = 0; i < Ec; ++i) {
        if (num_traits<T>::to_double(c.pmacro[i]) <= 0) continue;
        for (std::size_t s : c.MS[i]) c.pmicro[s] = T(c.p[s] / c.pmacro[i]);
    }

    // computeMacroRate: the micro rates out of the block, weighted by the
    // conditional probability of sitting in each of its micro-stages.
    c.macro_rate = Matrix<T>(Ec, Ec, zero);
    for (std::size_t i = 0; i < Ec; ++i)
        for (std::size_t j = 0; j < Ec; ++j)
            for (std::size_t mi : c.MS[i])
                for (std::size_t mj : c.MS[j]) c.macro_rate(i, j) += T(c.pmicro[mi] * c.E0(mi, mj));

    // newEmbweight: P(the previous macro-state was k | now entering e).
    c.prob_orig = Matrix<double>(Ec, Ec, 0.0);
    for (std::size_t x = 0; x < Ec; ++x) {
        double tot = 0.0;
        for (std::size_t h = 0; h < Ec; ++h)
            if (h != x)
                tot += num_traits<T>::to_double(c.pmacro[h]) *
                       num_traits<T>::to_double(c.macro_rate(h, x));
        if (!(tot > 0.0)) continue;
        for (std::size_t k = 0; k < Ec; ++k) {
            if (k == x) continue;
            c.prob_orig(k, x) = num_traits<T>::to_double(c.pmacro[k]) *
                                num_traits<T>::to_double(c.macro_rate(k, x)) / tot;
        }
    }

    // The macro networks: the first micro-stage's structure, its rates replaced
    // by the pmicro-weighted averages over the block.
    c.env = std::make_shared<Environment<T>>(e0.name() + "-compressed", Ec);
    for (std::size_t i = 0; i < Ec; ++i) {
        const std::size_t first = c.MS[i][0];
        qn::NetworkStruct<T> sn = e0.stage(first).model;
        const std::size_t M = sn.nstations, K = sn.nclasses;
        for (std::size_t s : c.MS[i])
            if (e0.stage(s).model.nstations != M || e0.stage(s).model.nclasses != K)
                throw InputError(
                    "SolverENV compression: macro-state " + std::to_string(i + 1) +
                    " merges stages with different stations or classes, whose rates cannot be "
                    "averaged entrywise");
        for (std::size_t m = 0; m < M; ++m) {
            const lang::NodeType nt = sn.stations[m].nodetype;
            // Only a Queue or a Delay has a service rate to average; a Source
            // carries an arrival process and a Join an infinite rate, and the
            // reference skips both.
            if (nt != lang::NodeType::Queue && nt != lang::NodeType::Delay) continue;
            for (std::size_t k = 0; k < K; ++k) {
                T acc = zero;
                for (std::size_t s : c.MS[i])
                    acc += T(c.pmicro[s] * e0.stage(s).model.rates(m, k));
                if (num_traits<T>::to_double(acc) > 0) sn.service[m][k] = lang::Distrib<T>::exp_rate(acc);
            }
        }
        sn.refresh_struct();
        c.env->set_stage(i, e0.stage(first).name + "+", e0.stage(first).type, sn);
    }

    for (std::size_t i = 0; i < Ec; ++i) {
        bool any_out = false;
        for (std::size_t j = 0; j < Ec; ++j) {
            if (i == j) continue;  // a self-loop would lengthen the holding time
            if (!(num_traits<T>::to_double(c.macro_rate(i, j)) > 0)) continue;
            // The reset policy of the representative micro arc. Two micro arcs
            // folding into one macro arc may carry DIFFERENT resets and a
            // std::function cannot be compared, so disagreement is detectable
            // only in whether a reset is present at all; that much is refused,
            // and beyond it the representative stands.
            const std::size_t fi = c.MS[i][0], fj = c.MS[j][0];
            const bool want = static_cast<bool>(e0.arc(fi, fj).reset);
            for (std::size_t a : c.MS[i])
                for (std::size_t b : c.MS[j])
                    if (e0.arc(a, b).enabled && static_cast<bool>(e0.arc(a, b).reset) != want)
                        throw UnsupportedError(
                            "SolverENV compression: the arcs folding into the macro transition " +
                            std::to_string(i + 1) + " -> " + std::to_string(j + 1) +
                            " do not agree on whether a reset policy applies, and one macro arc "
                            "can carry only one; split the partition so that reset policies are "
                            "uniform within it");
            c.env->add_transition(i, j, lang::Distrib<T>::exp_rate(c.macro_rate(i, j)),
                                  e0.arc(fi, fj).reset);
            any_out = true;
        }
        if (!any_out)
            throw InputError(
                "SolverENV compression: macro-state " + std::to_string(i + 1) +
                " has no outgoing transition, so the compressed environment is absorbing; the "
                "partition merged a whole recurrent class into one block");
    }
    return c;
}

/**
 * `probEnv = pMacro` and `probOrig = newEmbweight`, the two quantities
 * `applyCompression` overwrites on the environment.
 *
 * ORDER MATTERS: `SolverEnv<T>`'s constructor calls `Environment::init()`,
 * which recomputes both from the macro arcs, so this must be applied AFTER the
 * solver is constructed and BEFORE `solve()` is called. `solver_env_meanfield`
 * below does exactly that, and is the reason to prefer it over wiring the two
 * calls by hand.
 *
 * The two are consistent rather than contradictory: for an exponential
 * environment `Environment::init()` derives probEnv as the stationary law of
 * the macro generator, and aggregating a chain by its exact conditional
 * distributions reproduces the block sums of the original stationary law
 * exactly. So this overwrite replaces one estimate of the same quantity by
 * another, and the gap between them is a second reading of the decomposition
 * error alongside eps.
 */
template <class T>
void env_apply_macro_probabilities(Environment<T>& e, const EnvCompression<T>& c) {
    const std::size_t Ec = c.MS.size();
    if (e.nstages() != Ec)
        throw InputError(
            "SolverENV compression: the macro probabilities do not match the environment they "
            "are being applied to");
    e.prob_env.assign(Ec, 0.0);
    for (std::size_t i = 0; i < Ec; ++i) e.prob_env[i] = num_traits<T>::to_double(c.pmacro[i]);
    e.prob_orig = c.prob_orig;
}

/** A mean-field solve, with the compression that produced it. */
/** `aggregateCacheMeanfield_`'s output: one hit/miss vector per Cache node. */
template <class T>
struct CacheBlendResult {
    std::vector<std::size_t> nodes;         ///< 0-based Cache node indices of stage 1
    std::vector<std::vector<T> > hitprob;   ///< [cache][class], NaN where no flow
    std::vector<std::vector<T> > missprob;  ///< [cache][class]
};

template <class T>
struct EnvMeanfieldSolution {
    EnvSolution avg;              ///< what SolverEnv reported
    EnvCompression<T> compression;
    bool compressed = false;      ///< false when the solve ran on the original stages
    /**
     * `aggregateCacheMeanfield_`: environment-blended hit and miss probabilities
     * per Cache node. Empty when the model holds no Cache. The reference writes
     * these onto the stage-one node objects with `setResultHitProb`; this port
     * has no model-object layer at this level, so they ride in the result.
     */
    CacheBlendResult<T> cache;
};

namespace meanfield_detail {

/** 0-based Cache node indices of a stage, in `find(nodetype == Cache)` order. */
template <class T>
std::vector<std::size_t> cache_nodes_of(const qn::NetworkStruct<T>& sn) {
    std::vector<std::size_t> out;
    for (std::size_t ind = 0; ind < sn.nodes.size(); ++ind)
        if (sn.nodes[ind].nodetype == lang::NodeType::Cache) out.push_back(ind);
    return out;
}

/**
 * `aggregateCacheMeanfield_`: the environment-blended hit and miss ratios of
 * every Cache, by a mean-field fixed point over the caches' OWN occupancies.
 *
 * THIS IS A SECOND FIXED POINT, nested beside the queue-length one, and it
 * exists because the two carry different objects. `SolverEnv`'s coupling hands
 * each stage the marginal mean queue lengths its predecessors left; a cache's
 * state is its per-item occupancy, which no queue length encodes. So this sweep
 * integrates each stage's cache drift over its sojourn from an entry occupancy
 * that mixes its predecessors' exit occupancies by `prob_orig`, exactly the way
 * the queue-length handoff mixes means, and iterates the pair to convergence.
 *
 * THE RMF DRIFT RUNS IN PER-REQUEST TIME, NOT REAL TIME, which is the trap here.
 * `solver_fld_cacheqn_tran` returns a grid in units of cache requests, so the
 * holding-time CDF cannot be evaluated on it directly: the grid is divided by
 * the cache's total arrival rate first (`treal = t / Lam`). Skipping that
 * weights every stage as though its cache saw the same request rate, which
 * silently favours the slow stages.
 *
 * THE RATIO IS TAKEN AFTER THE BLEND, as in the state-vector twin: a per-stage
 * ratio weighted by prob_env averages ratios, which is not the ratio the
 * environment exhibits unless every stage carries the same total rate.
 */
template <class T>
CacheBlendResult<T> aggregate_cache_meanfield(Environment<T>& e, const EnvOptions& o) {
    CacheBlendResult<T> out;
    const std::size_t E = e.nstages();
    if (E == 0) return out;
    // The blend is over the Cache NODES of the stage networks; a layered
    // environment has none, and reading an absent NetworkStruct would report an
    // empty node list as though the model had been examined.
    if (e.has_lqn_stages()) return out;
    const qn::NetworkStruct<T>& sn1 = e.stage(0).model;
    out.nodes = cache_nodes_of(sn1);
    if (out.nodes.empty()) return out;
    const std::size_t nc = out.nodes.size(), K = sn1.nclasses;

    // Only a FLUID stage exposes the RMF cache transient this reads; the
    // reference returns without writing anything when any stage is not one,
    // rather than blending what it has.
    if (o.stage_solver != "fluid") {
        out.nodes.clear();
        return out;
    }

    // A finite window per stage, falling back to a few mean holding times when
    // the inner solver left the timespan open, as the reference does.
    std::vector<double> tend(E, 0.0);
    for (std::size_t s = 0; s < E; ++s) {
        double t1 = o.timespan_end;
        if (!std::isfinite(t1) || !(t1 > 0.0)) t1 = 20.0 * mam::map_mean(e.hold_time[s].map());
        tend[s] = t1;
    }

    std::vector<std::vector<std::vector<T> > > entry(E, std::vector<std::vector<T> >(nc));
    std::vector<std::vector<fluid::FluidCacheqnTranCache<T> > > tran(E);
    std::vector<std::vector<std::vector<T> > > wmass(E, std::vector<std::vector<T> >(nc));
    std::vector<std::vector<std::vector<T> > > exit_occ(E, std::vector<std::vector<T> >(nc));
    std::vector<T> prev_flat;
    const int sweeps = o.iter_max > 0 ? o.iter_max : 1;

    for (int sweep = 0; sweep < sweeps; ++sweep) {
        for (std::size_t s = 0; s < E; ++s) {
            fluid::FluidOptions fo;
            fo.method = "rmf";
            tran[s] = fluid::solver_fld_cacheqn_tran(e.stage(s).model, fo, 0.0, tend[s], entry[s]);
            for (std::size_t c = 0; c < nc; ++c) {
                std::size_t idx = tran[s].size();
                for (std::size_t q = 0; q < tran[s].size(); ++q)
                    if (tran[s][q].node == out.nodes[c]) idx = q;
                if (idx == tran[s].size()) continue;
                const fluid::FluidCacheqnTranCache<T>& tc = tran[s][idx];
                double Lam = 0.0;
                for (std::size_t k = 0; k < tc.arate.size(); ++k)
                    Lam += num_traits<T>::to_double(tc.arate[k]);
                if (!(Lam > 0.0)) continue;

                // PER-REQUEST TIME -> REAL TIME before the CDF is evaluated.
                std::vector<double> treal(tc.t.size(), 0.0);
                for (std::size_t j = 0; j < tc.t.size(); ++j)
                    treal[j] = num_traits<T>::to_double(tc.t[j]) / Lam;
                const std::vector<double> F = mam::map_cdf(e.hold_time[s].map(), treal);
                std::vector<T> w(treal.size(), num_traits<T>::from_int(0));
                for (std::size_t j = 1; j < treal.size(); ++j)
                    w[j] = num_traits<T>::from_double(F[j] - F[j - 1]);
                wmass[s][c] = w;

                double sw = 0.0;
                for (std::size_t j = 0; j < w.size(); ++j) sw += num_traits<T>::to_double(w[j]);
                if (!(sw > 0.0) || tc.xocc.rows() == 0) continue;
                std::vector<T> xo(tc.xocc.rows(), num_traits<T>::from_int(0));
                for (std::size_t a = 0; a < tc.xocc.rows(); ++a) {
                    T acc = num_traits<T>::from_int(0);
                    for (std::size_t j = 0; j < w.size() && j < tc.xocc.cols(); ++j)
                        acc = T(acc + T(tc.xocc(a, j) * w[j]));
                    xo[a] = T(acc / num_traits<T>::from_double(sw));
                }
                exit_occ[s][c] = xo;
            }
        }

        // Each stage's entry occupancy is its predecessors' exits, by prob_orig.
        std::vector<std::vector<std::vector<T> > > next(E, std::vector<std::vector<T> >(nc));
        for (std::size_t s = 0; s < E; ++s)
            for (std::size_t c = 0; c < nc; ++c) {
                std::vector<T> acc;
                for (std::size_t h = 0; h < E; ++h) {
                    const double po = e.prob_orig(h, s);
                    if (!(po > 0.0) || exit_occ[h][c].empty()) continue;
                    if (acc.empty()) acc.assign(exit_occ[h][c].size(), num_traits<T>::from_int(0));
                    if (acc.size() != exit_occ[h][c].size()) continue;
                    for (std::size_t a = 0; a < acc.size(); ++a)
                        acc[a] = T(acc[a] + T(num_traits<T>::from_double(po) * exit_occ[h][c][a]));
                }
                next[s][c] = acc;
            }

        std::vector<T> flat;
        for (std::size_t s = 0; s < E; ++s)
            for (std::size_t c = 0; c < nc; ++c)
                flat.insert(flat.end(), next[s][c].begin(), next[s][c].end());
        entry = next;
        if (!prev_flat.empty() && prev_flat.size() == flat.size()) {
            double dmax = 0.0;
            for (std::size_t a = 0; a < flat.size(); ++a)
                dmax = std::max(dmax, std::fabs(num_traits<T>::to_double(flat[a]) -
                                                num_traits<T>::to_double(prev_flat[a])));
            prev_flat = flat;
            if (dmax < o.iter_tol) break;
        } else {
            prev_flat = flat;
        }
    }

    const T nan = num_traits<T>::from_double(std::numeric_limits<double>::quiet_NaN());
    out.hitprob.assign(nc, std::vector<T>(K, nan));
    out.missprob.assign(nc, std::vector<T>(K, nan));
    for (std::size_t c = 0; c < nc; ++c) {
        std::vector<T> hitT(K, num_traits<T>::from_int(0)), missT(K, num_traits<T>::from_int(0));
        for (std::size_t s = 0; s < E; ++s) {
            if (wmass[s][c].empty()) continue;
            double sw = 0.0;
            bool ok = true;
            for (std::size_t j = 0; j < wmass[s][c].size(); ++j) {
                const double v = num_traits<T>::to_double(wmass[s][c][j]);
                if (!std::isfinite(v)) ok = false;
                sw += v;
            }
            if (!ok || !(sw > 0.0)) continue;
            std::size_t idx = tran[s].size();
            for (std::size_t q = 0; q < tran[s].size(); ++q)
                if (tran[s][q].node == out.nodes[c]) idx = q;
            if (idx == tran[s].size()) continue;
            const fluid::FluidCacheqnTranCache<T>& tc = tran[s][idx];
            const T pe = num_traits<T>::from_double(e.prob_env[s]);
            for (std::size_t k = 0; k < K && k < tc.arate.size(); ++k) {
                if (!(num_traits<T>::to_double(tc.arate[k]) > 0.0)) continue;
                T hbar = num_traits<T>::from_int(0), mbar = num_traits<T>::from_int(0);
                for (std::size_t j = 0; j < wmass[s][c].size() && j < tc.hitprob_t.cols(); ++j) {
                    hbar = T(hbar + T(tc.hitprob_t(k, j) * wmass[s][c][j]));
                    mbar = T(mbar + T(tc.missprob_t(k, j) * wmass[s][c][j]));
                }
                const T swT = num_traits<T>::from_double(sw);
                hitT[k] = T(hitT[k] + T(T(pe * tc.arate[k]) * T(hbar / swT)));
                missT[k] = T(missT[k] + T(T(pe * tc.arate[k]) * T(mbar / swT)));
            }
        }
        for (std::size_t k = 0; k < K; ++k) {
            const T tot = T(hitT[k] + missT[k]);
            if (num_traits<T>::to_double(tot) > 0.0) {
                out.hitprob[c][k] = T(hitT[k] / tot);
                out.missprob[c][k] = T(missT[k] / tot);
            }
        }
    }
    return out;
}

/** A Cache in a stage is only servable by the fluid backend; name the other case. */
template <class T>
void refuse_cache_stages(const Environment<T>& e, const EnvOptions& o) {
    if (o.stage_solver == "fluid") return;
    for (std::size_t s = 0; s < e.nstages(); ++s) {
        // A layered stage carries no NetworkStruct at all, so there is no node
        // table to scan; a Cache inside an LQN is a CacheTask, which lives in
        // its host's LAYER and is SolverLN's to serve.
        if (e.is_lqn(s)) continue;
        const qn::NetworkStruct<T>& sn = e.stage(s).model;
        for (std::size_t ind = 0; ind < sn.nodes.size(); ++ind)
            if (sn.nodes[ind].nodetype == lang::NodeType::Cache)
                throw UnsupportedError(
                    "SolverENV meanfield: stage " + std::to_string(s + 1) +
                    " holds a Cache, whose environment-blended hit and miss ratios come from "
                    "aggregateCacheMeanfield_ and its per-sweep solver_fld_cacheqn_tran RMF "
                    "transient; that transient exists only for a FLUID stage solver, and this "
                    "ensemble runs '" +
                    o.stage_solver + "' stages");
    }
}

}  // namespace meanfield_detail

/** The mean-field solve on the original stages, with no compression. */
template <class T>
EnvMeanfieldSolution<T> solver_env_meanfield(Environment<T>& e, const EnvOptions& o) {
    meanfield_detail::refuse_cache_stages(e, o);
    EnvMeanfieldSolution<T> out;
    SolverEnv<T> s(e, o);
    out.avg = s.solve();
    out.cache = meanfield_detail::aggregate_cache_meanfield(e, o);
    return out;
}

/**
 * The compressed mean-field solve: aggregate the environment, then run the
 * mean-field fixed point over the macro-states.
 *
 * The compressed environment is kept alive by the returned structure, which is
 * what `SolverEnv` held a reference to; reading `.avg` out of the result and
 * discarding the rest is safe, but the compression it came from travels with it
 * so that eps and epsMAX can be checked against the numbers they produced.
 */
template <class T>
EnvMeanfieldSolution<T> solver_env_meanfield(Environment<T>& e, const EnvOptions& o,
                                             const EnvCompressOptions& c) {
    meanfield_detail::refuse_cache_stages(e, o);
    EnvMeanfieldSolution<T> out;
    out.compression = env_compress(e, c);
    out.compressed = true;
    SolverEnv<T> s(*out.compression.env, o);
    env_apply_macro_probabilities(*out.compression.env, out.compression);
    out.avg = s.solve();
    out.cache = meanfield_detail::aggregate_cache_meanfield(*out.compression.env, o);
    return out;
}

}  // namespace env
}  // namespace line

#endif  // LINE_SOLVERS_ENV_SOLVER_ENV_MEANFIELD_H
