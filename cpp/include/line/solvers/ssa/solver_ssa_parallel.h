/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_SSA_SOLVER_SSA_PARALLEL_H
#define LINE_SOLVERS_SSA_SOLVER_SSA_PARALLEL_H

/**
 * SolverSSA, the `para` / `parallel` method: a port of
 * `solver_ssa_analyzer_parallel.m`.
 *
 * REPLICATION IS NOT A LONGER RUN. This is the whole point of the method and
 * the one thing a reader must take from this file. The budget `opt.samples` is
 * SPLIT across R independent replicas of `ceil(samples/R)` firings each, not
 * multiplied by R; each replica is a complete run of the serial engine from the
 * same initial state on its own random stream, and each produces its own time
 * average. What is returned is the MEAN OF THE R REPLICA ESTIMATES, and the
 * quantity that says how good it is, is the spread BETWEEN those R numbers --
 * not the number of firings behind them. Concatenating the R sample paths and
 * time-averaging the whole would give a similar point estimate and an error bar
 * too small by a factor that grows with the autocorrelation of the path, which
 * is precisely the mistake replication exists to avoid. So `avg` is the replica
 * mean, `*_sem` is the standard error OF THAT MEAN computed across replicas
 * with R-1 degrees of freedom, and both are reported together.
 *
 * WORKER-COUNT INVARIANCE. Replica r (0-based) is seeded `base_seed + r` and
 * given `ceil(samples/R)` firings, so its result is a function of
 * `(base_seed, samples, R)` alone. The reference's header dwells on this
 * because its earlier `spmd` implementation both divided the budget by, and
 * seeded from, the RUNTIME number of workers, so the same script answered
 * differently on two machines. Nothing here can depend on a worker count for a
 * second reason, below.
 *
 * THERE IS NO THREADING, AND THAT IS A DELIBERATE CHOICE, NOT A GAP. The
 * replicas are run one after another in the loop below. Three reasons, in order
 * of weight: (i) the parallelism is over REPLICAS and each replica's answer is
 * pinned by its seed, so the returned numbers are bit-for-bit the same whether
 * the loop is serial or spread over cores -- the concurrency is an
 * implementation detail of how long the call takes and is not part of the
 * answer, which is exactly the property the reference had to work to recover;
 * (ii) nothing else in this port starts a thread and the build declares no
 * threading dependency, so introducing one here would be a build-system change
 * made for a wall-clock gain in a solver whose cost is already the user's
 * chosen sample budget; (iii) `SsaSerialEngine` holds `sn` by const reference
 * and mutates nothing shared, so the loop is trivially parallelizable later by
 * whoever wants to pay for the dependency -- the estimator does not change.
 *
 * WHAT THE CACHE WRITE-BACK AVERAGES. The reference averages each Cache node's
 * realized hit and miss probabilities over the replicas with a fixed 1/R
 * weight. That is kept, NaN and all: a replica in which a cache saw no reads
 * reports 0/0 there, and the reference propagates it into the average rather
 * than dropping the replica, which is the honest outcome -- the estimate really
 * is undefined when a replica contributes no reads.
 *
 * DOUBLE ONLY, for the reason `solver_ssa_serial.h` and `solver_ssa_nrm.h` both
 * give: the sample path comes out of exponential clocks drawn as `-log(u)/rate`
 * and the answer's error is Monte Carlo error, not rounding.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <string>
#include <type_traits>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/solvers/ssa/solver_ssa_serial.h"
#include "line/solvers/ssa/ssa_event_cache.h"
#include "line/solvers/ssa/ssa_types.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace ssa {

/**
 * The replicated engine's knobs: the serial engine's, plus the two the
 * reference reads from `options.config` on this path.
 *
 * `nreplicas` defaults to 8, which is `SolverOptions('SSA')`'s own default and
 * not a number chosen here; `eventcache` defaults to false, which is what the
 * reference sets for every `options.lang` except `java`.
 */
struct SsaParallelOptions : SsaSerialOptions {
    /** `options.config.nreplicas`: R, the FIXED number of independent replicas. */
    std::size_t nreplicas = 8;
    /** `options.config.eventcache`: memoize `after_event` inside each replica. */
    bool eventcache = false;
};

/**
 * What the replicated analyzer returns.
 *
 * The per-replica tables are kept, not just their mean, because the mean alone
 * cannot be audited: a caller who wants a different confidence level, a
 * different combination rule or a look at whether one replica is an outlier
 * needs the R numbers. What is NOT kept is the R sample paths -- each is
 * `samples_per_replica` long and the estimator never looks at them again.
 */
template <class T>
struct SsaParallelSolution {
    /** The replica MEAN; `method` = "parallel". */
    SsaSolution avg;
    /** The R per-replica estimates, in replica order. */
    std::vector<SsaSolution> replica;
    /** The stream each replica ran on: `base_seed + r`. */
    std::vector<unsigned long> seed;

    std::size_t nreplicas = 0;
    /** `ceil(samples / R)`: the firings EACH replica performed. */
    std::size_t samples_per_replica = 0;
    /** `opt.samples`: the budget asked for, which `R * samples_per_replica` rounds up. */
    std::size_t samples_requested = 0;
    unsigned long base_seed = 0;

    /**
     * Standard error of the replica mean, per metric: `s / sqrt(R)` with `s`
     * the sample standard deviation ACROSS the R replica estimates.
     *
     * NOT in the reference, which returns the mean alone. It is here because a
     * replicated estimate whose between-replica spread is discarded is
     * indistinguishable from a single long run, and the whole reason to
     * replicate is that the two have different error bars. NaN when R = 1,
     * where no spread is observable -- which is the correct report, not zero.
     */
    Matrix<double> QN_sem, UN_sem, RN_sem, TN_sem;
    std::vector<double> XN_sem, CN_sem;

    /** The Cache write-back of the reference's final loop, averaged over replicas. */
    std::vector<SsaCacheRatio> cache;
};

namespace parallel_detail {

/** Sample mean and standard error of the mean of `x`, both NaN when it is short. */
inline void mean_sem(const std::vector<double>& x, double& mean, double& sem) {
    const std::size_t R = x.size();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (R == 0) {
        mean = nan;
        sem = nan;
        return;
    }
    double s = 0.0;
    for (std::size_t r = 0; r < R; ++r) s += x[r];
    mean = s / static_cast<double>(R);
    if (R < 2) {
        // One replica has no observable spread. Reporting 0 would claim the
        // estimate is exact, which is the failure mode this field exists to
        // prevent, so it reports "unknown" instead.
        sem = nan;
        return;
    }
    double ss = 0.0;
    for (std::size_t r = 0; r < R; ++r) ss += (x[r] - mean) * (x[r] - mean);
    sem = std::sqrt(ss / static_cast<double>(R - 1)) / std::sqrt(static_cast<double>(R));
}

/** Elementwise `mean_sem` over the same cell of every replica's matrix. */
inline void reduce_matrix(const std::vector<SsaSolution>& rep, Matrix<double> SsaSolution::*m,
                          Matrix<double>& mean, Matrix<double>& sem) {
    const std::size_t R = rep.size();
    const std::size_t nr = R ? (rep[0].*m).rows() : 0;
    const std::size_t nc = R ? (rep[0].*m).cols() : 0;
    mean = Matrix<double>(nr, nc, 0.0);
    sem = Matrix<double>(nr, nc, 0.0);
    std::vector<double> col(R, 0.0);
    for (std::size_t i = 0; i < nr; ++i)
        for (std::size_t j = 0; j < nc; ++j) {
            for (std::size_t r = 0; r < R; ++r) col[r] = (rep[r].*m)(i, j);
            double mu = 0.0, se = 0.0;
            mean_sem(col, mu, se);
            mean(i, j) = mu;
            sem(i, j) = se;
        }
}

/** Elementwise `mean_sem` over the same entry of every replica's vector. */
inline void reduce_vector(const std::vector<SsaSolution>& rep, std::vector<double> SsaSolution::*v,
                          std::vector<double>& mean, std::vector<double>& sem) {
    const std::size_t R = rep.size();
    const std::size_t n = R ? (rep[0].*v).size() : 0;
    mean.assign(n, 0.0);
    sem.assign(n, 0.0);
    std::vector<double> col(R, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t r = 0; r < R; ++r) col[r] = (rep[r].*v)[i];
        mean_sem(col, mean[i], sem[i]);
    }
}

}  // namespace parallel_detail

/**
 * `solver_ssa_analyzer_parallel.m`: run R replicas of the serial engine and
 * combine their estimates.
 *
 * The reference's `run_replica` subfunction is `solver_ssa_serial_analyzer`
 * here: the two compute the same per-station table from one sample path, and
 * factoring the replication away from the estimator is what makes it obvious
 * that the combination below touches only the R finished numbers.
 */
template <class T>
SsaParallelSolution<T> solver_ssa_parallel_analyzer(const qn::NetworkStruct<T>& sn,
                                                    const SsaParallelOptions& opt) {
    // `if constexpr`, not a run-time test: the replica reaches `map_mean` and
    // the logarithm of a uniform, so a Rational instantiation would fail to
    // COMPILE rather than refuse. The gate keeps the body uninstantiated.
    if constexpr (!std::is_same<T, double>::value) {
        (void)sn;
        (void)opt;
        throw UnsupportedError(
            "solver_ssa_parallel: an SSA sample path is generated from exponential clocks, which "
            "are logarithms of uniform draws; there is no exact value to compute and a wider "
            "float carries no information the Monte Carlo error does not swamp. Rerun with "
            "--arith double");
    } else {
        // The replica holding time is drawn as -mean*log(u), so the backend must
        // have a logarithm at all. The gate above keeps a backend without one
        // from reaching here; the assert names the reason rather than letting
        // the failure surface inside the uniform draw.
        static_assert(num_traits<T>::has_transcendental,
                      "solver_ssa_parallel: each replica is an SSA sample path generated from "
                      "exponential clocks drawn as -mean*log(u), which needs transcendental "
                      "arithmetic");

        if (opt.eventcache)
            throw UnsupportedError(
                "SolverSSA(method='parallel'): options.config.eventcache is set. The memo itself "
                "is ported (ssa::SsaEventCache, which the reference builds per replica via "
                "EventCache.create) and is verified against recomputation, but SsaSerialEngine "
                "still calls the free qn::after_event and takes no cache argument, so setting the "
                "flag here would advertise a memo the replicas never consult. Wire the engine's "
                "enabled() through SsaEventCache::after_event first");

        const std::size_t R = opt.nreplicas > 0 ? opt.nreplicas : 1;
        // ceil(samples/R): the budget is SPLIT, so the total firings performed
        // is R*ceil(samples/R), which rounds the request UP and never down --
        // a replica of zero firings would have no estimate to contribute.
        std::size_t per = (opt.samples + R - 1) / R;
        if (per == 0) per = 1;  // a replica of zero firings has no estimate to contribute

        SsaParallelSolution<T> out;
        out.nreplicas = R;
        out.samples_per_replica = per;
        out.samples_requested = opt.samples;
        out.base_seed = opt.seed;
        out.replica.reserve(R);
        out.seed.reserve(R);

        // Only the finished ESTIMATE of each replica is kept, never its sample
        // path: the combination below looks at the R tables and at nothing
        // else, and holding R paths at once would make the peak memory grow
        // with a budget the estimator does not read.
        std::vector<std::vector<SsaCacheRatio> > rep_cache;
        rep_cache.reserve(R);
        for (std::size_t r = 0; r < R; ++r) {
            // Sliced to the base on purpose: `nreplicas` and `eventcache` are
            // the replicATION's knobs and mean nothing inside a replica, while
            // everything the serial engine reads (cutoff, warmupfrac,
            // state_max) is carried through unchanged, as the reference's
            // `repoptions = laboptions` carries it.
            SsaSerialOptions ropt = opt;
            ropt.method = "serial";
            ropt.samples = per;
            ropt.verbose = false;  // the reference's VerboseLevel.SILENT per replica
            // Replica r's stream is fixed by r alone, which is what makes the
            // combined answer independent of the order the replicas run in and
            // of how many run at once.
            ropt.seed = opt.seed + static_cast<unsigned long>(r);
            const SsaSerialSolution<T> one = solver_ssa_serial_analyzer(sn, ropt);
            out.replica.push_back(one.avg);
            rep_cache.push_back(one.cache);
            out.seed.push_back(ropt.seed);
        }

        SsaSolution& a = out.avg;
        parallel_detail::reduce_matrix(out.replica, &SsaSolution::QN, a.QN, out.QN_sem);
        parallel_detail::reduce_matrix(out.replica, &SsaSolution::UN, a.UN, out.UN_sem);
        parallel_detail::reduce_matrix(out.replica, &SsaSolution::RN, a.RN, out.RN_sem);
        parallel_detail::reduce_matrix(out.replica, &SsaSolution::TN, a.TN, out.TN_sem);
        parallel_detail::reduce_vector(out.replica, &SsaSolution::XN, a.XN, out.XN_sem);
        parallel_detail::reduce_vector(out.replica, &SsaSolution::CN, a.CN, out.CN_sem);

        // The reference averages RN and CN with the same 1/R weight it uses for
        // the rest rather than recomputing them from the averaged QN and XN.
        // The two differ, because a ratio of means is not a mean of ratios, and
        // it is the reference's estimator that is reported here.

        a.method = "parallel";
        // The FIRINGS SPENT, which is not the effective sample size of the
        // estimate: that is R, the number of independent numbers averaged.
        // Anyone reading this field as a precision must read `*_sem` instead.
        a.samples = 0;
        a.simulated_time = 0.0;
        for (std::size_t r = 0; r < R; ++r) {
            a.samples += out.replica[r].samples;
            a.simulated_time += out.replica[r].simulated_time;
        }

        // The Cache write-back. `hitclass` is a property of the struct, not of a
        // replica, so the reference's `length(...hitclass) >= k` test is
        // replica-independent and is taken once from `sn`: a class outside it
        // keeps the zero the reference initializes and never adds to.
        const std::size_t K = sn.nclasses;
        for (typename std::map<std::size_t, qn::CacheParam<T> >::const_iterator ci =
                 sn.nodeparam.begin();
             ci != sn.nodeparam.end(); ++ci) {
            const std::size_t ind = ci->first;
            if (ind == 0 || ind > sn.nodes.size()) continue;
            if (sn.nodes[ind - 1].nodetype != lang::NodeType::Cache) continue;
            SsaCacheRatio cr;
            cr.node = ind;
            cr.hitprob.assign(K, 0.0);
            cr.missprob.assign(K, 0.0);
            cr.residt.assign(K, std::numeric_limits<double>::quiet_NaN());
            std::vector<double> dly(K, 0.0);
            bool any_delayed = false;
            for (std::size_t k = 1; k <= K; ++k) {
                if (ci->second.hitclass.size() < k || ci->second.missclass.size() < k) continue;
                double h = 0.0, m = 0.0, d = 0.0;
                for (std::size_t r = 0; r < R; ++r)
                    for (std::size_t c = 0; c < rep_cache[r].size(); ++c) {
                        if (rep_cache[r][c].node != ind) continue;
                        h += rep_cache[r][c].hitprob[k - 1] / static_cast<double>(R);
                        m += rep_cache[r][c].missprob[k - 1] / static_cast<double>(R);
                        // A replica that saw no merge leaves the field empty
                        // rather than reporting a zero share it never measured.
                        if (!rep_cache[r][c].delayedprob.empty()) {
                            d += rep_cache[r][c].delayedprob[k - 1] / static_cast<double>(R);
                            any_delayed = true;
                        }
                    }
                cr.hitprob[k - 1] = h;
                cr.missprob[k - 1] = m;
                dly[k - 1] = d;
            }
            if (any_delayed) cr.delayedprob = dly;
            out.cache.push_back(cr);
        }
        return out;
    }
}

/**
 * The `para` / `parallel` entry of `solver_ssa_analyzer.m`.
 *
 * The reference reaches this only after its NRM eligibility gate has declined
 * the model (`solver_ssa_analyzer.m` lines 143-157 run the NRM instead when it
 * is eligible, because one fast exact-enough run beats replicated simulation).
 * That preference belongs to the dispatcher and is not duplicated here, so this
 * entry always replicates the serial engine; asking it for `serial` says so by
 * name rather than quietly answering with one replica, whose error bar is a
 * factor sqrt(R) wider than the one requested.
 */
template <class T>
SsaParallelSolution<T> solver_ssa_parallel(const qn::NetworkStruct<T>& sn,
                                           const SsaParallelOptions& opt) {
    const std::string& m = opt.method;
    if (m == "para" || m == "parallel" || m == "default")
        return solver_ssa_parallel_analyzer(sn, opt);
    if (m == "serial" || m == "ssa")
        throw UnsupportedError(
            "SolverSSA(parallel): the '" + m +
            "' method is ONE run of the serial engine, not the mean of " +
            std::to_string(opt.nreplicas) +
            " independent replicas, and the two report the same quantity at different variances. "
            "Call solver_ssa_serial for it");
    throw UnsupportedError("SolverSSA(parallel): '" + m +
                           "' is not a method this entry accepts; it implements 'para' and "
                           "'parallel' and the 'default' alias that reaches them");
}

}  // namespace ssa
}  // namespace line

#endif  // LINE_SOLVERS_SSA_SOLVER_SSA_PARALLEL_H
