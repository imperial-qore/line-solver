/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_NC_NC_TYPES_H
#define LINE_SOLVERS_NC_NC_TYPES_H

/**
 * Controls and result shape shared by the normalizing-constant analyzers.
 *
 * The result is `mva::MvaSolution<T>`, not a type of its own: SolverNC returns
 * the same [Q,U,R,T,C,X,lG] tuple as SolverMVA, and sharing the struct is what
 * lets the fork-join fixed point (`fj_driver.h`) drive `nc_dispatch` as its
 * inner solve without a shim.
 */

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/solvers/mva/mva_types.h"
#include "line/solvers/cache_metrics.h"
#include "line/util/matrix.h"

namespace line {
namespace nc {

/** Controls, defaulting to `SolverOptions('NC')` in the reference. */
struct NcSolverOptions {
    std::string method = "default";
    double tol = 1e-4;       ///< options.tol
    double iter_tol = 1e-4;  ///< options.iter_tol, the eta stopping test
    int iter_max = 1000;
    /**
     * `options.config.highvar`. SolverOptions('NC') overrides the global
     * default and sets this to 'interp', so a non-product-form FCFS station is
     * rescaled by the WSC 2020 interpolation and the analyzer iterates. Setting
     * it to 'default' makes the outer loop run exactly once.
     */
    std::string highvar = "interp";
    /**
     * `options.config.fork_join`: which fork-join arm the shared fixed point
     * takes on a model with a Fork. `default`/`mmt`/`fjt` is the MMT transform,
     * `ht`/`heidelberger-trivedi` the Heidelberger-Trivedi one. See
     * `MvaOptions::fork_join`; the driver is the same function for both solvers.
     */
    std::string fork_join = "default";
    /**
     * `options.config.multiserver`: how a finite multiserver station is
     * represented. `default` keeps the historical dispatch (Seidmann on
     * `default`, the exact mu(n)=min(n,c) lattice on `exact`/`is`/`panald`,
     * and the lattice on the 2-station Delay+multiserver shape); `seidmann`
     * forces Seidmann everywhere, `exact`/`lld` forces the lattice everywhere
     * it is admissible. The SolverMVA values this solver has no counterpart for
     * (`softmin`, `conway`, `krzesinski`, `suri`, `erlang`) fall back to
     * `default` rather than throwing: one options object is commonly reused
     * across solvers. See `MvaOptions::multiserver`.
     */
    std::string multiserver = "default";
    std::size_t samples = 100000;  ///< options.samples, read by the estimators
    unsigned long seed = 23000;    ///< options.seed, read by the estimators
    /**
     * `options.config.mcmc_batches`: the non-overlapping batches the
     * Chen-O'Cinneide estimator splits its run into for the batch-means
     * confidence intervals. 30 is Schmeiser (1982), the count used in the
     * tables of the paper.
     */
    /**
     * `options.config.aghq_nodes`: nodes per simplex direction of the adaptive
     * Gauss-Hermite rule. q = 1 reproduces pfqn_le; the rule costs q^(M-1)
     * evaluations, so the default stays small.
     */
    std::size_t aghq_nodes = 3;
    std::size_t mcmc_batches = 30;
    /**
     * `options.config.mcmc_burnin`: the warm-up fraction that estimator
     * discards before it starts accumulating. The paper discards none and
     * ignores the initialization bias.
     */
    double mcmc_burnin = 0.1;
    /**
     * `options.config.algorithm` for `getCdfRespT`: 'exact' selects
     * `pfqn_stdf`, 'rd' the heuristic `pfqn_stdf_heur`. The reference defaults
     * this field into the options struct when it is absent.
     */
    std::string cdf_algorithm = "exact";
    /** `options.config.mem_tol` / `mem_maxiter`, read by the MEM fixed point. */
    double mem_tol = 1e-6;
    long mem_maxiter = 1000;
    /**
     * Whether the model this solve came from has a Fork, which the model handed
     * to the analyzer no longer does. Same role as `MvaOptions::base_has_fork`.
     */
    bool base_has_fork = false;
    /**
     * `options.config.slotted`. Runs the solve on a discrete time scale, the
     * same switch SolverLDES uses. It routes to `solver_nc_dt` and REFUSES a
     * model outside the discrete-time product form rather than falling back to
     * the continuous-time analyzers: the continuous-time answer to a slotted
     * question is a different number, not a worse one.
     */
    bool slotted = false;
    /** `options.config.slotlength`, the slot length in model time units. */
    double slotlength = 1.0;
};

/** The `[Q,U,R,T,C,X,lG]` of the reference, plus the algorithm that ran. */
template <class T>
struct NcSolution {
    mva::MvaSolution<T> sol;
    std::string actualmethod;
    Matrix<T> STeff;  ///< the service times of the last pass, MATLAB's STeff
    /**
     * Set only by the integrated cacheqn branch: the converged struct whose
     * routing carries the ACTUAL hit/miss probabilities, from which the runner
     * derives ArvR and ResidT. Empty for every other model. Same role as
     * `mva::DispatchResult::refreshed_struct`.
     */
    std::shared_ptr<qn::NetworkStruct<T>> refreshed_struct;
    /**
     * The reference's warning text, verbatim, empty when it did not warn. Same
     * role as `mva::AvgResult::warning`, which `solver_nc_run_analyzer` copies it into:
     * a statement about HOW the numbers were produced. The cache branch uses it
     * for the two cost-cap conditions -- a promotion path the caps block, and a
     * method that has no cost-capped counterpart and was switched -- because the
     * port has no `line_warning` channel and a FLAG no caller reads is the same
     * as saying nothing.
     */
    std::string warning;
    /**
     * (h) mean storage cost held by each cache list, K_j = sum_i sigma_i pi_ij.
     * EMPTY for every model without item sizes, which is every model but a cache
     * carrying `sn.nodeparam.itemsize` (ton21cache Sec. IX).
     */
    std::vector<T> listcost;
    /**
     * What the cache branches observed, EMPTY on a model with no Cache node.
     *
     * It rides here rather than being recomputed by the caller because the
     * analyzers already have it and throwing it away would mean solving the
     * cache twice to answer `getAvgCacheTable` -- and the second solve, being a
     * different call, would be free to disagree with the first.
     */
    solvers::CacheMetrics<T> cache;
};

}  // namespace nc
}  // namespace line

#endif  // LINE_SOLVERS_NC_NC_TYPES_H
