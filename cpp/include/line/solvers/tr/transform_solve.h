// Copyright (c) 2012-2026, QORE Lab, Imperial College London
// All rights reserved.
#ifndef LINE_SOLVERS_TR_TRANSFORM_SOLVE_H
#define LINE_SOLVERS_TR_TRANSFORM_SOLVE_H

/**
 * @file
 * Solver-agnostic driver of a model TRANSFORMATION, the sibling of
 * `solvers/mva/fj_driver.h`.
 *
 * A transformation rewrites the model into one or more subproblems, solves
 * those with a REAL solver, and maps the metrics back onto the original classes
 * and stations:
 *
 *     expand -> solve -> lift
 *
 * THE INNER SOLVER IS THE OUTER SOLVER. The caller supplies `inner_solve` as a
 * template parameter, exactly as `fj_driver.h` takes its `InnerSolve`, so a
 * transformation written once serves every engine rather than the one it was
 * first written for. `solver_ctmc_chain_aggregation` used to call
 * `solver_ctmc_run_analyzer` directly, which is what kept a transform with
 * nothing CTMC-specific in it out of reach of MVA, NC and Fluid.
 *
 * WHAT "ANY SOLVER" MEANS HERE, AND WHAT IT DOES NOT. The inner solve is a
 * TEMPLATED ANALYZER, not the `solvers/solver.h` facade, which is a
 * double-only interface whose only accessor is `avg_table()`. That reaches
 * MVA, NC, CTMC and Fluid; it does NOT reach JMT or LDES, which are subprocess
 * wrappers behind that facade.
 *
 * NO METHOD NAME TABLE HERE, deliberately. MATLAB, the JAR and python resolve
 * `options.config.transform` through a method name table because their options carry
 * a string-keyed config map. C++ selects the transform at the call site from a
 * typed option (`CtmcOptions::chain_aggregation`), so a runtime method name table
 * would add a failure mode without adding reach.
 *
 * Mirrors MATLAB `@NetworkSolver/transformSolve.m` with
 * `solver_tr_chains_analyzer.m`, the JAR `jline.solvers.tr.TransformSolve` with
 * `ChainsStrategy`, and python `line_solver/solvers/transform_driver.py`.
 */

#include <string>

#include <algorithm>
#include <cmath>
#include <vector>

#include "line/api/pfqn/pfqn_bk.h"
#include "line/api/sn/sn_aggregate_chains.h"
#include "line/api/sn/sn_remove_class.h"
#include "line/lang/dist_scale_rate.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mva/sn_chain.h"
#include "line/solvers/mva/solver_mva_runner.h"

namespace line {
namespace tr {

/**
 * Chain aggregation: collapse every chain onto a single class, solve, and map
 * the chain metrics back onto the classes.
 *
 * `api::sn_aggregate_chains` builds the collapsed model and
 * `mva::sn_deaggregate_chain_results` maps its metrics back through alpha, the
 * per-station share of the chain's visits each class carries.
 *
 * WHAT IS TRADED. Exactness on a non-product-form model: one aggregate service
 * law, fitted to the alpha-weighted first two moments, replaces the per-class
 * ones. On a product-form model the chain IS the unit MVA and convolution
 * already solve in, so the answer is exact and the state space is the smaller
 * one.
 *
 * SINGLE PASS: one solve of the aggregate determines the answer, so there is no
 * sweep and no convergence test.
 *
 * @param sn          the original struct
 * @param inner_solve solves the aggregated struct, normally with the caller's
 *                    own analyzer bound to its own options
 * @param method      the caller's method name, for the reported compound name
 */
template <class T, class InnerSolve>
mva::AvgResult<T> transform_solve_chains(const qn::NetworkStruct<T>& sn,
                                         InnerSolve inner_solve,
                                         const std::string& method) {
    // expand
    api::ChainAggregationResult<T> agg = api::sn_aggregate_chains(sn);

    // solve
    const mva::AvgResult<T> chain = inner_solve(agg.model.get_struct());

    // lift. The deaggregation reads the ORIGINAL struct's chain demands, which
    // is where alpha, Lchain, STchain and Vchain all come from; the transform's
    // own copy of them is the same object and is left to the caller.
    const mva::ChainDemands<T> d = mva::sn_get_demands_chain(sn);
    const mva::ClassResults<T> cls = mva::sn_deaggregate_chain_results(
        sn, d, chain.QN, chain.UN, chain.RN, chain.TN, chain.XN);

    mva::AvgResult<T> out;
    out.QN = cls.Q;
    out.UN = cls.U;
    out.RN = cls.R;
    out.TN = cls.Tp;
    out.CN = cls.C;
    out.XN = cls.X;
    out.AN = mva::sn_get_arvr_from_tput(sn, out.TN);
    out.WN = mva::sn_get_residt_from_respt(sn, out.RN);
    out.method = method;
    out.actualmethod = method + "/chainaggr";
    return out;
}


/**
 * LOAD CONCEALMENT (Birman-Kogan Algorithm 2) as a transformation, and the
 * first ITERATED one.
 *
 * Chain `l` is solved on its own against the residual capacity the others leave
 * it, `A_i = 1 - sum_{k!=l} L(i,k) X_k`, so it sees the concealed demand
 * `L(i,l)/A_i`. Sweeping the chains in GAUSS-SEIDEL order -- publishing each
 * chain's throughput the moment it is known, so chain `l+1` of the same sweep
 * already sees it -- and iterating to a fixed point is the algorithm.
 *
 * WHAT THIS ADDS OVER THE KERNEL. `pfqn_bklc` solves each single-chain
 * subproblem on a DEMAND VECTOR with the inner solve hard-wired to MVA or the
 * uniform expansion. Here the subproblem is a real single-class struct, so the
 * inner solve is whichever analyzer the caller bound, which is what makes the
 * concealment approximation measurable rather than merely asserted.
 *
 * IT IS NOT A STRICTLY BETTER LC: the kernel sees only `L`, while the chain
 * aggregation refits the chain service law to two moments.
 *
 * THE TOLERANCE IS FIXED AT 1e-10, a PARITY requirement and not a knob: a
 * looser one stops the sweep at a different iteration in each codebase.
 *
 * Mirrors MATLAB `solver_tr_lc_analyzer.m`, python `_lc_strategy` and the JAR
 * `jline.solvers.tr.LcStrategy`.
 *
 * @param sn          the original struct
 * @param inner_solve solves one concealed single-class struct
 * @param method      the caller's method name, for the reported compound name
 * @param iter_max    sweep cap; the kernel's own default is 1000
 */
template <class T, class InnerSolve>
mva::AvgResult<T> transform_solve_lc(const qn::NetworkStruct<T>& sn,
                                     InnerSolve inner_solve,
                                     const std::string& method,
                                     std::size_t iter_max = 1000) {
    const double kTol = 1e-10;
    const double kFineTol = 1e-12;

    // expand: chain aggregation, then one single-class struct per chain
    api::ChainAggregationResult<T> agg = api::sn_aggregate_chains(sn);
    const qn::NetworkStruct<T> snChain = agg.model.get_struct();
    const mva::ChainDemands<T> dem = mva::sn_get_demands_chain(snChain);
    const std::size_t M = snChain.nstations;
    const std::size_t R = dem.Lchain.cols();
    if (snChain.nclasses != R)
        throw InputError("transform_solve_lc: the chain-aggregated model has a class count that "
                         "differs from its chain count; load concealment needs one class per chain");

    // The concealment slows QUEUEING stations only: a delay holds no queue, so
    // zeroing its rows makes A come out as exactly 1 there.
    std::vector<bool> is_delay(M, false);
    Matrix<T> L(M, R, num_traits<T>::from_int(0));
    std::vector<double> Z(R, 0.0), N(R, 0.0);
    for (std::size_t i = 0; i < M; ++i)
        is_delay[i] = snChain.stations[i].sched == qn::SchedStrategy::INF;
    for (std::size_t r = 0; r < R; ++r) {
        N[r] = dem.Nchain[r];
        for (std::size_t i = 0; i < M; ++i) {
            if (is_delay[i])
                Z[r] += num_traits<T>::to_double(dem.Lchain(i, r));
            else
                L(i, r) = dem.Lchain(i, r);
        }
    }

    // One single-class struct per chain, built ONCE and re-concealed in place.
    std::vector<qn::NetworkStruct<T>> subs;
    std::vector<std::vector<lang::Distrib<T>>> base(R);
    subs.reserve(R);
    for (std::size_t l = 0; l < R; ++l) {
        qn::NetworkStruct<T> m = snChain;
        // Descending, so removing a higher class never shifts the index of the
        // one being kept or of any class still to be removed.
        for (std::size_t k = R; k >= 1; --k)
            if (k != l + 1) m = api::sn_remove_class(m, k);
        subs.push_back(m);
        base[l].resize(M);
        for (std::size_t i = 0; i < M; ++i)
            if (!is_delay[i]) base[l][i] = m.service[i][0];
    }

    // Step 1 of Algorithm 2, reproducing the kernel's own seed: the saddle point
    // utilizations of Corollary 1, the N/(Z+sum L) fallback and the 1/max L
    // capacity clamp. Seeding identically is what keeps the SWEEP COUNT, and not
    // only the fixed point, comparable with the kernel and the other codebases.
    std::vector<double> X(R, 0.0);
    // GUARDED AT COMPILE TIME. pfqn_bk static_asserts on
    // num_traits<T>::has_transcendental (its saddle point expansion needs log),
    // so at T = Rational the call must not be instantiated at all; an `if` alone
    // would still compile the branch and fire the assert. Rational then starts
    // from the closed-form bound below, which costs a few extra sweeps but
    // reaches the same fixed point.
    if constexpr (num_traits<T>::has_transcendental) {
        std::vector<T> Nv(R), Zv(R);
        for (std::size_t r = 0; r < R; ++r) {
            Nv[r] = num_traits<T>::from_double(N[r]);
            Zv[r] = num_traits<T>::from_double(Z[r]);
        }
        try {
            const pfqn::BkResult<T> bk = pfqn::pfqn_bk(L, Nv, Zv);
            if (bk.X.size() == R)
                for (std::size_t r = 0; r < R; ++r) X[r] = num_traits<T>::to_double(bk.X[r]);
        } catch (const std::exception&) {
            // the seed is a starting point, not an answer: an unusable saddle
            // point falls through to the closed-form bound below
        }
    }
    for (std::size_t r = 0; r < R; ++r) {
        if (!std::isfinite(X[r]) || X[r] < 0) X[r] = 0.0;
        double sum = Z[r], cap = 0.0;
        for (std::size_t i = 0; i < M; ++i) {
            const double l = num_traits<T>::to_double(L(i, r));
            sum += l;
            cap = std::max(cap, l);
        }
        if (X[r] == 0.0 && N[r] > 0 && sum > 0) X[r] = N[r] / sum;
        if (cap > 0) X[r] = std::min(X[r], 1.0 / cap);
    }

    auto conceal_all = [&]() {
        for (std::size_t l = 0; l < R; ++l)
            for (std::size_t i = 0; i < M; ++i) {
                if (is_delay[i]) continue;
                double busy = 0.0;
                for (std::size_t k = 0; k < R; ++k)
                    if (k != l) busy += num_traits<T>::to_double(L(i, k)) * X[k];
                const double a = std::max(1.0 - busy, kFineTol);
                // dist_scale_rate multiplies the RATE, so a factor of A_i divides
                // the mean by A_i: exactly the concealed demand L(i,l)/A_i.
                subs[l].set_service(i + 1, 1, lang::dist_scale_rate(base[l][i],
                                                                    num_traits<T>::from_double(a)));
            }
    };
    conceal_all();

    std::vector<mva::AvgResult<T>> res(R);
    std::vector<double> Xold(R, -1.0);
    std::size_t iters = 0;
    for (std::size_t it = 1; it <= std::max<std::size_t>(1, iter_max); ++it) {
        iters = it;
        for (std::size_t l = 0; l < R; ++l) {
            res[l] = inner_solve(subs[l]);
            // GAUSS-SEIDEL: publish chain l now, so chain l+1 of this same sweep
            // already sees it.
            X[l] = res[l].XN.empty() ? 0.0 : num_traits<T>::to_double(res[l].XN[0]);
            if (!(X[l] >= 0)) X[l] = 0.0;
            conceal_all();
        }
        double diff = 0.0, scale = 1.0;
        for (std::size_t r = 0; r < R; ++r) {
            diff = std::max(diff, std::abs(X[r] - Xold[r]));
            scale = std::max(scale, std::abs(X[r]));
        }
        Xold = X;
        if (diff <= kTol * scale) break;
    }

    // lift: reassemble the chain table, then deaggregate onto the classes
    Matrix<T> Q(M, R, num_traits<T>::from_int(0)), U(M, R, num_traits<T>::from_int(0));
    Matrix<T> Rr(M, R, num_traits<T>::from_int(0)), Tp(M, R, num_traits<T>::from_int(0));
    std::vector<T> Xc(R, num_traits<T>::from_int(0));
    for (std::size_t l = 0; l < R; ++l) {
        for (std::size_t i = 0; i < M && i < res[l].QN.rows(); ++i) {
            Q(i, l) = res[l].QN(i, 0);
            U(i, l) = res[l].UN(i, 0);
            Rr(i, l) = res[l].RN(i, 0);
            Tp(i, l) = res[l].TN(i, 0);
        }
        Xc[l] = num_traits<T>::from_double(X[l]);
    }

    mva::AvgResult<T> out;
    if (sn.nchains >= sn.nclasses) {
        // ONE CLASS PER CHAIN: the chain answer already IS the class answer, so
        // the lift is a re-indexing rather than a deaggregation.
        out.QN = Q; out.UN = U; out.RN = Rr; out.TN = Tp; out.XN = Xc;
        out.CN.assign(R, num_traits<T>::from_int(0));
        for (std::size_t l = 0; l < R; ++l) {
            T acc = num_traits<T>::from_int(0);
            for (std::size_t i = 0; i < M; ++i) acc = acc + Rr(i, l);
            out.CN[l] = acc;
        }
    } else {
        const mva::ChainDemands<T> d = mva::sn_get_demands_chain(sn);
        const mva::ClassResults<T> cls =
            mva::sn_deaggregate_chain_results(sn, d, Q, U, Rr, Tp, Xc);
        out.QN = cls.Q; out.UN = cls.U; out.RN = cls.R;
        out.TN = cls.Tp; out.CN = cls.C; out.XN = cls.X;
    }
    out.AN = mva::sn_get_arvr_from_tput(sn, out.TN);
    out.WN = mva::sn_get_residt_from_respt(sn, out.RN);
    out.method = method;
    out.actualmethod = method + "/lc";
    out.iter = static_cast<int>(iters);
    return out;
}

}  // namespace tr
}  // namespace line

#endif  // LINE_SOLVERS_TR_TRANSFORM_SOLVE_H
