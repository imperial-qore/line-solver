/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_LN_LQN_ANALYZERS_H
#define LINE_SOLVERS_LN_LQN_ANALYZERS_H

/**
 * The @@SolverLN methods that solver_ln.h does not carry.
 *
 * solver_ln.h ports the fixed point itself (construct, buildLayers,
 * buildLayersRecursive, init, initInterlock, converged, analyze, post,
 * updateMetricsDefault, updatePopulations, updateThinkTimes, updateLayers,
 * updateRoutingProbabilities, getEntryServiceMatrix, getEnsembleAvg), and
 * lqn_helpers.h carries lqn_fwd_rendezvous and the lqn_overtake_markov chain.
 * What is here is the remainder of matlab/src/solvers/LN/@@SolverLN:
 *
 *   overtake_prob        - the reduced three-state CTMC
 *   overtake_prob_markov - the input mapping onto the LQNS overtaking chain
 *   convergedStoch       - the Robbins-Monro / Polyak-Ruppert controller
 *
 * plus thin wrappers, below, naming the three @@SolverLN entry points whose work
 * lives on the solver itself: getSensitivityTable, getTranAvg and getCdfRespT.
 * updateMetricsMomentBased (the `moment3` method) is likewise a member of
 * SolverLN, because it replaces updateMetricsDefault inside the iteration and
 * writes the same state.
 *
 * WHAT SolverLN CALLS AND WHAT IT DOES NOT. `lqn_overtake_prob_markov` IS wired
 * in: solver_ln.h computes the `servt_ph1` / `servt_ph2` split and calls it from
 * update_metrics, where updateMetricsDefault.m:313-351 calls
 * overtake_prob_markov. It is still written as a free function over quantities
 * the caller supplies, and the tested phase residence `xj` is still a PARAMETER
 * rather than something derived here, because the split belongs to the solver's
 * iterate and inventing it from a phase-1 model would be a guess dressed as a
 * result. `overtake_prob`, the reduced three-state CTMC, is the alternative the
 * reference keeps beside it and no path selects.
 *
 * `LnStochController` IS wired in, since `layer_solver = 'ssa'` gave the port a
 * stochastic layer engine: `SolverLN::iterate` builds one whenever that engine
 * is selected, drives `relax_omega` from it and reports the Polyak-Ruppert
 * average as the final iterate. It is NOT used for any deterministic engine --
 * the two convergence tests both rewrite `results`, so only one may be in
 * force.
 *
 * INCLUDE ORDER. This header needs `LayerResult` and `SolverLN` complete, so it
 * is parsed after solver_ln.h, which includes it at its foot; the declaration
 * solver_ln.h needs stands near the top of that file. Do not turn the include
 * below into a forward declaration -- LnStochController copies LayerResults.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/lang/lqn/lqn_struct.h"
#include "line/num/number.h"
#include "line/solvers/ln/lqn_helpers.h"
#include "line/solvers/ln/solver_ln.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace ln {

using lang::CallType;
using lang::Distrib;
using lang::GlobalConstants;
using lqn::LqnStruct;

// ---------------------------------------------------------------------------
// overtake_prob
// ---------------------------------------------------------------------------

/** Stationary law of the three-state overtaking chain, in its own order. */
template <class T>
struct OvertakeCtmcState {
    T idle, phase1, phase2;
};

/**
 * Stationary law of the reduced overtaking chain of overtake_prob.m.
 *
 * The chain is the cycle idle -> phase 1 -> phase 2 -> idle, with rates lambda,
 * 1/S1 and 1/S2. A cycle visits every state exactly once per traversal, so the
 * stationary probabilities are proportional to the mean holding times
 * (1/lambda, S1, S2); the reference reaches the same numbers by solving the
 * augmented singular system in least squares, which is the same answer arrived
 * at less directly and only for a floating T. Scaling the ratios by lambda
 * removes the reciprocal, so nothing here divides by a rate.
 *
 * The caller must have established lambda > 0; at lambda = 0 the chain is
 * absorbed in `idle` and has no unique stationary law.
 */
template <class T>
OvertakeCtmcState<T> lqn_overtake_ctmc(const T& S1, const T& S2, const T& lambda) {
    const T one = num_traits<T>::from_int(1);
    const T den = T(one + lambda * (S1 + S2));
    OvertakeCtmcState<T> pi;
    pi.idle = T(one / den);
    pi.phase1 = T(lambda * S1 / den);
    pi.phase2 = T(lambda * S2 / den);
    return pi;
}

/**
 * Probability that an arrival at an entry finds the server in phase 2.
 *
 * Port of @@SolverLN/overtake_prob.m. S1 and S2 are the entry's phase-1 and
 * phase-2 service times, lambda the arrival rate at the entry (the reference
 * falls back to the task throughput when the entry's own is not yet resolved,
 * which is the caller's choice to make) and `mult` the task multiplicity.
 *
 * PASTA is what licenses reading the answer off the time-stationary law: the
 * arrival stream is Poisson, so arrivals see time averages and the probability
 * an arrival finds phase 2 is the probability the chain is in phase 2.
 *
 * The multi-server branch is an APPROXIMATION in the reference and stays one
 * here (overtake_prob.m, lines 77 to 91): the phase-2 fraction of a busy server
 * scaled by the utilization, saturating at the bare phase-2 fraction once the
 * offered load reaches one server's worth. It is not a c-server chain and does
 * not converge to one.
 */
template <class T>
T lqn_overtake_prob(const T& S1, const T& S2, const T& lambda, double mult) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    // Each of these makes the phase-2 interval unobservable rather than small,
    // so the chain is not merely near-degenerate, it does not exist.
    if (num_traits<T>::to_double(S2) < GlobalConstants::FineTol ||
        num_traits<T>::to_double(S1) < GlobalConstants::FineTol ||
        num_traits<T>::to_double(lambda) < GlobalConstants::FineTol)
        return zero;

    auto clamp01 = [&](const T& x) { return x < zero ? zero : (x > one ? one : x); };

    if (mult == 1.0) return clamp01(lqn_overtake_ctmc(S1, S2, lambda).phase2);

    // An infinite server is the c -> Inf limit of the same expression: the load
    // per server vanishes, so an arrival never meets a busy one.
    if (!std::isfinite(mult)) return zero;
    const T phase_frac = T(S2 / (S1 + S2));
    const T rho = T(lambda * (S1 + S2) / num_traits<T>::from_double(mult));
    if (rho >= one) return clamp01(phase_frac);
    return clamp01(T(phase_frac * rho));
}

// ---------------------------------------------------------------------------
// overtake_prob_markov
// ---------------------------------------------------------------------------

namespace detail {

/** Entry that owns an activity, via actsof; 0 when the activity has no entry. */
template <class T>
std::size_t ln_entry_of_activity(const LqnStruct<T>& lqn, std::size_t aidx) {
    for (std::size_t e = 1; e <= lqn.nentries; ++e) {
        const std::size_t eabs = lqn.eshift + e;
        for (std::size_t a : lqn.actsof[eabs])
            if (a == aidx) return eabs;
    }
    return 0;
}

}  // namespace detail

/**
 * Overtaking probability at a server entry, through the LQNS phased-server
 * chain rather than the reduced CTMC above.
 *
 * Port of @@SolverLN/overtake_prob_markov.m: the input-mapping layer that turns
 * the LayeredNetworkStruct plus the current fixed-point iterate into the
 * per-client-phase slice parameters lqn_overtake_markov consumes, one client
 * entry at a time, summing the contributions and truncating to 1 as LQNS does
 * in Markov_Phased_Server::PrOT_e.
 *
 * `xj` is the tested server phase's residence time, `self.servt_ph2(eidx)` in
 * the reference. It is a parameter because the C++ SolverLN has no phase split
 * to read it from; see the file header.
 *
 * `servt` and `tput` are indexed by element 1..nidx and `callresidt` by call
 * 1..ncalls, which is exactly the layout SolverLN::state_servt, state_tput and
 * state_callresidt return.
 *
 * A trap worth naming: the reference's own comment says the caller activities
 * are the SYNCHRONOUS ones, but the selection it writes is unfiltered, so an
 * asynchronous or forwarding caller into the same entry also contributes a
 * client entry. That is reproduced here; the per-phase call scan below does
 * filter to SYNC, which is where the distinction actually bites.
 */
template <class T>
T lqn_overtake_prob_markov(const LqnStruct<T>& lqn, const std::vector<T>& servt,
                           const std::vector<T>& callresidt, const std::vector<T>& tput,
                           std::size_t eidx, const T& xj) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    if (servt.size() <= lqn.nidx || tput.size() <= lqn.nidx ||
        callresidt.size() <= lqn.ncalls)
        throw InputError(
            "lqn_overtake_prob_markov: servt and tput must be indexed 1..nidx and callresidt "
            "1..ncalls, as SolverLN's state accessors return them");
    if (eidx <= lqn.eshift || eidx > lqn.eshift + lqn.nentries)
        throw InputError("lqn_overtake_prob_markov: eidx is not an entry index");

    if (!(num_traits<T>::to_double(xj) > GlobalConstants::FineTol)) return zero;
    const std::size_t server_tidx = lqn.parent[eidx];

    // client entries that reach this server entry, in first-caller order
    std::vector<std::size_t> caller_entries;
    for (std::size_t c = 1; c <= lqn.ncalls; ++c) {
        if (lqn.callpair_dst[c] != eidx) continue;
        const std::size_t ceidx = detail::ln_entry_of_activity(lqn, lqn.callpair_src[c]);
        if (ceidx == 0) continue;
        bool known = false;
        for (std::size_t v : caller_entries)
            if (v == ceidx) known = true;
        if (!known) caller_entries.push_back(ceidx);
    }
    if (caller_entries.empty()) return zero;

    T prOt = zero;
    for (std::size_t ceidx : caller_entries) {
        const std::size_t ctidx = lqn.parent[ceidx];
        const std::vector<std::size_t>& acts = lqn.actsof[ceidx];
        if (acts.empty()) continue;

        int maxPhaseA = 1;
        for (std::size_t aidx : acts) {
            const std::size_t a = aidx - lqn.ashift;
            if (a >= 1 && a <= lqn.nacts && lqn.actphase[a] > maxPhaseA) maxPhaseA = lqn.actphase[a];
        }
        const std::size_t nStates = std::size_t(maxPhaseA) + 1;

        Matrix<T> clientPhases(nStates, 5, zero);
        clientPhases(0, 0) = one;  // the think slice is never chopped by a call
        const Distrib<T>& z = lqn.think[ctidx];
        if (!z.disabled && z.mean > zero) clientPhases(0, 1) = z.mean;

        std::vector<T> y_aj(nStates, zero);
        for (int p = 1; p <= maxPhaseA; ++p) {
            T nSlices = one, service = zero, y_ij = zero, y_ik = zero, tk_num = zero;
            for (std::size_t aidx : acts) {
                const std::size_t a = aidx - lqn.ashift;
                if (a < 1 || a > lqn.nacts || lqn.actphase[a] != p) continue;
                service = T(service + servt[aidx]);
                for (std::size_t c : lqn.callsof[aidx]) {
                    if (lqn.calltype[c] != CallType::SYNC) continue;
                    const T y = lqn.callproc_mean[c];
                    if (y == zero) continue;
                    // a rendezvous cuts the phase into one more slice, which is
                    // what makes the client interruptible partway through it
                    nSlices = T(nSlices + y);
                    if (lqn.parent[lqn.callpair_dst[c]] == server_tidx) {
                        y_ij = T(y_ij + y);
                    } else {
                        y_ik = T(y_ik + y);
                        tk_num = T(tk_num + y * callresidt[c]);
                    }
                }
            }
            const T t_k = y_ik > zero ? T(tk_num / y_ik) : zero;
            const std::size_t row = std::size_t(p);
            clientPhases(row, 0) = nSlices;
            clientPhases(row, 1) = service;
            clientPhases(row, 2) = y_ij;
            clientPhases(row, 3) = y_ik;
            clientPhases(row, 4) = t_k;
            y_aj[row] = y_ij;
            y_aj[0] = T(y_aj[0] + y_ij);
        }
        if (y_aj[0] == zero) continue;  // this client never reaches the server task

        T prVisit = one;
        if (num_traits<T>::to_double(tput[ctidx]) > GlobalConstants::FineTol &&
            num_traits<T>::to_double(tput[ceidx]) > GlobalConstants::FineTol)
            prVisit = T(tput[ceidx] / tput[ctidx]);

        prOt = T(prOt + lqn_overtake_markov(clientPhases, prVisit, xj, y_aj));
    }

    // Contributions of independent client entries are summed, so nothing stops
    // the total exceeding one; LQNS truncates rather than renormalizing.
    if (prOt < zero) return zero;
    return prOt > one ? one : prOt;
}

// ---------------------------------------------------------------------------
// convergedStoch
// ---------------------------------------------------------------------------

// LnStochConfig itself lives in solver_ln.h: `SolverLN::iterate` names it by
// value, and a non-dependent name must be complete where the template is
// parsed, not where it is instantiated. Only the controller below is deferred.

namespace detail {

/**
 * Running mean prev + (raw - prev)/k, tolerant of a NaN on either side.
 *
 * A stochastic layer solver can return NaN for a metric it did not estimate
 * (an idle class in a short simulation run), and a single such sample would
 * otherwise poison the average for every remaining iteration.
 */
template <class T>
Matrix<T> ln_polyak(const Matrix<T>& prev, const Matrix<T>& raw, long k) {
    Matrix<T> out(prev.rows(), prev.cols());
    const T kk = num_traits<T>::from_int(k);
    for (std::size_t i = 0; i < prev.rows(); ++i)
        for (std::size_t j = 0; j < prev.cols(); ++j) {
            T m = T(prev(i, j) + (raw(i, j) - prev(i, j)) / kk);
            if (std::isnan(num_traits<T>::to_double(m))) m = raw(i, j);
            if (std::isnan(num_traits<T>::to_double(m))) m = prev(i, j);
            out(i, j) = m;
        }
    return out;
}

template <class T>
std::vector<T> ln_polyak_vec(const std::vector<T>& prev, const std::vector<T>& raw, long k) {
    std::vector<T> out(prev.size());
    const T kk = num_traits<T>::from_int(k);
    for (std::size_t i = 0; i < prev.size(); ++i) {
        T m = T(prev[i] + (raw[i] - prev[i]) / kk);
        if (std::isnan(num_traits<T>::to_double(m))) m = raw[i];
        if (std::isnan(num_traits<T>::to_double(m))) m = prev[i];
        out[i] = m;
    }
    return out;
}

}  // namespace detail

/**
 * Convergence controller for an ensemble whose layers are solved by a NOISY
 * method (simulation, or Monte Carlo normalizing constants).
 *
 * Port of @@SolverLN/convergedStoch.m. The deterministic test in
 * SolverLN::converged cannot terminate against noise: the successive-difference
 * error is bounded below by the standard error of the layer estimates, and the
 * layer-reset confirmation step only resamples that noise. This replaces it
 * with a stochastic approximation scheme:
 *
 *   1. Burn-in. Plain Picard at the relaxation init chose, to get near the
 *      fixed point fast while the noise still does not matter.
 *   2. Robbins-Monro. `relax_omega` then decays as a0/k^alpha, so the iterate
 *      converges almost surely under the contraction assumption the
 *      deterministic iteration already makes plus zero-mean bounded-variance
 *      noise (Robbins and Monro, 1951). The caller must actually APPLY
 *      relax_omega to the fed-forward iterate for any of this to hold.
 *   3. Polyak-Ruppert. Running averages of the layer results and of the
 *      reported iterate, which give the optimal O(1/sqrt(k)) rate and make the
 *      answer insensitive to a0 (Polyak and Juditsky, 1992).
 *   4. Stopping on the DRIFT OF THE AVERAGE, not of the iterate. That drift
 *      decays like 1/k even under persistent noise, so the test terminates, and
 *      it self-calibrates: noisier layers hold the drift above tolerance longer
 *      and buy themselves more averaging.
 *
 * The reference averages QN, UN, RN, TN, AN and WN; LayerResult carries no AN,
 * so five fields are averaged here and the sixth is not silently invented.
 */
template <class T>
class LnStochController {
public:
    explicit LnStochController(const LnStochConfig& cfg)
        : cfg_(cfg), omega_(cfg.relax_burnin), err_(1, 0.0) {}

    /**
     * Fold iteration `it` in and say whether the iteration may stop. `it` counts
     * from 1 and must advance by one per call. `layer_jobs` is the total closed
     * population of each layer, which normalizes the drift so that layers of
     * very different size contribute comparably.
     */
    bool update(int it, const std::vector<LayerResult<T>>& latest,
                const std::vector<double>& layer_jobs, const std::vector<T>& servt,
                const std::vector<T>& residt) {
        if (it < 1) return false;
        if (err_.size() <= std::size_t(it)) err_.resize(std::size_t(it) + 1, 0.0);

        // Scheduled one iteration ahead, as in the reference: the step this sets
        // is the one the NEXT updateMetrics applies.
        if (it >= cfg_.burnin)
            omega_ = std::min(1.0, cfg_.a0 / std::pow(std::max(1.0, double(it - cfg_.burnin + 1)),
                                                      cfg_.alpha));

        if (it <= cfg_.burnin) {
            err_[std::size_t(it)] = std::numeric_limits<double>::infinity();
            return false;
        }
        if (start_ < 0) start_ = it;

        const long k = k_ + 1;
        double err = 0.0;
        if (k == 1) {
            avg_ = latest;
        } else {
            for (std::size_t e = 0; e < latest.size() && e < avg_.size(); ++e) {
                const LayerResult<T> prev = avg_[e];
                avg_[e].QN = detail::ln_polyak(prev.QN, latest[e].QN, k);
                avg_[e].UN = detail::ln_polyak(prev.UN, latest[e].UN, k);
                avg_[e].RN = detail::ln_polyak(prev.RN, latest[e].RN, k);
                avg_[e].TN = detail::ln_polyak(prev.TN, latest[e].TN, k);
                avg_[e].WN = detail::ln_polyak(prev.WN, latest[e].WN, k);
                const double N = e < layer_jobs.size() ? layer_jobs[e] : 0.0;
                if (N > 0.0) {
                    double dmax = 0.0;
                    for (std::size_t i = 0; i < prev.QN.rows(); ++i)
                        for (std::size_t j = 0; j < prev.QN.cols(); ++j) {
                            const double d = std::abs(num_traits<T>::to_double(avg_[e].QN(i, j)) -
                                                      num_traits<T>::to_double(prev.QN(i, j)));
                            if (!std::isnan(d) && d > dmax) dmax = d;
                        }
                    err += dmax / N;
                }
            }
        }
        k_ = k;

        if (k == 1) {
            servt_avg_ = servt;
            residt_avg_ = residt;
        } else {
            servt_avg_ = detail::ln_polyak_vec(servt_avg_, servt, k);
            residt_avg_ = detail::ln_polyak_vec(residt_avg_, residt, k);
        }

        err_[std::size_t(it)] = err;
        if (k <= cfg_.conseq) return false;
        for (long w = 0; w < cfg_.conseq; ++w)
            if (!(err_[std::size_t(it - w)] < cfg_.iter_tol)) return false;
        return true;
    }

    double relax_omega() const { return omega_; }
    /** Per-iteration drift, 1-based; slot 0 is unused. */
    const std::vector<double>& iteration_error() const { return err_; }
    long averaging_count() const { return k_; }
    /** Iteration at which averaging started, -1 while still in burn-in. */
    long averaging_start() const { return start_; }
    const std::vector<LayerResult<T>>& averaged_results() const { return avg_; }
    const std::vector<T>& averaged_servt() const { return servt_avg_; }
    const std::vector<T>& averaged_residt() const { return residt_avg_; }

private:
    LnStochConfig cfg_;
    double omega_;
    std::vector<double> err_;
    long k_ = 0;
    long start_ = -1;
    std::vector<LayerResult<T>> avg_;
    std::vector<T> servt_avg_, residt_avg_;
};

// ---------------------------------------------------------------------------
// the remaining @@SolverLN entry points, as free functions over the solver
// ---------------------------------------------------------------------------

/**
 * @@SolverLN/getSensitivityTable.m: solve the ensemble, then concatenate each
 * LAYER solver's own sensitivity table under a leading Layer column.
 *
 * The work is `SolverLN::get_sensitivity_table`, which has to be a member -- it
 * perturbs each layer in place and re-enters `solve_layer` for that layer, so it
 * needs the fork views, the region routing and the cache refresh that only the
 * solver holds. This wrapper exists so that the operation is reachable under
 * the name the reference gives it.
 */
template <class T>
LnSensTable<T> lqn_sensitivity_table(SolverLN<T>& solver, const sens::SensOptions& opt) {
    return solver.get_sensitivity_table(opt);
}

/**
 * @@SolverLN/getTranAvg.m: the block-diagonal aggregate transient over the LQN
 * layers, in whichever coupling `LnOptions::ln_transient` names.
 */
template <class T>
LnTranSolution lqn_tran_avg(SolverLN<T>& solver) {
    return solver.get_tran_avg();
}

/**
 * @@SolverLN/getCdfRespT.m: the per-entry response-time distribution, which
 * only the `moment3` method produces.
 */
template <class T>
std::vector<LnCdf> lqn_cdf_respt(SolverLN<T>& solver) {
    return solver.get_cdf_respt();
}

}  // namespace ln
}  // namespace line

#endif  // LINE_SOLVERS_LN_LQN_ANALYZERS_H
