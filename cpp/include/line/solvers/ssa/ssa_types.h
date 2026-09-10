/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_SSA_SSA_TYPES_H
#define LINE_SOLVERS_SSA_SSA_TYPES_H

/**
 * Controls, results and the random source of SolverSSA.
 *
 * WHICH RANDOM STREAM, AND WHAT THAT COSTS. An SSA answer is a function of the
 * random stream, so a seed-fixed golden belongs to one implementation only.
 * MATLAB draws from `rand` (its own Mersenne Twister wrapper); the JAR and
 * native Python share an MT19937. This port uses MT19937 too, but it does NOT
 * claim stream compatibility with any of them: the ORDER in which the sample
 * path consumes draws is part of the algorithm, and no two of the four
 * implementations consume in the same order (the JAR, for instance, spends a
 * draw on a single-phase entry-phase pick where this port does not). So a
 * cross-codebase check against this engine is STATISTICAL, never exact --
 * which is what `_kb/14-cpp-multiprecision.md` records and what the tests
 * assert.
 */

#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include "line/util/matrix.h"

namespace line {
namespace ssa {

/** Helpers shared by the SSA engines that do not belong to a single one. */
namespace ssa_detail {

/**
 * Add the START/PREEMPT counts of one successor row to the per-state rate
 * accumulators, weighted by the rate of the arc that carries them. A node that
 * is not a station contributes nothing: only a station has a server to seize.
 */
template <class T, class Outcome>
inline void add_tag_rates(std::vector<std::vector<double>>& start,
                          std::vector<std::vector<double>>& preempt,
                          const T& sn, std::size_t node, const Outcome& oc, std::size_t row,
                          double w) {
    if (!(w != 0) || node == 0 || node > sn.nodes.size()) return;
    const std::size_t isf = sn.stateful_index(node);
    if (isf == 0 || isf > start.size()) return;
    if (row < oc.start.size())
        for (std::size_t j = 0; j < oc.start[row].size(); ++j) {
            const std::size_t cls = oc.start[row][j];
            if (cls >= 1 && cls <= start[isf - 1].size()) start[isf - 1][cls - 1] += w;
        }
    if (row < oc.preempt.size())
        for (std::size_t j = 0; j < oc.preempt[row].size(); ++j) {
            const std::size_t cls = oc.preempt[row][j];
            if (cls >= 1 && cls <= preempt[isf - 1].size()) preempt[isf - 1][cls - 1] += w;
        }
}

}  // namespace ssa_detail

/** Controls, defaulting to `SolverOptions('SSA')` in the reference. */
struct SsaOptions {
    /** `default` and `nrm` both select the Next Reaction Method here. */
    std::string method = "default";
    /** Reaction firings to simulate; `options.samples` in the reference. */
    std::size_t samples = 10000;
    /** `options.seed`; LINE's own default is 23000. */
    unsigned long seed = 23000;
    bool verbose = false;
    /**
     * `options.config.warmupfrac`: the leading fraction of the path discarded
     * before the means are taken.
     *
     * DECLARED HERE, NOT ON `SsaSerialOptions`, although only the serial and
     * replicated engines read one: `ssa_serial_options` builds the engine's
     * knobs by copying the BASE slice of the caller's, so a field on the
     * derived struct is unreachable from a caller and kept its default however
     * the request was spelled -- which is what made `--warmupfrac` a knob the
     * C++ CLI could not offer at all.
     */
    double warmupfrac = 0.0;
    /**
     * `options.config.state_space_gen` of `solver_ssa_analyzer_nrm.m`: which of
     * the two NRM engines runs. `none` and `default` take the plain engine,
     * which integrates the metrics along the path; anything else takes the
     * tabulating one, which records the states it visits and forms the means as
     * `pi * A` over them. The two answer the same question by different routes,
     * which is what makes them testable against each other.
     */
    std::string state_space_gen = "default";
};

/** What the analyzer returns, in the same shape as the MVA and fluid results. */
struct SsaSolution {
    Matrix<double> QN, UN, RN, TN;
    std::vector<double> CN, XN;
    /**
     * The DERIVED rates, (nstations x nclasses): how often per unit time a
     * class-r service STARTS at station i, and how often a class-r job in
     * service is PUSHED BACK into the buffer there. At a lossless station with
     * no in-service abandonment StartN == TN + PreemptN, up to simulation
     * error; SolverCTMC reports the exact value of the same quantity.
     */
    Matrix<double> StartN, PreemptN;
    /** The concrete algorithm, as the reference's `method`. */
    std::string method = "nrm";
    /** Simulated time the metrics are averaged over; the reference's `totalTime`. */
    double simulated_time = 0.0;
    /** Reaction firings actually performed. */
    std::size_t samples = 0;
};

/**
 * The uniform source, MATLAB's `rand`.
 *
 * A 53-bit uniform assembled from two MT19937 words, the construction numpy
 * and the JAR's RandomManager both use, shifted by half an ulp so the value is
 * strictly inside (0,1). The shift is not cosmetic: `-log(u)` is the
 * exponential clock of every reaction and an exact zero would make it
 * infinite, which the run loop reads as a deadlock.
 */
class SsaRng {
public:
    explicit SsaRng(unsigned long seed) : g_(static_cast<std::uint_fast32_t>(seed)) {}

    double uniform() {
        const std::uint64_t a = g_() >> 5, b = g_() >> 6;
        return (static_cast<double>(a) * 67108864.0 + static_cast<double>(b) + 0.5) /
               9007199254740992.0;
    }

    /** Uniform index in [0, n), the reference's `1 + floor(rand*n)`. */
    std::size_t index(std::size_t n) {
        if (n == 0) return 0;
        const std::size_t k = static_cast<std::size_t>(uniform() * static_cast<double>(n));
        return k < n ? k : n - 1;
    }

    /**
     * Index drawn from the unnormalized nonnegative weights `p`, the
     * reference's `drawFromDist`: an all-zero weight vector yields index 0.
     */
    std::size_t draw(const std::vector<double>& p) {
        double tot = 0.0;
        for (double x : p) tot += x;
        if (!(tot > 0.0)) return 0;
        const double u = uniform() * tot;
        double acc = 0.0;
        for (std::size_t i = 0; i < p.size(); ++i) {
            acc += p[i];
            if (acc > u) return i;
        }
        return p.size() - 1;
    }

private:
    std::mt19937 g_;
};

}  // namespace ssa
}  // namespace line

#endif  // LINE_SOLVERS_SSA_SSA_TYPES_H
