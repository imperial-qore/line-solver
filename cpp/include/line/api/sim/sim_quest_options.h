/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SIM_SIM_QUEST_OPTIONS_H
#define LINE_API_SIM_SIM_QUEST_OPTIONS_H

/**
 * Options of the QUEST procedures, with the published FQUEST defaults.
 *
 * Port of matlab/src/api/sim/sim_quest_options.m. The defaults b0 = 50, m0 = 500,
 * s = [32 24 16 10], beta = 0.30, eta = 0.2 and theta = 2.3 are the ones the
 * article reports after its own experimentation: b0 = 50 gives the warmup
 * randomness test enough power, 32 batches suffice to estimate the variance
 * parameter while fewer than 10 make the interval unreliable, and the decaying
 * warmup significance keeps the batch size from growing so far that truncation
 * eats a short sample. With these values the fourth warmup iteration runs at
 * beta*exp(-0.2*3^2.3) = 0.025.
 *
 * The reference also rejects unknown field names, which a struct cannot carry;
 * what remains here is the admissibility check, and it is not cosmetic. In
 * particular s must be STRICTLY DECREASING: sim_fquest walks it as a
 * batch-count ladder that only ever steps down, so a non-monotone s would let a
 * later stage re-test a batch count an earlier stage had already rejected and
 * the loop would no longer terminate at the ladder's end.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/util/error.h"

namespace line {
namespace sim {

/** Procedure constants shared by sim_fquest and sim_firquest. */
struct QuestOptions {
    long b0 = 50;                        ///< Initial batch count for the warmup stage
    long m0 = 500;                       ///< Initial batch size for the warmup stage
    std::vector<long> s{32, 24, 16, 10}; ///< Descending batch counts for the test stages
    double beta = 0.30;                  ///< Significance level of the stage tests
    double eta = 0.2;                    ///< Decay coefficient of the warmup significance
    double theta = 2.3;                  ///< Decay exponent of the warmup significance
    double weight = std::sqrt(12.0);     ///< Constant STS weight function
    bool force = true;                   ///< Deliver a heuristic interval when a test fails
};

/**
 * Validates an option set and returns it.
 *
 * @param options the constants to check
 * @return the same constants, once every admissibility condition holds
 */
inline QuestOptions sim_quest_options(const QuestOptions& options = QuestOptions()) {
    QuestOptions opt = options;
    if (opt.b0 < 3) throw InputError("sim_quest_options: b0 must be an integer >= 3");
    if (opt.m0 < 1) throw InputError("sim_quest_options: m0 must be a positive integer");
    if (opt.s.empty())
        throw InputError("sim_quest_options: s must be a nonempty vector of positive integers");
    for (std::size_t i = 0; i < opt.s.size(); ++i)
        if (opt.s[i] < 1)
            throw InputError("sim_quest_options: s must be a nonempty vector of positive integers");
    for (std::size_t i = 1; i < opt.s.size(); ++i)
        if (opt.s[i] >= opt.s[i - 1])
            throw InputError("sim_quest_options: s must be strictly decreasing");
    if (!(opt.beta > 0.0) || !(opt.beta < 1.0))
        throw InputError("sim_quest_options: beta must be a real scalar in (0,1)");
    if (!(opt.eta >= 0.0))
        throw InputError("sim_quest_options: eta must be a nonnegative real scalar");
    if (!(opt.theta > 0.0))
        throw InputError("sim_quest_options: theta must be a positive real scalar");
    if (opt.weight == 0.0)
        throw InputError("sim_quest_options: weight must be a nonzero real scalar");
    return opt;
}

}  // namespace sim
}  // namespace line

#endif  // LINE_API_SIM_SIM_QUEST_OPTIONS_H
