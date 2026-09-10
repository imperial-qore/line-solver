#pragma once
/**
 * @file fluid_conservation_guard.h
 * @brief Detects a moment-closure trajectory that has left the model.
 *
 * WHY IT EXISTS. The moment-closure drift can leave the simplex: on a station
 * where min(n,c) is not the identity the Gaussian correction to the per-class
 * share can drive a coordinate negative, and since the drift is conservative
 * another grows to match. In MATLAB `odeset('NonNegative')` projects the
 * ACCEPTED step, so the excursion is CLAMPED rather than reported -- which
 * injects mass, collapses the step size, and leaves the window never returning.
 * One MATLAB suite run sat in `test_CQN_Cox_CS_9` for 3h16m, and the 2026-08-27
 * run was killed after `test11_interlock_lqnx` had held the suite for 100
 * minutes, taking every block after it down with it.
 *
 * THIS PORT CANNOT HANG THE SAME WAY, and the difference is worth stating
 * rather than papering over: `solver_fluid.h` clamps the state AFTER each
 * window, not inside the integration, so the same divergence surfaces as a
 * finished window holding a state that is not a solution -- silently, with
 * every later window integrated from it. That is what this makes loud. Python
 * DOES reproduce the reference's semantics (`_integrate_nonnegative` clips and
 * restarts), so its guard runs per accepted step and halts the window.
 *
 * THE TEST IS AN EXACT INVARIANT, not a heuristic bound on time or magnitude.
 * The drift conserves the population of every CLOSED CHAIN exactly, so any
 * deviation is a divergence and nothing else. The tolerance is a generous
 * fraction of that population rather than a numerical tolerance: the
 * integrator's own error is ~1e-4 relative, while the documented excursion
 * reaches 5.2e4 against a true population of 0.05. A closed model whose
 * population has moved by TOL is no longer solving the model, whatever it is
 * converging to.
 *
 * THE CHAIN IS THE CONSERVED UNIT, NOT THE CLASS, and the difference is the
 * whole correctness of this check. A class population is what that class
 * STARTS with; class switching then moves jobs between the classes of one
 * chain, so only the chain total is invariant. Watching classes instead
 * condemns every class-switching model out of hand -- measured on
 * `cqn_twoclass_hyperl` (313 of 447 accepted states), on `init_state_ps` (286
 * of 310) and on every one of the 162 fluid layers an LQN builds under the
 * `srvn.cs` encoding, where the chain sum never moved at all. A cache model is
 * the same story with the hit/miss classes.
 *
 * A wall-clock budget would have caught the same thing and was rejected: it
 * makes the answer depend on how busy the host is, so the same model would fall
 * back on one machine and not on another. This invariant is deterministic.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/solvers/fluid/fluid_nonhyperbolic.h"
#include "line/solvers/fluid/fluid_odes.h"

namespace line {
namespace fluid {

/** Relative population drift that counts as having left the model. */
inline constexpr double kFluidConservationTol = 0.1;

// The "is the closure active" predicate is FluidClosure::gaussian(), which
// already exists: the first pass of the closure runs at sigma2 = 0 -- that pass
// IS the first-order solve -- and only the later ones can diverge. Gating on it
// keeps the guard off `closing` and `matrix`, which are also the fallback
// ladder's own rungs and must not be sent down a fallback by their own
// watchdog.

/**
 * The classes of each chain, as 0-based column indices.
 *
 * `sn.inchain` holds them 1-based, one vector per chain. A struct declaring no
 * chains degrades to one chain per class, which is the safe reading rather
 * than a guess: with no chain map there is no class switching to merge
 * classes, so each class IS its own conserved unit.
 */
template <typename T>
inline std::vector<std::vector<std::size_t>> fluid_chain_partition(
        const qn::NetworkStruct<T>& sn, std::size_t K) {
    std::vector<std::vector<std::size_t>> out;
    if (sn.inchain.empty()) {
        out.reserve(K);
        for (std::size_t r = 0; r < K; ++r) out.push_back({r});
        return out;
    }
    out.reserve(sn.inchain.size());
    for (const std::vector<std::size_t>& chain : sn.inchain) {
        std::vector<std::size_t> members;
        members.reserve(chain.size());
        for (std::size_t r1 : chain)
            if (r1 >= 1 && r1 - 1 < K) members.push_back(r1 - 1);
        out.push_back(members);
    }
    return out;
}

/**
 * The closed chain whose conserved population has drifted past `tol`, or -1.
 *
 * @param sn  the network struct, for the chain membership and populations
 * @param L   the state layout, which fixes each (station,class) block
 * @param x   a state vector of length `L.nstates`
 * @param tol relative deviation that counts as having left the model
 */
template <typename T>
int fluid_conservation_violation(const qn::NetworkStruct<T>& sn,
                                 const FluidLayout& L,
                                 const std::vector<double>& x,
                                 double tol = kFluidConservationTol) {
    const std::size_t M = L.qidx.size();
    if (M == 0) return -1;
    const std::size_t K = L.qidx[0].size();
    const std::vector<std::vector<std::size_t>> chains = fluid_chain_partition(sn, K);
    for (std::size_t c = 0; c < chains.size(); ++c) {
        double target = 0.0;
        for (std::size_t r : chains[c])
            target += static_cast<double>(sn.classes[r].population);
        if (!std::isfinite(target) || target <= 0.0) {
            continue;  // open, or absent: no conserved population to check
        }
        double mass = 0.0;
        bool any = false;
        for (std::size_t r : chains[c]) {
            for (std::size_t i = 0; i < M; ++i) {
                const std::size_t n = L.kic[i][r];
                if (n == 0) continue;
                const std::size_t base = L.qidx[i][r];
                if (base + n > x.size()) continue;  // a different state vector
                for (std::size_t k = 0; k < n; ++k) mass += x[base + k];
                any = true;
            }
        }
        if (!any) continue;
        if (std::fabs(mass - target) > tol * std::max(1.0, target)) {
            return static_cast<int>(c);
        }
    }
    return -1;
}

}  // namespace fluid
}  // namespace line
