#ifndef LINE_API_MC_CTMC_MEMORY_GATE_H
#define LINE_API_MC_CTMC_MEMORY_GATE_H

/**
 * @file
 * Host-aware memory pre-gate for SolverCTMC.
 *
 * Port of `matlab/src/api/mc/ctmc_memory_gate.m`, mirrored by the JAR
 * (`MemoryGuard.gate`) and native Python (`ctmc_memory_gate`).
 *
 * THE POLICY IS REFUSAL, NOT A WARNING. When the predicted peak exceeds the
 * safe budget the solve is REFUSED with an error naming the two numbers, and
 * only an explicit `force` downgrades that to a warning. A warning the caller
 * can walk past is the wrong contract here: the failure mode it precedes is the
 * OOM killer taking the whole process, which produces no diagnostic at all and,
 * observed on 2026-08-01, took down five sessions and the terminal hosting
 * them. Refusing costs the caller a re-run with another solver; not refusing
 * costs it everything running on the machine.
 *
 * MDD AND CFTP ARE NEVER GATED. Both are served by their own analyzers and
 * return before any explicit-generator path, so the gate is not on their route.
 * That exemption is load-bearing rather than incidental: `mdd` holds the
 * reachable set in a decision diagram whose size is governed by the diagram's
 * compression, not by the state count, so it routinely solves models this
 * estimator scores at exp(200). Gating on explicit size would clamp exactly the
 * method that exists to beat it. Keep any new non-explicit method on the same
 * side of the gate.
 *
 * The predictor is the power law of the reference, `bytes = alpha * N^beta`,
 * evaluated in log space. MATLAB and Python calibrate alpha/beta once per host
 * by profiling a sparse LU and cache the fit; this port uses the reference's
 * FALLBACK coefficients unconditionally, which is the same branch those two
 * take when calibration cannot run. That is deliberate: the fallback is the
 * conservative end of the fitted range, it is deterministic across hosts and
 * runs, and a gate that silently changes its verdict after a background
 * profiling step is worse than one that is slightly pessimistic.
 */

#include <cmath>
#include <sstream>
#include <string>

#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#endif

namespace line {
namespace mc {

/// 8 bytes of value plus 8 amortized for the index, per stored nonzero.
constexpr double CTMC_BYTES_PER_NZ = 16.0;
/// Fraction of available memory the solver may target.
constexpr double CTMC_DEFAULT_SAFETY_FRACTION = 0.6;
/// Fallback power-law coefficients, identical to MATLAB and Python.
constexpr double CTMC_FALLBACK_ALPHA = CTMC_BYTES_PER_NZ * 8.0;
constexpr double CTMC_FALLBACK_BETA = 1.3;
/// Conservative available-memory default when the host probe fails.
constexpr double CTMC_FALLBACK_AVAIL_BYTES = 1.0 * 1024.0 * 1024.0 * 1024.0;

/** The gate verdict, plus the message the caller reports either way. */
struct CtmcGateResult {
    bool ok = true;
    std::string message;
    double predicted_gb = 0.0;
    double budget_gb = 0.0;
};

/**
 * Available physical memory in bytes.
 *
 * Never throws: a probe failure returns the conservative fallback, because a
 * gate that cannot measure the host must not thereby become permissive.
 */
inline double ctmc_available_memory_bytes() {
#if defined(_SC_AVPHYS_PAGES) && defined(_SC_PAGE_SIZE)
    const long pages = ::sysconf(_SC_AVPHYS_PAGES);
    const long page_size = ::sysconf(_SC_PAGE_SIZE);
    if (pages > 0 && page_size > 0)
        return static_cast<double>(pages) * static_cast<double>(page_size);
#endif
    return CTMC_FALLBACK_AVAIL_BYTES;
}

/**
 * Decide whether a state space of log-size `log_nstates` can be solved here.
 *
 * @param log_nstates    natural log of the worst-case state count
 * @param force          true downgrades a refusal to a warning
 * @param safety_fraction fraction of available memory the solve may target
 */
inline CtmcGateResult ctmc_memory_gate(double log_nstates, bool force = false,
                                       double safety_fraction = CTMC_DEFAULT_SAFETY_FRACTION) {
    CtmcGateResult res;
    const double avail = ctmc_available_memory_bytes();
    const double budget = safety_fraction * avail;

    // Compared in LOG space: the predicted byte count of an intractable model
    // overflows a double long before it exceeds the budget, and inf > x is a
    // comparison that has already lost the margin it was meant to report.
    const double log_pred = std::log(CTMC_FALLBACK_ALPHA) + CTMC_FALLBACK_BETA * log_nstates;
    const double log_budget = std::log(budget > 1.0 ? budget : 1.0);

    res.predicted_gb = std::exp(log_pred < 700.0 ? log_pred : 700.0) / (1024.0 * 1024.0 * 1024.0);
    res.budget_gb = budget / (1024.0 * 1024.0 * 1024.0);
    if (log_pred <= log_budget) return res;

    // An intractable model predicts a byte count with hundreds of digits, and
    // printing it in full says nothing a magnitude does not. Above a petabyte
    // the figure is reported as a power of ten.
    std::ostringstream pred;
    const double log10_gb = (log_pred - std::log(1024.0 * 1024.0 * 1024.0)) / std::log(10.0);
    if (log10_gb > 6.0) {
        pred.setf(std::ios::fixed);
        pred.precision(1);
        pred << "1e" << log10_gb;
    } else {
        pred.setf(std::ios::fixed);
        pred.precision(2);
        pred << res.predicted_gb;
    }

    std::ostringstream os;
    os.setf(std::ios::fixed);
    os.precision(2);
    os << "CTMC predicted peak memory ~" << pred.str() << " GB exceeds the safe budget ~"
       << res.budget_gb << " GB (" << static_cast<int>(100.0 * safety_fraction + 0.5) << "% of "
       << avail / (1024.0 * 1024.0 * 1024.0)
       << " GB available). Reduce the state space (e.g. lower 'cutoff'), use the 'mdd' method, "
          "use another solver (MVA/NC/FLD), or set force=true to override.";
    res.message = os.str();
    res.ok = force;
    return res;
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_MEMORY_GATE_H
