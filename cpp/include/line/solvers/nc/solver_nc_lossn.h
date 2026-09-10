/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_NC_SOLVER_NC_LOSSN_H
#define LINE_SOLVERS_NC_SOLVER_NC_LOSSN_H

/**
 * Port of `solver_nc_lossn_analyzer.m`: the open LOSS NETWORK, which is a
 * Source, ONE multiclass Delay sitting inside a Finite Capacity Region under a
 * DROP rule, and a Sink.
 *
 * WHAT THE ANALYZER ACTUALLY COMPUTES. There is no queueing: the single station
 * is an infinite server, so a job that is admitted never waits and leaves after
 * one service time. The only question is which arrivals are ADMITTED, and the
 * region answers it through a linear admission rule A n <= C on the per-class
 * occupancy vector n. Everything else -- carried throughput, mean population,
 * response time -- follows from the per-class blocking probability by Little's
 * law at the infinite server.
 *
 * WHERE THE ROWS OF A COME FROM, in the order the simulation engines test them:
 * the global job cap, the memory budget weighted by the per-class sizes, the
 * per-class job caps, and any explicit linear constraint. A row left unbounded
 * (the -1 sentinel) is DROPPED rather than given a surrogate capacity, because
 * a large finite surrogate would report a small but non-zero blocking where the
 * truth is none. A region that declares no bounded row at all has no admission
 * rule and is refused rather than solved as an unconstrained delay.
 *
 * Every row is a function of the per-class occupancy of the REGION only, which
 * is all the FiniteCapacityRegion API can express, so no row can distinguish
 * the stations inside the region. That is why a second member station would
 * change nothing in A and yet would break the Little's-law recovery below,
 * which attributes the whole carried load to one infinite server; a region with
 * more than one member station is therefore refused rather than collapsed.
 *
 * THE THREE METHODS, ALL PORTED.
 *
 *   exact / ms  `lossn_manjunath`, the Manjunath-Sikdar contour-integral transform,
 *               and the default on integral constraints as in the reference.
 *               EXACT: g(C) is a residue, i.e. a coefficient of a truncated
 *               multivariate power series, so the answer carries neither an
 *               iteration tolerance nor a sampling error. It is the only one of
 *               the three that runs under exact arithmetic, since every step is
 *               rational and the metrics are ratios. Cost is the product of
 *               (C_j+1) over the simultaneously live rows, so it is exact but
 *               not unconditionally cheap.
 *   erlangfp    `lossn_erlangfp`, the Erlang fixed-point (reduced-load)
 *               approximation. IT IS AN APPROXIMATION: it assumes the links
 *               block INDEPENDENTLY, which is false whenever two rows of A
 *               share a class, and it is exact only in the single-row case
 *               where the assumption is vacuous and the fixed point collapses
 *               to Erlang's loss formula. Never compare it to an exact solver
 *               at a tight tolerance on a multi-row region.
 *   mci         `lossn_mci`, Ross-Wang importance sampling. Unbiased, with a
 *               confidence interval, and it carries a normalizing constant.
 *               The method to reach for when the exact transform's live-grid
 *               product is prohibitive.
 *
 *   rec         `lossn_rec`, MDD-rec: the same constant as the exact sum over
 *               the admissible set, obtained by one memoised walk of the
 *               decision diagram holding it. It places no integrality demand on
 *               A or C, which is why it -- and not 'erlangfp', which this port
 *               refuses there -- is what 'default' takes on a FRACTIONAL region.
 *               Before it existed a fractional region had no exact route here at
 *               all, only the Monte Carlo 'mci' whose answer is a random
 *               variable.
 *
 * WHY 'erlangfp' NEEDS INTEGRAL A AND C AND 'mci' DOES NOT. The ported
 * `lossn_erlangfp` raises (1-E_i) to the power A(i,r) through `num_pow_int`,
 * which takes an unsigned exponent, and calls `erlang_b` with an `int`
 * capacity; a fractional entry would be TRUNCATED silently, solving a different
 * region. The reference evaluates both through `factln`, i.e. a gamma function,
 * so it accepts fractional arguments -- which is exactly why its 'default'
 * falls back to 'erlangfp' when the region declares fractional class sizes. The
 * port cannot follow that fallback and refuses it by name. `lossn_mci` sums
 * over an integer lattice of states but compares in real arithmetic, so it has
 * no such restriction and is the method to use on a fractional region.
 * `lossn_manjunath` needs integral rows for a different reason -- the residue argument
 * counts whole units of capacity -- so an EXPLICIT 'exact' on a fractional
 * region is refused by `lossn_manjunath` itself rather than downgraded here, exactly
 * as the reference does; only 'default' falls back.
 *
 * ARITHMETIC. There is no blanket transcendental gate: under exact arithmetic
 * the 'exact' method answers exactly and the other two refuse by name, which is
 * strictly more useful than refusing the whole analyzer. `lG` is a double in
 * every arithmetic, because it is a logarithm.
 */

#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "line/api/da/da_fpi.h"
#include "line/api/lossn/lossn_erlangfp.h"
#include "line/api/lossn/lossn_mci.h"
#include "line/api/lossn/lossn_manjunath.h"
#include "line/api/lossn/lossn_rec.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/nc/nc_types.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace nc {

/** What the loss-network analyzer returns beyond the usual metric table. */
template <class T>
struct NcLossnSolution {
    NcSolution<T> sol;
    Matrix<T> A;          ///< (J x K) rows of the admission rule A n <= C
    std::vector<T> Cvec;  ///< (J) the right-hand side
    std::vector<T> nu;    ///< (K) offered load per class
    std::vector<T> Loss;  ///< (K) blocking probability per class
    std::vector<T> E;     ///< (J) per-row blocking, empty outside the fixed point
};

namespace detail {

/**
 * The 1-based stations of region `f`.
 *
 * `members` is the authoritative flag and is read first. A struct assembled
 * without it -- which is what MATLAB's `sn.region{f}` is, a capacity matrix and
 * nothing else -- is read back the reference's way instead, from the -1
 * sentinel: a station whose row carries any bounded entry is a member. The two
 * agree on every region the builder produces; they differ only on a member
 * station left wholly unbounded, which the sentinel test cannot see.
 */
template <class T>
std::vector<std::size_t> lossn_region_members(const qn::NetworkStruct<T>& sn, std::size_t f) {
    const typename qn::NetworkStruct<T>::Region& rg = sn.regions[f];
    bool anyFlag = false;
    for (bool b : rg.members)
        if (b) anyFlag = true;

    std::vector<std::size_t> out;
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        bool in = false;
        if (anyFlag) {
            in = i < rg.members.size() && rg.members[i];
        } else if (i < rg.cap.size()) {
            for (double c : rg.cap[i])
                if (c >= 0.0) in = true;
            if (i < rg.maxmem.size() && rg.maxmem[i] >= 0.0) in = true;
        }
        if (in) out.push_back(i + 1);
    }
    return out;
}

/**
 * Port of the local `lossn_region_constraints` of the reference: the rows of
 * A n <= C for region `f` at member station `st` (1-based).
 */
template <class T>
void lossn_region_constraints(const qn::NetworkStruct<T>& sn, std::size_t f, std::size_t st,
                              std::size_t K, Matrix<T>& A, std::vector<T>& Cvec) {
    const typename qn::NetworkStruct<T>::Region& rg = sn.regions[f];
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::vector<double>& caps = rg.cap[st - 1];

    std::vector<std::vector<T>> rows;
    std::vector<T> rhs;

    // Global job cap: sum_r n_r <= globalMaxJobs.
    if (K < caps.size() && caps[K] >= 0.0) {
        rows.push_back(std::vector<T>(K, one));
        rhs.push_back(num_traits<T>::from_double(caps[K]));
    }

    // Memory budget: sum_r size_r n_r <= globalMaxMemory. The class sizes are
    // the row, so this is the one row that can legitimately be fractional.
    if (st - 1 < rg.maxmem.size() && rg.maxmem[st - 1] >= 0.0) {
        std::vector<T> row(K, one);
        for (std::size_t r = 0; r < K && r < rg.size.size(); ++r) row[r] = rg.size[r];
        rows.push_back(row);
        rhs.push_back(num_traits<T>::from_double(rg.maxmem[st - 1]));
    }

    // Per-class job caps, already folded with the per-class memory caps by
    // `add_region`, so there is no separate per-class memory row.
    for (std::size_t r = 0; r < K && r < caps.size(); ++r) {
        if (caps[r] < 0.0) continue;
        std::vector<T> row(K, zero);
        row[r] = one;
        rows.push_back(row);
        rhs.push_back(num_traits<T>::from_double(caps[r]));
    }

    // Explicit linear constraints from FiniteCapacityRegion.setConstraint.
    for (std::size_t k = 0; k < rg.lincon_A.rows(); ++k) {
        std::vector<T> row(K, zero);
        for (std::size_t r = 0; r < K && r < rg.lincon_A.cols(); ++r) row[r] = rg.lincon_A(k, r);
        rows.push_back(row);
        rhs.push_back(k < rg.lincon_b.size() ? rg.lincon_b[k] : zero);
    }

    if (rows.empty())
        throw UnsupportedError(
            "solver_nc_lossn_analyzer: the finite capacity region declares no bounded constraint, "
            "so it admits every arrival and is not a loss network; give it a global job cap, a "
            "memory budget, a per-class cap or an explicit linear constraint");

    A = Matrix<T>(rows.size(), K, zero);
    for (std::size_t j = 0; j < rows.size(); ++j)
        for (std::size_t r = 0; r < K; ++r) A(j, r) = rows[j][r];
    Cvec = rhs;
}

/** `options.method` split on '.' and '/' and lowercased, as `nc` does elsewhere. */
inline std::vector<std::string> lossn_tokens(const std::string& method) {
    std::vector<std::string> toks;
    std::string tok;
    for (char ch : method) {
        if (ch == '.' || ch == '/') {
            toks.push_back(tok);
            tok.clear();
        } else {
            tok += static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
        }
    }
    toks.push_back(tok);
    return toks;
}

inline bool lossn_has_token(const std::vector<std::string>& toks, const char* what) {
    for (const std::string& t : toks)
        if (t == what) return true;
    return false;
}

}  // namespace detail

/**
 * True when the model has the SHAPE of a loss network -- open, one region, one
 * member station, that station an infinite server -- whatever admission rule
 * the region applies.
 *
 * SEPARATE FROM THE RULE TEST ON PURPOSE. A shape test that also demanded DROP
 * would make a WAITQ region on this very shape indistinguishable, to the
 * caller, from a model with no region at all, and the runner would then answer
 * it as an unconstrained network. The reference splits the two questions the
 * same way: `runAnalyzer.m:307-327` first recognises the shape, and only then
 * branches on the rule, raising on WAITQ rather than falling through.
 *
 * The member test is `isinf(nservers)`, not `nodetype == Delay`, to match the
 * reference. An infinite-server Queue has the shape and must reach the analyzer
 * so its by-name refusal fires; excluding it here would silently hide it.
 */
template <class T>
bool nc_has_lossn_shape(const qn::NetworkStruct<T>& sn) {
    if (sn.regions.size() != 1) return false;
    const std::vector<std::size_t> mem = detail::lossn_region_members(sn, 0);
    if (mem.size() != 1) return false;
    const std::size_t st = mem[0];
    if (st == 0 || st > sn.stations.size()) return false;
    if (std::isfinite(sn.stations[st - 1].nservers)) return false;

    for (const qn::JobClass& c : sn.classes)
        if (std::isfinite(c.population) && c.population > 0.0) return false;
    return true;
}

/**
 * True when the model is a loss network: the shape above, with EVERY class
 * dropped at the region.
 *
 * The DROP rule is what makes it a LOSS network rather than a blocking one. A
 * WAITQ region holds the arrival back instead of discarding it, which is a
 * queueing phenomenon the Erlang model has no state for, so it must not be
 * routed here.
 *
 * ALL classes, not merely one: the reference tests `all(regionrule(1,:) ==
 * DROP)`. A region that discards one class and holds another back is a mixed
 * system whose blocked class occupies the region while it waits, so the
 * per-class loss probabilities the Erlang fixed point returns would not be the
 * ones the model implies.
 */
template <class T>
bool nc_is_lossn_model(const qn::NetworkStruct<T>& sn) {
    if (!nc_has_lossn_shape(sn)) return false;
    if (sn.regions[0].rule.size() < sn.nclasses) return false;
    for (std::size_t r = 0; r < sn.nclasses; ++r)
        if (sn.regions[0].rule[r] != lang::DropStrategy::DROP) return false;
    return true;
}

/**
 * Port of `solver_nc_lossn_analyzer.m`.
 *
 * @param sn  the refreshed struct; must satisfy `nc_is_lossn_model`
 * @param opt solver controls; `method` selects erlangfp / mci
 */
template <class T>
NcLossnSolution<T> solver_nc_lossn_analyzer(const qn::NetworkStruct<T>& sn,
                                            const NcSolverOptions& opt) {
    NcLossnSolution<T> out;
    {
        const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
        const std::size_t K = sn.nclasses, M = sn.nstations;

        // 1. The delay inside the region.
        const std::vector<std::size_t> mem = detail::lossn_region_members(sn, 0);
        if (sn.regions.empty())
            throw UnsupportedError(
                "solver_nc_lossn_analyzer: the model declares no finite capacity region, so it is "
                "not a loss network");
        if (mem.size() != 1)
            throw UnsupportedError(
                "solver_nc_lossn_analyzer: the finite capacity region holds " +
                std::to_string(mem.size()) +
                " stations; the admission rule can only see the occupancy of the region as a "
                "whole, so the carried load could not be attributed to a station");
        const std::size_t delayIdx = mem[0];
        if (sn.stations[delayIdx - 1].nodetype != qn::NodeType::Delay)
            throw UnsupportedError(
                "solver_nc_lossn_analyzer: the station inside the finite capacity region is not a "
                "Delay; a loss network holds admitted jobs at an infinite server, and a queueing "
                "station would make the response time depend on the population");

        // 2. Offered load. A route carries nu_r = arrival rate times mean
        // holding time INSIDE the region, i.e. the visit ratio at the delay over
        // its service rate. The bare arrival rate would be right only for unit
        // mean service times.
        out.nu.assign(K, zero);
        std::vector<T> lambda(K, zero), mu(K, zero);
        const std::size_t dsf = sn.stateful_of_station(delayIdx);
        for (std::size_t r = 0; r < K; ++r) {
            const std::size_t src = sn.classes[r].refstat;
            if (src == 0 || src > M) continue;
            if (!sn.disabled[src - 1][r]) lambda[r] = sn.rates(src - 1, r);
            // A class the delay never serves never enters the region: it offers
            // no load and is reported blocked with probability zero, which is
            // what a rate read out of a disabled pair could not express.
            if (sn.disabled[delayIdx - 1][r]) continue;
            mu[r] = sn.rates(delayIdx - 1, r);
            if (mu[r] == zero)
                throw UnsupportedError(
                    "solver_nc_lossn_analyzer: class " + std::to_string(r + 1) +
                    " is enabled at the delay with a zero service rate, so its mean holding time "
                    "in the region is unbounded and its offered load is undefined");

            T V = one;
            const std::size_t rsf = sn.stateful_of_station(src);
            for (std::size_t c = 0; c < sn.nchains; ++c) {
                if (!sn.chains[c][r]) continue;
                const T vref = sn.visits[c](rsf - 1, r);
                if (vref > zero) V = T(sn.visits[c](dsf - 1, r) / vref);
                break;
            }
            out.nu[r] = T(lambda[r] * V / mu[r]);
        }

        // 3. The admission rule.
        detail::lossn_region_constraints(sn, 0, delayIdx, K, out.A, out.Cvec);
        const std::size_t J = out.Cvec.size();

        // 4. Method selection, following the reference's precedence: an explicit
        // 'mci' wins, then an explicit 'erlangfp', then the exact transform.
        const std::vector<std::string> toks = detail::lossn_tokens(opt.method);
        bool integral = true;
        for (std::size_t j = 0; j < J; ++j) {
            const double c = num_traits<T>::to_double(out.Cvec[j]);
            if (std::fabs(c - std::round(c)) > 1e-9) integral = false;
            for (std::size_t r = 0; r < K; ++r) {
                const double a = num_traits<T>::to_double(out.A(j, r));
                if (std::fabs(a - std::round(a)) > 1e-9) integral = false;
            }
        }

        std::string chosen;
        if (detail::lossn_has_token(toks, "mci")) {
            chosen = "mci";
        } else if (detail::lossn_has_token(toks, "erlangfp")) {
            chosen = "erlangfp";
        } else if (detail::lossn_has_token(toks, "rec")) {
            chosen = "rec";
        } else if (detail::lossn_has_token(toks, "exact") ||
                   detail::lossn_has_token(toks, "manjunath") ||
                   detail::lossn_has_token(toks, "ms")) {
            // An EXPLICIT request for the transform is honoured even on a
            // fractional region, where `lossn_manjunath` refuses it by name; silently
            // answering with the Erlang approximation instead would report an
            // approximation under the name of an exact method.
            chosen = "exact";
        } else {
            // Only 'default' falls back, and only on a fractional region -- to
            // MDD-rec, which is exact there, rather than to 'erlangfp', which
            // this port refuses on a fractional region.
            chosen = integral ? "exact" : "rec";
        }

        double lG = std::numeric_limits<double>::quiet_NaN();
        std::vector<T> QLen(K, zero);
        out.Loss.assign(K, zero);
        int niter = 0;
        std::string method;

        if (chosen == "rec") {
            const lossn::LossnRecResult<T> rr = lossn::lossn_rec<T>(out.nu, out.A, out.Cvec);
            QLen = rr.QLen;
            out.Loss = rr.Loss;
            lG = rr.lG;
            niter = rr.iterations;
            method = "lossn.rec";
        } else if (chosen == "exact") {
            const lossn::LossnManjunathResult<T> mr =
                lossn::lossn_manjunath<T>(out.nu, out.A, out.Cvec, lossn::LossnManjunathOptions());
            QLen = mr.QLen;
            out.Loss = mr.Loss;
            lG = mr.lG;
            // The transform is direct, so the iteration count is 1 by
            // definition, as in the reference.
            niter = static_cast<int>(mr.iterations);
            method = "lossn.exact";
        } else if (chosen == "mci") {
            if constexpr (!num_traits<T>::has_transcendental) {
                throw UnsupportedError(
                    "solver_nc_lossn_analyzer: the 'mci' method forms its importance weights in "
                    "log space and reports a confidence interval, so it is unavailable under exact "
                    "arithmetic -- and meaningless there in any case, since the estimate is a "
                    "random variable. Use 'exact', the Manjunath-Sikdar transform, which is "
                    "rational throughout");
            } else {
                lossn::LossnMciOptions<T> mciopt;
                mciopt.samples = opt.samples;
                const lossn::LossnMciResult<T> mr =
                    lossn::lossn_mci<T>(out.nu, out.A, out.Cvec, mciopt,
                                        static_cast<std::uint64_t>(opt.seed));
                QLen = mr.QLen;
                out.Loss = mr.Loss;
                lG = mr.lG;
                // The sampler does not iterate; the reference reports the
                // realised sample count in the slot the other methods use for
                // the iteration count, so a caller can tell how much work
                // produced the estimate.
                niter = static_cast<int>(mr.nsamples);
                method = "lossn.mci";
            }
        } else {
            if constexpr (!num_traits<T>::has_transcendental) {
                throw UnsupportedError(
                    "solver_nc_lossn_analyzer: the 'erlangfp' method evaluates Erlang's loss "
                    "formula through logs and stops on a tolerance, so it is unavailable under "
                    "exact arithmetic and would not be exact there anyway. Use 'exact', the "
                    "Manjunath-Sikdar transform");
            } else {
                if (!integral)
                    throw UnsupportedError(
                        "solver_nc_lossn_analyzer: the 'erlangfp' method as ported takes integer "
                        "circuit requirements and an integer link capacity, and this region "
                        "declares fractional ones (a memory budget with fractional class sizes, or "
                        "an explicit linear constraint). The reference evaluates Erlang B through "
                        "factln and accepts them; the port would truncate them and solve a "
                        "different region. Use 'mci', which compares in real arithmetic");
                std::vector<int> Cint(J, 0);
                for (std::size_t j = 0; j < J; ++j)
                    Cint[j] = static_cast<int>(std::lround(num_traits<T>::to_double(out.Cvec[j])));
                // The reference hardwires these inside lossn_erlangfp rather
                // than reading options.iter_tol, and stops on a non-finite
                // increment, which the shared driver does only when asked.
                da::FpiOptions fpopt;
                fpopt.nanstop = true;
                const lossn::ErlangFpResult<T> fr =
                    lossn::lossn_erlangfp<T>(out.nu, out.A, Cint, fpopt);
                QLen = fr.QLen;
                out.Loss = fr.Loss;
                out.E = fr.E;
                niter = static_cast<int>(fr.iterations);
                method = "lossn.erlangfp";
            }
        }

        // 5. The metric table. QLen is the CARRIED load E[n_r], so the carried
        // throughput follows from Little's law at the infinite server.
        out.sol.sol.Q = Matrix<T>(M, K, zero);
        out.sol.sol.U = Matrix<T>(M, K, zero);
        out.sol.sol.R = Matrix<T>(M, K, zero);
        out.sol.sol.Tp = Matrix<T>(M, K, zero);
        out.sol.sol.C.assign(K, zero);
        out.sol.sol.X.assign(K, zero);
        for (std::size_t r = 0; r < K; ++r) {
            const std::size_t src = sn.classes[r].refstat;
            const T Xc = T(lambda[r] * (one - out.Loss[r]));
            out.sol.sol.X[r] = Xc;
            out.sol.sol.Tp(delayIdx - 1, r) = Xc;
            // The source emits the ACCEPTED (post-drop) rate, so that the
            // routing-based arrival rate `sn_get_arvr_from_tput` derives is
            // non-zero at the delay and agrees with the flow-conserving
            // departure throughput and with the rate SolverJMT simulates.
            if (src >= 1 && src <= M) out.sol.sol.Tp(src - 1, r) = Xc;
            out.sol.sol.Q(delayIdx - 1, r) = QLen[r];
            // An infinite server never queues, so the response time is the
            // service time and the "utilization" is the mean number of busy
            // servers, which is the population itself.
            if (mu[r] != zero) out.sol.sol.R(delayIdx - 1, r) = T(one / mu[r]);
            out.sol.sol.U(delayIdx - 1, r) = QLen[r];
        }
        out.sol.sol.lG = lG;
        out.sol.sol.iter = niter;
        out.sol.sol.method = method;
        out.sol.actualmethod = method;
        return out;
    }
}

}  // namespace nc
}  // namespace line

#endif  // LINE_SOLVERS_NC_SOLVER_NC_LOSSN_H
