/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The remaining `@@SolverCTMC` accessors: `getGenerator` / `getInfGen`,
 * `getStateSpace` / `getStateSpaceAggr` and the `getTranProb*` family. The
 * symbolic getters live in solver_ctmc_symbolic.h.
 *
 * WHAT AN ACCESSOR HERE IS INDEXED BY. Every quantity below is read off the
 * chain `solver_ctmc_analyzer` SOLVED, which is the enumerated space after the
 * DROP regions have censored it and, on a reducible encoding, restricted to one
 * weakly connected component. The reference's `getGenerator` and
 * `getStateSpace` instead call `solver_ctmc` and `State.spaceGenerator`
 * directly and so report the pre-restriction space, which on such a model has
 * more rows than its own `pi` vector does. The states that differ carry zero
 * stationary mass and no other accessor in this port can index them, so they
 * are dropped here rather than handed back as rows nothing else accepts.
 *
 * WHY THE FILTRATION IS AN OUTPUT AND NOT A DERIVATION. `eventFilt` is not
 * recoverable from Q: assembling the generator adds every synchronization's
 * contribution into one entry, and no inspection of the sum says which event
 * put what there. It has to be recorded while Q is built, which is what
 * `CtmcOptions::keep_filtration` turns on, so the accessors that return it
 * force that flag on rather than reporting an empty filtration.
 */
#ifndef LINE_SOLVERS_CTMC_SOLVER_CTMC_GETTERS_H
#define LINE_SOLVERS_CTMC_SOLVER_CTMC_GETTERS_H

#include <cmath>
#include <cstddef>
#include <limits>
#include <set>
#include <string>
#include <vector>

#include "line/api/mc/ctmc_passage.h"
#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/state.h"
#include "line/lang/qn/state_events.h"
#include "line/solvers/ctmc/solver_ctmc.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_prob.h"
#include "line/solvers/ctmc/solver_ctmc_transient.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

#include "line/api/sim/sim_runlength.h"

namespace line {
namespace ctmc {

/** `[infGen, eventFilt, ev]` of `@@SolverCTMC/getGenerator.m`. */
template <class T>
struct CtmcGenerator {
    Matrix<T> Q;                     ///< the infinitesimal generator
    std::vector<NetState<T>> space;  ///< row i of Q is space[i]
    /**
     * `eventFilt`: `filt[a]` holds only what synchronization `sync[a]`
     * contributed, so `sum_a filt[a]` is the off-diagonal part of Q. The two
     * vectors are indexed alike; that pairing is the whole content of the
     * filtration and is why `sync` travels with it.
     */
    std::vector<Matrix<T>> filt;
    std::vector<Sync<T>> sync;  ///< `ev`, the reference's `sn.sync`
    /**
     * The DERIVED START/PREEMPT filtrations, indexed [station-1][class-1]. They
     * are NOT part of `filt` and are NOT paired with `sync`: a START rides on an
     * arc `filt` already carries, so a caller summing `filt` must not see them.
     */
    std::vector<std::vector<Matrix<T>>> start_filt, preempt_filt;
};

/** `[stateSpace, localStateSpace]` of `@@SolverCTMC/getStateSpace.m`. */
template <class T>
struct CtmcStateSpace {
    /**
     * The enumerated states. `space[s].local[f]` is already the reference's
     * `localStateSpace{f}` row: this port keeps a state as its per-node blocks
     * and never flattens it, so the local decomposition costs nothing to
     * report and cannot disagree with the flat form.
     */
    std::vector<NetState<T>> space;
    Matrix<T> flat;                       ///< the blocks concatenated, as MATLAB returns them
    std::vector<std::size_t> node_width;  ///< column width of each stateful node's block
    /**
     * `localStateSpace{f}`: one matrix per stateful node, its DISTINCT local
     * rows in first-appearance order.
     *
     * WITHOUT THIS THE SECOND RETURN VALUE CANNOT BE FORMED. `getStateSpace.m`
     * under `lang='cpp'` returns `{stateSpace}` -- a single cell holding the
     * whole flat space -- because the local decomposition was not on the wire,
     * and a caller that indexes `localStateSpace{f}` per node then reads the
     * global space for every node. `node_width` alone is not enough: it says
     * where a block ENDS, not which rows a node admits.
     *
     * IT IS THE ENUMERATED CHAIN'S DECOMPOSITION, not `spaceGenerator`'s raw
     * per-node enumeration. The reference's `qnc.space{f}` is built before the
     * cartesian product and can hold a local row no global state uses; every
     * row here appears in `space`. On a model whose lattice admits every
     * combination the two coincide, and where they differ this one is the
     * decomposition of the chain that was actually solved.
     */
    std::vector<Matrix<T>> local;
};

/** The time-dependent answer of one `getTranProb*` query. */
template <class T>
struct CtmcTranProb {
    /**
     * The reference returns `Pi_t = [t, pi_t]`, one matrix with time glued on
     * as column 1. The two are kept apart here because they are not the same
     * quantity and concatenating them forces every consumer to know that
     * column 1 is not a probability.
     */
    std::vector<T> t;
    Matrix<T> pit;     ///< (ntimes x nstates) occupancy over the solved chain
    Matrix<T> labels;  ///< (nstates x width) the state descriptor the query asked for
};

/**
 * `SolverCTMC.getStartRate` and `getPreemptRate`: the DERIVED rates the
 * START/PREEMPT filtration reduces to.
 *
 * StartN(i,r) is how often per unit time a class-r service STARTS at station i,
 * and PreemptN(i,r) how often a class-r job in service is pushed back into the
 * buffer there. Both are computed by `solver_ctmc_avg_from_pi` already; what was
 * missing was any way to ask for them, because `solver_ctmc_run_analyzer` returns an
 * `mva::AvgResult` and that shape has no slot for a derived rate. The reference
 * exposes them as plain getters, and so do these.
 *
 * They need the filtration, so the solve is repeated with `keep_filtration` set
 * rather than read off a result that may not carry it: a caller who already has
 * one passes the `CtmcSolution` overload instead and pays nothing.
 *
 * At a lossless station with no in-service abandonment
 * StartN == TN + PreemptN, which is the identity to check them against.
 */
template <class T>
Matrix<T> ctmc_get_start_rate(const NetworkStruct<T>& sn, const CtmcOptions& opt) {
    CtmcOptions o = opt;
    o.keep_filtration = true;
    return solver_ctmc_analyzer(sn, o).avg.StartN;
}

/** As above, for a caller whose solution already carries the filtration. */
template <class T>
Matrix<T> ctmc_get_start_rate(const NetworkStruct<T>&, const CtmcSolution<T>& d) {
    return d.avg.StartN;
}

/** `SolverCTMC.getPreemptRate`; see {@link ctmc_get_start_rate}. */
template <class T>
Matrix<T> ctmc_get_preempt_rate(const NetworkStruct<T>& sn, const CtmcOptions& opt) {
    CtmcOptions o = opt;
    o.keep_filtration = true;
    return solver_ctmc_analyzer(sn, o).avg.PreemptN;
}

/** As above, for a caller whose solution already carries the filtration. */
template <class T>
Matrix<T> ctmc_get_preempt_rate(const NetworkStruct<T>&, const CtmcSolution<T>& d) {
    return d.avg.PreemptN;
}

namespace getters_detail {

/** Per-node block widths, which every state shares by construction. */
template <class T>
std::vector<std::size_t> node_widths(const std::vector<NetState<T>>& space) {
    std::vector<std::size_t> w;
    if (space.empty()) return w;
    for (std::size_t f = 0; f < space[0].local.size(); ++f) w.push_back(space[0].local[f].size());
    return w;
}

/** Concatenate the per-node blocks of every state into one matrix. */
template <class T>
Matrix<T> flatten(const std::vector<NetState<T>>& space) {
    const std::vector<std::size_t> w = node_widths(space);
    std::size_t total = 0;
    for (std::size_t f = 0; f < w.size(); ++f) total += w[f];
    Matrix<T> out(space.size(), total, num_traits<T>::from_int(0));
    for (std::size_t s = 0; s < space.size(); ++s) {
        std::size_t c = 0;
        for (std::size_t f = 0; f < space[s].local.size(); ++f)
            for (std::size_t j = 0; j < space[s].local[f].size(); ++j) out(s, c++) = space[s].local[f][j];
    }
    return out;
}

/**
 * The DISTINCT local rows of each stateful node, in first-appearance order.
 *
 * First-appearance and not sorted, because the reference's `qnc.space{f}` is
 * built in enumeration order and a caller that pairs a local index with a row
 * of `space` must see the same order on both sides.
 */
template <class T>
std::vector<Matrix<T>> local_spaces(const std::vector<NetState<T>>& space) {
    std::vector<Matrix<T>> out;
    if (space.empty()) return out;
    const std::size_t NF = space[0].local.size();
    for (std::size_t f = 0; f < NF; ++f) {
        std::vector<std::vector<T>> rows;
        std::set<std::vector<double>> seen;
        for (std::size_t s = 0; s < space.size(); ++s) {
            if (space[s].local.size() <= f) continue;
            const std::vector<T>& r = space[s].local[f];
            std::vector<double> key(r.size());
            for (std::size_t j = 0; j < r.size(); ++j) key[j] = num_traits<T>::to_double(r[j]);
            if (!seen.insert(key).second) continue;
            rows.push_back(r);
        }
        const std::size_t w = rows.empty() ? 0 : rows[0].size();
        Matrix<T> m(rows.size(), w, num_traits<T>::from_int(0));
        for (std::size_t i = 0; i < rows.size(); ++i)
            for (std::size_t j = 0; j < w && j < rows[i].size(); ++j) m(i, j) = rows[i][j];
        out.push_back(m);
    }
    return out;
}

/**
 * The stateful index of a node, refusing a node that carries no state.
 *
 * A stateless node -- a Router, a ClassSwitch -- has no block in any state, so
 * there is nothing to slice out for it and the query is a caller error rather
 * than an empty answer.
 */
template <class T>
std::size_t stateful_or_throw(const NetworkStruct<T>& sn, std::size_t ind, const char* what) {
    if (ind < 1 || ind > sn.nodes.size())
        throw InputError(std::string(what) + ": the node index is out of range");
    const std::size_t isf = sn.stateful_index(ind);
    if (isf == 0)
        throw InputError(std::string(what) +
                         ": node " + std::to_string(ind) +
                         " is stateless and holds no block of the network state");
    return isf;
}

/**
 * The timespan gate of the `getTranProb*` family.
 *
 * The reference refuses an infinite horizon by name because there is nothing to
 * integrate to: `pi(t)` on [0, Inf) is the stationary vector, which is what
 * `getProb*` already answers.
 */
template <class T>
void check_timespan(const T& t0, const T& t1, const char* what) {
    const double a = num_traits<T>::to_double(t0), b = num_traits<T>::to_double(t1);
    if (!std::isfinite(a) || !std::isfinite(b) || !(b > a))
        throw InputError(std::string(what) +
                         " requires a finite timespan [t0, t1] with t1 > t0; for the limit as t "
                         "grows use the stationary family (getProb, getProbAggr, getProbSys, "
                         "getProbSysAggr)");
}

}  // namespace getters_detail

/**
 * Port of `@@SolverCTMC/getGenerator.m`: the generator, its event filtration and
 * the synchronization list the filtration is indexed by.
 *
 * The synchronization list is rebuilt rather than carried, `refresh_sync` being
 * a pure function of the struct: it returns the same list, in the same order,
 * that indexed `filt` when the generator was assembled.
 */
template <class T>
CtmcGenerator<T> ctmc_get_generator(const NetworkStruct<T>& sn, const CtmcSolution<T>& d) {
    if (d.chain.filt.empty())
        throw InputError(
            "ctmc_get_generator: the solution carries no event filtration, which cannot be "
            "recovered from Q because its entries have already summed every synchronization's "
            "contribution; re-solve with CtmcOptions::keep_filtration = true");
    CtmcGenerator<T> g;
    g.Q = d.chain.Q;
    g.space = d.chain.space;
    g.filt = d.chain.filt;
    g.start_filt = d.chain.start_filt;
    g.preempt_filt = d.chain.preempt_filt;
    g.sync = refresh_sync(sn);
    return g;
}

/**
 * As above, solving the chain first.
 *
 * `keep_filtration` is forced on because the filtration is half the answer and
 * cannot be reconstructed afterwards; a caller who only wants Q should read
 * `CtmcSolution::chain` instead and not pay for one n x n matrix per event.
 */
template <class T>
CtmcGenerator<T> ctmc_get_generator(const NetworkStruct<T>& sn, const CtmcOptions& opt) {
    CtmcOptions o = opt;
    o.keep_filtration = true;
    return ctmc_get_generator(sn, solver_ctmc_analyzer(sn, o));
}

/** `@@SolverCTMC/getInfGen.m`, a pure alias of `getGenerator` in the reference. */
template <class T>
CtmcGenerator<T> ctmc_get_infgen(const NetworkStruct<T>& sn, const CtmcSolution<T>& d) {
    return ctmc_get_generator(sn, d);
}

/** `@@SolverCTMC/getInfGen.m`, solving the chain first. */
template <class T>
CtmcGenerator<T> ctmc_get_infgen(const NetworkStruct<T>& sn, const CtmcOptions& opt) {
    return ctmc_get_generator(sn, opt);
}

/**
 * Port of `@@SolverCTMC/getStateSpace.m`.
 *
 * The reference derives `localStateSpace` by cutting the flat matrix at the
 * width of each node's own space, a split that can only be got right by
 * carrying those widths alongside. Here the state IS the split, so `flat` is
 * the derived form and the widths are reported for a caller comparing columns
 * against MATLAB.
 */
template <class T>
CtmcStateSpace<T> ctmc_get_state_space(const NetworkStruct<T>&, const CtmcSolution<T>& d) {
    CtmcStateSpace<T> s;
    s.space = d.chain.space;
    s.flat = getters_detail::flatten(d.chain.space);
    s.node_width = getters_detail::node_widths(d.chain.space);
    s.local = getters_detail::local_spaces(d.chain.space);
    return s;
}

/** As above, solving the chain first. */
template <class T>
CtmcStateSpace<T> ctmc_get_state_space(const NetworkStruct<T>& sn, const CtmcOptions& opt) {
    return ctmc_get_state_space(sn, solver_ctmc_analyzer(sn, opt));
}

/** The answer of `@@SolverCTMC/getCdfFirstPassT.m`: the [F(t), t] curve with
 * its grid, density and resolved state sets. */
struct CtmcFirstPassage {
    std::vector<double> t;            ///< the grid, 1000 points to the horizon
    std::vector<double> F;            ///< CDF at t, clamped to [0, 1]
    std::vector<double> f;            ///< density at t
    std::vector<std::size_t> source;  ///< resolved 0-based rows; empty = conditional stationary
    std::vector<std::size_t> target;  ///< resolved 0-based rows
};

namespace getters_detail {

/**
 * A state set given as 1-based ROW INDICES into the state space or as matrices
 * of state rows, resolved to 0-based row indices against `flat`. An
 * unrecognised row is an error rather than a silent drop, since a passage into
 * a state that is not in the space is not a slow passage but an undefined one.
 */
template <class T>
std::vector<std::size_t> resolve_state_set(const Matrix<double>& S, const Matrix<T>& flat,
                                           const char* name) {
    std::set<std::size_t> idx;
    if (S.rows() == 0 || S.cols() == 0) return std::vector<std::size_t>();
    const std::size_t n = flat.rows();
    bool is_index_vector = (S.rows() == 1 || S.cols() == 1);
    if (is_index_vector) {
        for (std::size_t i = 0; i < S.rows() && is_index_vector; ++i)
            for (std::size_t j = 0; j < S.cols(); ++j) {
                const double v = S(i, j);
                if (v != std::floor(v) || v < 1 || v > static_cast<double>(n)) {
                    is_index_vector = false;
                    break;
                }
            }
    }
    if (is_index_vector) {
        for (std::size_t i = 0; i < S.rows(); ++i)
            for (std::size_t j = 0; j < S.cols(); ++j)
                idx.insert(static_cast<std::size_t>(S(i, j)) - 1);
    } else {
        if (S.cols() != flat.cols())
            throw InputError(std::string("getCdfFirstPassT: a state row in set ") + name +
                             " has " + std::to_string(S.cols()) + " columns where the state "
                             "space has " + std::to_string(flat.cols()));
        for (std::size_t i = 0; i < S.rows(); ++i) {
            bool found = false;
            for (std::size_t r = 0; r < n && !found; ++r) {
                bool eq = true;
                for (std::size_t j = 0; j < S.cols() && eq; ++j)
                    if (num_traits<T>::to_double(flat(r, j)) != S(i, j)) eq = false;
                if (eq) {
                    idx.insert(r);
                    found = true;
                }
            }
            if (!found)
                throw InputError(std::string("A state given in set ") + name +
                                 " is not in the state space.");
        }
    }
    return std::vector<std::size_t>(idx.begin(), idx.end());
}

}  // namespace getters_detail

/**
 * Port of `@@SolverCTMC/getCdfFirstPassT.m`: the distribution of the FIRST
 * PASSAGE TIME from state set A into state set B, on the CTMC underlying the
 * model. An empty A starts from the conditional stationary law on the
 * complement of B.
 *
 * THIS IS NOT getCdfRespT. That getter times a tagged job between an arrival
 * at a station and its departure, through the event filtration; this one times
 * the chain between two sets of states the caller names, and answers questions
 * the filtration cannot express -- the writer cycle time of a readers-writers
 * model, the time to fill a buffer, the time to leave a degraded region.
 *
 * Reference: P. G. Harrison and W. J. Knottenbelt, "Passage Time Distributions
 * in Large Markov Chains", 2002.
 */
template <class T>
CtmcFirstPassage ctmc_cdf_firstpasst(const NetworkStruct<T>&, const CtmcSolution<T>& d,
                                     const Matrix<double>& A, const Matrix<double>& B,
                                     const std::string& method = "expm") {
    static_assert(num_traits<T>::has_transcendental,
                  "getCdfFirstPassT takes a matrix exponential and therefore requires an "
                  "arithmetic with transcendental functions");
    const Matrix<T>& Q = d.chain.Q;
    const std::size_t n = Q.rows();
    const Matrix<T> flat = getters_detail::flatten(d.chain.space);

    CtmcFirstPassage out;
    out.target = getters_detail::resolve_state_set(B, flat, "B");
    if (out.target.empty())
        throw InputError(
            "The target state set B is empty: a first passage time into no state is undefined.");
    out.source = getters_detail::resolve_state_set(A, flat, "A");

    std::vector<T> pi0;
    if (!out.source.empty()) {
        pi0.assign(n, num_traits<T>::from_int(0));
        const T w = num_traits<T>::from_int(1) / num_traits<T>::from_int(
                        static_cast<long>(out.source.size()));
        for (std::size_t idx : out.source) pi0[idx] = w;
    }

    // The horizon is chosen the way the response-time getter chooses it: 100
    // events at the slowest rate in the chain.
    double min_rate = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            const double v = std::abs(num_traits<T>::to_double(Q(i, j)));
            if (v > GlobalConstants::FineTol && v < min_rate) min_rate = v;
        }
    const double thor = std::abs(100.0 / min_rate);
    std::vector<double> tset(1000);
    for (std::size_t i = 0; i < tset.size(); ++i)
        tset[i] = thor * static_cast<double>(i) / static_cast<double>(tset.size() - 1);

    const mc::PassageCurve<T> curve = mc::ctmc_passage_time(Q, pi0, out.target, tset, method);
    out.t = tset;
    out.F = curve.F;
    out.f = curve.f;
    return out;
}

/** As above, solving the chain first. */
template <class T>
CtmcFirstPassage ctmc_cdf_firstpasst(const NetworkStruct<T>& sn, const CtmcOptions& opt,
                                     const Matrix<double>& A, const Matrix<double>& B,
                                     const std::string& method = "expm") {
    return ctmc_cdf_firstpasst(sn, solver_ctmc_analyzer(sn, opt), A, B, method);
}

/** The answer of `@@SolverCTMC/getFirstPassTMoments.m`. */
template <class T>
struct CtmcFirstPassageMoments {
    std::vector<T> m;   ///< (nmax) moments for a passage started uniformly in A
    Matrix<T> mall;     ///< (nstates x nmax), one row per starting state
    std::vector<std::size_t> source;  ///< resolved 0-based rows; empty = conditional stationary
    std::vector<std::size_t> target;  ///< resolved 0-based rows
};

/**
 * Port of `@@SolverCTMC/getFirstPassTMoments.m`: moments of order 1..nmax of
 * the first passage time from state set A into state set B.
 *
 * NO TRANSFORM INVERSION AND NO TIME GRID ARE INVOLVED. The moments come from
 * Eq. 3 of Harrison and Knottenbelt (2002) -- one linear solve per order -- so
 * they are EXACT and are not limited by the horizon a CDF would have to be
 * truncated at. That is the whole reason this getter exists beside
 * `ctmc_cdf_firstpasst`: the variance or the skewness of a passage time costs
 * nmax solves here and a numerical integration of a truncated curve there.
 *
 * `mall` is zero on B and infinite where B cannot be reached, as the reference
 * reports it. A and B name states as in `ctmc_cdf_firstpasst`.
 */
template <class T>
CtmcFirstPassageMoments<T> ctmc_firstpasst_moments(const NetworkStruct<T>&,
                                                   const CtmcSolution<T>& d,
                                                   const Matrix<double>& A,
                                                   const Matrix<double>& B,
                                                   std::size_t nmax = 3) {
    static_assert(num_traits<T>::has_transcendental,
                  "getFirstPassTMoments marks an unreachable target with an infinity and "
                  "therefore requires an arithmetic that has one");
    const Matrix<T>& Q = d.chain.Q;
    const std::size_t n = Q.rows();
    const Matrix<T> flat = getters_detail::flatten(d.chain.space);

    CtmcFirstPassageMoments<T> out;
    out.target = getters_detail::resolve_state_set(B, flat, "B");
    if (out.target.empty())
        throw InputError(
            "The target state set B is empty: a first passage time into no state is undefined.");
    out.source = getters_detail::resolve_state_set(A, flat, "A");

    std::vector<T> pi0;
    if (!out.source.empty()) {
        pi0.assign(n, num_traits<T>::from_int(0));
        const T w = num_traits<T>::from_int(1) / num_traits<T>::from_int(
                        static_cast<long>(out.source.size()));
        for (std::size_t idx : out.source) pi0[idx] = w;
    }

    const mc::PassageMoments<T> pm = mc::ctmc_passage_moments(Q, pi0, out.target, nmax);
    out.m = pm.m;
    out.mall = pm.mall;
    return out;
}

/** As above, solving the chain first. */
template <class T>
CtmcFirstPassageMoments<T> ctmc_firstpasst_moments(const NetworkStruct<T>& sn,
                                                   const CtmcOptions& opt,
                                                   const Matrix<double>& A,
                                                   const Matrix<double>& B,
                                                   std::size_t nmax = 3) {
    return ctmc_firstpasst_moments(sn, solver_ctmc_analyzer(sn, opt), A, B, nmax);
}

/**
 * Port of `@@SolverCTMC/getStateSpaceAggr.m`: the per-(station, class) job
 * counts of every state, in column block order `(ist-1)*K + k`.
 *
 * The reference returns `[]` with a warning when the model has not been solved,
 * since its copy is a by-product cached by a previous run. There is no such
 * cache here: the aggregate is a function of the state space alone and is
 * recomputed, so the accessor either answers or throws.
 */
template <class T>
Matrix<T> ctmc_get_state_space_aggr(const NetworkStruct<T>& sn, const CtmcOptions& opt) {
    return ctmc_state_space_aggr(sn, solver_ctmc_analyzer(sn, opt).chain.space);
}

/** As above, for a caller who has already solved the chain. */
template <class T>
Matrix<T> ctmc_get_state_space_aggr(const NetworkStruct<T>& sn, const CtmcSolution<T>& d) {
    return ctmc_state_space_aggr(sn, d.chain.space);
}

namespace getters_detail {

/**
 * The per-class marginal of one node's block, for every state.
 *
 * `prob_detail::marginal_of` is the single decoder of a local row into class
 * counts -- it knows the phase blocks of a station and the leading count
 * columns of a Cache -- and is reused rather than re-derived, so an aggregated
 * transient query and an aggregated stationary one cannot disagree about what
 * a state holds.
 *
 * A SOURCE ROW IS ZERO, not Inf. `to_marginal` reports an infinite reservoir
 * for an EXT station, which describes the encoding rather than a queue length;
 * `ctmc_state_space_aggr` zeroes it for the same reason, and a label matrix
 * that disagreed with it would make the system and per-node aggregates
 * inconsistent.
 */
template <class T>
Matrix<T> node_marginal_labels(const NetworkStruct<T>& sn, std::size_t ind, std::size_t isf,
                               const std::vector<NetState<T>>& space) {
    const std::size_t K = sn.nclasses;
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> out(space.size(), K, zero);
    const std::size_t ist = sn.nodes[ind - 1].station;
    if (ist != 0 && sn.stations[ist - 1].nodetype == NodeType::Source) return out;
    for (std::size_t s = 0; s < space.size(); ++s) {
        const std::vector<T> m = prob_detail::marginal_of(sn, ind, space[s].local[isf - 1]);
        for (std::size_t k = 0; k < K && k < m.size(); ++k) out(s, k) = m[k];
    }
    return out;
}

/**
 * Refuse a (struct, trajectory) pair that cannot belong to one model.
 *
 * The overloads taking a ready `CtmcTransient` cannot verify in general that it
 * was integrated from the struct they are handed, and they should not try: the
 * point of those overloads is that four queries share one integration. This
 * catches the one mismatch that would read past the end of a state instead of
 * merely answering about the wrong model.
 */
template <class T>
void check_pair(std::size_t isf, const std::vector<NetState<T>>& space, const char* what) {
    if (!space.empty() && isf > space[0].local.size())
        throw InputError(std::string(what) +
                         ": the transient solution has fewer stateful nodes than the model it "
                         "was queried with, so it was not integrated from that model");
}

/** The block of one node, for every state: the reference's `SSnode`. */
template <class T>
Matrix<T> node_labels(std::size_t isf, const std::vector<NetState<T>>& space) {
    const std::size_t w = space.empty() ? 0 : space[0].local[isf - 1].size();
    Matrix<T> out(space.size(), w, num_traits<T>::from_int(0));
    for (std::size_t s = 0; s < space.size(); ++s)
        for (std::size_t j = 0; j < w; ++j) out(s, j) = space[s].local[isf - 1][j];
    return out;
}

/** The time grid and occupancy every `getTranProb*` query shares. */
template <class T>
CtmcTranProb<T> tran_common(const CtmcTransient<T>& tr) {
    CtmcTranProb<T> p;
    p.t = tr.t;
    p.pit = tr.pit;
    return p;
}

}  // namespace getters_detail

/**
 * Port of `@@SolverCTMC/getTranProb.m`: pi(t) over the whole chain, labelled by
 * one node's local state.
 *
 * IT IS NOT A PER-STATE TRANSIENT PROBABILITY, despite the name's symmetry with
 * `getProb`. The reference returns the FULL occupancy vector together with the
 * node's slice of the state space, leaving the caller to sum the rows sharing
 * the local state it cares about; that is a strictly richer answer than one
 * marginal and is reproduced as such.
 *
 * @param ind 1-based node index
 * @param sn the refreshed network struct
 * @param tr transient solution whose pi(t) is being labelled
 */
template <class T>
CtmcTranProb<T> ctmc_get_tran_prob(const NetworkStruct<T>& sn, const CtmcTransient<T>& tr,
                                   std::size_t ind) {
    assert_phase_type_states(sn, "getTranProb");
    const std::size_t isf = getters_detail::stateful_or_throw(sn, ind, "getTranProb");
    getters_detail::check_pair(isf, tr.chain.chain.space, "getTranProb");
    CtmcTranProb<T> p = getters_detail::tran_common(tr);
    p.labels = getters_detail::node_labels(isf, tr.chain.chain.space);
    return p;
}

/** As above, integrating the forward equation first. */
template <class T>
CtmcTranProb<T> ctmc_get_tran_prob(const NetworkStruct<T>& sn, const CtmcOptions& opt,
                                   std::size_t ind, const T& t0, const T& t1) {
    static_assert(num_traits<T>::has_transcendental,
                  "getTranProb integrates the forward equation, whose adaptive step controller "
                  "is transcendental; use --arith double or real");
    assert_phase_type_states(sn, "getTranProb");
    getters_detail::check_timespan(t0, t1, "getTranProb");
    return ctmc_get_tran_prob(sn, solver_ctmc_transient_analyzer(sn, opt, t0, t1), ind);
}

/**
 * Port of `@@SolverCTMC/getTranProbAggr.m`: pi(t), labelled by one node's
 * per-class job counts.
 *
 * THE REFERENCE SLICES THE AGGREGATE BY NODE INDEX, `SSa(:, (jnd-1)*K+1 :
 * jnd*K)`, while `ctmc_ssg` writes that matrix in STATION blocks
 * `(ist-1)*K+1 : ist*K`. The two indices coincide only when every node is a
 * station, so on a model carrying a Router or a ClassSwitch the reference reads
 * the wrong block. The labels are decoded from the node's own state here
 * instead of sliced out of a station-indexed matrix, which sidesteps the
 * mismatch and extends to a stateful non-station -- a Cache -- that no station
 * block describes at all.
 *
 * @param ind 1-based node index
 * @param sn the refreshed network struct
 * @param tr transient solution whose pi(t) is being labelled
 */
template <class T>
CtmcTranProb<T> ctmc_get_tran_prob_aggr(const NetworkStruct<T>& sn, const CtmcTransient<T>& tr,
                                        std::size_t ind) {
    assert_phase_type_states(sn, "getTranProbAggr");
    const std::size_t isf = getters_detail::stateful_or_throw(sn, ind, "getTranProbAggr");
    getters_detail::check_pair(isf, tr.chain.chain.space, "getTranProbAggr");
    CtmcTranProb<T> p = getters_detail::tran_common(tr);
    p.labels = getters_detail::node_marginal_labels(sn, ind, isf, tr.chain.chain.space);
    return p;
}

/** As above, integrating the forward equation first. */
template <class T>
CtmcTranProb<T> ctmc_get_tran_prob_aggr(const NetworkStruct<T>& sn, const CtmcOptions& opt,
                                        std::size_t ind, const T& t0, const T& t1) {
    static_assert(num_traits<T>::has_transcendental,
                  "getTranProbAggr integrates the forward equation, whose adaptive step "
                  "controller is transcendental; use --arith double or real");
    assert_phase_type_states(sn, "getTranProbAggr");
    getters_detail::check_timespan(t0, t1, "getTranProbAggr");
    return ctmc_get_tran_prob_aggr(sn, solver_ctmc_transient_analyzer(sn, opt, t0, t1), ind);
}

/**
 * Port of `@@SolverCTMC/getTranProbSys.m`: pi(t), labelled by the whole network
 * state with its phases.
 */
template <class T>
CtmcTranProb<T> ctmc_get_tran_prob_sys(const NetworkStruct<T>& sn, const CtmcTransient<T>& tr) {
    assert_phase_type_states(sn, "getTranProbSys");
    CtmcTranProb<T> p = getters_detail::tran_common(tr);
    p.labels = getters_detail::flatten(tr.chain.chain.space);
    return p;
}

/** As above, integrating the forward equation first. */
template <class T>
CtmcTranProb<T> ctmc_get_tran_prob_sys(const NetworkStruct<T>& sn, const CtmcOptions& opt,
                                       const T& t0, const T& t1) {
    static_assert(num_traits<T>::has_transcendental,
                  "getTranProbSys integrates the forward equation, whose adaptive step "
                  "controller is transcendental; use --arith double or real");
    assert_phase_type_states(sn, "getTranProbSys");
    getters_detail::check_timespan(t0, t1, "getTranProbSys");
    return ctmc_get_tran_prob_sys(sn, solver_ctmc_transient_analyzer(sn, opt, t0, t1));
}

/**
 * Port of `@@SolverCTMC/getTranProbSysAggr.m`: pi(t), labelled by the network's
 * per-(station, class) job counts.
 *
 * The labels are `ctmc_state_space_aggr`, the same matrix the transient
 * analyzer already integrates Q(t) and U(t) against, so a caller summing these
 * rows by hand reproduces its `QNt` exactly.
 */
template <class T>
CtmcTranProb<T> ctmc_get_tran_prob_sys_aggr(const NetworkStruct<T>& sn,
                                            const CtmcTransient<T>& tr) {
    assert_phase_type_states(sn, "getTranProbSysAggr");
    CtmcTranProb<T> p = getters_detail::tran_common(tr);
    p.labels = ctmc_state_space_aggr(sn, tr.chain.chain.space);
    return p;
}

/** As above, integrating the forward equation first. */
template <class T>
CtmcTranProb<T> ctmc_get_tran_prob_sys_aggr(const NetworkStruct<T>& sn, const CtmcOptions& opt,
                                            const T& t0, const T& t1) {
    static_assert(num_traits<T>::has_transcendental,
                  "getTranProbSysAggr integrates the forward equation, whose adaptive step "
                  "controller is transcendental; use --arith double or real");
    assert_phase_type_states(sn, "getTranProbSysAggr");
    getters_detail::check_timespan(t0, t1, "getTranProbSysAggr");
    return ctmc_get_tran_prob_sys_aggr(sn, solver_ctmc_transient_analyzer(sn, opt, t0, t1));
}

/*
 * `@@SolverCTMC/getSymbolicGenerator.m` and `getSymbolicSolution.m` are NOT
 * here: they live in solver_ctmc_symbolic.h as `ctmc_symbolic_generator` and
 * `ctmc_symbolic_solution`.
 *
 * They were refused by name until 2026-07-31, on the grounds that this port had
 * no computer algebra. Half of that was never true -- the generator is LINEAR in
 * the event symbols, so it is one numeric filtration per event and needs no
 * algebra to assemble -- and the other half stopped being true when `api/sym`
 * gained a client for the line-sage-rest service the reference itself uses.
 * They are kept in their own header so that this one, which every CTMC accessor
 * includes, does not drag in a socket-using HTTP client.
 */

/*
 * `@@SolverCTMC/getSjrnT.m` and `sjrnT.m` are aliases of `getCdfRespT` in the
 * reference and are NOT defined here: an alias belongs beside the entry point
 * it forwards to, and the CTMC response-time CDF lives in its own header (see
 * `solver_mam_get_sjrn_t` in solver_mam_runner.h and `solver_nc_sjrnt` in
 * solver_nc_cdf.h for the two established precedents). Defining it here would
 * make this header depend on that one for nothing but a forwarding call.
 */

/**
 * Port of `@@SolverCTMC/getAsymptoticVariance.m`: the asymptotic variance of the
 * time-average of a reward along a sample path of this model's CTMC.
 *
 * WHAT IT IS FOR. A simulation estimate of a steady-state mean has a standard
 * error that shrinks like sqrt(sigma^2/t), where sigma^2 is NOT the stationary
 * variance of the reward but its ASYMPTOTIC variance, which also carries the
 * autocorrelation of the path. That number is what says how long a run has to
 * be, and `sim_runlength` turns it into a run length for a target precision. It
 * cannot be guessed from the stationary variance: on M/M/1 the two differ by a
 * factor that blows up like (1-rho)^-2.
 *
 * The reward is a function of the state ROW, evaluated on the state space the
 * generator was built from.
 */
template <class T>
sim::AsymVarResult<T> ctmc_get_asymptotic_variance(
    const NetworkStruct<T>& sn, const CtmcOptions& opt,
    const std::function<T(const NetState<T>&)>& reward) {
    const CtmcGenerator<T> g = ctmc_get_generator(sn, opt);
    std::vector<T> f;
    f.reserve(g.space.size());
    for (std::size_t i = 0; i < g.space.size(); ++i) f.push_back(reward(g.space[i]));
    return sim::sim_asymvar_ctmc<T>(g.Q, f);
}

}  // namespace ctmc
}  // namespace line

#endif  // LINE_SOLVERS_CTMC_SOLVER_CTMC_GETTERS_H
