/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The SolverCTMC probability family: `solver_ctmc_joint`, `_jointaggr`,
 * `_marg`, `_margaggr`, and the gate `@@SolverCTMC/assertPhaseTypeStates` puts
 * in front of all four.
 *
 * THE DISTINCTION THE FOUR NAMES ENCODE, because it is easy to mix up:
 *
 *   joint      P(the whole network is in exactly this state), phases included
 *   jointaggr  P(the whole network holds exactly these per-class counts),
 *              summed over every phase and buffer arrangement that realizes them
 *   marg       per STATION, P(that station is in exactly this local state)
 *   margaggr   per STATION, P(that station holds exactly these per-class counts)
 *
 * `joint` and `marg` are per-state answers and `jointaggr` and `margaggr` are
 * aggregates of them, which is precisely why the ME gate below applies to all
 * four and not only to the first two: an aggregate over PHASES is a probability
 * under a matrix-exponential, but an aggregate over STATES sharing a marginal
 * is not, because the sum still runs over signed terms.
 */
#ifndef LINE_SOLVERS_CTMC_SOLVER_CTMC_PROB_H
#define LINE_SOLVERS_CTMC_SOLVER_CTMC_PROB_H

#include <cstddef>
#include <string>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/state.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace ctmc {

/**
 * Port of `@@SolverCTMC/assertPhaseTypeStates`: refuse a query whose answer
 * would be a per-state probability under a matrix-exponential process.
 *
 * A matrix-exponential embeds in the generator with POSITIVE off-diagonal
 * entries in D0 (equivalently a signed entry vector), so the stationary vector
 * is a SIGNED measure: only its aggregates over each phase block are
 * probabilities. Mean measures stay exact, being linear in that vector, but a
 * per-state or transient answer is not a probability at all, and uniformization
 * -- a Poisson mixture of powers of I + Q/lambda -- diverges on a signed
 * generator. Such queries are refused rather than answered with a number that
 * looks like a probability.
 *
 * THE TEST IS ON THE MATRICES, not on a flag. MATLAB carries `sn.isph`, which
 * this NetworkStruct has no counterpart for; the property it records is exactly
 * "D0 has no positive off-diagonal and D1 is non-negative", so that is what is
 * checked here rather than a field being invented to hold the answer.
 */
template <class T>
void assert_phase_type_states(const NetworkStruct<T>& sn, const std::string& what) {
    // Same matrix test the avg reduction gates its clamp on, so one predicate
    // decides both "refuse a per-state answer" and "keep the negative mass".
    if (!ctmc_all_phasetype(sn))
        throw UnsupportedError(
            what +
            " is unavailable: the model has a matrix-exponential (ME) service or arrival "
            "process, so the stationary vector of the generator is a signed measure and "
            "per-state probabilities and uniformization-based transients do not exist. "
            "Mean measures (getAvg, getAvgTable) remain exact");
}

namespace prob_detail {

/** Left-pad a local row with zeros to the width the enumerated space uses. */
template <class T>
std::vector<T> align(const std::vector<T>& row, std::size_t width) {
    if (row.size() >= width) return row;
    std::vector<T> out(width - row.size(), num_traits<T>::from_int(0));
    out.insert(out.end(), row.begin(), row.end());
    return out;
}

/** The per-class marginal of one stateful node's local row, or empty. */
template <class T>
std::vector<T> marginal_of(const NetworkStruct<T>& sn, std::size_t ind,
                           const std::vector<T>& row) {
    const std::size_t R = sn.nclasses;
    const std::size_t ist = sn.nodes[ind - 1].station;
    if (ist == 0) {
        // A stateful non-station -- a Cache -- carries one count column per
        // class ahead of its own local block, which is the marginal itself.
        return std::vector<T>(row.begin(), row.begin() + std::min(R, row.size()));
    }
    std::vector<std::size_t> ph(R, 1), shift(R, 0);
    std::size_t w = 0;
    for (std::size_t r = 0; r < R; ++r) {
        ph[r] = sn.phasessz_of(ist, r + 1);
        shift[r] = w;
        w += ph[r];
    }
    return qn::to_marginal(sn, ist, row, ph, shift, sn.nvars_of(ind)).nir;
}

}  // namespace prob_detail

/**
 * Port of `solver_ctmc_joint`: P(the network is in exactly `state`).
 *
 * Returns zero when the state is not in the enumerated space, which is the
 * honest answer for a state the encoding cannot represent or the dynamics
 * cannot reach -- not an error, because `findrows` returning nothing is how the
 * reference reports the same thing.
 */
template <class T>
T solver_ctmc_joint(const NetworkStruct<T>& sn, const CtmcSolution<T>& d,
                    const NetState<T>& state) {
    assert_phase_type_states(sn, "getProbSys");
    NetState<T> q = state;
    for (std::size_t f = 0; f < q.local.size() && f < d.chain.space[0].local.size(); ++f)
        q.local[f] = prob_detail::align(q.local[f], d.chain.space[0].local[f].size());
    const std::vector<double> key = ctmc_detail::state_key(q);
    for (std::size_t s = 0; s < d.chain.space.size(); ++s)
        if (ctmc_detail::state_key(d.chain.space[s]) == key) return d.pi[s];
    return num_traits<T>::from_int(0);
}

/**
 * Port of `solver_ctmc_jointaggr`: P(the network holds exactly these per-class
 * counts), summed over every phase and buffer arrangement that realizes them.
 */
template <class T>
T solver_ctmc_jointaggr(const NetworkStruct<T>& sn, const CtmcSolution<T>& d,
                        const NetState<T>& state) {
    assert_phase_type_states(sn, "getProbSysAggr");
    const std::vector<std::size_t>& sfn = sn.stateful_nodes;
    std::vector<std::vector<T>> want(sfn.size());
    for (std::size_t f = 0; f < sfn.size(); ++f)
        want[f] = prob_detail::marginal_of(sn, sfn[f], state.local[f]);

    T acc = num_traits<T>::from_int(0);
    for (std::size_t s = 0; s < d.chain.space.size(); ++s) {
        bool same = true;
        for (std::size_t f = 0; f < sfn.size() && same; ++f) {
            const std::vector<T> got =
                prob_detail::marginal_of(sn, sfn[f], d.chain.space[s].local[f]);
            if (got.size() != want[f].size()) { same = false; break; }
            for (std::size_t r = 0; r < got.size(); ++r)
                if (num_traits<T>::to_double(got[r]) !=
                    num_traits<T>::to_double(want[f][r])) { same = false; break; }
        }
        if (same) acc += d.pi[s];
    }
    return acc;
}

/**
 * Port of `solver_ctmc_marg`: per STATION, P(that station is in exactly its
 * local slice of `state`), marginalized over every other node.
 *
 * @return one entry per station, indexed 0-based
 */
template <class T>
std::vector<T> solver_ctmc_marg(const NetworkStruct<T>& sn, const CtmcSolution<T>& d,
                                const NetState<T>& state) {
    assert_phase_type_states(sn, "getProb");
    const T zero = num_traits<T>::from_int(0);
    std::vector<T> out(sn.nstations, zero);
    for (std::size_t ist = 1; ist <= sn.nstations; ++ist) {
        const std::size_t isf = sn.stateful_of_station(ist);
        if (isf == 0) continue;
        const std::vector<T> want =
            prob_detail::align(state.local[isf - 1], d.chain.space[0].local[isf - 1].size());
        for (std::size_t s = 0; s < d.chain.space.size(); ++s)
            if (d.chain.space[s].local[isf - 1] == want) out[ist - 1] += d.pi[s];
    }
    return out;
}

/**
 * Port of `solver_ctmc_margaggr`: per STATION, P(that station holds exactly
 * these per-class counts).
 *
 * This is the aggregate `solver_ctmc_marg` is the refinement of: it sums the
 * per-state probabilities of every local state sharing the marginal, so on a
 * single-phase model with no buffer ordering the two coincide.
 */
template <class T>
std::vector<T> solver_ctmc_margaggr(const NetworkStruct<T>& sn, const CtmcSolution<T>& d,
                                    const NetState<T>& state) {
    assert_phase_type_states(sn, "getProbAggr");
    const T zero = num_traits<T>::from_int(0);
    std::vector<T> out(sn.nstations, zero);
    for (std::size_t ist = 1; ist <= sn.nstations; ++ist) {
        const std::size_t isf = sn.stateful_of_station(ist);
        const std::size_t ind = sn.node_of_station(ist);
        if (isf == 0) continue;
        const std::vector<T> want = prob_detail::marginal_of(sn, ind, state.local[isf - 1]);
        for (std::size_t s = 0; s < d.chain.space.size(); ++s) {
            const std::vector<T> got =
                prob_detail::marginal_of(sn, ind, d.chain.space[s].local[isf - 1]);
            bool same = got.size() == want.size();
            for (std::size_t r = 0; r < got.size() && same; ++r)
                if (num_traits<T>::to_double(got[r]) != num_traits<T>::to_double(want[r]))
                    same = false;
            if (same) out[ist - 1] += d.pi[s];
        }
    }
    return out;
}

/**
 * Port of `solver_ctmc_ratecomplement`: the long-run rate of an action as seen
 * from each TANGIBLE state, given the action's rate filter `D`.
 *
 * Vanishing states are removed from the generator by stochastic complementation,
 * so an action that fires only in vanishing states -- a fork firing, a join
 * departure, the firing of an immediate SPN mode -- would be lost if its rate
 * were read off the tangible rows alone. The rate observed from tangible state
 * s is the direct exit rate via the action plus the expected number of firings
 * along the vanishing chain entered from s:
 *
 *   r = D(nonimm,:)*1 + Q12*(-Q22)^-1*(D(imm,:)*1)
 *
 * NOT YET REACHED BY THIS PORT'S GENERATOR, which gives an immediate transition
 * the reference's ~1e8 rate rather than eliminating it, so there are no
 * vanishing states to complement out. It is ported at its own signature so the
 * elimination can be added without re-deriving the correction.
 *
 * @param nonimm 0-based tangible row indices, in the order they appear in Q11
 * @param imm    0-based vanishing row indices
 * @param Q12    tangible-to-vanishing block
 * @param Q22    vanishing-to-vanishing block
 * @param D generator whose immediate states are being eliminated
 */
template <class T>
std::vector<T> solver_ctmc_ratecomplement(const Matrix<T>& D,
                                          const std::vector<std::size_t>& nonimm,
                                          const std::vector<std::size_t>& imm,
                                          const Matrix<T>& Q12, const Matrix<T>& Q22) {
    const T zero = num_traits<T>::from_int(0);
    std::vector<T> r(nonimm.size(), zero);
    for (std::size_t a = 0; a < nonimm.size(); ++a)
        for (std::size_t j = 0; j < D.cols(); ++j) r[a] += D(nonimm[a], j);
    if (imm.empty()) return r;

    std::vector<T> b(imm.size(), zero);
    for (std::size_t a = 0; a < imm.size(); ++a)
        for (std::size_t j = 0; j < D.cols(); ++j) b[a] += D(imm[a], j);
    Matrix<T> negQ22(Q22.rows(), Q22.cols());
    for (std::size_t a = 0; a < Q22.rows(); ++a)
        for (std::size_t c = 0; c < Q22.cols(); ++c) negQ22(a, c) = T(-Q22(a, c));
    const std::vector<T> x = solve(negQ22, b);
    for (std::size_t a = 0; a < nonimm.size() && a < Q12.rows(); ++a)
        for (std::size_t c = 0; c < Q12.cols() && c < x.size(); ++c) r[a] += T(Q12(a, c) * x[c]);
    return r;
}

}  // namespace ctmc
}  // namespace line

#endif  // LINE_SOLVERS_CTMC_SOLVER_CTMC_PROB_H
