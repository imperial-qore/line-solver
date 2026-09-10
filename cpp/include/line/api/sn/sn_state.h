/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_STATE_H
#define LINE_API_SN_SN_STATE_H

/**
 * Ports of matlab/src/api/sn/sn_get_state_aggr.m and sn_is_state_valid.m.
 *
 * THE STATE ARGUMENT. The reference reads `sn.state`, a cell of one row per
 * STATEFUL node that the model layer writes at `initDefault` time. This struct
 * carries only the pieces of a declared state that cannot be derived (a Place's
 * initial marking, a warm cache's contents, a declared prior over a declared
 * space), and the solvers build the rest. Both functions therefore take the
 * state explicitly, which is also what makes them usable on a candidate state
 * rather than only on the one the struct happens to hold.
 *
 * WHAT VALIDITY MEANS. Four things, in the reference's order: no class is
 * present at a station that does not serve it; no station holds more jobs of a
 * class than its capacity allows; the number of jobs IN SERVICE does not exceed
 * the server count at a discipline that cannot overlap service, and never
 * exceeds the number present; and every closed chain holds exactly its
 * population. The last is the one that catches a hand-written initial state,
 * and it is a relative test at CoarseTol because a fractional population is
 * legal.
 *
 * ARITHMETIC: field. Counting and comparison.
 */

#include <cmath>
#include <cstddef>
#include <map>
#include <utility>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/state.h"
#include "line/lang/qn/state_events.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace api {

/**
 * Port of sn_get_state_aggr: the per-class job counts of each stateful node's
 * state row, with the phase and buffer encoding aggregated away.
 *
 * @param state state[isf] is the row of the (isf+1)-th stateful node
 * @return out[isf] is that node's (ni, nir) pair
 */
template <class T>
std::vector<std::pair<T, std::vector<T>>> sn_get_state_aggr(
    const qn::NetworkStruct<T>& sn, const std::vector<std::vector<T>>& state) {
    std::vector<std::pair<T, std::vector<T>>> out;
    out.reserve(state.size());
    for (std::size_t isf = 0; isf < state.size(); ++isf) {
        const std::size_t ind = sn.stateful_nodes[isf];
        out.push_back(qn::to_marginal_aggr(sn, ind, state[isf]));
    }
    return out;
}

/**
 * The (nstations x nclasses) per-class job counts of the model's OWN state.
 *
 * This is what the reference reads as `State.toMarginal(sn, ist, state{isf})`
 * inside every `getProb*` getter, and it is a DECLARED quantity rather than a
 * derived one: `setState` and `initFromMarginal` put the jobs where the caller
 * asked, and the probability of "the model's state" is the probability of that
 * placement. The writers emit the row as a one-row `stateSpace` with a `[1]`
 * `statePrior` and the reader stores the pair, so it reaches here.
 *
 * A STATION THAT DECLARES NOTHING FALLS BACK TO THE DEFAULT MARKING -- every
 * closed class's jobs at its reference station -- which is what `initDefault`
 * would have put there, so a model that was never initialized answers exactly as
 * it did before any of this existed. The fallback is per station and not
 * all-or-nothing, matching the reference, where `setState` on one node leaves
 * the others at whatever they held.
 */
template <class T>
Matrix<T> sn_declared_marginal(const qn::NetworkStruct<T>& sn) {
    const std::size_t M = sn.nstations, R = sn.nclasses;
    Matrix<T> nir(M, R, num_traits<T>::from_int(0));
    std::vector<bool> declared(M, false);
    for (std::size_t i = 0; i < M; ++i) {
        const std::size_t ind = sn.station_to_node[i];  // 1-based node index
        const typename std::map<std::size_t, Matrix<T>>::const_iterator ss =
            sn.statespace.find(ind);
        if (ss == sn.statespace.end() || ss->second.rows() != 1 || ss->second.cols() == 0) continue;
        std::vector<T> row(ss->second.cols());
        for (std::size_t c = 0; c < ss->second.cols(); ++c) row[c] = ss->second(0, c);
        // A row the encoding cannot decode is left to the default marking rather
        // than reported as zeros: zeros are a state, and the wrong one.
        std::pair<T, std::vector<T>> m;
        try {
            m = qn::to_marginal_aggr(sn, ind, row);
        } catch (const Error&) {
            continue;
        }
        if (m.second.size() != R) continue;
        for (std::size_t r = 0; r < R; ++r) nir(i, r) = m.second[r];
        declared[i] = true;
    }
    for (std::size_t r = 0; r < R; ++r) {
        const double pop = sn.classes[r].population;
        if (!std::isfinite(pop)) continue;  // an open class has no initial jobs
        const std::size_t rs = sn.classes[r].refstat;  // 1-based station
        if (rs >= 1 && rs <= M && !declared[rs - 1])
            nir(rs - 1, r) = num_traits<T>::from_double(pop);
    }
    return nir;
}

/**
 * Port of `State.isValid`, which is the whole body of sn_is_state_valid once
 * the marginals are formed.
 *
 * @param n (nstations x nclasses) jobs present
 * @param s (nstations x nclasses) jobs in service
 */
template <class T>
bool sn_state_counts_valid(const qn::NetworkStruct<T>& sn, const Matrix<T>& n, const Matrix<T>& s) {
    const std::size_t M = sn.nstations, R = sn.nclasses;
    if (n.rows() == 0 && s.rows() != 0) return false;
    for (std::size_t i = 0; i < M; ++i) {
        const std::size_t nd = sn.station_to_node[i];
        const bool is_place = nd != 0 && sn.nodes[nd - 1].nodetype == qn::NodeType::Place;
        for (std::size_t r = 0; r < R; ++r) {
            if (!is_place && !sn.disabled.empty() && sn.disabled[i][r] &&
                num_traits<T>::to_double(n(i, r)) > 0.0)
                return false;
            if (i < sn.classcap.size() && r < sn.classcap[i].size() &&
                num_traits<T>::to_double(n(i, r)) > sn.classcap[i][r])
                return false;
        }
    }
    if (s.rows() != 0) {
        for (std::size_t i = 0; i < M; ++i) {
            const std::size_t nd = sn.station_to_node[i];
            if (nd != 0 && sn.nodes[nd - 1].nodetype == qn::NodeType::Place) continue;
            if (!(sn.stations[i].nservers > 0.0)) continue;
            double sums = 0.0;
            for (std::size_t r = 0; r < R; ++r) sums += num_traits<T>::to_double(s(i, r));
            if (sums > sn.stations[i].nservers) {
                switch (sn.stations[i].sched) {
                    case qn::SchedStrategy::FCFS:
                    case qn::SchedStrategy::SIRO:
                    case qn::SchedStrategy::LCFS:
                    case qn::SchedStrategy::HOL:
                    case qn::SchedStrategy::POLLING: return false;
                    default: break;
                }
            }
            // the reference's `any(n<s)` is over the WHOLE pair of matrices,
            // not this station's row; reproduced, because a state that puts a
            // job in service where none is present is invalid wherever it is
            for (std::size_t a = 0; a < M; ++a)
                for (std::size_t r = 0; r < R; ++r)
                    if (num_traits<T>::to_double(n(a, r)) < num_traits<T>::to_double(s(a, r)))
                        return false;
        }
    }
    for (std::size_t c = 0; c < sn.nchains; ++c) {
        double njobs_chain = 0.0;
        bool open = false;
        for (std::size_t r = 0; r < R; ++r) {
            if (!sn.chains[c][r]) continue;
            if (std::isinf(sn.classes[r].population)) open = true;
            njobs_chain += sn.classes[r].population;
        }
        if (open || std::isinf(njobs_chain)) continue;
        double statejobs = 0.0;
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r)
                if (sn.chains[c][r]) statejobs += num_traits<T>::to_double(n(i, r));
        if (!(statejobs > 0.0)) return false;
        if (std::fabs(1.0 - njobs_chain / statejobs) > qn::GlobalConstants::CoarseTol)
            return false;
    }
    return true;
}

/**
 * Port of sn_is_state_valid: form the station marginals of `state` and test
 * them.
 *
 * A station whose row carries several candidate states is answered on its FIRST
 * row, which is what the reference does after warning: a validity question
 * about a set of states has no single answer.
 */
template <class T>
bool sn_is_state_valid(const qn::NetworkStruct<T>& sn, const std::vector<std::vector<T>>& state) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = sn.nstations, R = sn.nclasses;
    Matrix<T> n(M, R, zero), s(M, R, zero);
    for (std::size_t ist = 1; ist <= M; ++ist) {
        const std::size_t isf = sn.stateful_of_station(ist);
        if (isf == 0 || isf > state.size()) return false;
        const std::size_t ind = sn.stateful_nodes[isf - 1];
        const qn::RowLayout<T> L = qn::row_layout(sn, ind, state[isf - 1].size());
        std::vector<std::size_t> ph(R, 1), shift(R, 0);
        for (std::size_t r = 0; r < R; ++r) {
            ph[r] = L.K[r];
            shift[r] = L.Ks[r];
        }
        const qn::Marginal<T> m = qn::to_marginal(sn, ist, state[isf - 1], ph, shift, L.nvar);
        for (std::size_t r = 0; r < R; ++r) {
            n(ist - 1, r) = m.nir[r];
            s(ist - 1, r) = m.sir[r];
        }
    }
    return sn_state_counts_valid(sn, n, s);
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_STATE_H
