/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_CTMC_SOLVER_CTMC_MAPC_H
#define LINE_SOLVERS_CTMC_SOLVER_CTMC_MAPC_H

/**
 * @file
 * @ingroup line_solvers
 * Port of `solver_ctmc_mapc.m`: the pair form of a MAP (and MMPP2) service at an FCFS
 * station, so the CTMC follows the law JMT and both LDES engines sample.
 *
 * Those engines keep ONE sampler per (station, class): every draw starts in the phase the
 * previous draw ENDED in, draws chained in service-start order, idle periods included. With
 * c > 1 servers the next start can happen while earlier draws are still running, so the
 * landing phase of a draw must be known when it starts. With V = (-D0)^-1 D1 a draw started
 * in h ends in j w.p. V(h,j); a busy server is a PAIR (i,j) that moves i->k at
 * D0(i,k)V(k,j)/V(i,j) and completes at D1(i,j)/V(i,j). The class memory variable is the
 * landing of the most recently STARTED draw: a start from h enters (h,j) w.p. V(h,j) and sets
 * h := j, while moves and completions keep it. For a renewal MAP the chain reduces in law to
 * PH/c. The pair law is stored as the lifted MAP D0p (conditioned moves) and
 * D1p((i,j),(j,j')) = D1(i,j)/V(i,j)*V(j,j'), equivalent in law to the original; the
 * bookkeeping read by the state handlers goes to `sn.ctmcmapc`.
 *
 * As in the reference, only a multiserver station is rewritten. A single server keeps its law and
 * gets a `single` entry: the memory is the phase the last completion left the MAP in, which the
 * handlers here did not carry before, so every draw restarted from the entry law. Pairs at c = 1
 * would be exact too but multiply the space, most of it unreachable, for no gain in the law.
 */

#include <cmath>
#include <cstddef>
#include <map>
#include <utility>
#include <vector>

#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/state_events.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace ctmc {

namespace mapc_detail {

template <class T>
bool wants_mapc(const qn::NetworkStruct<T>& sn, std::size_t ist, std::size_t r) {
    const qn::Station<T>& st = sn.stations[ist - 1];
    if (st.sched != lang::SchedStrategy::FCFS || !std::isfinite(st.nservers) || st.nservers < 1)
        return false;
    if (sn.mapc_of(ist, r) != nullptr) return false;
    const lang::Distrib<T>& d = sn.service[ist - 1][r - 1];
    if (d.disabled) return false;
    if (d.type != lang::ProcessType::MAP && d.type != lang::ProcessType::MMPP2) return false;
    return d.D0.rows() > 0 && d.D1.rows() == d.D0.rows();
}

/** Rewrite class r at station ist into pair form (c > 1) or mark it carried (c = 1); a private copy. */
template <class T>
void mapc_class(qn::NetworkStruct<T>& sn, std::size_t ist, std::size_t r) {
    lang::Distrib<T>& d = sn.service[ist - 1][r - 1];
    const T zero = num_traits<T>::from_int(0);
    const std::size_t p = d.D0.rows();
    if (sn.stations[ist - 1].nservers == 1) {
        qn::CtmcMapc<T> mc1;
        mc1.p = p;
        mc1.single = true;
        sn.ctmcmapc[std::make_pair(ist, r)] = mc1;
        return;
    }
    Matrix<T> nD0(p, p, zero);
    for (std::size_t i = 0; i < p; ++i)
        for (std::size_t j = 0; j < p; ++j) nD0(i, j) = T(zero - d.D0(i, j));
    Matrix<T> V = matmul(inverse(nD0), d.D1);
    qn::CtmcMapc<T> mc;
    mc.p = p;
    std::vector<std::vector<long>> tp(p, std::vector<long>(p, -1));
    for (std::size_t i = 0; i < p; ++i)
        for (std::size_t j = 0; j < p; ++j) {
            if (std::fabs(num_traits<T>::to_double(V(i, j))) < 1e-14) V(i, j) = zero;
            if (num_traits<T>::to_double(V(i, j)) > 0) {
                tp[i][j] = static_cast<long>(mc.pairs.size());
                mc.pairs.push_back(std::make_pair(i, j));
            }
        }
    const std::size_t np = mc.pairs.size();
    Matrix<T> D0p(np, np, zero), D1p(np, np, zero);
    mc.done.assign(np, zero);
    for (std::size_t t = 0; t < np; ++t) {
        const std::size_t i = mc.pairs[t].first, j = mc.pairs[t].second;
        mc.done[t] = T(d.D1(i, j) / V(i, j));
        D0p(t, t) = d.D0(i, i);
        for (std::size_t k = 0; k < p; ++k)
            if (k != i && num_traits<T>::to_double(d.D0(i, k)) != 0 && tp[k][j] >= 0)
                D0p(t, static_cast<std::size_t>(tp[k][j])) = T(d.D0(i, k) * V(k, j) / V(i, j));
        for (std::size_t jn = 0; jn < p; ++jn)
            if (tp[j][jn] >= 0) D1p(t, static_cast<std::size_t>(tp[j][jn])) = T(mc.done[t] * V(j, jn));
    }
    mc.V = V;
    d.D0 = D0p;  // same law, so type, mean and scv are kept
    d.D1 = D1p;
    sn.ctmcmapc[std::make_pair(ist, r)] = mc;
}

}  // namespace mapc_detail

/**
 * Rewrite a copy of `in` into `out` when it has a MAP service at an FCFS station; false,
 * leaving `out` untouched, when there is nothing to rewrite. Idempotent. A declared initial
 * state at a rewritten station is carried as its marginal and rebuilt in the pair layout.
 */
template <class T>
bool ctmc_mapc(const qn::NetworkStruct<T>& in, qn::NetworkStruct<T>& out) {
    const std::size_t R = in.nclasses;
    bool any = false;
    for (std::size_t ist = 1; ist <= in.stations.size() && !any; ++ist)
        for (std::size_t r = 1; r <= R && !any; ++r) any = mapc_detail::wants_mapc(in, ist, r);
    if (!any) return false;
    out = in;
    for (std::size_t ind = 1; ind <= in.nodes.size(); ++ind) {
        const std::size_t ist = in.nodes[ind - 1].station;
        if (ist == 0) continue;
        bool changed = false;
        for (std::size_t r = 1; r <= R; ++r)
            if (mapc_detail::wants_mapc(in, ist, r)) {
                mapc_detail::mapc_class(out, ist, r);
                changed = true;
            }
        if (!changed) continue;
        const auto sp = in.statespace.find(ind);
        if (sp == in.statespace.end() || sp->second.rows() == 0) continue;
        if (sp->second.rows() > 1)
            throw UnsupportedError("SolverCTMC cannot place a distribution over initial states on "
                                   "the MAP service of station '" + in.nodes[ind - 1].name + "'.");
        std::vector<T> row(sp->second.cols());
        for (std::size_t c = 0; c < sp->second.cols(); ++c) row[c] = sp->second(0, c);
        const std::pair<T, std::vector<T>> mg = qn::to_marginal_aggr(in, ind, row);
        Matrix<T> m(1, R, num_traits<T>::from_int(0));
        for (std::size_t r = 0; r < R; ++r) m(0, r) = mg.second[r];
        out.statespace[ind] = m;
    }
    return true;
}

}  // namespace ctmc
}  // namespace line

#endif  // LINE_SOLVERS_CTMC_SOLVER_CTMC_MAPC_H
