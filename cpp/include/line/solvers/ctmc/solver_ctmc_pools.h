/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_CTMC_SOLVER_CTMC_POOLS_H
#define LINE_SOLVERS_CTMC_SOLVER_CTMC_POOLS_H

/**
 * @file
 * @ingroup line_solvers
 * Port of `solver_ctmc_pools.m`: rewrite every heterogeneous-server station
 * (`Queue.addServerType`) into the form the CTMC state space enumerates exactly,
 * following the LDES semantics.
 *
 * A class-r job in service at a pooled station occupies one server of ONE
 * compatible pool t and is served by that pool's law, or by the station's own
 * law for r when the pool declares none. The per-class service process becomes a
 * block-diagonal phase-type law with one block per compatible pool, in ascending
 * pool order, so a phase of the class-r server block names both the pool and the
 * service phase. The bookkeeping read by `qn::after_event_station_pool` and
 * `qn::from_marginal_pool` goes to `sn.ctmcpool`. Under ALIS and FAIRNESS, when
 * some class has two or more compatible pools, the global rotating pool order is
 * part of the state: its index into `perms` is the station's single local
 * variable (nvars column 2R+1). The PHASE synchronizations a multi-phase pool law
 * needs follow from `refresh_sync`, which reads the rewritten laws.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/state_events.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace ctmc {

namespace pools_detail {

/** A renewal phase-type law is what a pool server can carry; anything else is refused by name. */
template <class T>
void require_ph(const lang::Distrib<T>& d, const mam::Map<T>& m, const std::string& name,
                const std::string& cname, const std::string& what) {
    const lang::ProcessType pt = d.type;
    if (pt == lang::ProcessType::MAP || pt == lang::ProcessType::MMPP2 ||
        pt == lang::ProcessType::MMAP || pt == lang::ProcessType::ME ||
        pt == lang::ProcessType::RAP)
        throw UnsupportedError(
            "SolverCTMC serves a heterogeneous server pool with a renewal phase-type law; class '" +
            cname + "' at station '" + name + "' has a " + std::string(lang::process_to_text(pt)) +
            " law in " + what + ".");
    const std::size_t n = m.D0.rows();
    const std::vector<T> a = mam::map_pie(m);
    double asum = 0, d1n = 0, err = 0;
    for (std::size_t j = 0; j < n; ++j) asum += num_traits<T>::to_double(a[j]);
    asum = std::max(asum, lang::GlobalConstants::FineTol);
    bool bad = false;
    for (std::size_t i = 0; i < n; ++i) {
        double ex = 0;
        for (std::size_t j = 0; j < n; ++j) ex -= num_traits<T>::to_double(m.D0(i, j));
        if (ex < -lang::GlobalConstants::FineTol) bad = true;
        for (std::size_t j = 0; j < n; ++j) {
            if (i != j && num_traits<T>::to_double(m.D0(i, j)) < -lang::GlobalConstants::FineTol) bad = true;
            const double d1 = num_traits<T>::to_double(m.D1(i, j));
            d1n += std::fabs(d1);
            err += std::fabs(d1 - ex * num_traits<T>::to_double(a[j]) / asum);
        }
    }
    if (bad || err > lang::GlobalConstants::CoarseTol * std::max(1.0, d1n))
        throw UnsupportedError(
            "SolverCTMC serves a heterogeneous server pool with a renewal phase-type law; class '" +
            cname + "' at station '" + name + "' has a correlated or matrix-exponential law in " +
            what + ".");
}

/** The features a pooled station cannot be combined with, `sub_refuse_features` of the reference. */
template <class T>
void refuse_features(const qn::NetworkStruct<T>& sn, std::size_t ind, std::size_t ist) {
    const std::size_t R = sn.nclasses;
    const qn::Station<T>& st = sn.stations[ist - 1];
    std::string why;
    bool lld = false;
    for (std::size_t j = 0; j < st.lldscaling.size(); ++j)
        if (num_traits<T>::to_double(st.lldscaling[j]) != 1.0) lld = true;
    bool retrial = false;
    const auto rit = sn.retrialparam.find(ist);
    if (rit != sn.retrialparam.end())
        for (std::size_t r = 0; r < rit->second.retrial_proc.size(); ++r)
            if (!rit->second.retrial_proc[r].disabled) retrial = true;
    bool balk = false, reneg = false, immf = false, reply = false, par = false, rr = false,
         sig = false;
    for (std::size_t r = 0; r < st.balking.size(); ++r)
        if (st.balking[r].strategy != lang::BalkingStrategy::NONE) balk = true;
    for (std::size_t r = 0; r < st.impatience.size(); ++r)
        if (st.impatience[r] != lang::ImpatienceType::NONE) reneg = true;
    if (sn.immfeed.size() >= ist)
        for (std::size_t r = 0; r < sn.immfeed[ist - 1].size(); ++r)
            if (sn.immfeed[ist - 1][r]) immf = true;
    if (sn.replyblock.size() >= ind)
        for (std::size_t r = 0; r < sn.replyblock[ind - 1].size(); ++r)
            if (sn.replyblock[ind - 1][r]) reply = true;
    for (std::size_t r = 0; r < st.server_parallelism.size(); ++r)
        if (st.server_parallelism[r] > 1) par = true;
    const std::vector<lang::RoutingStrategy>& rt = sn.nodes[ind - 1].routing;
    for (std::size_t r = 0; r < rt.size(); ++r)
        if (rt[r] == lang::RoutingStrategy::RROBIN || rt[r] == lang::RoutingStrategy::WRROBIN)
            rr = true;
    for (std::size_t s = 0; s < sn.issignal.size() && s < R && !sig; ++s) {
        if (!sn.issignal[s]) continue;
        const std::size_t col = (ind - 1) * R + s;
        for (std::size_t i = 0; i < sn.rtnodes.rows() && col < sn.rtnodes.cols(); ++i)
            if (num_traits<T>::to_double(sn.rtnodes(i, col)) > 0) { sig = true; break; }
    }
    if (lld) why = "load-dependent service";
    else if (st.cdscaling) why = "class-dependent service";
    else if (st.jdscaling) why = "joint-dependent service";
    else if (static_cast<bool>(sn.gdscaling)) why = "global dependence";
    else if (sn.breakdownparam.count(ist)) why = "server breakdowns";
    else if (retrial) why = "retrial";
    else if (balk) why = "balking";
    else if (reneg) why = "reneging";
    else if (sn.isbasblocking.size() >= ind && sn.isbasblocking[ind - 1]) why = "BAS blocking";
    else if (immf) why = "immediate feedback";
    else if (reply) why = "synchronous calls";
    else if (par) why = "server parallelism";
    else if (rr) why = "round-robin routing";
    else if (sig) why = "signals";
    if (!why.empty())
        throw UnsupportedError("SolverCTMC does not combine heterogeneous server pools with " + why +
                               " (station '" + sn.nodes[ind - 1].name + "').");
}

/** Rewrite one station; the struct is the solver's private copy. */
template <class T>
void pool_station(qn::NetworkStruct<T>& sn, std::size_t ind, std::size_t ist) {
    const std::size_t R = sn.nclasses;
    qn::Station<T>& st = sn.stations[ist - 1];
    const std::string& name = sn.nodes[ind - 1].name;
    const std::size_t NT = st.server_types.size();
    switch (st.sched) {
        case lang::SchedStrategy::FCFS:
        case lang::SchedStrategy::HOL:
        case lang::SchedStrategy::LCFS:
        case lang::SchedStrategy::LCFSPRIO:
        case lang::SchedStrategy::SIRO:
            break;
        default:
            throw UnsupportedError(
                "SolverCTMC supports heterogeneous server pools under FCFS, HOL, FCFSPRIO, LCFS, "
                "LCFSPRIO and SIRO; station '" + name + "' uses " +
                std::string(lang::sched_to_text(st.sched)) + ".");
    }
    refuse_features(sn, ind, ist);
    if (sn.nvars_of(ind) > 0)
        throw UnsupportedError(
            "Station '" + name + "' combines heterogeneous server pools with another feature that "
            "owns a local variable (polling, BAS blocking, breakdown, a modulated process); "
            "SolverCTMC cannot represent both.");
    qn::CtmcPool<T> pool;
    pool.ntypes = NT;
    pool.policy = st.hetero_policy;
    pool.count.assign(NT, 0.0);
    pool.compat.assign(NT, std::vector<bool>(R, false));
    for (std::size_t t = 0; t < NT; ++t) {
        pool.count[t] = st.server_types[t].count;
        const std::vector<bool>& c = st.server_types[t].compatible;
        for (std::size_t r = 0; r < R; ++r) pool.compat[t][r] = c.empty() || (r < c.size() && c[r]);
    }
    pool.pools.assign(R, std::vector<std::size_t>());
    pool.off.assign(R, std::vector<std::size_t>());
    pool.len.assign(R, std::vector<std::size_t>());
    pool.alpha.assign(R, std::vector<std::vector<T>>());
    pool.exit.assign(R, std::vector<std::vector<T>>());
    pool.D0.assign(R, std::vector<Matrix<T>>());
    pool.fsfrate.assign(NT, std::vector<double>(R, 0.0));
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t r = 0; r < R; ++r) {
        const lang::Distrib<T>& base = sn.service[ist - 1][r];
        const bool baseOk = !base.disabled && base.type != lang::ProcessType::DISABLED;
        bool hasPoolLaw = false;
        for (std::size_t t = 0; t < NT; ++t) {
            const std::vector<lang::Distrib<T>>& sv = st.server_types[t].service;
            hasPoolLaw = hasPoolLaw || (pool.compat[t][r] && r < sv.size() && !sv[r].disabled);
        }
        if (!baseOk && !hasPoolLaw) continue;  // the class is not served here
        for (std::size_t t = 0; t < NT; ++t)
            if (pool.compat[t][r]) pool.pools[r].push_back(t + 1);
        const std::string& cname = sn.classes[r].name;
        if (pool.pools[r].empty())
            throw UnsupportedError("Station '" + name +
                                   "' declares no server pool compatible with class '" + cname +
                                   "', so a job of that class would wait forever.");
        mam::Map<T> basemap;
        if (baseOk) {
            basemap = lang::dist_to_map(base);
            require_ph(base, basemap, name, cname, "its default service");
        }
        std::vector<mam::Map<T>> blocks;
        std::size_t shift = 0;
        for (std::size_t k = 0; k < pool.pools[r].size(); ++k) {
            const std::size_t t = pool.pools[r][k];
            const std::vector<lang::Distrib<T>>& sv = st.server_types[t - 1].service;
            mam::Map<T> law;
            if (r < sv.size() && !sv[r].disabled) {
                law = lang::dist_to_map(sv[r]);
                require_ph(sv[r], law, name, cname, "pool '" + st.server_types[t - 1].name + "'");
            } else {
                if (!baseOk)
                    throw UnsupportedError(
                        "Server pool '" + st.server_types[t - 1].name + "' of station '" + name +
                        "' accepts class '" + cname + "' but declares no law for it, and the "
                        "station's own service for the class is disabled.");
                law = basemap;
            }
            const std::size_t n = law.D0.rows();
            std::vector<T> a = mam::map_pie(law);
            T as = zero;
            for (std::size_t j = 0; j < n; ++j) as += a[j];
            for (std::size_t j = 0; j < n; ++j) a[j] = T(a[j] / as);
            std::vector<T> ex(n, zero);
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < n; ++j) ex[i] -= law.D0(i, j);
            pool.off[r].push_back(shift);
            pool.len[r].push_back(n);
            pool.alpha[r].push_back(a);
            pool.exit[r].push_back(ex);
            pool.D0[r].push_back(law.D0);
            pool.fsfrate[t - 1][r] = 1.0 / num_traits<T>::to_double(mam::map_mean(law));
            blocks.push_back(law);
            shift += n;
        }
        Matrix<T> D0x(shift, shift, zero), D1x(shift, shift, zero);
        std::vector<T> exx(shift, zero);
        for (std::size_t k = 0; k < blocks.size(); ++k) {
            const std::size_t o = pool.off[r][k];
            for (std::size_t i = 0; i < pool.len[r][k]; ++i) {
                exx[o + i] = pool.exit[r][k][i];
                for (std::size_t j = 0; j < pool.len[r][k]; ++j) D0x(o + i, o + j) = blocks[k].D0(i, j);
            }
        }
        // a completion restarts in the first block's entry law; which block a job takes
        // is decided by the event handler, so this D1 only fixes the renewal form
        for (std::size_t i = 0; i < shift; ++i)
            for (std::size_t j = 0; j < pool.len[r][0]; ++j) D1x(i, j) = T(exx[i] * pool.alpha[r][0][j]);
        lang::Distrib<T> d;
        d.type = lang::ProcessType::PH;
        d.disabled = false;
        d.D0 = D0x;
        d.D1 = D1x;
        lang::dist_refresh_moments(d);
        sn.service[ist - 1][r] = d;
        if (sn.disabled.size() >= ist && sn.disabled[ist - 1].size() > r) sn.disabled[ist - 1][r] = false;
    }
    std::vector<std::size_t> ncls(NT, 0);
    for (std::size_t t = 0; t < NT; ++t)
        for (std::size_t r = 0; r < R; ++r) ncls[t] += pool.compat[t][r] ? 1 : 0;
    pool.alfsorder.resize(NT);
    for (std::size_t t = 0; t < NT; ++t) pool.alfsorder[t] = t + 1;
    std::stable_sort(pool.alfsorder.begin(), pool.alfsorder.end(),
                     [&](std::size_t a, std::size_t b) { return ncls[a - 1] < ncls[b - 1]; });
    bool multi = false;
    for (std::size_t r = 0; r < R; ++r) multi = multi || pool.pools[r].size() >= 2;
    pool.rotate = (pool.policy == lang::HeteroSchedPolicy::ALIS ||
                   pool.policy == lang::HeteroSchedPolicy::FAIRNESS) && multi;
    if (pool.rotate) {
        std::vector<std::size_t> p(NT);
        for (std::size_t t = 0; t < NT; ++t) p[t] = t + 1;
        do {
            pool.perms.push_back(p);  // lexicographic, identity first
        } while (std::next_permutation(p.begin(), p.end()));
        sn.nvars[ind - 1][2 * R] = 1;
    }
    double total = 0;
    for (std::size_t t = 0; t < NT; ++t) total += pool.count[t];
    st.nservers = total;  // the pools are the server bank, as in LDES
    sn.ctmcpool[ind] = pool;
}

}  // namespace pools_detail

/** True when the struct has a heterogeneous-server station that `ctmc_pools` rewrites. */
template <class T>
bool ctmc_has_pools(const qn::NetworkStruct<T>& sn) {
    for (std::size_t ind = 1; ind <= sn.nodes.size(); ++ind) {
        const std::size_t ist = sn.nodes[ind - 1].station;
        if (ist == 0 || sn.stations[ist - 1].server_types.empty() || sn.ctmcpool.count(ind)) continue;
        const lang::SchedStrategy s = sn.stations[ist - 1].sched;
        // PAS/OI model heterogeneous compatible servers through the OI rank rate, not pools
        if (s == lang::SchedStrategy::PAS || s == lang::SchedStrategy::OI) continue;
        return true;
    }
    return false;
}

/**
 * Rewrite a copy of `in` into `out` when it has pooled stations; false, leaving
 * `out` untouched, when there is nothing to rewrite. Idempotent: a station that
 * already carries `ctmcpool` is skipped. A declared initial state at a pooled
 * station is carried as its marginal and rebuilt in the pooled layout.
 */
template <class T>
bool ctmc_pools(const qn::NetworkStruct<T>& in, qn::NetworkStruct<T>& out) {
    if (!ctmc_has_pools(in)) return false;
    if (in.isfjaugmented)
        throw UnsupportedError(
            "SolverCTMC does not combine heterogeneous server pools with fork-join: the tag "
            "augmentation adds sibling classes the pool compatibility does not name.");
    out = in;
    const std::size_t R = in.nclasses;
    for (std::size_t ind = 1; ind <= in.nodes.size(); ++ind) {
        const std::size_t ist = in.nodes[ind - 1].station;
        if (ist == 0 || in.stations[ist - 1].server_types.empty() || in.ctmcpool.count(ind)) continue;
        const lang::SchedStrategy s = in.stations[ist - 1].sched;
        if (s == lang::SchedStrategy::PAS || s == lang::SchedStrategy::OI) continue;
        pools_detail::pool_station(out, ind, ist);
        const auto sp = in.statespace.find(ind);
        if (sp == in.statespace.end() || sp->second.rows() == 0) continue;
        if (sp->second.rows() > 1)
            throw UnsupportedError("SolverCTMC cannot place a distribution over initial states on "
                                   "the heterogeneous server pools of station '" +
                                   in.nodes[ind - 1].name + "'.");
        std::vector<T> row(sp->second.cols());
        for (std::size_t c = 0; c < sp->second.cols(); ++c) row[c] = sp->second(0, c);
        const std::pair<T, std::vector<T>> mg = qn::to_marginal_aggr(in, ind, row);
        Matrix<T> m(1, R, num_traits<T>::from_int(0));
        for (std::size_t r = 0; r < R; ++r) m(0, r) = mg.second[r];
        out.statespace[ind] = m;
    }
    return true;
}

/** Utilization at a pooled station: the mean number of busy servers of the class over the bank. */
template <class T>
bool ctmc_pool_station(const qn::NetworkStruct<T>& sn, std::size_t ist) {
    const std::size_t ind = sn.node_of_station(ist);
    return ind != 0 && !sn.ctmcpool.empty() && sn.ctmcpool.count(ind) != 0;
}

}  // namespace ctmc
}  // namespace line

#endif  // LINE_SOLVERS_CTMC_SOLVER_CTMC_POOLS_H
