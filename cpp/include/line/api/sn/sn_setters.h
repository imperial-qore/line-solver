/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_SETTERS_H
#define LINE_API_SN_SN_SETTERS_H

/**
 * Ports of the sn_set_* family and of sn_refresh_process_fields
 * (matlab/src/api/sn).
 *
 * These edit an ALREADY REFRESHED struct in place. The model layer is the
 * normal way to change a network, but a sensitivity sweep, a design-of-
 * experiments driver or an optimiser wants to move one number and re-solve
 * without rebuilding and re-refreshing the whole struct -- rebuilding would
 * also renumber chains and visits, which is exactly what such a sweep must
 * hold fixed.
 *
 * WHAT `refresh` MEANS HERE. Only the caller knows whether the edit invalidated
 * the derived fields, so each setter takes the flag rather than deciding: a
 * rate or SCV change invalidates the process representation, a population
 * change invalidates the visits, and everything else invalidates nothing.
 *
 * sn_refresh_process_fields rebuilds the (D0, D1) representation from the first
 * two moments alone, by the reference's rule: SCV one is exponential, SCV below
 * one is Erlang-ceil(1/SCV), SCV above one is a two-phase hyperexponential, and
 * a hyperexponential that comes out infeasible falls back to exponential. A
 * non-positive, infinite or undefined rate is left alone -- there is no
 * representation to build.
 *
 * ARITHMETIC: field for everything but map_hyperexp, which needs a square root
 * and therefore refuses under exact arithmetic; an SCV above one is the only
 * path that reaches it.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace api {

/**
 * Port of sn_refresh_process_fields: rebuild `sn.service[ist][r]` from the
 * (rate, SCV) pair currently in `sn.rates` and `sn.scv`.
 *
 * Indices are 1-based, as everywhere in this api layer.
 */
template <class T>
void sn_refresh_process_fields(qn::NetworkStruct<T>& sn, std::size_t ist, std::size_t r) {
    const double rate = num_traits<T>::to_double(sn.rates(ist - 1, r - 1));
    if (!std::isfinite(rate) || !(rate > 0.0)) return;
    const T mean = T(num_traits<T>::from_int(1) / sn.rates(ist - 1, r - 1));
    const double scv = num_traits<T>::to_double(sn.scv(ist - 1, r - 1));
    lang::Distrib<T>& d = sn.service[ist - 1][r - 1];
    mam::Map<T> m;
    lang::ProcessType pt = lang::ProcessType::EXP;
    if (std::isnan(scv) || std::fabs(scv - 1.0) < 1e-10) {
        m = mam::map_exponential(sn.rates(ist - 1, r - 1));
    } else if (scv < 1.0) {
        const unsigned k = static_cast<unsigned>(std::max(1.0, std::ceil(1.0 / scv)));
        m = mam::map_erlang(mean, k);
        pt = lang::ProcessType::ERLANG;
    } else {
        bool ok = true;
        try {
            m = mam::map_hyperexp(mean, sn.scv(ist - 1, r - 1),
                                  num_traits<T>::from_rational(99, 100));
            pt = lang::ProcessType::HYPEREXP;
        } catch (const Error&) {
            ok = false;
        }
        if (!ok || m.D0.rows() == 0) {
            m = mam::map_exponential(sn.rates(ist - 1, r - 1));
            pt = lang::ProcessType::EXP;
        }
    }
    d.type = pt;
    d.disabled = false;
    d.D0 = m.D0;
    d.D1 = m.D1;
    d.mean = mean;
    d.scv = sn.scv(ist - 1, r - 1);
    sn.disabled[ist - 1][r - 1] = false;
}

/** Port of sn_set_service: write a (rate, SCV) pair at one (station, class). */
template <class T>
void sn_set_service(qn::NetworkStruct<T>& sn, std::size_t ist, std::size_t r, const T& rate,
                    const T& scv, bool auto_refresh = false) {
    sn.rates(ist - 1, r - 1) = rate;
    sn.scv(ist - 1, r - 1) = scv;
    if (auto_refresh) sn_refresh_process_fields(sn, ist, r);
}

/** Port of sn_set_arrival: the same, at whichever station is the Source. */
template <class T>
void sn_set_arrival(qn::NetworkStruct<T>& sn, std::size_t r, const T& rate, const T& scv,
                    bool auto_refresh = false) {
    std::size_t ist = 0;
    for (std::size_t a = 0; a < sn.nodes.size(); ++a)
        if (sn.nodes[a].nodetype == qn::NodeType::Source) ist = sn.nodes[a].station;
    if (ist == 0) throw InputError("sn_set_arrival: no Source station found in network");
    sn_set_service(sn, ist, r, rate, scv, auto_refresh);
}

/**
 * Port of sn_set_service_batch: write a whole (rate, SCV) table, skipping the
 * entries the caller left undefined.
 *
 * MATLAB marks "leave this one alone" with NaN; this port takes an explicit
 * mask, because a struct instantiated at Rational has no NaN to mark it with.
 */
template <class T>
void sn_set_service_batch(qn::NetworkStruct<T>& sn, const Matrix<T>& rates, const Matrix<T>& scvs,
                          const std::vector<std::vector<bool>>& set_rate,
                          const std::vector<std::vector<bool>>& set_scv,
                          bool auto_refresh = false) {
    std::vector<std::pair<std::size_t, std::size_t>> touched;
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t r = 0; r < sn.nclasses; ++r) {
            if (i < set_rate.size() && r < set_rate[i].size() && set_rate[i][r]) {
                sn.rates(i, r) = rates(i, r);
                touched.push_back(std::make_pair(i + 1, r + 1));
            }
            if (i < set_scv.size() && r < set_scv[i].size() && set_scv[i][r])
                sn.scv(i, r) = scvs(i, r);
        }
    if (!auto_refresh) return;
    for (std::size_t k = 0; k < touched.size(); ++k)
        sn_refresh_process_fields(sn, touched[k].first, touched[k].second);
}

/**
 * Port of sn_set_population: change a class population and the closed total.
 *
 * `auto_refresh` re-solves the visit ratios, which the reference does because a
 * class that becomes open (or closed) changes which chains are closed and
 * therefore how the visits are normalised.
 */
template <class T>
void sn_set_population(qn::NetworkStruct<T>& sn, std::size_t r, double njobs,
                       bool auto_refresh = false) {
    sn.classes[r - 1].population = njobs;
    if (auto_refresh) sn.refresh_chains();
}

/** Port of sn_set_servers. */
template <class T>
void sn_set_servers(qn::NetworkStruct<T>& sn, std::size_t ist, double nservers) {
    sn.stations[ist - 1].nservers = nservers;
}

/** Port of sn_set_priority. */
template <class T>
void sn_set_priority(qn::NetworkStruct<T>& sn, std::size_t r, int priority) {
    sn.classes[r - 1].prio = priority;
}

/**
 * Port of sn_set_fork_fanout.
 *
 * The reference writes `sn.nodeparam{f}.fanOut`, a field no MATLAB solver reads
 * back; the quantity this struct carries and the fork-join machinery does read
 * is the per-branch multiplicity `tasks_per_link`, so that is what is written.
 */
template <class T>
void sn_set_fork_fanout(qn::NetworkStruct<T>& sn, std::size_t fork_node, double fanout) {
    if (fork_node == 0 || fork_node > sn.nodes.size() ||
        sn.nodes[fork_node - 1].nodetype != qn::NodeType::Fork)
        throw InputError("sn_set_fork_fanout: node " + std::to_string(fork_node) +
                         " is not a Fork node");
    sn.nodes[fork_node - 1].tasks_per_link = fanout;
}

/** Port of sn_set_routing: replace the class-expanded stateful routing wholesale. */
template <class T>
void sn_set_routing(qn::NetworkStruct<T>& sn, const Matrix<T>& rt, bool auto_refresh = false) {
    sn.rt = rt;
    if (auto_refresh) sn.refresh_chains();
}

/** Port of sn_set_routing_prob: one entry of the class-expanded stateful routing. */
template <class T>
void sn_set_routing_prob(qn::NetworkStruct<T>& sn, std::size_t from_stateful,
                         std::size_t from_class, std::size_t to_stateful, std::size_t to_class,
                         const T& prob, bool auto_refresh = false) {
    const std::size_t K = sn.nclasses;
    sn.rt((from_stateful - 1) * K + (from_class - 1), (to_stateful - 1) * K + (to_class - 1)) =
        prob;
    if (auto_refresh) sn.refresh_chains();
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_SETTERS_H
