/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Port of the event half of MATLAB's `+State` package: the successor states an
 * event produces at one node, with their rates and probabilities. This is what
 * turns the enumerated state space of `state.h` into a generator.
 *
 * THE ACTIVE / PASSIVE CONVENTION. An event is ACTIVE at the node that
 * schedules it and PASSIVE at the node that receives it -- a DEP at one station
 * IS the ARV at the next. Only the active half knows a rate, so the passive
 * half returns the sentinel -1 and the generator assembly substitutes the
 * active rate. The sentinel is unambiguous because a rate is never negative.
 *
 * WHY PROBABILITIES ARE SEPARATE FROM RATES. One event can have several
 * successors: an entering job picks its service phase, a signal picks the job
 * it removes, a random-order queue picks whom to serve. The rate belongs to the
 * event and the probability to the choice, so a single (state, event) pair
 * yields a ROW of successors and the generator entry is rate * prob.
 */
#ifndef LINE_LANG_QN_STATE_EVENTS_H
#define LINE_LANG_QN_STATE_EVENTS_H

#include <cmath>
#include <cstddef>
#include <limits>
#include <algorithm>
#include <functional>
#include <utility>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/api/mam/map_moment.h"
#include "line/lang/qn/polling_info.h"
#include "line/lang/qn/state.h"
#include "line/util/error.h"

namespace line {
namespace qn {

using lang::EventType;

/**
 * What one event produces at one node: the successor rows, their rates and
 * their probabilities, all three the same length.
 */
template <class T>
struct EventOutcome {
    std::vector<std::vector<T>> space;  ///< successor local state rows
    std::vector<T> rate;                ///< per-row rate, -1 on a passive half
    std::vector<T> prob;                ///< per-row probability of the choice
    /// START annotation: the 1-based classes that BEGIN or RESUME holding a
    /// server on each successor row. An instantaneous tag on the arc the row
    /// already carries, never an event of its own, so no rate, probability or
    /// state depends on it. Usually empty; kept as a list so one arc can start
    /// several jobs (a region release cascade does).
    std::vector<std::vector<std::size_t>> start;
    /// PREEMPT annotation: the 1-based classes pushed back into the buffer.
    std::vector<std::vector<std::size_t>> preempt;
    bool empty() const { return space.empty(); }
};

/**
 * Tag the successor row just appended to OUT: START_CLS begins service on it
 * and PREEMPT_CLS is displaced by it, either 0 for none. The tag vectors are
 * grown to match `space`, so only the arcs that carry a tag need to say
 * anything and every other row is an empty list.
 */
template <typename T>
inline void tag_last(EventOutcome<T>& out, std::size_t start_cls, std::size_t preempt_cls) {
    if (out.space.empty()) return;
    out.start.resize(out.space.size());
    out.preempt.resize(out.space.size());
    if (start_cls) out.start.back().push_back(start_cls);
    if (preempt_cls) out.preempt.back().push_back(preempt_cls);
}

/**
 * Bring the tag vectors up to one entry per successor, so a caller can index
 * them exactly like `space`.
 */
template <typename T>
inline void pad_tags(EventOutcome<T>& out) {
    out.start.resize(out.space.size());
    out.preempt.resize(out.space.size());
}

/**
 * Port of `State.toMarginalAggr`: the job counts of one node's state row,
 * without the per-phase detail `to_marginal` also computes.
 *
 * It is NOT simply a projection of `to_marginal`: it accepts stateful
 * non-stations (whose row is a plain per-class count ahead of the local
 * variables), and it leaves the preemptive families out of its buffer switch,
 * so for those it reports only the jobs IN SERVICE. That asymmetry is the
 * reference's, and it is load-bearing -- the arrival branch uses this to test
 * for room, where counting a preempted job twice would refuse a valid arrival.
 *
 * @param sn      the network struct
 * @param ind     NODE index (1-based), as in the reference
 * @param state_i the node's state row
 * @return (ni, nir): total jobs and jobs per class
 */
template <class T>
std::pair<T, std::vector<T>> to_marginal_aggr(const NetworkStruct<T>& sn, std::size_t ind,
                                              const std::vector<T>& state_i) {
    const std::size_t R = sn.nclasses;
    const T zero = num_traits<T>::from_int(0);
    std::vector<T> nir(R, zero);
    if (ind == 0 || ind > sn.nodes.size())
        throw InputError("to_marginal_aggr: node index is out of range");
    const NodeDef& nd = sn.nodes[ind - 1];
    const std::size_t ist = nd.station;
    const std::size_t nvar = sn.nvars_of(ind);

    // A Join of an FJ-augmented struct is a station whose row is a bare per-class
    // count: it holds jobs waiting to synchronize, with no buffer/phase split and
    // no service at all, so the station path below would read its counts as
    // buffer tags.
    if (sn.isfjaugmented && nd.nodetype == NodeType::Join && state_i.size() >= R) {
        T ni = zero;
        for (std::size_t r = 0; r < R; ++r) {
            nir[r] = state_i[state_i.size() - R + r];
            ni += nir[r];
        }
        return std::make_pair(ni, nir);
    }

    // A stateful non-station carries a per-class count ahead of its local
    // variables. A node with nothing but bookkeeping (a Router's round-robin
    // pointer) still has to report R zeros, so callers can index nir[r].
    if (ist == 0) {
        const std::size_t bufw = state_i.size() > nvar ? state_i.size() - nvar : 0;
        for (std::size_t r = 0; r < R && r < bufw; ++r) nir[r] = state_i[r];
        T ni = zero;
        for (std::size_t r = 0; r < R; ++r) ni += nir[r];
        return std::make_pair(ni, nir);
    }

    // A Source reports zero, not Inf: its jobs are external, and the EXT
    // sentinel `to_marginal` returns describes the encoding, not a count that
    // an arrival branch could compare against a capacity.
    if (nd.nodetype == NodeType::Source) return std::make_pair(zero, nir);

    std::vector<std::size_t> K(R, 1), Ks(R, 0);
    std::size_t srvw = 0;
    for (std::size_t r = 0; r < R; ++r) {
        K[r] = sn.phasessz_of(ist, r + 1);
        Ks[r] = srvw;
        srvw += K[r];
    }
    if (state_i.size() < nvar + srvw)
        throw InputError("to_marginal_aggr: state row is narrower than its server block");
    const std::size_t srv0 = state_i.size() - nvar - srvw;

    for (std::size_t r = 0; r < R; ++r)
        for (std::size_t k = 0; k < K[r]; ++k) nir[r] += state_i[srv0 + Ks[r] + k];

    const SchedStrategy sched = sn.stations[ist - 1].sched;
    if (sched == SchedStrategy::EXT) {
        // Reached only by an EXT station that is not a Source node, since the
        // Source returns zero above. It carries the same clamp as `to_marginal`
        // and for the same reason: an exact type has no infinity, and building
        // one from a double Inf throws rather than saturating.
        const T ext = num_traits<T>::is_exact
                          ? num_traits<T>::from_double(GlobalConstants::MaxInt)
                          : num_traits<T>::from_double(
                                std::numeric_limits<double>::infinity());
        for (std::size_t r = 0; r < R; ++r) nir[r] = ext;
    } else if (state_detail::buffer_is_class_tag(sched)) {
        for (std::size_t r = 0; r < R; ++r) {
            const T tag = num_traits<T>::from_int(static_cast<long>(r + 1));
            for (std::size_t b = 0; b < srv0; ++b)
                if (state_i[b] == tag) nir[r] += num_traits<T>::from_int(1);
        }
    } else if (state_detail::buffer_is_tag_phase_pairs(sched)) {
        // Only the EVEN positions are class tags; the odd ones record the phase
        // each preempted job was interrupted in. Without this arm a preempted
        // job was invisible here and nir counted the server alone, unlike
        // `to_marginal`, which has carried the paired decode all along.
        if (srv0 > 1)
            for (std::size_t r = 0; r < R; ++r) {
                const T tag = num_traits<T>::from_int(static_cast<long>(r + 1));
                for (std::size_t b = 0; b < srv0; b += 2)
                    if (state_i[b] == tag) nir[r] += num_traits<T>::from_int(1);
            }
    } else if (state_detail::buffer_is_per_class_count(sched)) {
        for (std::size_t r = 0; r < R && r < srv0; ++r) nir[r] += state_i[r];
    } else if (sched == SchedStrategy::PAS || sched == SchedStrategy::OI) {
        // A PAS / OI row is the ORDERED JOB LIST and nothing else: entry `b` is
        // the 1-based class of the job in position b, 0 for an empty slot, which
        // is the encoding `after_event_station_pas` reads and writes. It has no
        // buffer/server split at all, so the server-block sum taken above
        // counted the LAST LIST POSITION as a phase occupancy -- it is discarded
        // here and the whole row is scanned instead.
        //
        // WITHOUT THIS BRANCH the station reports a queue length of about zero
        // while the jobs are demonstrably in it, in EVERY solver that reduces a
        // state through this function (SolverCTMC's `solver_ctmc_avg_from_pi`
        // and SolverSSA's serial analyzer both do), and the capacity filter in
        // `after_event_station_arv` compares that zero against the station's
        // bound. A wrong NUMBER, never an error.
        for (std::size_t r = 0; r < R; ++r) nir[r] = zero;
        const std::size_t w = state_i.size() > nvar ? state_i.size() - nvar : 0;
        for (std::size_t b = 0; b < w; ++b) {
            const double v = num_traits<T>::to_double(state_i[b]);
            const long tag = static_cast<long>(v + 0.5);
            if (tag >= 1 && static_cast<std::size_t>(tag) <= R)
                nir[tag - 1] += num_traits<T>::from_int(1);
        }
    }

    // A disabled class holds no jobs whatever the row says. A Place is exempt:
    // its tokens are not services, so it has no rate to be disabled.
    if (nd.nodetype != NodeType::Place)
        for (std::size_t r = 0; r < R; ++r)
            if (sn.disabled[ist - 1][r]) nir[r] = zero;

    T ni = zero;
    for (std::size_t r = 0; r < R; ++r) ni += nir[r];
    return std::make_pair(ni, nir);
}

/**
 * Port of `State.isPhysicalCapacity`: true when the bound at (ist, class) is a
 * PHYSICAL capacity rather than a state-space CUTOFF on an open class.
 *
 * The distinction decides what a refused arrival means. The producer's
 * capacity arguments have the cutoff folded in -- `solver_ssa` overwrites
 * cap/classcap with min(cutoff, physical) -- so at a cutoff boundary they read
 * finite even with no physical cap. Treating that as physical would turn a
 * state-space TRUNCATION into a self-loop loss, which reports a wrong arrival
 * rate and perturbs the sample path.
 *
 * The in-producer signal is the DROP RULE: `refreshCapacity` sets a non-WAITQ
 * rule exactly when the capacity is physical, and a cutoff-bounded open class
 * keeps the WAITQ default.
 */
template <class T>
bool is_physical_capacity(const NetworkStruct<T>& sn, std::size_t ist, std::size_t cls) {
    if (sn.droprule.size() < ist || sn.droprule[ist - 1].size() < cls) return false;
    const DropStrategy dr = sn.droprule[ist - 1][cls - 1];
    return dr != DropStrategy::WAITQ && static_cast<int>(dr) != 0;
}

/**
 * Port of `State.arrivalIsLost`: true when an arrival that finds no room is
 * LOST, false when it must BLOCK the upstream instead. Every refusal path
 * branches on this, and the two outcomes are encoded differently:
 *
 *   LOST    -> leave the state UNCHANGED, a self-loop. The event still fires,
 *              so the OFFERED job reaches the arrival-rate statistic and the
 *              loss appears as ArvR - Tput. A self-loop cancels on the
 *              generator diagonal, so it cannot move the stationary law.
 *   BLOCKED -> return NO rows. That disables the upstream departure until room
 *              frees, which is what the become-blocked edge tests for.
 *
 * The rule is the CLASS TYPE, not the drop rule: a closed network's population
 * is a defining invariant, so a closed job can never be dropped. An explicit
 * BAS/BBS/RSRD rule asks for blocking for any class.
 */
template <class T>
bool arrival_is_lost(const NetworkStruct<T>& sn, std::size_t ist, std::size_t cls) {
    if (sn.droprule.size() >= ist && sn.droprule[ist - 1].size() >= cls) {
        const DropStrategy dr = sn.droprule[ist - 1][cls - 1];
        if (dr == DropStrategy::BAS || dr == DropStrategy::BBS || dr == DropStrategy::RSRD)
            return false;  // the user asked for blocking explicitly
    }
    // The rule above sees only THIS station's declaration. Under the upstream
    // declaration form the BAS rule sits on the blocking station, not on the
    // destination where the refusal happens, so the destination side is recorded
    // separately by `refresh_bas_blocking`. Without this branch an open class
    // refused here would be declared lost, the become-blocked edge would never
    // fire, and the blocking station would behave as if its destination were
    // unbounded.
    if (sn.isbasdestination.size() >= ist && sn.isbasdestination[ist - 1].size() >= cls &&
        sn.isbasdestination[ist - 1][cls - 1])
        return false;
    // Open -> lost, closed -> blocked. This decides only what happens once the
    // arrival has already been refused, never whether it is refused.
    const double nj = sn.njobs()[cls - 1];
    return !std::isfinite(nj);
}

/** How a station's state row splits into [buffer | server | local vars]. */
template <class T>
struct RowLayout {
    std::vector<std::size_t> K;   ///< phases per class
    std::vector<std::size_t> Ks;  ///< offset of class r's phase block
    std::size_t srvw = 0;         ///< total server width
    std::size_t nvar = 0;         ///< local-variable width
    std::size_t bufw = 0;         ///< buffer width, the only discipline-dependent part
};

template <class T>
RowLayout<T> row_layout(const NetworkStruct<T>& sn, std::size_t ind, std::size_t width) {
    const std::size_t R = sn.nclasses;
    const std::size_t ist = sn.nodes[ind - 1].station;
    RowLayout<T> L;
    L.K.assign(R, 1);
    L.Ks.assign(R, 0);
    for (std::size_t r = 0; r < R; ++r) {
        L.K[r] = sn.phasessz_of(ist, r + 1);
        L.Ks[r] = L.srvw;
        L.srvw += L.K[r];
    }
    L.nvar = sn.nvars_of(ind);
    if (width < L.nvar + L.srvw)
        throw InputError("after_event_station: the state row is narrower than its server block");
    L.bufw = width - L.nvar - L.srvw;
    return L;
}

/**
 * The entry-phase distribution `pie{ist}{class}`: which phase a service STARTS
 * in. This is `map_pie`, the equilibrium embedded at DEPARTURE instants, and
 * NOT `map_prob`, the time-stationary law of D0+D1 -- the two differ whenever
 * the process is not exponential (for Erlang-2, entry is [1,0] while the
 * time-stationary law is [0.5,0.5]).
 *
 * A Place has no service process, so its "phases" carry no rate; the reference
 * falls back on a uniform choice there rather than leaving the vector NaN.
 */
template <class T>
std::vector<T> entry_phase_dist(const NetworkStruct<T>& sn, std::size_t ist, std::size_t cls) {
    const std::size_t nph = sn.phases_of(ist, cls);
    const lang::Distrib<T>& d = sn.service[ist - 1][cls - 1];
    std::vector<T> pie;
    if (d.D0.rows() == nph && d.D1.rows() == nph && nph > 0 && !d.disabled) {
        mam::Map<T> m;
        m.D0 = d.D0;
        m.D1 = d.D1;
        try {
            pie = mam::map_pie(m);
        } catch (const Error&) {
            pie.clear();  // a zero-rate process has no entry law; fall back below
        }
    }
    bool ok = pie.size() == nph;
    if (ok) {
        T s = num_traits<T>::from_int(0);
        for (std::size_t k = 0; k < nph; ++k) s += pie[k];
        ok = num_traits<T>::to_double(s) > 0 && std::isfinite(num_traits<T>::to_double(s));
    }
    if (!ok) {
        pie.assign(nph, num_traits<T>::from_int(0));
        if (nph > 0) {
            const T u = T(num_traits<T>::from_int(1) / num_traits<T>::from_int(
                                                          static_cast<long>(nph)));
            for (std::size_t k = 0; k < nph; ++k) pie[k] = u;
        }
    }
    return pie;
}

/** Where node `ind` keeps its reply-block counters inside the local vars. */
struct ReplyBlockInfo {
    std::vector<std::size_t> classes;  ///< 1-based calling classes holding a block
    std::vector<std::size_t> slot;     ///< slot[r-1] = 0-based column, or npos
    std::size_t width = 0;
};

/** Defined below; the polling branches of ARV, DEP and SWITCH use these. */
template <class T>
void polling_get(const PollingInfo<T>& pi, const std::vector<T>& var, std::size_t srvclass,
                 std::size_t& pos, std::size_t& swk, long& ctr);
template <class T>
std::vector<T> polling_set(const PollingInfo<T>& pi, std::vector<T> var, std::size_t pos,
                           std::size_t swk, long ctr);
template <class T>
void polling_next(const PollingInfo<T>& pi, std::size_t pos, const std::vector<long>& nbuf,
                  std::size_t R, bool arrived, std::size_t& q, int& mode, long& budget);
template <class T>
void polling_land(const NetworkStruct<T>& sn, std::size_t ist, const PollingInfo<T>& pi,
                  std::size_t q, int mode, long budget, const std::vector<T>& buf,
                  const std::vector<T>& srv, const std::vector<T>& var, const RowLayout<T>& L,
                  std::vector<std::vector<T>>& rows, std::vector<T>& probs);

/** Defined below; the reply block subtracts held servers in the ARV branch. */
template <class T>
double reply_blocked(const NetworkStruct<T>& sn, std::size_t ind, const std::vector<T>& var);

/** Defined below; the departure branch records a server held for a reply. */
template <class T>
ReplyBlockInfo reply_block_info(const NetworkStruct<T>& sn, std::size_t ind);

/** Defined below; the ARV and DEP branches divert to it before any slicing. */
template <class T>
EventOutcome<T> after_event_station_pas(const NetworkStruct<T>& sn, std::size_t ind,
                                        const std::vector<T>& inspace, EventType event,
                                        std::size_t cls);

/**
 * Port of the ARV branch of `State.afterEventStation`: an arriving class-`cls`
 * job joins node `ind`, whose local state is `inspace`.
 *
 * The event is PASSIVE -- the upstream departure sets the rate -- so every row
 * returned carries the -1 sentinel. It is nevertheless the branch with the most
 * successors, because the entering job chooses its service phase, and under the
 * preemptive disciplines it also chooses which job to displace.
 *
 * ONE ROW IN, MANY ROWS OUT. The reference threads a whole matrix of input
 * rows through this handler and partitions them with logical masks
 * (`idle_srv`, `all_busy_srv`). Every caller in the CTMC and SSA paths passes a
 * single row, so those masks degenerate to a branch, which is what this port
 * writes. The successor set is identical.
 */
template <class T>
EventOutcome<T> after_event_station_arv(const NetworkStruct<T>& sn, std::size_t ind,
                                        const std::vector<T>& inspace, std::size_t cls) {
    const std::size_t R = sn.nclasses;
    const std::size_t ist = sn.nodes[ind - 1].station;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T minus_one = num_traits<T>::from_int(-1);
    EventOutcome<T> out;
    if (ist == 0) throw InputError("after_event_station_arv: node is not a station");
    const SchedStrategy sched = sn.stations[ist - 1].sched;
    // A pass-and-swap station has no buffer/server split at all, so it must be
    // diverted before any slicing happens.
    if (sched == SchedStrategy::PAS || sched == SchedStrategy::OI)
        return after_event_station_pas(sn, ind, inspace, EventType::ARV, cls);
    const RowLayout<T> L = row_layout(sn, ind, inspace.size());

    // A Place holds a marking, not a service facility: it has no servers and no
    // phases, so an arriving token only increments the class marking. Running
    // it through the scheduling branches below would write into a phase slot
    // the row does not have, widening the state so it no longer matches the
    // enumerated space -- and the arrival would then be silently dropped.
    if (sn.nodes[ind - 1].nodetype == NodeType::Place) {
        std::vector<T> row = inspace;
        const double cap = sn.classcap[ist - 1][cls - 1];
        if (num_traits<T>::to_double(inspace[cls - 1]) < cap) {
            row[cls - 1] += one;
            out.space.push_back(row);
            out.rate.push_back(minus_one);
            out.prob.push_back(one);
        } else {
            // Place full: the arrival is blocked and lost, with no state change.
            out.space.push_back(row);
            out.rate.push_back(minus_one);
            out.prob.push_back(zero);
        }
        return out;
    }

    const std::pair<T, std::vector<T>> mg = to_marginal_aggr(sn, ind, inspace);
    const double ni = num_traits<T>::to_double(mg.first);
    const double nir_c = num_traits<T>::to_double(mg.second[cls - 1]);
    const double cap_i = sn.cap[ist - 1];
    const double ccap = sn.classcap[ist - 1][cls - 1];
    const double S = sn.stations[ist - 1].nservers;
    const std::vector<T> pentry = entry_phase_dist(sn, ist, cls);

    for (std::size_t kentry = 0; kentry < L.K[cls - 1]; ++kentry) {
        std::vector<T> buf(inspace.begin(), inspace.begin() + L.bufw);
        std::vector<T> srv(inspace.begin() + L.bufw, inspace.begin() + L.bufw + L.srvw);
        std::vector<T> var(inspace.begin() + L.bufw + L.srvw, inspace.end());
        std::vector<std::vector<T>> cand;  // (buf, srv, var) triples, flattened below
        std::vector<T> cand_prob;
        // START/PREEMPT tag of each candidate: the class that takes a server on
        // that row and the class it displaces, 0 for neither. Filtered with the
        // rows themselves at the capacity gate below.
        std::vector<std::size_t> cand_start, cand_preempt;

        double occ = 0;  // jobs currently in service
        for (std::size_t j = 0; j < srv.size(); ++j) occ += num_traits<T>::to_double(srv[j]);

        if (sched == SchedStrategy::EXT) {
            // A Source accepts a virtual arrival from the Sink for any open
            // class: the reservoir is unbounded, so the state does not move.
            if (!std::isfinite(sn.njobs()[cls - 1])) {
                out.space.push_back(inspace);
                out.rate.push_back(zero);
                out.prob.push_back(one);
                return out;
            }
            continue;
        }

        if (sched == SchedStrategy::PS || sched == SchedStrategy::INF ||
            sched == SchedStrategy::DPS || sched == SchedStrategy::GPS ||
            sched == SchedStrategy::PSPRIO || sched == SchedStrategy::DPSPRIO ||
            sched == SchedStrategy::GPSPRIO || sched == SchedStrategy::LPS) {
            // Every job is in service at once, so the arrival never queues.
            const std::size_t col = L.Ks[cls - 1] + kentry;
            std::size_t started = 0;
            if (num_traits<T>::to_double(srv[col]) < ccap) {
                srv[col] += one;
                started = cls;  // the job enters service at once
                cand_prob.push_back(pentry[kentry]);
            } else {
                cand_prob.push_back(zero);
            }
            std::vector<T> row = buf;
            row.insert(row.end(), srv.begin(), srv.end());
            row.insert(row.end(), var.begin(), var.end());
            cand.push_back(row);
            cand_start.push_back(started);
            cand_preempt.push_back(0);
        } else if (sched == SchedStrategy::POLLING) {
            // The CONTROLLER decides who is served, not the arrival: a job
            // joins its class buffer and waits for the server to walk to it,
            // even when the facility is idle, because the server is then in a
            // switchover. The one exception is a PARKED server, which only
            // arises with an empty station and immediate switchovers: it
            // reaches the arriving job in zero time and opens a visit at once.
            const PollingInfo<T> pinfo = polling_info(sn, ind);
            std::size_t srvclass = 0;
            for (std::size_t r = 1; r <= R; ++r) {
                double tot = 0;
                for (std::size_t p = 0; p < L.K[r - 1]; ++p)
                    tot += num_traits<T>::to_double(srv[L.Ks[r - 1] + p]);
                if (tot > 0) { srvclass = r; break; }
            }
            std::size_t pos = 0, swk = 0;
            long ctr = 0;
            polling_get(pinfo, var, srvclass, pos, swk, ctr);
            if (srvclass == 0 && swk == 0) {
                std::vector<long> nbuf(R, 0);
                for (std::size_t r = 0; r < R && r < L.bufw; ++r)
                    nbuf[r] = static_cast<long>(num_traits<T>::to_double(buf[r]));
                nbuf[cls - 1] += 1;
                std::size_t q = 0;
                int mode = 0;
                long budget = 0;
                polling_next(pinfo, pos, nbuf, R, true, q, mode, budget);
                srv[L.Ks[cls - 1] + kentry] += one;
                var = polling_set(pinfo, var, q, 0, budget);
                cand_start.push_back(cls);  // a parked server takes it at once
            } else {
                buf[cls - 1] += one;
                cand_start.push_back(0);
            }
            cand_preempt.push_back(0);
            cand_prob.push_back(pentry[kentry]);
            std::vector<T> row = buf;
            row.insert(row.end(), srv.begin(), srv.end());
            row.insert(row.end(), var.begin(), var.end());
            cand.push_back(row);
        } else if (sched == SchedStrategy::SIRO || sched == SchedStrategy::SEPT ||
                   sched == SchedStrategy::LEPT) {
            // Test the SERVER occupancy, not the total count: the two agree in
            // work-conserving states, but an immediate-feedback self-loop
            // transiently leaves an idle server with a non-empty buffer, and
            // the fed-back job must re-enter the vacated server.
            if (occ < S) {
                srv[L.Ks[cls - 1] + kentry] += one;
                cand_start.push_back(cls);
            } else {
                buf[cls - 1] += one;
                cand_start.push_back(0);
            }
            cand_preempt.push_back(0);
            cand_prob.push_back(pentry[kentry]);
            std::vector<T> row = buf;
            row.insert(row.end(), srv.begin(), srv.end());
            row.insert(row.end(), var.begin(), var.end());
            cand.push_back(row);
        } else if (state_detail::buffer_is_class_tag(sched)) {
            // ORDERED BUFFER. An idle server takes the job; otherwise it goes
            // to the first empty slot, and if there is none the arrival is
            // refused -- which is a LOSS or a BLOCK, never a silent drop.
            // Servers held for a pending REPLY are NOT available to an
            // arriving job, so a job may already have to wait while the raw
            // occupancy is below the server count. Zero for every model
            // without reply signals.
            const double seff = S - reply_blocked(sn, ind, var);
            if (occ < seff) {
                srv[L.Ks[cls - 1] + kentry] += one;
                std::vector<T> row = buf;
                row.insert(row.end(), srv.begin(), srv.end());
                row.insert(row.end(), var.begin(), var.end());
                cand.push_back(row);
                cand_prob.push_back(pentry[kentry]);
                cand_start.push_back(cls);
                cand_preempt.push_back(0);
            } else {
                std::size_t slot = 0;  // 1-based index of the LAST empty slot
                for (std::size_t b = 0; b < L.bufw; ++b)
                    if (num_traits<T>::to_double(buf[b]) == 0) slot = b + 1;
                // A structurally free column is not enough: the CAPACITY must
                // also permit the placement. Gating on the capacity only for a
                // PHYSICAL bound keeps a state-space cutoff a truncation --
                // firing the gate at a cutoff would turn it into a self-loop.
                const bool has_room = is_physical_capacity(sn, ist, cls)
                                          ? (ni < cap_i && nir_c < ccap)
                                          : true;
                if (slot > 0 && has_room) {
                    buf[slot - 1] = num_traits<T>::from_int(static_cast<long>(cls));
                    std::vector<T> row = buf;
                    row.insert(row.end(), srv.begin(), srv.end());
                    row.insert(row.end(), var.begin(), var.end());
                    cand.push_back(row);
                    cand_prob.push_back(pentry[kentry]);
                    cand_start.push_back(0);  // the job waits in the buffer
                    cand_preempt.push_back(0);
                } else if (arrival_is_lost(sn, ist, cls) &&
                           is_physical_capacity(sn, ist, cls)) {
                    // LOST: keep the row unchanged, so the event still fires
                    // and the OFFERED job reaches the arrival-rate statistic.
                    // A self-loop cancels on the generator diagonal, so the
                    // stationary law cannot move.
                    //
                    // ONLY A DECLARED BUFFER MAY LOSE A JOB. `slot == 0` also
                    // fires when the row is merely as wide as the STATE-SPACE
                    // CUTOFF let it be, and there the reference emits no
                    // successor at all -- the state above the cutoff is absent,
                    // not refused. Emitting the self-loop there costs nothing on
                    // the generator (the diagonal absorbs it) and everything on
                    // the RATES, which count the loop as a departure of the
                    // upstream Source: on mqn_multiserver_fcfs 20 such loops
                    // moved Source Tput from 0.24763 to 0.26040, an arrival rate
                    // no job ever carried, and left ArvR above Tput at a station
                    // that drops nothing. `has_room` above already draws exactly
                    // this line for the capacity test.
                    cand.push_back(inspace);
                    cand_prob.push_back(pentry[kentry]);
                    cand_start.push_back(0);  // the job is lost: it starts nothing
                    cand_preempt.push_back(0);
                }
                // BLOCKED: emit nothing, which disables the upstream departure
                // until room frees. That absence is what the become-blocked
                // edge tests for.
            }
        } else if (state_detail::buffer_is_tag_phase_pairs(sched)) {
            // PREEMPTIVE FAMILY. The buffer holds [class, phase] pairs, because
            // a displaced job must remember where it was interrupted. An idle
            // server simply takes the arrival; a busy one forces a CHOICE of
            // victim, so the event has one successor per (class, phase) in
            // service, weighted by that server's share of the occupancy.
            if (occ < S) {
                srv[L.Ks[cls - 1] + kentry] += one;
                std::vector<T> row = buf;
                row.insert(row.end(), srv.begin(), srv.end());
                row.insert(row.end(), var.begin(), var.end());
                cand.push_back(row);
                cand_prob.push_back(pentry[kentry]);
                cand_start.push_back(cls);
                cand_preempt.push_back(0);
            } else {
                // Priority-awareness is a property of the DECLARED policy, never
                // of the data: inferring it from the class priorities turned
                // plain LCFSPR/FCFSPR into something that is neither the base
                // policy nor the PRIO variant.
                const bool prio_aware = sched == SchedStrategy::FCFSPRPRIO ||
                                        sched == SchedStrategy::FCFSPIPRIO ||
                                        sched == SchedStrategy::LCFSPRPRIO ||
                                        sched == SchedStrategy::LCFSPIPRIO;
                const bool lcfs_family = sched == SchedStrategy::LCFSPRPRIO ||
                                         sched == SchedStrategy::LCFSPIPRIO;
                // PR resumes the victim in the phase it held; PI restarts it
                // from the entry phase. That single value is the whole
                // difference between the two families in this branch.
                const bool resume = sched == SchedStrategy::LCFSPR ||
                                    sched == SchedStrategy::LCFSPRPRIO ||
                                    sched == SchedStrategy::FCFSPR ||
                                    sched == SchedStrategy::FCFSPRPRIO;
                bool can_preempt_any = false;
                for (std::size_t cp = 1; cp <= R; ++cp) {
                    if (prio_aware) {
                        // Across priority groups a strictly higher-priority
                        // arrival preempts. WITHIN a group the base discipline
                        // decides: LCFS-PR keeps the NEWEST job in service, so
                        // an equal-priority arrival preempts; FCFS-PR never
                        // lets an arrival preempt.
                        const int pa = sn.classes[cls - 1].prio;
                        const int pv = sn.classes[cp - 1].prio;
                        if (lcfs_family ? (pa > pv) : (pa >= pv)) continue;
                    }
                    for (std::size_t pp = 0; pp < L.K[cp - 1]; ++pp) {
                        const std::size_t vcol = L.Ks[cp - 1] + pp;
                        const double busy = num_traits<T>::to_double(srv[vcol]);
                        if (busy <= 0) continue;
                        can_preempt_any = true;
                        std::vector<T> b2 = buf, s2 = srv;
                        s2[vcol] -= one;
                        s2[L.Ks[cls - 1] + kentry] += one;
                        // Rightmost empty PAIR, which is where a displaced job
                        // is stored; the class column is one left of the zero
                        // the scan finds.
                        std::size_t slot = 0;
                        for (std::size_t b = 0; b < L.bufw; ++b)
                            if (num_traits<T>::to_double(b2[b]) == 0) slot = b;
                        if (slot == 0) continue;  // no room to hold the victim
                        b2[slot - 1] = num_traits<T>::from_int(static_cast<long>(cp));
                        b2[slot] = resume ? num_traits<T>::from_int(static_cast<long>(pp + 1))
                                          : one;
                        std::vector<T> row = b2;
                        row.insert(row.end(), s2.begin(), s2.end());
                        row.insert(row.end(), var.begin(), var.end());
                        cand.push_back(row);
                        // the displaced job leaves the server and the arriving
                        // one takes it, on the same arc
                        cand_start.push_back(cls);
                        cand_preempt.push_back(cp);
                        // The victim is drawn uniformly among the jobs in
                        // service, so its share of the occupancy weights the
                        // successor.
                        cand_prob.push_back(T(pentry[kentry] * num_traits<T>::from_double(busy / occ)));
                    }
                }
                // Every busy server holds a job this arrival may not displace,
                // so the job WAITS instead. Without this the loop above emits
                // nothing, the arrival transition does not exist at all, and
                // the class can never enter a busy station -- its queue is
                // then silently understated. It is stored as a (class, entry
                // phase) pair exactly as a preempted job is, so promotion
                // resumes it from that phase.
                if (prio_aware && !can_preempt_any) {
                    std::size_t slot = 0;
                    for (std::size_t b = 0; b < L.bufw; ++b)
                        if (num_traits<T>::to_double(buf[b]) == 0) slot = b;
                    if (slot > 0) {
                        std::vector<T> b2 = buf;
                        b2[slot - 1] = num_traits<T>::from_int(static_cast<long>(cls));
                        b2[slot] = num_traits<T>::from_int(static_cast<long>(kentry + 1));
                        std::vector<T> row = b2;
                        row.insert(row.end(), srv.begin(), srv.end());
                        row.insert(row.end(), var.begin(), var.end());
                        cand.push_back(row);
                        cand_prob.push_back(pentry[kentry]);
                        cand_start.push_back(0);  // it preempts nothing and waits
                        cand_preempt.push_back(0);
                    }
                }
            }
        } else {
            throw UnsupportedError(
                std::string("after_event_station_arv: the ") + lang::sched_to_text(sched) +
                " discipline is not ported yet");
        }

        // The capacity filter of the reference: drop any successor that would
        // exceed the station or class bound.
        for (std::size_t c = 0; c < cand.size(); ++c) {
            const std::pair<T, std::vector<T>> og = to_marginal_aggr(sn, ind, cand[c]);
            if (num_traits<T>::to_double(og.second[cls - 1]) > ccap) continue;
            if (num_traits<T>::to_double(og.first) > cap_i) continue;
            out.space.push_back(cand[c]);
            out.rate.push_back(minus_one);
            out.prob.push_back(cand_prob[c]);
            // the capacity gate drops rows, so the tags are attached here
            tag_last(out, c < cand_start.size() ? cand_start[c] : 0,
                     c < cand_preempt.size() ? cand_preempt[c] : 0);
        }
    }
    pad_tags(out);
    return out;
}

/**
 * `sn.mu` and `sn.phi` for one (station, class), derived as MATLAB's
 * `Markovian.getMu` / `getPhi` derive them from the (D0, D1) pair:
 *
 *   mu(k)  = -D0(k,k)                  total exit rate of phase k
 *   phi(k) = sum_j D1(k,j) / -D0(k,k)  probability the exit COMPLETES service
 *
 * so mu*phi is the departure rate and mu*(1-phi) the phase-advance rate. An
 * Immediate process has D0(1,1) = 0 and takes phi = 1, since every exit of a
 * zero-duration service is a completion.
 */
template <class T>
std::pair<std::vector<T>, std::vector<T>> phase_rates(const NetworkStruct<T>& sn,
                                                      std::size_t ist, std::size_t cls) {
    const lang::Distrib<T>& d = sn.service[ist - 1][cls - 1];
    const std::size_t n = sn.phases_of(ist, cls);
    std::vector<T> mu(n, num_traits<T>::from_int(0)), phi(n, num_traits<T>::from_int(1));
    if (d.D0.rows() != n || d.D1.rows() != n) return std::make_pair(mu, phi);
    for (std::size_t k = 0; k < n; ++k) {
        const T dk = T(-d.D0(k, k));
        mu[k] = dk;
        if (num_traits<T>::to_double(dk) == 0) {
            phi[k] = num_traits<T>::from_int(1);  // Immediate: every exit completes
            continue;
        }
        T s = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < n; ++j) s += d.D1(k, j);
        phi[k] = T(s / dk);
    }
    return std::make_pair(mu, phi);
}

/** The limited-load-dependent multiplier at population `n`, 1 when unset. */
template <class T>
T lld_factor(const NetworkStruct<T>& sn, std::size_t ist, double n) {
    const std::vector<T>& s = sn.stations[ist - 1].lldscaling;
    if (s.empty()) return num_traits<T>::from_int(1);
    // Beyond the declared levels the scaling holds at its last value, which is
    // what "limited" load-dependence means: the curve is flat past the limit.
    if (!std::isfinite(n) || n >= static_cast<double>(s.size())) return s.back();
    if (n < 1) return num_traits<T>::from_int(1);
    return s[static_cast<std::size_t>(n) - 1];
}

/**
 * Port of `State.cdclassfactor`: the class-dependence multiplier of a class-`cls`
 * rate at the per-class population `nir`.
 *
 * `cdscaling` maps a 1 x R population vector to the R dimensionless scalings
 * beta_r(n); the component of the class whose service is firing is the factor.
 * A handle returning a SCALAR is the neutral case and yields 1 for every class,
 * which is why the read is clamped to the vector's last entry rather than
 * indexed blindly -- the reference's `v(min(class, numel(v)))`.
 *
 * `jdscaling` is folded in multiplicatively here, exactly as
 * `State.afterEventInit` folds eta_i into the effective per-station handle.
 * Doing it at the point of use rather than by rewriting the struct keeps the two
 * fields distinguishable for the product-form tests elsewhere.
 */
template <class T>
T cd_factor(const NetworkStruct<T>& sn, std::size_t ist, const std::vector<T>& nir,
            std::size_t cls) {
    const T one = num_traits<T>::from_int(1);
    const Station<T>& st = sn.stations[ist - 1];
    if (!st.cdscaling && !st.jdscaling) return one;
    T f = one;
    const CdScaling<T>* handles[2] = {&st.cdscaling, &st.jdscaling};
    for (std::size_t h = 0; h < 2; ++h) {
        const CdScaling<T>& fun = *handles[h];
        if (!fun) continue;
        const std::vector<T> v = fun(nir);
        if (v.empty())
            throw InputError(
                "cd_factor: the class-dependence map returned an empty scaling vector");
        f = T(f * v[std::min(cls, v.size()) - 1]);
    }
    return f;
}

/**
 * The population the *PRIO disciplines actually share the server among.
 *
 * While every job fits in a server nobody is waiting and precedence is moot, so
 * this is the plain marginal. Once the station saturates only the MOST URGENT
 * group present is served (a LOWER `classprio` value is more urgent in LINE),
 * and both the share and the load-dependent lookup are taken over that group
 * alone -- `nirprio` / `niprio` in `State.afterEventStation`.
 */
template <class T>
struct PrioPop {
    std::vector<T> nir;      ///< `nirprio` when masked, the plain marginal otherwise
    double ni = 0;           ///< `niprio` when masked, the plain total otherwise
    bool masked = false;     ///< whether the saturated-station mask was applied
    bool served = true;      ///< false when `cls` is not in the most urgent group
};

/** Compute the *PRIO effective population; a no-op for every other discipline. */
template <class T>
PrioPop<T> prio_pop(const NetworkStruct<T>& sn, std::size_t ist, const Marginal<T>& m,
                    std::size_t cls, double ni, double S) {
    const std::size_t R = sn.nclasses;
    const SchedStrategy sched = sn.stations[ist - 1].sched;
    PrioPop<T> p;
    p.nir = m.nir;
    p.ni = ni;
    const bool prio_aware = sched == SchedStrategy::PSPRIO ||
                            sched == SchedStrategy::DPSPRIO ||
                            sched == SchedStrategy::GPSPRIO;
    if (!prio_aware || ni <= S) return p;
    double best = std::numeric_limits<double>::infinity();
    for (std::size_t r = 0; r < R; ++r)
        if (num_traits<T>::to_double(m.nir[r]) > 0)
            best = std::min(best, static_cast<double>(sn.classes[r].prio));
    if (static_cast<double>(sn.classes[cls - 1].prio) != best) {
        p.served = false;
        return p;
    }
    p.masked = true;
    p.ni = 0;
    for (std::size_t r = 0; r < R; ++r) {
        if (static_cast<double>(sn.classes[r].prio) != best)
            p.nir[r] = num_traits<T>::from_int(0);
        p.ni += num_traits<T>::to_double(p.nir[r]);
    }
    return p;
}

/** Defined below; DEP and PHASE must share one definition of the share. */
template <class T>
T service_share(const NetworkStruct<T>& sn, std::size_t ist, const Marginal<T>& m,
                std::size_t cls, double ni, double S);

/**
 * Port of the DEP branch of `State.afterEventStation`: a class-`cls` job
 * completes service at station `ind`.
 *
 * This is the ACTIVE half, so unlike ARV it carries real rates, and the rate is
 * where the disciplines actually differ -- the state update is nearly the same
 * for all of them. Two rate conventions appear, and they are not
 * interchangeable:
 *
 *   mu(k)*phi(k)*kir  -- the processor-sharing family, where the completion
 *                        rate is the phase rate times the share of the server
 *   D1(k,kdest)*kir   -- the FCFS family, where the MAP matrix already encodes
 *                        both the completion and the phase the NEXT service
 *                        starts in, so the destination phase is enumerated
 */
template <class T>
EventOutcome<T> after_event_station_dep(const NetworkStruct<T>& sn, std::size_t ind,
                                        const std::vector<T>& inspace, std::size_t cls,
                                        bool no_promote = false) {
    const std::size_t R = sn.nclasses;
    const std::size_t ist = sn.nodes[ind - 1].station;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    EventOutcome<T> out;
    if (ist == 0) throw InputError("after_event_station_dep: node is not a station");
    const SchedStrategy sched = sn.stations[ist - 1].sched;
    if (sched == SchedStrategy::PAS || sched == SchedStrategy::OI)
        return after_event_station_pas(sn, ind, inspace, EventType::DEP, cls);
    const RowLayout<T> L = row_layout(sn, ind, inspace.size());
    const double S = sn.stations[ist - 1].nservers;

    std::vector<std::size_t> ph(R, 1), shift(R, 0);
    for (std::size_t r = 0; r < R; ++r) {
        ph[r] = L.K[r];
        shift[r] = L.Ks[r];
    }
    const Marginal<T> m = to_marginal(sn, ist, inspace, ph, shift, L.nvar);
    if (num_traits<T>::to_double(m.sir[cls - 1]) <= 0) return out;  // nothing to depart

    const std::pair<std::vector<T>, std::vector<T>> mp = phase_rates(sn, ist, cls);
    const std::vector<T>& mu = mp.first;
    const std::vector<T>& phi = mp.second;
    const lang::Distrib<T>& d = sn.service[ist - 1][cls - 1];

    double ni = 0;
    for (std::size_t r = 0; r < R; ++r) ni += num_traits<T>::to_double(m.nir[r]);
    // Both dependence factors multiply every rate this branch emits, so they are
    // folded into ONE scalar named `lld` for the arithmetic below. Their
    // ARGUMENTS differ at a saturated *PRIO station, and the reference is not
    // uniform about it: the load-dependent lookup takes the priority-masked
    // total `niprio` for all three *PRIO disciplines, while the class-dependence
    // handle takes the masked vector `nirprio` only for DPSPRIO and GPSPRIO and
    // the plain marginal for PSPRIO (afterEventStation.m:709-711 vs :746-748 and
    // :808-810). That asymmetry is reproduced rather than smoothed: smoothing it
    // would move every reported metric on a model that has both.
    const PrioPop<T> pp = prio_pop(sn, ist, m, cls, ni, S);
    const bool cd_takes_prio =
        pp.masked && (sched == SchedStrategy::DPSPRIO || sched == SchedStrategy::GPSPRIO);
    const T lld = T(lld_factor(sn, ist, pp.ni) *
                    cd_factor(sn, ist, cd_takes_prio ? pp.nir : m.nir, cls));

    // A retrial station does NOT promote from the orbit on completion: an
    // orbiting job re-enters only through a RETRY at the retrial rate. The
    // same suppression serves an immediate-feedback self-loop, where the
    // departing job holds the server for its own re-arrival.
    bool suppress_promote = no_promote;
    {
        const typename std::map<std::size_t, RetrialParam<T>>::const_iterator rit =
            sn.retrialparam.find(ist);
        if (rit != sn.retrialparam.end())
            for (std::size_t r = 0; r < rit->second.retrial_proc.size(); ++r)
                if (!rit->second.retrial_proc[r].disabled) { suppress_promote = true; break; }
    }

    for (std::size_t k = 0; k < L.K[cls - 1]; ++k) {
        std::vector<T> buf(inspace.begin(), inspace.begin() + L.bufw);
        std::vector<T> srv(inspace.begin() + L.bufw, inspace.begin() + L.bufw + L.srvw);
        std::vector<T> var(inspace.begin() + L.bufw + L.srvw, inspace.end());
        const std::size_t col = L.Ks[cls - 1] + k;
        if (num_traits<T>::to_double(srv[col]) <= 0) continue;
        const T kir = m.kir[cls - 1][k];

        if (sched == SchedStrategy::EXT) {
            // A Source EMITS an arrival. Its reservoir is unbounded, so no job
            // count changes; what moves is the MODULATING PHASE, from k to
            // kentry at rate D1(k, kentry). For an exponential source that is
            // the single entry lambda and the state is unchanged, which is why
            // the arrival stream is memoryless; for a MAP it is exactly the
            // correlation the source is there to produce.
            if (!std::isfinite(sn.njobs()[cls - 1])) {
                for (std::size_t ke = 0; ke < L.K[cls - 1]; ++ke) {
                    const T arv = d.D1(k, ke);
                    if (num_traits<T>::to_double(arv) <= 0) continue;
                    std::vector<T> row = inspace;
                    row[L.bufw + L.Ks[cls - 1] + k] -= one;
                    row[L.bufw + L.Ks[cls - 1] + ke] += one;
                    out.space.push_back(row);
                    out.rate.push_back(T(lld * arv));
                    out.prob.push_back(one);
                }
            }
        } else if (sched == SchedStrategy::INF || sched == SchedStrategy::PS ||
            sched == SchedStrategy::LPS || sched == SchedStrategy::DPS ||
            sched == SchedStrategy::GPS || sched == SchedStrategy::PSPRIO ||
            sched == SchedStrategy::DPSPRIO || sched == SchedStrategy::GPSPRIO) {
            srv[col] -= one;
            // The same share PHASE uses: a completion and an internal phase
            // advance are driven by the identical fraction of the server, so
            // they must never be computed two different ways.
            const T rate = T(mu[k] * phi[k] * kir * service_share(sn, ist, m, cls, ni, S));
            std::vector<T> row = buf;
            row.insert(row.end(), srv.begin(), srv.end());
            row.insert(row.end(), var.begin(), var.end());
            out.space.push_back(row);
            out.rate.push_back(T(lld * rate));
            out.prob.push_back(one);
        } else if (state_detail::buffer_is_class_tag(sched) && sched != SchedStrategy::FCFS) {
            // HOL, LCFS and LCFSPRIO share FCFS's ordered buffer but NOT its
            // rate convention: the completion rate is mu*phi*kir, summed over
            // destination phases rather than enumerated. What separates the
            // three is only WHICH waiting job is promoted.
            const bool has_waiting = ni > S && !suppress_promote;
            const T rate = T(mu[k] * phi[k] * kir);
            srv[col] -= one;
            if (!has_waiting) {
                std::vector<T> row = buf;
                row.insert(row.end(), srv.begin(), srv.end());
                row.insert(row.end(), var.begin(), var.end());
                out.space.push_back(row);
                out.rate.push_back(T(lld * rate));
                out.prob.push_back(one);
                continue;
            }
            // Position of the job that starts service, 0-based; L.bufw = none.
            std::size_t pos = L.bufw;
            if (sched == SchedStrategy::LCFS) {
                // Plain LCFS is NOT priority-aware: it always promotes the most
                // recent arrival. Arrivals fill the rightmost empty slot, so
                // the newest job is the FIRST nonzero column. Branching here on
                // the class priorities would silently turn every LCFS station
                // with distinct priorities into an LCFSPRIO one.
                for (std::size_t b = 0; b < L.bufw; ++b)
                    if (num_traits<T>::to_double(buf[b]) != 0) { pos = b; break; }
            } else {
                // HOL and LCFSPRIO serve the highest-priority waiting group
                // first; in LINE a LOWER classprio value is more urgent. Within
                // the group HOL takes the oldest job (rightmost) and LCFSPRIO
                // the newest (leftmost), which is how each relates to its
                // non-priority base discipline.
                double best = std::numeric_limits<double>::infinity();
                for (std::size_t b = 0; b < L.bufw; ++b) {
                    const double v = num_traits<T>::to_double(buf[b]);
                    if (v <= 0) continue;
                    const double p = sn.classes[static_cast<std::size_t>(v) - 1].prio;
                    if (p < best) best = p;
                }
                if (std::isfinite(best))
                    for (std::size_t b = 0; b < L.bufw; ++b) {
                        const double v = num_traits<T>::to_double(buf[b]);
                        if (v <= 0) continue;
                        if (sn.classes[static_cast<std::size_t>(v) - 1].prio != best) continue;
                        pos = b;
                        if (sched == SchedStrategy::LCFSPRIO) break;  // leftmost
                    }
            }
            if (pos == L.bufw) continue;
            const std::size_t hc = static_cast<std::size_t>(num_traits<T>::to_double(buf[pos]));
            std::vector<T> b2 = buf;
            if (sched == SchedStrategy::LCFS) {
                // LCFS clears the slot IN PLACE, leaving a hole; the next
                // arrival refills it, since arrivals seek the rightmost empty
                // slot. The priority variants instead close the gap.
                b2[pos] = zero;
            } else {
                for (std::size_t b = pos; b > 0; --b) b2[b] = b2[b - 1];
                b2[0] = zero;
            }
            const std::vector<T> pentry = entry_phase_dist(sn, ist, hc);
            for (std::size_t ke = 0; ke < L.K[hc - 1]; ++ke) {
                std::vector<T> s3 = srv;
                s3[L.Ks[hc - 1] + ke] += one;
                std::vector<T> row = b2;
                row.insert(row.end(), s3.begin(), s3.end());
                row.insert(row.end(), var.begin(), var.end());
                out.space.push_back(row);
                out.rate.push_back(T(lld * rate * pentry[ke]));
                out.prob.push_back(one);
                tag_last(out, hc, 0);  // the promoted job takes the freed server
            }
        } else if (state_detail::buffer_is_class_tag(sched)) {
            // FCFS. D1(k, kdest) is the completion rate that leaves the process
            // in phase kdest, so the destination phase has to be enumerated
            // rather than summed away.
            //
            // A synchronous call: this departing job KEEPS its server until its
            // REPLY returns here, so the server is not handed to a waiting job;
            // it is recorded as held in the reply block instead. Servers already
            // held that way are likewise unavailable, so a job can be waiting
            // while the raw occupancy is below the server count.
            std::vector<T> var2 = var;
            bool holds_reply = false;
            if (sn.replyblock.size() >= ind && sn.replyblock[ind - 1].size() >= cls &&
                sn.replyblock[ind - 1][cls - 1]) {
                const ReplyBlockInfo ri = reply_block_info(sn, ind);
                const std::size_t sl = ri.slot[cls - 1];
                if (sl != static_cast<std::size_t>(-1) && sl < var2.size()) {
                    var2[sl] += one;
                    holds_reply = true;
                }
            }
            const double nb = reply_blocked(sn, ind, var);
            const bool has_waiting = ni > (S - nb) && !suppress_promote && !holds_reply;
            for (std::size_t kd = 0; kd < L.K[cls - 1]; ++kd) {
                const T rate = T(d.D1(k, kd) * kir);
                std::vector<T> s2 = srv;
                s2[col] -= one;
                if (!has_waiting) {
                    std::vector<T> row = buf;
                    row.insert(row.end(), s2.begin(), s2.end());
                    row.insert(row.end(), var2.begin(), var2.end());
                    out.space.push_back(row);
                    out.rate.push_back(T(lld * rate));
                    out.prob.push_back(one);
                    continue;
                }
                // Promote the head of the buffer. The head is the LAST column:
                // the buffer shifts RIGHT as jobs join, so the oldest job sits
                // at the end -- which is what makes the discipline first-come.
                const double headv = num_traits<T>::to_double(buf[L.bufw - 1]);
                if (headv <= 0) continue;
                const std::size_t hc = static_cast<std::size_t>(headv);
                std::vector<T> b2(L.bufw, zero);
                for (std::size_t b = 1; b < L.bufw; ++b) b2[b] = buf[b - 1];
                const std::vector<T> pentry = entry_phase_dist(sn, ist, hc);
                for (std::size_t ke = 0; ke < L.K[hc - 1]; ++ke) {
                    std::vector<T> s3 = s2;
                    s3[L.Ks[hc - 1] + ke] += one;
                    std::vector<T> row = b2;
                    row.insert(row.end(), s3.begin(), s3.end());
                    row.insert(row.end(), var2.begin(), var2.end());
                    const T r3 = T(lld * rate * pentry[ke]);
                    out.space.push_back(row);
                    out.rate.push_back(r3);
                    // A branch that cannot happen carries probability zero, not
                    // one: a PH whose entry vector does not reach every phase
                    // would otherwise contribute phantom departures.
                    out.prob.push_back(num_traits<T>::to_double(r3) == 0 ? zero : one);
                    tag_last(out, hc, 0);  // the head of the buffer takes the freed server
                }
            }
        } else if (sched == SchedStrategy::POLLING) {
            // A completion ends the VISIT unless the discipline still allows
            // another job of the same class; when it ends, the server walks the
            // cyclic order to the next tangible controller state.
            const PollingInfo<T> pinfo = polling_info(sn, ind);
            const T rate = T(mu[k] * phi[k] * kir);
            if (num_traits<T>::to_double(rate) <= 0) continue;
            std::size_t pos = 0, swk = 0;
            long ctr = 0;
            polling_get(pinfo, var, cls, pos, swk, ctr);
            // No job can complete while the server is walking.
            if (swk != 0) continue;
            std::vector<long> nbuf(R, 0);
            for (std::size_t r = 0; r < R && r < L.bufw; ++r)
                nbuf[r] = static_cast<long>(num_traits<T>::to_double(buf[r]));
            srv[col] -= one;
            long ctrnext = 0;
            bool goon = false;
            switch (pinfo.ptype) {
                case lang::PollingType::EXHAUSTIVE:
                    ctrnext = 0;
                    goon = nbuf[cls - 1] > 0;  // the visit ends when it drains
                    break;
                case lang::PollingType::GATED:
                    ctrnext = ctr - 1;  // one of the gated jobs completed
                    goon = ctrnext > 0;
                    break;
                case lang::PollingType::KLIMITED:
                    ctrnext = ctr - 1;  // one of the K permitted services used
                    goon = ctrnext > 0 && nbuf[cls - 1] > 0;
                    break;
                case lang::PollingType::DECREMENTING:
                    ctrnext = ctr;  // the target level is fixed for the visit
                    goon = nbuf[cls - 1] > ctr;
                    break;
            }
            std::size_t q = cls;
            int mode = 1;
            long budget = ctrnext;
            if (!goon) polling_next(pinfo, cls, nbuf, R, false, q, mode, budget);
            std::vector<std::vector<T>> rows;
            std::vector<T> probs;
            polling_land(sn, ist, pinfo, q, mode, budget, buf, srv, var, L, rows, probs);
            for (std::size_t j = 0; j < rows.size(); ++j) {
                out.space.push_back(rows[j]);
                out.rate.push_back(T(lld * rate * probs[j]));
                out.prob.push_back(one);
                // mode 1 opens a visit, pulling a waiting class-q job into the
                // server; a switchover or a park starts nobody
                tag_last(out, mode == 1 ? q : 0, 0);
            }
        } else if (state_detail::buffer_is_tag_phase_pairs(sched)) {
            // THE PREEMPT-RESUME / PREEMPT-INDEPENDENT FAMILY, all eight members.
            // Every arm of the reference carries the same rate law, mu*phi*kir
            // summed over destination phases (afterEventStation.m:1045, :1081,
            // :1118, :1150, :1188, :1234, :1286, :1332); what separates them is
            // WHICH waiting job is promoted and in WHICH phase it restarts.
            const T rate = T(mu[k] * phi[k] * kir);
            srv[col] -= one;
            const bool prio_aware = sched == SchedStrategy::FCFSPRPRIO ||
                                    sched == SchedStrategy::FCFSPIPRIO ||
                                    sched == SchedStrategy::LCFSPRPRIO ||
                                    sched == SchedStrategy::LCFSPIPRIO;
            const bool lcfs = sched == SchedStrategy::LCFSPR ||
                              sched == SchedStrategy::LCFSPI ||
                              sched == SchedStrategy::LCFSPRPRIO ||
                              sched == SchedStrategy::LCFSPIPRIO;
            // PR resumes the promoted job in the phase it was interrupted in;
            // PI discards that phase and restarts from the entry distribution.
            // That single value is the whole PR-vs-PI difference.
            const bool resume = sched == SchedStrategy::LCFSPR ||
                                sched == SchedStrategy::LCFSPRPRIO ||
                                sched == SchedStrategy::FCFSPR ||
                                sched == SchedStrategy::FCFSPRPRIO;
            // Class column, 0-based and hence EVEN; L.bufw means "nobody waits".
            std::size_t pos = L.bufw;
            if (ni > S && !suppress_promote && L.bufw >= 2) {
                double best = std::numeric_limits<double>::infinity();
                if (prio_aware)
                    // Only the class columns carry a class tag: reading the
                    // phase columns too would let a phase index masquerade as a
                    // class and win the group, which is what the reference's
                    // FCFSPIPRIO/LCFSPIPRIO arms do (:1290, :1336) while its
                    // PR-PRIO arms correctly restrict to `class_cols` (:1194).
                    for (std::size_t b = 0; b + 1 < L.bufw; b += 2) {
                        const double v = num_traits<T>::to_double(buf[b]);
                        if (v <= 0) continue;
                        const double p = sn.classes[static_cast<std::size_t>(v) - 1].prio;
                        if (p < best) best = p;
                    }
                for (std::size_t b = 0; b + 1 < L.bufw; b += 2) {
                    const double v = num_traits<T>::to_double(buf[b]);
                    if (v <= 0) continue;
                    if (prio_aware && sn.classes[static_cast<std::size_t>(v) - 1].prio != best)
                        continue;
                    pos = b;
                    // The buffer is newest-first, so LCFS takes the FIRST such
                    // pair (:1052 colfirstnnz) and FCFS the LAST (:1121
                    // colLastNnz). LCFS does NOT take the rightmost slot.
                    if (lcfs) break;
                }
            }
            if (pos == L.bufw) {
                std::vector<T> row = buf;
                row.insert(row.end(), srv.begin(), srv.end());
                row.insert(row.end(), var.begin(), var.end());
                out.space.push_back(row);
                out.rate.push_back(T(lld * rate));
                out.prob.push_back(one);
                continue;
            }
            const std::size_t hc = static_cast<std::size_t>(num_traits<T>::to_double(buf[pos]));
            const std::size_t kst =
                static_cast<std::size_t>(num_traits<T>::to_double(buf[pos + 1]));
            if (hc == 0 || hc > R || kst == 0 || kst > L.K[hc - 1])
                throw InputError("after_event_station_dep: station '" +
                                 sn.stations[ist - 1].name +
                                 "' holds a waiting job whose [class, phase] pair is malformed");
            // Close the hole by padding a whole EMPTY PAIR on the left, keeping
            // the buffer right-aligned as `from_marginal` enumerates it and as
            // afterEventStation.m:1267-1274 requires. Removing the pair in place
            // (:1055, :1125) or padding one slot on each side (:1222) leaves a
            // layout the enumerator never emits, so the successor is unreachable
            // and the generator becomes reducible.
            std::vector<T> b2(L.bufw, zero);
            for (std::size_t b = 2; b <= pos + 1; ++b) b2[b] = buf[b - 2];
            for (std::size_t b = pos + 2; b < L.bufw; ++b) b2[b] = buf[b];
            if (resume) {
                std::vector<T> s3 = srv;
                s3[L.Ks[hc - 1] + kst - 1] += one;
                std::vector<T> row = b2;
                row.insert(row.end(), s3.begin(), s3.end());
                row.insert(row.end(), var.begin(), var.end());
                out.space.push_back(row);
                out.rate.push_back(T(lld * rate));
                out.prob.push_back(one);
                tag_last(out, hc, 0);  // the promoted job resumes on the freed server
            } else {
                const std::vector<T> pentry = entry_phase_dist(sn, ist, hc);
                for (std::size_t ke = 0; ke < L.K[hc - 1]; ++ke) {
                    std::vector<T> s3 = srv;
                    s3[L.Ks[hc - 1] + ke] += one;
                    std::vector<T> row = b2;
                    row.insert(row.end(), s3.begin(), s3.end());
                    row.insert(row.end(), var.begin(), var.end());
                    out.space.push_back(row);
                    out.rate.push_back(T(lld * rate * pentry[ke]));
                    out.prob.push_back(one);
                    tag_last(out, hc, 0);
                }
            }
        } else if (state_detail::buffer_is_per_class_count(sched)) {
            // A per-class-count buffer promotes by COUNT, since no order is
            // recorded. SIRO picks the next job at random, so every waiting
            // class is a distinct successor with its share of the queue.
            srv[col] -= one;
            const T rate = T(mu[k] * phi[k] * kir);
            double waiting = 0;
            for (std::size_t r = 0; r < R; ++r) waiting += num_traits<T>::to_double(buf[r]);
            if (waiting <= 0 || suppress_promote) {
                std::vector<T> row = buf;
                row.insert(row.end(), srv.begin(), srv.end());
                row.insert(row.end(), var.begin(), var.end());
                out.space.push_back(row);
                out.rate.push_back(T(lld * rate));
                out.prob.push_back(one);
                continue;
            }
            for (std::size_t r = 1; r <= R; ++r) {
                const double nb = num_traits<T>::to_double(buf[r - 1]);
                if (nb <= 0) continue;
                const std::vector<T> pentry = entry_phase_dist(sn, ist, r);
                for (std::size_t ke = 0; ke < L.K[r - 1]; ++ke) {
                    std::vector<T> b2 = buf, s2 = srv;
                    b2[r - 1] -= one;
                    s2[L.Ks[r - 1] + ke] += one;
                    std::vector<T> row = b2;
                    row.insert(row.end(), s2.begin(), s2.end());
                    row.insert(row.end(), var.begin(), var.end());
                    out.space.push_back(row);
                    out.rate.push_back(T(lld * rate));
                    out.prob.push_back(T(num_traits<T>::from_double(nb / waiting) * pentry[ke]));
                    tag_last(out, r, 0);  // the drawn waiting class takes the server
                }
            }
        } else {
            throw UnsupportedError(
                std::string("after_event_station_dep: the ") + lang::sched_to_text(sched) +
                " discipline is not ported yet");
        }
    }
    pad_tags(out);
    return out;
}

/**
 * The fraction of the station's capacity a class-`cls` job in phase `k`
 * receives, which is the only part of the rate the discipline decides.
 *
 * PHASE and DEP share this factor exactly: an internal phase transition is
 * driven by the same server share as a completion, which is why a job under PS
 * advances through its phases more slowly when the station is busy. Only the
 * MATRIX differs -- D0(k,kdest) for a phase advance, D1 for a completion.
 */
template <class T>
T service_share(const NetworkStruct<T>& sn, std::size_t ist, const Marginal<T>& m,
                std::size_t cls, double ni, double S) {
    const std::size_t R = sn.nclasses;
    const SchedStrategy sched = sn.stations[ist - 1].sched;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    // The *PRIO variants behave as their base discipline while every job fits
    // in a server -- with n <= c nobody is waiting, so precedence is moot. Once
    // the station saturates, only the MOST URGENT group present is served, and
    // the sharing is computed among that group alone. `prio_pop` is the single
    // definition of that group, shared with the load-dependent lookup.
    const PrioPop<T> p = prio_pop(sn, ist, m, cls, ni, S);
    if (!p.served) return zero;  // a lower-priority class gets no service at all
    const std::vector<T>& nir = p.nir;
    const double nieff = p.ni;

    if (sched == SchedStrategy::PS || sched == SchedStrategy::LPS ||
        sched == SchedStrategy::PSPRIO)
        return nieff > 0 ? num_traits<T>::from_double(std::min(nieff, S) / nieff) : zero;
    if (sched == SchedStrategy::DPS || sched == SchedStrategy::GPS ||
        sched == SchedStrategy::DPSPRIO || sched == SchedStrategy::GPSPRIO) {
        if (S > 1)
            throw UnsupportedError(
                "state_events: multi-server DPS/GPS stations are not supported");
        const bool dps = sched == SchedStrategy::DPS || sched == SchedStrategy::DPSPRIO;
        const std::vector<T>& w = sn.stations[ist - 1].schedparam;
        T wsum = zero;
        for (std::size_t r = 0; r < R; ++r) wsum += w[r];
        T denom = zero;
        for (std::size_t r = 0; r < R; ++r) {
            // DPS weights by the job COUNTS; GPS shares between the classes
            // PRESENT, so each contributes at most one however many it holds.
            const T nr = dps ? nir[r] : (num_traits<T>::to_double(nir[r]) > 0 ? one : zero);
            denom += T(w[r] / wsum * nr);
        }
        const T nc = nir[cls - 1];
        if (num_traits<T>::to_double(denom) == 0 || num_traits<T>::to_double(nc) == 0)
            return zero;
        const T sh = T((w[cls - 1] / wsum) / denom);
        return dps ? sh : T(sh / nc);
    }
    // INF and the queueing disciplines: a job in service holds a whole server.
    return one;
}

/**
 * Port of the PHASE branch of `State.afterEventStation`: service advances a
 * phase WITHOUT completing.
 *
 * The rate is D0(k, kdest), the off-diagonal of the hidden generator, times the
 * same server share a completion gets. Keeping PHASE and DEP on one share is
 * what makes a phase-type service slow down consistently under contention; a
 * phase advance at full speed under PS would shorten the effective service.
 */
template <class T>
EventOutcome<T> after_event_station_phase(const NetworkStruct<T>& sn, std::size_t ind,
                                          const std::vector<T>& inspace, std::size_t cls) {
    const std::size_t R = sn.nclasses;
    const std::size_t ist = sn.nodes[ind - 1].station;
    const T one = num_traits<T>::from_int(1);
    EventOutcome<T> out;
    if (ist == 0) throw InputError("after_event_station_phase: node is not a station");
    const RowLayout<T> L = row_layout(sn, ind, inspace.size());
    const double S = sn.stations[ist - 1].nservers;

    std::vector<std::size_t> ph(R, 1), shift(R, 0);
    for (std::size_t r = 0; r < R; ++r) {
        ph[r] = L.K[r];
        shift[r] = L.Ks[r];
    }
    const Marginal<T> m = to_marginal(sn, ist, inspace, ph, shift, L.nvar);
    if (num_traits<T>::to_double(m.nir[cls - 1]) <= 0) return out;

    double ni = 0;
    for (std::size_t r = 0; r < R; ++r) ni += num_traits<T>::to_double(m.nir[r]);
    // Unlike DEP, the PHASE branch takes the UNMASKED population for both
    // factors even at a saturated *PRIO station (afterEventStation.m:1688): the
    // masking there is a property of the completion rate, not of the lookup, and
    // `service_share` already zeroes a non-urgent class's advance.
    const T lld = T(lld_factor(sn, ist, ni) * cd_factor(sn, ist, m.nir, cls));
    const T share = service_share(sn, ist, m, cls, ni, S);
    const lang::Distrib<T>& d = sn.service[ist - 1][cls - 1];
    if (d.D0.rows() != L.K[cls - 1]) return out;

    for (std::size_t k = 0; k < L.K[cls - 1]; ++k) {
        if (num_traits<T>::to_double(inspace[L.bufw + L.Ks[cls - 1] + k]) <= 0) continue;
        for (std::size_t kd = 0; kd < L.K[cls - 1]; ++kd) {
            if (kd == k) continue;  // the diagonal is the exit rate, not a move
            std::vector<T> row = inspace;
            row[L.bufw + L.Ks[cls - 1] + k] -= one;
            row[L.bufw + L.Ks[cls - 1] + kd] += one;
            out.space.push_back(row);
            out.rate.push_back(T(lld * d.D0(k, kd) * m.kir[cls - 1][k] * share));
            out.prob.push_back(one);
        }
    }
    return out;
}

/**
 * Port of the RENEGE branch: a WAITING class-`cls` job abandons the queue.
 *
 * Patience is exponential, so every waiting job abandons at the same rate and
 * the aggregate out of this state is (waiting count) * mu. Which job leaves is
 * therefore immaterial -- waiting jobs are exchangeable under memoryless
 * patience -- so the reference removes the first tagged slot and re-pads a zero
 * on the left, keeping the buffer in the right-aligned form the arrival handler
 * expects.
 */
template <class T>
EventOutcome<T> after_event_station_renege(const NetworkStruct<T>& sn, std::size_t ind,
                                           const std::vector<T>& inspace, std::size_t cls,
                                           const T& impatience_mu) {
    const std::size_t R = sn.nclasses;
    const std::size_t ist = sn.nodes[ind - 1].station;
    EventOutcome<T> out;
    if (ist == 0) throw InputError("after_event_station_renege: node is not a station");
    const RowLayout<T> L = row_layout(sn, ind, inspace.size());

    std::vector<std::size_t> ph(R, 1), shift(R, 0);
    for (std::size_t r = 0; r < R; ++r) {
        ph[r] = L.K[r];
        shift[r] = L.Ks[r];
    }
    const Marginal<T> m = to_marginal(sn, ist, inspace, ph, shift, L.nvar);
    const double waiting = num_traits<T>::to_double(m.nir[cls - 1]) -
                           num_traits<T>::to_double(m.sir[cls - 1]);
    if (waiting <= 0) return out;

    const T tag = num_traits<T>::from_int(static_cast<long>(cls));
    std::size_t slot = L.bufw;
    for (std::size_t b = 0; b < L.bufw; ++b)
        if (inspace[b] == tag) { slot = b; break; }
    if (slot == L.bufw) return out;

    std::vector<T> row = inspace;
    for (std::size_t b = slot; b > 0; --b) row[b] = row[b - 1];
    row[0] = num_traits<T>::from_int(0);
    out.space.push_back(row);
    out.rate.push_back(T(num_traits<T>::from_double(waiting) * impatience_mu));
    out.prob.push_back(num_traits<T>::from_int(1));
    return out;
}

/**
 * Port of the RETRY branch: an ORBITING class-`cls` job retries entry.
 *
 * The retry succeeds only when a server is free; otherwise the job stays in
 * orbit and the event is not generated at all, which is exactly what
 * distinguishes a retrial queue from a queue whose buffer is called an orbit.
 *
 * @param constant_policy CONSTANT retrial: one controller retries for the whole
 *        orbit, so the rate does NOT scale with the orbit size. Under the
 *        default LINEAR policy every orbiting job carries its own timer and the
 *        aggregate rate is (orbit size) * mu.
 * @param sn the refreshed network struct
 * @param ind index of the station the event fires at
 * @param inspace the state the event is applied to
 * @param cls class of the retrying job
 * @param retrial_mu retrial rate of that class
 */
template <class T>
EventOutcome<T> after_event_station_retry(const NetworkStruct<T>& sn, std::size_t ind,
                                          const std::vector<T>& inspace, std::size_t cls,
                                          const T& retrial_mu, bool constant_policy = false) {
    const std::size_t R = sn.nclasses;
    const std::size_t ist = sn.nodes[ind - 1].station;
    const T one = num_traits<T>::from_int(1);
    EventOutcome<T> out;
    if (ist == 0) throw InputError("after_event_station_retry: node is not a station");
    const RowLayout<T> L = row_layout(sn, ind, inspace.size());
    const double S = sn.stations[ist - 1].nservers;

    std::vector<std::size_t> ph(R, 1), shift(R, 0);
    for (std::size_t r = 0; r < R; ++r) {
        ph[r] = L.K[r];
        shift[r] = L.Ks[r];
    }
    const Marginal<T> m = to_marginal(sn, ist, inspace, ph, shift, L.nvar);
    const double orbit = num_traits<T>::to_double(m.nir[cls - 1]) -
                         num_traits<T>::to_double(m.sir[cls - 1]);
    double occ = 0;
    for (std::size_t j = 0; j < L.srvw; ++j) occ += num_traits<T>::to_double(inspace[L.bufw + j]);
    if (orbit <= 0 || occ >= S) return out;

    const T tag = num_traits<T>::from_int(static_cast<long>(cls));
    std::size_t slot = L.bufw;
    for (std::size_t b = 0; b < L.bufw; ++b)
        if (inspace[b] == tag) { slot = b; break; }
    if (slot == L.bufw) return out;

    const std::vector<T> pentry = entry_phase_dist(sn, ist, cls);
    const T agg = constant_policy ? retrial_mu
                                  : T(num_traits<T>::from_double(orbit) * retrial_mu);
    for (std::size_t ke = 0; ke < L.K[cls - 1]; ++ke) {
        if (num_traits<T>::to_double(pentry[ke]) <= 0) continue;
        std::vector<T> row = inspace;
        for (std::size_t b = slot; b > 0; --b) row[b] = row[b - 1];
        row[0] = num_traits<T>::from_int(0);
        row[L.bufw + L.Ks[cls - 1] + ke] += one;
        out.space.push_back(row);
        out.rate.push_back(T(agg * pentry[ke]));
        out.prob.push_back(one);
        // A successful retry is the only way into the server at a retrial
        // station (its departures never promote from the orbit), so it carries
        // the START the startRate == TN + preemptRate identity needs.
        tag_last(out, cls, 0);
    }
    pad_tags(out);
    return out;
}

/**
 * Port of the FAILURE and REPAIR branches: the server goes down, or comes back.
 *
 * Only the STATUS column moves, which is the trailing local variable. Jobs in
 * service are NOT lost: service is memoryless here, so an interrupted job
 * resumes on repair with no state to remember. The passive half of the
 * synchronization is LOCAL, so no job moves anywhere in the network either.
 *
 * @param up   true for REPAIR (0 -> 1), false for FAILURE (1 -> 0)
 * @param mu   breakdownMu for a failure, repairMu for a repair
 * @param sn the refreshed network struct
 * @param ind index of the station the event fires at
 * @param inspace the state the event is applied to
 */
template <class T>
EventOutcome<T> after_event_station_breakdown(const NetworkStruct<T>& sn, std::size_t ind,
                                              const std::vector<T>& inspace, bool up,
                                              const T& mu) {
    EventOutcome<T> out;
    if (inspace.empty()) return out;
    (void)sn;
    (void)ind;
    const double status = num_traits<T>::to_double(inspace.back());
    // A failure needs an UP server and a repair a DOWN one; anything else is
    // not an admissible transition and must produce no edge at all.
    if ((up && status != 0) || (!up && status != 1)) return out;
    std::vector<T> row = inspace;
    row.back() = num_traits<T>::from_int(up ? 1 : 0);
    out.space.push_back(row);
    out.rate.push_back(mu);
    out.prob.push_back(num_traits<T>::from_int(1));
    return out;
}

/**
 * Port of `State.afterEventFork`: an event at a STATEFUL Fork node.
 *
 * The fork's state is a plain per-class count of PARENT jobs momentarily held
 * between their arrival and the firing. An arrival buffers one; a DEPARTURE DOES
 * NOT EXIST, because the multi-branch emission is atomic across several nodes and
 * cannot be decomposed into a departure here plus an arrival there --
 * `refresh_sync` emits no DEP sync for a Fork and `after_fj_event` fires instead.
 *
 * The arrival rate is left UNSET (the reference writes -1) because an ARV is the
 * PASSIVE half of a synchronization: the rate belongs to the active departure
 * upstream.
 */
template <class T>
EventOutcome<T> after_event_fork(const NetworkStruct<T>& sn, std::size_t ind,
                                 const std::vector<T>& inspace, EventType event,
                                 std::size_t cls) {
    const std::size_t R = sn.nclasses;
    EventOutcome<T> out;
    if (event != EventType::ARV) return out;
    if (inspace.size() < R)
        throw InputError("after_event_fork: the Fork state row at node '" + sn.nodes[ind - 1].name +
                         "' is narrower than the class set");
    std::vector<T> row = inspace;
    // The counts are the LAST R columns, which is where `from_marginal_node`
    // puts them and where `after_fj_event` reads them.
    row[row.size() - R + cls - 1] += num_traits<T>::from_int(1);
    out.space.push_back(row);
    out.rate.push_back(num_traits<T>::from_int(-1));
    out.prob.push_back(num_traits<T>::from_int(1));
    return out;
}

/**
 * Port of `State.afterEventJoin`: an event at a Join node of an FJ-augmented
 * struct.
 *
 * The join's state is a plain per-class count vector of BUFFERED jobs, so it
 * bypasses the buffer/server/local slicing every other station takes -- a join
 * performs no service, it performs a rendezvous.
 *
 * ARV buffers the arriving job or sibling. DEP in an ORIGINAL class r fires when
 * either a plain (never-forked) class-r job is buffered, or some tag has its full
 * required sibling multiset present; the firing consumes the siblings of the
 * LOWEST complete tag, which is the same canonical choice the fork's allocation
 * makes and is what keeps the two in step. DEP in an AUXILIARY class is refused
 * outright: a sibling never departs on its own, it is consumed by the parent's
 * firing, and letting it depart would release a job the fork never emitted.
 */
template <class T>
EventOutcome<T> after_event_join(const NetworkStruct<T>& sn, std::size_t ind,
                                 const std::vector<T>& inspace, EventType event,
                                 std::size_t cls) {
    const std::size_t R = sn.nclasses;
    const T one = num_traits<T>::from_int(1);
    EventOutcome<T> out;
    if (inspace.size() < R)
        throw InputError("after_event_join: the Join state row at node '" + sn.nodes[ind - 1].name +
                         "' is narrower than the class set");
    const std::size_t off = inspace.size() - R;

    if (event == EventType::ARV) {
        std::vector<T> row = inspace;
        row[off + cls - 1] += one;
        out.space.push_back(row);
        out.rate.push_back(num_traits<T>::from_int(-1));
        out.prob.push_back(one);
        return out;
    }
    if (event != EventType::DEP) return out;

    const typename std::map<std::size_t, FjJoinParam>::const_iterator jit =
        sn.fjjoinparam.find(ind);
    const FjJoinParam* fjp = jit == sn.fjjoinparam.end() ? 0 : &jit->second;

    if (fjp) {
        for (std::size_t x = 0; x < fjp->origclasses.size(); ++x) {
            const std::map<std::size_t, std::vector<std::vector<std::size_t>>>::const_iterator ait =
                fjp->auxmatrix.find(fjp->origclasses[x]);
            if (ait == fjp->auxmatrix.end()) continue;
            for (std::size_t b = 0; b < ait->second.size(); ++b)
                for (std::size_t t = 0; t < ait->second[b].size(); ++t)
                    if (ait->second[b][t] == cls) return out;  // an auxiliary class
        }
    }

    const T imm = num_traits<T>::from_double(lang::GlobalConstants::Immediate);
    if (num_traits<T>::to_double(inspace[off + cls - 1]) > 0) {
        // A plain job that never went through the fork: it passes straight
        // through, since a join is a rendezvous only for siblings.
        std::vector<T> row = inspace;
        row[off + cls - 1] -= one;
        out.space.push_back(row);
        out.rate.push_back(imm);
        out.prob.push_back(one);
        return out;
    }
    if (!fjp) return out;
    const std::map<std::size_t, std::vector<std::vector<std::size_t>>>::const_iterator ait =
        fjp->auxmatrix.find(cls);
    const std::map<std::size_t, std::vector<std::size_t>>::const_iterator rit =
        fjp->required.find(cls);
    if (ait == fjp->auxmatrix.end() || rit == fjp->required.end()) return out;
    const std::vector<std::vector<std::size_t>>& aux = ait->second;
    const std::vector<std::size_t>& req = rit->second;
    if (aux.empty()) return out;
    const std::size_t B = aux.size(), Tt = aux[0].size();
    for (std::size_t t = 0; t < Tt; ++t) {
        bool complete = true;
        for (std::size_t b = 0; b < B && complete; ++b) {
            const std::size_t a = aux[b][t];
            const double need = b < req.size() ? static_cast<double>(req[b]) : 1.0;
            if (num_traits<T>::to_double(inspace[off + a - 1]) < need) complete = false;
        }
        if (!complete) continue;
        std::vector<T> row = inspace;
        for (std::size_t b = 0; b < B; ++b) {
            const std::size_t a = aux[b][t];
            const long need = b < req.size() ? static_cast<long>(req[b]) : 1;
            row[off + a - 1] -= num_traits<T>::from_int(need);
        }
        out.space.push_back(row);
        out.rate.push_back(imm);
        out.prob.push_back(one);
        return out;  // the LOWEST complete tag only; the rest are permutations
    }
    return out;
}

/**
 * Port of `State.afterEventStation`'s dispatch: the successors of one event at
 * one station.
 *
 * The event-specific rates that are not derivable from `sn` alone -- patience,
 * retrial and breakdown -- are passed in, because the reference reads them from
 * fields (`impatienceMu`, `retrialMu`, `breakdownMu`) that this port has not
 * yet grown. Every rate that IS derivable is computed from the struct.
 */
template <class T>
void rr_advance_row(const NetworkStruct<T>& sn, std::size_t ind, std::size_t cls,
                    std::vector<std::vector<T>>& rows);

/**
 * Port of `State.afterEventRouter`: a Router holds a job for the instant it
 * takes to decide where it goes.
 *
 * The row is [per-class counts | local vars], with no buffer and no phase: a
 * Router serves nothing, so there is no service state to carry. An ARRIVAL adds
 * the job at an UNSPECIFIED rate (-1), which is the reference's marker for a
 * passive half whose rate the active half sets; a DEPARTURE removes it at the
 * Immediate rate and advances the dispatch pointer, so the router never holds a
 * job for a positive length of time.
 */
template <class T>
EventOutcome<T> after_event_router(const NetworkStruct<T>& sn, std::size_t ind,
                                   const std::vector<T>& inspace, EventType event,
                                   std::size_t cls) {
    EventOutcome<T> out;
    const std::size_t R = sn.nclasses;
    if (inspace.size() < R) return out;
    const T one = num_traits<T>::from_int(1);
    if (event == EventType::ARV) {
        std::vector<T> row = inspace;
        row[cls - 1] = T(row[cls - 1] + one);
        out.space.push_back(row);
        // Passive action: the rate is the active half's, not this node's.
        out.rate.push_back(num_traits<T>::from_int(-1));
        out.prob.push_back(one);
        return out;
    }
    if (event == EventType::DEP) {
        if (!(num_traits<T>::to_double(inspace[cls - 1]) > 0)) return out;
        std::vector<T> row = inspace;
        row[cls - 1] = T(row[cls - 1] - one);
        out.space.push_back(row);
        out.rate.push_back(num_traits<T>::from_double(lang::GlobalConstants::Immediate));
        out.prob.push_back(one);
        rr_advance_row(sn, ind, cls, out.space);
        return out;
    }
    return out;
}

/**
 * Advance the round-robin dispatch pointer of (IND, CLS) in every successor row.
 *
 * The local-variable block is the TAIL of a state row, so the pointer is located
 * from the right; `rr_var_slot` gives its 1-based index inside that block. A
 * no-op wherever the pair does not dispatch round-robin.
 */
template <class T>
void rr_advance_row(const NetworkStruct<T>& sn, std::size_t ind, std::size_t cls,
                    std::vector<std::vector<T>>& rows) {
    if (sn.rr_var_slot(ind, cls) == 0) return;
    const std::size_t w = sn.nvars_of(ind);
    if (w == 0) return;
    for (std::size_t i = 0; i < rows.size(); ++i) {
        if (rows[i].size() < w) continue;
        std::vector<T> var(rows[i].end() - w, rows[i].end());
        sn.rr_advance(ind, cls, var);
        std::copy(var.begin(), var.end(), rows[i].end() - w);
    }
}

template <class T>
EventOutcome<T> after_event_station(const NetworkStruct<T>& sn, std::size_t ind,
                                    const std::vector<T>& inspace, EventType event,
                                    std::size_t cls, bool no_promote = false,
                                    const T& aux_rate = num_traits<T>::from_int(0)) {
    switch (event) {
        case EventType::ARV:
            // A signal class arriving is not an arrival at all: it removes
            // jobs and is annihilated, so it never reaches the scheduling
            // branches. A REPLY signal is the exception -- it completes a
            // synchronous call and then joins as an ordinary job.
            if (sn.issignal.size() >= cls && sn.issignal[cls - 1]) {
                if (!(sn.signaltype.size() >= cls &&
                      sn.signaltype[cls - 1] == lang::SignalType::REPLY))
                    return after_event_station_signal(sn, ind, inspace, cls);
                // A REPLY takes the reply path only at a station that actually
                // holds a block for it; elsewhere it is a plain job class and
                // falls through to ordinary arrival handling.
                if (reply_block_info(sn, ind).width > 0)
                    return after_event_station_reply(sn, ind, inspace, cls);
            }
            return after_event_station_arv(sn, ind, inspace, cls);
        case EventType::DEP: {
            EventOutcome<T> out = after_event_station_dep(sn, ind, inspace, cls, no_promote);
            // ROUND-ROBIN DISPATCH advances on every completion, which is what
            // makes the next destination deterministic; the generator then reads
            // the pointer OUT OF THIS SUCCESSOR to pick the link. Without the
            // advance the pointer is a frozen coordinate and every job takes the
            // same link, which is random routing with the wrong support rather
            // than round robin. Reference: `afterEventStation.m:632-658`.
            rr_advance_row(sn, ind, cls, out.space);
            // TRUE BAS, the departure half. When the marker is already set the
            // front job has COMPLETED and is being held, so this DEP is not a
            // service completion at all -- it is the instant transfer of that
            // held job downstream, which the generator only offers when the
            // destination has room. It therefore fires at the Immediate rate and
            // CLEARS the marker; the complementary become-blocked edge (0 -> 1) is
            // added by the generator, the only place that can see the
            // destination's occupancy.
            //
            // Gate on `isbasblocking`, not on this station's own drop rule: under
            // the destination declaration form the rule is not here.
            if (ind <= sn.isbasblocking.size() && sn.isbasblocking[ind - 1] &&
                !inspace.empty() && num_traits<T>::to_double(inspace.back()) == 1 &&
                !out.space.empty()) {
                // 1e7 VERBATIM, not GlobalConstants::Immediate (1e8). The rate is
                // large but FINITE, so the blocked states keep a proportional
                // share of the stationary mass and the analyzer's queue-length
                // shift reads it; using 1e8 would divide that share by ten and
                // move every reported queue length. The reference hardcodes 1e7.
                const T imm = num_traits<T>::from_double(1e7);
                for (std::size_t i = 0; i < out.space.size(); ++i) {
                    out.space[i].back() = num_traits<T>::from_int(0);
                    out.rate[i] = imm;
                }
            }
            return out;
        }
        case EventType::PHASE:
            return after_event_station_phase(sn, ind, inspace, cls);
        case EventType::RENEGE:
            return after_event_station_renege(sn, ind, inspace, cls, aux_rate);
        case EventType::RETRY:
            return after_event_station_retry(sn, ind, inspace, cls, aux_rate);
        case EventType::SWITCH:
            return after_event_station_switch(sn, ind, inspace, cls);
        case EventType::FAILURE:
            return after_event_station_breakdown(sn, ind, inspace, false, aux_rate);
        case EventType::REPAIR:
            return after_event_station_breakdown(sn, ind, inspace, true, aux_rate);
        case EventType::LOCAL:
            return EventOutcome<T>();  // a dummy event moves nothing
        default:
            throw UnsupportedError(std::string("after_event_station: the ") +
                                   lang::event_to_text(event) +
                                   " event is not ported yet");
    }
}

/**
 * Port of `State.afterEventTransition`, the PHASE arm: one running server of
 * the given mode advances its firing phase (`cls` is interpreted as the MODE,
 * as in the reference). ENABLE and FIRE are global events handled by
 * `after_global_event`, so they return an empty outcome here, matching the
 * reference's no-op arms.
 *
 * RATE. The MATLAB body multiplies the phase-k move rate by BOTH
 * kir(:,mode,k) and nir(mode) (afterEventTransition.m:38-40), but nir is the
 * sum of kir over the phases, so the extra factor counts the running servers
 * twice; the JAR and the native python carry D0(k,kdest) * kir alone, and
 * this port follows them.
 *
 * The row layout is the one `after_global_event` slices:
 * [buf(nmodes) | srv(sum fK) | fired(nmodes) | var].
 */
template <class T>
EventOutcome<T> after_event_transition(const NetworkStruct<T>& sn, std::size_t ind,
                                       const std::vector<T>& inspace, EventType event,
                                       std::size_t mode) {
    EventOutcome<T> out;
    if (event != EventType::PHASE) return out;  // ENABLE / FIRE are global
    const typename std::map<std::size_t, TransitionParam<T>>::const_iterator it =
        sn.transparam.find(ind);
    if (it == sn.transparam.end())
        throw InputError("after_event_transition: node has no TransitionParam");
    const TransitionParam<T>& tp = it->second;
    if (mode == 0 || mode > tp.nmodes)
        throw InputError("after_event_transition: mode index is out of range");
    const T one = num_traits<T>::from_int(1);

    std::vector<std::size_t> fK(tp.nmodes, 1), fKs(tp.nmodes, 0);
    std::size_t tot = 0;
    for (std::size_t m = 0; m < tp.nmodes; ++m) {
        fK[m] = m < tp.firingphases.size() && tp.firingphases[m] > 0 ? tp.firingphases[m] : 1;
        fKs[m] = tot;
        tot += fK[m];
    }
    if (fK[mode - 1] <= 1) return out;  // a single phase has no internal move
    if (mode - 1 >= tp.firingproc.size() ||
        tp.firingproc[mode - 1].D0.rows() != fK[mode - 1])
        return out;
    const Matrix<T>& D0 = tp.firingproc[mode - 1].D0;

    for (std::size_t k = 0; k < fK[mode - 1]; ++k) {
        const std::size_t idx = tp.nmodes + fKs[mode - 1] + k;
        const T cnt = inspace[idx];
        if (!(num_traits<T>::to_double(cnt) > 0)) continue;
        for (std::size_t kd = 0; kd < fK[mode - 1]; ++kd) {
            if (kd == k) continue;  // the diagonal is the exit rate, not a move
            if (!(num_traits<T>::to_double(D0(k, kd)) > 0)) continue;
            std::vector<T> row = inspace;
            row[idx] -= one;
            row[tp.nmodes + fKs[mode - 1] + kd] += one;
            out.space.push_back(row);
            out.rate.push_back(T(D0(k, kd) * cnt));
            out.prob.push_back(one);
        }
    }
    return out;
}

/**
 * Port of `State.afterEvent`: the successors of one event at one NODE.
 *
 * The reference's body is mostly slicing -- it cuts `inspace` into buffer,
 * server and local-variable blocks and hands the pieces to the per-node-type
 * handler. This port slices inside each handler instead (`row_layout`), so what
 * remains here is the dispatch itself and the guards that precede it.
 *
 * A class the station does not accept short-circuits: `phases_of` is zero
 * there, and every downstream index into the server block would be out of
 * range. That guard is the reference's `K(class) == 0` test.
 *
 * `cls` IS A MODE, NOT A CLASS, on a Transition's PHASE action; see the guard.
 */
template <class T>
EventOutcome<T> after_event(const NetworkStruct<T>& sn, std::size_t ind,
                            const std::vector<T>& inspace, EventType event, std::size_t cls,
                            bool no_promote = false,
                            const T& aux_rate = num_traits<T>::from_int(0)) {
    if (ind == 0 || ind > sn.nodes.size())
        throw InputError("after_event: node index is out of range");
    const NodeDef& nd = sn.nodes[ind - 1];
    // THE `cls` SLOT IS NOT ALWAYS A CLASS. On a Transition's PHASE action it
    // carries the MODE, which `refresh_sync` puts there deliberately ("one
    // server phase-change action per MODE, not per class") and
    // `after_event_transition` reads back as such. Bounding it by `nclasses`
    // refused every mode past the class count: spn_basic_closed is one class
    // and three modes, so modes 2 and 3 raised "class index is out of range"
    // and no closed SPN with more modes than classes could be walked at all.
    // The mode bound belongs to the handler, which already applies it against
    // `tp.nmodes`, so only the non-Transition case is checked here.
    const bool cls_is_mode = (nd.nodetype == NodeType::Transition && event == EventType::PHASE);
    if (!cls_is_mode && (cls == 0 || cls > sn.nclasses))
        throw InputError("after_event: class index is out of range");

    // A Join of an FJ-augmented struct IS a station, but its state is a bare
    // per-class count vector: it performs a rendezvous, not a service, so it must
    // bypass the buffer/server slicing before `after_event_station` sees it.
    if (sn.isfjaugmented && nd.nodetype == NodeType::Join)
        return after_event_join(sn, ind, inspace, event, cls);

    if (nd.station != 0) {
        // A class with no service process at this station cannot be involved
        // in any event here.
        if (sn.phases_of(nd.station, cls) == 0) return EventOutcome<T>();
        return after_event_station(sn, ind, inspace, event, cls, no_promote, aux_rate);
    }
    if (!nd.stateful) return EventOutcome<T>();  // a stateless node holds nothing

    if (nd.nodetype == NodeType::Cache) return after_event_cache(sn, ind, inspace, event, cls);
    if (nd.nodetype == NodeType::Fork) return after_event_fork(sn, ind, inspace, event, cls);
    if (nd.nodetype == NodeType::Transition)
        return after_event_transition(sn, ind, inspace, event, cls);
    if (nd.nodetype == NodeType::Router) return after_event_router(sn, ind, inspace, event, cls);

    throw UnsupportedError("after_event: events at stateful non-station node '" + nd.name +
                           "' are not ported yet");
}

namespace signal_detail {

/** Merge duplicate destinations so the generator sees one entry per state. */
template <class T>
void merge_states(std::vector<std::vector<T>>& sp, std::vector<T>& pr) {
    std::vector<std::vector<T>> us;
    std::vector<T> up;
    for (std::size_t i = 0; i < sp.size(); ++i) {
        std::size_t at = us.size();
        for (std::size_t j = 0; j < us.size(); ++j)
            if (us[j] == sp[i]) { at = j; break; }
        if (at == us.size()) {
            us.push_back(sp[i]);
            up.push_back(pr[i]);
        } else {
            up[at] += pr[i];
        }
    }
    sp.swap(us);
    pr.swap(up);
}

}  // namespace signal_detail

/**
 * Port of `State.signalBatchPMF`: the batch size a negative signal removes.
 *
 * The pmf is CLIPPED at the eligible population: an oversized batch empties
 * the station rather than driving the queue negative, so the whole tail
 * P(B >= n) lumps onto "remove all n". That is the same clipping LDES applies
 * with min(B,n) and the tail term MAM uses.
 */
template <class T>
std::pair<std::vector<std::size_t>, std::vector<T>> signal_batch_pmf(
    const NetworkStruct<T>& sn, std::size_t cls, std::size_t ntot) {
    std::vector<std::size_t> kv;
    std::vector<T> kp;
    if (sn.signalremdist.size() < cls || sn.signalremdist[cls - 1].empty()) {
        kv.push_back(1);
        kp.push_back(num_traits<T>::from_int(1));
        return std::make_pair(kv, kp);
    }
    const std::vector<T>& d = sn.signalremdist[cls - 1];
    T head_sum = num_traits<T>::from_int(0);
    for (std::size_t b = 0; b < ntot; ++b) {
        const T p = b < d.size() ? d[b] : num_traits<T>::from_int(0);
        if (num_traits<T>::to_double(p) > 0) {
            kv.push_back(b);
            kp.push_back(p);
        }
        head_sum += p;
    }
    const double tail = 1.0 - num_traits<T>::to_double(head_sum);
    if (tail > 0) {
        kv.push_back(ntot);
        kp.push_back(num_traits<T>::from_double(tail));
    }
    if (kv.empty()) {
        kv.push_back(1);
        kp.push_back(num_traits<T>::from_int(1));
        return std::make_pair(kv, kp);
    }
    T tot = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < kp.size(); ++i) tot += kp[i];
    if (num_traits<T>::to_double(tot) > 0)
        for (std::size_t i = 0; i < kp.size(); ++i) kp[i] = T(kp[i] / tot);
    return std::make_pair(kv, kp);
}

/**
 * Port of `State.afterEventStationSignal`: a G-network signal arrives.
 *
 * A signal NEVER joins the station. It removes jobs already there and is
 * annihilated, so the event is passive throughout and the successors differ
 * only in which victims were taken.
 *
 * Victim selection has two tiers, and conflating them is the trap: FCFS and
 * LCFS rank by AGE, which only an ordered buffer records, so at a per-class
 * count buffer an age policy degenerates to a uniform draw. They also drain
 * the waiting line completely before touching a server, whereas RANDOM draws
 * uniformly across waiting and in-service jobs alike.
 */
template <class T>
EventOutcome<T> after_event_station_signal(const NetworkStruct<T>& sn, std::size_t ind,
                                           const std::vector<T>& inspace, std::size_t cls) {
    const std::size_t R = sn.nclasses;
    const std::size_t ist = sn.nodes[ind - 1].station;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T minus_one = num_traits<T>::from_int(-1);
    EventOutcome<T> out;
    if (ist == 0) throw InputError("after_event_station_signal: node is not a station");
    const RowLayout<T> L = row_layout(sn, ind, inspace.size());
    const SchedStrategy sched = sn.stations[ist - 1].sched;
    const double S = sn.stations[ist - 1].nservers;

    std::vector<std::size_t> ph(R, 1), shift(R, 0);
    for (std::size_t r = 0; r < R; ++r) {
        ph[r] = L.K[r];
        shift[r] = L.Ks[r];
    }
    const Marginal<T> m = to_marginal(sn, ist, inspace, ph, shift, L.nvar);

    // A CATASTROPHE empties the station outright, ignoring the batch pmf: a
    // catastrophe removes every job by definition.
    if (sn.signaltype.size() >= cls && sn.signaltype[cls - 1] == lang::SignalType::CATASTROPHE) {
        std::vector<T> row(L.bufw + L.srvw, zero);
        row.insert(row.end(), inspace.end() - L.nvar, inspace.end());
        out.space.push_back(row);
        out.rate.push_back(minus_one);
        out.prob.push_back(one);
        return out;
    }

    // Eligible victim classes: the declared target, or every non-signal class
    // for the classic untargeted negative customer.
    std::vector<std::size_t> tgt;
    const std::size_t declared = sn.signaltarget.size() >= cls ? sn.signaltarget[cls - 1] : 0;
    if (declared >= 1) {
        tgt.push_back(declared);
    } else {
        for (std::size_t r = 1; r <= R; ++r)
            if (sn.issignal.size() < r || !sn.issignal[r - 1]) tgt.push_back(r);
    }
    std::vector<std::size_t> elig;
    std::size_t ntot = 0;
    for (std::size_t i = 0; i < tgt.size(); ++i) {
        const double nr = num_traits<T>::to_double(m.nir[tgt[i] - 1]);
        if (nr > 0) {
            elig.push_back(tgt[i]);
            ntot += static_cast<std::size_t>(nr);
        }
    }
    if (elig.empty()) {
        // No victim: the signal simply vanishes, leaving the state unchanged.
        out.space.push_back(inspace);
        out.rate.push_back(minus_one);
        out.prob.push_back(one);
        return out;
    }

    const lang::RemovalPolicy policy = sn.signalrempolicy.size() >= cls
                                           ? sn.signalrempolicy[cls - 1]
                                           : lang::RemovalPolicy::RANDOM;
    const bool ordered = state_detail::buffer_is_class_tag(sched) ||
                         state_detail::buffer_is_tag_phase_pairs(sched);
    const bool paired = state_detail::buffer_is_tag_phase_pairs(sched);
    const bool counted = state_detail::buffer_is_per_class_count(sched);
    const std::pair<std::vector<std::size_t>, std::vector<T>> pmf =
        signal_batch_pmf(sn, cls, ntot);

    std::vector<std::vector<T>> acc;
    std::vector<T> accp;
    for (std::size_t ik = 0; ik < pmf.first.size(); ++ik) {
        if (num_traits<T>::to_double(pmf.second[ik]) <= 0) continue;
        std::vector<std::vector<T>> cur(1, inspace);
        std::vector<T> curp(1, one);
        // Remove one job at a time: sequential draws without replacement give
        // a uniform choice of the removed SUBSET.
        for (std::size_t step = 0; step < pmf.first[ik]; ++step) {
            std::vector<std::vector<T>> nxt;
            std::vector<T> nxtp;
            for (std::size_t row = 0; row < cur.size(); ++row) {
                std::vector<T> buf(cur[row].begin(), cur[row].begin() + L.bufw);
                std::vector<T> srv(cur[row].begin() + L.bufw,
                                   cur[row].begin() + L.bufw + L.srvw);
                const std::vector<T> var(cur[row].begin() + L.bufw + L.srvw, cur[row].end());

                // Waiting victims, as (position, class, multiplicity).
                std::vector<std::size_t> wpos, wcls, wwt;
                if (ordered) {
                    for (std::size_t b = 0; b < L.bufw; b += paired ? 2 : 1) {
                        const double v = num_traits<T>::to_double(buf[b]);
                        if (v <= 0) continue;
                        bool ok = false;
                        for (std::size_t e = 0; e < elig.size(); ++e)
                            if (elig[e] == static_cast<std::size_t>(v)) ok = true;
                        if (!ok) continue;
                        wpos.push_back(b);
                        wcls.push_back(static_cast<std::size_t>(v));
                        wwt.push_back(1);
                    }
                } else if (counted) {
                    for (std::size_t e = 0; e < elig.size(); ++e) {
                        const std::size_t r = elig[e];
                        if (r > L.bufw) continue;
                        const double v = num_traits<T>::to_double(buf[r - 1]);
                        if (v <= 0) continue;
                        wpos.push_back(r - 1);
                        wcls.push_back(r);
                        wwt.push_back(static_cast<std::size_t>(v));
                    }
                }
                // In-service victims, per (class, phase).
                std::vector<std::size_t> scls, sph, scnt;
                for (std::size_t e = 0; e < elig.size(); ++e) {
                    const std::size_t r = elig[e];
                    for (std::size_t p = 0; p < L.K[r - 1]; ++p) {
                        const double v = num_traits<T>::to_double(srv[L.Ks[r - 1] + p]);
                        if (v <= 0) continue;
                        scls.push_back(r);
                        sph.push_back(p);
                        scnt.push_back(static_cast<std::size_t>(v));
                    }
                }
                std::size_t nwait = 0, nsrv = 0;
                for (std::size_t i = 0; i < wwt.size(); ++i) nwait += wwt[i];
                for (std::size_t i = 0; i < scnt.size(); ++i) nsrv += scnt[i];
                if (nwait == 0 && nsrv == 0) {
                    // Already drained: nothing left for this step to remove.
                    nxt.push_back(cur[row]);
                    nxtp.push_back(curp[row]);
                    continue;
                }

                std::vector<std::vector<T>> sp;
                std::vector<T> pr;
                const bool age_ordered =
                    ordered && (policy == lang::RemovalPolicy::FCFS ||
                                policy == lang::RemovalPolicy::LCFS);
                if (age_ordered && nwait > 0) {
                    // The head of line is the LAST occupied slot, the newest
                    // arrival the first: the buffer is right-aligned.
                    std::size_t pick = 0;
                    for (std::size_t i = 0; i < wpos.size(); ++i)
                        if (policy == lang::RemovalPolicy::FCFS ? wpos[i] > wpos[pick]
                                                                : wpos[i] < wpos[pick])
                            pick = i;
                    std::vector<T> b2 = buf;
                    b2.erase(b2.begin() + wpos[pick], b2.begin() + wpos[pick] + (paired ? 2 : 1));
                    b2.insert(b2.begin(), paired ? 2 : 1, zero);
                    std::vector<T> nr = b2;
                    nr.insert(nr.end(), srv.begin(), srv.end());
                    nr.insert(nr.end(), var.begin(), var.end());
                    sp.push_back(nr);
                    pr.push_back(one);
                } else {
                    // RANDOM draws over everything present; FCFS/LCFS at a
                    // count buffer drain the waiting line first, and only reach
                    // the servers once it is empty.
                    const std::size_t total = policy == lang::RemovalPolicy::RANDOM
                                                  ? nwait + nsrv
                                                  : (nwait > 0 ? nwait : nsrv);
                    for (std::size_t i = 0; i < wpos.size(); ++i) {
                        std::vector<T> b2 = buf;
                        if (counted) {
                            b2[wcls[i] - 1] -= one;
                        } else {
                            b2.erase(b2.begin() + wpos[i],
                                     b2.begin() + wpos[i] + (paired ? 2 : 1));
                            b2.insert(b2.begin(), paired ? 2 : 1, zero);
                        }
                        std::vector<T> nr = b2;
                        nr.insert(nr.end(), srv.begin(), srv.end());
                        nr.insert(nr.end(), var.begin(), var.end());
                        sp.push_back(nr);
                        pr.push_back(num_traits<T>::from_double(static_cast<double>(wwt[i]) /
                                                                static_cast<double>(total)));
                    }
                    if (policy == lang::RemovalPolicy::RANDOM || nwait == 0) {
                        for (std::size_t i = 0; i < scls.size(); ++i) {
                            std::vector<T> b2 = buf, s2 = srv;
                            s2[L.Ks[scls[i] - 1] + sph[i]] -= one;
                            // The freed server pulls in the head of line, where
                            // the station keeps one at all.
                            double occ = 0;
                            for (std::size_t j = 0; j < s2.size(); ++j)
                                occ += num_traits<T>::to_double(s2[j]);
                            if (L.bufw > 0 && occ < S) {
                                if (ordered) {
                                    std::size_t hp = L.bufw;
                                    for (std::size_t b = 0; b < L.bufw; b += paired ? 2 : 1)
                                        if (num_traits<T>::to_double(b2[b]) > 0) hp = b;
                                    if (hp != L.bufw) {
                                        const std::size_t pc = static_cast<std::size_t>(
                                            num_traits<T>::to_double(b2[hp]));
                                        std::size_t pp = 0;
                                        if (paired) {
                                            const double v =
                                                num_traits<T>::to_double(b2[hp + 1]);
                                            pp = v >= 1 ? static_cast<std::size_t>(v) - 1 : 0;
                                        }
                                        b2.erase(b2.begin() + hp,
                                                 b2.begin() + hp + (paired ? 2 : 1));
                                        b2.insert(b2.begin(), paired ? 2 : 1, zero);
                                        s2[L.Ks[pc - 1] + pp] += one;
                                    }
                                } else if (counted) {
                                    // A count buffer carries no order, so the
                                    // lowest-indexed waiting class is promoted
                                    // to keep the map single-valued; the actual
                                    // service order is resolved by the rates.
                                    for (std::size_t r = 1; r <= R && r <= L.bufw; ++r)
                                        if (num_traits<T>::to_double(b2[r - 1]) > 0) {
                                            b2[r - 1] -= one;
                                            s2[L.Ks[r - 1]] += one;
                                            break;
                                        }
                                }
                            }
                            std::vector<T> nr = b2;
                            nr.insert(nr.end(), s2.begin(), s2.end());
                            nr.insert(nr.end(), var.begin(), var.end());
                            sp.push_back(nr);
                            pr.push_back(num_traits<T>::from_double(
                                static_cast<double>(scnt[i]) / static_cast<double>(total)));
                        }
                    }
                }
                signal_detail::merge_states(sp, pr);
                for (std::size_t i = 0; i < sp.size(); ++i) {
                    nxt.push_back(sp[i]);
                    nxtp.push_back(T(curp[row] * pr[i]));
                }
            }
            signal_detail::merge_states(nxt, nxtp);
            cur.swap(nxt);
            curp.swap(nxtp);
        }
        for (std::size_t i = 0; i < cur.size(); ++i) {
            acc.push_back(cur[i]);
            accp.push_back(T(pmf.second[ik] * curp[i]));
        }
    }
    signal_detail::merge_states(acc, accp);
    out.space.swap(acc);
    out.prob.swap(accp);
    out.rate.assign(out.space.size(), minus_one);
    return out;
}

/**
 * Port of `State.passAndSwap`: the transition a service completion triggers at
 * a pass-and-swap station (Dorsman and Gardner 2024, Sect. 2.3).
 *
 * The completing job scans FORWARD from its own position for the first job it
 * may swap with per the graph G, takes that job's place and ejects it; the
 * ejected job repeats the scan. The chain ends at a job with no swappable
 * successor, and THAT job departs -- which is why the departing class is in
 * general not the class whose service completed.
 *
 * @param c 0-based list of 1-based class indices, oldest first
 * @param p 0-based position whose service token completed
 * @param G class-compatibility graph, G[a][b] true when a may swap with b
 * @return (the list after the transition, the 1-based departing class)
 */
template <class T>
std::pair<std::vector<std::size_t>, std::size_t> pass_and_swap(
    const std::vector<std::size_t>& c, std::size_t p,
    const std::vector<std::vector<bool>>& G) {
    const std::size_t n = c.size();
    if (p >= n) throw InputError("pass_and_swap: position is out of range for the state");
    std::vector<std::size_t> chain(1, p);
    std::size_t moving = c[p], cur = p;
    for (;;) {
        std::size_t q = n;
        for (std::size_t j = cur + 1; j < n; ++j)
            if (moving - 1 < G.size() && c[j] - 1 < G[moving - 1].size() &&
                G[moving - 1][c[j] - 1]) {
                q = j;
                break;
            }
        if (q == n) break;  // no swappable successor: this job departs
        chain.push_back(q);
        moving = c[q];
        cur = q;
    }
    const std::size_t dep = c[chain.back()];
    // Shift classes one step along the chain; the last is overwritten because
    // it departed, and the head-of-chain slot is then removed.
    std::vector<std::size_t> cnew = c;
    for (std::size_t i = 0; i + 1 < chain.size(); ++i) cnew[chain[i + 1]] = c[chain[i]];
    cnew.erase(cnew.begin() + chain[0]);
    return std::make_pair(cnew, dep);
}

/**
 * Port of `State.afterEventStationPAS`: events at a pass-and-swap station.
 *
 * The state here is NOT the [buffer | server] split every other discipline
 * uses: it is the full ordered list of class indices, left-aligned and zero
 * padded, with no server block at all. Service is governed by the rate
 * function mu(c) rather than by a per-class rate, so the whole notion of "in
 * service" is replaced by a token at each position.
 */
/**
 * Per-position service rate increments Delta_mu(c1..cp) = mu(c1..cp) - mu(c1..c_{p-1}).
 */
template <class T, class F>
inline std::vector<double> pas_increments(const F& mu_fun, const std::vector<std::size_t>& c) {
    std::vector<double> inc(c.size(), 0.0);
    double mu_prev = 0.0;
    for (std::size_t p = 0; p < c.size(); ++p) {
        const std::vector<std::size_t> prefix(c.begin(), c.begin() + p + 1);
        const double mu_cur = num_traits<T>::to_double(mu_fun(prefix));
        inc[p] = mu_cur - mu_prev;
        mu_prev = mu_cur;
    }
    return inc;
}

/**
 * Tag the successor just appended to OUT with the PAS positions that started
 * service on it: those of CNEW that are served (Delta_mu > 0) and were not
 * served in COLD. Under a swap the tag follows the POSITION rather than the job
 * identity, since pass-and-swap redefines which job holds a position.
 */
template <class T, class F>
inline void pas_tag_started(EventOutcome<T>& out, const F& mu_fun,
                            const std::vector<std::size_t>& cold,
                            const std::vector<std::size_t>& cnew) {
    if (!mu_fun || out.space.empty()) return;
    const std::vector<double> inc_new = pas_increments<T>(mu_fun, cnew);
    const std::vector<double> inc_old = pas_increments<T>(mu_fun, cold);
    for (std::size_t p = 0; p < inc_new.size(); ++p) {
        if (inc_new[p] <= 0) continue;
        if (p < inc_old.size() && inc_old[p] > 0) continue;  // already served
        tag_last(out, cnew[p], 0);
    }
}

template <class T>
EventOutcome<T> after_event_station_pas(const NetworkStruct<T>& sn, std::size_t ind,
                                        const std::vector<T>& inspace, EventType event,
                                        std::size_t cls) {
    const std::size_t ist = sn.nodes[ind - 1].station;
    const T one = num_traits<T>::from_int(1), zero = num_traits<T>::from_int(0);
    EventOutcome<T> out;
    if (ist == 0) throw InputError("after_event_station_pas: node is not a station");
    const typename std::map<std::size_t, typename NetworkStruct<T>::PasParam>::const_iterator it =
        sn.pasparam.find(ist);
    if (it == sn.pasparam.end() || !it->second.svc_rate_fun)
        throw InputError(
            "after_event_station_pas: the station has no service rate function mu(c); set one "
            "with set_pas");

    const std::size_t V = sn.nvars_of(ind);
    const std::size_t W = inspace.size() - V;
    std::vector<std::size_t> c;
    for (std::size_t i = 0; i < W; ++i) {
        const double v = num_traits<T>::to_double(inspace[i]);
        if (v > 0) c.push_back(static_cast<std::size_t>(v));
    }
    const std::vector<T> var(inspace.begin() + W, inspace.end());
    const double cap = sn.cap[ist - 1];

    if (event == EventType::ARV) {
        // The arrival joins at the BACK of the list, which is what records the
        // order the rate function is a function of.
        if (static_cast<double>(c.size()) >= cap) return out;  // full: lost
        std::vector<std::size_t> nc = c;
        nc.push_back(cls);
        // NO SLOT IN THE ENCODING IS A BLOCK, NOT A LOSS. The list occupies one
        // row position per job, so a row of width W holds W jobs; emitting a
        // successor here and letting `row.resize(W, zero)` cut the list back to
        // W would DESTROY the job that did not fit while reporting a state that
        // looks exactly like the pre-arrival one -- a customer silently gone
        // from a closed network, with no error anywhere. Returning no successor
        // disables the upstream departure instead, which is what every other
        // discipline's capacity filter does and what the enumerated CTMC space
        // does at its own width boundary.
        if (nc.size() > W) return out;
        std::vector<T> row;
        for (std::size_t i = 0; i < nc.size(); ++i)
            row.push_back(num_traits<T>::from_int(static_cast<long>(nc[i])));
        row.resize(W, zero);
        row.insert(row.end(), var.begin(), var.end());
        out.space.push_back(row);
        out.rate.push_back(num_traits<T>::from_int(-1));
        out.prob.push_back(one);
        // A PAS station has one clock for the whole station and no server to
        // hold, so "in service" means "at a position whose rate increment
        // Delta_mu is positive": a job starts exactly when a position goes from
        // a zero increment to a positive one. With mu(c) = 1 only the head is
        // served (M/M/1), with mu(c) = |c| every position is (M/M/inf), and the
        // rule reproduces both.
        pas_tag_started<T>(out, it->second.svc_rate_fun, c, nc);
        pad_tags(out);
        return out;
    }
    if (event != EventType::DEP) return out;  // PAS service is exponential: no PHASE

    // Each position holds a service token firing at the INCREMENT of mu over
    // the prefix ending there. Summing the increments telescopes to mu(c), so
    // the station's total service rate is exactly the rate function.
    T mu_prev = zero;
    for (std::size_t p = 0; p < c.size(); ++p) {
        const std::vector<std::size_t> prefix(c.begin(), c.begin() + p + 1);
        const T mu_cur = it->second.svc_rate_fun(prefix);
        const T ratep = T(mu_cur - mu_prev);
        mu_prev = mu_cur;
        if (num_traits<T>::to_double(ratep) <= 0) continue;  // position unserved
        const std::pair<std::vector<std::size_t>, std::size_t> ps =
            pass_and_swap<T>(c, p, it->second.swap_graph);
        // The completing token need not eject its own class, so only the
        // positions whose chain ends in THIS class contribute to its departure.
        if (ps.second != cls) continue;
        std::vector<T> row;
        for (std::size_t i = 0; i < ps.first.size(); ++i)
            row.push_back(num_traits<T>::from_int(static_cast<long>(ps.first[i])));
        row.resize(W, zero);
        row.insert(row.end(), var.begin(), var.end());
        out.space.push_back(row);
        out.rate.push_back(ratep);
        out.prob.push_back(one);
        pas_tag_started<T>(out, it->second.svc_rate_fun, c, ps.first);
    }
    pad_tags(out);
    return out;
}

/**
 * Port of `State.replyBlockInfo`: the layout of the reply block.
 *
 * The block trails the modulation, routing and shared-node columns of nvars,
 * occupying columns 2R+1+r. Appending is deliberate -- every existing nvars
 * reader keeps its indices, and the columns stay zero-width for models without
 * reply signals, so no other model changes state width.
 */
template <class T>
ReplyBlockInfo reply_block_info(const NetworkStruct<T>& sn, std::size_t ind) {
    const std::size_t R = sn.nclasses;
    ReplyBlockInfo ri;
    ri.slot.assign(R, static_cast<std::size_t>(-1));
    if (sn.nvars.size() < ind || sn.nvars[ind - 1].size() < 3 * R + 1) return ri;
    std::size_t pos = 0;
    for (std::size_t j = 0; j < 2 * R + 1; ++j) pos += sn.nvars[ind - 1][j];
    for (std::size_t r = 1; r <= R; ++r)
        if (sn.nvars[ind - 1][2 * R + r] > 0) {
            ri.slot[r - 1] = pos;
            ri.classes.push_back(r);
            ++pos;
            ++ri.width;
        }
    return ri;
}

/** How many servers node `ind` is holding for pending replies, from its vars. */
template <class T>
double reply_blocked(const NetworkStruct<T>& sn, std::size_t ind, const std::vector<T>& var) {
    const ReplyBlockInfo ri = reply_block_info(sn, ind);
    double nb = 0;
    for (std::size_t i = 0; i < ri.classes.size(); ++i) {
        const std::size_t s = ri.slot[ri.classes[i] - 1];
        // The slot is an index into the local-variable block, which is the
        // TAIL of the row, so it is offset from the start of `var`.
        if (s != static_cast<std::size_t>(-1) && s < var.size())
            nb += num_traits<T>::to_double(var[s]);
    }
    return nb;
}

/**
 * Port of `State.afterEventStationReply`: a REPLY signal completes a
 * synchronous call at the station holding the server for it.
 *
 * A REPLY is not a negative customer. It releases one held server and then
 * JOINS as an ordinary job carrying the call result onward, so unlike
 * NEGATIVE or CATASTROPHE it is not annihilated.
 *
 * PASS-THROUGH is the subtle part: the released server is taken by the reply
 * ITSELF, never by a waiting job. The reply is work this station already paid
 * for, so queueing it behind the residents both misreports its residence and
 * steals capacity. Its service is typically Immediate, so the server is handed
 * straight back and the ordinary departure path then promotes the head of
 * line -- which also keeps the occupancy within the server count, unlike
 * admitting the reply on top of a promoted job.
 */
template <class T>
EventOutcome<T> after_event_station_reply(const NetworkStruct<T>& sn, std::size_t ind,
                                          const std::vector<T>& inspace, std::size_t cls) {
    const std::size_t R = sn.nclasses;
    const std::size_t ist = sn.nodes[ind - 1].station;
    const T one = num_traits<T>::from_int(1), zero = num_traits<T>::from_int(0);
    EventOutcome<T> out;
    if (ist == 0) throw InputError("after_event_station_reply: node is not a station");
    const RowLayout<T> L = row_layout(sn, ind, inspace.size());
    const ReplyBlockInfo ri = reply_block_info(sn, ind);
    const double S = sn.stations[ist - 1].nservers;

    // The calling class this reply releases: the one whose expected reply IS
    // this class and which holds a block here.
    std::size_t callclass = 0;
    for (std::size_t i = 0; i < ri.classes.size(); ++i) {
        const std::size_t r = ri.classes[i];
        if (sn.syncreply.size() >= r && sn.syncreply[r - 1] == cls) {
            callclass = r;
            break;
        }
    }

    std::vector<T> buf(inspace.begin(), inspace.begin() + L.bufw);
    std::vector<T> srv(inspace.begin() + L.bufw, inspace.begin() + L.bufw + L.srvw);
    std::vector<T> var(inspace.begin() + L.bufw + L.srvw, inspace.end());

    if (callclass > 0) {
        const std::size_t s = ri.slot[callclass - 1];
        if (s != static_cast<std::size_t>(-1) && s < var.size() &&
            num_traits<T>::to_double(var[s]) > 0)
            var[s] -= one;
    }

    // The reply joins: into a free server, enumerating its entry phase, or at
    // the tail of the buffer. Servers still held for OTHER pending replies are
    // not available, which is what `Seff` subtracts.
    const double seff = S - reply_blocked(sn, ind, var);
    double occ = 0;
    for (std::size_t j = 0; j < srv.size(); ++j) occ += num_traits<T>::to_double(srv[j]);
    if (occ < seff) {
        const std::vector<T> pentry = entry_phase_dist(sn, ist, cls);
        for (std::size_t ke = 0; ke < L.K[cls - 1]; ++ke) {
            if (num_traits<T>::to_double(pentry[ke]) <= 0) continue;
            std::vector<T> s2 = srv;
            s2[L.Ks[cls - 1] + ke] += one;
            std::vector<T> row = buf;
            row.insert(row.end(), s2.begin(), s2.end());
            row.insert(row.end(), var.begin(), var.end());
            out.space.push_back(row);
            out.rate.push_back(num_traits<T>::from_int(-1));
            out.prob.push_back(pentry[ke]);
        }
        return out;
    }
    // Every available server is busy: queue at the tail, which for a
    // right-aligned buffer is the LAST empty slot.
    std::vector<T> b2 = buf;
    std::size_t slot = b2.size();
    for (std::size_t b = 0; b < b2.size(); ++b)
        if (num_traits<T>::to_double(b2[b]) == 0) slot = b;
    if (slot == b2.size()) {
        b2.insert(b2.begin(), zero);
        slot = 0;
    }
    b2[slot] = num_traits<T>::from_int(static_cast<long>(cls));
    std::vector<T> row = b2;
    row.insert(row.end(), srv.begin(), srv.end());
    row.insert(row.end(), var.begin(), var.end());
    out.space.push_back(row);
    out.rate.push_back(num_traits<T>::from_int(-1));
    out.prob.push_back(one);
    return out;
}


/** Read (pos, swk, ctr) out of the local-variable block. */
template <class T>
void polling_get(const PollingInfo<T>& pi, const std::vector<T>& var, std::size_t srvclass,
                 std::size_t& pos, std::size_t& swk, long& ctr) {
    pos = pi.ipos != static_cast<std::size_t>(-1)
              ? static_cast<std::size_t>(num_traits<T>::to_double(var[pi.off + pi.ipos]))
              : (srvclass > 0 ? srvclass : 1);
    swk = pi.iswk != static_cast<std::size_t>(-1)
              ? static_cast<std::size_t>(num_traits<T>::to_double(var[pi.off + pi.iswk]))
              : 0;
    ctr = pi.ictr != static_cast<std::size_t>(-1)
              ? static_cast<long>(num_traits<T>::to_double(var[pi.off + pi.ictr]))
              : 0;
}

/** Write (pos, swk, ctr) back into the local-variable block. */
template <class T>
std::vector<T> polling_set(const PollingInfo<T>& pi, std::vector<T> var, std::size_t pos,
                           std::size_t swk, long ctr) {
    if (pi.ipos != static_cast<std::size_t>(-1))
        var[pi.off + pi.ipos] = num_traits<T>::from_int(static_cast<long>(pos));
    if (pi.iswk != static_cast<std::size_t>(-1))
        var[pi.off + pi.iswk] = num_traits<T>::from_int(static_cast<long>(swk));
    if (pi.ictr != static_cast<std::size_t>(-1))
        var[pi.off + pi.ictr] = num_traits<T>::from_int(ctr);
    return var;
}

/** Port of `State.pollingBudget`: how many services this visit may perform. */
template <class T>
long polling_budget(const PollingInfo<T>& pi, long nbufq) {
    switch (pi.ptype) {
        case lang::PollingType::EXHAUSTIVE: return 0;  // unused: drains instead
        case lang::PollingType::GATED: return nbufq;   // exactly those found
        case lang::PollingType::KLIMITED: return static_cast<long>(pi.pk);
        case lang::PollingType::DECREMENTING: return nbufq - 1;
        default: throw InputError("polling_budget: unsupported polling type");
    }
}

/**
 * Port of `State.pollingNext`: where the server goes from buffer `pos`.
 *
 * Returns mode 1 to open a visit at q, 2 to enter the switchover into q, and 0
 * to PARK. Parking is reachable only with every switchover immediate, where a
 * server completing a full lap without finding work would otherwise cycle in
 * zero time forever.
 */
template <class T>
void polling_next(const PollingInfo<T>& pi, std::size_t pos, const std::vector<long>& nbuf,
                  std::size_t R, bool arrived, std::size_t& q, int& mode, long& budget) {
    if (arrived && pi.polled[pos - 1] && nbuf[pos - 1] > 0) {
        // The switchover into pos is already paid for, so a visit starts here.
        q = pos;
        mode = 1;
        budget = polling_budget(pi, nbuf[pos - 1]);
        return;
    }
    std::size_t p = pos;
    for (std::size_t step = 0; step < R; ++step) {  // a full lap, ending at pos
        p = p % R + 1;
        if (!pi.polled[p - 1]) continue;
        if (pi.has_sw[p - 1]) {
            q = p;
            mode = 2;
            budget = 0;
            return;
        }
        if (nbuf[p - 1] > 0) {
            q = p;
            mode = 1;
            budget = polling_budget(pi, nbuf[p - 1]);
            return;
        }
    }
    q = pos;
    mode = 0;
    budget = 0;
}

/** Port of `State.pollingLand`: the states the walk lands in, with weights. */
template <class T>
void polling_land(const NetworkStruct<T>& sn, std::size_t ist, const PollingInfo<T>& pi,
                  std::size_t q, int mode, long budget, const std::vector<T>& buf,
                  const std::vector<T>& srv, const std::vector<T>& var, const RowLayout<T>& L,
                  std::vector<std::vector<T>>& rows, std::vector<T>& probs) {
    const T one = num_traits<T>::from_int(1);
    if (mode == 1) {
        // Open or continue a visit at q: pull a waiting class-q job in.
        std::vector<T> b2 = buf;
        b2[q - 1] -= one;
        const std::vector<T> pentry = entry_phase_dist(sn, ist, q);
        for (std::size_t ke = 0; ke < L.K[q - 1]; ++ke) {
            if (num_traits<T>::to_double(pentry[ke]) <= 0) continue;
            std::vector<T> s2 = srv;
            s2[L.Ks[q - 1] + ke] += one;
            std::vector<T> row = b2;
            row.insert(row.end(), s2.begin(), s2.end());
            const std::vector<T> v2 = polling_set(pi, var, q, 0, budget);
            row.insert(row.end(), v2.begin(), v2.end());
            rows.push_back(row);
            probs.push_back(pentry[ke]);
        }
    } else if (mode == 2) {
        // Enter the switchover into q: the facility stays EMPTY while walking.
        for (std::size_t ke = 0; ke < pi.ksw[q - 1]; ++ke) {
            if (num_traits<T>::to_double(pi.sw_pie[q - 1][ke]) <= 0) continue;
            std::vector<T> row = buf;
            row.insert(row.end(), srv.begin(), srv.end());
            const std::vector<T> v2 = polling_set(pi, var, q, ke + 1, 0);
            row.insert(row.end(), v2.begin(), v2.end());
            rows.push_back(row);
            probs.push_back(pi.sw_pie[q - 1][ke]);
        }
    } else {
        // Park: held until the next arrival wakes the server.
        std::vector<T> row = buf;
        row.insert(row.end(), srv.begin(), srv.end());
        const std::vector<T> v2 = polling_set(pi, var, q, 0, 0);
        row.insert(row.end(), v2.begin(), v2.end());
        rows.push_back(row);
        probs.push_back(one);
    }
}


/**
 * Port of the SWITCH branch: a polling server advances its switchover timer.
 *
 * Unlike PHASE, which carries only the internal transitions of a phase-type and
 * leaves the absorption to DEP, this event carries BOTH -- a completed
 * switchover moves no job, so there is no departure to attach the absorption
 * to. It is therefore emitted even for a single-phase switchover, where it
 * consists of the absorption alone.
 */
template <class T>
EventOutcome<T> after_event_station_switch(const NetworkStruct<T>& sn, std::size_t ind,
                                           const std::vector<T>& inspace, std::size_t cls) {
    const std::size_t R = sn.nclasses;
    const std::size_t ist = sn.nodes[ind - 1].station;
    const T one = num_traits<T>::from_int(1);
    EventOutcome<T> out;
    if (ist == 0) throw InputError("after_event_station_switch: node is not a station");
    const PollingInfo<T> pinfo = polling_info(sn, ind);
    if (!pinfo.valid || !pinfo.has_sw[cls - 1]) return out;
    const RowLayout<T> L = row_layout(sn, ind, inspace.size());

    const std::vector<T> buf(inspace.begin(), inspace.begin() + L.bufw);
    const std::vector<T> srv(inspace.begin() + L.bufw, inspace.begin() + L.bufw + L.srvw);
    const std::vector<T> var(inspace.begin() + L.bufw + L.srvw, inspace.end());

    std::size_t pos = 0, swk = 0;
    long ctr = 0;
    polling_get(pinfo, var, 0, pos, swk, ctr);
    // The server must actually be inside the switchover into this buffer.
    if (pos != cls || swk == 0) return out;

    // Internal transitions of the switchover phase-type.
    for (std::size_t kd = 0; kd < pinfo.ksw[cls - 1]; ++kd) {
        if (kd + 1 == swk) continue;
        const T r0 = pinfo.sw_d0[cls - 1](swk - 1, kd);
        if (num_traits<T>::to_double(r0) <= 0) continue;
        std::vector<T> row = buf;
        row.insert(row.end(), srv.begin(), srv.end());
        const std::vector<T> v2 = polling_set(pinfo, var, cls, kd + 1, 0);
        row.insert(row.end(), v2.begin(), v2.end());
        out.space.push_back(row);
        out.rate.push_back(r0);
        out.prob.push_back(one);
    }
    // Absorption: the server arrives and either opens a visit or walks on.
    T rate = num_traits<T>::from_int(0);
    for (std::size_t j = 0; j < pinfo.ksw[cls - 1]; ++j) rate += pinfo.sw_d1[cls - 1](swk - 1, j);
    if (num_traits<T>::to_double(rate) <= 0) return out;
    std::vector<long> nbuf(R, 0);
    for (std::size_t r = 0; r < R && r < L.bufw; ++r)
        nbuf[r] = static_cast<long>(num_traits<T>::to_double(buf[r]));
    std::size_t q = 0;
    int mode = 0;
    long budget = 0;
    polling_next(pinfo, cls, nbuf, R, true, q, mode, budget);
    std::vector<std::vector<T>> rows;
    std::vector<T> probs;
    polling_land(sn, ist, pinfo, q, mode, budget, buf, srv, var, L, rows, probs);
    for (std::size_t j = 0; j < rows.size(); ++j) {
        // A switchover completing over an empty buffer starts the next leg at
        // once, and when that leg re-enters the SAME phase of the same
        // switchover the landing state IS the departure state. Such a self-loop
        // is not a transition: emitting it would inflate the row's exit rate.
        if (rows[j] == inspace) continue;
        out.space.push_back(rows[j]);
        out.rate.push_back(T(rate * probs[j]));
        out.prob.push_back(one);
        // A completed switchover that opens a visit pulls a waiting class-q job
        // into the server, so it starts service just as an ARV or a DEP
        // promotion does. This is the one service start a polling station
        // reaches through neither, and leaving it untagged would break
        // startRate == TN + preemptRate there for no reason other than the name
        // of the carrier event.
        tag_last(out, mode == 1 ? q : 0, 0);
    }
    pad_tags(out);
    return out;
}


/**
 * Port of `State.afterEventCache`: events at a Cache node.
 *
 * A Cache is stateful but is NOT a station, so its row is [per-class counts |
 * cache contents | retrieval bitmap] with no buffer or server block. The
 * contents region holds one column per cached slot, laid out list by list;
 * `cpos(i,j)` is position j of list i.
 *
 * A READ is INSTANTANEOUS: every branch fires at `GlobalConstants::Immediate`,
 * because the read is a routing decision rather than a service. The job enters
 * in its read class and leaves in the hit or miss class, so the transition
 * both moves the job between classes and rewrites the cache contents.
 */
template <class T>
EventOutcome<T> after_event_cache(const NetworkStruct<T>& sn, std::size_t ind,
                                  const std::vector<T>& inspace, EventType event,
                                  std::size_t cls) {
    const std::size_t R = sn.nclasses;
    const T one = num_traits<T>::from_int(1), zero = num_traits<T>::from_int(0);
    const T imm = num_traits<T>::from_double(GlobalConstants::Immediate);
    EventOutcome<T> out;
    const typename std::map<std::size_t, CacheParam<T>>::const_iterator ci = sn.nodeparam.find(ind);
    if (ci == sn.nodeparam.end()) throw InputError("after_event_cache: node has no CacheParam");
    const CacheParam<T>& cp = ci->second;
    const std::size_t h = cp.itemcap.size();
    const std::size_t n = cp.nitems;
    std::size_t tcc = 0;  // total cache capacity, the width of the contents region
    for (std::size_t i = 0; i < h; ++i)
        if (cp.itemcap[i] > 0) tcc += static_cast<std::size_t>(cp.itemcap[i]);
    // cpos(i,j): position j (1-based) of list i (1-based) in the contents region.
    const std::vector<int>& m = cp.itemcap;
    struct Cpos {
        const std::vector<int>& m;
        std::size_t operator()(std::size_t i, std::size_t j) const {
            std::size_t base = 0;
            for (std::size_t t = 0; t + 1 < i; ++t) base += static_cast<std::size_t>(m[t]);
            return base + j - 1;
        }
    } cpos{m};

    std::vector<T> srv(inspace.begin(), inspace.begin() + R);
    std::vector<T> var(inspace.begin() + R, inspace.end());

    if (event == EventType::ARV) {
        srv[cls - 1] += one;
        std::vector<T> row = srv;
        row.insert(row.end(), var.begin(), var.end());
        out.space.push_back(row);
        out.rate.push_back(num_traits<T>::from_int(-1));  // passive
        out.prob.push_back(one);
        return out;
    }

    if (event == EventType::DEP) {
        if (num_traits<T>::to_double(srv[cls - 1]) <= 0) return out;
        // A departure only moves the job: the occupancy bit of a fetch was
        // already set by the READ that began it, exactly as the reference's DEP
        // branch does nothing beyond decrementing the class.
        srv[cls - 1] -= one;
        std::vector<T> row = srv;
        row.insert(row.end(), var.begin(), var.end());
        out.space.push_back(row);
        out.rate.push_back(imm);  // the departure is instantaneous
        out.prob.push_back(one);
        return out;
    }

    if (event != EventType::READ) return out;

    // A READ needs exactly one job present, and it must be of the reading class.
    double tot = 0;
    for (std::size_t r = 0; r < R; ++r) tot += num_traits<T>::to_double(srv[r]);
    if (num_traits<T>::to_double(srv[cls - 1]) <= 0 || tot != 1) return out;
    if (cls - 1 >= cp.pread.size() || cp.pread[cls - 1].empty()) return out;
    const std::vector<T>& p = cp.pread[cls - 1];

    // The delayed-hit retrieval system. Block A (n columns) marks the items in
    // flight; block B (one column per retrieval class) counts the secondary
    // requests merged onto those fetches. Its width is read off the row rather
    // than off the parameters, because a struct whose space predates the block
    // still has to be walkable.
    const bool retr = cp.retrieval_capacity > 0 && !cp.retrieval_classes.empty();
    std::vector<std::size_t> rc_list, rc_items, rc_orig;
    if (retr) cache_retrieval_class_map(cp, rc_list, rc_items, rc_orig);
    const std::size_t block_b = tcc + n;
    std::size_t width_b = var.size() > block_b ? var.size() - block_b : 0;
    if (width_b != rc_list.size()) width_b = 0;
    // -1 is the sample path's unbounded merge; an exact solver enumerates block
    // B only up to the level it declared, and a merge past that level is refused
    // rather than folded into a state the space does not hold.
    const double maxpend = cp.max_pending_retrieval < 0
                               ? std::numeric_limits<double>::infinity()
                               : static_cast<double>(cp.max_pending_retrieval);
    bool from_retrieval = false;
    for (std::size_t j = 0; j < rc_list.size(); ++j)
        if (rc_list[j] == cls) from_retrieval = true;

    for (std::size_t k = 1; k <= n; ++k) {
        if (k - 1 >= p.size() || num_traits<T>::to_double(p[k - 1]) <= 0) continue;
        std::vector<T> srv_e = srv;
        srv_e[cls - 1] -= one;
        // The item is searched ONLY in the contents region; the trailing
        // retrieval slots are a different namespace.
        std::size_t posk = 0;
        for (std::size_t c = 0; c < tcc && c < var.size(); ++c)
            if (static_cast<std::size_t>(num_traits<T>::to_double(var[c])) == k) {
                posk = c + 1;
                break;
            }
        // A RETURNING RETRIEVAL always completes the miss that started it, so it
        // takes the miss branch even where the enumeration produced a (then
        // unreachable) state holding item k. A retrieval class has no hit class
        // to switch into either.
        if (from_retrieval) posk = 0;
        const Matrix<T>& ac = cp.accost.empty() || cls - 1 >= cp.accost.size() ||
                                      k - 1 >= cp.accost[cls - 1].size()
                                  ? Matrix<T>()
                                  : cp.accost[cls - 1][k - 1];
        const bool have_ac = ac.rows() >= h + 1 && ac.cols() >= h + 1;

        if (posk == 0) {
            // CACHE MISS, or one leg of a retrieval. A fetch is in flight iff
            // block A's bit for item k is set.
            const bool in_flight =
                retr && tcc + k <= var.size() && num_traits<T>::to_double(var[tcc + k - 1]) != 0;
            std::size_t r_class = 0;
            if (retr && k - 1 < cp.retrieval_classes.size() &&
                cls - 1 < cp.retrieval_classes[k - 1].size())
                r_class = cp.retrieval_classes[k - 1][cls - 1];
            // A returning retrieval not recorded in the bitmap is an unreachable
            // artifact of the enumeration; it is not continued.
            if (from_retrieval && !in_flight) continue;

            if (!from_retrieval && r_class != 0) {
                // BEGIN a fetch, or MERGE onto one already running. The merge is
                // the delayed hit: it adds nothing to the cache's server block,
                // and is held in block B until the fetch completes and releases
                // it in the hit class of the class that issued it.
                std::vector<T> srv_b = srv_e;
                std::vector<T> var_b = var;
                if (!in_flight) {
                    srv_b[r_class - 1] += one;
                    if (tcc + k <= var_b.size()) var_b[tcc + k - 1] = one;
                } else {
                    if (width_b == 0) continue;
                    std::size_t bslot = width_b;
                    for (std::size_t j = 0; j < rc_list.size(); ++j)
                        if (rc_list[j] == r_class) { bslot = j; break; }
                    if (bslot >= width_b) continue;
                    double pend = 0;
                    for (std::size_t j = 0; j < width_b; ++j)
                        pend += num_traits<T>::to_double(var_b[block_b + j]);
                    if (pend >= maxpend) continue;  // beyond the truncation level
                    var_b[block_b + bslot] += one;
                }
                std::vector<T> row = srv_b;
                row.insert(row.end(), var_b.begin(), var_b.end());
                out.space.push_back(row);
                out.rate.push_back(T(p[k - 1] * imm));
                out.prob.push_back(one);
                continue;
            }

            // The fetch is complete (or there is no retrieval system): the job
            // leaves in the miss class, the item may be admitted into one of the
            // lists, and every request merged onto this fetch is released in the
            // SAME transition as a delayed hit.
            if (cls - 1 >= cp.missclass.size() || cp.missclass[cls - 1] == 0) continue;
            std::vector<T> srv_m = srv_e;
            srv_m[cp.missclass[cls - 1] - 1] += one;
            std::vector<T> var_m = var;
            if (tcc + k <= var_m.size()) var_m[tcc + k - 1] = zero;  // retrieval done
            for (std::size_t j = 0; j < width_b; ++j) {
                if (rc_items[j] != k) continue;
                const double held = num_traits<T>::to_double(var_m[block_b + j]);
                if (held <= 0) continue;
                const std::size_t oc = rc_orig[j];
                if (oc - 1 >= cp.hitclass.size() || cp.hitclass[oc - 1] == 0) continue;
                srv_m[cp.hitclass[oc - 1] - 1] += var_m[block_b + j];
                var_m[block_b + j] = zero;
            }

            // Column 1 of the access cost is the REJECT branch: the item passes
            // through without being cached at all.
            const T rej = have_ac ? ac(0, 0) : zero;
            if (num_traits<T>::to_double(rej) > 0) {
                std::vector<T> row = srv_m;
                row.insert(row.end(), var_m.begin(), var_m.end());
                out.space.push_back(row);
                out.rate.push_back(T(rej * p[k - 1] * imm));
                out.prob.push_back(one);
            }
            for (std::size_t l = 1; l <= h; ++l) {
                const T w = have_ac ? ac(0, l) : (l == 1 ? one : zero);
                if (num_traits<T>::to_double(w) <= 0) continue;
                if (m[l - 1] <= 0) continue;
                const std::size_t ml = static_cast<std::size_t>(m[l - 1]);
                if (cp.replacestrat == lang::ReplacementStrategy::RR) {
                    // Random replacement: the item lands uniformly in any slot.
                    for (std::size_t rr = 1; rr <= ml; ++rr) {
                        std::vector<T> vp = var_m;
                        vp[cpos(l, rr)] = num_traits<T>::from_int(static_cast<long>(k));
                        std::vector<T> row = srv_m;
                        row.insert(row.end(), vp.begin(), vp.end());
                        out.space.push_back(row);
                        out.rate.push_back(T(w * p[k - 1] /
                                             num_traits<T>::from_int(static_cast<long>(ml)) * imm));
                        out.prob.push_back(one);
                    }
                } else {
                    // The ordered families insert at the HEAD, shifting the list
                    // down by one and evicting its tail.
                    std::vector<T> vp = var_m;
                    for (std::size_t j = ml; j >= 2; --j) vp[cpos(l, j)] = var_m[cpos(l, j - 1)];
                    vp[cpos(l, 1)] = num_traits<T>::from_int(static_cast<long>(k));
                    T rate = T(w * p[k - 1] * imm);
                    // q-LRU admits a miss only with probability q; the rest
                    // passes through uncached.
                    if (cp.replacestrat == lang::ReplacementStrategy::QLRU) {
                        const T q = cp.qlru;
                        if (num_traits<T>::to_double(q) < 1) {
                            std::vector<T> row0 = srv_m;
                            row0.insert(row0.end(), var_m.begin(), var_m.end());
                            out.space.push_back(row0);
                            out.rate.push_back(T(rate * T(one - q)));
                            out.prob.push_back(one);
                        }
                        rate = T(rate * q);
                    }
                    if (num_traits<T>::to_double(rate) <= 0) continue;
                    std::vector<T> row = srv_m;
                    row.insert(row.end(), vp.begin(), vp.end());
                    out.space.push_back(row);
                    out.rate.push_back(rate);
                    out.prob.push_back(one);
                }
            }
        } else {
            // CACHE HIT: the job leaves in the hit class. Which list it was
            // found in decides how the contents move.
            if (cls - 1 >= cp.hitclass.size() || cp.hitclass[cls - 1] == 0) continue;
            std::vector<T> srv_h = srv_e;
            srv_h[cp.hitclass[cls - 1] - 1] += one;
            std::size_t li = 1, acc = 0;
            for (std::size_t t = 0; t < h; ++t) {
                const std::size_t mt = m[t] > 0 ? static_cast<std::size_t>(m[t]) : 0;
                if (posk <= acc + mt) { li = t + 1; break; }
                acc += mt;
            }
            const std::size_t j = posk - acc;
            if (li < h) {
                // A HIT BELOW THE TERMINAL LIST PROMOTES, and that is what makes
                // an h-list cache more than h caches side by side. Row `li` of
                // the access cost routes the item from list li into list
                // inew >= li, exactly as `afterEventCache.m` does over `inew =
                // i:h`; an empty accost is the reference default
                // `diag(ones(1,h),1)` with a 1 in the bottom-right, the linear
                // cache that moves the item one list up. Omitting this branch
                // froze every list above the first at its initial contents: on
                // m=[2,1] only 12 of the 60 configurations stayed reachable, and
                // the CTMC hit ratio of cache_compare_replc came out 4.9% low.
                for (std::size_t inew = li; inew <= h; ++inew) {
                    if (m[inew - 1] <= 0) continue;
                    const std::size_t mn = static_cast<std::size_t>(m[inew - 1]);
                    const T w = have_ac ? ac(li, inew) : (inew == li + 1 ? one : zero);
                    if (num_traits<T>::to_double(w) <= 0) continue;
                    if (cp.replacestrat == lang::ReplacementStrategy::RR) {
                        // Random replacement swaps with a uniformly drawn slot.
                        for (std::size_t r = 1; r <= mn; ++r) {
                            std::vector<T> vp = var;
                            vp[cpos(li, j)] = var[cpos(inew, r)];
                            vp[cpos(inew, r)] = num_traits<T>::from_int(static_cast<long>(k));
                            std::vector<T> row = srv_h;
                            row.insert(row.end(), vp.begin(), vp.end());
                            out.space.push_back(row);
                            out.rate.push_back(
                                T(w * p[k - 1] /
                                  num_traits<T>::from_int(static_cast<long>(mn)) * imm));
                            out.prob.push_back(one);
                        }
                        continue;
                    }
                    // Every read below is of the UNMODIFIED row, so the three
                    // moves compose in the reference's order even where the
                    // source list and the target list are the same one.
                    std::vector<T> vp = var;
                    // The LRU family closes the gap in list li; FIFO orders by
                    // insertion and so leaves list li otherwise untouched.
                    const bool ordered = cp.replacestrat != lang::ReplacementStrategy::FIFO;
                    if (ordered)
                        for (std::size_t t = j; t >= 2; --t) vp[cpos(li, t)] = var[cpos(li, t - 1)];
                    // The tail evicted from the target list takes the slot the
                    // promoted item vacated.
                    vp[cpos(li, ordered ? 1 : j)] = var[cpos(inew, mn)];
                    for (std::size_t t = mn; t >= 2; --t) vp[cpos(inew, t)] = var[cpos(inew, t - 1)];
                    vp[cpos(inew, 1)] = num_traits<T>::from_int(static_cast<long>(k));
                    std::vector<T> row = srv_h;
                    row.insert(row.end(), vp.begin(), vp.end());
                    out.space.push_back(row);
                    out.rate.push_back(T(w * p[k - 1] * imm));
                    out.prob.push_back(one);
                }
            } else if (cp.replacestrat == lang::ReplacementStrategy::RR ||
                       cp.replacestrat == lang::ReplacementStrategy::FIFO ||
                       cp.replacestrat == lang::ReplacementStrategy::SFIFO) {
                // A hit in the terminal list does not reorder these: FIFO orders
                // by INSERTION, and random replacement has no order to disturb.
                std::vector<T> row = srv_h;
                row.insert(row.end(), var.begin(), var.end());
                out.space.push_back(row);
                out.rate.push_back(T(p[k - 1] * imm));
                out.prob.push_back(one);
            } else {
                // LRU and its relatives promote the hit item to the head of its
                // list, which is the whole content of "recently used".
                std::vector<T> vp = var;
                for (std::size_t t = j; t >= 2; --t) vp[cpos(li, t)] = var[cpos(li, t - 1)];
                vp[cpos(li, 1)] = var[cpos(li, j)];
                std::vector<T> row = srv_h;
                row.insert(row.end(), vp.begin(), vp.end());
                out.space.push_back(row);
                out.rate.push_back(T(p[k - 1] * imm));
                out.prob.push_back(one);
            }
        }
    }
    return out;
}


/** One half of a GLOBAL synchronization: a mode event at a node. */
template <class T>
struct ModeEvent {
    EventType event = EventType::LOCAL;
    std::size_t node = 0;   ///< 1-based node index (a Transition, or a place)
    std::size_t mode = 0;   ///< 1-based mode index
    std::size_t cls = 1;    ///< 1-based class the arc moves
    T weight = num_traits<T>::from_int(1);  ///< arc multiplicity
};

/**
 * A GLOBAL synchronization: an SPN mode event and the place arcs it drives.
 *
 * Unlike an ordinary Sync, which pairs ONE active with ONE passive, a firing
 * touches every input and output place at once -- that atomicity is what makes
 * a Petri net transition a transition. PRE passives consume, POST produce, and
 * LOCAL passives are read-only (an inhibitor place, whose marking is tested but
 * never moved).
 */
template <class T>
struct GlobalSync {
    ModeEvent<T> active;
    std::vector<ModeEvent<T>> passive;
};

/**
 * Port of `MNetwork.refreshGlobalSync`: the ENABLE and FIRE synchronizations.
 *
 * An inhibiting place enters as a LOCAL passive rather than a PRE, and only
 * when it is not already an enabling or firing place: its marking is read for
 * the inhibition test but no token crosses the arc.
 */
template <class T>
std::vector<GlobalSync<T>> refresh_global_sync(const NetworkStruct<T>& sn) {
    std::vector<GlobalSync<T>> gsync;
    const T one = num_traits<T>::from_int(1);
    for (std::size_t ind = 1; ind <= sn.nodes.size(); ++ind) {
        if (sn.nodes[ind - 1].nodetype != NodeType::Transition) continue;
        const typename std::map<std::size_t, TransitionParam<T>>::const_iterator it =
            sn.transparam.find(ind);
        if (it == sn.transparam.end()) continue;
        const TransitionParam<T>& tp = it->second;
        for (int pass = 0; pass < 2; ++pass) {
            for (std::size_t m = 1; m <= tp.nmodes; ++m) {
                // ONE ENTRY PER (place, class) ARC, which is what makes a
                // multiclass net a different net from the class-summed one: a
                // Class2 token at a place must not satisfy a Class1 pre-arc, so
                // the pair travels into the passive rather than the place alone.
                std::vector<std::pair<std::size_t, std::size_t>> enab, fire, inhib;
                const Matrix<T>& en = tp.enabling[m - 1];
                const Matrix<T>& fi = tp.firing[m - 1];
                const Matrix<T>& ih = tp.inhibiting[m - 1];
                for (std::size_t q = 0; q < en.rows(); ++q)
                    for (std::size_t r = 0; r < en.cols(); ++r)
                        if (num_traits<T>::to_double(en(q, r)) > 0)
                            enab.push_back(std::make_pair(q + 1, r + 1));
                for (std::size_t q = 0; q < fi.rows(); ++q)
                    for (std::size_t r = 0; r < fi.cols(); ++r)
                        if (num_traits<T>::to_double(fi(q, r)) > 0)
                            fire.push_back(std::make_pair(q + 1, r + 1));
                for (std::size_t q = 0; q < ih.rows(); ++q)
                    for (std::size_t r = 0; r < ih.cols(); ++r) {
                        if (std::isinf(num_traits<T>::to_double(ih(q, r)))) continue;
                        // The de-duplication is against the PLACE, not the
                        // (place, class) pair: the passive's job is to bring the
                        // place's marginal into the outcome, and one copy of a
                        // place carries every class of it.
                        bool dup = false;
                        for (std::size_t i = 0; i < enab.size(); ++i)
                            dup = dup || enab[i].first == q + 1;
                        for (std::size_t i = 0; i < fire.size(); ++i)
                            dup = dup || fire[i].first == q + 1;
                        for (std::size_t i = 0; i < inhib.size(); ++i)
                            dup = dup || inhib[i].first == q + 1;
                        if (!dup) inhib.push_back(std::make_pair(q + 1, r + 1));
                    }
                GlobalSync<T> g;
                g.active.event = pass == 0 ? EventType::ENABLE : EventType::FIRE;
                g.active.node = ind;
                g.active.mode = m;
                if (pass == 0) {
                    // An ENABLE only READS the markings, so every passive is
                    // LOCAL: enabling it is a test, not a token movement. One
                    // per PLACE here, since a read of a place reads every class.
                    std::vector<std::size_t> seen;
                    auto once = [&](std::size_t q, std::size_t r) {
                        for (std::size_t i = 0; i < seen.size(); ++i)
                            if (seen[i] == q) return;
                        seen.push_back(q);
                        g.passive.push_back(ModeEvent<T>{EventType::LOCAL, q, m, r, one});
                    };
                    for (std::size_t i = 0; i < enab.size(); ++i) once(enab[i].first, enab[i].second);
                    for (std::size_t i = 0; i < inhib.size(); ++i)
                        once(inhib[i].first, inhib[i].second);
                } else {
                    for (std::size_t i = 0; i < enab.size(); ++i)
                        g.passive.push_back(ModeEvent<T>{EventType::PRE, enab[i].first, m,
                                                         enab[i].second,
                                                         en(enab[i].first - 1, enab[i].second - 1)});
                    for (std::size_t i = 0; i < fire.size(); ++i)
                        g.passive.push_back(ModeEvent<T>{EventType::POST, fire[i].first, m,
                                                         fire[i].second,
                                                         fi(fire[i].first - 1, fire[i].second - 1)});
                    for (std::size_t i = 0; i < inhib.size(); ++i)
                        g.passive.push_back(ModeEvent<T>{EventType::LOCAL, inhib[i].first, m,
                                                         inhib[i].second, one});
                }
                gsync.push_back(g);
            }
        }
    }
    return gsync;
}

/** What one global event produces: a whole network state per outcome. */
template <class T>
struct GlobalOutcome {
    std::vector<NetState<T>> space;
    std::vector<T> rate, prob;
    /**
     * True where the outcome is a firing COMPLETION, i.e. one that applied the
     * PRE/POST updates. Callers must NOT re-derive this from the markings: a
     * transition whose firing returns exactly what its enabling consumed
     * leaves every marking invariant yet still completed.
     */
    std::vector<bool> completion;
    bool empty() const { return space.empty(); }
};

/**
 * Port of `State.afterGlobalEvent`: an SPN mode ENABLEs or FIREs.
 *
 * This is the one handler that rewrites SEVERAL nodes at once, because a
 * firing is atomic across all its arcs. The Transition's own row records how
 * many servers of each mode are idle, running (and in which firing phase), and
 * have just fired; the places are rewritten through the PRE and POST passives.
 */
template <class T>
GlobalOutcome<T> after_global_event(const NetworkStruct<T>& sn, const NetState<T>& glspace,
                                    const GlobalSync<T>& gl) {
    const std::size_t R = sn.nclasses;
    const std::size_t ind = gl.active.node;
    const std::size_t mode = gl.active.mode;
    const T one = num_traits<T>::from_int(1), zero = num_traits<T>::from_int(0);
    const T imm = num_traits<T>::from_double(GlobalConstants::Immediate);
    GlobalOutcome<T> out;
    const std::size_t isf = sn.stateful_index(ind);
    if (isf == 0) throw InputError("after_global_event: the transition is not stateful");
    const typename std::map<std::size_t, TransitionParam<T>>::const_iterator it =
        sn.transparam.find(ind);
    if (it == sn.transparam.end()) throw InputError("after_global_event: node has no TransitionParam");
    const TransitionParam<T>& tp = it->second;

    std::vector<std::size_t> fK(tp.nmodes, 1), fKs(tp.nmodes, 0);
    std::size_t tot = 0;
    for (std::size_t m = 0; m < tp.nmodes; ++m) {
        fK[m] = m < tp.firingphases.size() && tp.firingphases[m] > 0 ? tp.firingphases[m] : 1;
        fKs[m] = tot;
        tot += fK[m];
    }
    const std::vector<T>& row = glspace.local[isf - 1];
    std::vector<T> buf(row.begin(), row.begin() + tp.nmodes);
    std::vector<T> srv(row.begin() + tp.nmodes, row.begin() + tp.nmodes + tot);
    std::vector<T> fired(row.begin() + tp.nmodes + tot,
                         row.begin() + 2 * tp.nmodes + tot);
    const std::vector<T> var(row.begin() + 2 * tp.nmodes + tot, row.end());

    // The marking of every place this mode reads, node-indexed.
    std::vector<std::vector<T>> ep(sn.nodes.size() + 1, std::vector<T>(R, zero));
    for (std::size_t j = 0; j < gl.passive.size(); ++j) {
        const std::size_t pn = gl.passive[j].node;
        const std::size_t pisf = sn.stateful_index(pn);
        if (pisf == 0) continue;
        const std::pair<T, std::vector<T>> mg = to_marginal_aggr(sn, pn, glspace.local[pisf - 1]);
        ep[pn] = mg.second;
    }

    // The enabling DEGREE: how many concurrent firings the marking supports.
    // An inhibitor arc disables the mode outright once its threshold is met.
    //
    // EVERY TEST IS PER (place, class), the elementwise comparison
    // `afterGlobalEvent.m:85` makes against `enabling_m`. Summing the marking
    // over classes first, as this port did until 2026-08-12, let a Class2 token
    // satisfy a Class1 pre-arc: on a net where Mode1 needs two Class1 tokens at
    // P1 and Mode2 one Class2 token there, the summed test fires Mode1 off a
    // marking that holds no Class1 token at all.
    const Matrix<T>& en_m = tp.enabling[mode - 1];
    const Matrix<T>& ih_m = tp.inhibiting[mode - 1];
    bool inhibited = false;
    for (std::size_t q = 0; q < ih_m.rows(); ++q)
        for (std::size_t r = 0; r < ih_m.cols() && r < R; ++r) {
            const double thr = num_traits<T>::to_double(ih_m(q, r));
            if (std::isinf(thr)) continue;
            if (num_traits<T>::to_double(ep[q + 1][r]) >= thr) inhibited = true;
        }
    bool under = false;
    for (std::size_t q = 0; q < en_m.rows(); ++q)
        for (std::size_t r = 0; r < en_m.cols() && r < R; ++r) {
            const double need = num_traits<T>::to_double(en_m(q, r));
            if (need <= 0) continue;
            if (num_traits<T>::to_double(ep[q + 1][r]) < need) under = true;
        }
    long mark_degree = 0;
    if (!inhibited && !under) {
        long d = 1;
        for (;;) {
            bool ok = true;
            for (std::size_t q = 0; q < en_m.rows() && ok; ++q)
                for (std::size_t r = 0; r < en_m.cols() && r < R && ok; ++r) {
                    const double need = num_traits<T>::to_double(en_m(q, r)) * d;
                    if (need <= 0) continue;
                    if (num_traits<T>::to_double(ep[q + 1][r]) < need) ok = false;
                }
            if (!ok) break;
            ++d;
        }
        mark_degree = d - 1;
    }
    const double svm = mode - 1 < tp.nmodeservers.size() ? tp.nmodeservers[mode - 1] : 1.0;
    const long nsrv = std::isfinite(svm) ? static_cast<long>(svm)
                                         : static_cast<long>(GlobalConstants::MaxInt);

    if (gl.active.event == EventType::ENABLE) {
        long running = 0;
        for (std::size_t k = 0; k < fK[mode - 1]; ++k)
            running += static_cast<long>(num_traits<T>::to_double(srv[fKs[mode - 1] + k]));
        if (inhibited || under) {
            // Disabled: every server of this mode returns to the idle pool.
            std::vector<T> b2 = buf, s2 = srv;
            b2[mode - 1] = num_traits<T>::from_int(nsrv);
            for (std::size_t k = 0; k < fK[mode - 1]; ++k) s2[fKs[mode - 1] + k] = zero;
            std::vector<T> nr = b2;
            nr.insert(nr.end(), s2.begin(), s2.end());
            nr.insert(nr.end(), fired.begin(), fired.end());
            nr.insert(nr.end(), var.begin(), var.end());
            if (nr == row) return out;  // already disabled: not a transition
            NetState<T> ns = glspace;
            ns.local[isf - 1] = nr;
            out.space.push_back(ns);
            out.rate.push_back(imm);
            out.prob.push_back(one);
            out.completion.push_back(false);
            return out;
        }
        const long want = std::min(mark_degree, nsrv);
        if (running == want) return out;  // nothing to do
        if (running < want) {
            // Start servers, distributing them over the firing phases by the
            // entry law; the multinomial weight is the probability of that
            // split.
            const long nadd = want - running;
            std::vector<T> pe(fK[mode - 1], zero);
            if (mode - 1 < tp.firingproc.size() && tp.firingproc[mode - 1].D0.rows() ==
                                                       static_cast<std::size_t>(fK[mode - 1])) {
                mam::Map<T> mp;
                mp.D0 = tp.firingproc[mode - 1].D0;
                mp.D1 = tp.firingproc[mode - 1].D1;
                const std::vector<T> pv = mam::map_pie(mp);
                for (std::size_t k = 0; k < fK[mode - 1]; ++k) pe[k] = pv[k];
            } else {
                pe[0] = one;
            }
            // Enumerate the splits of nadd over the phases.
            std::vector<std::vector<long>> combs;
            std::vector<long> cur(fK[mode - 1], 0);
            std::function<void(std::size_t, long)> rec = [&](std::size_t k, long left) {
                if (k + 1 == fK[mode - 1]) {
                    cur[k] = left;
                    combs.push_back(cur);
                    return;
                }
                for (long v = left; v >= 0; --v) {
                    cur[k] = v;
                    rec(k + 1, left - v);
                }
            };
            rec(0, nadd);
            for (std::size_t i = 0; i < combs.size(); ++i) {
                std::vector<T> b2 = buf, s2 = srv;
                b2[mode - 1] -= num_traits<T>::from_int(nadd);
                double logp = std::lgamma(static_cast<double>(nadd) + 1.0);
                bool zeroprob = false;
                for (std::size_t k = 0; k < fK[mode - 1]; ++k) {
                    s2[fKs[mode - 1] + k] += num_traits<T>::from_int(combs[i][k]);
                    const double pk = num_traits<T>::to_double(pe[k]);
                    if (pk > 0)
                        logp += combs[i][k] * std::log(pk) -
                                std::lgamma(static_cast<double>(combs[i][k]) + 1.0);
                    else if (combs[i][k] > 0)
                        zeroprob = true;
                }
                std::vector<T> nr = b2;
                nr.insert(nr.end(), s2.begin(), s2.end());
                nr.insert(nr.end(), fired.begin(), fired.end());
                nr.insert(nr.end(), var.begin(), var.end());
                NetState<T> ns = glspace;
                ns.local[isf - 1] = nr;
                out.space.push_back(ns);
                out.rate.push_back(imm);
                out.prob.push_back(zeroprob ? zero : num_traits<T>::from_double(std::exp(logp)));
                out.completion.push_back(false);
            }
            return out;
        }
        // Stop the surplus servers, chosen uniformly across the phases: the
        // weight is the multivariate hypergeometric probability of that choice.
        const long ndiff = running - want;
        std::vector<long> sv(fK[mode - 1], 0);
        for (std::size_t k = 0; k < fK[mode - 1]; ++k)
            sv[k] = static_cast<long>(num_traits<T>::to_double(srv[fKs[mode - 1] + k]));
        std::vector<std::vector<long>> combs;
        std::vector<long> cur(fK[mode - 1], 0);
        std::function<void(std::size_t, long)> rec = [&](std::size_t k, long left) {
            if (k + 1 == fK[mode - 1]) {
                if (left > sv[k]) return;
                cur[k] = left;
                combs.push_back(cur);
                return;
            }
            for (long v = std::min(left, sv[k]); v >= 0; --v) {
                cur[k] = v;
                rec(k + 1, left - v);
            }
        };
        rec(0, ndiff);
        std::vector<double> w(combs.size(), 0.0);
        double wmax = -1e300;
        for (std::size_t i = 0; i < combs.size(); ++i) {
            double lw = 0;
            for (std::size_t k = 0; k < fK[mode - 1]; ++k)
                lw += std::lgamma(static_cast<double>(sv[k]) + 1.0) -
                      std::lgamma(static_cast<double>(combs[i][k]) + 1.0) -
                      std::lgamma(static_cast<double>(sv[k] - combs[i][k]) + 1.0);
            w[i] = lw;
            wmax = std::max(wmax, lw);
        }
        double wsum = 0;
        for (std::size_t i = 0; i < w.size(); ++i) {
            w[i] = std::exp(w[i] - wmax);
            wsum += w[i];
        }
        for (std::size_t i = 0; i < combs.size(); ++i) {
            std::vector<T> b2 = buf, s2 = srv;
            for (std::size_t k = 0; k < fK[mode - 1]; ++k)
                s2[fKs[mode - 1] + k] = num_traits<T>::from_int(sv[k] - combs[i][k]);
            b2[mode - 1] += num_traits<T>::from_int(ndiff);
            std::vector<T> nr = b2;
            nr.insert(nr.end(), s2.begin(), s2.end());
            nr.insert(nr.end(), fired.begin(), fired.end());
            nr.insert(nr.end(), var.begin(), var.end());
            NetState<T> ns = glspace;
            ns.local[isf - 1] = nr;
            out.space.push_back(ns);
            out.rate.push_back(imm);
            out.prob.push_back(num_traits<T>::from_double(wsum > 0 ? w[i] / wsum : 0.0));
            out.completion.push_back(false);
        }
        return out;
    }

    if (gl.active.event != EventType::FIRE) return out;

    const bool immediate_mode = mode - 1 < tp.timing.size() &&
                                tp.timing[mode - 1] == lang::TimingStrategy::IMMEDIATE;
    const T fw = immediate_mode && mode - 1 < tp.fireweight.size() ? tp.fireweight[mode - 1] : one;
    const long en_degree = inhibited ? 0 : std::min(mark_degree, nsrv);
    const long imm_servers = immediate_mode ? std::min(mark_degree, nsrv) : 0;

    for (std::size_t k = 0; k < fK[mode - 1]; ++k) {
        const double in_k = num_traits<T>::to_double(srv[fKs[mode - 1] + k]);
        const bool fires = immediate_mode ? (k == 0 && imm_servers >= 1)
                                          : (in_k > 0 && en_degree >= 1);
        if (!fires) continue;
        T rate = zero;
        if (immediate_mode) {
            rate = T(imm * fw * num_traits<T>::from_int(imm_servers));
        } else {
            T d1sum = zero;
            if (mode - 1 < tp.firingproc.size())
                for (std::size_t j = 0; j < tp.firingproc[mode - 1].D1.cols(); ++j)
                    d1sum += tp.firingproc[mode - 1].D1(k, j);
            rate = T(d1sum * num_traits<T>::from_double(in_k));
            // A marking-dependent firing rate g_mode(marking) is exact here,
            // because the CTMC evaluates it per enumerated state.
            if (mode - 1 < tp.firingdep.size() && tp.firingdep[mode - 1]) {
                std::vector<T> mk;
                for (std::size_t q = 1; q <= sn.nodes.size(); ++q) {
                    T s2 = zero;
                    for (std::size_t r = 0; r < R; ++r) s2 += ep[q][r];
                    mk.push_back(s2);
                }
                rate = T(rate * tp.firingdep[mode - 1](mk));
            }
        }
        if (num_traits<T>::to_double(rate) <= 0) continue;

        std::vector<T> b2 = buf, s2 = srv;
        if (in_k > 0) {
            s2[fKs[mode - 1] + k] -= one;  // the firing server leaves execution
            b2[mode - 1] += one;           // and returns to the idle pool
        }
        NetState<T> ns = glspace;
        std::vector<T> nr = b2;
        nr.insert(nr.end(), s2.begin(), s2.end());
        nr.insert(nr.end(), fired.begin(), fired.end());
        nr.insert(nr.end(), var.begin(), var.end());
        ns.local[isf - 1] = nr;

        // The arcs fire ATOMICALLY with the mode: PRE consumes from every
        // input place and POST produces into every output place, in one
        // transition. Splitting them would let the net occupy a state in which
        // the tokens have left one place and not arrived at the other.
        for (std::size_t j = 0; j < gl.passive.size(); ++j) {
            const ModeEvent<T>& pe2 = gl.passive[j];
            if (pe2.event != EventType::PRE && pe2.event != EventType::POST) continue;
            const std::size_t pisf = sn.stateful_index(pe2.node);
            if (pisf == 0) continue;
            std::vector<T>& prow = ns.local[pisf - 1];
            const double wgt = num_traits<T>::to_double(pe2.weight);
            const std::size_t pist = sn.nodes[pe2.node - 1].station;
            const SchedStrategy psched =
                pist != 0 ? sn.stations[pist - 1].sched : SchedStrategy::INF;
            const std::size_t c = pe2.cls;
            if (pe2.event == EventType::PRE) {
                if (state_detail::buffer_is_class_tag(psched)) {
                    // An ordered buffer: consume from the head end, matching
                    // the discipline rather than a count.
                    long left = static_cast<long>(wgt);
                    const T tag = num_traits<T>::from_int(static_cast<long>(c));
                    if (psched == SchedStrategy::LCFS) {
                        for (std::size_t b = 0; b < prow.size() && left > 0; ++b)
                            if (prow[b] == tag) { prow[b] = zero; --left; }
                    } else {
                        for (std::size_t b = prow.size(); b-- > 0 && left > 0;)
                            if (prow[b] == tag) { prow[b] = zero; --left; }
                    }
                } else if (prow.size() > R) {
                    // A Place with a [count | server] split: drain the server
                    // slot into the count, as the reference does.
                    const double totc = num_traits<T>::to_double(prow[c - 1]) +
                                        num_traits<T>::to_double(prow[R + c - 1]);
                    prow[c - 1] = num_traits<T>::from_double(totc - wgt);
                    prow[R + c - 1] = zero;
                } else if (c - 1 < prow.size()) {
                    prow[c - 1] -= num_traits<T>::from_double(wgt);
                }
            } else {
                if (state_detail::buffer_is_class_tag(psched)) {
                    for (long q = 0; q < static_cast<long>(wgt); ++q)
                        prow.insert(prow.begin(), num_traits<T>::from_int(static_cast<long>(c)));
                } else if (c - 1 < prow.size()) {
                    prow[c - 1] += num_traits<T>::from_double(wgt);
                }
            }
        }
        out.space.push_back(ns);
        out.rate.push_back(rate);
        out.prob.push_back(one);
        out.completion.push_back(true);
    }
    return out;
}

/** One half of a synchronization: an event at a node, in a class. */
template <class T>
struct SyncEvent {
    EventType event = EventType::LOCAL;
    std::size_t node = 0;   ///< 1-based node index, or `local` for the dummy
    std::size_t cls = 0;    ///< 1-based class index
    T prob = num_traits<T>::from_int(1);  ///< routing probability, passive half
    /**
     * The routing probability is a FUNCTION of the network state, so `prob`
     * holds only the state-independent placeholder and the generator must read
     * `rt_state(sn, state)(rt_row, rt_col)` instead.
     *
     * The reference decides this per ACTIVE NODE (`sn.isstatedep(node_a,3)`) and
     * then calls whatever `prob` holds, which throws where a node routes one
     * class state-dependently and another by probability. Deciding it per
     * SYNCHRONIZATION agrees wherever the reference runs at all, and does not
     * throw where it would.
     */
    bool statedep = false;
    std::size_t rt_row = 0;  ///< (isf-1)*nclasses + (r-1) of the active half
    std::size_t rt_col = 0;  ///< (jsf-1)*nclasses + (s-1) of this passive half
};

/**
 * One synchronization: an ACTIVE event and the PASSIVE event it drives.
 *
 * Every transition of the CTMC is one of these. The active half sets the rate;
 * the passive half is where the job lands, weighted by the routing probability.
 * A LOCAL passive half means the active event moves no job out of its node --
 * a phase change, a reneging job leaving the system, a server failing.
 */
template <class T>
struct Sync {
    SyncEvent<T> active, passive;
};

/**
 * Port of `MNetwork.refreshSync`: the synchronization list.
 *
 * The ORDER matters as much as the content: the generator adds rates in this
 * order, and while addition is commutative in exact arithmetic it is not in
 * floating point, so a reordered list perturbs the last digits of every
 * reported metric.
 *
 * @param impatience_classes (station x class) true where reneging is declared
 * @param breakdown_nodes    1-based node indices with a server breakdown
 * @param sn the refreshed network struct
 */
template <class T>
std::vector<Sync<T>> refresh_sync(
    const NetworkStruct<T>& sn,
    const std::vector<std::vector<bool>>& impatience_classes = std::vector<std::vector<bool>>(),
    const std::vector<std::size_t>& breakdown_nodes = std::vector<std::size_t>()) {
    const std::size_t R = sn.nclasses;
    const std::size_t local = sn.nodes.size() + 1;  // the dummy passive node
    const T one = num_traits<T>::from_int(1);
    std::vector<Sync<T>> sync;

    for (std::size_t ind = 1; ind <= sn.nodes.size(); ++ind) {
        const NodeDef& nd = sn.nodes[ind - 1];
        const std::size_t ist = nd.station;
        for (std::size_t r = 1; r <= R; ++r) {
            // A phase-change action exists only for a multi-phase service:
            // with one phase there is no internal transition to make.
            if (ist != 0 && sn.phases_of(ist, r) > 1) {
                Sync<T> s;
                s.active = SyncEvent<T>{EventType::PHASE, ind, r, one};
                s.passive = SyncEvent<T>{EventType::LOCAL, local, r, one};
                sync.push_back(s);
            }
            if (ist != 0 && impatience_classes.size() >= ist &&
                impatience_classes[ist - 1].size() >= r && impatience_classes[ist - 1][r - 1]) {
                Sync<T> s;
                s.active = SyncEvent<T>{EventType::RENEGE, ind, r, one};
                s.passive = SyncEvent<T>{EventType::LOCAL, local, r, one};
                sync.push_back(s);
            }
            if (ist != 0) {
                const typename std::map<std::size_t, RetrialParam<T>>::const_iterator rit =
                    sn.retrialparam.find(ist);
                if (rit != sn.retrialparam.end() && rit->second.retrial_proc.size() >= r &&
                    !rit->second.retrial_proc[r - 1].disabled) {
                    Sync<T> s;
                    s.active = SyncEvent<T>{EventType::RETRY, ind, r, one};
                    s.passive = SyncEvent<T>{EventType::LOCAL, local, r, one};
                    sync.push_back(s);
                }
            }
            // Failure and repair are properties of the SERVER, not of a class,
            // so exactly one pair is emitted per station rather than one per
            // class -- hence the r == 1 guard.
            if (ist != 0 && r == 1) {
                bool has_bd = false;
                for (std::size_t b = 0; b < breakdown_nodes.size(); ++b)
                    if (breakdown_nodes[b] == ind) { has_bd = true; break; }
                if (has_bd) {
                    Sync<T> f, rp;
                    f.active = SyncEvent<T>{EventType::FAILURE, ind, r, one};
                    f.passive = SyncEvent<T>{EventType::LOCAL, local, r, one};
                    sync.push_back(f);
                    rp.active = SyncEvent<T>{EventType::REPAIR, ind, r, one};
                    rp.passive = SyncEvent<T>{EventType::LOCAL, local, r, one};
                    sync.push_back(rp);
                }
            }
            // A polling station needs a SWITCH action per buffer whose
            // entering switchover is a real (non-immediate) walk.
            if (ist != 0 && sn.stations[ist - 1].sched == SchedStrategy::POLLING) {
                const PollingInfo<T> pinfo = polling_info(sn, ind);
                if (pinfo.valid && pinfo.has_sw[r - 1]) {
                    Sync<T> s2;
                    s2.active = SyncEvent<T>{EventType::SWITCH, ind, r, one};
                    s2.passive = SyncEvent<T>{EventType::LOCAL, local, r, one};
                    sync.push_back(s2);
                }
            }
            if (!nd.stateful) continue;
            // A stateful Fork emits no departure sync: the atomic multi-branch
            // emission is a fork firing synchronization instead.
            if (nd.nodetype == NodeType::Fork) continue;

            // A CACHE READ IS ITS OWN ACTION, not a routing decision. The read
            // consults the contents, rewrites them and switches the job into the
            // hit or the miss class, so it moves no job between nodes and its
            // passive half is the dummy.
            //
            // THE READ CLASS THEREFORE EMITS NO DEPARTURE HERE. `refresh_routing`
            // resolves its unresolved cache split to a uniform half-half so the
            // VISIT equations have a number; the reference leaves the same two
            // entries NaN, and `ceil(NaN) > 0` is false, so no departure sync is
            // built from them. Reproducing the half-half as a synchronization
            // would make the sample path decide hit against miss by a coin
            // instead of by reading the cache, which is what it did.
            bool cache_read_class = false;
            if (nd.nodetype == NodeType::Cache) {
                const typename std::map<std::size_t, CacheParam<T>>::const_iterator ci =
                    sn.nodeparam.find(ind);
                if (ci != sn.nodeparam.end()) {
                    if (r - 1 < ci->second.pread.size() && !ci->second.pread[r - 1].empty()) {
                        Sync<T> s;
                        s.active = SyncEvent<T>{EventType::READ, ind, r, one};
                        s.passive = SyncEvent<T>{EventType::READ, local, r, one};
                        sync.push_back(s);
                    }
                    cache_read_class =
                        r - 1 < ci->second.hitclass.size() && ci->second.hitclass[r - 1] != 0;
                }
            }
            // A stateful Transition emits one server phase-change action per
            // MODE, not per class (the JAR gates the same way on the first
            // class); its cls slot carries the mode, as `after_event_transition`
            // reads it.
            if (nd.nodetype == NodeType::Transition && r == 1) {
                const typename std::map<std::size_t, TransitionParam<T>>::const_iterator ti =
                    sn.transparam.find(ind);
                if (ti != sn.transparam.end()) {
                    for (std::size_t m = 1; m <= ti->second.nmodes; ++m) {
                        Sync<T> s;
                        s.active = SyncEvent<T>{EventType::PHASE, ind, m, one};
                        s.passive = SyncEvent<T>{EventType::LOCAL, local, m, one};
                        sync.push_back(s);
                    }
                }
            }
            if (cache_read_class) continue;

            const std::size_t isf = sn.stateful_index(ind);
            for (std::size_t jnd = 1; jnd <= sn.nodes.size(); ++jnd) {
                if (!sn.nodes[jnd - 1].stateful) continue;
                const std::size_t jsf = sn.stateful_index(jnd);
                for (std::size_t s = 1; s <= R; ++s) {
                    const T p = sn.rt((isf - 1) * R + (r - 1), (jsf - 1) * R + (s - 1));
                    if (num_traits<T>::to_double(p) <= 0) continue;
                    Sync<T> ns;
                    ns.active = SyncEvent<T>{EventType::DEP, ind, r, one};
                    ns.passive = SyncEvent<T>{EventType::ARV, jnd, s, p};
                    // SDR reaches `rt` as the uniform placeholder `refresh_routing`
                    // writes, so the SUPPORT of the mask is right and the values
                    // are not: the pair exists exactly where a link does, and the
                    // generator replaces the probability state by state. This is
                    // the reference's `rtmask = rtfun(emptystate, emptystate)`,
                    // which likewise keeps every connected pair.
                    if (nd.routing.size() >= s && nd.routing[s - 1] == RoutingStrategy::SDR) {
                        ns.passive.statedep = true;
                        ns.passive.rt_row = (isf - 1) * R + (r - 1);
                        ns.passive.rt_col = (jsf - 1) * R + (s - 1);
                    }
                    sync.push_back(ns);
                }
            }
        }
    }
    return sync;
}

/**
 * Port of `State.afterFJEvent`: fire ONE entry of the fork firing list.
 *
 * A fork firing is atomic across several nodes -- it consumes the parent at the
 * fork and places one sibling at each branch head in the same instant -- so
 * unlike every ordinary transition it cannot be decomposed into an active half
 * and a passive half. It therefore takes the whole network state, exactly as an
 * SPN global synchronization does.
 *
 * THE TWO ENABLING CONDITIONS.
 *
 *   1. The fork holds at least one class-r parent.
 *
 *   2. This entry's tag is the LOWEST FREE tag for this (fork, class). A tag is
 *      free when its auxiliary classes have zero occupancy NETWORK-WIDE, which is
 *      why the test scans every stateful node and not just the branches: a
 *      sibling in transit is still outstanding. Without the canonical choice
 *      every firing would produce one successor per free tag, all of them
 *      relabellings of each other, and the chain would carry a factorial number
 *      of duplicate states.
 *
 * The emission is applied SEQUENTIALLY over a growing outcome list rather than
 * branch-by-branch into one state, because that is what handles the three cases
 * a single pass would get wrong: two branches sharing a head node, `weight > 1`
 * repeated emissions on one branch, and a non-exponential sibling service whose
 * phase-entry mixture makes one arrival into several outcomes.
 */
template <class T>
GlobalOutcome<T> after_fj_event(const NetworkStruct<T>& sn, const FjSync<T>& e,
                                const NetState<T>& gl) {
    const std::size_t R = sn.nclasses;
    GlobalOutcome<T> out;
    const std::size_t isf_f = sn.stateful_index(e.fork);
    if (isf_f == 0) return out;
    const std::vector<T>& fs = gl.local[isf_f - 1];
    if (fs.size() < R) return out;
    if (num_traits<T>::to_double(fs[fs.size() - R + e.cls - 1]) < 1) return out;

    // Tag occupancy, network-wide.
    const std::size_t B = e.auxall.size();
    if (B == 0 || e.tag == 0 || e.tag > e.auxall[0].size()) return out;
    const std::size_t Tt = e.auxall[0].size();
    std::vector<double> nglobal(R, 0.0);
    for (std::size_t isf = 1; isf <= sn.stateful_nodes.size(); ++isf) {
        const std::pair<T, std::vector<T>> mg =
            to_marginal_aggr(sn, sn.stateful_nodes[isf - 1], gl.local[isf - 1]);
        for (std::size_t r = 0; r < R; ++r) {
            const double v = num_traits<T>::to_double(mg.second[r]);
            // A Source encodes its reservoir as an infinite marginal; adding it
            // would make every tag look occupied. An auxiliary class is never
            // generated by a Source, so skipping the non-finite entries cannot
            // hide a real sibling.
            if (std::isfinite(v)) nglobal[r] += v;
        }
    }
    std::vector<double> occ(Tt, 0.0);
    for (std::size_t t = 0; t < Tt; ++t)
        for (std::size_t b = 0; b < B; ++b) occ[t] += nglobal[e.auxall[b][t] - 1];
    if (occ[e.tag - 1] > 0) return out;                       // this tag is in use
    for (std::size_t t = 0; t + 1 < e.tag; ++t)
        if (occ[t] == 0) return out;                          // a lower tag is free

    NetState<T> seed = gl;
    seed.local[isf_f - 1][seed.local[isf_f - 1].size() - R + e.cls - 1] -=
        num_traits<T>::from_int(1);

    std::vector<NetState<T>> partials(1, seed);
    std::vector<T> partprob(1, num_traits<T>::from_int(1));
    // The emission list, one entry per sibling. `weightlink` is filled only when
    // the fork sends different counts down different links; the interleave below
    // reproduces `repmat(1:B,1,w)` exactly in the uniform case, so a plain fork
    // walks the order it always did.
    std::vector<std::size_t> emissions;
    if (e.weightlink.empty()) {
        const std::size_t w = e.weight == 0 ? 1 : e.weight;
        for (std::size_t rep = 0; rep < w; ++rep)
            for (std::size_t b = 0; b < B; ++b) emissions.push_back(b);
    } else {
        std::size_t wmax = 0;
        for (std::size_t b = 0; b < e.weightlink.size(); ++b)
            if (e.weightlink[b] > wmax) wmax = e.weightlink[b];
        for (std::size_t rep = 1; rep <= wmax; ++rep)
            for (std::size_t b = 0; b < B; ++b)
                if (b < e.weightlink.size() && e.weightlink[b] >= rep) emissions.push_back(b);
    }
    for (std::size_t ei = 0; ei < emissions.size(); ++ei) {
        const std::size_t b = emissions[ei];
        const std::size_t bh = e.branchheads[b];
        const std::size_t isf_b = sn.stateful_index(bh);
        if (isf_b == 0) return GlobalOutcome<T>();
        const std::size_t a = e.auxclasses[b];
        std::vector<NetState<T>> nextp;
        std::vector<T> nextq;
        for (std::size_t pp = 0; pp < partials.size(); ++pp) {
            const EventOutcome<T> arv =
                after_event(sn, bh, partials[pp].local[isf_b - 1], EventType::ARV, a);
            if (arv.space.empty()) return GlobalOutcome<T>();  // blocked: no firing
            for (std::size_t io = 0; io < arv.space.size(); ++io) {
                NetState<T> ns = partials[pp];
                ns.local[isf_b - 1] = arv.space[io];
                nextp.push_back(ns);
                nextq.push_back(io < arv.prob.size() ? T(partprob[pp] * arv.prob[io])
                                                     : partprob[pp]);
            }
        }
        partials.swap(nextp);
        partprob.swap(nextq);
    }

    out.space = partials;
    for (std::size_t i = 0; i < partials.size(); ++i) {
        out.rate.push_back(num_traits<T>::from_double(lang::GlobalConstants::Immediate));
        out.prob.push_back(T(partprob[i] * e.prob));
        out.completion.push_back(true);
    }
    return out;
}

}  // namespace qn
}  // namespace line

#endif  // LINE_LANG_QN_STATE_EVENTS_H
