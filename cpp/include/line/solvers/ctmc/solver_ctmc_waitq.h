/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Port of `solver_ctmc_fcr_waitq.m`: the reachability-built generator of a model
 * whose finite capacity region applies WAITQ.
 *
 * WHY THIS IS A SEPARATE GENERATOR AND NOT A FILTER. Under DROP a refused job is
 * destroyed, so the chain never occupies a forbidden state and censoring the
 * enumerated space IS the censored chain -- that is `solver_ctmc_fcr.h`. Under
 * WAITQ the refused job LEAVES its upstream station and parks in a per-region
 * FIFO of (class, destination) tokens that sits outside every station. That FIFO
 * is state, it is owned by the region rather than by any node, and no filter on
 * the per-node state space can represent it. The state here is therefore
 * augmented as [per-node states, buf_1, ..., buf_F] and the transition relation
 * is rebuilt around it.
 *
 * THE FOUR RULES THAT MAKE WAITQ WHAT IT IS, all from the JMT reference:
 *
 *   1. RELEASE IS STRICTLY HEAD-OF-LINE. After every transition that frees
 *      capacity, tokens leave in FIFO order and a head that still does not fit
 *      blocks the whole queue behind it, even where a token further back would
 *      fit. Releasing the first token that fits instead would be a different
 *      discipline with a different mean.
 *   2. A FRESH ARRIVAL IS NOT QUEUED BEHIND THE FIFO. It is admitted whenever the
 *      constraints permit, so it overtakes a head stuck on a different
 *      constraint. The FIFO gates only the jobs already in it.
 *   3. THE RELEASE CASCADE IS PART OF THE SAME TRANSITION. Freeing one slot can
 *      release a token whose admission frees another slot, and so on; the chain
 *      jumps straight to the settled state. Splitting the cascade into separate
 *      transitions would invent intermediate states with a residence time the
 *      model does not have.
 *   4. PARKED JOBS ARE IN NO STATION AND IN NO REGION. Station queue lengths
 *      exclude them, which is the JMT report convention; they are visible only
 *      through `ctmc_waitq_parked`, and any population accounting has to add
 *      them back by hand.
 *
 * A CLASS SWITCH INSIDE ONE REGION IS THE SUBTLE CASE. Such a hop leaves the
 * region occupancy unchanged in total but moves one job between the per-class
 * counts, so it can violate a per-class cap that the pre-transition state
 * satisfied. It cannot be tested before the transition either, because the
 * departure frees the old class's slot first. It is therefore deferred: the
 * re-entry is resolved after the release cascade has settled, and parked at the
 * TAIL if it still does not fit, since it was refused after every token already
 * in the queue.
 */
#ifndef LINE_SOLVERS_CTMC_SOLVER_CTMC_WAITQ_H
#define LINE_SOLVERS_CTMC_SOLVER_CTMC_WAITQ_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <deque>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "line/api/mc/ctmc_solve.h"
#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/solver_feature_sets.h"
#include "line/lang/qn/state.h"
#include "line/lang/qn/state_events.h"
#include "line/solvers/ctmc/solver_ctmc.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_fcr.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace ctmc {

/**
 * One augmented state: the network state, plus the token FIFO of every region.
 *
 * A token is `(dest_node - 1) * K + cls`, both 1-based, which is the reference's
 * encoding: a parked job knows the class it will enter in and the node it was on
 * its way to when it was refused.
 */
template <class T>
struct WaitqState {
    NetState<T> net;
    std::vector<std::vector<std::size_t>> buf;
};

/** The chain the WAITQ walk produces, alongside the FIFOs its states carry. */
template <class T>
struct WaitqResult {
    CtmcResult<T> chain;                                     ///< Q, the net halves, rates, filt
    std::vector<std::vector<std::vector<std::size_t>>> buf;  ///< per state, per region
};

namespace waitq_detail {

/**
 * The admission constraints of one region, flattened out of the Region struct
 * once so the inner loop is not re-deriving them at every candidate state.
 *
 * The predicate this feeds is the one `ctmc_region_admissible` applies to a
 * whole state; it is restated per region because WAITQ has to ask about ONE
 * region at a time, on a candidate population vector that no state realizes yet.
 */
template <class T>
struct RegionCaps {
    std::vector<bool> member;   ///< per station, 0-based
    std::vector<double> ccap;   ///< per class, -1 = unbounded
    double gcap = -1.0;
    double memcap = -1.0;
    std::vector<double> size;   ///< per class memory footprint
    Matrix<T> A;                ///< the linear constraint A x <= b
    std::vector<T> b;
    std::vector<bool> iswaitq;  ///< per class; false = this class is DROP here
    std::size_t lmax = 0;       ///< FIFO length bound
};

template <class T>
std::vector<RegionCaps<T>> extract_caps(const NetworkStruct<T>& sn) {
    const std::size_t M = sn.stations.size(), K = sn.nclasses;
    std::vector<RegionCaps<T>> out;
    out.reserve(sn.regions.size());
    for (std::size_t f = 0; f < sn.regions.size(); ++f) {
        const typename NetworkStruct<T>::Region& rg = sn.regions[f];
        RegionCaps<T> c;
        c.member.assign(M, false);
        c.ccap.assign(K, -1.0);
        c.size.assign(K, 1.0);
        c.iswaitq.assign(K, false);
        for (std::size_t i = 0; i < M; ++i) {
            if (i >= rg.members.size() || !rg.members[i]) continue;
            c.member[i] = true;
            // The caps are replicated on every member row, so taking the
            // tightest is both correct and safe against a hand-built struct.
            for (std::size_t k = 0; k < K && k < rg.cap[i].size(); ++k)
                if (rg.cap[i][k] != -1.0)
                    c.ccap[k] = c.ccap[k] == -1.0 ? rg.cap[i][k] : std::min(c.ccap[k], rg.cap[i][k]);
            if (K < rg.cap[i].size() && rg.cap[i][K] != -1.0)
                c.gcap = c.gcap == -1.0 ? rg.cap[i][K] : std::min(c.gcap, rg.cap[i][K]);
            if (i < rg.maxmem.size() && rg.maxmem[i] != -1.0)
                c.memcap = c.memcap == -1.0 ? rg.maxmem[i] : std::min(c.memcap, rg.maxmem[i]);
        }
        for (std::size_t k = 0; k < K; ++k) {
            if (k < rg.size.size()) c.size[k] = num_traits<T>::to_double(rg.size[k]);
            c.iswaitq[k] = k < rg.rule.size() && rg.rule[k] != DropStrategy::DROP;
        }
        c.A = rg.lincon_A;
        c.b = rg.lincon_b;
        out.push_back(c);
    }
    return out;
}

/** The reference's `violates`: does this per-class population break region f. */
template <class T>
bool violates(const RegionCaps<T>& c, const std::vector<double>& x) {
    double total = 0, memory = 0;
    for (std::size_t k = 0; k < x.size(); ++k) {
        total += x[k];
        memory += x[k] * c.size[k];
        if (c.ccap[k] != -1.0 && x[k] > c.ccap[k] + 1e-9) return true;
    }
    if (c.gcap != -1.0 && total > c.gcap + 1e-9) return true;
    if (c.memcap != -1.0 && memory > c.memcap + 1e-9) return true;
    for (std::size_t row = 0; row < c.A.rows() && row < c.b.size(); ++row) {
        double lhs = 0;
        for (std::size_t k = 0; k < x.size() && k < c.A.cols(); ++k)
            lhs += num_traits<T>::to_double(c.A(row, k)) * x[k];
        if (lhs > num_traits<T>::to_double(c.b[row]) + 1e-9) return true;
    }
    return false;
}

/**
 * The reference's `regionAggr`: the per-class population a region currently
 * holds, summed over its member stations.
 *
 * `to_marginal_aggr` rather than `to_marginal`, for the reason the reference's
 * arrival branch uses it: it leaves a preempted job out of the count, and
 * counting one twice would refuse an admission the region has room for.
 */
template <class T>
std::vector<double> region_aggr(const NetworkStruct<T>& sn, const RegionCaps<T>& c,
                                const NetState<T>& net) {
    const std::size_t K = sn.nclasses;
    std::vector<double> x(K, 0.0);
    for (std::size_t i = 0; i < c.member.size(); ++i) {
        if (!c.member[i]) continue;
        const std::size_t ind = sn.node_of_station(i + 1);
        const std::size_t isf = sn.stateful_index(ind);
        if (isf == 0) continue;
        const std::pair<T, std::vector<T>> mg = qn::to_marginal_aggr(sn, ind, net.local[isf - 1]);
        for (std::size_t k = 0; k < K && k < mg.second.size(); ++k) {
            const double v = num_traits<T>::to_double(mg.second[k]);
            // A Source reports an infinite reservoir, which describes the
            // encoding and not an occupancy; adding it would make every region
            // permanently full.
            if (std::isfinite(v)) x[k] += v;
        }
    }
    return x;
}

/** The flattened key of an augmented state, for exact index lookup. */
template <class T>
std::vector<double> waitq_key(const WaitqState<T>& ws) {
    std::vector<double> key = ctmc_detail::state_key(ws.net);
    for (std::size_t f = 0; f < ws.buf.size(); ++f) {
        key.push_back(-3.0);  // a separator no node block and no token emits
        for (std::size_t j = 0; j < ws.buf[f].size(); ++j)
            key.push_back(static_cast<double>(ws.buf[f][j]));
    }
    return key;
}

/**
 * The FIFO length bound of each region: at most every job that could be parked
 * at once.
 *
 * A CLOSED CLASS IS BOUNDED BY ITS WHOLE CHAIN, not by its own population, since
 * a job may switch into the class before being refused. An open class has no
 * population to bound it and takes the state-space cutoff instead, which is what
 * makes the FIFO finite at all -- and what makes the walk TRUNCATE rather than
 * refuse when the bound is reached, exactly as the reference does.
 */
template <class T>
void resolve_lmax(const NetworkStruct<T>& sn, const std::vector<std::size_t>& cutoff,
                  std::vector<RegionCaps<T>>& caps) {
    const std::size_t K = sn.nclasses;
    const std::vector<double> N = sn.njobs();
    std::vector<std::size_t> tokbound(K, 0);
    for (std::size_t r = 0; r < K; ++r) {
        bool any = false;
        for (std::size_t f = 0; f < caps.size(); ++f)
            if (caps[f].iswaitq[r]) any = true;
        if (!any) continue;
        double chainpop = 0;
        bool finite = true;
        for (std::size_t c = 0; c < sn.nchains; ++c) {
            bool holds = false;
            for (std::size_t a = 0; a < sn.inchain[c].size(); ++a)
                if (sn.inchain[c][a] == r + 1) holds = true;
            if (!holds) continue;
            for (std::size_t a = 0; a < sn.inchain[c].size(); ++a) {
                const double nj = N[sn.inchain[c][a] - 1];
                if (std::isfinite(nj))
                    chainpop += nj;
                else
                    finite = false;
            }
        }
        tokbound[r] = finite ? static_cast<std::size_t>(chainpop)
                             : (r < cutoff.size() ? cutoff[r] : 0);
    }
    for (std::size_t f = 0; f < caps.size(); ++f) {
        std::size_t s = 0;
        for (std::size_t r = 0; r < K; ++r)
            if (caps[f].iswaitq[r]) s += tokbound[r];
        caps[f].lmax = s;
    }
}

/**
 * The widest local row stateful node `ind` admits, which is the width every
 * state of that node carries in the ENUMERATED space.
 *
 * WHY A WALK NEEDS THIS AND THE ENUMERATION DOES NOT. `space_generator` builds
 * each node's rows for every marginal and RIGHT-ALIGNS the narrow ones into the
 * widest it saw (the reference's `fromMarginalBounds`), so an FCFS station whose
 * buffer is one slot wide when empty carries the buffer of its FULLEST marginal
 * in every state. A walk seeded from `default_init_state` gets the row at the
 * NATURAL width of the initial marginal instead -- for an empty FCFS station a
 * one-slot buffer -- and nothing ever widens it: the arrival branch here places a
 * job in an existing empty slot and refuses when there is none, where MATLAB's
 * `afterEventStation` prepends a column first. The station then holds nservers+1
 * jobs and no more; for a CLOSED class the refusal is a BLOCK rather than a loss,
 * so no successor is emitted at all and the upstream departure never fires.
 *
 * MEASURED on Think -> Q1 -> Q2 with three jobs and an UNBOUNDED region, which
 * must reproduce the region-less solution exactly: X came out 0.935049 against
 * 0.962806. That is not a perturbation, it is the exact answer for a different
 * model -- the same network truncated to two jobs per queue, which is what a
 * one-slot buffer in front of one server means.
 *
 * The width of a node's row is a property of that node alone: for an ordered
 * (class-tag) buffer it grows with the TOTAL jobs held and not with how they
 * split across classes, and for every other encoding it is constant. So one call
 * to `from_marginal_node` at a maximal admissible marginal settles it. The
 * descent to a smaller total is not defensive padding: `from_marginal_node`
 * returns NO rows for a marginal the station cannot hold, and the joint bound
 * (`cap` against the sum of the per-class bounds) can be met by a marginal that
 * some other constraint inside the handler still rejects.
 *
 * `ssa::serial_detail::max_row_width` is this same computation, which the serial
 * simulator needs for this same reason -- it walks the encoding too. The two are
 * stated twice because a CTMC generator has no business including the
 * simulator's header, and they must be changed together; the one place that
 * would hold a single copy is `State` itself.
 */
template <class T>
std::size_t max_row_width(const NetworkStruct<T>& sn, std::size_t ind,
                          const std::vector<std::size_t>& cutoff) {
    const std::size_t R = sn.nclasses;
    const std::size_t ist = sn.nodes[ind - 1].station;
    std::vector<std::size_t> ph(R, 1);
    if (ist != 0)
        for (std::size_t r = 0; r < R; ++r) ph[r] = sn.phasessz_of(ist, r + 1);

    // A Source and a stateful non-station (a Cache, a Join) hold a per-class
    // count and their local variables, both of fixed width, so the empty
    // marginal already gives the final width.
    if (ist == 0 || sn.stations[ist - 1].nodetype == NodeType::Source) {
        std::vector<T> row;
        if (!qn::from_marginal_node_first(sn, ind, std::vector<std::size_t>(R, 0), ph, row))
            throw UnsupportedError("SolverCTMC: node '" + sn.nodes[ind - 1].name +
                                   "' admits no state at all, so the WAITQ walk cannot start");
        return row.size();
    }

    // The per-class bound: the class population when closed, the state-space
    // cutoff when open, never above the station's own per-class capacity.
    const std::vector<double> N = sn.njobs();
    std::vector<std::size_t> bound(R, 0);
    for (std::size_t r = 0; r < R; ++r) {
        double b = std::isfinite(N[r]) ? N[r]
                                       : static_cast<double>(r < cutoff.size() ? cutoff[r] : 0);
        const double cc = sn.classcap[ist - 1][r];
        if (cc < b) b = cc;
        bound[r] = b > 0 ? static_cast<std::size_t>(b) : 0;
    }
    std::size_t total = 0;
    for (std::size_t r = 0; r < R; ++r) total += bound[r];
    const double tcap = sn.cap[ist - 1];
    if (std::isfinite(tcap) && tcap < static_cast<double>(total))
        total = tcap > 0 ? static_cast<std::size_t>(tcap) : 0;

    for (std::size_t t = total + 1; t-- > 0;) {
        std::vector<std::size_t> n(R, 0);
        std::size_t left = t;
        for (std::size_t r = 0; r < R && left > 0; ++r) {
            n[r] = std::min(left, bound[r]);
            left -= n[r];
        }
        if (left > 0) continue;  // this total does not fit the per-class bounds
        const std::vector<std::vector<T>> rows = qn::from_marginal_node(sn, ind, n, ph);
        if (!rows.empty()) return rows[0].size();
    }
    throw UnsupportedError("SolverCTMC: station '" + sn.stations[ist - 1].name +
                           "' admits no state at all, so the WAITQ walk cannot start");
}

/**
 * The initial network state, padded to the encoding width of every node.
 *
 * The LEFT pad is not a convention chosen here: it is what `space_generator`
 * does to the narrow rows, and every slicer in `state_events.h` measures its
 * blocks from the RIGHT-hand end, so a right pad would shift the server block
 * and silently decode the wrong queue.
 */
template <class T>
NetState<T> wide_init_state(const NetworkStruct<T>& sn,
                            const std::vector<std::size_t>& cutoff) {
    NetState<T> init;
    if (!analyzer_detail::default_init_state(sn, init))
        throw UnsupportedError(
            "SolverCTMC: the model's initial state admits no state at some node; check the class "
            "populations against their reference stations");
    for (std::size_t f = 0; f < sn.stateful_nodes.size() && f < init.local.size(); ++f) {
        const std::size_t w = max_row_width(sn, sn.stateful_nodes[f], cutoff);
        if (init.local[f].size() >= w) continue;
        init.local[f].insert(init.local[f].begin(), w - init.local[f].size(),
                             num_traits<T>::from_int(0));
    }
    return init;
}

/**
 * One tentative successor, before the release cascade has run on it.
 *
 * `arv` and `dep_isf` carry the per-class rate bookkeeping through the cascade
 * rather than letting the caller add it: the cascade can DROP a branch (a full
 * FIFO), and a rate charged before that happens would count a transition that
 * never fired.
 */
template <class T>
struct Emission {
    NetState<T> net;
    std::vector<std::vector<std::size_t>> buf;
    T w = num_traits<T>::from_int(0);
    std::vector<std::pair<std::size_t, std::size_t>> arv;  ///< (stateful index, class), 1-based
    std::size_t dep_isf = 0, dep_cls = 0;                  ///< 0 when the active half is not a DEP
    std::size_t pf = 0;                                    ///< deferred re-entry region, 1-based
    std::size_t pcls = 0, pdest = 0;
    bool pwaitq = false;
};

/** One item of the release cascade: a tentative state and how it got there. */
template <class T>
struct Work {
    NetState<T> net;
    std::vector<std::vector<std::size_t>> buf;
    T prob;
    std::vector<std::pair<std::size_t, std::size_t>> arv;
    std::size_t dep_isf = 0, dep_cls = 0;
    std::size_t pf = 0, pcls = 0, pdest = 0;
    bool pwaitq = false;
};

/**
 * One SETTLED successor of an augmented state: where the transition goes, what
 * it weighs, and the rate statistics it carries.
 *
 * "Settled" means the release cascade has run to its fixed point, so a caller
 * never sees an intermediate state with a token that could still have left. The
 * weight is the rate times the probability of this particular leaf of the
 * cascade tree, which is why one synchronization can appear several times.
 */
template <class T>
struct Successor {
    std::size_t sync = 0;  ///< index into the `sync` list the enumeration was given
    WaitqState<T> next;
    T w = num_traits<T>::from_int(0);
    std::vector<std::pair<std::size_t, std::size_t>> arv;  ///< (stateful, class), 1-based
    std::size_t dep_isf = 0, dep_cls = 0;                  ///< 0 when nothing departed
};

/**
 * Every transition enabled in one augmented state, with its settled successor.
 *
 * THIS IS THE WAITQ TRANSITION RELATION, and it is written once. The chain
 * builder walks it breadth-first to fill a generator; the SSA serial engine
 * draws one of the returned moves at each step. Two copies of the four rules in
 * this file's header would be two chances for the region discipline to differ
 * between the exact solver and the simulator while both looked correct.
 *
 * `sync` and `caps` are passed in rather than rebuilt because both are loop
 * invariants of either caller, and `resolve_lmax` has to have run on `caps`
 * already -- an unbounded FIFO admits every token and the walk would not
 * terminate.
 */
template <class T>
void waitq_successors(const NetworkStruct<T>& sn, const std::vector<Sync<T>>& sync,
                      const std::vector<RegionCaps<T>>& caps, const WaitqState<T>& ws,
                      std::vector<Successor<T>>& out) {
    const std::size_t K = sn.nclasses;
    const std::size_t F = caps.size();
    const std::size_t local = sn.nodes.size() + 1;
    const NetState<T>& net = ws.net;
    const std::vector<std::vector<std::size_t>>& buf = ws.buf;
    out.clear();

    std::vector<std::vector<double>> xf(F);
    for (std::size_t f = 0; f < F; ++f) xf[f] = region_aggr(sn, caps[f], net);

    // State-dependent routing, evaluated once at the state the transitions leave
    // from, exactly as the ordinary generator does. A region and an SDR
    // subnetwork are independent constructs, so this relation must carry eq. (10)
    // too -- reading the uniform placeholder in `sn.rt` here would answer a
    // uniformly-routed model under an SDR name at every model that has both.
    Matrix<T> rt_now;
    const bool sdr = sn.has_sdr_routing();
    if (sdr) rt_now = qn::rt_state(sn, net.local);

    /*
     * The reference's `emit`: run the release cascade to a fixed point, then
     * record one successor per settled leaf.
     *
     * A WORK LIST and not a recursion because a release can have several
     * outcomes -- an arrival that splits over the phase it starts service in --
     * so the cascade is a tree, every leaf is a distinct successor of the SAME
     * transition, and each carries its own share of the probability.
     */
    const auto emit = [&](std::size_t a, const Emission<T>& e) {
        if (num_traits<T>::to_double(e.w) <= 0) return;
        std::deque<Work<T>> work;
        Work<T> w0;
        w0.net = e.net;
        w0.buf = e.buf;
        w0.prob = num_traits<T>::from_int(1);
        w0.arv = e.arv;
        w0.dep_isf = e.dep_isf;
        w0.dep_cls = e.dep_cls;
        w0.pf = e.pf;
        w0.pcls = e.pcls;
        w0.pdest = e.pdest;
        w0.pwaitq = e.pwaitq;
        work.push_back(w0);

        while (!work.empty()) {
            const Work<T> it = work.front();
            work.pop_front();

            bool progressed = false;
            for (std::size_t f = 0; f < F && !progressed; ++f) {
                if (it.buf[f].empty()) continue;
                const std::size_t tok = it.buf[f][0];
                const std::size_t dest = (tok - 1) / K + 1;
                const std::size_t r = (tok - 1) % K + 1;
                std::vector<double> x = region_aggr(sn, caps[f], it.net);
                x[r - 1] += 1.0;
                // Head of line: this region's FIFO stays blocked even where a
                // token behind the head would fit.
                if (violates(caps[f], x)) continue;
                const std::size_t isf_d = sn.stateful_index(dest);
                if (isf_d == 0) continue;
                const qn::EventOutcome<T> od =
                    qn::after_event(sn, dest, it.net.local[isf_d - 1], EventType::ARV, r);
                if (od.space.empty()) continue;
                for (std::size_t id = 0; id < od.space.size(); ++id) {
                    if (num_traits<T>::to_double(od.prob[id]) <= 0) continue;
                    Work<T> nx = it;
                    nx.net.local[isf_d - 1] = od.space[id];
                    nx.buf[f].erase(nx.buf[f].begin());
                    nx.prob = T(it.prob * od.prob[id]);
                    nx.arv.push_back(std::make_pair(isf_d, r));
                    work.push_back(nx);
                }
                progressed = true;
            }
            if (progressed) continue;

            if (it.pf != 0) {
                // The cascade has settled, so the deferred class switch can be
                // tested against a region occupancy that already reflects the
                // departure which freed the old class's slot.
                const std::size_t f = it.pf - 1;
                std::vector<double> x = region_aggr(sn, caps[f], it.net);
                x[it.pcls - 1] += 1.0;
                Work<T> nx = it;
                nx.pf = 0;
                if (violates(caps[f], x)) {
                    // A DROP class loses the switching job outright; a WAITQ one
                    // joins the tail.
                    if (it.pwaitq) {
                        if (nx.buf[f].size() >= caps[f].lmax) continue;
                        nx.buf[f].push_back((it.pdest - 1) * K + it.pcls);
                    }
                    work.push_back(nx);
                    continue;
                }
                const std::size_t isf_d = sn.stateful_index(it.pdest);
                bool admitted = false;
                if (isf_d != 0) {
                    const qn::EventOutcome<T> od = qn::after_event(
                        sn, it.pdest, it.net.local[isf_d - 1], EventType::ARV, it.pcls);
                    for (std::size_t id = 0; id < od.space.size(); ++id) {
                        if (num_traits<T>::to_double(od.prob[id]) <= 0) continue;
                        Work<T> ny = nx;
                        ny.net.local[isf_d - 1] = od.space[id];
                        ny.prob = T(it.prob * od.prob[id]);
                        ny.arv.push_back(std::make_pair(isf_d, it.pcls));
                        work.push_back(ny);
                        admitted = true;
                    }
                }
                if (!admitted) {
                    // The region had room but the destination's own state does
                    // not admit the job -- its station capacity, say -- so it
                    // parks rather than vanishing.
                    if (nx.buf[f].size() >= caps[f].lmax) continue;
                    nx.buf[f].push_back((it.pdest - 1) * K + it.pcls);
                    work.push_back(nx);
                }
                continue;
            }

            const T contrib = T(e.w * it.prob);
            if (num_traits<T>::to_double(contrib) <= 0) continue;
            Successor<T> s;
            s.sync = a;
            s.next.net = it.net;
            s.next.buf = it.buf;
            s.w = contrib;
            s.arv = it.arv;
            s.dep_isf = it.dep_isf;
            s.dep_cls = it.dep_cls;
            out.push_back(s);
        }
    };

    for (std::size_t a = 0; a < sync.size(); ++a) {
        const Sync<T>& sy = sync[a];
        const std::size_t node_a = sy.active.node;
        const std::size_t isf_a = sn.stateful_index(node_a);
        if (isf_a == 0) continue;
        const std::size_t cls_a = sy.active.cls;
        const std::size_t stat_a = sn.nodes[node_a - 1].station;

        const qn::EventOutcome<T> oa =
            qn::after_event(sn, node_a, net.local[isf_a - 1], sy.active.event, cls_a);
        for (std::size_t ia = 0; ia < oa.space.size(); ++ia) {
            const T rate = oa.rate[ia];
            if (num_traits<T>::to_double(rate) <= 0) continue;

            const std::size_t node_p = sy.passive.node;
            if (node_p == local) {
                Emission<T> e;
                e.net = net;
                e.net.local[isf_a - 1] = oa.space[ia];
                e.buf = buf;
                e.w = rate;
                emit(a, e);
                continue;
            }
            const std::size_t isf_p = sn.stateful_index(node_p);
            if (isf_p == 0) continue;
            const std::size_t cls_p = sy.passive.cls;
            const std::size_t stat_p = sn.nodes[node_p - 1].station;
            const T w = T(rate * (sy.passive.statedep
                                      ? rt_now(sy.passive.rt_row, sy.passive.rt_col)
                                      : sy.passive.prob));
            if (num_traits<T>::to_double(w) <= 0) continue;
            const std::size_t dep_isf = sy.active.event == EventType::DEP ? isf_a : 0;

            // REGION ENTRY is the only place the caps are tested on the way
            // in: the passive station is inside the region and the active
            // node is not, so this transition raises the region's occupancy.
            std::size_t blockedf = 0, droppedf = 0;
            if (sy.passive.event == EventType::ARV && stat_p != 0) {
                for (std::size_t f = 0; f < F; ++f) {
                    if (!caps[f].member[stat_p - 1]) continue;
                    if (stat_a != 0 && caps[f].member[stat_a - 1]) continue;
                    std::vector<double> xn = xf[f];
                    xn[cls_p - 1] += 1.0;
                    if (!violates(caps[f], xn)) continue;
                    if (caps[f].iswaitq[cls_p - 1])
                        blockedf = f + 1;
                    else
                        droppedf = f + 1;
                    break;
                }
            }

            if (droppedf != 0) {
                // DROP: the refused job is destroyed, so only the active
                // half applies. It still departed, which is why the loss
                // shows up downstream as ArvR - Tput.
                Emission<T> e;
                e.net = net;
                e.net.local[isf_a - 1] = oa.space[ia];
                e.buf = buf;
                e.w = w;
                e.dep_isf = dep_isf;
                e.dep_cls = cls_a;
                emit(a, e);
                continue;
            }

            // A class switch between two members of the same region: see the
            // file header. Deferred, not tested here.
            std::size_t switchf = 0;
            if (blockedf == 0 && sy.passive.event == EventType::ARV && cls_p != cls_a &&
                stat_a != 0 && stat_p != 0)
                for (std::size_t f = 0; f < F; ++f)
                    if (caps[f].member[stat_a - 1] && caps[f].member[stat_p - 1]) {
                        switchf = f + 1;
                        break;
                    }

            if (switchf != 0) {
                Emission<T> e;
                e.net = net;
                e.net.local[isf_a - 1] = oa.space[ia];
                e.buf = buf;
                e.w = w;
                e.dep_isf = dep_isf;
                e.dep_cls = cls_a;
                e.pf = switchf;
                e.pcls = cls_p;
                e.pdest = node_p;
                e.pwaitq = caps[switchf - 1].iswaitq[cls_p - 1];
                emit(a, e);
                continue;
            }

            if (blockedf != 0) {
                // The FIFO bound is a TRUNCATION of an open class's
                // unbounded queue, not an error: the reference drops the
                // transition, and so does this. It is the same kind of
                // approximation the state-space cutoff already is.
                if (buf[blockedf - 1].size() >= caps[blockedf - 1].lmax) continue;
                Emission<T> e;
                e.net = net;
                e.net.local[isf_a - 1] = oa.space[ia];
                e.buf = buf;
                e.buf[blockedf - 1].push_back((node_p - 1) * K + cls_p);
                e.w = w;
                e.dep_isf = dep_isf;
                e.dep_cls = cls_a;
                emit(a, e);
                continue;
            }

            const std::vector<T>& psrc = node_p == node_a ? oa.space[ia] : net.local[isf_p - 1];
            const qn::EventOutcome<T> op = qn::after_event(sn, node_p, psrc, sy.passive.event, cls_p);
            // No successor at the passive node is a BLOCKED arrival, the
            // reference's true-BAS branch. Every model that declares true
            // blocking has already been refused by name, so what remains is
            // a destination refusing the job outside any region, and the
            // transition simply does not fire.
            for (std::size_t ip = 0; ip < op.space.size(); ++ip) {
                if (num_traits<T>::to_double(op.prob[ip]) <= 0) continue;
                Emission<T> e;
                e.net = net;
                e.net.local[isf_a - 1] = oa.space[ia];
                e.net.local[isf_p - 1] = op.space[ip];
                e.buf = buf;
                e.w = T(w * op.prob[ip]);
                e.dep_isf = dep_isf;
                e.dep_cls = cls_a;
                if (sy.passive.event == EventType::ARV)
                    e.arv.push_back(std::make_pair(isf_p, cls_p));
                emit(a, e);
            }
        }
    }
}

}  // namespace waitq_detail

/** True when the model declares a region that applies anything other than DROP. */
template <class T>
bool ctmc_has_waitq_region(const NetworkStruct<T>& sn) {
    for (std::size_t f = 0; f < sn.regions.size(); ++f)
        for (std::size_t r = 0; r < sn.regions[f].rule.size(); ++r)
            if (sn.regions[f].rule[r] != DropStrategy::DROP) return true;
    return false;
}

/**
 * The combinations the reference gates, plus the two this port cannot represent.
 *
 * The first three are the reference's own: each needs a semantics for what a
 * parked token means that the reference declines to define. The last two are
 * this port's, and both are about state that does not exist here rather than
 * about semantics.
 */
template <class T>
void ctmc_check_waitq_support(const NetworkStruct<T>& sn) {
    if (!sn.transparam.empty())
        throw UnsupportedError(
            "SolverCTMC: a WAITQ finite capacity region is not supported together with stochastic "
            "Petri net transitions; a firing is atomic across its arcs and has no single "
            "destination to park a refused token against");
    if (sn.has_fork() || !sn.fj.empty())
        throw UnsupportedError(
            "SolverCTMC: a WAITQ finite capacity region is not supported together with fork-join; "
            "a forked task refused entry would park without its siblings, and the join has no rule "
            "for a sibling that never arrived");
    for (std::size_t i = 0; i < sn.nodes.size(); ++i) {
        const qn::NodeDef& nd = sn.nodes[i];
        for (std::size_t r = 0; r < nd.routing.size(); ++r) {
            const lang::RoutingStrategy rs = nd.routing[r];
            if (rs == lang::RoutingStrategy::PROB || rs == lang::RoutingStrategy::RAND ||
                rs == lang::RoutingStrategy::DISABLED)
                continue;
            throw UnsupportedError(
                "SolverCTMC: a WAITQ finite capacity region is not supported together with the "
                "state-dependent routing at node '" +
                nd.name +
                "'; the destination a token parks against is fixed when the job is refused, and a "
                "state-dependent choice would have to be re-decided on release");
        }
    }
    for (std::size_t f = 0; f < sn.regions.size(); ++f)
        for (std::size_t r = 0; r < sn.regions[f].rule.size(); ++r) {
            const DropStrategy d = sn.regions[f].rule[r];
            if (d == DropStrategy::DROP || d == DropStrategy::WAITQ) continue;
            throw UnsupportedError(
                "SolverCTMC: finite capacity region " + std::to_string(f + 1) +
                " applies a blocking rule other than DROP or WAITQ to class " +
                std::to_string(r + 1) +
                "; BAS, BBS and RSRD hold the job AT ITS SERVER, which needs a blocked-server "
                "marker column in the node state that this port does not allocate");
        }
    // TRUE BAS BLOCKING IS ORTHOGONAL TO THE REGION RULE and is refused on its
    // own terms. The reference holds the completing job at its server by setting
    // the last column of the station's state row; `refresh_local_vars` here
    // reserves that shared column for a Cache, a breakdown or a polling
    // controller and never for a BAS marker, so there is no column to set and no
    // handler that would release it. Voiding the departure instead is NOT
    // equivalent: it frees the server to re-serve the same job, which for
    // exponential service redraws the completion and understates throughput.
    for (std::size_t i = 0; i < sn.droprule.size() && i < sn.stations.size(); ++i)
        for (std::size_t r = 0; r < sn.droprule[i].size(); ++r) {
            const DropStrategy d = sn.droprule[i][r];
            if (d != DropStrategy::BAS && d != DropStrategy::BBS && d != DropStrategy::RSRD)
                continue;
            throw UnsupportedError(
                "SolverCTMC: station '" + sn.stations[i].name +
                "' applies true blocking (BAS/BBS/RSRD) to class " + std::to_string(r + 1) +
                ", which holds the completing job at its server; this port allocates no "
                "blocked-server marker in the node state, and voiding the departure instead would "
                "let the server re-serve the same job and understate throughput");
        }
}

/**
 * Port of the reachability walk of `solver_ctmc_fcr_waitq.m`.
 *
 * The space is WALKED and not enumerated, for the reason an SPN is: no
 * population marginal produces a state whose region FIFO is non-empty, so a
 * lattice enumeration would emit only the empty-buffer states and every blocking
 * transition would land outside the space.
 */
template <class T>
WaitqResult<T> solver_ctmc_waitq(const NetworkStruct<T>& sn, const CtmcOptions& opt) {
    ctmc_check_waitq_support(sn);

    const std::size_t K = sn.nclasses;
    const std::size_t NF = sn.stateful_nodes.size();
    const T zero = num_traits<T>::from_int(0);

    const std::vector<Sync<T>> sync = refresh_sync(sn);
    std::vector<waitq_detail::RegionCaps<T>> caps = waitq_detail::extract_caps(sn);
    const std::size_t F = caps.size();
    const std::vector<std::size_t> cutoff = analyzer_detail::resolve_cutoff(sn, opt);
    waitq_detail::resolve_lmax(sn, cutoff, caps);

    // AT THE ENCODING WIDTH, not at the width of the initial marginal: the rows
    // never grow during the walk, so a station seeded with the narrow buffer of
    // its empty state can never hold more than nservers+1 jobs. See
    // `max_row_width`.
    WaitqState<T> init;
    init.net = waitq_detail::wide_init_state(sn, cutoff);
    init.buf.assign(F, std::vector<std::size_t>());
    for (std::size_t f = 0; f < F; ++f)
        if (waitq_detail::violates(caps[f], waitq_detail::region_aggr(sn, caps[f], init.net)))
            throw InputError(
                "SolverCTMC: the initial state violates the constraints of finite capacity region " +
                std::to_string(f + 1) + "; the region cannot hold the model's initial population");

    std::vector<WaitqState<T>> states;
    states.push_back(init);
    std::map<std::vector<double>, std::size_t> index;
    index[waitq_detail::waitq_key(init)] = 0;
    std::deque<std::size_t> frontier;
    frontier.push_back(0);

    // Triplets rather than a matrix, because the state count is not known until
    // the walk has finished and a growing dense matrix would be recopied at
    // every doubling.
    std::vector<std::size_t> ta, ti, tj;
    std::vector<T> tv;
    std::vector<std::vector<std::vector<T>>> arv, dep;
    arv.push_back(std::vector<std::vector<T>>(NF, std::vector<T>(K, zero)));
    dep.push_back(std::vector<std::vector<T>>(NF, std::vector<T>(K, zero)));

    // Index an augmented state, appending it and its rate blocks when new.
    const auto state_index = [&](const WaitqState<T>& ws) -> std::size_t {
        const std::vector<double> key = waitq_detail::waitq_key(ws);
        const typename std::map<std::vector<double>, std::size_t>::const_iterator it =
            index.find(key);
        if (it != index.end()) return it->second;
        if (states.size() >= opt.state_max)
            throw UnsupportedError(
                "SolverCTMC: the WAITQ state space exceeds the cap of " +
                std::to_string(opt.state_max) +
                " states; the region FIFO multiplies the space by its own occupancy");
        states.push_back(ws);
        arv.push_back(std::vector<std::vector<T>>(NF, std::vector<T>(K, zero)));
        dep.push_back(std::vector<std::vector<T>>(NF, std::vector<T>(K, zero)));
        index[key] = states.size() - 1;
        frontier.push_back(states.size() - 1);
        return states.size() - 1;
    };

    // The walk. `frontier` is consumed from the front, so states are expanded in
    // discovery order and the numbering matches the reference's own
    // breadth-first construction. The transition relation itself is
    // `waitq_successors`, shared with the SSA serial engine so that the region
    // discipline cannot differ between the exact solver and the simulator.
    std::vector<waitq_detail::Successor<T>> succ;
    while (!frontier.empty()) {
        const std::size_t s = frontier.front();
        frontier.pop_front();
        const WaitqState<T> ws = states[s];  // by value: `states` grows below
        waitq_detail::waitq_successors(sn, sync, caps, ws, succ);
        for (std::size_t e = 0; e < succ.size(); ++e) {
            const waitq_detail::Successor<T>& su = succ[e];
            const std::size_t dst = state_index(su.next);
            ta.push_back(su.sync);
            ti.push_back(s);
            tj.push_back(dst);
            tv.push_back(su.w);
            if (su.dep_isf != 0) dep[s][su.dep_isf - 1][su.dep_cls - 1] += su.w;
            for (std::size_t q = 0; q < su.arv.size(); ++q)
                arv[s][su.arv[q].first - 1][su.arv[q].second - 1] += su.w;
        }
    }

    const std::size_t n = states.size();
    WaitqResult<T> res;
    res.chain.Q = Matrix<T>(n, n, zero);
    res.chain.space.reserve(n);
    res.buf.reserve(n);
    for (std::size_t s = 0; s < n; ++s) {
        res.chain.space.push_back(states[s].net);
        res.buf.push_back(states[s].buf);
    }
    res.chain.arv_rates = arv;
    res.chain.dep_rates = dep;
    if (opt.keep_filtration) res.chain.filt.assign(sync.size(), Matrix<T>(n, n, zero));
    for (std::size_t e = 0; e < tv.size(); ++e) {
        res.chain.Q(ti[e], tj[e]) += tv[e];
        if (opt.keep_filtration) res.chain.filt[ta[e]](ti[e], tj[e]) += tv[e];
    }
    // A refused arrival that leaves every node unchanged is a self-loop, and the
    // diagonal has to cancel it exactly as it does in the default generator.
    make_infgen(res.chain.Q);
    return res;
}

/**
 * The mean number of parked jobs per class, over a stationary law.
 *
 * IT IS IN NO QLen. A parked job is in no station and in no region, so
 * `solver_ctmc_avg_from_pi` cannot see it and the model's population is
 * conserved only once this is added back. That is the JMT report convention,
 * not an omission.
 */
template <class T>
std::vector<T> ctmc_waitq_parked(const NetworkStruct<T>& sn, const WaitqResult<T>& r,
                                 const std::vector<T>& pi) {
    const std::size_t K = sn.nclasses;
    std::vector<T> out(K, num_traits<T>::from_int(0));
    for (std::size_t s = 0; s < r.buf.size() && s < pi.size(); ++s)
        for (std::size_t f = 0; f < r.buf[s].size(); ++f)
            for (std::size_t j = 0; j < r.buf[s][f].size(); ++j) {
                const std::size_t cls = (r.buf[s][f][j] - 1) % K + 1;
                out[cls - 1] += pi[s];
            }
    return out;
}

/** A solved WAITQ model: the usual CTMC solution, plus what the FIFOs hold. */
template <class T>
struct WaitqSolution {
    CtmcSolution<T> sol;
    std::vector<std::vector<std::vector<std::size_t>>> buf;
    std::vector<T> parked;  ///< mean parked jobs per class
};

/**
 * Build the WAITQ chain, solve it, and map it onto the same means every other
 * CTMC path reports.
 *
 * NO WEAKLY-CONNECTED-COMPONENT STEP, unlike `solver_ctmc_analyzer`. That step
 * exists because the lattice enumeration emits states the dynamics cannot reach;
 * this walk starts at the initial state and applies the same handlers the
 * generator does, so every state it holds is reachable by construction and
 * restricting to a component could only remove states the model does occupy.
 */
template <class T>
WaitqSolution<T> solver_ctmc_waitq_analyzer(const NetworkStruct<T>& sn, const CtmcOptions& opt) {
    check_method(opt.method);
    WaitqSolution<T> out;
    const WaitqResult<T> r = solver_ctmc_waitq(sn, opt);
    out.sol.cutoff = analyzer_detail::resolve_cutoff(sn, opt);
    out.sol.pi = mc::ctmc_solve(r.chain.Q);
    out.sol.avg = solver_ctmc_avg_from_pi(sn, r.chain, out.sol.pi);
    out.sol.chain = r.chain;
    out.sol.actualmethod = opt.method;
    out.buf = r.buf;
    out.parked = ctmc_waitq_parked(sn, r, out.sol.pi);
    return out;
}

/** A CTMC solve routed to whichever path the model's region rules require. */
template <class T>
struct CtmcAnySolution {
    CtmcSolution<T> sol;
    std::vector<std::vector<std::vector<std::size_t>>> buf;  ///< empty off the WAITQ path
    std::vector<T> parked;  ///< mean parked jobs per class; empty off the WAITQ path
    bool waitq = false;
};

/**
 * The entry point a caller who does not know which path a model needs should
 * use: pick the WAITQ walk when a region asks for anything other than DROP, and
 * the lattice analyzer otherwise.
 *
 * IT IS A SEPARATE ENTRY POINT AND NOT A CHANGE TO `solver_ctmc_run_analyzer`. That
 * function is the DEFAULT path, and a WAITQ region reaching it is a defect it
 * must keep reporting: its generator carries no token buffer, so it would solve
 * the region as DROP. The dispatch is therefore placed here, where both paths
 * are in scope, and the default path keeps refusing by name.
 *
 * `ctmc_check_support` runs on both branches. `ctmc_check_waitq_support` covers
 * fork-join and true blocking on its own and needs no entry for class or joint
 * dependence: this generator never scales a rate itself, it walks on whatever
 * `EventOutcome::rate` the shared `qn::after_event` -> `after_event_station`
 * dispatcher hands back (state_events.h:1591, :1675), and `cd_factor` is
 * folded into that rate INSIDE the dispatcher (state_events.h:831-832, :1298)
 * for every caller alike. `solver_ctmc.h`'s default generator calls the same
 * `qn::after_event`, so the two paths scale identically by construction, not
 * by coincidence. GREPPING THIS FILE FOR `cd_factor` FINDS NOTHING BY DESIGN --
 * that is not a missing call, it is this generator having no rate arithmetic
 * of its own to put one in. Trace the dispatcher, not the grep, before
 * concluding otherwise. MATLAB has the identical shape for the identical
 * reason: `solver_ctmc_fcr_waitq.m` computes every rate through
 * `State.afterEventHashed`, which calls the same `State.afterEvent` the
 * default generator uses, so it too carries no cdscaling/jdscaling text of
 * its own.
 */
template <class T>
CtmcAnySolution<T> solver_ctmc_analyzer_any(const NetworkStruct<T>& sn, const CtmcOptions& opt) {
    ctmc_check_support(sn);
    // The WAITQ branch short-circuits below and would otherwise never see the
    // gate that solver_ctmc_analyzer applies on the default path.
    qn::feature_gate("SolverCTMC", qn::ctmc_feature_set(opt.method), sn);
    CtmcAnySolution<T> out;
    if (!ctmc_has_waitq_region(sn)) {
        out.sol = solver_ctmc_analyzer(sn, opt);
        return out;
    }
    const WaitqSolution<T> w = solver_ctmc_waitq_analyzer(sn, opt);
    out.sol = w.sol;
    out.buf = w.buf;
    out.parked = w.parked;
    out.waitq = true;
    return out;
}

/** Solve on whichever path applies and format, mirroring `solver_ctmc_run_analyzer`. */
template <class T>
mva::AvgResult<T> solver_ctmc_run_analyzer_any(const NetworkStruct<T>& sn, const CtmcOptions& opt) {
    return solver_ctmc_avg_table(sn, solver_ctmc_analyzer_any(sn, opt).sol, opt.method);
}

}  // namespace ctmc
}  // namespace line

#endif  // LINE_SOLVERS_CTMC_SOLVER_CTMC_WAITQ_H
