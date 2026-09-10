/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Port of the MATLAB `+State` package: the encoding that turns a station's
 * state row into marginal job counts.
 *
 * WHAT A STATE ROW IS. For a station the row is [buffer | server | vars]:
 * `nvars(ind,:)` local-variable slots at the end, `sum(phasesz)` server slots
 * before them (one per class-phase pair, laid out by `phaseshift`), and
 * whatever remains at the front is the buffer. The buffer encoding is NOT
 * uniform across disciplines -- FCFS/LCFS store a CLASS TAG per waiting
 * position, so class-r jobs are counted by matching the tag, while SIRO,
 * POLLING, SEPT, LEPT and SRPT store a per-class COUNT in column r. Reading
 * one as the other silently produces a plausible number, which is why the
 * discipline switch below is transcribed case by case rather than collapsed.
 *
 * THE EXT SENTINEL. `to_marginal` returns nir = +Inf for a Source, and that is
 * deliberate: a Source is an infinite reservoir and the value describes the
 * ENCODING, not a queue length. A caller that treats it as a queue length gets
 * Inf, which is exactly the defect that reached SolverCTMC's averagers in
 * MATLAB (see `_kb/07-cross-language-parity.md`, the Source-row section). Any
 * consumer must test `is_source` before using nir, not clamp the Inf away.
 */
#ifndef LINE_LANG_QN_STATE_H
#define LINE_LANG_QN_STATE_H

#include <algorithm>
#include <map>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/api/pfqn/pfqn_comb_common.h"
#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/polling_info.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace qn {

/** What `State.toMarginal` returns for one station and one state row. */
template <class T>
struct Marginal {
    T ni = T();                            ///< total jobs in the station
    std::vector<T> nir;                    ///< jobs per class
    std::vector<T> sir;                    ///< jobs in service per class
    std::vector<std::vector<T>> kir;       ///< jobs in service per class and phase
};

namespace state_detail {

/** True when the discipline stores a per-class COUNT in buffer column r. */
inline bool buffer_is_per_class_count(SchedStrategy s) {
    switch (s) {
        case SchedStrategy::SIRO:
        case SchedStrategy::POLLING:
        case SchedStrategy::SEPT:
        case SchedStrategy::LEPT:
        case SchedStrategy::SRPT:
        case SchedStrategy::SRPTPRIO:
            return true;
        default:
            return false;
    }
}

/**
 * True when the discipline stores a CLASS TAG per waiting position. HOL is
 * MATLAB's FCFSPRIO, so the two are the same enumerator and appear once.
 */
inline bool buffer_is_class_tag(SchedStrategy s) {
    switch (s) {
        case SchedStrategy::FCFS:
        case SchedStrategy::HOL:
        case SchedStrategy::LCFS:
        case SchedStrategy::LCFSPRIO:
            return true;
        default:
            return false;
    }
}

/**
 * True when the buffer interleaves [class, phase] pairs, so only the even
 * positions carry class tags. The preemptive families do this to remember the
 * phase every preempted job was interrupted in.
 */
inline bool buffer_is_tag_phase_pairs(SchedStrategy s) {
    switch (s) {
        // The whole preempt-resume / preempt-independent family, exactly as
        // `toMarginal` groups it. PI restarts a preempted job from its entry
        // phase and PR resumes it in place, but BOTH must record the phase per
        // waiting job, so they share the paired encoding.
        case SchedStrategy::FCFSPI:
        case SchedStrategy::FCFSPIPRIO:
        case SchedStrategy::FCFSPR:
        case SchedStrategy::FCFSPRPRIO:
        case SchedStrategy::LCFSPI:
        case SchedStrategy::LCFSPIPRIO:
        case SchedStrategy::LCFSPR:
        case SchedStrategy::LCFSPRPRIO:
            return true;
        default:
            return false;
    }
}

}  // namespace state_detail

/**
 * Port of `State.toMarginal` for a STATION, one state row at a time.
 *
 * @param sn        the network struct
 * @param ist       station index (1-based, as elsewhere in NetworkStruct)
 * @param state_i   the station's state row
 * @param phasesz   per-class phase counts
 * @param phaseshift per-class offset into the server block
 * @param nvar      width of the trailing local-variable block, 0 when absent
 * @return the marginal counts; nir is +Inf for every class at a Source
 */
template <class T>
Marginal<T> to_marginal(const NetworkStruct<T>& sn, std::size_t ist,
                        const std::vector<T>& state_i, const std::vector<std::size_t>& phasesz,
                        const std::vector<std::size_t>& phaseshift, std::size_t nvar = 0) {
    const std::size_t R = sn.nclasses;
    if (ist == 0 || ist > sn.stations.size())
        throw InputError("to_marginal: station index " + std::to_string(ist) + " is out of range");

    const T zero = num_traits<T>::from_int(0);
    Marginal<T> m;
    m.nir.assign(R, zero);
    m.sir.assign(R, zero);
    std::size_t maxph = 1;
    for (std::size_t r = 0; r < R; ++r) maxph = std::max(maxph, phasesz[r]);
    m.kir.assign(R, std::vector<T>(maxph, zero));

    // A Join of an FJ-augmented struct: the row is a bare per-class count of jobs
    // WAITING to synchronize. Nothing is in service and nothing is in a phase,
    // which is why `sir` and `kir` stay zero rather than mirroring `nir`.
    const std::size_t jnd = sn.node_of_station(ist);
    if (sn.isfjaugmented && jnd != 0 && sn.nodes[jnd - 1].nodetype == NodeType::Join &&
        state_i.size() >= R) {
        for (std::size_t r = 0; r < R; ++r) m.nir[r] = state_i[state_i.size() - R + r];
        return m;
    }

    // AN ORDER-INDEPENDENT ROW HAS NO SERVER BLOCK. PAS and OI encode the
    // station as the ORDERED LIST of the class indices present, one slot per
    // job the buffer can hold, so its width is the capacity and not
    // sum(phasesz). Slicing a server block off it read the last list positions
    // as phase occupancies, and where the capacity is SMALLER than the class
    // count the width test below fired outright: `to_marginal: state row is
    // narrower than the server block it declares` on pas_compatibility_5class
    // (five classes, capacity three). The arm further down decodes the list; it
    // needs neither `srv0` nor the per-phase sum.
    const SchedStrategy sched = sn.stations[ist - 1].sched;
    const bool ordered_list = (sched == SchedStrategy::PAS || sched == SchedStrategy::OI);

    // [buffer | server | vars]: slice from the RIGHT, since only the buffer
    // width varies with the discipline.
    std::size_t srvw = 0;
    for (std::size_t r = 0; r < R; ++r) srvw += phasesz[r];
    if (!ordered_list && state_i.size() < nvar + srvw)
        throw InputError("to_marginal: state row is narrower than the server block it declares");
    const std::size_t srv0 = ordered_list ? 0 : state_i.size() - nvar - srvw;
    const std::size_t bufw = srv0;

    if (!ordered_list) {
        for (std::size_t r = 0; r < R; ++r) {
            for (std::size_t k = 0; k < phasesz[r]; ++k) {
                const T v = state_i[srv0 + phaseshift[r] + k];
                m.kir[r][k] = v;
                m.sir[r] += v;
            }
        }
    }

    if (sched == SchedStrategy::EXT) {
        // Infinite reservoir: a statement about the encoding, not a queue
        // length. Consumers must branch on the station being a Source.
        //
        // AN EXACT TYPE HAS NO INFINITY. Rational is a field of quotients of
        // integers, so building one from a double Inf throws outright ("Cannot
        // convert a non-finite number to an integer") -- and the sentinel is
        // built for every Source row whether or not anybody reads it, which
        // killed the whole exact run inside `after_event_station_dep`
        // (state_events.h): it calls this on the Source before its own EXT
        // branch ever looks at sir/kir. Clamp to MaxInt there, the same clamp
        // `from_marginal_node` applies to an infinite server count for the same
        // reason -- still absurd as a queue length, so a consumer that forgot
        // its `is_source` test is still visibly wrong rather than plausibly
        // wrong, but representable. Types that DO carry an infinity keep it, so
        // the double path is untouched.
        const T ext = num_traits<T>::is_exact
                          ? num_traits<T>::from_double(GlobalConstants::MaxInt)
                          : num_traits<T>::from_double(
                                std::numeric_limits<double>::infinity());
        for (std::size_t r = 0; r < R; ++r) m.nir[r] = ext;
    } else if (sched == SchedStrategy::PAS || sched == SchedStrategy::OI) {
        // `State.toMarginal`'s PAS arm (toMarginal.m:105-130), which this port
        // did not have: the row is the ORDERED LIST of class indices and there
        // is no server block to slice at all, so the generic sum above read the
        // last list POSITION as a phase occupancy. A station holding jobs then
        // reported a queue length of about zero in every consumer of this
        // decode -- `ctmc_state_space_aggr` and with it SolverCTMC's averages,
        // SolverSSA's serial analyzer, and the reward and transient analyzers.
        //
        // IN SERVICE IS NOT "IN THE SERVER" HERE but "receiving a positive rate
        // increment": Delta mu over the prefix ending at that position, exactly
        // as the reference computes it. That is what makes an order-independent
        // station's utilization a quantity rather than a slot count.
        for (std::size_t r = 0; r < R; ++r) {
            m.nir[r] = zero;
            m.sir[r] = zero;
            for (std::size_t k = 0; k < maxph; ++k) m.kir[r][k] = zero;
        }
        const std::size_t w = state_i.size() > nvar ? state_i.size() - nvar : 0;
        std::vector<std::size_t> clist;
        for (std::size_t b = 0; b < w; ++b) {
            const long tag = static_cast<long>(num_traits<T>::to_double(state_i[b]) + 0.5);
            if (tag >= 1 && static_cast<std::size_t>(tag) <= R) {
                m.nir[tag - 1] += num_traits<T>::from_int(1);
                clist.push_back(static_cast<std::size_t>(tag));
            }
        }
        const typename std::map<std::size_t,
                                typename NetworkStruct<T>::PasParam>::const_iterator pit =
            sn.pasparam.find(ist);
        if (pit != sn.pasparam.end() && pit->second.svc_rate_fun) {
            T muprev = zero;
            for (std::size_t p = 0; p < clist.size(); ++p) {
                const std::vector<std::size_t> prefix(clist.begin(), clist.begin() + p + 1);
                const T mucur = pit->second.svc_rate_fun(prefix);
                if (num_traits<T>::to_double(T(mucur - muprev)) > 0) {
                    m.sir[clist[p] - 1] += num_traits<T>::from_int(1);
                    m.kir[clist[p] - 1][0] += num_traits<T>::from_int(1);
                }
                muprev = mucur;
            }
        } else {
            // No rate function declared: the reference falls back to "every job
            // present is in service", which is the OI reading of an empty swap
            // graph.
            m.sir = m.nir;
            for (std::size_t r = 0; r < R; ++r) m.kir[r][0] = m.nir[r];
        }
    } else if (state_detail::buffer_is_class_tag(sched)) {
        for (std::size_t r = 0; r < R; ++r) {
            T waiting = zero;
            for (std::size_t b = 0; b < bufw; ++b)
                if (state_i[b] == num_traits<T>::from_int(static_cast<long>(r + 1))) waiting += num_traits<T>::from_int(1);
            m.nir[r] = T(m.sir[r] + waiting);
        }
    } else if (state_detail::buffer_is_tag_phase_pairs(sched)) {
        if (bufw > 1) {
            for (std::size_t r = 0; r < R; ++r) {
                T waiting = zero;
                for (std::size_t b = 0; b < bufw; b += 2)  // even positions are the class tags
                    if (state_i[b] == num_traits<T>::from_int(static_cast<long>(r + 1))) waiting += num_traits<T>::from_int(1);
                m.nir[r] = T(m.sir[r] + waiting);
            }
        } else {
            m.nir = m.sir;
        }
    } else if (state_detail::buffer_is_per_class_count(sched)) {
        for (std::size_t r = 0; r < R; ++r)
            m.nir[r] = bufw >= R ? T(m.sir[r] + state_i[r]) : m.sir[r];
    } else {
        // INF, PS, DPS, GPS and the rest: everything present is in service.
        m.nir = m.sir;
    }

    // A Place is a token container: the buffer/server split in its encoding is
    // an artifact of a transition FIRE relocating surviving tokens, so fold the
    // buffer slot back in for the INF-family disciplines it uses. Without this
    // the token count collapses to the oscillating server slot while the
    // dynamics stay correct -- the measurement moves, the model does not.
    if (sn.stations[ist - 1].nodetype == NodeType::Place) {
        switch (sched) {
            case SchedStrategy::INF:
            case SchedStrategy::PS:
            case SchedStrategy::PSPRIO:
            case SchedStrategy::DPS:
            case SchedStrategy::DPSPRIO:
            case SchedStrategy::GPS:
            case SchedStrategy::GPSPRIO:
            case SchedStrategy::LPS:
                if (bufw >= R)
                    for (std::size_t r = 0; r < R; ++r) m.nir[r] = T(m.sir[r] + state_i[r]);
                break;
            default:
                break;
        }
        for (std::size_t r = 0; r < R; ++r)
            if (sn.disabled[ist - 1][r]) {
                for (std::size_t k = 0; k < phasesz[r]; ++k) m.kir[r][k] = zero;
                m.sir[r] = zero;
            }
    } else {
        for (std::size_t r = 0; r < R; ++r)
            if (sn.disabled[ist - 1][r]) {
                m.nir[r] = zero;
                for (std::size_t k = 0; k < phasesz[r]; ++k) m.kir[r][k] = zero;
                m.sir[r] = zero;
            }
    }

    m.ni = zero;
    for (std::size_t r = 0; r < R; ++r) m.ni += m.nir[r];
    return m;
}

/**
 * Total jobs held by every STATION at one network state, indexed by station.
 *
 * A Source is an infinite reservoir whose `to_marginal` reports +Inf, which
 * describes the encoding and not a queue length, so it is reported as zero here
 * exactly as `ctmc_state_space_aggr` does. A station with no stateful node
 * holds nothing.
 */
template <class T>
std::vector<double> station_populations(const NetworkStruct<T>& sn,
                                        const std::vector<std::vector<T>>& local) {
    const std::size_t M = sn.stations.size(), K = sn.nclasses;
    std::vector<double> n(M, 0.0);
    for (std::size_t ist = 1; ist <= M; ++ist) {
        const std::size_t isf = sn.stateful_of_station(ist);
        if (isf == 0 || isf > local.size()) continue;
        if (sn.stations[ist - 1].nodetype == NodeType::Source) continue;
        std::vector<std::size_t> ph(K, 1), shift(K, 0);
        std::size_t w = 0;
        for (std::size_t k = 0; k < K; ++k) {
            ph[k] = sn.phasessz_of(ist, k + 1);
            shift[k] = w;
            w += ph[k];
        }
        const Marginal<T> m =
            to_marginal(sn, ist, local[isf - 1], ph, shift, sn.nvars_of(sn.node_of_station(ist)));
        n[ist - 1] = num_traits<T>::to_double(m.ni);
    }
    return n;
}

/**
 * Port of `sn.rtfun`: the routing over the stateful nodes AT ONE STATE.
 *
 * `rt` is a constant matrix because almost every routing strategy is; SDR is
 * not, and eq. (10) of Krzesinski (1987) makes the split out of the entry
 * centre a function of the branch and subnetwork populations. This rebuilds the
 * node-level table with the SDR rows re-evaluated at `local` and eliminates the
 * stateless nodes through the SAME stochastic complement `rt` is built with, so
 * a model with a Router between the centres is complemented identically either
 * way.
 *
 * The residual mass returns the customer to the departure centre -- the busy
 * form of waiting of Sec. 2.5 -- and a branch closed by its own population bound
 * simply receives nothing, because `pfqn_sdrprob` returns zero there.
 *
 * `sub_sdr` returns zero off the class diagonal: SDR does not switch class, so
 * only the (r, r) block of the entry row carries mass.
 *
 * @param sn    the refreshed network struct
 * @param local per-stateful-node state rows, the `NetState::local` of the state
 * @return the (nstateful * nclasses) square routing table at that state
 */
template <class T>
Matrix<T> rt_state(const NetworkStruct<T>& sn, const std::vector<std::vector<T>>& local) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t K = sn.nclasses, I = sn.nodes.size();
    Matrix<T> full = sn.rtnodes;
    if (full.rows() != I * K)
        throw InputError("rt_state: the node-level routing has not been refreshed");
    if (sn.sdr.empty()) return sn.stoch_comp_stateful(full, K);

    const pfqn::SdrCoeff co = pfqn::pfqn_sdrcoeff(sn.sdr);
    const std::vector<double> Pb = pfqn::pfqn_sdrprob(co, station_populations(sn, local));
    const double Ped = pfqn::pfqn_sdrped(Pb);

    // The destination split out of the entry centre, by NODE: a branch is
    // entered at its own entry centre, and the denied customer goes to d.
    std::vector<double> p(I, 0.0);
    for (std::size_t b = 1; b < sn.sdr_nodes.branch.size(); ++b)
        p[sn.sdr_nodes.entryOf[b]] += Pb[b];
    p[sn.sdr_nodes.departure] += Ped;

    for (std::size_t ind = 1; ind <= I; ++ind) {
        const NodeDef& nd = sn.nodes[ind - 1];
        for (std::size_t r = 1; r <= K && r <= nd.routing.size(); ++r) {
            if (nd.routing[r - 1] != RoutingStrategy::SDR) continue;
            for (std::size_t jnd = 1; jnd <= I; ++jnd)
                for (std::size_t s = 1; s <= K; ++s)
                    full((ind - 1) * K + (r - 1), (jnd - 1) * K + (s - 1)) =
                        s == r ? num_traits<T>::from_double(p[jnd - 1]) : zero;
        }
    }
    return sn.stoch_comp_stateful(full, K);
}

/**
 * Port of `State.cartesian`: pair every row of `a` with every row of `b`.
 *
 * An empty `a` is the identity, as in the reference, so a fold over classes can
 * start from nothing. Row order is a's outer, b's inner -- the same order
 * `fromMarginal` relies on when it appends the server block after the buffer.
 */
template <class T>
std::vector<std::vector<T>> cartesian(const std::vector<std::vector<T>>& a,
                                      const std::vector<std::vector<T>>& b) {
    if (a.empty()) return b;
    if (b.empty()) return a;
    std::vector<std::vector<T>> out;
    out.reserve(a.size() * b.size());
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < b.size(); ++j) {
            std::vector<T> row = a[i];
            row.insert(row.end(), b[j].begin(), b[j].end());
            out.push_back(row);
        }
    return out;
}

/**
 * Port of `State.spaceClosedSingle`: the ways to place `n` jobs over `m`
 * phases. `m == 0` yields NO rows, not one empty row -- a class with no service
 * process contributes nothing, and returning an empty row instead would let it
 * multiply the product by one and silently survive the fold.
 */
template <class T>
std::vector<std::vector<T>> space_closed_single(std::size_t m, std::size_t n) {
    std::vector<std::vector<T>> out;
    if (m == 0) return out;
    const std::vector<std::vector<int>> rows =
        pfqn::multichoose_rows(static_cast<int>(m), static_cast<int>(n));
    out.reserve(rows.size());
    for (std::size_t i = 0; i < rows.size(); ++i) {
        std::vector<T> r;
        r.reserve(rows[i].size());
        for (std::size_t j = 0; j < rows[i].size(); ++j)
            r.push_back(num_traits<T>::from_int(rows[i][j]));
        out.push_back(r);
    }
    return out;
}

/**
 * `space_closed_single` with a PER-SLOT bound, the reference's
 * `spaceClosedSingle(M, N, caps)`.
 *
 * The bound is applied BEFORE recursing, which is the whole point: a slot of
 * capacity 0 takes only the zero item, so a class that can occupy one node out
 * of M enumerates M rows rather than binomial(n+M-1, M-1). The unbounded form
 * generates the full lattice and leaves the caller to reject the impossible
 * rows one at a time, at one `from_marginal` call each.
 *
 * `caps` is one entry per slot; an entry of `-1` means unbounded.
 */
template <class T>
void space_closed_single_capped_rec(std::size_t m, long n, const std::vector<long>& caps,
                                    std::size_t off, std::vector<T>& row,
                                    std::vector<std::vector<T>>& out) {
    if (m == 0) {
        if (n == 0) out.push_back(row);
        return;
    }
    long room = 0;
    bool unbounded = false;
    for (std::size_t k = off; k < off + m && k < caps.size(); ++k) {
        if (caps[k] < 0) { unbounded = true; break; }
        room += caps[k];
    }
    if (!unbounded && n > room) return;
    const long here = off < caps.size() ? caps[off] : -1;
    const long hi = here < 0 ? n : std::min<long>(n, here);
    for (long i = 0; i <= hi; ++i) {
        row.push_back(num_traits<T>::from_int(i));
        space_closed_single_capped_rec(m - 1, n - i, caps, off + 1, row, out);
        row.pop_back();
    }
}

/** @see space_closed_single_capped_rec */
template <class T>
std::vector<std::vector<T>> space_closed_single_capped(std::size_t m, std::size_t n,
                                                       const std::vector<long>& caps) {
    std::vector<std::vector<T>> out;
    if (m == 0) return out;
    std::vector<T> row;
    row.reserve(m);
    space_closed_single_capped_rec<T>(m, static_cast<long>(n), caps, 0, row, out);
    return out;
}

/**
 * Port of `matlab/util/multiset_perms.m` on an ASCENDING multiset, ROW ORDER
 * INCLUDED. (Until 2026-08-19 the reference was the vendored `uniqueperms.m`,
 * removed for want of a license grant; the replacement reproduces its listing
 * exactly, so this port is unchanged.)
 *
 * The row order is not an aesthetic choice: row 0 is what `default_init_state`
 * takes as the initial state, so it decides which communicating class a chain
 * made reducible by the swap graph settles in. The reference's own ordering is
 * inconsistent between its two branches and is mirrored here rather than
 * normalised: an all-distinct multiset goes through `perms`, which is REVERSE
 * lexicographic (descending first), while a multiset with a repeat recurses
 * over its unique values in ASCENDING order. Both are reproduced.
 *
 * The caller supplies `vec` sorted ascending, which is what the class-major
 * build in the PAS arm of `from_marginal_core` produces; the recursion erases
 * one element and preserves that order.
 */
inline std::vector<std::vector<std::size_t> > pas_multiset_perms(
    const std::vector<std::size_t>& vec) {
    std::vector<std::vector<std::size_t> > pu;
    if (vec.empty()) return pu;
    std::vector<std::size_t> uvec = vec;
    std::sort(uvec.begin(), uvec.end());
    uvec.erase(std::unique(uvec.begin(), uvec.end()), uvec.end());
    if (uvec.size() == 1) {
        pu.push_back(vec);
        return pu;
    }
    if (uvec.size() == vec.size()) {
        std::vector<std::size_t> p(vec.rbegin(), vec.rend());
        do {
            pu.push_back(p);
        } while (std::prev_permutation(p.begin(), p.end()));
        return pu;
    }
    for (std::size_t i = 0; i < uvec.size(); ++i) {
        std::vector<std::size_t> v = vec;
        for (std::size_t j = 0; j < v.size(); ++j)
            if (v[j] == uvec[i]) {
                v.erase(v.begin() + static_cast<std::ptrdiff_t>(j));
                break;
            }
        const std::vector<std::vector<std::size_t> > tmp = pas_multiset_perms(v);
        for (std::size_t t = 0; t < tmp.size(); ++t) {
            std::vector<std::size_t> row(1, uvec[i]);
            row.insert(row.end(), tmp[t].begin(), tmp[t].end());
            pu.push_back(row);
        }
    }
    return pu;
}

/**
 * Port of `State.fromMarginal` for the station families CTMC enumerates:
 * every local state in which station `ist` holds exactly `n[r]` class-r jobs.
 *
 * COVERED: Queue, Delay, Source and Place under the disciplines whose buffer is
 * either absent (INF/PS/DPS/GPS) or a per-class count. The row is laid out
 * [buffer | server] to match `to_marginal`, and the server block is the
 * cartesian fold of `space_closed_single(phases[r], n[r])` over classes.
 *
 * Retrial stations are covered here too, by the (in-service, orbit) split, and
 * PAS/OI by the ordered class-index list that `to_marginal` already decoded.
 * Transition nodes are NOT: a Transition is stateful but is not a station, so
 * it never carries a station index -- `from_marginal_node` handles it, which is
 * also where the reference puts it.
 *
 * REFUSED BY NAME, because a wrong guess is indistinguishable from a correct
 * one downstream: a station carrying a MAP/MMPP2 arrival process under a
 * discipline whose encoding does not carry the modulating phase.
 */
template <class T>
std::vector<std::vector<T>> from_marginal_core(const NetworkStruct<T>& sn, std::size_t ist,
                                               const std::vector<std::size_t>& n,
                                               const std::vector<std::size_t>& phases) {
    const std::size_t R = sn.nclasses;
    if (ist == 0 || ist > sn.stations.size())
        throw InputError("from_marginal: station index " + std::to_string(ist) + " is out of range");
    if (n.size() != R || phases.size() != R)
        throw InputError("from_marginal: n and phases must have one entry per class");

    const Station<T>& st = sn.stations[ist - 1];
    const SchedStrategy sched = st.sched;

    // The reference's MAP guard. A MAP/MMPP2 arrival process modulates between
    // phases, and only the FCFS encoding carries that phase in the state; under
    // any other discipline the enumeration would drop the modulating chain and
    // return states that cannot represent the process. A Source is exempt: its
    // phase block IS the modulating chain, which the EXT branch builds.
    for (std::size_t r = 0; r < R; ++r) {
        const ProcessType pt = sn.procid(ist, r + 1);
        if ((pt == ProcessType::MAP || pt == ProcessType::MMPP2) &&
            sched != SchedStrategy::FCFS && st.nodetype != NodeType::Source)
            throw UnsupportedError(
                "from_marginal: a MAP/MMPP2 process at a non-FCFS station is not supported; "
                "only the FCFS encoding carries the modulating phase");
    }

    // RETRIAL: no waiting line. An arrival that finds every server busy joins an
    // ORBIT and re-attempts at the retrial rate, so the state is the
    // (in-service, orbit) SPLIT and every admissible split is a distinct state
    // -- including the idle-server ones, which an ordinary queue cannot occupy
    // while jobs wait. A completion does NOT promote from the orbit.
    // The reference's own test is `any(~cellfun(@@isempty, sn.retrialProc(ist,:)))`:
    // an entry exists only once some class actually has a retrial process, so a
    // present-but-empty record must NOT switch the encoding.
    bool is_retrial = false;
    {
        const typename std::map<std::size_t, RetrialParam<T> >::const_iterator rit =
            sn.retrialparam.find(ist);
        if (rit != sn.retrialparam.end())
            for (std::size_t r = 0; r < rit->second.retrial_proc.size(); ++r)
                if (!rit->second.retrial_proc[r].disabled) { is_retrial = true; break; }
    }
    if (is_retrial) {
        std::size_t rr = R;  // first class actually present
        for (std::size_t r = 0; r < R; ++r)
            if (n[r] > 0) { rr = r; break; }
        if (rr == R) {  // empty station: one idle row, no orbit slots
            std::size_t w = 0;
            for (std::size_t r = 0; r < R; ++r) w += phases[r];
            return std::vector<std::vector<T>>{std::vector<T>(w, num_traits<T>::from_int(0))};
        }
        const double S = st.nservers;
        const std::size_t maxsrv =
            std::isfinite(S) ? std::min(n[rr], static_cast<std::size_t>(S)) : n[rr];
        const std::size_t maxorbit = n[rr];
        std::vector<std::vector<T>> res3;
        for (std::size_t csrv = 0; csrv <= maxsrv; ++csrv) {
            const std::size_t orbit = n[rr] - csrv;
            std::vector<T> buf(maxorbit - orbit, num_traits<T>::from_int(0));
            buf.insert(buf.end(), orbit, num_traits<T>::from_int(static_cast<long>(rr + 1)));
            std::vector<std::vector<T>> srv3;
            bool ok3 = true;
            for (std::size_t cls = 0; cls < R; ++cls) {
                const std::size_t want = cls == rr ? csrv : 0;
                const std::vector<std::vector<T>> sc = space_closed_single<T>(phases[cls], want);
                if (sc.empty() && want > 0) { ok3 = false; break; }
                srv3 = cartesian(srv3, sc);
            }
            if (!ok3) continue;
            for (std::size_t i3 = 0; i3 < srv3.size(); ++i3) {
                std::vector<T> row = buf;
                row.insert(row.end(), srv3[i3].begin(), srv3[i3].end());
                res3.push_back(row);
            }
        }
        return res3;
    }

    // The capacity gate of the reference: a Source has no finite buffer, every
    // other station refuses a marginal it cannot hold, returning NO rows rather
    // than an unreachable one.
    std::vector<std::vector<T>> out;
    if (sched != SchedStrategy::EXT)
        for (std::size_t r = 0; r < R; ++r)
            // Compare in DOUBLE: classcap is +Inf for an uncapped class, and
            // casting Inf to size_t is undefined behaviour -- it made this gate
            // reject every marginal and from_marginal returned no states at all.
            if (sn.classcap.size() >= ist && r < sn.classcap[ist - 1].size() &&
                static_cast<double>(n[r]) > sn.classcap[ist - 1][r])
                return out;

    if (sched == SchedStrategy::PAS || sched == SchedStrategy::OI) {
        // `fromMarginal.m:457-475`. THE PAS/OI LOCAL STATE IS AN ORDERED LIST,
        // not a count: the row is the class index of the job in each of the
        // `sn.cap(ist)` positions, left-aligned and zero-padded. There is no
        // buffer/server split, no phase block, and no server count, which is
        // why this arm returns before the cartesian fold below.
        //
        // Without it a PAS station fell through to the count-shaped default,
        // one column per class, while `to_marginal`'s PAS arm decoded the same
        // row as a list -- the two encodings disagreed and nothing errored. A
        // count of 2 in column 1 decoded as one class-1 job in position 0 and
        // whatever class index 2 names in position 1, so queue length at a PAS
        // station read near zero in every consumer of the decode: SolverSSA's
        // serial analyzer AND SolverCTMC's averages, which share this walk.
        //
        // REFUSE AN INFINITE CAPACITY, as the reference does. The list encoding
        // has no width without one and substituting a default would silently
        // truncate the state space instead of reporting that the model is
        // underspecified.
        const double W = ist <= sn.cap.size() ? sn.cap[ist - 1]
                                              : std::numeric_limits<double>::infinity();
        if (!std::isfinite(W))
            throw InputError(
                "from_marginal: PAS stations require finite capacity for state-space generation");
        const std::size_t w = static_cast<std::size_t>(W);
        std::size_t tot = 0;
        for (std::size_t r = 0; r < R; ++r) tot += n[r];
        if (tot == 0)
            return std::vector<std::vector<T> >(
                1, std::vector<T>(w, num_traits<T>::from_int(0)));
        if (tot > w) return out;  // infeasible: exceeds the station's total capacity
        std::vector<std::size_t> vi;
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t j = 0; j < n[r]; ++j) vi.push_back(r + 1);
        const std::vector<std::vector<std::size_t> > mi = pas_multiset_perms(vi);
        out.reserve(mi.size());
        for (std::size_t i = 0; i < mi.size(); ++i) {
            std::vector<T> row;
            row.reserve(w);
            for (std::size_t j = 0; j < mi[i].size(); ++j)
                row.push_back(num_traits<T>::from_int(static_cast<long>(mi[i][j])));
            row.resize(w, num_traits<T>::from_int(0));
            out.push_back(row);
        }
        return out;
    }

    // `space_closed_single` is empty exactly when the class has no phase, so
    // this is the unreachability test the eager fold below used to perform,
    // separated from the fold itself. THE FOLD IS DEFERRED to its only reader,
    // the INF/PS arm: it places every job of every class in a phase, which is
    // C(n+K-1, K-1) rows for K phases, and the ordered-buffer arm right below
    // never reads it -- it takes si[r] from the permutation tail. Building it
    // there anyway cost 8.5e8 discarded rows on an M/ME/1 with an order-11
    // service and cutoff 30 (test_cme_distribution.py::test_ctmc_mg1_is_exact),
    // which reads as a hang rather than as an error.
    //
    // THE TEST IS ON THE CONTENT, NOT ON THE WIDTH. `phases` is the WIDTH
    // vector (`phasessz_of`, one column per class even where the class has no
    // process), so a disabled class still occupies a column and the emptiness
    // has to be asked of `phases_of` directly.
    for (std::size_t r = 0; r < R; ++r)
        if (sn.phases_of(ist, r + 1) == 0 && n[r] > 0)
            return out;  // no service process, yet jobs demanded: unreachable

    if (state_detail::buffer_is_class_tag(sched) || state_detail::buffer_is_tag_phase_pairs(sched)) {
        // ORDERED BUFFER, as `fromMarginal`'s FCFS/LCFS branch builds it: the
        // waiting positions record WHICH class occupies each slot, so the state
        // is a permutation of the multiset {r repeated n[r] times}, not a count.
        // The last S entries of each permutation are the jobs in service and
        // the rest is the buffer, which is why the marginal alone does not
        // determine the state and a count-shaped buffer would be wrong.
        // The preempt family stores [class, phase] PAIRS, so one permutation of
        // the waiting classes yields one row per assignment of an interruption
        // phase to each waiting job, exactly as `fromMarginal` interleaves
        // `mi_buf` with `bkstate` (fromMarginal.m:310-325).
        const bool paired = state_detail::buffer_is_tag_phase_pairs(sched);
        std::vector<std::size_t> vi;
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t j = 0; j < n[r]; ++j) vi.push_back(r + 1);
        const double S = st.nservers;
        const std::size_t nsrv =
            std::isfinite(S) ? static_cast<std::size_t>(S) : vi.size();
        if (vi.empty()) {
            // Empty station: one state, all-zero buffer and idle servers. The
            // preempt-resume buffer is [class, phase] pairs, so its empty width
            // is even -- a one-column buffer there would misalign to_marginal.
            const std::size_t bw = state_detail::buffer_is_tag_phase_pairs(sched) ? 2u : 1u;
            std::vector<T> row(bw, num_traits<T>::from_int(0));
            std::size_t srvw2 = 0;
            for (std::size_t r = 0; r < R; ++r) srvw2 += phases[r];
            row.insert(row.end(), srvw2, num_traits<T>::from_int(0));
            return std::vector<std::vector<T>>{row};
        }
        // Descending first, walked down with prev_permutation: the SAME set of
        // permutations as an ascending next_permutation walk, but row 0 is the
        // descending buffer, which is the row MATLAB State.fromMarginalAndStarted
        // returns. Anything taking rows[0] as the initial state (default_init_state)
        // then agrees with the reference, and on a chain made reducible by
        // non-overtaking routing that is what picks the closed communicating class.
        std::sort(vi.begin(), vi.end(), [](std::size_t a, std::size_t b) { return a > b; });
        std::vector<std::vector<T>> res2;
        do {
            // Split the permutation: the tail is in service, the head waits.
            const std::size_t insrv = std::min(nsrv, vi.size());
            std::vector<std::size_t> si(R, 0);
            for (std::size_t j = vi.size() - insrv; j < vi.size(); ++j) si[vi[j] - 1] += 1;
            std::vector<std::vector<T>> kst;
            bool ok2 = true;
            for (std::size_t r = 0; r < R; ++r) {
                const std::vector<std::vector<T>> sr = space_closed_single<T>(phases[r], si[r]);
                if (sr.empty() && si[r] > 0) { ok2 = false; break; }
                kst = cartesian(kst, sr);
            }
            if (!ok2) continue;
            std::vector<std::size_t> wait;
            for (std::size_t j = 0; j + insrv < vi.size(); ++j) wait.push_back(vi[j]);
            std::vector<std::vector<T>> bufs;
            if (!paired) {
                std::vector<T> b;
                for (std::size_t j = 0; j < wait.size(); ++j)
                    b.push_back(num_traits<T>::from_int(static_cast<long>(wait[j])));
                if (b.empty()) b.push_back(num_traits<T>::from_int(0));
                bufs.push_back(b);
            } else if (wait.empty()) {
                // An empty paired buffer is still two columns wide, so that a
                // narrower row left-padded to the widest one keeps its parity.
                bufs.push_back(std::vector<T>(2, num_traits<T>::from_int(0)));
            } else {
                bufs.push_back(std::vector<T>());
                for (std::size_t j = 0; j < wait.size(); ++j) {
                    std::vector<std::vector<T>> next;
                    for (std::size_t b = 0; b < bufs.size(); ++b)
                        for (std::size_t p = 0; p < phases[wait[j] - 1]; ++p) {
                            std::vector<T> row = bufs[b];
                            row.push_back(num_traits<T>::from_int(static_cast<long>(wait[j])));
                            row.push_back(num_traits<T>::from_int(static_cast<long>(p + 1)));
                            next.push_back(row);
                        }
                    bufs.swap(next);
                }
            }
            for (std::size_t bi = 0; bi < bufs.size(); ++bi)
                for (std::size_t i2 = 0; i2 < kst.size(); ++i2) {
                    std::vector<T> row = bufs[bi];
                    row.insert(row.end(), kst[i2].begin(), kst[i2].end());
                    res2.push_back(row);
                }
        } while (std::prev_permutation(vi.begin(), vi.end()));
        return res2;
    }

    if (!state_detail::buffer_is_per_class_count(sched)) {
        // The INF/PS families track only jobs in the servers, so every job of
        // every class carries a phase and this fold IS the local state space.
        std::vector<std::vector<T>> srv;
        for (std::size_t r = 0; r < R; ++r)
            srv = cartesian(srv, space_closed_single<T>(phases[r], n[r]));
        return srv;
    }

    // Per-class-count buffer: everything beyond the servers waits. The split
    // between in-service and waiting is NOT determined by the marginal, so it is
    // enumerated here.
    //
    // THE SPLIT MUST DRIVE THE SERVER BLOCK, NOT BE RECOVERED FROM IT. The `srv`
    // fold above places all n[r] jobs of every class, because that is what the
    // INF and PS families need. An earlier version of this branch built the
    // server block the same way and then tried to read the in-service count back
    // out of it -- but that count is identically n[r], the very number just
    // placed there, so the buffer came out identically zero and every marginal
    // whose total exceeded the server count was discarded. A POLLING station
    // then held at most `nservers` jobs and never queued one: a two-class model
    // with lambda 0.2/0.3 and mu 2 reported the throughput of an M/M/1/1,
    // 0.16/0.24 against the correct 0.20/0.30, with QLen equal to Util in every
    // state because nothing ever waited. The buffered states were not merely
    // improbable, they were absent, so the arrival handler's successors had
    // nowhere to land and the generator dropped those edges.
    //
    // The other two branches of this function already do it this way -- the
    // class-tag branch takes si[r] from the permutation tail, the retrial branch
    // takes `want` from the in-service count -- and they are the two that worked.
    //
    // HOW MANY JOBS ARE IN SERVICE IS A PROPERTY OF THE DISCIPLINE. The
    // reference splits this family in two, and the split is exactly the
    // work-conservation of the station:
    //
    //   POLLING has its own `fromMarginal` case, and it emits the EMPTY-facility
    //   row -- all n waiting, nothing served -- alongside the one-job-in-service
    //   rows. A polling server may sit idle with a backlog because it is walking
    //   between buffers or parked, so that row is a state it genuinely occupies.
    //   Its single facility also means at most one job is ever in service.
    //
    //   SIRO, SEPT, LEPT and the SRPT pair are work conserving, and the
    //   reference enumerates EXACTLY min(S, sum n) in service (`multichoosecon`).
    //   Emitting the partially-idle rows for them would widen the space past
    //   MATLAB's with states no transition enters.
    //
    // So the bound is `sum s == maxsrv`, relaxed to `sum s >= 0` for POLLING.
    const double S = st.nservers;
    if (sched == SchedStrategy::POLLING && S != 1.0)
        throw UnsupportedError(
            "from_marginal: a polling station must have exactly one server; the controller "
            "encoding pins the visit to a single class at a time, so a multi-server polling "
            "station cannot be represented (the reference refuses it the same way)");
    const std::size_t nsrv =
        std::isfinite(S) ? static_cast<std::size_t>(S) : static_cast<std::size_t>(-1);
    std::size_t ntot = 0;
    for (std::size_t r = 0; r < R; ++r) ntot += n[r];
    const std::size_t maxsrv = std::min(nsrv, ntot);
    const std::size_t minsrv = sched == SchedStrategy::POLLING ? 0 : maxsrv;

    std::vector<std::vector<T>> res;
    std::vector<std::size_t> sv(R, 0);
    for (;;) {
        std::size_t stot = 0;
        for (std::size_t r = 0; r < R; ++r) stot += sv[r];
        if (stot >= minsrv && stot <= maxsrv) {
            std::vector<std::vector<T>> srv_s;
            bool ok = true;
            for (std::size_t r = 0; r < R && ok; ++r) {
                const std::vector<std::vector<T>> sr = space_closed_single<T>(phases[r], sv[r]);
                if (sr.empty() && sv[r] > 0) ok = false;
                else srv_s = cartesian(srv_s, sr);
            }
            if (ok) {
                std::vector<T> buf(R, num_traits<T>::from_int(0));
                for (std::size_t r = 0; r < R; ++r)
                    buf[r] = num_traits<T>::from_int(static_cast<long>(n[r] - sv[r]));
                for (std::size_t i = 0; i < srv_s.size(); ++i) {
                    std::vector<T> row = buf;
                    row.insert(row.end(), srv_s[i].begin(), srv_s[i].end());
                    res.push_back(row);
                }
            }
        }
        // Mixed-radix increment of the per-class in-service counts.
        std::size_t r = 0;
        for (; r < R; ++r) {
            if (sv[r] < n[r]) { ++sv[r]; break; }
            sv[r] = 0;
        }
        if (r == R) break;
    }
    return res;
}

/**
 * Port of `State.pollingBlocks` + `State.pollingProject`: every controller
 * configuration compatible with one (buffer, server) row.
 *
 * ENUMERATING PER ROW rather than taking a blind cartesian product is what
 * keeps the space tight and the chain irreducible. pos is PINNED to the class
 * in service; a switchover excludes a busy facility; and a park excludes a
 * non-empty station, because with only immediate switchovers the server would
 * have reached the waiting work in zero time. A cartesian product would admit
 * "serving buffer 1 while a class-2 job holds the server", which no transition
 * can enter or leave consistently with the marginals.
 *
 * @param srvclass 1-based class holding the single service facility, 0 if idle
 * @param nbuf     per-class waiting counts
 * @param pi the polling configuration of the station
 * @return one row per configuration, PROJECTED onto the columns `pinfo`
 *         materializes; the elided ones are reconstructible from the rest of
 *         the state (`polling_get`), so keeping them would split each state
 *         into copies no observation can tell apart
 */
template <class T>
std::vector<std::vector<T>> polling_blocks(const PollingInfo<T>& pi, std::size_t srvclass,
                                           const std::vector<std::size_t>& nbuf) {
    const std::size_t R = pi.polled.size();
    std::vector<std::vector<long>> trips;  // full [pos, swk, ctr]
    if (srvclass > 0) {
        // A buffer outside the cyclic order can hold no job at all.
        if (!pi.polled[srvclass - 1]) return std::vector<std::vector<T>>();
        std::vector<long> ctrset;
        switch (pi.ptype) {
            case lang::PollingType::EXHAUSTIVE:
                ctrset.push_back(0);
                break;
            case lang::PollingType::GATED:
                // ctr counts the jobs admitted at the polling instant that have
                // not completed, the one in service included, so ctr >= 1 and
                // the ctr-1 still uncompleted ones all wait in the buffer.
                for (long c = 1; c <= static_cast<long>(nbuf[srvclass - 1]) + 1; ++c)
                    ctrset.push_back(c);
                break;
            case lang::PollingType::KLIMITED:
                for (long c = 1; c <= static_cast<long>(pi.pk); ++c) ctrset.push_back(c);
                break;
            case lang::PollingType::DECREMENTING:
                // ctr is the population the visit is driving the class down to.
                for (long c = 0; c <= static_cast<long>(nbuf[srvclass - 1]); ++c)
                    ctrset.push_back(c);
                break;
        }
        for (std::size_t i = 0; i < ctrset.size(); ++i)
            trips.push_back(std::vector<long>{static_cast<long>(srvclass), 0, ctrset[i]});
    } else {
        bool anysw = false;
        for (std::size_t q = 1; q <= R; ++q) {
            if (!pi.has_sw[q - 1]) continue;
            anysw = true;
            // Every phase of a non-immediate switchover is dwelt in, and jobs
            // may wait meanwhile: this is exactly what makes a polling station
            // non-work-conserving.
            for (std::size_t k = 1; k <= pi.ksw[q - 1]; ++k)
                trips.push_back(std::vector<long>{static_cast<long>(q), static_cast<long>(k), 0});
        }
        std::size_t total = 0;
        for (std::size_t r = 0; r < nbuf.size(); ++r) total += nbuf[r];
        if (!anysw && total == 0) {
            std::size_t first = 0;
            for (std::size_t r = 1; r <= R; ++r)
                if (pi.polled[r - 1]) { first = r; break; }
            trips.push_back(std::vector<long>{static_cast<long>(first), 0, 0});
        }
    }

    const std::size_t npos = static_cast<std::size_t>(-1);
    std::vector<std::vector<T>> out;
    for (std::size_t i = 0; i < trips.size(); ++i) {
        std::vector<T> row;
        if (pi.ipos != npos) row.push_back(num_traits<T>::from_int(trips[i][0]));
        if (pi.iswk != npos) row.push_back(num_traits<T>::from_int(trips[i][1]));
        if (pi.ictr != npos) row.push_back(num_traits<T>::from_int(trips[i][2]));
        out.push_back(row);
    }
    return out;
}

/**
 * Append the trailing local-variable block to every row the core builders
 * produce, which is what makes a state row as wide as `nvars_of` declares.
 *
 * Shared by `from_marginal` and `from_marginal_and_started`: the two differ in
 * how the (buffer, server) split is chosen, never in what follows it, so the
 * controller, the BAS marker and the reply counters are appended once here.
 *
 * WITHOUT THIS THE ROWS ARE ONE BLOCK TOO NARROW and every slicer that takes
 * `nvar` clear of the right-hand end reads the server block from the wrong
 * offset. Measured on a two-class EXHAUSTIVE polling station: the controller
 * was absent from the enumerated space entirely and the station behaved as a
 * capacity-one queue, reporting Tput 0.16/0.24 against MATLAB's 0.20/0.30.
 *
 * The polling controller is enumerated per row (`polling_blocks`), and so is the
 * ROUND-ROBIN POINTER: which link the next job takes is a coordinate of the
 * state, so a space that emitted it at zero would carry one configuration of a
 * dispatcher that has several, and the chain could never leave it. RROBIN
 * enumerates its outlinks and WRROBIN the positions of its weighted cycle,
 * exactly as `State.fromMarginal`'s `sub_routevars` does.
 *
 * The modulating phase and the REPLY blocked-server counters are still emitted
 * at zero: exact for a model that declares none of them, and width-correct
 * rather than enumerated for one that does; see `14-cpp-multiprecision.md`.
 */
template <class T>
std::vector<std::vector<T>> append_local_vars(const NetworkStruct<T>& sn, std::size_t ist,
                                              std::vector<std::vector<T>> rows,
                                              const std::vector<std::size_t>& n,
                                              const std::vector<std::size_t>& phases) {
    if (rows.empty()) return rows;
    const std::size_t ind = sn.node_of_station(ist);
    // No node means no local-variable block to append; a Layer submodel carries
    // stations without a node map at all.
    if (ind == 0) return rows;
    const std::size_t width = sn.nvars_of(ind);
    if (width == 0) return rows;

    const std::size_t R = sn.nclasses;
    const PollingInfo<T> pi = polling_info(sn, ind);
    const std::size_t pw = pi.valid ? pi.width : 0;
    const T zero = num_traits<T>::from_int(0);

    // The true-BAS blocked marker takes the shared node-block column, so it is
    // MUTUALLY EXCLUSIVE with the polling controller and is enumerated instead of
    // it. Only a NON-EMPTY station can hold a blocked job -- the marker says "the
    // front job here has completed and is waiting for room" -- so an empty
    // marginal carries the single value 0 rather than both.
    const bool bas = ind <= sn.isbasblocking.size() && sn.isbasblocking[ind - 1];
    std::size_t ntot = 0;
    for (std::size_t r = 0; r < n.size(); ++r) ntot += n[r];

    // The dispatch pointers, one column per class that routes round-robin, in
    // CLASS ORDER and immediately after the modulating phases: that is where
    // `rr_var_slot` reads them and where the reference appends them.
    std::size_t rrw = 0;
    if (ind <= sn.nvars.size())
        for (std::size_t r = 0; r < R && R + r < sn.nvars[ind - 1].size(); ++r)
            rrw += sn.nvars[ind - 1][R + r];
    std::vector<std::vector<T>> ptrsets(1, std::vector<T>());
    for (std::size_t r = 1; r <= R; ++r) {
        if (sn.rr_var_slot(ind, r) == 0) continue;
        std::vector<T> vals;
        if (sn.nodes[ind - 1].routing[r - 1] == RoutingStrategy::RROBIN) {
            const std::vector<std::size_t> ol = sn.rr_outlinks(ind, r);
            for (std::size_t d = 0; d < ol.size(); ++d)
                vals.push_back(num_traits<T>::from_int(static_cast<long>(ol[d])));
        } else {
            const std::vector<std::size_t> cy = sn.rr_weighted_outlinks(ind, r);
            for (std::size_t d = 0; d < cy.size(); ++d)
                vals.push_back(num_traits<T>::from_int(static_cast<long>(d + 1)));
        }
        if (vals.empty()) vals.push_back(zero);
        std::vector<std::vector<T>> grown;
        for (std::size_t g = 0; g < ptrsets.size(); ++g)
            for (std::size_t v = 0; v < vals.size(); ++v) {
                std::vector<T> row = ptrsets[g];
                row.push_back(vals[v]);
                grown.push_back(row);
            }
        ptrsets.swap(grown);
    }

    std::vector<std::vector<T>> out;
    for (std::size_t i = 0; i < rows.size(); ++i) {
        // The controller depends on WHICH class holds the service facility and
        // on how many jobs wait, so recover both from the row just built.
        std::vector<std::size_t> nbuf(R, 0);
        std::size_t srvclass = 0;
        if (pw) {
            std::size_t srvw = 0;
            for (std::size_t r = 0; r < R; ++r) srvw += phases[r];
            const std::size_t srv0 = rows[i].size() - srvw;
            std::size_t off = 0;
            for (std::size_t r = 0; r < R; ++r) {
                std::size_t c = 0;
                for (std::size_t k = 0; k < phases[r]; ++k)
                    c += static_cast<std::size_t>(
                        num_traits<T>::to_double(rows[i][srv0 + off + k]));
                off += phases[r];
                if (c > 0 && srvclass == 0) srvclass = r + 1;
                // The per-class-count buffer POLLING uses puts the waiting
                // count in column r, ahead of the server block.
                if (srv0 >= R) nbuf[r] = static_cast<std::size_t>(
                    num_traits<T>::to_double(rows[i][r]));
            }
        }
        const std::vector<std::vector<T>> blocks =
            pw ? polling_blocks(pi, srvclass, nbuf)
               : std::vector<std::vector<T>>(1, std::vector<T>());
        // No admissible controller configuration means the (buffer, server)
        // split itself is unoccupiable, so the row is DROPPED rather than
        // emitted without a controller.
        // ORDER IS THE ENCODING: nvars is [modulation | routing | node block |
        // reply], so the controller sits between the routing pointers and the
        // reply counters, not at the end. `polling_info` computes its offset
        // the same way, and the two must agree or every read is shifted.
        const std::size_t head = pi.valid ? pi.off : width - pw;
        for (std::size_t b = 0; b < blocks.size(); ++b) {
        for (std::size_t g = 0; g < ptrsets.size(); ++g) {
            std::vector<T> row = rows[i];
            row.insert(row.end(), head - rrw, zero);  // modulating phases
            row.insert(row.end(), ptrsets[g].begin(), ptrsets[g].end());  // dispatch pointers
            row.insert(row.end(), blocks[b].begin(), blocks[b].end());
            row.insert(row.end(), width - head - pw, zero);  // reply counters
            if (!bas) {
                out.push_back(row);
                continue;
            }
            // The marker is the LAST column, which is where the departure
            // handler and the generator's become-blocked edge both read it.
            // `refresh_bas_blocking` refuses BAS together with a reply block for
            // exactly this reason: the reply counters would trail it.
            for (std::size_t v = 0; v <= (ntot > 0 ? 1u : 0u); ++v) {
                std::vector<T> r2 = row;
                r2.back() = num_traits<T>::from_int(static_cast<long>(v));
                out.push_back(r2);
            }
        }
        }
    }
    return out;
}

template <class T>
std::vector<std::vector<T>> from_marginal(const NetworkStruct<T>& sn, std::size_t ist,
                                          const std::vector<std::size_t>& n,
                                          const std::vector<std::size_t>& phases) {
    const std::size_t jnd0 = sn.node_of_station(ist);
    // A Join of an FJ-augmented struct: its state is the per-class count vector
    // itself, DETERMINED by the marginal rather than enumerated from it. There is
    // no buffer order to choose and no phase to be in, so exactly one row exists.
    if (sn.isfjaugmented && jnd0 != 0 && sn.nodes[jnd0 - 1].nodetype == NodeType::Join) {
        std::vector<T> row(sn.nclasses, num_traits<T>::from_int(0));
        for (std::size_t r = 0; r < sn.nclasses && r < n.size(); ++r)
            row[r] = num_traits<T>::from_int(static_cast<long>(n[r]));
        return std::vector<std::vector<T>>(1, row);
    }
    // SYNCHRONOUS CALL (REPLY signal), `State.fromMarginal:44-90`. A node with a
    // reply block holds one server per job that has left for its callee and is
    // waiting for the reply, and THAT COUNT IS NOT DERIVABLE FROM THE MARGINAL:
    // the job is at the callee, not here. So the held counts are ENUMERATED and
    // the rest of the state is built with the REMAINING servers -- with b held,
    // only S-b jobs can be in service, a configuration the plain enumeration
    // never produces. Without this the counter column exists but is zero in
    // every row, the departure that would hold a server has no successor in the
    // space, and the chain dead-ends: measured on the `test_ctmc_reply` model,
    // 11 states with the whole population frozen at the caller and every rate
    // zero, against 21 enumerated and X = 0.47059 in MATLAB, the JAR and Python.
    if (jnd0 != 0 && sn.replyblock.size() >= jnd0) {
        std::vector<std::size_t> rclasses;
        for (std::size_t r = 0; r < sn.replyblock[jnd0 - 1].size(); ++r)
            if (sn.replyblock[jnd0 - 1][r]) rclasses.push_back(r + 1);
        if (!rclasses.empty()) {
            const double Sd = sn.stations[ist - 1].nservers;
            const std::size_t S =
                std::isfinite(Sd) ? static_cast<std::size_t>(Sd) : static_cast<std::size_t>(0);
            // The recursion runs on a struct with the block CLEARED, which is
            // what terminates it: the copy takes this branch no further.
            NetworkStruct<T> snb = sn;
            snb.replyblock[jnd0 - 1].assign(snb.replyblock[jnd0 - 1].size(), false);
            if (snb.nvars.size() >= jnd0)
                for (std::size_t r = 1; r <= sn.nclasses && 2 * sn.nclasses + r < snb.nvars[jnd0 - 1].size();
                     ++r)
                    snb.nvars[jnd0 - 1][2 * sn.nclasses + r] = 0;
            std::vector<std::vector<std::size_t>> bspace(1, std::vector<std::size_t>());
            for (std::size_t i = 0; i < rclasses.size(); ++i) {
                std::vector<std::vector<std::size_t>> next;
                for (std::size_t j = 0; j < bspace.size(); ++j)
                    for (std::size_t v = 0; v <= S; ++v) {
                        std::vector<std::size_t> row = bspace[j];
                        row.push_back(v);
                        next.push_back(row);
                    }
                bspace.swap(next);
            }
            std::vector<std::vector<std::vector<T>>> subs(bspace.size());
            std::size_t maxw = 0;
            for (std::size_t bi = 0; bi < bspace.size(); ++bi) {
                std::size_t tot = 0;
                for (std::size_t i = 0; i < bspace[bi].size(); ++i) tot += bspace[bi][i];
                if (tot > S) continue;
                snb.stations[ist - 1].nservers = static_cast<double>(S - tot);
                subs[bi] = from_marginal(snb, ist, n, phases);
                for (std::size_t i = 0; i < subs[bi].size(); ++i)
                    maxw = std::max(maxw, subs[bi][i].size());
            }
            // A held server pushes a job into the buffer, so the sub-spaces come
            // out at DIFFERENT buffer widths. The buffer is RIGHT-aligned, empty
            // slots padding the left, so widen the narrow rows on the left.
            std::vector<std::vector<T>> out;
            for (std::size_t bi = 0; bi < bspace.size(); ++bi) {
                for (std::size_t i = 0; i < subs[bi].size(); ++i) {
                    std::vector<T> row;
                    row.reserve(maxw + bspace[bi].size());
                    if (subs[bi][i].size() < maxw)
                        row.assign(maxw - subs[bi][i].size(), num_traits<T>::from_int(0));
                    row.insert(row.end(), subs[bi][i].begin(), subs[bi][i].end());
                    for (std::size_t j = 0; j < bspace[bi].size(); ++j)
                        row.push_back(num_traits<T>::from_int(static_cast<long>(bspace[bi][j])));
                    out.push_back(row);
                }
            }
            std::sort(out.begin(), out.end(), [](const std::vector<T>& a, const std::vector<T>& b) {
                for (std::size_t i = 0; i < a.size() && i < b.size(); ++i) {
                    const double av = num_traits<T>::to_double(a[i]);
                    const double bv = num_traits<T>::to_double(b[i]);
                    if (av != bv) return av < bv;
                }
                return a.size() < b.size();
            });
            out.erase(std::unique(out.begin(), out.end(),
                                  [](const std::vector<T>& a, const std::vector<T>& b) {
                                      if (a.size() != b.size()) return false;
                                      for (std::size_t i = 0; i < a.size(); ++i)
                                          if (num_traits<T>::to_double(a[i]) !=
                                              num_traits<T>::to_double(b[i]))
                                              return false;
                                      return true;
                                  }),
                      out.end());
            return out;
        }
    }
    return append_local_vars(sn, ist, from_marginal_core(sn, ist, n, phases), n, phases);
}

/**
 * Port of `State.fromMarginalAndStarted`: ONE state realizing both a per-class
 * occupancy `n` and a per-class STARTED count `s`.
 *
 * HOW IT DIFFERS FROM `from_marginal_core`, which is the reason it is a separate
 * builder rather than a filter over it. `from_marginal` ENUMERATES every
 * (buffer, server) split consistent with a marginal, because the marginal alone
 * does not determine which jobs hold the servers. Here the split is GIVEN: `s`
 * says how many jobs of each class are in service, so exactly one split is
 * meant and the buffer contents follow as `n - s`. Filtering the enumeration
 * would be both quadratic and wrong at an empty station, where the reference
 * emits a row of a prescribed width rather than selecting one.
 *
 * EVERY STARTED JOB IS PLACED IN PHASE ONE, which is what makes the result an
 * INITIAL state rather than a member of the stationary space: a job that has
 * just started service has not advanced through its phase-type law yet.
 *
 * ONE BUFFER ORDERING IS EMITTED under the ordered-buffer disciplines, the
 * descending-sorted one. The reference enumerates every permutation and then
 * keeps the lexicographic maximum through its trailing unique/flip, so building
 * that row directly is the same answer without the factorial.
 *
 * A SHARED SERVER IGNORES `s` ENTIRELY. Under INF/PS/DPS/GPS/LPS and their
 * priority variants every job present is in a server, so the state holds `n`,
 * not `s`; writing `s` there would lose the queued jobs. That is not a
 * simplification but the reference's own branch, and it was a real defect in
 * the Python twin until 2026-08-09.
 *
 * PAS/OI ignores `s` for a different reason: its local state is the ordered
 * class-index list of the jobs present, which carries no service split at all,
 * so this delegates to the marginal builder there.
 *
 * @param sn     the network struct
 * @param ist    station index (1-based)
 * @param n      per-class occupancy
 * @param s      per-class started count; must satisfy s[r] <= n[r]
 * @param phases per-class phase counts
 */
template <class T>
std::vector<std::vector<T>> from_marginal_and_started_core(const NetworkStruct<T>& sn,
                                                           std::size_t ist,
                                                           const std::vector<std::size_t>& n,
                                                           const std::vector<std::size_t>& s,
                                                           const std::vector<std::size_t>& phases) {
    const std::size_t R = sn.nclasses;
    if (ist == 0 || ist > sn.stations.size())
        throw InputError("from_marginal_and_started: station index " + std::to_string(ist) +
                         " is out of range");
    if (n.size() != R || s.size() != R || phases.size() != R)
        throw InputError("from_marginal_and_started: n, s and phases must have one entry per class");

    const Station<T>& st = sn.stations[ist - 1];
    const SchedStrategy sched = st.sched;
    std::vector<std::vector<T>> out;

    std::size_t ntot = 0, stot = 0;
    for (std::size_t r = 0; r < R; ++r) {
        if (s[r] > n[r]) return out;  // more started than present: no such state
        ntot += n[r];
        stot += s[r];
    }

    // The reference's two pre-switch guards. Both are "no such state" rather
    // than an error, since the caller may be sweeping a lattice.
    if (ist <= sn.classcap.size())
        for (std::size_t r = 0; r < R && r < sn.classcap[ist - 1].size(); ++r)
            if (static_cast<double>(n[r]) > sn.classcap[ist - 1][r]) return out;
    const double S = st.nservers;
    if (S > 0.0 && std::isfinite(S) && static_cast<double>(stot) > S) return out;

    // PAS/OI: the started counts are immaterial to the list encoding.
    if (sched == SchedStrategy::PAS || sched == SchedStrategy::OI)
        return from_marginal_core(sn, ist, n, phases);

    // A Source generates rather than holds: its phase block IS the arrival
    // process, seeded in phase one for every class it serves, and the leading
    // column is the infinite population the reference marks with Inf.
    if (sched == SchedStrategy::EXT || st.nodetype == NodeType::Source) {
        // ONE BLOCK PER CLASS, INCLUDING A CLASS THE SOURCE DOES NOT GENERATE:
        // `phases` is the width vector, so a disabled arrival contributes one
        // always-zero column rather than nothing. That is the reference's row
        // (`[Inf 1 0 0 0 0 0]` for a Source generating the first of six
        // classes), and a narrower one cannot be decoded by anyone who built it
        // from `sn.phasessz` -- see `NetworkStruct::phasessz_of`.
        std::vector<std::vector<T>> srv;
        for (std::size_t r = 0; r < R; ++r) {
            if (phases[r] == 0) continue;
            std::vector<T> init(phases[r], num_traits<T>::from_int(0));
            if (r < sn.classes.size() && std::isinf(sn.classes[r].population) &&
                !sn.disabled[ist - 1][r])
                init[0] = num_traits<T>::from_int(1);
            srv = cartesian(srv, std::vector<std::vector<T>>(1, init));
        }
        if (srv.empty()) srv.push_back(std::vector<T>());
        for (std::size_t i = 0; i < srv.size(); ++i) {
            std::vector<T> row(1, num_traits<T>::from_double(
                                     std::numeric_limits<double>::infinity()));
            row.insert(row.end(), srv[i].begin(), srv[i].end());
            out.push_back(row);
        }
        return out;
    }

    /** The server block: `cnt[r]` jobs of class r, all of them in phase one. */
    const auto phase_one_block = [&](const std::vector<std::size_t>& cnt,
                                     bool* ok) -> std::vector<T> {
        std::vector<T> blk;
        *ok = true;
        for (std::size_t r = 0; r < R; ++r) {
            // No service process, yet jobs demanded there: unreachable. An
            // empty factor must annihilate the fold, not be absorbed by it.
            // The COLUMN is still emitted -- `phases` is the width vector and a
            // disabled class holds an always-zero column, as the reference does.
            if (sn.phases_of(ist, r + 1) == 0 && cnt[r] > 0) *ok = false;
            if (phases[r] == 0) continue;
            blk.push_back(num_traits<T>::from_int(static_cast<long>(cnt[r])));
            blk.insert(blk.end(), phases[r] - 1, num_traits<T>::from_int(0));
        }
        return blk;
    };
    std::size_t srvw = 0;
    for (std::size_t r = 0; r < R; ++r) srvw += phases[r];

    if (state_detail::buffer_is_per_class_count(sched)) {
        // UNORDERED buffer: one waiting COUNT per class ahead of the server
        // block. Below the server count every job is in service and the buffer
        // is empty, which is the reference's own special case.
        const bool all_in_service = std::isfinite(S) && static_cast<double>(ntot) <= S;
        std::vector<std::size_t> insrv(R, 0), wait(R, 0);
        for (std::size_t r = 0; r < R; ++r) {
            insrv[r] = all_in_service ? n[r] : s[r];
            wait[r] = n[r] - insrv[r];
        }
        bool ok = false;
        const std::vector<T> blk = phase_one_block(insrv, &ok);
        if (!ok) return out;
        std::vector<T> row;
        for (std::size_t r = 0; r < R; ++r)
            row.push_back(num_traits<T>::from_int(static_cast<long>(wait[r])));
        row.insert(row.end(), blk.begin(), blk.end());
        out.push_back(row);
        return out;
    }

    const bool paired = state_detail::buffer_is_tag_phase_pairs(sched);
    if (state_detail::buffer_is_class_tag(sched) || paired) {
        // ORDERED buffer: the waiting positions carry class tags, and under the
        // preempt family a [class, phase] PAIR per position. The empty buffer
        // keeps the width the decoder expects -- one column, or two when
        // paired, since an odd-width paired buffer is half a pair and
        // `to_marginal` would read every slot shifted.
        const std::size_t bw = paired ? 2u : 1u;
        if (ntot == 0) {
            std::vector<T> row(bw + srvw, num_traits<T>::from_int(0));
            out.push_back(row);
            return out;
        }
        std::vector<std::size_t> inbuf;
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t j = 0; j < n[r] - s[r]; ++j) inbuf.push_back(r + 1);
        // Descending, which is the lexicographic maximum the reference's
        // trailing unique/flip leaves standing.
        std::sort(inbuf.begin(), inbuf.end(),
                  [](std::size_t a, std::size_t b) { return a > b; });

        bool ok = false;
        const std::vector<T> blk = phase_one_block(s, &ok);
        if (!ok) return out;

        std::vector<T> row;
        if (inbuf.empty()) {
            row.assign(bw, num_traits<T>::from_int(0));
        } else {
            for (std::size_t j = 0; j < inbuf.size(); ++j) {
                row.push_back(num_traits<T>::from_int(static_cast<long>(inbuf[j])));
                // A preempted job is recorded in its LAST service phase, as the
                // reference builds the lexicographic maximum.
                if (paired)
                    row.push_back(num_traits<T>::from_int(
                        static_cast<long>(phases[inbuf[j] - 1])));
            }
        }
        row.insert(row.end(), blk.begin(), blk.end());
        out.push_back(row);
        return out;
    }

    // SHARED SERVER (INF/PS/DPS/GPS/LPS and the priority variants): no buffer,
    // and every job present is in a server, so the block holds n and not s.
    bool ok = false;
    const std::vector<T> blk = phase_one_block(n, &ok);
    if (!ok) return out;
    out.push_back(blk);
    return out;
}

/**
 * `State.fromMarginalAndStarted` with the trailing local-variable block
 * appended, i.e. the counterpart of `from_marginal` for a prescribed service
 * split. A Join of an FJ-augmented struct carries its per-class counts and no
 * split, exactly as in `from_marginal`.
 */
template <class T>
std::vector<std::vector<T>> from_marginal_and_started(const NetworkStruct<T>& sn, std::size_t ist,
                                                      const std::vector<std::size_t>& n,
                                                      const std::vector<std::size_t>& s,
                                                      const std::vector<std::size_t>& phases) {
    const std::size_t jnd0 = sn.node_of_station(ist);
    if (sn.isfjaugmented && jnd0 != 0 && sn.nodes[jnd0 - 1].nodetype == NodeType::Join) {
        std::vector<T> row(sn.nclasses, num_traits<T>::from_int(0));
        for (std::size_t r = 0; r < sn.nclasses && r < n.size(); ++r)
            row[r] = num_traits<T>::from_int(static_cast<long>(n[r]));
        return std::vector<std::vector<T>>(1, row);
    }
    return append_local_vars(sn, ist, from_marginal_and_started_core(sn, ist, n, s, phases), n,
                             phases);
}


/**
 * Port of `State.fromMarginal` at its OWN signature: the reference indexes by
 * NODE, not by station, and derives the station internally. That distinction is
 * load-bearing, not cosmetic -- a Transition is stateful but is NOT a station,
 * so it carries no station index and the station-indexed overload can never
 * reach it. Routing every caller through here is what makes the Transition
 * branch below live rather than dead code.
 *
 * A node that is not stateful holds no jobs and contributes no local state.
 */
template <class T>
std::vector<std::vector<T>> from_marginal_node(const NetworkStruct<T>& sn, std::size_t ind,
                                               const std::vector<std::size_t>& n,
                                               const std::vector<std::size_t>& phases) {
    if (ind == 0 || ind > sn.nodes.size())
        throw InputError("from_marginal_node: node index " + std::to_string(ind) +
                         " is out of range");
    const NodeDef& nd = sn.nodes[ind - 1];
    if (nd.station != 0) return from_marginal(sn, nd.station, n, phases);
    if (!nd.stateful) return std::vector<std::vector<T>>();

    const std::size_t R = sn.nclasses;
    if (n.size() != R)
        throw InputError("from_marginal_node: n must have one entry per class");

    if (nd.nodetype == NodeType::Transition) {
        // PER-MODE state, not per-class: [free servers | firing phases | fired].
        // The reference emits the all-idle row -- every mode's servers free, no
        // firing in progress -- and lets the event handlers walk from there, so
        // the marginal n plays no part. An infinite server count is clamped to
        // MaxInt because the row is a COUNT vector and Inf is not a count.
        const typename std::map<std::size_t, TransitionParam<T> >::const_iterator it =
            sn.transparam.find(ind);
        if (it == sn.transparam.end())
            throw InputError(
                "from_marginal_node: transition node has no TransitionParam; build it with "
                "add_transition");
        const TransitionParam<T>& tp = it->second;
        std::vector<T> row;
        for (std::size_t mm = 0; mm < tp.nmodes; ++mm) {
            const double sv = mm < tp.nmodeservers.size() ? tp.nmodeservers[mm] : 1.0;
            row.push_back(num_traits<T>::from_double(
                std::isfinite(sv) ? sv : static_cast<double>(GlobalConstants::MaxInt)));
        }
        std::size_t fph = 0;
        for (std::size_t mm = 0; mm < tp.nmodes; ++mm)
            fph += mm < tp.firingphases.size() && tp.firingphases[mm] > 0 ? tp.firingphases[mm] : 1;
        row.insert(row.end(), fph, num_traits<T>::from_int(0));
        row.insert(row.end(), tp.nmodes, num_traits<T>::from_int(0));
        return std::vector<std::vector<T> >{row};
    }

    // A CACHE HOLDS ONE READING JOB, which is what makes the read a decision
    // rather than a queue: `spaceGeneratorNodes` caps every class there at 1 and
    // the node at 1 job in total. Without the cap the lattice admits a cache
    // holding two reads, and since `after_event_cache` fires a READ only on a
    // node holding exactly one job, that state is ABSORBING -- on an open model
    // the whole stationary mass ends up in it and every cache rate reads zero.
    //
    // A completing fetch is the one exception: it releases every merged
    // secondary request into the hit class in ONE immediate transition, so the
    // hit classes momentarily hold up to `max_pending_retrieval` extra jobs.
    if (nd.nodetype == NodeType::Cache) {
        const typename std::map<std::size_t, CacheParam<T> >::const_iterator cc =
            sn.nodeparam.find(ind);
        if (cc != sn.nodeparam.end()) {
            const std::size_t pend =
                cc->second.retrieval_capacity > 0 && cc->second.max_pending_retrieval > 0
                    ? static_cast<std::size_t>(cc->second.max_pending_retrieval)
                    : 0;
            std::size_t tot = 0;
            for (std::size_t r = 0; r < R; ++r) {
                bool is_hit = false;
                for (std::size_t u = 0; u < cc->second.hitclass.size(); ++u)
                    if (cc->second.hitclass[u] == r + 1) is_hit = true;
                if (n[r] > (is_hit ? 1 + pend : 1)) return std::vector<std::vector<T> >();
                tot += n[r];
            }
            if (tot > 1 + pend) return std::vector<std::vector<T> >();
        }
    }

    // Any other stateful non-station (a Cache, a stateful ClassSwitch): the
    // reference folds `spaceClosedSingle(1, n(r))` per class, i.e. one column
    // per class holding its count, with no phase and no buffer.
    std::vector<std::vector<T> > acc;
    for (std::size_t r = 0; r < R; ++r) {
        const std::vector<std::vector<T> > sr = space_closed_single<T>(1, n[r]);
        if (sr.empty() && n[r] > 0) return std::vector<std::vector<T> >();
        acc = cartesian(acc, sr);
    }

    // A DISPATCHER CARRIES ITS POINTER IN THE STATE, for the same reason a cache
    // carries its contents: which link the next job takes is not derivable from
    // the marginal. RROBIN enumerates its outlinks and WRROBIN the positions of
    // its weighted cycle, exactly as `sub_routevars` does, and the columns sit
    // in CLASS ORDER right after the per-class counts, where `rr_var_slot` reads
    // them. Emitting the node one column narrower than `nvars_of` declares makes
    // every successor a state the space does not contain, so the node becomes
    // ABSORBING and the whole chain stalls with the jobs sitting in it.
    {
        std::vector<std::vector<T> > ptrsets(1, std::vector<T>());
        bool any = false;
        for (std::size_t r = 1; r <= R; ++r) {
            if (sn.rr_var_slot(ind, r) == 0) continue;
            any = true;
            std::vector<T> vals;
            if (sn.nodes[ind - 1].routing[r - 1] == RoutingStrategy::RROBIN) {
                const std::vector<std::size_t> ol = sn.rr_outlinks(ind, r);
                for (std::size_t u = 0; u < ol.size(); ++u)
                    vals.push_back(num_traits<T>::from_int(static_cast<long>(ol[u])));
            } else {
                const std::vector<std::size_t> cy = sn.rr_weighted_outlinks(ind, r);
                for (std::size_t u = 0; u < cy.size(); ++u)
                    vals.push_back(num_traits<T>::from_int(static_cast<long>(u + 1)));
            }
            if (vals.empty()) vals.push_back(num_traits<T>::from_int(0));
            std::vector<std::vector<T> > grown;
            for (std::size_t g = 0; g < ptrsets.size(); ++g)
                for (std::size_t v = 0; v < vals.size(); ++v) {
                    std::vector<T> row = ptrsets[g];
                    row.push_back(vals[v]);
                    grown.push_back(row);
                }
            ptrsets.swap(grown);
        }
        if (any) acc = cartesian(acc, ptrsets);
    }

    // A CACHE CARRIES ITS CONTENTS IN THE STATE, and they are not derivable
    // from any marginal: the row is [per-class counts | contents | occupancy],
    // exactly as `after_event_cache` slices it. Enumerating only the counts
    // left the contents region absent, so every READ landed on a state one
    // block wider than anything enumerated and the space came out EMPTY.
    //
    // A contents slot holds an ITEM INDEX and an item sits in at most one slot --
    // a cache holding two copies of one item is not a state the replacement
    // policies can reach or leave.
    //
    // THE ENUMERATED CACHE IS FULL, which is what `State.spaceCache` builds:
    // every ordered placement of `tcc` DISTINCT items, with no empty slot. Under
    // every replacement policy a cache that has once been filled never empties
    // again, so a partly empty configuration is transient and enumerating it
    // only spreads stationary mass over states the model leaves for good. It is
    // load-bearing beyond tidiness for a retrieval cache: block A may be set
    // only for an item that is NOT cached, and the count of such items is the
    // retrieval-system capacity `nitems - tcc`, which is a number only if the
    // cached set is full. A cache with more slots than items can never fill, so
    // the empty slot survives there.
    if (nd.nodetype == NodeType::Cache) {
        const typename std::map<std::size_t, CacheParam<T> >::const_iterator ci =
            sn.nodeparam.find(ind);
        if (ci == sn.nodeparam.end()) return acc;
        const CacheParam<T>& cp = ci->second;
        std::size_t tcc = 0;
        for (std::size_t u = 0; u < cp.itemcap.size(); ++u)
            if (cp.itemcap[u] > 0) tcc += static_cast<std::size_t>(cp.itemcap[u]);

        const std::size_t first_item = tcc <= cp.nitems ? 1 : 0;
        std::vector<std::vector<T> > contents(1, std::vector<T>());
        for (std::size_t slot = 0; slot < tcc; ++slot) {
            std::vector<std::vector<T> > grown;
            for (std::size_t c = 0; c < contents.size(); ++c)
                for (std::size_t item = first_item; item <= cp.nitems; ++item) {
                    bool dup = false;
                    if (item > 0)
                        for (std::size_t j = 0; j < contents[c].size(); ++j)
                            if (num_traits<T>::to_double(contents[c][j]) ==
                                static_cast<double>(item)) { dup = true; break; }
                    if (dup) continue;
                    std::vector<T> row = contents[c];
                    row.push_back(num_traits<T>::from_int(static_cast<long>(item)));
                    grown.push_back(row);
                }
            contents.swap(grown);
        }
        // Block A, the delayed-hit occupancy bitmap: one bit per item, set while
        // a retrieval for it is in flight.
        if (cp.retrieval_capacity > 0)
            for (std::size_t i = 0; i < cp.nitems; ++i) {
                std::vector<std::vector<T> > grown;
                for (std::size_t c = 0; c < contents.size(); ++c)
                    for (int b = 0; b <= 1; ++b) {
                        std::vector<T> row = contents[c];
                        row.push_back(num_traits<T>::from_int(b));
                        grown.push_back(row);
                    }
                contents.swap(grown);
            }
        // Block B, the merged secondary requests: one count per retrieval class,
        // nonzero only where that class's item is in flight (block A is set) and
        // summing to at most the declared truncation level. `State.spaceCache`
        // enumerates the same compositions; -1 means the caller is walking a
        // sample path rather than enumerating, so no pending state is generated.
        if (cp.retrieval_capacity > 0) {
            std::vector<std::size_t> rcl, rci, rco;
            cache_retrieval_class_map(cp, rcl, rci, rco);
            if (!rcl.empty()) {
                const long maxpend =
                    cp.max_pending_retrieval > 0 ? cp.max_pending_retrieval : 0;
                std::vector<std::vector<T> > grown;
                for (std::size_t c = 0; c < contents.size(); ++c) {
                    // which slots may carry a count in THIS row: the ones whose
                    // item is currently being fetched
                    std::vector<std::size_t> active;
                    for (std::size_t j = 0; j < rcl.size(); ++j)
                        if (num_traits<T>::to_double(contents[c][tcc + rci[j] - 1]) != 0)
                            active.push_back(j);
                    std::vector<std::vector<long> > pend(1, std::vector<long>(rcl.size(), 0));
                    for (std::size_t a = 0; a < active.size(); ++a) {
                        std::vector<std::vector<long> > next;
                        for (std::size_t q = 0; q < pend.size(); ++q) {
                            long used = 0;
                            for (std::size_t j = 0; j < rcl.size(); ++j) used += pend[q][j];
                            for (long v = 0; v + used <= maxpend; ++v) {
                                std::vector<long> row = pend[q];
                                row[active[a]] = v;
                                next.push_back(row);
                            }
                        }
                        pend.swap(next);
                    }
                    for (std::size_t q = 0; q < pend.size(); ++q) {
                        std::vector<T> row = contents[c];
                        for (std::size_t j = 0; j < rcl.size(); ++j)
                            row.push_back(num_traits<T>::from_int(pend[q][j]));
                        grown.push_back(row);
                    }
                }
                contents.swap(grown);
            }
        }
        std::vector<std::vector<T> > joint = cartesian(acc, contents);
        if (cp.retrieval_capacity <= 0) return joint;

        // THE LOCALLY INVALID ROWS OF A RETRIEVAL CACHE, `spaceGeneratorNodes.m`
        // lines 206-280. The cartesian product above enumerates combinations the
        // dynamics can never reach, and they are not harmless: a row holding a
        // returning retrieval job BESIDE another job disables the READ (which
        // needs the node to hold exactly one) while leaving the departure
        // enabled, so the fetch bounces back to its queue without completing.
        // That inflated the retrieval flow by 20% on retrieval_simple while the
        // hit and miss shares, being ratios, stayed exactly right -- the kind of
        // defect only a flow identity catches.
        std::vector<bool> is_hit(R, false), is_miss(R, false);
        for (std::size_t u = 0; u < cp.hitclass.size(); ++u)
            if (cp.hitclass[u] >= 1 && cp.hitclass[u] <= R) is_hit[cp.hitclass[u] - 1] = true;
        for (std::size_t u = 0; u < cp.missclass.size(); ++u)
            if (cp.missclass[u] >= 1 && cp.missclass[u] <= R) is_miss[cp.missclass[u] - 1] = true;
        const long maxpend = cp.max_pending_retrieval > 0 ? cp.max_pending_retrieval : 0;
        std::vector<std::vector<T> > kept;
        for (std::size_t row = 0; row < joint.size(); ++row) {
            const std::vector<T>& st = joint[row];
            double tot = 0, other = 0, miss = 0;
            for (std::size_t r = 0; r < R; ++r) {
                const double v = num_traits<T>::to_double(st[r]);
                tot += v;
                if (is_miss[r]) miss += v;
                else if (!is_hit[r]) other += v;
            }
            // ONE JOB READS AT A TIME. The single exception is the state a
            // completing fetch lands in: the miss-class job that was the fetch,
            // together with the delayed hits it released in the same transition.
            if (!(tot <= 1 || (other == 0 && miss <= 1 && tot <= 1 + maxpend))) continue;
            bool ok = true;
            long nbits = 0;
            for (std::size_t i = 0; i < cp.nitems && ok; ++i) {
                if (num_traits<T>::to_double(st[R + tcc + i]) == 0) continue;
                ++nbits;
                // An item cannot be cached and in flight at the same time.
                for (std::size_t c = 0; c < tcc; ++c)
                    if (static_cast<std::size_t>(num_traits<T>::to_double(st[R + c])) == i + 1)
                        ok = false;
            }
            if (!ok || nbits > cp.retrieval_capacity) continue;
            kept.push_back(st);
        }
        return kept;
    }
    return acc;
}

/**
 * Port of `State.fromMarginalAndStarted` at its OWN signature, which indexes by
 * NODE rather than by station, mirroring `from_marginal_node`.
 *
 * A Petri-net element is REFUSED BY NAME, as the reference refuses it: a
 * Transition's local state is per MODE, not per class, and a Place holds
 * tokens, so neither has a notion of a started job to prescribe. A node that is
 * not stateful holds no jobs and contributes no local state.
 */
template <class T>
std::vector<std::vector<T>> from_marginal_node_and_started(
    const NetworkStruct<T>& sn, std::size_t ind, const std::vector<std::size_t>& n,
    const std::vector<std::size_t>& s, const std::vector<std::size_t>& phases) {
    if (ind == 0 || ind > sn.nodes.size())
        throw InputError("from_marginal_node_and_started: node index " + std::to_string(ind) +
                         " is out of range");
    const NodeDef& nd = sn.nodes[ind - 1];
    if (nd.nodetype == NodeType::Transition || nd.nodetype == NodeType::Place)
        throw UnsupportedError(
            "from_marginal_node_and_started cannot be used on Petri net elements");
    if (nd.station != 0) return from_marginal_and_started(sn, nd.station, n, s, phases);
    if (!nd.stateful) return std::vector<std::vector<T>>();
    // Every other stateful non-station carries no service split, so the
    // marginal builder answers it: a Cache holds its reading job and its
    // contents, a stateful ClassSwitch one count per class.
    return from_marginal_node(sn, ind, n, phases);
}

/**
 * Port of `State.fromMarg`: the state space with a given TOTAL queue length.
 *
 * This is the class-summed counterpart of `from_marginal_node`. Where that one
 * fixes how many jobs of EACH class the node holds, this one fixes only how
 * many it holds ALTOGETHER and returns the union over every class split of
 * `ntot` the node can hold.
 *
 * A class DISABLED at the station has `classcap` 0 and is excluded from the
 * split enumeration up front rather than after the fact. Asking
 * `from_marginal_node` for a job of such a class yields an EMPTY local space,
 * and `space_closed_single(0,1)` returning no rows empties the whole cartesian
 * product, so the job would silently disappear -- the trap documented at
 * `space_closed_single` above.
 *
 * `ntot == 0` has the single empty split, which `from_marginal_node` answers
 * with the per-discipline empty state; the width is NOT re-derived here.
 *
 * The buffer is RIGHT-aligned, so sub-spaces of different width are padded on
 * the LEFT before they are stacked, exactly as the reference does for the
 * reply-block sub-spaces. Rows are then uniqued and reversed, which puts the
 * empty state first and the states with jobs in phase 1 earlier.
 *
 * The twin over BOTH totals is `from_marg_node_started` below.
 */
template <class T>
std::vector<std::vector<T>> from_marg_node(const NetworkStruct<T>& sn, std::size_t ind,
                                           std::size_t ntot,
                                           const std::vector<std::size_t>& phases) {
    if (ind == 0 || ind > sn.nodes.size())
        throw InputError("from_marg_node: node index " + std::to_string(ind) + " is out of range");
    const NodeDef& nd = sn.nodes[ind - 1];
    const std::size_t R = sn.nclasses;

    // Per-class capacity of the node, unbounded when it is not a station.
    std::vector<double> ccap(R, std::numeric_limits<double>::infinity());
    if (nd.station != 0 && sn.classcap.size() >= nd.station) {
        const std::vector<double>& row = sn.classcap[nd.station - 1];
        for (std::size_t r = 0; r < R && r < row.size(); ++r) ccap[r] = row[r];
    }

    std::vector<std::vector<int>> nset;
    if (ntot == 0) {
        nset.push_back(std::vector<int>(R, 0));
    } else {
        const std::vector<std::vector<int>> all =
            pfqn::multichoose_rows(static_cast<int>(R), static_cast<int>(ntot));
        for (std::size_t j = 0; j < all.size(); ++j) {
            bool ok = true;
            for (std::size_t r = 0; r < R && ok; ++r)
                if (static_cast<double>(all[j][r]) > ccap[r]) ok = false;
            if (ok) nset.push_back(all[j]);
        }
    }

    std::vector<std::vector<T>> space;
    std::size_t maxw = 0;
    std::vector<std::vector<std::vector<T>>> subspaces;
    for (std::size_t j = 0; j < nset.size(); ++j) {
        std::vector<std::size_t> nj(R, 0);
        for (std::size_t r = 0; r < R; ++r) nj[r] = static_cast<std::size_t>(nset[j][r]);
        const std::vector<std::vector<T>> sj = from_marginal_node(sn, ind, nj, phases);
        if (sj.empty()) continue;
        for (std::size_t a = 0; a < sj.size(); ++a) maxw = std::max(maxw, sj[a].size());
        subspaces.push_back(sj);
    }
    for (std::size_t j = 0; j < subspaces.size(); ++j)
        for (std::size_t a = 0; a < subspaces[j].size(); ++a) {
            std::vector<T> row = subspaces[j][a];
            if (row.size() < maxw)
                row.insert(row.begin(), maxw - row.size(), num_traits<T>::from_int(0));
            space.push_back(row);
        }
    if (space.empty()) return space;

    std::sort(space.begin(), space.end());
    space.erase(std::unique(space.begin(), space.end()), space.end());
    std::reverse(space.begin(), space.end());
    return space;
}

namespace state_detail {

/**
 * Port of `multichoosecon.m`: the ways to draw `S` units from the availability
 * vector `n`.
 *
 * Unlike `multichoose_rows`, the count drawn from category r is capped by
 * `n[r]`, so the enumeration never proposes a job of a class the station does
 * not hold -- which is what keeps `from_marg_node_started` from asking the
 * per-class builder for a state that cannot exist.
 */
inline std::vector<std::vector<std::size_t> > multichoosecon(const std::vector<std::size_t>& n,
                                                             std::size_t S) {
    const std::size_t R = n.size();
    std::vector<std::vector<std::size_t> > out;
    if (R == 0) return out;
    if (S == 0) {
        out.push_back(std::vector<std::size_t>(R, 0));
        return out;
    }
    if (S == 1) {
        for (std::size_t i = 0; i < R; ++i)
            if (n[i] > 0) {
                std::vector<std::size_t> row(R, 0);
                row[i] = 1;
                out.push_back(row);
            }
        return out;
    }
    for (std::size_t i = 0; i < R; ++i) {
        if (n[i] == 0) continue;
        std::vector<std::size_t> n1 = n;
        --n1[i];
        const std::vector<std::vector<std::size_t> > tail = multichoosecon(n1, S - 1);
        for (std::size_t k = 0; k < tail.size(); ++k) {
            std::vector<std::size_t> row = tail[k];
            ++row[i];
            out.push_back(row);
        }
    }
    return out;
}

}  // namespace state_detail

/**
 * Port of `State.fromMargAndStarted`: the states with a given TOTAL queue
 * length AND a given TOTAL number of started jobs.
 *
 * Where `from_marginal_node_and_started` takes one per-class occupancy and one
 * per-class started vector and builds ONE row, this takes only the two totals
 * and returns the union of that row over every pair consistent with them:
 * `sum(n) == ntot`, `sum(s) == stot`, and `s <= n` elementwise.
 *
 * The started counts are drawn with `multichoosecon` from the jobs actually
 * present rather than from `multichoose_rows` over all classes, so a split is
 * never proposed that puts more of a class in service than the station holds.
 * Classes disabled at the station are excluded through `classcap` for the
 * reason documented on `from_marg_node`: an empty local space is ABSORBED by
 * the cartesian product instead of annihilating it, so the job would silently
 * disappear rather than the state being rejected.
 *
 * Widths are unified on the LEFT and the rows uniqued and reversed, exactly as
 * `from_marg_node` does, since the buffer is right-aligned.
 */
template <class T>
std::vector<std::vector<T>> from_marg_node_started(const NetworkStruct<T>& sn, std::size_t ind,
                                                   std::size_t ntot, std::size_t stot,
                                                   const std::vector<std::size_t>& phases) {
    if (ind == 0 || ind > sn.nodes.size())
        throw InputError("from_marg_node_started: node index " + std::to_string(ind) +
                         " is out of range");
    std::vector<std::vector<T>> space;
    if (stot > ntot) return space;

    const NodeDef& nd = sn.nodes[ind - 1];
    const std::size_t R = sn.nclasses;

    std::vector<double> ccap(R, std::numeric_limits<double>::infinity());
    if (nd.station != 0 && sn.classcap.size() >= nd.station) {
        const std::vector<double>& row = sn.classcap[nd.station - 1];
        for (std::size_t r = 0; r < R && r < row.size(); ++r) ccap[r] = row[r];
    }

    std::vector<std::vector<int>> nset;
    if (ntot == 0) {
        nset.push_back(std::vector<int>(R, 0));
    } else {
        const std::vector<std::vector<int>> all =
            pfqn::multichoose_rows(static_cast<int>(R), static_cast<int>(ntot));
        for (std::size_t j = 0; j < all.size(); ++j) {
            bool ok = true;
            for (std::size_t r = 0; r < R && ok; ++r)
                if (static_cast<double>(all[j][r]) > ccap[r]) ok = false;
            if (ok) nset.push_back(all[j]);
        }
    }

    std::size_t maxw = 0;
    std::vector<std::vector<std::vector<T>>> subspaces;
    for (std::size_t j = 0; j < nset.size(); ++j) {
        std::vector<std::size_t> nj(R, 0);
        for (std::size_t r = 0; r < R; ++r) nj[r] = static_cast<std::size_t>(nset[j][r]);
        const std::vector<std::vector<std::size_t> > sset =
            state_detail::multichoosecon(nj, stot);
        for (std::size_t k = 0; k < sset.size(); ++k) {
            const std::vector<std::vector<T>> sjk =
                from_marginal_node_and_started(sn, ind, nj, sset[k], phases);
            if (sjk.empty()) continue;
            for (std::size_t a = 0; a < sjk.size(); ++a) maxw = std::max(maxw, sjk[a].size());
            subspaces.push_back(sjk);
        }
    }
    for (std::size_t j = 0; j < subspaces.size(); ++j)
        for (std::size_t a = 0; a < subspaces[j].size(); ++a) {
            std::vector<T> row = subspaces[j][a];
            if (row.size() < maxw)
                row.insert(row.begin(), maxw - row.size(), num_traits<T>::from_int(0));
            space.push_back(row);
        }
    if (space.empty()) return space;

    std::sort(space.begin(), space.end());
    space.erase(std::unique(space.begin(), space.end()), space.end());
    std::reverse(space.begin(), space.end());
    return space;
}

/**
 * The FIRST row `from_marginal_node` emits, BUILT rather than enumerated.
 *
 * Every caller that seeds a sample path takes row 0 and discards the rest, and
 * for a Cache the rest is the ordered placement of distinct items into its
 * slots -- (nitems+1)*nitems*(nitems-1)*... rows. Materializing that to read one
 * row's contents is what took SolverSSA on tut06_cache_lru_zipf (1000 items, 50
 * slots) to 12 GB and an OOM kill of the whole host; the peak is invariant in
 * `samples` precisely because the walk itself is lazy and only the SEED costs
 * this. A simulator must not pay an exact solver's enumeration to start.
 *
 * Row 0 holds items 1..tcc in slot order, by construction and not by choice:
 * `cartesian` emits `a[0] ++ b[0]`, the contents loop offers the lowest
 * admissible item first at every slot and its duplicate test rejects the ones
 * already placed, and the delayed-hit occupancy bits start clear. So the row is
 * the per-class counts, the first `tcc` items, then zeros -- bit for bit what
 * the enumeration returned, which is what lets the CTMC find its declared
 * initial state in the space. A cache with more slots than items cannot fill
 * and starts empty, as the enumeration also has it. Every other node kind
 * delegates, its first row already being cheap.
 *
 * @return false when the node admits no state at all, as an empty return does
 */
template <class T>
bool from_marginal_node_first(const NetworkStruct<T>& sn, std::size_t ind,
                              const std::vector<std::size_t>& n,
                              const std::vector<std::size_t>& phases, std::vector<T>& out) {
    if (ind == 0 || ind > sn.nodes.size())
        throw InputError("from_marginal_node_first: node index " + std::to_string(ind) +
                         " is out of range");
    const NodeDef& nd = sn.nodes[ind - 1];
    if (nd.station != 0 || !nd.stateful || nd.nodetype != NodeType::Cache) {
        const std::vector<std::vector<T> > rows = from_marginal_node(sn, ind, n, phases);
        if (rows.empty()) return false;
        out = rows[0];
        return true;
    }

    const std::size_t R = sn.nclasses;
    if (n.size() != R)
        throw InputError("from_marginal_node_first: n must have one entry per class");
    out.clear();
    for (std::size_t r = 0; r < R; ++r) {
        const std::vector<std::vector<T> > sr = space_closed_single<T>(1, n[r]);
        // an absent placement drops the class's column, exactly as the fold does
        if (sr.empty()) {
            if (n[r] > 0) return false;
            continue;
        }
        out.insert(out.end(), sr[0].begin(), sr[0].end());
    }

    const typename std::map<std::size_t, CacheParam<T> >::const_iterator ci =
        sn.nodeparam.find(ind);
    if (ci == sn.nodeparam.end()) return true;
    const CacheParam<T>& cp = ci->second;
    std::size_t tcc = 0;
    for (std::size_t u = 0; u < cp.itemcap.size(); ++u)
        if (cp.itemcap[u] > 0) tcc += static_cast<std::size_t>(cp.itemcap[u]);
    for (std::size_t slot = 0; slot < tcc; ++slot)
        out.push_back(num_traits<T>::from_int(
            tcc <= cp.nitems ? static_cast<long>(slot + 1) : 0L));
    if (cp.retrieval_capacity > 0) {
        out.insert(out.end(), cp.nitems, num_traits<T>::from_int(0));
        std::vector<std::size_t> rcl, rci, rco;
        cache_retrieval_class_map(cp, rcl, rci, rco);
        out.insert(out.end(), rcl.size(), num_traits<T>::from_int(0));
    }
    return true;
}

/**
 * Port of `State.spaceClosedMulti`: how N[r] class-r jobs distribute over M
 * stateful nodes, for every class, as the cartesian fold of the per-class
 * placements. Column block r holds class r's counts, M columns wide.
 */
template <class T>
std::vector<std::vector<T>> space_closed_multi(
    std::size_t M, const std::vector<std::size_t>& N,
    const std::vector<std::vector<long>>& caps = std::vector<std::vector<long>>()) {
    std::vector<std::vector<T>> ss;
    for (std::size_t r = 0; r < N.size(); ++r) {
        const std::vector<std::vector<T>> sr =
            r < caps.size() && !caps[r].empty() ? space_closed_single_capped<T>(M, N[r], caps[r])
                                                : space_closed_single<T>(M, N[r]);
        if (sr.empty()) return std::vector<std::vector<T>>();
        ss = r == 0 ? sr : cartesian(ss, sr);
    }
    return ss;
}

/**
 * Port of `State.spaceClosedMultiCS`: the same, but a CHAIN's population is
 * shared among its classes, so the split between them is itself enumerated.
 *
 * Class switching moves a job between classes of one chain, so only the CHAIN
 * total is invariant. Enumerating over per-class populations alone would fix a
 * split the model does not fix, and drop every state reachable by a switch.
 *
 * @param M      number of stateful nodes (sources excluded by the caller)
 * @param N      per-class populations
 * @param chains chains[c][r] true when class r belongs to chain c
 */
template <class T>
std::vector<std::vector<T>> space_closed_multi_cs(
    std::size_t M, const std::vector<std::size_t>& N,
    const std::vector<std::vector<bool>>& chains,
    const std::vector<std::vector<long>>& caps = std::vector<std::vector<long>>()) {
    const std::size_t C = chains.size();
    const std::size_t R = N.size();
    // Per chain, the ways its total splits across the classes it contains.
    std::vector<std::vector<std::vector<int>>> chainInitPos(C);
    std::vector<std::size_t> inchain_sz(C, 0);
    for (std::size_t c = 0; c < C; ++c) {
        std::size_t tot = 0, k = 0;
        for (std::size_t r = 0; r < R; ++r)
            if (chains[c][r]) { tot += N[r]; ++k; }
        inchain_sz[c] = k;
        chainInitPos[c] = k == 0 ? std::vector<std::vector<int>>{std::vector<int>()}
                                 : pfqn::multichoose_rows(static_cast<int>(k), static_cast<int>(tot));
    }
    std::vector<std::vector<T>> ss;
    std::vector<std::size_t> v(C, 0);
    for (;;) {
        std::vector<std::size_t> subN;
        subN.reserve(R);
        for (std::size_t c = 0; c < C; ++c) {
            const std::vector<int>& row = chainInitPos[c][v[c]];
            for (std::size_t j = 0; j < row.size(); ++j)
                subN.push_back(static_cast<std::size_t>(row[j]));
        }
        // THE SPLIT IS SCATTERED BACK ONTO GLOBAL CLASS POSITIONS before the
        // per-slot caps are applied: `subN` is built chain by chain, so its
        // order is the concatenation of the chains, not the class order `caps`
        // is indexed by.
        std::vector<std::size_t> subNg(R, 0);
        {
            std::size_t k = 0;
            for (std::size_t c = 0; c < C; ++c)
                for (std::size_t r = 0; r < R; ++r)
                    if (chains[c][r]) subNg[r] = subN[k++];
        }
        const std::vector<std::vector<T>> blk =
            caps.empty() ? space_closed_multi<T>(M, subN)
                         : space_closed_multi<T>(M, subNg, caps);
        ss.insert(ss.end(), blk.begin(), blk.end());
        // odometer over the per-chain splits
        std::size_t c = 0;
        for (; c < C; ++c) {
            if (++v[c] < chainInitPos[c].size()) break;
            v[c] = 0;
        }
        if (c == C) break;
    }
    return ss;
}

/** One network state: the per-stateful-node local rows it is composed of. */
template <class T>
struct NetState {
    std::vector<std::vector<T>> local;  ///< local[isf] is that node's state row
};

/**
 * Port of `State.initialOccupancy`: the class-r jobs node `ind` holds in the
 * DECLARED initial state, or 0 when there is none.
 *
 * Only a Place answers nonzero. Every other node type encodes its local state
 * differently, so column r there is not an occupancy -- and a station whose
 * visit ratio is zero is genuinely never entered, while a Place with no input
 * arc is a TRANSIENT state of the chain rather than an absent one.
 */
template <class T>
std::size_t state_initial_occupancy(const NetworkStruct<T>& sn, std::size_t ind, std::size_t r) {
    if (ind == 0 || ind > sn.nodes.size()) return 0;
    if (sn.nodes[ind - 1].nodetype != NodeType::Place) return 0;
    const typename std::map<std::size_t, std::vector<T>>::const_iterator it =
        sn.initmarking.find(ind);
    if (it == sn.initmarking.end() || r >= it->second.size()) return 0;
    const double v = num_traits<T>::to_double(it->second[r]);
    return v > 0 ? static_cast<std::size_t>(v) : 0;
}

/**
 * Port of the `capacityc` table of `State.spaceGeneratorNodes`: the largest
 * class-r marginal node `ind` may hold in the enumerated space.
 *
 * WITHOUT THIS TABLE THE LATTICE IS PLACED BLIND. `space_closed_multi_cs`
 * spreads a class's population over EVERY stateful node, including the ones the
 * class never visits, so a mixed model enumerates open jobs at a closed-only
 * station and closed jobs at an open-only one. Those states are unreachable, but
 * they are not free: they enlarge the generator, and on `mqn_multiserver_fcfs`
 * (1904 states against the reference's 1304) they moved the Source's departure
 * rate from 0.24763 to 0.26040 -- an arrival rate the model never offers.
 *
 * The bound is the reference's, in its order: a class that does not visit the
 * node is capped at its DECLARED initial occupancy (nonzero only for an SPN
 * Place, which may hold tokens at time zero without ever being re-entered), a
 * disabled service at 0, an open class at its cutoff and a closed one at its
 * chain population, both cut down by the station's per-class buffer and then by
 * any finite-capacity region the station belongs to.
 *
 * @param sn      the network struct
 * @param cutoff  per-class population bound for open classes
 * @param cutoff_mat  optional (nstations x nclasses) bound, the reference's
 *                    matrix `options.cutoff`; overrides `cutoff` per station
 * @return (nnodes x nclasses) capacities, 0-based in both indices
 */
template <class T>
std::vector<std::vector<std::size_t>> space_capacity_c(
    const NetworkStruct<T>& sn, const std::vector<std::size_t>& cutoff,
    const std::vector<std::vector<std::size_t>>& cutoff_mat =
        std::vector<std::vector<std::size_t>>()) {
    const std::size_t R = sn.nclasses, N = sn.nodes.size();
    std::vector<std::vector<std::size_t>> cap(N, std::vector<std::size_t>(R, 0));
    const std::vector<double> njobs = sn.njobs();

    // maxPending mirrors `State.spaceGeneratorNodes`: a completing fetch
    // releases every merged secondary request in one immediate transition.
    std::size_t maxpend = 0;
    for (std::size_t r = 0; r < R; ++r)
        if (!std::isfinite(njobs[r])) maxpend = std::max(maxpend, cutoff[r]);
    if (maxpend > 0) --maxpend;

    for (std::size_t ind = 1; ind <= N; ++ind) {
        const NodeDef& nd = sn.nodes[ind - 1];
        if (!nd.stateful) continue;
        const std::size_t ist = nd.station;
        const std::size_t isf = sn.stateful_index(ind);
        if (ist != 0 && nd.nodetype != NodeType::Source) {
            for (std::size_t r = 0; r < R; ++r) {
                std::size_t c = R;  // r's chain, or R when the struct carries none
                for (std::size_t cc = 0; cc < sn.chains.size(); ++cc)
                    if (r < sn.chains[cc].size() && sn.chains[cc][r]) { c = cc; break; }
                const bool novisit = c < sn.visits.size() && sn.visits[c].rows() != 0 &&
                                     isf != 0 && isf - 1 < sn.visits[c].rows() &&
                                     num_traits<T>::to_double(sn.visits[c](isf - 1, r)) == 0;
                if (novisit) {
                    cap[ind - 1][r] = state_initial_occupancy(sn, ind, r);
                    continue;
                }
                if (nd.nodetype != NodeType::Place && ist - 1 < sn.disabled.size() &&
                    r < sn.disabled[ist - 1].size() && sn.disabled[ist - 1][r]) {
                    cap[ind - 1][r] = 0;
                    continue;
                }
                double b;
                if (!std::isfinite(njobs[r])) {
                    // `capacityc(ind,r) = min(cutoff(ist,r), classcap(ist,r))`:
                    // the reference's cutoff is a (station x class) MATRIX
                    // wherever a model needs a different truncation per station,
                    // and only the uniform case collapses to one number.
                    b = static_cast<double>(
                        (!cutoff_mat.empty() && ist - 1 < cutoff_mat.size() &&
                         r < cutoff_mat[ist - 1].size())
                            ? cutoff_mat[ist - 1][r]
                            : cutoff[r]);
                } else {
                    b = 0.0;  // the whole CHAIN's population may sit here
                    for (std::size_t k = 0; k < R; ++k)
                        if (c < sn.chains.size() && k < sn.chains[c].size() && sn.chains[c][k] &&
                            std::isfinite(njobs[k]))
                            b += njobs[k];
                }
                if (ist - 1 < sn.classcap.size() && r < sn.classcap[ist - 1].size())
                    b = std::min(b, sn.classcap[ist - 1][r]);
                if (b > 0)
                    for (std::size_t f = 0; f < sn.regions.size(); ++f) {
                        const typename NetworkStruct<T>::Region& rg = sn.regions[f];
                        if (ist - 1 >= rg.cap.size()) continue;
                        if (r < rg.cap[ist - 1].size() && rg.cap[ist - 1][r] >= 0)
                            b = std::min(b, rg.cap[ist - 1][r]);
                        if (R < rg.cap[ist - 1].size() && rg.cap[ist - 1][R] >= 0)
                            b = std::min(b, rg.cap[ist - 1][R]);
                    }
                if (!(b > 0)) b = 0;
                cap[ind - 1][r] = std::isfinite(b) ? static_cast<std::size_t>(b)
                                                   : static_cast<std::size_t>(-1);
            }
            continue;
        }
        switch (nd.nodetype) {
            case NodeType::Cache: {
                for (std::size_t r = 0; r < R; ++r) cap[ind - 1][r] = 1;
                const typename std::map<std::size_t, CacheParam<T>>::const_iterator ci =
                    sn.nodeparam.find(ind);
                if (ci != sn.nodeparam.end() && ci->second.retrieval_capacity > 0)
                    for (std::size_t r = 0; r < R; ++r)
                        for (std::size_t k = 0; k < ci->second.hitclass.size(); ++k)
                            if (ci->second.hitclass[k] == r + 1) cap[ind - 1][r] = 1 + maxpend;
                break;
            }
            case NodeType::Router:
                for (std::size_t r = 0; r < R; ++r) {
                    std::size_t c = R;
                    for (std::size_t cc = 0; cc < sn.chains.size(); ++cc)
                        if (r < sn.chains[cc].size() && sn.chains[cc][r]) { c = cc; break; }
                    cap[ind - 1][r] =
                        (c < sn.nodevisits.size() && sn.nodevisits[c].rows() != 0 &&
                         ind - 1 < sn.nodevisits[c].rows() &&
                         num_traits<T>::to_double(sn.nodevisits[c](ind - 1, r)) > 0)
                            ? 1
                            : 0;
                }
                break;
            default:
                // A Transition holds no class-indexed job, and a Source is an
                // infinite reservoir the lattice never places into.
                for (std::size_t r = 0; r < R; ++r)
                    cap[ind - 1][r] = static_cast<std::size_t>(-1);
                break;
        }
    }
    return cap;
}

/**
 * Port of `State.spaceGenerator`: every network state, reachable or not.
 *
 * The population lattice is walked with `space_closed_multi_cs`, and each
 * lattice row is expanded per node through `from_marginal`; the network states
 * are the cartesian product of those local spaces. An open class has no finite
 * population, so `cutoff` bounds it -- WITHOUT a cutoff the lattice is infinite
 * and the reference errors rather than truncating silently.
 *
 * @param sn      the network struct
 * @param cutoff  per-class population bound for open classes
 * @param maxst   refuse beyond this many states (the reference's ctmc_max_states)
 * @param cutoff_mat  optional (nstations x nclasses) bound, the reference's
 *                    matrix `options.cutoff`
 */
template <class T>
std::vector<NetState<T>> space_generator(
    const NetworkStruct<T>& sn, const std::vector<std::size_t>& cutoff,
    std::size_t maxst = 3000000,
    const std::vector<std::vector<std::size_t>>& cutoff_mat =
        std::vector<std::vector<std::size_t>>()) {
    const std::size_t R = sn.nclasses;
    if (cutoff.size() != R)
        throw InputError("space_generator: cutoff must have one entry per class");

    // Closed classes carry their own population; open ones are bounded by the
    // cutoff. A missing cutoff on an open class is an error in the reference,
    // because a silently truncated lattice yields a state space that looks
    // complete and is not.
    // `njobs()` is a population COUNT, always a plain double; it is not carried
    // in T, so reading it as a vector<T> fails to compile on every backend but
    // double -- which is how this survived until SolverCTMC was instantiated at
    // higher precision.
    const std::vector<double> njobs = sn.njobs();
    std::vector<std::size_t> Np(R, 0);
    for (std::size_t r = 0; r < R; ++r) {
        const double nj = njobs[r];
        if (std::isfinite(nj)) {
            Np[r] = static_cast<std::size_t>(nj);
        } else {
            if (cutoff[r] == 0)
                throw InputError(
                    "space_generator: class " + std::to_string(r) +
                    " is open, so its population is unbounded; supply a cutoff for it");
            Np[r] = cutoff[r];
        }
    }

    // The per-node bound the reference places the lattice under. An open class's
    // lattice height is then the LARGEST capacity any node grants it, not the raw
    // cutoff: levels above that produce only rows every node refuses.
    const std::vector<std::vector<std::size_t>> capc = space_capacity_c(sn, cutoff, cutoff_mat);
    const std::size_t unbounded = static_cast<std::size_t>(-1);
    for (std::size_t r = 0; r < R; ++r) {
        if (std::isfinite(njobs[r])) continue;
        std::size_t hi = 0;
        bool any_unbounded = false;
        for (std::size_t ind = 1; ind <= sn.nodes.size(); ++ind) {
            if (!sn.nodes[ind - 1].stateful || sn.nodes[ind - 1].nodetype == NodeType::Source)
                continue;
            if (capc[ind - 1][r] == unbounded) any_unbounded = true;
            else hi = std::max(hi, capc[ind - 1][r]);
        }
        if (!any_unbounded) Np[r] = std::min(Np[r], hi);
    }

    std::vector<std::vector<bool>> chains = sn.chains;
    if (chains.empty()) {  // no class switching: each class is its own chain
        chains.assign(R, std::vector<bool>(R, false));
        for (std::size_t r = 0; r < R; ++r) chains[r][r] = true;
    }

    // The lattice is over the STATEFUL NODES, not the stations: a Transition or
    // a Cache is stateful without being a station, and iterating stations drops
    // its local block entirely. Sources hold no jobs, so they take no column.
    const std::vector<std::size_t>& sfn = sn.stateful_nodes;
    const std::size_t NF = sfn.size();
    std::vector<std::size_t> lat_col(NF, static_cast<std::size_t>(-1));
    std::size_t Mp = 0;
    for (std::size_t k = 0; k < NF; ++k)
        if (sn.nodes[sfn[k] - 1].nodetype != NodeType::Source) lat_col[k] = Mp++;

    // The population LATTICE, not just its top row. An open class has no fixed
    // population, so every level from 0 to its cutoff is a distinct set of
    // states; only a CLOSED class is pinned to its own N. Enumerating the top
    // row alone yields a chain that can never empty -- for an M/M/1/K it left
    // exactly one state, and a one-state generator is trivially valid, which is
    // why the closed-model test could not detect this.
    std::vector<bool> is_open(R, false);
    for (std::size_t r = 0; r < R; ++r)
        is_open[r] = !std::isfinite(njobs[r]);

    // THE LATTICE IS CAPACITY-BOUND, not just the per-node local rows. `capc`
    // already zeroes a (node, class) pair the class never visits, and the
    // composition below drops the rows that breach it -- but it drops them ONE
    // AT A TIME, after `space_closed_multi_cs` has spread every class over every
    // slot. On the class-switching chain Source->Q1(A)->Q2(B)->Q3(C) that is
    // (cutoff+1)^(3*3) candidates for (cutoff+1)^3 reachable states, and the
    // rejection is not free: each candidate costs a `from_marginal`.
    // `lat_caps[class][slot]` is `capc` in the SLOT order the lattice row uses.
    std::vector<std::vector<long>> lat_caps(R, std::vector<long>(Mp, -1));
    {
        bool ok = Mp > 0;
        for (std::size_t k = 0; k < NF && ok; ++k) {
            if (lat_col[k] == static_cast<std::size_t>(-1)) continue;
            const std::size_t ind = sfn[k];
            for (std::size_t r = 0; r < R; ++r) {
                const std::size_t cap = capc[ind - 1][r];
                lat_caps[r][lat_col[k]] =
                    cap == unbounded ? -1 : static_cast<long>(cap);
            }
        }
        if (!ok) lat_caps.clear();
    }
    std::vector<std::vector<T>> pos;
    std::vector<std::size_t> nlev(R, 0);
    for (;;) {
        bool admissible = true;
        for (std::size_t r = 0; r < R; ++r)
            if (!is_open[r] && nlev[r] != Np[r]) { admissible = false; break; }
        if (admissible) {
            const std::vector<std::vector<T>> part =
                space_closed_multi_cs<T>(Mp, nlev, chains, lat_caps);
            pos.insert(pos.end(), part.begin(), part.end());
            if (pos.size() > maxst)
                throw UnsupportedError(
                    "space_generator: the population lattice exceeds the cap of " +
                    std::to_string(maxst) + " states; raise it or use another solver");
        }
        // Mixed-radix increment over the per-class levels, the reference's pprod.
        std::size_t r = 0;
        for (; r < R; ++r) {
            if (nlev[r] < Np[r]) { ++nlev[r]; break; }
            nlev[r] = 0;
        }
        if (r == R) break;
    }
    // Distinct levels can yield the same lattice row when a class is absent
    // from it, so drop the duplicates the reference removes with `unique`.
    std::sort(pos.begin(), pos.end());
    pos.erase(std::unique(pos.begin(), pos.end()), pos.end());
    // Per node and per lattice row, the local rows -- collected BEFORE any
    // composition because their width is not uniform. A marginal of 3 jobs at
    // an FCFS station needs two buffer columns where a marginal of 2 needs one,
    // and the reference unifies the two by RIGHT-ALIGNING the narrower row into
    // the widest (`fromMarginalBounds`), padding zeros on the left. Composing
    // before padding leaves an arrival unable to match its own successor,
    // because the successor is a wider vector and no index lookup can find it.
    std::vector<std::vector<std::vector<std::vector<T>>>> allper(pos.size());
    std::vector<bool> row_ok(pos.size(), true);
    std::vector<std::size_t> maxw(NF, 0);
    std::vector<NetState<T>> out;
    for (std::size_t j = 0; j < pos.size(); ++j) {
        std::vector<std::vector<std::vector<T>>> per(NF);
        bool ok = true;
        for (std::size_t k = 0; k < NF && ok; ++k) {
            const std::size_t ind = sfn[k];
            const std::size_t ist = sn.nodes[ind - 1].station;
            std::vector<std::size_t> ph(R, 1), nmarg(R, 0);
            // phases_of takes a 1-BASED class index, as elsewhere in NetworkStruct.
            if (ist != 0)
                for (std::size_t r = 0; r < R; ++r) ph[r] = sn.phasessz_of(ist, r + 1);
            if (lat_col[k] == static_cast<std::size_t>(-1)) {
                // A Source is an infinite reservoir: one job per class in
                // service, no marginal of its own -- but ONLY for the classes it
                // actually generates. A class the Source does not generate has a
                // disabled arrival, so phases_of is 0, and demanding one job of
                // it makes space_closed_single(0, 1) return nothing: the Source
                // yields no local row and the WHOLE state space comes out empty.
                // Invisible on a single-class open model, because its one class
                // is generated; a Source-Cache-Sink model, whose Hit and Miss
                // classes exist only downstream, hits it immediately.
                for (std::size_t r = 0; r < R; ++r)
                    if (!sn.disabled[ist - 1][r]) nmarg[r] = 1;
            } else {
                for (std::size_t r = 0; r < R; ++r)
                    nmarg[r] = static_cast<std::size_t>(
                        num_traits<T>::to_double(pos[j][r * Mp + lat_col[k]]));
                // `any(stateMarg_i > capacityc(ind,:))` of the reference: the node
                // yields NO row, so the whole lattice row is dropped by the
                // cartesian product below.
                for (std::size_t r = 0; r < R && ok; ++r)
                    if (capc[ind - 1][r] != static_cast<std::size_t>(-1) &&
                        nmarg[r] > capc[ind - 1][r])
                        ok = false;
                if (!ok) break;
            }
            per[k] = from_marginal_node(sn, ind, nmarg, ph);
            if (per[k].empty()) ok = false;
            for (std::size_t b = 0; b < per[k].size(); ++b)
                maxw[k] = std::max(maxw[k], per[k][b].size());
        }
        row_ok[j] = ok;
        allper[j].swap(per);
    }

    for (std::size_t j = 0; j < pos.size(); ++j) {
        if (!row_ok[j]) continue;
        std::vector<std::vector<std::vector<T>>>& per = allper[j];
        // Right-align every row into the node's widest, as the reference does.
        for (std::size_t k = 0; k < NF; ++k)
            for (std::size_t b = 0; b < per[k].size(); ++b)
                if (per[k][b].size() < maxw[k])
                    per[k][b].insert(per[k][b].begin(), maxw[k] - per[k][b].size(),
                                     num_traits<T>::from_int(0));
        // Cartesian product across nodes, carrying the per-node rows.
        std::vector<NetState<T>> acc(1);
        for (std::size_t i = 0; i < NF; ++i) {
            std::vector<NetState<T>> next;
            for (std::size_t a = 0; a < acc.size(); ++a)
                for (std::size_t b = 0; b < per[i].size(); ++b) {
                    NetState<T> ns = acc[a];
                    ns.local.push_back(per[i][b]);
                    next.push_back(ns);
                }
            acc.swap(next);
            if (out.size() + acc.size() > maxst)
                throw UnsupportedError(
                    "space_generator: the state space exceeds the cap of " +
                    std::to_string(maxst) + " states; raise it or use another solver");
        }
        out.insert(out.end(), acc.begin(), acc.end());
    }
    // The same network state can be produced by two lattice rows once the
    // widths are unified, so drop duplicates as the reference's `unique` does.
    std::vector<std::vector<std::vector<T>>> seen;
    std::vector<NetState<T>> uniq;
    for (std::size_t i = 0; i < out.size(); ++i) {
        bool dup = false;
        for (std::size_t p = 0; p < seen.size() && !dup; ++p)
            if (seen[p] == out[i].local) dup = true;
        if (dup) continue;
        seen.push_back(out[i].local);
        uniq.push_back(out[i]);
    }
    return uniq;
}

}  // namespace qn
}  // namespace line

#endif  // LINE_LANG_QN_STATE_H
