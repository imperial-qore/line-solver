/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_TO_QRF_BLOCKING_H
#define LINE_API_SN_SN_TO_QRF_BLOCKING_H

/**
 * The QRF BAS blocking tables (f, MR, BB, MM, ZZ, MM1), derived from an sn.
 *
 * Port of matlab/src/api/sn/sn_to_qrf_capacity.m and
 * matlab/src/api/sn/sn_to_qrf_blocking.m.
 *
 * `qrf_bas` describes a Blocking-After-Service network by a finite-capacity
 * queue f and an enumeration of the BLOCKING CONFIGURATIONS reachable behind
 * it. Everything in that enumeration is implied by the model, so it is derived
 * here rather than demanded from the caller; `options.config.qrf_params`
 * remains an explicit override.
 *
 * The tables, and the constraint that reads each one in `qrf_bas`:
 *
 *   f         the ONE finite-capacity queue. The formulation carries a scalar f
 *             (ZERO4/ZERO7/ZERO8, THM30, THM3I, THM3L all index it), so a model
 *             with two binding buffers is refused here.
 *   F(i)      min(buffer size, N) for every queue; N where the buffer is
 *             unbounded, since no queue can hold more than the population.
 *   BB(m,i)   1 iff queue i is blocked in configuration m.
 *   ZZ(m)     the blocking depth of configuration m.
 *   MM(m,0)   head of the FIFO blocking order: the queue that takes the slot
 *             when f completes. `qrf_bas` reads ONLY the first column.
 *   MM1(m,j)  index of the configuration reached from m when j becomes blocked.
 *             Read by THM3L alone, at depth ZM-1.
 *
 * THREE INVARIANTS, each a correctness condition rather than a convention:
 *
 *   1. Configuration 1 MUST be the empty one. ZERO4 iterates `m = 2:MR` and
 *      ZERO5/ZERO7/ZERO8 test `m >= 2` to mean "some queue is blocked".
 *   2. ZM = max(ZZ) MUST be the reachable maximum. `qrf_bas` recomputes ZM from
 *      ZZ and closes the depth ladder there, so a truncated enumeration excises
 *      states the real chain visits and the polytope stops containing the true
 *      distribution -- the bound stops bounding. The size guard therefore
 *      REFUSES; it never truncates.
 *   3. Blocking APPENDS at the tail: a queue that becomes blocked joins behind
 *      those already waiting, so MM1's successor is the configuration with j
 *      appended, and the head MM(m,0) names never moves.
 *
 * THE ENUMERATION IS THE FULL ORDERED ONE, and it has to be. A (set, head)
 * collapse looks sound -- the LP reads configurations only through BB, ZZ,
 * MM(:,0) and MM1, and both objective and readout sum over m -- and it would
 * shrink MR from sum_z P(B,z) to 1 + sum_z C(B,z)*z. It was tried and it is
 * WRONG. Merging the depth-ZM configurations that share a set and a head makes
 * several THM3L rows, one per depth-(ZM-1) predecessor, reference the SAME
 * merged successor block. That is extra coupling the fine system does not have,
 * so the collapsed polytope is strictly SMALLER, not a projection of the fine
 * one, and it can cut off the true distribution. Measured on a 4-station model
 * with three feeders (B=3, ZM=3, MR 13 collapsed vs 16 full), the collapse
 * reported upper bounds of 0.681/0.979/0.768 where the full enumeration gives
 * 0.709/0.982/0.800: tighter, from a coarser state space, which is the
 * signature of a cut that is not valid.
 *
 * So MR is factorial in the number of feeders B, and the size guard is what
 * keeps that honest: it REFUSES an oversized instance rather than trimming the
 * enumeration, because trimming is the same unsound cut by another name.
 *
 * WHO CAN BE BLOCKED is read from `sn.isbasblocking`, not from `sn.droprule`:
 * LINE accepts the BAS declaration on the upstream station or on the full
 * destination, and reading droprule at the capped station sees only the second
 * (BUG-83).
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "line/api/sn/sn_get_buffer_size.h"
#include "line/lang/qn/network_struct.h"

namespace line {
namespace sn {

/**
 * Variable-count ceiling of the derived LP. A guard, not a tuning knob: the
 * enumeration cannot be truncated (invariant 2), so an oversized model is
 * refused rather than approximated.
 */
const double kQrfDefaultMaxVars = 5e5;

/** Per-station occupancy bounds, and which of them bind. */
struct QrfCapacity {
    std::vector<int> F;         ///< (M) occupancy bound of each station, in jobs
    std::vector<bool> binding;  ///< (M) true where the buffer can refuse a job
    std::string msg;            ///< empty on success
};

/**
 * F is an OCCUPANCY BOUND, not a declared capacity: the station's buffer where
 * that buffer BINDS, and the population N everywhere else, since no queue of a
 * closed model can hold more than N jobs. Binding is decided by
 * `sn_get_buffer_size`, the single place in LINE that makes that call.
 *
 * Both QRF blocking bounds need this. `qrf.bas` needs it beside the blocking
 * tables; `qrf.rsrd` needs it ALONE, since its PBB constraint reads only which
 * queues can be full and it carries no blocking tables at all.
 */
template <class T>
QrfCapacity sn_to_qrf_capacity(const qn::NetworkStruct<T>& sn) {
    QrfCapacity out;
    const std::size_t M = sn.nstations;
    out.F.assign(M, 0);
    out.binding.assign(M, false);

    double Nd = 0.0;
    const std::vector<double> njobs = sn.njobs();
    for (std::size_t r = 0; r < njobs.size(); ++r) Nd += njobs[r];
    if (!(Nd >= 1.0) || std::isinf(Nd)) {
        out.msg = "the QRF bounds need a closed model with a finite population.";
        return out;
    }
    const int N = static_cast<int>(Nd + 0.5);

    for (std::size_t i = 0; i < M; ++i) {
        const double b = sn_get_buffer_size(sn, i + 1);  // 1-based
        out.binding[i] = std::isfinite(b);
        out.F[i] = (!std::isfinite(b) || b > static_cast<double>(N)) ? N
                                                                    : static_cast<int>(b + 0.5);
        if (out.F[i] < 1) {
            std::ostringstream os;
            os << "station " << (i + 1) << " has capacity " << out.F[i]
               << ": the QRF bounds need every queue to be able to hold at least one job.";
            out.msg = os.str();
            return out;
        }
    }
    return out;
}

/** The derived blocking tables, in the reference's 1-based queue indexing. */
struct QrfBlocking {
    int f = 1;                            ///< finite-capacity queue, 1-based
    std::vector<int> F;                   ///< (M) occupancy bounds
    int MR = 1;                           ///< number of blocking configurations
    std::vector<std::vector<int> > BB;    ///< (MR x M) blocking state
    std::vector<std::vector<int> > MM;    ///< (MR x .) blocking order, 1-based, 0 = absent
    std::vector<std::vector<int> > MM1;   ///< (MR x M) successor map, 1-based, 0 = absent
    std::vector<int> ZZ;                  ///< (MR) blocked count per configuration
    int ZM = 0;                           ///< maximum reachable blocking depth
    std::vector<int> blockers;            ///< 1-based stations that can be blocked behind f
    std::string msg;                      ///< empty on success
};

namespace detail {

/** The one-configuration table for a model in which no blocking is reachable. */
inline QrfBlocking qrf_empty_blocking(const std::vector<int>& F, int f_one_based,
                                      std::size_t M) {
    QrfBlocking b;
    b.f = f_one_based;
    b.F = F;
    b.MR = 1;
    b.BB.assign(1, std::vector<int>(M, 0));
    b.MM.assign(1, std::vector<int>(2, 0));
    b.MM1.assign(1, std::vector<int>(M, 0));
    b.ZZ.assign(1, 0);
    b.ZM = 0;
    return b;
}

/**
 * 0-based station indices that hold a completed job when f is full. A blocker
 * must route into f, must not be f, and must not be an infinite server (which
 * has a server per job and cannot be held). BAS itself is read from the BUG-83
 * marker, with a structural fallback for an sn built without it.
 */
template <class T>
std::vector<std::size_t> qrf_blockers(const qn::NetworkStruct<T>& sn, std::size_t f,
                                      std::size_t M) {
    std::vector<bool> declared(M, false);
    bool any = false;
    for (std::size_t i = 0; i < M; ++i) {
        const std::size_t ind = (i + 1 <= sn.station_to_node.size()) ? sn.station_to_node[i] : 0;
        if (ind >= 1 && ind <= sn.isbasblocking.size() && sn.isbasblocking[ind - 1]) {
            declared[i] = true;
            any = true;
        }
    }
    if (!any && !sn.droprule.empty()) {
        // Fallback: BAS declared on the upstream station or on the full
        // destination, the same two forms refresh_bas_blocking resolves.
        bool dest_bas = false;
        if (sn.droprule.size() > f)
            for (std::size_t r = 0; r < sn.droprule[f].size(); ++r)
                if (sn.droprule[f][r] == qn::DropStrategy::BAS) dest_bas = true;
        for (std::size_t i = 0; i < M; ++i) {
            if (i == f || i >= sn.droprule.size()) continue;
            bool here_bas = false;
            for (std::size_t r = 0; r < sn.droprule[i].size(); ++r)
                if (sn.droprule[i][r] == qn::DropStrategy::BAS) here_bas = true;
            if (dest_bas || here_bas) declared[i] = true;
        }
    }

    const std::size_t R = sn.nclasses;
    std::vector<std::size_t> out;
    for (std::size_t i = 0; i < M; ++i) {
        if (i == f || !declared[i]) continue;
        if (sn.stations[i].sched == qn::SchedStrategy::INF) continue;
        // Summed over class pairs so the test survives a multiclass sn, even
        // though the QRF gate upstream admits one class only.
        bool routes = false;
        for (std::size_t r = 0; r < R && !routes; ++r)
            for (std::size_t s = 0; s < R; ++s)
                if (sn.rt(i * R + r, f * R + s) > 0) {
                    routes = true;
                    break;
                }
        if (routes) out.push_back(i);
    }
    std::sort(out.begin(), out.end());
    return out;
}

/**
 * Every ordered sequence of distinct blockers up to length ZM, first entry the
 * head. Deterministic order -- depth ascending, then subsets lexicographic by
 * ascending station index, then the orders of each subset sorted -- so every
 * codebase emits identical tables. The empty configuration is first
 * (invariant 1).
 */
inline std::vector<std::vector<std::size_t> > qrf_enumerate_permutations(
    const std::vector<std::size_t>& blockers, int ZM) {
    std::vector<std::vector<std::size_t> > cfg;
    cfg.push_back(std::vector<std::size_t>());
    const std::size_t nb = blockers.size();
    for (int z = 1; z <= ZM; ++z) {
        std::vector<bool> pick(nb, false);
        for (std::size_t i = 0; i < static_cast<std::size_t>(z) && i < nb; ++i) pick[i] = true;
        // iterate subsets in lexicographic order of ascending index
        std::vector<std::vector<std::size_t> > subsets;
        std::vector<std::size_t> idx(z, 0);
        for (int d = 0; d < z; ++d) idx[d] = d;
        while (true) {
            std::vector<std::size_t> members;
            for (int d = 0; d < z; ++d) members.push_back(blockers[idx[d]]);
            subsets.push_back(members);
            int d = z - 1;
            while (d >= 0 && idx[d] == nb - z + d) --d;
            if (d < 0) break;
            ++idx[d];
            for (int e = d + 1; e < z; ++e) idx[e] = idx[e - 1] + 1;
        }
        for (std::size_t s = 0; s < subsets.size(); ++s) {
            std::vector<std::size_t> order = subsets[s];  // already ascending
            do {
                cfg.push_back(order);
            } while (std::next_permutation(order.begin(), order.end()));
        }
    }
    return cfg;
}

/**
 * Identity of a configuration: the whole blocking order, since that is what
 * distinguishes configurations in the enumeration `qrf_bas` is entitled to.
 */
inline std::string qrf_cfg_key(const std::vector<std::size_t>& seq) {
    std::ostringstream os;
    for (std::size_t i = 0; i < seq.size(); ++i) os << seq[i] << ',';
    return os.str();
}

/** Total service phases across stations, which sizes the QRF variable space. */
template <class T>
int qrf_total_phases(const std::vector<std::pair<Matrix<T>, Matrix<T> > >& MAPs) {
    int total = 0;
    for (std::size_t i = 0; i < MAPs.size(); ++i)
        total += std::max<int>(1, static_cast<int>(MAPs[i].first.rows()));
    return total;
}

}  // namespace detail

/**
 * @param sn       network structure
 * @param Ktot     total service phases, which sizes the LP together with MR and N
 * @param max_vars variable-count ceiling; pass kQrfDefaultMaxVars for the default
 * @return the derived tables, or a QrfBlocking carrying a non-empty msg
 */
template <class T>
QrfBlocking sn_to_qrf_blocking(const qn::NetworkStruct<T>& sn, int Ktot,
                               double max_vars = kQrfDefaultMaxVars) {
    QrfBlocking out;
    const std::size_t M = sn.nstations;

    double Nd = 0.0;
    const std::vector<double> njobs = sn.njobs();
    for (std::size_t r = 0; r < njobs.size(); ++r) Nd += njobs[r];
    const int N = static_cast<int>(Nd + 0.5);

    const QrfCapacity cap = sn_to_qrf_capacity(sn);
    if (!cap.msg.empty()) {
        out.msg = cap.msg;
        return out;
    }

    std::vector<std::size_t> fcand;
    for (std::size_t i = 0; i < M; ++i)
        if (cap.binding[i]) fcand.push_back(i);

    if (fcand.empty()) {
        // No binding buffer: callers gate on sn_has_blocking first, so this is
        // a defensive branch rather than a normal path.
        return detail::qrf_empty_blocking(cap.F, 1, M);
    }
    if (fcand.size() > 1) {
        std::ostringstream os;
        os << "'qrf.bas' models a single finite-capacity queue (its f is a scalar), but "
           << fcand.size() << " stations have a binding buffer: ";
        for (std::size_t c = 0; c < fcand.size(); ++c) {
            if (c) os << ", ";
            os << sn.stations[fcand[c]].name;
        }
        os << ". Use 'qrf.rsrd', whose PBB constraint sums over every full queue and therefore "
              "admits several, or cap only one station.";
        out.msg = os.str();
        return out;
    }
    const std::size_t f = fcand[0];

    const std::vector<std::size_t> blockers = detail::qrf_blockers(sn, f, M);

    // Blocking needs f at capacity plus one held job per blocked queue, so the
    // population caps the depth as tightly as the feeder count does.
    int ZM = std::min<int>(static_cast<int>(blockers.size()), N - cap.F[f]);
    if (ZM < 0) ZM = 0;
    if (ZM == 0) {
        QrfBlocking b = detail::qrf_empty_blocking(cap.F, static_cast<int>(f) + 1, M);
        for (std::size_t i = 0; i < blockers.size(); ++i)
            b.blockers.push_back(static_cast<int>(blockers[i]) + 1);
        return b;
    }

    const std::vector<std::vector<std::size_t> > cfg =
        detail::qrf_enumerate_permutations(blockers, ZM);
    const int MR = static_cast<int>(cfg.size());

    // Size guard: refuse, never truncate (invariant 2).
    const double n_vars =
        static_cast<double>(MR) * (N + 1) * (N + 1) * Ktot * Ktot + Ktot;
    if (n_vars > max_vars) {
        std::ostringstream os;
        os << "the QRF BAS linear program for this model would carry " << n_vars
           << " variables (MR=" << MR << " blocking configurations, N=" << N << ", " << Ktot
           << " service phases in total), above the qrf_maxvars limit of " << max_vars
           << ". The enumeration cannot be truncated -- a depth below the reachable maximum ZM="
           << ZM << " excises states the chain visits, and the result would no longer bound. "
           << "Reduce the population, the number of stations feeding " << sn.stations[f].name
           << ", or the phase counts; or raise the limit deliberately.";
        out.msg = os.str();
        return out;
    }

    out.f = static_cast<int>(f) + 1;
    out.F = cap.F;
    out.MR = MR;
    out.ZM = ZM;
    out.BB.assign(MR, std::vector<int>(M, 0));
    out.MM.assign(MR, std::vector<int>(std::max<std::size_t>(2, blockers.size()), 0));
    out.MM1.assign(MR, std::vector<int>(M, 0));
    out.ZZ.assign(MR, 0);
    for (std::size_t i = 0; i < blockers.size(); ++i)
        out.blockers.push_back(static_cast<int>(blockers[i]) + 1);

    std::map<std::string, int> index;
    for (int m = 0; m < MR; ++m) {
        const std::vector<std::size_t>& seq = cfg[m];
        out.ZZ[m] = static_cast<int>(seq.size());
        for (std::size_t z = 0; z < seq.size(); ++z) {
            out.BB[m][seq[z]] = 1;
            // only the first column is read; the rest records the full order
            out.MM[m][z] = static_cast<int>(seq[z]) + 1;
        }
        index[detail::qrf_cfg_key(seq)] = m;
    }
    for (int m = 0; m < MR; ++m) {
        if (out.ZZ[m] >= ZM) continue;  // THM3L reads MM1 only below ZM
        for (std::size_t b = 0; b < blockers.size(); ++b) {
            const std::size_t j = blockers[b];
            if (out.BB[m][j]) continue;
            std::vector<std::size_t> succ = cfg[m];
            succ.push_back(j);
            const std::map<std::string, int>::const_iterator it =
                index.find(detail::qrf_cfg_key(succ));
            if (it != index.end()) out.MM1[m][j] = it->second + 1;
        }
    }
    return out;
}

}  // namespace sn
}  // namespace line

#endif  // LINE_API_SN_SN_TO_QRF_BLOCKING_H
