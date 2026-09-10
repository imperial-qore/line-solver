/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Finite Capacity Regions in SolverCTMC: the DROP rule, as a filter on the
 * enumerated state space, and the gate that refuses WAITQ.
 *
 * WHY DROP IS A FILTER AND WAITQ IS NOT. Under DROP a job refused entry to a
 * full region is LOST, so the chain simply never occupies a state violating the
 * region's caps: removing those states and letting `make_infgen` re-close the
 * rows is exactly the censored chain, and it is what the reference does (its
 * `spaceGeneratorNodes` bound plus the post-generation filter). Under WAITQ the
 * job leaves its station and PARKS in a per-region FIFO outside every station,
 * to be released head-of-line as capacity frees; that queue is extra state the
 * region owns, so the state vector has to be augmented with one token buffer
 * per region and the transition relation rebuilt around it. The two are not
 * variants of one mechanism, which is why only the first is here and the second
 * is refused BY NAME rather than approximated by dropping.
 *
 * THE REGION BOUND ALSO BELONGS IN THE ENUMERATION, not only after it. The
 * automatic cutoff is region-blind and can far exceed any reachable population
 * -- cutoff 10 for a region capped at 4 -- and filtering only after generation
 * means enumerating an intractable space first. Filtering here is correct but
 * not sufficient for large models; see the note on `ctmc_region_class_cap`.
 */
#ifndef LINE_SOLVERS_CTMC_SOLVER_CTMC_FCR_H
#define LINE_SOLVERS_CTMC_SOLVER_CTMC_FCR_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/state.h"
#include "line/solvers/ctmc/solver_ctmc.h"
#include "line/util/error.h"

namespace line {
namespace ctmc {

using lang::DropStrategy;

/**
 * Refuse the region rules this port does not implement.
 *
 * WAITQ needs the per-region token FIFO described above; BAS, BBS and RSRD are
 * blocking rules whose held-job marker this generator does not carry either.
 * Each is named so a caller learns which rule stopped it rather than seeing a
 * region silently behave as DROP.
 */
template <class T>
void ctmc_check_region_rules(const NetworkStruct<T>& sn) {
    for (std::size_t f = 0; f < sn.regions.size(); ++f)
        for (std::size_t r = 0; r < sn.regions[f].rule.size(); ++r) {
            const DropStrategy d = sn.regions[f].rule[r];
            if (d == DropStrategy::DROP) continue;
            const char* nm = d == DropStrategy::WAITQ  ? "WAITQ"
                             : d == DropStrategy::BAS  ? "BAS"
                             : d == DropStrategy::BBS  ? "BBS"
                             : d == DropStrategy::RSRD ? "RSRD"
                                                       : "an unknown rule";
            throw UnsupportedError(
                std::string("SolverCTMC: finite capacity region ") + std::to_string(f + 1) +
                " applies " + nm + " to class " + std::to_string(r + 1) +
                ", which augments the state with a per-region waiting queue (WAITQ) or a held-job "
                "marker (BAS/BBS/RSRD); only DROP is ported, and it is a filter on the state "
                "space rather than extra state");
        }
}

/**
 * True when `nir` -- the per-(station, class) counts of one state, in
 * `(ist-1)*K + k` order -- satisfies every region.
 *
 * The caps are compared against the counts SUMMED OVER THE MEMBER STATIONS,
 * which is the whole point of a region: a per-station cap cannot express "at
 * most 6 jobs between these three stations". The reference's -1 sentinel means
 * unbounded and is skipped rather than compared.
 */
template <class T>
bool ctmc_region_admissible(const NetworkStruct<T>& sn, const std::vector<T>& nir) {
    const std::size_t M = sn.stations.size(), K = sn.nclasses;
    for (std::size_t f = 0; f < sn.regions.size(); ++f) {
        const typename NetworkStruct<T>::Region& rg = sn.regions[f];
        double total = 0, memory = 0;
        std::vector<double> per_class(K, 0.0);
        double gcap = -1.0, memcap = -1.0;
        std::vector<double> ccap(K, -1.0);
        bool any_member = false;
        for (std::size_t i = 0; i < M; ++i) {
            if (i >= rg.members.size() || !rg.members[i]) continue;
            any_member = true;
            for (std::size_t k = 0; k < K; ++k) {
                const double n = num_traits<T>::to_double(nir[i * K + k]);
                if (!std::isfinite(n)) continue;  // a Source's Inf sentinel
                per_class[k] += n;
                total += n;
                memory += n * num_traits<T>::to_double(rg.size[k]);
            }
            // The caps are replicated on every member row, so the first member
            // carries them; taking the tightest guards a hand-built struct.
            for (std::size_t k = 0; k < K; ++k)
                if (rg.cap[i][k] != -1.0)
                    ccap[k] = ccap[k] == -1.0 ? rg.cap[i][k] : std::min(ccap[k], rg.cap[i][k]);
            if (rg.cap[i][K] != -1.0)
                gcap = gcap == -1.0 ? rg.cap[i][K] : std::min(gcap, rg.cap[i][K]);
            if (rg.maxmem[i] != -1.0)
                memcap = memcap == -1.0 ? rg.maxmem[i] : std::min(memcap, rg.maxmem[i]);
        }
        if (!any_member) continue;
        if (gcap != -1.0 && total > gcap + 1e-9) return false;
        if (memcap != -1.0 && memory > memcap + 1e-9) return false;
        for (std::size_t k = 0; k < K; ++k)
            if (ccap[k] != -1.0 && per_class[k] > ccap[k] + 1e-9) return false;
        // The linear constraint A n <= b, evaluated on the region-wide counts.
        for (std::size_t row = 0; row < rg.lincon_A.rows() && row < rg.lincon_b.size(); ++row) {
            double lhs = 0;
            for (std::size_t k = 0; k < K && k < rg.lincon_A.cols(); ++k)
                lhs += num_traits<T>::to_double(rg.lincon_A(row, k)) * per_class[k];
            if (lhs > num_traits<T>::to_double(rg.lincon_b[row]) + 1e-9) return false;
        }
    }
    return true;
}

/**
 * The states of `space` a DROP region admits, in their original order.
 *
 * Order is preserved so the caller can restrict the arrival and departure rate
 * arrays with the same index set; reordering here would silently misalign them.
 */
template <class T>
std::vector<NetState<T>> ctmc_filter_regions(const NetworkStruct<T>& sn,
                                             const std::vector<NetState<T>>& space) {
    if (sn.regions.empty()) return space;
    ctmc_check_region_rules(sn);
    const Matrix<T> A = ctmc_state_space_aggr(sn, space);
    std::vector<NetState<T>> out;
    for (std::size_t s = 0; s < space.size(); ++s) {
        std::vector<T> nir(A.cols());
        for (std::size_t c = 0; c < A.cols(); ++c) nir[c] = A(s, c);
        if (ctmc_region_admissible(sn, nir)) out.push_back(space[s]);
    }
    if (out.empty())
        throw UnsupportedError(
            "SolverCTMC: no state satisfies the finite capacity regions; check that the region "
            "caps admit the model's population");
    return out;
}

/**
 * True where a class sits at a station inside a DROP region, per station.
 *
 * `solver_ctmc_avg_from_pi` needs it for the same reason it needs a finite
 * capacity: a job that can be dropped never entered service, so the offered
 * arrival rate is not what the server did and only the carried rate is a
 * utilization.
 */
template <class T>
std::vector<bool> ctmc_in_drop_region(const NetworkStruct<T>& sn) {
    std::vector<bool> in(sn.stations.size(), false);
    for (std::size_t f = 0; f < sn.regions.size(); ++f) {
        bool has_drop = false;
        for (std::size_t r = 0; r < sn.regions[f].rule.size(); ++r)
            if (sn.regions[f].rule[r] == DropStrategy::DROP) has_drop = true;
        if (!has_drop) continue;
        for (std::size_t i = 0; i < in.size() && i < sn.regions[f].members.size(); ++i)
            if (sn.regions[f].members[i]) in[i] = true;
    }
    return in;
}

}  // namespace ctmc
}  // namespace line

#endif  // LINE_SOLVERS_CTMC_SOLVER_CTMC_FCR_H
