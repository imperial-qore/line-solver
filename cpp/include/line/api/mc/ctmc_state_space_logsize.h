#ifndef LINE_API_MC_CTMC_STATE_SPACE_LOGSIZE_H
#define LINE_API_MC_CTMC_STATE_SPACE_LOGSIZE_H

/**
 * @file
 * Worst-case log-size of the CTMC state space induced by a NetworkStruct.
 *
 * Port of `matlab/src/api/mc/ctmc_state_space_logsize.m`. The estimate is the
 * product of four factors, summed in log space:
 *
 *   1. job placements: stars-and-bars C(n_k+M-1, M-1) per class over the
 *      stations that keep no ordered buffer, with an open class truncated at
 *      the cutoff;
 *   2. buffer orderings: a station outside the share family keeps the CLASS
 *      SEQUENCE of the jobs it holds, so with K>1 classes one occupancy vector
 *      is as many states as its sequences;
 *   3. service phases: the phase count raised to the number of jobs that can be
 *      in service concurrently at the station;
 *   4. routing state: one pointer over the outgoing links per (node, class)
 *      routed RROBIN or WRROBIN.
 *
 * It is computed in LOG space throughout because the quantity it exists to
 * detect overflows a double: the `intractableCTMC` fixture (8 PS stations,
 * N=400, Erlang-5) sits at exp(200), and a linear-space estimator would report
 * `inf` for everything above exp(709) and lose the ability to rank one
 * intractable model against another.
 *
 * This is the quantity fed to `ctmc_memory_gate`. It is separate from the gate
 * so a caller such as SolverAUTO can screen CTMC out of a ranking without
 * building the chain.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <limits>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"

namespace line {
namespace mc {



/** Largest (m_1..m_K) box the exact ordered-buffer DP will walk. */
inline constexpr double kOrderGridMax = 1.0e6;

/**
 * Log count of (placement, ordering) configurations over ALL order-preserving
 * stations at once, POPULATION CONSERVED. `caps_per[a][k]` bounds class k at
 * ordered station a, `cap_tot[a]` bounds the buffer TOTAL there (a finite
 * station capacity is a slot count, not a per-class bound), `njobs` is the
 * population to share out, and `m_rem` share stations take the leftovers.
 *
 * Cutoff truncates an OPEN class's population in the network exactly as the
 * plain stars-and-bars term treats it, so open classes are conserved too.
 */
inline double log_ordered_joint(const std::vector<std::vector<int>>& caps_per,
                                const std::vector<double>& cap_tot,
                                const std::vector<int>& njobs, std::size_t m_rem) {
    const std::size_t Kb = njobs.size();
    std::vector<std::size_t> dims(Kb);
    std::size_t nstate = 1;
    for (std::size_t k = 0; k < Kb; ++k) { dims[k] = njobs[k] + 1; nstate *= dims[k]; }
    const double NEG = -std::numeric_limits<double>::infinity();
    auto idx_of = [&](const std::vector<int>& v) {
        std::size_t ix = 0, mult = 1;
        for (std::size_t k = 0; k < Kb; ++k) { ix += v[k] * mult; mult *= dims[k]; }
        return ix;
    };
    auto sub_of = [&](std::size_t ix) {
        std::vector<int> v(Kb);
        for (std::size_t k = 0; k < Kb; ++k) { v[k] = static_cast<int>(ix % dims[k]); ix /= dims[k]; }
        return v;
    };
    std::vector<double> L(nstate, NEG);
    L[idx_of(njobs)] = 0.0;
    for (std::size_t a = 0; a < caps_per.size(); ++a) {
        std::vector<double> Ln(nstate, NEG);
        for (std::size_t si = 0; si < nstate; ++si) {
            if (!std::isfinite(L[si])) continue;
            const std::vector<int> rem = sub_of(si);
            std::vector<int> av(Kb);
            for (std::size_t k = 0; k < Kb; ++k) av[k] = std::min(caps_per[a][k], rem[k]);
            std::vector<int> m(Kb, 0);
            for (;;) {
                int t = 0;
                for (std::size_t k = 0; k < Kb; ++k) t += m[k];
                if (!(std::isfinite(cap_tot[a]) && t > cap_tot[a])) {
                    double v = L[si] + std::lgamma(static_cast<double>(t) + 1.0);
                    for (std::size_t k = 0; k < Kb; ++k) v -= std::lgamma(static_cast<double>(m[k]) + 1.0);
                    std::vector<int> nx(Kb);
                    for (std::size_t k = 0; k < Kb; ++k) nx[k] = rem[k] - m[k];
                    const std::size_t di = idx_of(nx);
                    if (!std::isfinite(Ln[di])) Ln[di] = v;
                    else { const double mx = std::max(Ln[di], v);
                           Ln[di] = mx + std::log(std::exp(Ln[di] - mx) + std::exp(v - mx)); }
                }
                std::size_t pos = 0;
                while (pos < Kb && m[pos] == av[pos]) { m[pos] = 0; ++pos; }
                if (pos == Kb) break;
                ++m[pos];
            }
        }
        L.swap(Ln);
    }
    std::vector<double> terms;
    for (std::size_t si = 0; si < nstate; ++si) {
        if (!std::isfinite(L[si])) continue;
        const std::vector<int> rem = sub_of(si);
        double v = L[si];
        if (m_rem >= 1) {
            for (std::size_t k = 0; k < Kb; ++k) {
                const double r = rem[k];
                v += std::lgamma(r + static_cast<double>(m_rem)) - std::lgamma(r + 1.0) -
                     std::lgamma(static_cast<double>(m_rem));
            }
        } else {
            bool leftover = false;
            for (std::size_t k = 0; k < Kb; ++k) if (rem[k] > 0) leftover = true;
            if (leftover) continue;
        }
        terms.push_back(v);
    }
    if (terms.empty()) return NEG;
    double top = terms[0];
    for (double v : terms) top = std::max(top, v);
    double acc = 0.0;
    for (double v : terms) acc += std::exp(v - top);
    return top + std::log(acc);
}

/** Options the estimator reads; only the cutoff matters. */
struct CtmcSizeOptions {
    double cutoff = -1.0;  ///< < 0 or non-finite = not given, take the solver default
};

/**
 * Worst-case log state-space size of `sn`.
 *
 * @param sn  the struct, ideally after `sn_nonmarkov_toph`, since a
 *            non-Markovian service becomes phases the raw struct does not carry
 * @return    natural log of the worst-case state count
 */
template <class T>
double ctmc_state_space_logsize(const qn::NetworkStruct<T>& sn,
                                const CtmcSizeOptions& opt = CtmcSizeOptions()) {
    using lang::RoutingStrategy;
    using lang::SchedStrategy;

    const std::size_t M = sn.nstations;
    const std::size_t K = sn.nclasses;
    if (M == 0 || K == 0) return 0.0;

    // The analyzer's own default for an open or mixed model, so the estimate is
    // taken against the space that would actually be built.
    double cutoff = opt.cutoff;
    if (!(cutoff > 0.0) || !std::isfinite(cutoff))
        cutoff = std::ceil(std::pow(6000.0, 1.0 / static_cast<double>(M * K)));

    // ORDERED BUFFERS, computed EXACTLY -- see the python twin for the measured
    // numbers. Omitting it under-priced gallery_mmap1_multiclass 311x; BOUNDING
    // it instead of computing it double counts and refused a working model.
    const auto is_share = [](SchedStrategy s) {
        return s == SchedStrategy::INF || s == SchedStrategy::PS || s == SchedStrategy::DPS ||
               s == SchedStrategy::GPS || s == SchedStrategy::PSPRIO ||
               s == SchedStrategy::DPSPRIO || s == SchedStrategy::GPSPRIO ||
               s == SchedStrategy::LPS;
    };
    std::vector<bool> is_buffered(K, true);
    std::size_t Kb = 0;
    for (std::size_t k = 0; k < K; ++k) {
        if (k < sn.issignal.size() && sn.issignal[k]) is_buffered[k] = false;
        if (is_buffered[k]) ++Kb;
    }
    std::size_t n_ord = 0;
    if (Kb > 1) {
        for (std::size_t i = 0; i < M; ++i) {
            const SchedStrategy s = sn.stations[i].sched;
            if (s == SchedStrategy::EXT || is_share(s)) continue;
            ++n_ord;
        }
    }

    double log_nstates = 0.0;
    std::vector<double> nk_eff(K, 0.0);
    const std::vector<double>& njobs = sn.njobs();
    for (std::size_t k = 0; k < K; ++k)
        nk_eff[k] = (k < njobs.size() && std::isfinite(njobs[k])) ? njobs[k] : cutoff;

    const auto place = [&](double nk, double ms) {
        return std::lgamma(1.0 + nk + ms - 1.0) - std::lgamma(1.0 + ms - 1.0) -
               std::lgamma(1.0 + nk);
    };
    // A ZERO per-class capacity means the class is DISABLED at that station, so it
    // never occupies a slot and the placement term must spread it over the
    // stations that admit it, not over all M. ld_whittle_bandwidth disables each
    // of its three PS routes for the other two classes; counting all M=4 priced it
    // at C(9,6)^3 = 592704 states, 7852 GB under the quadratic byte model, and the
    // gate refused a model whose true space is 7^3 = 343 and solves at once.
    const auto disabled_at = [&](std::size_t i, std::size_t k) {
        return i < sn.classcap.size() && k < sn.classcap[i].size() && sn.classcap[i][k] == 0.0;
    };
    const auto admitting_all = [&](std::size_t k) {
        std::size_t n = 0;
        for (std::size_t i = 0; i < M; ++i)
            if (!disabled_at(i, k)) ++n;
        return n < 1 ? std::size_t(1) : n;
    };
    const auto admitting_rem = [&](std::size_t k) {
        std::size_t n = 0;
        for (std::size_t i = 0; i < M; ++i) {
            const SchedStrategy sc = sn.stations[i].sched;
            if (!(sc == SchedStrategy::EXT || is_share(sc))) continue;
            if (!disabled_at(i, k)) ++n;
        }
        return n;
    };
    if (n_ord == 0) {
        for (std::size_t k = 0; k < K; ++k)
            log_nstates += place(nk_eff[k], static_cast<double>(admitting_all(k)));
    } else {
        const std::size_t m_rem = M - n_ord;
        for (std::size_t k = 0; k < K; ++k)
            if (!is_buffered[k])
                log_nstates += place(nk_eff[k], static_cast<double>(admitting_all(k)));
        std::vector<int> caps;
        double grid = 1.0;
        for (std::size_t k = 0; k < K; ++k)
            if (is_buffered[k]) {
                caps.push_back(static_cast<int>(std::floor(nk_eff[k])));
                grid *= (std::floor(nk_eff[k]) + 1.0);
            }
        if (grid <= kOrderGridMax) {
            std::vector<std::vector<int>> caps_per;
            std::vector<double> cap_tot;
            std::vector<int> njb;
            for (std::size_t k = 0; k < K; ++k)
                if (is_buffered[k]) njb.push_back(static_cast<int>(std::floor(nk_eff[k])));
            for (std::size_t i = 0; i < M; ++i) {
                const SchedStrategy sc = sn.stations[i].sched;
                if (sc == SchedStrategy::EXT || is_share(sc)) continue;
                std::vector<int> per;
                std::size_t bi = 0;
                for (std::size_t k = 0; k < K; ++k)
                    if (is_buffered[k]) {
                        int c = njb[bi];
                        if (i < sn.classcap.size() && k < sn.classcap[i].size() &&
                            std::isfinite(sn.classcap[i][k]))
                            c = std::min(c, static_cast<int>(std::floor(sn.classcap[i][k])));
                        per.push_back(c);
                        ++bi;
                    }
                caps_per.push_back(per);
                const double c = sn.stations[i].cap;
                cap_tot.push_back(std::isfinite(c) && c >= 0
                                      ? std::floor(c) : std::numeric_limits<double>::infinity());
            }
            log_nstates += log_ordered_joint(caps_per, cap_tot, njb, m_rem);
        } else {
            double total = 0.0;
            for (int c : caps) total += c;
            const double lkb = std::log(static_cast<double>(Kb));
            log_nstates += static_cast<double>(n_ord) *
                ((total + 1.0) * lkb - std::log(static_cast<double>(Kb) - 1.0) +
                 std::log1p(-std::exp(-(total + 1.0) * lkb)));
            for (std::size_t k = 0; k < K; ++k)
                if (is_buffered[k]) {
                    const std::size_t mk = admitting_rem(k);
                    if (mk >= 1) log_nstates += place(nk_eff[k], static_cast<double>(mk));
                }
        }
    }

    // A sharing discipline can hold every job in service at once; a queueing one
    // holds at most its server count.
    for (std::size_t i = 1; i <= M; ++i) {
        const SchedStrategy sched = sn.stations[i - 1].sched;
        const bool share = sched == SchedStrategy::INF || sched == SchedStrategy::PS ||
                           sched == SchedStrategy::DPS || sched == SchedStrategy::GPS ||
                           sched == SchedStrategy::PSPRIO || sched == SchedStrategy::DPSPRIO ||
                           sched == SchedStrategy::GPSPRIO || sched == SchedStrategy::LPS;
        for (std::size_t r = 1; r <= K; ++r) {
            const double p = static_cast<double>(sn.phasessz_of(i, r));
            if (!std::isfinite(p) || p <= 1.0) continue;
            double m;
            if (sched == SchedStrategy::EXT)
                m = 1.0;
            else if (share)
                m = nk_eff[r - 1];
            else
                m = std::min(nk_eff[r - 1], sn.stations[i - 1].nservers);
            if (!std::isfinite(m)) m = nk_eff[r - 1];
            log_nstates += std::lgamma(1.0 + m + p - 1.0) - std::lgamma(1.0 + p - 1.0) -
                           std::lgamma(1.0 + m);
        }
    }

    // Round-robin routing is stateful: the pointer over the outgoing links is
    // part of the state. The reference reads `sn.connmatrix`; this port has no
    // such field and takes the out-degree from `rtnodes`, the same graph after
    // the refresh resolved the routing strategies.
    for (std::size_t ind = 1; ind <= sn.nodes.size(); ++ind) {
        const std::vector<RoutingStrategy>& rt_i = sn.nodes[ind - 1].routing;
        std::size_t nrr = 0;
        for (std::size_t r = 0; r < K && r < rt_i.size(); ++r)
            if (rt_i[r] == RoutingStrategy::RROBIN || rt_i[r] == RoutingStrategy::WRROBIN) ++nrr;
        if (nrr == 0) continue;
        const std::size_t nout = sn.downstream_stations(ind).size();
        if (nout <= 1) continue;
        log_nstates += static_cast<double>(nrr) * std::log(static_cast<double>(nout));
    }

    return log_nstates;
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_STATE_SPACE_LOGSIZE_H
