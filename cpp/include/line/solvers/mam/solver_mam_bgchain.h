/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MAM_SOLVER_MAM_BGCHAIN_H
#define LINE_SOLVERS_MAM_SOLVER_MAM_BGCHAIN_H

/**
 * Port of `solver_mam_bgchain.m` and its three helpers: the analyzer that
 * treats the CLOSED classes as a background modulating chain and the OPEN
 * classes as matrix-analytic queues driven by it. A purely CLOSED model is the
 * degenerate case of the same construction -- with no open work to take a share
 * of the servers the chain alone answers, and it answers with the EXACT closed
 * CTMC at chain granularity -- so only a purely OPEN model is refused.
 *
 * WHY IT EXISTS. `dec.source` replaces a closed chain by a Poisson surrogate at
 * the current throughput iterate, which is exactly the part of a mixed model it
 * has least information about; measured against SolverCTMC that costs it 10-24%
 * on the queue lengths. But the closed population vector is the one part of a
 * mixed model whose state space is BOUNDED, so it can be solved exactly. This
 * method does that, and hands each open station a station-local Markovian
 * environment read off the chain, turning it into a level-dependent QBD whose
 * phase carries the number of closed jobs competing for its server.
 *
 *   1  background chain   the closed population vector over the stations the
 *                         closed classes visit          (mam_bgchain_ctmc)
 *   2  environment        that chain lumped onto the closed occupancy of ONE
 *                         station                       (mam_bgchain_env)
 *   3  open station       a MAP/PH/c queue modulated by that environment,
 *                         solved as a level-dependent QBD (mam_bgchain_station)
 *   4  fixed point        the capacity share feeds step 1 and closes
 *
 * TAGGED-CLASS ITERATION. Step 1 is a population process of dimension (closed
 * chains) x (stations), so its state space is exponential in the number of
 * closed chains R. The method therefore keeps ONE chain free at a time: the
 * tagged chain r is carried exactly, the other R-1 collapse into flow-equivalent
 * aggregate classes whose population is their total and whose service time and
 * routing at each station are their throughput-weighted means
 * (Chandy-Herzog-Woo). Each chain takes its turn as the tagged one and reads its
 * own metrics off the chain it is exact in; the open results are averaged over
 * the passes.
 *
 * HOW MUCH TO AGGREGATE is options.config.bgaggr, the number G of aggregate
 * classes; the background chain then carries 1 + G. G = 1 is the classic
 * tagged/aggregate pair and the default, so the chain stays two-class whatever R
 * is; G >= R-1 aggregates nothing, carries every closed chain exactly, and answers
 * in ONE pass instead of solving the same chain R times. Passing R reaches that, so
 * asking for no aggregation needs no magic value. The cost is the state space, the
 * product over the 1 + G classes of nchoosek(N_b + Mc - 1, Mc - 1), capped by
 * bgstates_max.
 *
 * WHICH CHAINS SHARE A GROUP is decided by similarity of per-station SERVICE
 * DEMAND. An aggregate carries the flow-weighted mean of its members' service times
 * and routing, so it is exact when they place the same demand at every station and
 * distorts in proportion to how far apart they are; grouping the demand-similar
 * chains together keeps the aggregation where it is harmless and away from the
 * chains it would misrepresent.
 *
 * THE EXCHANGED QUANTITY IS THE SHARE, NOT THE MEAN OCCUPANCY. The two halves
 * iterate on `cshare(i,e) = E[min(e+k,c) e/(e+k)]`, the mean number of servers
 * of station i that its e closed jobs hold, averaged over the open occupancy k.
 * Exchanging the mean open occupancy instead and rebuilding the share from it is
 * what a first cut does and it is WRONG: e/(e+k) is convex in k, so Jensen
 * biases the closed service rate down and the closed throughput with it. The
 * environment generator inside the QBD is level-dependent for the same reason.
 *
 * EXACTNESS, measured against SolverCTMC and exact MVA on mixed models of two to
 * four stations: PS or INF with ANY service law (exponential, Erlang, HyperExp,
 * Coxian), any number of servers, Poisson or MAP arrivals and one to four closed
 * chains agree to 4-5 significant digits, as does FCFS with class-INDEPENDENT
 * rates. FCFS with class-DEPENDENT rates keeps the closed queue lengths within
 * ~1% while the open queue length reads 14-20% low, the server being held here
 * in random order rather than head-of-line.
 *
 * PS IS INSENSITIVE to the service law beyond its mean, and the method honours
 * that rather than approximating it: at a PS station the open service is
 * replaced by the exponential of the same mean before the QBD is built. This QBD
 * tracks ONE service phase for the whole station, so carrying the phase-type
 * there makes the open queue length inherit the SCV-sensitivity of an M/PH/1
 * FCFS queue -- measured, a HyperExp of SCV 4 read 21% high where the exact
 * answer is the exponential one to five digits. At an FCFS station the service
 * law IS carried, collapsed into one phase-type process scaled by the share (the
 * same collapse `solver_mam_ldqbd` documents), and the background chain reads
 * only the MEAN closed service time.
 *
 * ARITHMETIC. The path runs `ctmc_solve` on the background chain and the
 * level-dependent QBD recursion on each station, the latter falling back to a
 * pseudo-inverse on a singular level, so it is gated on transcendental
 * arithmetic exactly as `solver_mam_ldqbd` is.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "line/api/mam/ldqbd.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mam/mmap_assemble.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/api/sn/sn_rt_stations.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/state.h"
#include "line/solvers/mam/mam_types.h"
#include "line/solvers/mva/mva_types.h"
#include "line/solvers/mva/sn_chain.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

namespace bgchain_detail {

/**
 * Block of one background class: the ways to place N jobs over the m stations
 * that class visits.
 *
 * This is `qn::space_closed_single`, the lattice primitive the CTMC solver
 * enumerates a closed population over, so the row order and the row count are
 * the reference's rather than this file's; only the support differs, being the
 * class's own stations rather than all of them.
 */
inline std::vector<std::vector<int>> closed_block(std::size_t m, int N) {
    // The primitive is templated on the numeric type of a state row and
    // `num_traits` is specialised for the field types only, so it is taken at
    // double and rounded back, as the JAR twin does off its Matrix.
    const std::vector<std::vector<double>> rows =
        qn::space_closed_single<double>(m, static_cast<std::size_t>(N));
    std::vector<std::vector<int>> out(rows.size(), std::vector<int>(m, 0));
    for (std::size_t r = 0; r < rows.size(); ++r)
        for (std::size_t t = 0; t < m; ++t)
            out[r][t] = static_cast<int>(std::lround(rows[r][t]));
    return out;
}

/**
 * Share of the server capacity that k open jobs hold when e closed jobs are also
 * present: min(k+e,c) busy servers times the open fraction k/(k+e).
 */
inline double open_share(double k, double e, double c) {
    const double tot = k + e;
    if (!(tot > 0.0)) return 0.0;
    return std::min(tot, c) * (k / tot);
}

/** The mirror image of open_share, for the closed jobs. */
inline double closed_share(double e, double k, double c) {
    const double tot = e + k;
    if (!(tot > 0.0)) return 0.0;
    return std::min(tot, c) * (e / tot);
}

/**
 * Group the columns of D into G clusters by per-station SERVICE DEMAND
 * (mam_bgchain_groups.m).
 *
 * WHY DEMAND IS THE RIGHT CRITERION. The aggregate that replaces a group carries
 * the flow-weighted mean of its members' service times and routing, so the group
 * aggregates EXACTLY when its members place the same demand at every station and
 * distorts both quantities in proportion to how far apart they are. The distance
 * is the symmetric relative L1 gap between the demand vectors,
 *
 *   dist(a,b) = sum_i |D[i][a] - D[i][b]| / ((sum_i D[i][a] + sum_i D[i][b])/2),
 *
 * scale-relative rather than absolute: it separates two chains whose demand
 * PROFILE across the stations differs and two whose profile agrees but whose
 * magnitude does not, and being dimensionless it groups a model the same way
 * whatever its time unit.
 *
 * WHY COMPLETE LINKAGE. Agglomerative from singletons, merging the pair of
 * clusters whose WORST member-to-member distance is smallest. The aggregation
 * error inside a group is driven by its worst mismatch, not its average one.
 *
 * DETERMINISM. Ties break on the lexicographically smallest pair of cluster
 * indices and the groups are relabelled by their smallest member, so the same
 * input gives the same grouping in MATLAB, the JAR, Python and C++.
 *
 * @param D (Mc x n) per-station demand, one column per chain
 * @param G number of groups wanted, clamped to [1, n]
 * @return group index in 0..G-1 of each chain
 */
inline std::vector<std::size_t> mam_bgchain_groups(const std::vector<std::vector<double>>& D,
                                                   std::size_t G) {
    const std::size_t Mc = D.size();
    const std::size_t n = (Mc == 0) ? 0 : D[0].size();
    std::vector<std::size_t> grp(n, 0);
    if (n == 0) return grp;
    if (G < 1) G = 1;
    if (G > n) G = n;

    std::vector<double> tot(n, 0.0);
    for (std::size_t c = 0; c < n; ++c)
        for (std::size_t i = 0; i < Mc; ++i) tot[c] += D[i][c];
    std::vector<std::vector<double>> dist(n, std::vector<double>(n, 0.0));
    for (std::size_t a = 0; a < n; ++a) {
        for (std::size_t b = a + 1; b < n; ++b) {
            const double den = (tot[a] + tot[b]) / 2.0;
            double d = 0.0;
            if (den > 1e-14) {
                double acc = 0.0;
                for (std::size_t i = 0; i < Mc; ++i) acc += std::fabs(D[i][a] - D[i][b]);
                d = acc / den;
            }
            dist[a][b] = d;
            dist[b][a] = d;
        }
    }

    std::vector<std::vector<std::size_t>> clusters(n);
    for (std::size_t a = 0; a < n; ++a) clusters[a].push_back(a);
    std::vector<bool> active(n, true);
    std::size_t nactive = n;
    while (nactive > G) {
        double best = std::numeric_limits<double>::infinity();
        long bp = -1, bq = -1;
        for (std::size_t p = 0; p < n; ++p) {
            if (!active[p]) continue;
            for (std::size_t q = p + 1; q < n; ++q) {
                if (!active[q]) continue;
                double d = 0.0;
                for (std::size_t x : clusters[p])
                    for (std::size_t y : clusters[q])
                        if (dist[x][y] > d) d = dist[x][y];
                if (d < best - 1e-14) {
                    best = d;
                    bp = static_cast<long>(p);
                    bq = static_cast<long>(q);
                }
            }
        }
        if (bp < 0) break;
        std::vector<std::size_t>& cp = clusters[static_cast<std::size_t>(bp)];
        std::vector<std::size_t>& cq = clusters[static_cast<std::size_t>(bq)];
        cp.insert(cp.end(), cq.begin(), cq.end());
        std::sort(cp.begin(), cp.end());
        cq.clear();
        active[static_cast<std::size_t>(bq)] = false;
        --nactive;
    }

    // Relabel by smallest member, so the group numbering is canonical
    std::vector<std::size_t> live;
    for (std::size_t p = 0; p < n; ++p)
        if (active[p]) live.push_back(p);
    std::sort(live.begin(), live.end(), [&](std::size_t a, std::size_t b) {
        return clusters[a][0] < clusters[b][0];
    });
    for (std::size_t g = 0; g < live.size(); ++g)
        for (std::size_t x : clusters[live[g]]) grp[x] = g;
    return grp;
}

}  // namespace bgchain_detail

/** The solved background modulating chain. */
template <class T>
struct BgchainCtmc {
    /** space[s][i][b]: class-b jobs held by station i in state s. */
    std::vector<std::vector<std::vector<int>>> space;
    /** totocc[s][i]: closed jobs of any background class held by station i. */
    std::vector<std::vector<int>> totocc;
    std::vector<T> pi;
    Matrix<T> Q;
    std::vector<std::vector<T>> QLen, Tput, Ubusy;  ///< (Mc x B)
    std::size_t nstates = 0;
};

/**
 * Build and solve the background chain (mam_bgchain_ctmc.m).
 *
 * A station holding e closed jobs serves background class b at rate
 * `n[i][b] / STb[i][b]` when it is an infinite server, and
 * `cshare[i][e] * (n[i][b]/e) / STb[i][b]` otherwise, splitting the capacity the
 * closed jobs hold over the background classes in proportion to their counts.
 * That is exact under PS and is the random-order surrogate under FCFS.
 *
 * @param Nb      population of each background class, size B in {1,2}
 * @param STb     (Mc x B) mean service time per station per background class
 * @param Pb      B row-stochastic (Mc x Mc) routing matrices
 * @param isinf_i whether each station of the support is an infinite server
 * @param nsrv    servers of each station of the support
 * @param cshare  (Mc x (Nmax+1)) mean servers the e closed jobs hold
 * @param supp    supp[i][b]: station i is on the route of background class b.
 *                NOT an optimization -- a chain that never visits a station
 *                cannot hold jobs there, and enumerating the union of every
 *                chain's stations puts probability on unreachable configurations
 *                that also ABSORB, because the chain's routing matrix has a zero
 *                row at an unvisited station which row-normalizes to a self-loop.
 *                The generator turns reducible and population conservation
 *                silently fails
 * @param states_max cap on the number of chain states
 */
template <class T>
BgchainCtmc<T> mam_bgchain_ctmc(const std::vector<int>& Nb,
                                const std::vector<std::vector<T>>& STb,
                                const std::vector<Matrix<T>>& Pb,
                                const std::vector<bool>& isinf_i,
                                const std::vector<double>& nsrv,
                                const std::vector<std::vector<double>>& cshare,
                                const std::vector<std::vector<bool>>& supp,
                                std::size_t states_max) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t Mc = STb.size();
    const std::size_t B = Nb.size();

    // Each class is enumerated over ITS OWN stations only; see the supp param.
    std::vector<std::vector<std::vector<int>>> sp(B);      // expanded to Mc columns
    std::vector<std::vector<std::vector<int>>> compb(B);   // reduced to the class's stations
    std::vector<std::vector<std::size_t>> idxb(B);
    std::vector<std::size_t> nst(B);
    std::size_t nstates = 1;
    for (std::size_t b = 0; b < B; ++b) {
        for (std::size_t i = 0; i < Mc; ++i)
            if (supp.empty() || supp[i][b]) idxb[b].push_back(i);
        if (idxb[b].empty()) {
            if (Nb[b] > 0)
                throw UnsupportedError("mam_bgchain_ctmc: background class " + std::to_string(b + 1) +
                                       " holds " + std::to_string(Nb[b]) +
                                       " jobs but visits no station");
            idxb[b].push_back(0);  // an empty class still needs one slot
        }
        compb[b] = bgchain_detail::closed_block(idxb[b].size(), Nb[b]);
        sp[b].assign(compb[b].size(), std::vector<int>(Mc, 0));
        for (std::size_t r = 0; r < compb[b].size(); ++r)
            for (std::size_t t = 0; t < idxb[b].size(); ++t)
                sp[b][r][idxb[b][t]] = compb[b][r][t];
        nst[b] = compb[b].size();
        nstates *= nst[b];
    }
    if (nstates > states_max)
        throw UnsupportedError(
            "mam_bgchain_ctmc: the background chain of this model has " + std::to_string(nstates) +
            " states, above the limit of " + std::to_string(states_max) +
            ". The chain enumerates the closed-class population vector over the " +
            std::to_string(Mc) +
            " stations the closed classes visit, so its size grows as nchoosek(N+Mc-1,Mc-1) per "
            "class, and it carries " + std::to_string(B) +
            " classes. Lower options.config.bgaggr to aggregate more of the closed chains into "
            "fewer classes, raise options.config.bgstates_max to solve it anyway, or reduce the "
            "closed populations.");

    // Per-class transition targets: tgt[b][s][i*Mc+j] is the class-b state
    // reached from s when one job moves from station i to station j, or -1.
    std::vector<std::vector<std::vector<long>>> tgt(B);
    for (std::size_t b = 0; b < B; ++b) {
        const std::size_t mb = idxb[b].size();
        std::map<std::vector<int>, long> index;
        for (std::size_t s = 0; s < nst[b]; ++s) index[compb[b][s]] = static_cast<long>(s);
        tgt[b].assign(nst[b], std::vector<long>(Mc * Mc, -1));
        for (std::size_t s = 0; s < nst[b]; ++s) {
            for (std::size_t ii = 0; ii < mb; ++ii) {
                if (compb[b][s][ii] == 0) continue;
                for (std::size_t jj = 0; jj < mb; ++jj) {
                    if (jj == ii) continue;
                    std::vector<int> cand = compb[b][s];
                    --cand[ii];
                    ++cand[jj];
                    const auto it = index.find(cand);
                    if (it != index.end())
                        tgt[b][s][idxb[b][ii] * Mc + idxb[b][jj]] = it->second;
                }
            }
        }
    }

    // Joint space, class 0 outermost.
    std::vector<std::size_t> strideb(B, 1);
    for (std::size_t b = 0; b < B; ++b) {
        std::size_t s = 1;
        for (std::size_t b2 = b + 1; b2 < B; ++b2) s *= nst[b2];
        strideb[b] = s;
    }
    BgchainCtmc<T> out;
    out.nstates = nstates;
    out.space.assign(nstates, std::vector<std::vector<int>>(Mc, std::vector<int>(B, 0)));
    out.totocc.assign(nstates, std::vector<int>(Mc, 0));
    std::vector<std::vector<std::size_t>> subidx(nstates, std::vector<std::size_t>(B, 0));
    for (std::size_t s = 0; s < nstates; ++s) {
        std::size_t rem = s;
        for (std::size_t bb = B; bb-- > 0;) {
            subidx[s][bb] = rem % nst[bb];
            rem /= nst[bb];
        }
        for (std::size_t b = 0; b < B; ++b) {
            const std::vector<int>& vec = sp[b][subidx[s][b]];
            for (std::size_t i = 0; i < Mc; ++i) {
                out.space[s][i][b] = vec[i];
                out.totocc[s][i] += vec[i];
            }
        }
    }

    std::vector<std::vector<T>> mu(Mc, std::vector<T>(B, zero));
    for (std::size_t b = 0; b < B; ++b)
        for (std::size_t i = 0; i < Mc; ++i)
            if (num_traits<T>::to_double(STb[i][b]) > 0.0)
                mu[i][b] = T(num_traits<T>::from_int(1) / STb[i][b]);

    Matrix<T> Q(nstates, nstates, zero);
    std::vector<std::vector<std::vector<T>>> rate_full(
        nstates, std::vector<std::vector<T>>(Mc, std::vector<T>(B, zero)));
    std::vector<std::vector<T>> cap_busy(nstates, std::vector<T>(Mc, zero));
    for (std::size_t s = 0; s < nstates; ++s) {
        for (std::size_t i = 0; i < Mc; ++i) {
            const int eclosed = out.totocc[s][i];
            if (eclosed == 0) continue;
            double held_d;
            if (isinf_i[i]) {
                held_d = static_cast<double>(eclosed);
            } else {
                const std::size_t ecap =
                    std::min<std::size_t>(static_cast<std::size_t>(eclosed), cshare[i].size() - 1);
                held_d = cshare[i][ecap];
            }
            if (!(held_d > 0.0)) continue;
            const T held = num_traits<T>::from_double(held_d);
            cap_busy[s][i] = held;
            for (std::size_t b = 0; b < B; ++b) {
                if (out.space[s][i][b] == 0 || mu[i][b] == zero) continue;
                const T frac = num_traits<T>::from_double(static_cast<double>(out.space[s][i][b]) /
                                                          static_cast<double>(eclosed));
                const T r = T(held * frac * mu[i][b]);
                rate_full[s][i][b] = r;
                for (std::size_t j = 0; j < Mc; ++j) {
                    if (j == i) continue;
                    if (!(num_traits<T>::to_double(Pb[b](i, j)) > 0.0)) continue;
                    const long tsub = tgt[b][subidx[s][b]][i * Mc + j];
                    if (tsub < 0) continue;
                    const std::size_t sdest =
                        s + (static_cast<std::size_t>(tsub) - subidx[s][b]) * strideb[b];
                    Q(s, sdest) += T(r * Pb[b](i, j));
                }
            }
        }
    }

    out.Q = mc::ctmc_makeinfgen(Q);
    if (nstates == 1) {
        out.pi.assign(1, num_traits<T>::from_int(1));
    } else {
        out.pi = mc::ctmc_solve(out.Q);
    }
    T tot = zero;
    for (T& v : out.pi) {
        if (num_traits<T>::to_double(v) < 0.0) v = zero;
        tot += v;
    }
    if (tot != zero)
        for (T& v : out.pi) v /= tot;

    out.QLen.assign(Mc, std::vector<T>(B, zero));
    out.Tput.assign(Mc, std::vector<T>(B, zero));
    out.Ubusy.assign(Mc, std::vector<T>(B, zero));
    for (std::size_t s = 0; s < nstates; ++s) {
        const T p = out.pi[s];
        if (p == zero) continue;
        for (std::size_t i = 0; i < Mc; ++i)
            for (std::size_t b = 0; b < B; ++b) {
                out.QLen[i][b] += T(p * num_traits<T>::from_int(out.space[s][i][b]));
                out.Tput[i][b] += T(p * rate_full[s][i][b]);
            }
    }
    for (std::size_t i = 0; i < Mc; ++i) {
        if (isinf_i[i]) {
            out.Ubusy[i] = out.QLen[i];
            continue;
        }
        for (std::size_t s = 0; s < nstates; ++s) {
            const T p = out.pi[s];
            const int occ = out.totocc[s][i];
            if (p == zero || occ == 0) continue;
            for (std::size_t b = 0; b < B; ++b) {
                const T share = num_traits<T>::from_double(
                    static_cast<double>(out.space[s][i][b]) / static_cast<double>(occ));
                out.Ubusy[i][b] +=
                    T(p * cap_busy[s][i] * share / num_traits<T>::from_double(nsrv[i]));
            }
        }
    }
    return out;
}

/** The environment one station sees: the background chain lumped onto its occupancy. */
template <class T>
struct BgchainEnv {
    Matrix<T> A;
    std::vector<T> phi;
    std::vector<int> esup;  ///< ascending
};

/**
 * Lump the background chain onto the closed occupancy of station i
 * (mam_bgchain_env.m).
 *
 * Station i does not observe the whole closed population vector, only how many
 * closed jobs compete with the open ones for its server. The lumped generator is
 * the stationary-weighted aggregation of Q over the level sets
 * {s : totocc(s,i) = e}, exact when the partition is lumpable in the
 * Kemeny-Snell sense and the standard exact-aggregation approximation otherwise.
 * The diagonal is rebuilt from the off-diagonal row sums, so the result is a
 * proper generator whatever the lumping error is. Environment states of zero
 * stationary probability are unreachable and are dropped.
 */
template <class T>
BgchainEnv<T> mam_bgchain_env(const BgchainCtmc<T>& bg, std::size_t i) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t n = bg.nstates;

    std::vector<int> levels;
    for (std::size_t s = 0; s < n; ++s) levels.push_back(bg.totocc[s][i]);
    std::sort(levels.begin(), levels.end());
    levels.erase(std::unique(levels.begin(), levels.end()), levels.end());

    std::vector<T> wAll(levels.size(), zero);
    for (std::size_t s = 0; s < n; ++s) {
        const std::size_t pos = static_cast<std::size_t>(
            std::lower_bound(levels.begin(), levels.end(), bg.totocc[s][i]) - levels.begin());
        wAll[pos] += bg.pi[s];
    }

    BgchainEnv<T> env;
    std::vector<T> w;
    for (std::size_t e = 0; e < levels.size(); ++e) {
        if (num_traits<T>::to_double(wAll[e]) > 1e-14) {
            env.esup.push_back(levels[e]);
            w.push_back(wAll[e]);
        }
    }
    if (env.esup.empty()) {
        // degenerate chain: the station never holds a closed job
        env.A = Matrix<T>(1, 1, zero);
        env.phi.assign(1, num_traits<T>::from_int(1));
        env.esup.assign(1, 0);
        return env;
    }

    const std::size_t me = env.esup.size();
    std::vector<long> lvl(n, -1);
    for (std::size_t s = 0; s < n; ++s) {
        const auto it = std::lower_bound(env.esup.begin(), env.esup.end(), bg.totocc[s][i]);
        if (it != env.esup.end() && *it == bg.totocc[s][i])
            lvl[s] = static_cast<long>(it - env.esup.begin());
    }

    env.A = Matrix<T>(me, me, zero);
    if (me > 1) {
        for (std::size_t s = 0; s < n; ++s) {
            if (lvl[s] < 0 || bg.pi[s] == zero) continue;
            for (std::size_t sp = 0; sp < n; ++sp) {
                if (sp == s || lvl[sp] < 0 || lvl[sp] == lvl[s]) continue;
                const T q = bg.Q(s, sp);
                if (q == zero) continue;
                env.A(static_cast<std::size_t>(lvl[s]), static_cast<std::size_t>(lvl[sp])) +=
                    T(bg.pi[s] * q);
            }
        }
        for (std::size_t e = 0; e < me; ++e) {
            T diag = zero;
            for (std::size_t ep = 0; ep < me; ++ep) {
                if (ep == e) continue;
                env.A(e, ep) /= w[e];
                diag += env.A(e, ep);
            }
            env.A(e, e) = T(-diag);
        }
    }
    T wsum = zero;
    for (const T& v : w) wsum += v;
    env.phi.assign(me, zero);
    for (std::size_t e = 0; e < me; ++e) env.phi[e] = T(w[e] / wsum);
    return env;
}

/** What one modulated station QBD returns. */
template <class T>
struct BgchainStation {
    T QLen, Util, Tput, ploss;
    std::vector<T> penv;
    std::vector<double> cshare;  ///< E[min(e+k,c) e/(e+k) | e], indexed by esup
    std::vector<int> esup;
};

/**
 * Solve the open classes of one station as a modulated level-dependent QBD
 * (mam_bgchain_station.m).
 *
 * Level = number of open jobs held by the station, phase = (arrival MAP phase,
 * environment state, service phase). With k open and e closed jobs present the
 * open aggregate completes at rate `min(k+e,c) k/(k+e)` times the phase-type
 * completion rate of one busy server: the dependence on k makes the QBD
 * level-dependent, the dependence on e makes it modulated.
 *
 * The environment is level-dependent too. A lumped transition that LOWERS the
 * closed occupancy is a closed completion here, so it carries the closed share
 * and level k rescales it by the ratio to the averaged share `gref` the chain
 * was built at; a transition that RAISES the occupancy is an arrival from
 * elsewhere and is left alone. Without this the closed jobs would drain at their
 * mean-field rate however long the open queue is, and the positive correlation
 * between the two occupancies would be lost.
 *
 * The level space is truncated at Kmax. An arrival at the top level is lost but
 * still advances the arrival phase, so the arrival process keeps its exact
 * marginal and autocorrelation and only the queue tail is cut.
 */
template <class T>
BgchainStation<T> mam_bgchain_station(const Matrix<T>& Da0, const Matrix<T>& Da1,
                                      const std::vector<T>& alpha_s, const Matrix<T>& Tsvc,
                                      const Matrix<T>& Ain, const std::vector<int>& esup_in,
                                      double nservers, const std::vector<double>& gref_in,
                                      std::size_t Kmax) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t ma = Da0.rows();
    const std::size_t me = esup_in.size();
    const std::size_t ms = alpha_s.size();
    if (Kmax < 1) Kmax = 1;

    // esup arrives ascending from mam_bgchain_env; gref is indexed to match.
    const std::vector<int>& esup = esup_in;
    const std::vector<double>& gref = gref_in;

    Matrix<T> tvec(ms, 1, zero);
    for (std::size_t i = 0; i < ms; ++i) {
        T s = zero;
        for (std::size_t j = 0; j < ms; ++j) s += Tsvc(i, j);
        tvec(i, 0) = T(-s);
    }
    Matrix<T> alphaRow(1, ms, zero);
    for (std::size_t j = 0; j < ms; ++j) alphaRow(0, j) = alpha_s[j];

    Matrix<T> Ime(me, me, zero), Ima(ma, ma, zero), Ims(ms, ms, zero);
    for (std::size_t i = 0; i < me; ++i) Ime(i, i) = one;
    for (std::size_t i = 0; i < ma; ++i) Ima(i, i) = one;
    for (std::size_t i = 0; i < ms; ++i) Ims(i, i) = one;

    // Split the environment into the closed departures from this station (which
    // the open level throttles) and the arrivals to it (which it does not).
    Matrix<T> Adown(me, me, zero), Aup(me, me, zero);
    for (std::size_t e = 0; e < me; ++e)
        for (std::size_t ep = 0; ep < me; ++ep) {
            if (ep < e) Adown(e, ep) = Ain(e, ep);
            else if (ep > e) Aup(e, ep) = Ain(e, ep);
        }

    auto env_at_level = [&](std::size_t k) {
        Matrix<T> Ak(me, me, zero);
        for (std::size_t e = 0; e < me; ++e) {
            const double g = bgchain_detail::closed_share(static_cast<double>(esup[e]),
                                                          static_cast<double>(k), nservers);
            const double ratio = (gref[e] > 0.0) ? g / gref[e] : 1.0;
            const T rt = num_traits<T>::from_double(ratio);
            T diag = zero;
            for (std::size_t ep = 0; ep < me; ++ep) {
                if (ep == e) continue;
                const T v = T(Aup(e, ep) + rt * Adown(e, ep));
                Ak(e, ep) = v;
                diag += v;
            }
            Ak(e, e) = T(-diag);
        }
        return Ak;
    };
    auto diagm = [&](const std::vector<T>& v) {
        Matrix<T> D(v.size(), v.size(), zero);
        for (std::size_t i = 0; i < v.size(); ++i) D(i, i) = v[i];
        return D;
    };

    // The C++ ldqbd indexes q0/q1/q2 BY LEVEL, with q2[0] unused, unlike the
    // reference's 1-based cell arrays.
    std::vector<Matrix<T>> Q0(Kmax), Q1(Kmax + 1), Q2(Kmax + 1);
    std::vector<std::vector<T>> phiae(Kmax);

    Q1[0] = qbd_detail::madd(kron(Da0, Ime), kron(Ima, env_at_level(0)));
    Q0[0] = kron(kron(Da1, Ime), alphaRow);

    const Matrix<T> Da0kron = kron(kron(Da0, Ime), Ims);
    const Matrix<T> Da1kron = kron(kron(Da1, Ime), Ims);
    Q2[0] = Matrix<T>(1, 1, zero);  // unused
    for (std::size_t k = 1; k <= Kmax; ++k) {
        std::vector<T> rep(ma * me, zero);
        for (std::size_t a = 0; a < ma; ++a)
            for (std::size_t e = 0; e < me; ++e)
                rep[a * me + e] = num_traits<T>::from_double(bgchain_detail::open_share(
                    static_cast<double>(k), static_cast<double>(esup[e]), nservers));
        phiae[k - 1] = rep;
        const Matrix<T> Ak = env_at_level(k);
        Q1[k] = qbd_detail::madd(
            qbd_detail::madd(Da0kron, kron(kron(Ima, Ak), Ims)), kron(diagm(rep), Tsvc));
        if (k < Kmax) Q0[k] = Da1kron;
        if (k == 1) {
            Q2[1] = kron(diagm(rep), tvec);
        } else {
            Q2[k] = kron(diagm(rep), matmul(tvec, alphaRow));
        }
    }
    // truncation: an arrival at the top level is lost, its phase transition is kept
    Q1[Kmax] = qbd_detail::madd(Q1[Kmax], Da1kron);

    const LdqbdResult<T> res = ldqbd(Q0, Q1, Q2);
    const std::vector<T>& plev = res.pi.pi;

    BgchainStation<T> out;
    out.esup = esup;
    out.QLen = zero;
    for (std::size_t k = 0; k <= Kmax; ++k)
        out.QLen += T(num_traits<T>::from_int(static_cast<int>(k)) * plev[k]);
    out.ploss = plev[Kmax];

    std::vector<T> penv(me, zero), gacc(me, zero);
    T util = zero, tput = zero;
    for (std::size_t k = 0; k <= Kmax; ++k) {
        const std::vector<T>& pk = res.pi.pi_level[k];
        std::vector<T> marg(me, zero);
        if (k == 0) {
            for (std::size_t a = 0; a < ma; ++a)
                for (std::size_t e = 0; e < me; ++e) marg[e] += pk[a * me + e];
        } else {
            for (std::size_t a = 0; a < ma; ++a)
                for (std::size_t e = 0; e < me; ++e) {
                    T block = zero, dep = zero;
                    for (std::size_t s = 0; s < ms; ++s) {
                        const T v = pk[(a * me + e) * ms + s];
                        block += v;
                        dep += T(v * tvec(s, 0));
                    }
                    marg[e] += block;
                    util += T(block * phiae[k - 1][a * me + e]);
                    tput += T(dep * phiae[k - 1][a * me + e]);
                }
        }
        for (std::size_t e = 0; e < me; ++e) {
            penv[e] += marg[e];
            gacc[e] += T(marg[e] * num_traits<T>::from_double(bgchain_detail::closed_share(
                                       static_cast<double>(esup[e]), static_cast<double>(k),
                                       nservers)));
        }
    }
    T psum = zero;
    for (const T& v : penv) psum += v;
    if (psum != zero)
        for (std::size_t e = 0; e < me; ++e) {
            penv[e] /= psum;
            gacc[e] /= psum;
        }
    out.penv = penv;
    out.cshare.assign(me, 0.0);
    for (std::size_t e = 0; e < me; ++e) {
        const double pe = num_traits<T>::to_double(penv[e]);
        out.cshare[e] = (pe > 1e-14) ? num_traits<T>::to_double(gacc[e]) / pe : 0.0;
    }
    out.Util = T(util / num_traits<T>::from_double(nservers));
    out.Tput = tput;
    return out;
}

/**
 * Port of `solver_mam_bgchain.m`.
 *
 * @param L   the refreshed struct; must be mixed (at least one open and one
 *            closed chain)
 * @param opt the MAM options; `cutoff` bounds the open level truncation,
 *            `bgstates_max` the background chain, `qbdphases_max` each station
 */
/**
 * Number of states of the background-chain CTMC solver_mam_bgchain would build
 * on this model, WITHOUT building it. Mirrors `mam_bgchain_states.m`.
 *
 * The size is what decides whether bgchain is affordable and mam_bgchain_ctmc
 * only discovers it after the partition is fixed, so the default-method chooser
 * needs it up front. The count follows the partition below: a pass carries the
 * tagged closed chain as background class 0 and the demand-similar groups of the
 * other closed chains as classes 1..G, each enumerating the compositions of its
 * population over the stations its members visit. Merging two chains onto the
 * UNION of their supports can raise the count as easily as lower it, so the
 * passes are enumerated rather than bounded and the largest returned: that is
 * the one mam_bgchain_ctmc would refuse. Returns 0 when bgchain does not apply
 * to the model at all.
 */
template <class T>
double bgchain_states(const qn::NetworkStruct<T>& L, const MamOptions& opt) {
    const std::size_t M = L.nstations, C = L.nchains;
    const mva::ChainDemands<T> dem = mva::sn_get_demands_chain(L);

    std::vector<std::size_t> closedChains;
    for (std::size_t c = 0; c < C; ++c) {
        bool open = false;
        for (std::size_t k : L.inchain[c])
            if (std::isinf(L.classes[k - 1].population)) open = true;
        if (!open && dem.Nchain[c] > 0.0) closedChains.push_back(c);
    }
    const std::size_t R = closedChains.size();
    if (R == 0) return 0.0;

    std::vector<std::size_t> cst;
    for (std::size_t i = 0; i < M; ++i) {
        bool visited = false;
        for (std::size_t c : closedChains)
            if (num_traits<T>::to_double(dem.Vchain(i, c)) > 1e-14) visited = true;
        if (visited) cst.push_back(i);
    }
    const std::size_t Mc = cst.size();
    if (Mc == 0) return 0.0;

    const std::size_t bgaggr_opt = (opt.bgaggr > 0) ? opt.bgaggr : 1;
    const std::size_t naggr =
        std::min(std::max<std::size_t>(bgaggr_opt, 1), std::max<std::size_t>(R - 1, 1));
    const bool noAggr = (R == 1) || (naggr >= R - 1);
    const std::size_t npass = noAggr ? 1 : R;

    // nchoosek in floating point, so a chain far above any usable size still compares
    auto binomial = [](double n, double k) {
        if (k < 0.0 || k > n) return 0.0;
        const double kk = std::min(k, n - k);
        double acc = 1.0;
        for (double i = 1.0; i <= kk; i += 1.0) acc = acc * (n - kk + i) / i;
        return acc;
    };

    double worst = 0.0;
    for (std::size_t pidx = 0; pidx < npass; ++pidx) {
        std::vector<std::vector<std::size_t>> members;
        if (noAggr) {
            for (std::size_t c : closedChains) members.push_back({c});
        } else {
            std::vector<std::size_t> others;
            for (std::size_t oi = 0; oi < R; ++oi)
                if (oi != pidx) others.push_back(closedChains[oi]);
            std::vector<std::vector<double>> D(Mc, std::vector<double>(others.size(), 0.0));
            for (std::size_t ii = 0; ii < Mc; ++ii)
                for (std::size_t oi = 0; oi < others.size(); ++oi)
                    D[ii][oi] = num_traits<T>::to_double(dem.Lchain(cst[ii], others[oi]));
            const std::vector<std::size_t> grp = bgchain_detail::mam_bgchain_groups(D, naggr);
            members.push_back({closedChains[pidx]});
            for (std::size_t g = 0; g < naggr; ++g) {
                std::vector<std::size_t> mem;
                for (std::size_t oi = 0; oi < others.size(); ++oi)
                    if (grp[oi] == g) mem.push_back(others[oi]);
                members.push_back(mem);
            }
        }

        double n = 1.0;
        for (const std::vector<std::size_t>& mem : members) {
            double Nb = 0.0;
            for (std::size_t o : mem) Nb += std::llround(dem.Nchain[o]);
            std::size_t m = 0;
            for (std::size_t ii = 0; ii < Mc; ++ii)
                for (std::size_t o : mem)
                    if (num_traits<T>::to_double(dem.Vchain(cst[ii], o)) > 1e-14) {
                        ++m;
                        break;
                    }
            if (m == 0) m = 1;   // an empty class still needs one slot to be indexed by
            n *= binomial(Nb + static_cast<double>(m) - 1.0, static_cast<double>(m) - 1.0);
            if (!std::isfinite(n)) return std::numeric_limits<double>::infinity();
        }
        worst = std::max(worst, n);
    }
    return worst;
}

template <class T>
mva::MvaSolution<T> solver_mam_bgchain(const qn::NetworkStruct<T>& L, const MamOptions& opt) {
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "solver_mam_bgchain: the level-dependent QBD recursion inverts a matrix per level and "
            "falls back to a pseudo-inverse when a level is singular, neither of which is exact "
            "arithmetic; rerun with --arith double or --arith real");
    } else {
    using lang::SchedStrategy;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations, K = L.nclasses, C = L.nchains;

    // ---- class-level service times ---------------------------------------
    Matrix<T> S(M, K, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k) {
            const double r = num_traits<T>::to_double(L.rates(i, k));
            if (std::isfinite(r) && r > 0.0) S(i, k) = T(one / L.rates(i, k));
        }

    const api::SnRtStations<T> rtv = api::sn_rt_stations(L);
    const Matrix<T>& rtst = rtv.rtst;
    const Matrix<T>& V = rtv.Vst;
    const mva::ChainDemands<T> dem = mva::sn_get_demands_chain(L);

    std::vector<bool> isopenchain(C, false);
    std::vector<std::size_t> openChains, closedChains;
    for (std::size_t c = 0; c < C; ++c) {
        bool open = false;
        for (std::size_t k : L.inchain[c])
            if (std::isinf(L.classes[k - 1].population)) open = true;
        isopenchain[c] = open;
        if (open) openChains.push_back(c);
        else if (dem.Nchain[c] > 0.0) closedChains.push_back(c);
    }
    const std::size_t R = closedChains.size();
    if (R == 0)
        throw UnsupportedError(
            "solver_mam_bgchain: the bgchain method requires at least one closed class: the "
            "background chain IS the closed population vector, so a purely open model has nothing "
            "to build it from. Use dec.source.");

    // ---- the support of the background chain -----------------------------
    std::vector<std::size_t> cst;
    for (std::size_t i = 0; i < M; ++i) {
        bool visited = false;
        for (std::size_t c : closedChains)
            if (num_traits<T>::to_double(dem.Vchain(i, c)) > 1e-14) visited = true;
        if (visited) cst.push_back(i);
    }
    const std::size_t Mc = cst.size();
    if (Mc == 0)
        throw UnsupportedError("solver_mam_bgchain: the closed classes of this model visit no "
                               "station");

    // ---- chain-level station routing, folding the class axis of rt --------
    std::vector<Matrix<T>> Pchain(C, Matrix<T>(M, M, zero));
    for (std::size_t c = 0; c < C; ++c) {
        Matrix<T> P(M, M, zero);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t k1 : L.inchain[c]) {
                const T a = dem.alpha(i, k1 - 1);
                if (!(num_traits<T>::to_double(a) > 0.0)) continue;
                for (std::size_t j = 0; j < M; ++j) {
                    T acc = zero;
                    for (std::size_t k2 : L.inchain[c])
                        acc += rtst(i * K + (k1 - 1), j * K + (k2 - 1));
                    if (acc != zero) P(i, j) += T(a * acc);
                }
            }
        Pchain[c] = P;
    }

    // ---- open arrival streams --------------------------------------------
    std::vector<double> lambdaChain(C, 0.0);
    std::map<std::size_t, Map<T>> chainArrival;
    for (std::size_t c : openChains) {
        const std::size_t isrc = L.classes[L.inchain[c][0] - 1].refstat;  // 1-based
        double lam = 0.0;
        std::vector<Mmap<T>> parts;
        for (std::size_t k : L.inchain[c]) {
            const double rk = num_traits<T>::to_double(L.rates(isrc - 1, k - 1));
            if (!std::isfinite(rk) || !(rk > 0.0)) continue;
            lam += rk;
            const Map<T> mk = lang::dist_to_map(L.service[isrc - 1][k - 1]);
            Mmap<T> mm;
            mm.D0 = mk.D0;
            mm.D1 = mk.D1;
            mm.Dc.assign(1, mk.D1);
            parts.push_back(mm);
        }
        lambdaChain[c] = lam;
        if (!parts.empty()) {
            const Mmap<T> sup =
                (parts.size() == 1) ? parts[0] : mmap_super_safe(parts, opt.space_max);
            chainArrival[c] = Map<T>{sup.D0, sup.D1};
        }
    }

    Matrix<T> lambdaOpen(M, K, zero);
    std::vector<bool> isopenclass(K, false);
    for (std::size_t c : openChains)
        for (std::size_t k : L.inchain[c]) {
            isopenclass[k - 1] = true;
            for (std::size_t i = 0; i < M; ++i)
                lambdaOpen(i, k - 1) = T(num_traits<T>::from_double(lambdaChain[c]) * V(i, k - 1));
        }

    // ---- fixed-point state -----------------------------------------------
    mva::MvaSolution<T> sol;
    sol.Q = Matrix<T>(M, K, zero);
    sol.U = Matrix<T>(M, K, zero);
    sol.R = Matrix<T>(M, K, zero);
    sol.Tp = Matrix<T>(M, K, zero);
    sol.C.assign(K, zero);
    sol.X.assign(K, zero);

    int Ntot = 0;
    for (std::size_t c : closedChains) Ntot += static_cast<int>(std::llround(dem.Nchain[c]));

    // cshare[i][e]: mean number of servers of station i that its e closed jobs
    // hold once the open work has taken its share. Starts at min(e,c).
    std::vector<std::vector<double>> cshare(M, std::vector<double>(Ntot + 1, 0.0));
    std::vector<double> nsrvAll(M, 1.0);
    std::vector<bool> isinfAll(M, false);
    for (std::size_t i = 0; i < M; ++i) {
        nsrvAll[i] = L.stations[i].nservers;
        isinfAll[i] = (L.stations[i].sched == SchedStrategy::INF);
        for (int e = 0; e <= Ntot; ++e)
            cshare[i][e] = std::min(static_cast<double>(e), nsrvAll[i]);
    }
    std::vector<double> Xclosed(C, 0.0);
    for (std::size_t c : closedChains) {
        double denom = 0.0;
        for (std::size_t i = 0; i < M; ++i) denom += num_traits<T>::to_double(dem.Lchain(i, c));
        if (denom > 0.0) Xclosed[c] = dem.Nchain[c] / denom;
    }

    const std::size_t states_max = (opt.bgstates_max > 0) ? opt.bgstates_max : 20000;
    const std::size_t phases_max = (opt.qbdphases_max > 0) ? opt.qbdphases_max : 500;

    // How many aggregate classes the background chain carries, and which chains
    // share each of them. opt.bgaggr is the number of AGGREGATE classes G: G = 1
    // is the classic tagged/aggregate pair, G >= R-1 aggregates nothing.
    const std::size_t bgaggr_opt = (opt.bgaggr > 0) ? opt.bgaggr : 1;
    const std::size_t naggr =
        std::min(std::max<std::size_t>(bgaggr_opt, 1), std::max<std::size_t>(R - 1, 1));
    // With nothing left to aggregate ONE background chain carries every closed
    // chain exactly, so the tagged loop would solve the same chain R times over.
    const bool noAggr = (R == 1) || (naggr >= R - 1);

    // The grouping is a property of the demands, not of the iterate, so it is
    // fixed once here rather than recomputed inside the fixed point.
    std::vector<std::vector<std::size_t>> othersOf(R), grpOf(R);
    for (std::size_t ridx = 0; ridx < R; ++ridx) {
        for (std::size_t oi = 0; oi < R; ++oi)
            if (oi != ridx) othersOf[ridx].push_back(closedChains[oi]);
        if (!noAggr && !othersOf[ridx].empty()) {
            std::vector<std::vector<double>> D(Mc, std::vector<double>(othersOf[ridx].size(), 0.0));
            for (std::size_t ii = 0; ii < Mc; ++ii)
                for (std::size_t oi = 0; oi < othersOf[ridx].size(); ++oi)
                    D[ii][oi] = num_traits<T>::to_double(dem.Lchain(cst[ii], othersOf[ridx][oi]));
            grpOf[ridx] = bgchain_detail::mam_bgchain_groups(D, naggr);
        }
    }
    const std::size_t npass = noAggr ? 1 : R;

    Matrix<T> TNprev(M, K, num_traits<T>::from_double(1e300));
    int totiter = 0;
    const double relax = 0.5;

    auto maxdiff = [&](const Matrix<T>& a, const Matrix<T>& b) {
        double m = 0.0;
        for (std::size_t i = 0; i < a.rows(); ++i)
            for (std::size_t j = 0; j < a.cols(); ++j)
                m = std::max(m, std::fabs(num_traits<T>::to_double(a(i, j)) -
                                          num_traits<T>::to_double(b(i, j))));
        return m;
    };

    while (maxdiff(sol.Tp, TNprev) > opt.tol && totiter < opt.iter_max) {
        ++totiter;
        TNprev = sol.Tp;

        std::vector<double> QopenAcc(M, 0.0);
        std::vector<std::vector<double>> cshareAcc(M, std::vector<double>(Ntot + 1, 0.0));

        for (std::size_t pidx = 0; pidx < npass; ++pidx) {
            const std::size_t r = closedChains[pidx];

            // Background classes: class 0 is the tagged chain, classes 1..G the
            // flow-equivalent aggregates of the demand-similar groups. With no
            // aggregation every closed chain is a class of its own, in chain order.
            std::vector<std::vector<std::size_t>> members;
            if (noAggr) {
                for (std::size_t c : closedChains) members.push_back({c});
            } else {
                members.push_back({r});
                for (std::size_t g = 0; g < naggr; ++g) {
                    std::vector<std::size_t> mem;
                    for (std::size_t oi = 0; oi < othersOf[pidx].size(); ++oi)
                        if (grpOf[pidx][oi] == g) mem.push_back(othersOf[pidx][oi]);
                    members.push_back(mem);
                }
            }
            const std::size_t B = members.size();

            std::vector<int> Nb(B, 0);
            std::vector<std::vector<T>> STb(Mc, std::vector<T>(B, zero));
            std::vector<Matrix<T>> Pb;
            // A class can only hold jobs at the stations its members visit; see
            // mam_bgchain_ctmc on why the union makes the chain reducible.
            std::vector<std::vector<bool>> suppb(Mc, std::vector<bool>(B, false));
            for (std::size_t b = 0; b < B; ++b) {
                const std::vector<std::size_t>& mem = members[b];
                for (std::size_t o : mem) {
                    Nb[b] += static_cast<int>(std::llround(dem.Nchain[o]));
                    for (std::size_t ii = 0; ii < Mc; ++ii)
                        if (num_traits<T>::to_double(dem.Vchain(cst[ii], o)) > 1e-14)
                            suppb[ii][b] = true;
                }
                if (mem.size() == 1) {
                    // a group of one is carried exactly: no mean to take
                    Matrix<T> P0(Mc, Mc, zero);
                    for (std::size_t ii = 0; ii < Mc; ++ii) {
                        STb[ii][b] = dem.STchain(cst[ii], mem[0]);
                        for (std::size_t jj = 0; jj < Mc; ++jj)
                            P0(ii, jj) = Pchain[mem[0]](cst[ii], cst[jj]);
                    }
                    Pb.push_back(P0);
                } else if (mem.empty()) {
                    Pb.push_back(Matrix<T>(Mc, Mc, zero));
                } else {
                    std::vector<std::vector<double>> w(Mc, std::vector<double>(mem.size(), 0.0));
                    for (std::size_t ii = 0; ii < Mc; ++ii) {
                        double rowsum = 0.0;
                        for (std::size_t oi = 0; oi < mem.size(); ++oi) {
                            w[ii][oi] = Xclosed[mem[oi]] *
                                        num_traits<T>::to_double(dem.Vchain(cst[ii], mem[oi]));
                            rowsum += w[ii][oi];
                        }
                        for (std::size_t oi = 0; oi < mem.size(); ++oi)
                            w[ii][oi] = (rowsum > 0.0) ? w[ii][oi] / rowsum
                                                       : 1.0 / static_cast<double>(mem.size());
                    }
                    Matrix<T> Pagg(Mc, Mc, zero);
                    for (std::size_t ii = 0; ii < Mc; ++ii) {
                        T st = zero;
                        for (std::size_t oi = 0; oi < mem.size(); ++oi) {
                            const T wv = num_traits<T>::from_double(w[ii][oi]);
                            st += T(wv * dem.STchain(cst[ii], mem[oi]));
                            for (std::size_t jj = 0; jj < Mc; ++jj)
                                Pagg(ii, jj) += T(wv * Pchain[mem[oi]](cst[ii], cst[jj]));
                        }
                        STb[ii][b] = st;
                    }
                    Pb.push_back(Pagg);
                }
            }
            // Row-normalize, leaving an all-zero row as a self-loop so the chain
            // stays a proper Markov chain on its support.
            for (std::size_t b = 0; b < B; ++b) {
                for (std::size_t ii = 0; ii < Mc; ++ii) {
                    T s = zero;
                    for (std::size_t jj = 0; jj < Mc; ++jj) s += Pb[b](ii, jj);
                    if (num_traits<T>::to_double(s) > 1e-14) {
                        for (std::size_t jj = 0; jj < Mc; ++jj) Pb[b](ii, jj) /= s;
                    } else {
                        for (std::size_t jj = 0; jj < Mc; ++jj) Pb[b](ii, jj) = zero;
                        Pb[b](ii, ii) = one;
                    }
                }
            }

            std::vector<bool> isinfC(Mc, false);
            std::vector<double> nsrvC(Mc, 1.0);
            std::vector<std::vector<double>> cshareC(Mc);
            for (std::size_t ii = 0; ii < Mc; ++ii) {
                isinfC[ii] = isinfAll[cst[ii]];
                nsrvC[ii] = nsrvAll[cst[ii]];
                cshareC[ii] = cshare[cst[ii]];
            }

            const BgchainCtmc<T> bg =
                mam_bgchain_ctmc(Nb, STb, Pb, isinfC, nsrvC, cshareC, suppb, states_max);

            // ---- closed-class metrics of every chain this pass carries EXACTLY:
            // the tagged one always, and every chain when nothing was aggregated.
            const std::size_t bmax = noAggr ? B : 1;
            for (std::size_t b = 0; b < bmax; ++b) {
                if (members[b].empty()) continue;
                const std::size_t rb = members[b][0];
                for (std::size_t k : L.inchain[rb])
                    for (std::size_t i = 0; i < M; ++i) {
                        sol.Q(i, k - 1) = zero;
                        sol.U(i, k - 1) = zero;
                        sol.R(i, k - 1) = zero;
                        sol.Tp(i, k - 1) = zero;
                    }
                for (std::size_t ii = 0; ii < Mc; ++ii) {
                    const std::size_t i = cst[ii];
                    for (std::size_t k : L.inchain[rb]) {
                        const T a = dem.alpha(i, k - 1);
                        if (!(num_traits<T>::to_double(a) > 0.0)) continue;
                        // THROUGHPUT splits by VISIT share, OCCUPANCY by DEMAND
                        // share. A chain queue divided by alpha alone gives every
                        // class of a station the same response time, impossible
                        // at a Delay where R must be the class service time; the
                        // weight is alpha*ST/STchain, as sn_deaggregate_chain_
                        // results applies it.
                        const T stc = dem.STchain(i, rb);
                        const T w = (num_traits<T>::to_double(stc) > 1e-14) ? T(a * S(i, k - 1) / stc) : a;
                        const T q = T(bg.QLen[ii][b] * w);
                        const T x = T(bg.Tput[ii][b] * a);
                        sol.Q(i, k - 1) = q;
                        sol.Tp(i, k - 1) = x;
                        sol.U(i, k - 1) = isinfC[ii] ? q : T(bg.Ubusy[ii][b] * w);
                        sol.R(i, k - 1) = (num_traits<T>::to_double(x) > 1e-14) ? T(q / x) : zero;
                    }
                }
                const std::size_t iref = L.classes[L.inchain[rb][0] - 1].refstat;  // 1-based
                T tputref = zero;
                for (std::size_t k : L.inchain[rb]) tputref += sol.Tp(iref - 1, k - 1);
                const double vref = num_traits<T>::to_double(dem.Vchain(iref - 1, rb));
                Xclosed[rb] = (vref > 1e-14) ? num_traits<T>::to_double(tputref) / vref
                                             : num_traits<T>::to_double(tputref);
            }

            std::vector<double> Uclosed(M, 0.0);
            for (std::size_t ii = 0; ii < Mc; ++ii) {
                double u = 0.0;
                for (std::size_t b = 0; b < B; ++b) u += num_traits<T>::to_double(bg.Ubusy[ii][b]);
                Uclosed[cst[ii]] = u;
            }

            // ---- one pass of the open side ---------------------------------
            for (std::size_t i = 0; i < M; ++i) {
                bool solved = false;
                const SchedStrategy sc = L.stations[i].sched;
                if (sc != SchedStrategy::EXT && sc != SchedStrategy::INF) {
                    std::vector<std::size_t> kopen;
                    for (std::size_t k = 0; k < K; ++k)
                        if (num_traits<T>::to_double(lambdaOpen(i, k)) > 1e-14) kopen.push_back(k);
                    if (!kopen.empty()) {
                        // aggregate open arrival MAP at the station
                        bool haveDa = false;
                        Map<T> Da;
                        for (std::size_t c : openChains) {
                            if (!(lambdaChain[c] > 1e-14)) continue;
                            const auto it = chainArrival.find(c);
                            if (it == chainArrival.end()) continue;
                            T vsum = zero;
                            for (std::size_t k : L.inchain[c]) vsum += V(i, k - 1);
                            const double rate_ic =
                                lambdaChain[c] * num_traits<T>::to_double(vsum);
                            if (!(rate_ic > 1e-14)) continue;
                            const Map<T> scaled =
                                map_scale(it->second, T(num_traits<T>::from_double(1.0 / rate_ic)));
                            if (!haveDa) {
                                Da = scaled;
                                haveDa = true;
                            } else {
                                std::vector<Mmap<T>> two(2);
                                two[0].D0 = Da.D0; two[0].D1 = Da.D1; two[0].Dc.assign(1, Da.D1);
                                two[1].D0 = scaled.D0; two[1].D1 = scaled.D1;
                                two[1].Dc.assign(1, scaled.D1);
                                const Mmap<T> sup = mmap_super_safe(two, opt.space_max);
                                Da = Map<T>{sup.D0, sup.D1};
                            }
                        }
                        if (haveDa) {
                            // arrival-weighted phase-type mixture of the open service laws
                            T lamtot = zero, svcwork = zero;
                            for (std::size_t k : kopen) {
                                lamtot += lambdaOpen(i, k);
                                svcwork += T(lambdaOpen(i, k) * S(i, k));
                            }
                            std::vector<std::vector<T>> pies;
                            std::vector<Matrix<T>> subgens;
                            std::size_t msTotal = 0;
                            // PROCESSOR SHARING IS INSENSITIVE to the service law beyond its mean, so
                            // carrying the phase-type representation at a PS station is not merely unnecessary,
                            // it is WRONG. This QBD tracks ONE service phase for the whole station, which makes
                            // the open queue length inherit the SCV-sensitivity of an M/PH/1 FCFS queue;
                            // measured against SolverCTMC, a HyperExp of SCV 4 then read 21% high where the
                            // exact answer is the exponential one to five digits. The exponential of the same
                            // mean is exact here, and it shrinks the QBD's phase count as a side effect.
                            const bool isPSstation = (sc == SchedStrategy::PS);
                            for (std::size_t k : kopen) {
                                Map<T> phk;
                                if (isPSstation) {
                                    const T mu = T(one / S(i, k));
                                    phk.D0 = Matrix<T>(1, 1, T(-mu));
                                    phk.D1 = Matrix<T>(1, 1, mu);
                                } else {
                                    phk = map_scale(lang::dist_to_map(L.service[i][k]), S(i, k));
                                }
                                std::vector<T> pik = map_pie(phk);
                                const T w = T(lambdaOpen(i, k) / lamtot);
                                for (T& v : pik) v = T(v * w);
                                pies.push_back(pik);
                                subgens.push_back(phk.D0);
                                msTotal += phk.D0.rows();
                            }
                            std::vector<T> alphaS(msTotal, zero);
                            Matrix<T> Tblk(msTotal, msTotal, zero);
                            std::size_t off = 0;
                            for (std::size_t idx = 0; idx < pies.size(); ++idx) {
                                const std::size_t n = subgens[idx].rows();
                                for (std::size_t a = 0; a < n; ++a) {
                                    alphaS[off + a] = pies[idx][a];
                                    for (std::size_t b = 0; b < n; ++b)
                                        Tblk(off + a, off + b) = subgens[idx](a, b);
                                }
                                off += n;
                            }

                            Matrix<T> Aenv(1, 1, zero);
                            std::vector<int> esup(1, 0);
                            const auto pos = std::find(cst.begin(), cst.end(), i);
                            if (pos != cst.end()) {
                                const BgchainEnv<T> env = mam_bgchain_env(
                                    bg, static_cast<std::size_t>(pos - cst.begin()));
                                Aenv = env.A;
                                esup = env.esup;
                            }

                            const std::size_t nphases = Da.D0.rows() * esup.size() * msTotal;
                            if (nphases > phases_max)
                                throw UnsupportedError(
                                    "solver_mam_bgchain: the modulated QBD of station " +
                                    std::to_string(i + 1) + " needs " + std::to_string(nphases) +
                                    " phases (" + std::to_string(Da.D0.rows()) + " arrival x " +
                                    std::to_string(esup.size()) + " environment x " +
                                    std::to_string(msTotal) + " service), above the limit of " +
                                    std::to_string(phases_max) +
                                    ". The environment axis is the closed population held by the "
                                    "station, so it grows with the closed population. Raise "
                                    "options.config.qbdphases_max, or reduce the closed population "
                                    "or the order of the arrival and service processes.");

                            std::vector<double> gref(esup.size(), 0.0);
                            for (std::size_t e = 0; e < esup.size(); ++e)
                                gref[e] = cshare[i][std::min<std::size_t>(
                                    static_cast<std::size_t>(esup[e]),
                                    cshare[i].size() - 1)];

                            // truncation level: the explicit cutoff, else enough
                            // levels for the geometric tail left by the closed
                            // traffic to be negligible
                            std::size_t Kmax;
                            if (opt.cutoff > 0) {
                                Kmax = std::max<std::size_t>(2, opt.cutoff);
                            } else {
                                const double free = std::max(1e-8, 1.0 - Uclosed[i]);
                                const double lam = num_traits<T>::to_double(lamtot);
                                const double Smix = num_traits<T>::to_double(svcwork) / lam;
                                double rho = lam * Smix / (nsrvAll[i] * free);
                                rho = std::min(std::max(rho, 1e-3), 1.0 - 1e-3);
                                const long kv =
                                    static_cast<long>(std::ceil(std::log(1e-8) / std::log(rho)));
                                Kmax = static_cast<std::size_t>(std::min<long>(
                                    std::max<long>(kv, 20), 200));
                            }

                            const BgchainStation<T> st =
                                mam_bgchain_station(Da.D0, Da.D1, alphaS, Tblk, Aenv, esup,
                                                    nsrvAll[i], gref, Kmax);
                            QopenAcc[i] += num_traits<T>::to_double(st.QLen);
                            // The QBD only saw the environment states the chain
                            // reaches; interpolate the rest so the next chain has
                            // a share wherever it may go, clipped to what a closed
                            // job can physically hold.
                            for (int e = 0; e <= Ntot; ++e) {
                                double g;
                                const std::size_t n = st.esup.size();
                                if (n == 1) {
                                    g = st.cshare[0];
                                } else if (e <= st.esup[0]) {
                                    const double slope = (st.cshare[1] - st.cshare[0]) /
                                                         (st.esup[1] - st.esup[0]);
                                    g = st.cshare[0] + slope * (e - st.esup[0]);
                                } else if (e >= st.esup[n - 1]) {
                                    const double slope = (st.cshare[n - 1] - st.cshare[n - 2]) /
                                                         (st.esup[n - 1] - st.esup[n - 2]);
                                    g = st.cshare[n - 1] + slope * (e - st.esup[n - 1]);
                                } else {
                                    std::size_t lo = 0;
                                    while (lo + 1 < n && st.esup[lo + 1] < e) ++lo;
                                    const double w = static_cast<double>(e - st.esup[lo]) /
                                                     static_cast<double>(st.esup[lo + 1] - st.esup[lo]);
                                    g = st.cshare[lo] + w * (st.cshare[lo + 1] - st.cshare[lo]);
                                }
                                cshareAcc[i][e] +=
                                    std::min(std::max(g, 0.0),
                                             std::min(static_cast<double>(e), nsrvAll[i]));
                            }
                            solved = true;
                        }
                    }
                }
                if (!solved) {
                    // a station with no open queue keeps the share it had
                    for (int e = 0; e <= Ntot; ++e) cshareAcc[i][e] += cshare[i][e];
                }
            }
        }

        std::vector<double> Qopen(M, 0.0);
        for (std::size_t i = 0; i < M; ++i) {
            Qopen[i] = QopenAcc[i] / static_cast<double>(npass);
            for (int e = 0; e <= Ntot; ++e)
                cshare[i][e] = (1.0 - relax) * cshare[i][e] +
                               relax * (cshareAcc[i][e] / static_cast<double>(npass));
        }

        // ---- open-class metrics from the aggregate station results --------
        for (std::size_t i = 0; i < M; ++i) {
            std::vector<std::size_t> kopen;
            for (std::size_t k = 0; k < K; ++k)
                if (isopenclass[k] && num_traits<T>::to_double(lambdaOpen(i, k)) > 1e-14)
                    kopen.push_back(k);
            const SchedStrategy sc = L.stations[i].sched;
            if (kopen.empty()) {
                for (std::size_t k = 0; k < K; ++k)
                    if (isopenclass[k]) {
                        sol.Tp(i, k) = lambdaOpen(i, k);
                        sol.Q(i, k) = zero;
                        sol.U(i, k) = zero;
                        sol.R(i, k) = zero;
                    }
                continue;
            }
            T lamtot = zero, work = zero;
            for (std::size_t k : kopen) {
                lamtot += lambdaOpen(i, k);
                work += T(lambdaOpen(i, k) * S(i, k));
            }
            const T Smix = T(work / lamtot);
            for (std::size_t k : kopen) {
                sol.Tp(i, k) = lambdaOpen(i, k);
                if (sc == SchedStrategy::EXT) {
                    sol.Q(i, k) = zero;
                    sol.U(i, k) = zero;
                    sol.R(i, k) = zero;
                } else if (sc == SchedStrategy::INF) {
                    sol.R(i, k) = S(i, k);
                    sol.Q(i, k) = T(lambdaOpen(i, k) * S(i, k));
                    sol.U(i, k) = sol.Q(i, k);
                } else {
                    const T Rtot = T(num_traits<T>::from_double(Qopen[i]) / lamtot);
                    T rk;
                    if (sc == SchedStrategy::PS) {
                        // processor sharing: residence scales with the demand
                        rk = T(Rtot * S(i, k) / Smix);
                    } else {
                        // FCFS and its variants: the wait is class-blind, the
                        // service time is not
                        const T cand = T(Rtot - Smix + S(i, k));
                        rk = (num_traits<T>::to_double(cand) > num_traits<T>::to_double(S(i, k)))
                                 ? cand
                                 : S(i, k);
                    }
                    sol.R(i, k) = rk;
                    sol.Q(i, k) = T(lambdaOpen(i, k) * rk);
                    // Utilization Law: a c-server station holds TN*S/c
                    sol.U(i, k) =
                        T(lambdaOpen(i, k) * S(i, k) / num_traits<T>::from_double(nsrvAll[i]));
                }
            }
        }
    }

    for (std::size_t c = 0; c < C; ++c)
        for (std::size_t k : L.inchain[c])
            sol.X[k - 1] = num_traits<T>::from_double(isopenchain[c] ? lambdaChain[c] : Xclosed[c]);
    for (std::size_t k = 0; k < K; ++k) {
        T acc = zero;
        for (std::size_t i = 0; i < M; ++i) acc += sol.R(i, k);
        sol.C[k] = acc;
    }
    sol.iter = totiter;
    return sol;
    }
}

}  // namespace mam
}  // namespace line

#endif  // LINE_SOLVERS_MAM_SOLVER_MAM_BGCHAIN_H
