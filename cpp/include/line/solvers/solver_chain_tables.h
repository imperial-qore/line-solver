/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_SOLVER_CHAIN_TABLES_H
#define LINE_SOLVERS_SOLVER_CHAIN_TABLES_H

/**
 * The CHAIN-level and SYSTEM-level views of a solved model.
 *
 * Ports of `@@NetworkSolver/getAvgSys.m`, `getAvgChain.m` and
 * `getAvgNodeChain.m` with the tables built on top of them
 * (`getAvgSysTable`, `getAvgChainTable`, `getAvgNodeChainTable`).
 *
 * A CHAIN IS NOT A CLASS AND THE AGGREGATION IS NOT A SUM FOR EVERY METRIC,
 * which is the whole reason these live apart from the AvgTable. Queue lengths,
 * utilizations, arrival rates, throughputs and residence times are ADDITIVE
 * over the classes of a chain, so their chain value is the row sum. Response
 * time is NOT: a chain's response time at a station is the per-visit time
 * averaged over the classes with the VISIT SHARE alpha as the weight, because a
 * job of the chain arrives as one class or another in proportion to how often
 * that class visits. Summing it instead would report the total time a job would
 * spend if it were every class at once.
 *
 * The system view is a third thing again. `getAvgSys` returns one response time
 * and one throughput PER CHAIN, both measured at the chain's reference station:
 * the throughput is the completing flow INTO that station, read off the routing
 * matrix, and the response time is the cycle time -- Little's law on a closed
 * chain, the visit-weighted sum of the per-class system times on an open one.
 *
 * INDEXING: these functions take a solved `mva::AvgResult` (station x class)
 * and the struct it was solved from, and never re-solve. That keeps them usable
 * from any solver whose runner returns an AvgResult, which is what the reference
 * means by putting them on @@NetworkSolver rather than on one solver.
 *
 * ONE DEVIATION FROM THE REFERENCE, STATED: `getAvgSys.m` indexes `sn.rt` and,
 * in one branch, `sn.visits{c}` with a STATION index, while both matrices are
 * indexed by STATEFUL node. The two orders coincide on every model whose
 * stateful nodes are all stations -- which is every model the reference is
 * exercised on -- and differ as soon as one is not (a Cache, a Logger). This
 * port uses `stateful_of_station`, so it agrees with the reference wherever the
 * reference is self-consistent and is correct where it is not.
 */

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/solvers/mva/sn_chain.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace solvers {

namespace chain_detail {

/**
 * Mean of the maximum of independent exponentials, by inclusion-exclusion.
 *
 * With one exponential of rate lambda_i per parallel branch, the expected time
 * until ALL have finished is
 *
 *     E[max] = sum_{k=1..n} (-1)^(k-1) sum_{|S|=k} 1 / (sum_{i in S} lambda_i),
 *
 * which is what `getAvgSys.m` and `pathsCS.m` both spell out with `nchoosek`.
 * The rates come from the branch response times as lambda_i = 1/r_i, so the
 * whole construction reads the parallel section as a race between exponentials
 * whose means are the measured branch times. That is an ASSUMPTION about the
 * branch laws, not a measurement of them, and it is the reference's -- a branch
 * whose time is far from exponential is the one case where this quantity is not
 * the synchronization delay it is reported as.
 *
 * Enumerated over bitmasks rather than by `nchoosek`, which materializes every
 * combination as a matrix row. The cap is on the branch count, not on the
 * subset count: 2^n terms alternate in sign and cancel catastrophically well
 * before memory becomes the issue, and a fork with more than 20 parallel
 * branches is a modelling question, not a numerical one.
 */
template <class T>
T exp_max_mean(const std::vector<T>& branch_times) {
    const std::size_t n = branch_times.size();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (n == 0) return zero;
    if (n > 20)
        throw UnsupportedError(
            "getAvgSys: the fork-join section has " + std::to_string(n) +
            " parallel paths, and the reference's synchronization delay is an "
            "inclusion-exclusion sum of 2^n alternating terms; past 20 paths that sum has lost "
            "every significant digit to cancellation, so it is refused rather than reported");

    std::vector<T> lambda(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double ri = num_traits<T>::to_double(branch_times[i]);
        if (!(ri > 0.0))
            throw NumericError(
                "getAvgSys: a path from the fork to the join has non-positive total response "
                "time, so its exponential rate is not defined; the branch carries no station "
                "with a finite response time");
        lambda[i] = T(one / branch_times[i]);
    }

    T d0 = zero;
    const unsigned long long total = 1ull << n;
    for (unsigned long long mask = 1; mask < total; ++mask) {
        T rate = zero;
        std::size_t bits = 0;
        for (std::size_t i = 0; i < n; ++i)
            if (mask & (1ull << i)) {
                rate = T(rate + lambda[i]);
                ++bits;
            }
        const T term = T(one / rate);
        if (bits % 2 == 1)
            d0 = T(d0 + term);
        else
            d0 = T(d0 - term);
    }
    return d0;
}

/** What one `pathsCS` walk returns. */
template <class T>
struct PathsResult {
    std::vector<T> times;               ///< total response time of each fork-to-join path
    std::vector<std::size_t> stations;  ///< 1-based stations lying on those paths
};

/**
 * Port of `ModelAdapter.pathsCS`: enumerate every path from `cur` to `stop` and
 * total the response time along each.
 *
 * THE ROUTING MATRIX IS `rtnodes` HERE AND `rtorig` IN THE REFERENCE, and the
 * indexing differs with it: `cell2mat(getLinkedRoutingMatrix)` is CLASS-major,
 * addressed as `(class-1)*nnodes + node`, while `rtnodes` is NODE-major,
 * `(node-1)*nclasses + class`. Only the SUPPORT of the matrix is read -- the
 * walk follows nonzero successors and never multiplies by a probability -- and
 * the two describe the same graph once the refresh has resolved class
 * switching, so the enumerated path SET is the same. The order in which the
 * paths come out differs, and nothing downstream depends on it: `exp_max_mean`
 * is symmetric in its argument.
 *
 * A NESTED FORK IS COLLAPSED BEFORE THE WALK CONTINUES, which is why this
 * mutates `RN`: the inner join's response time is set to the inner section's
 * own synchronization delay and the inner branch stations are zeroed, so the
 * outer walk crosses the inner section as a single station. The `RN(join) == 0`
 * test is what makes that happen once rather than once per outer path.
 *
 * A CYCLE IS REFUSED, NOT WALKED. The reference has no guard and recurses until
 * MATLAB runs out of stack; a routing loop inside a fork-join section has no
 * finite path set, so there is nothing to enumerate and saying so is the only
 * available answer.
 */
template <class T>
PathsResult<T> paths_cs(const qn::NetworkStruct<T>& sn, const Matrix<T>& P, std::size_t cur,
                        std::size_t stop, std::size_t cls, Matrix<T>& RN, const T& elapsed,
                        const std::vector<std::size_t>& acc_stations,
                        std::vector<std::pair<std::size_t, std::size_t>>& on_path) {
    const std::size_t K = sn.nclasses, N = sn.nodes.size();
    const T zero = num_traits<T>::from_int(0);

    PathsResult<T> out;
    if (cur == stop) {
        out.times.push_back(elapsed);
        out.stations = acc_stations;
        return out;
    }
    for (std::size_t i = 0; i < on_path.size(); ++i)
        if (on_path[i].first == cur && on_path[i].second == cls)
            throw UnsupportedError(
                "getAvgSys: the fork-join section of node '" + sn.nodes[cur - 1].name +
                "' contains a routing cycle, so the set of paths from the fork to the join is "
                "infinite and the synchronization delay has no inclusion-exclusion form");
    on_path.push_back(std::make_pair(cur, cls));

    T here = zero;
    std::vector<std::size_t> stations = acc_stations;
    const std::size_t ist = sn.nodes[cur - 1].station;
    if (ist != 0) {
        here = RN(ist - 1, cls - 1);
        stations.push_back(ist);
    }

    const std::size_t row = (cur - 1) * K + (cls - 1);
    for (std::size_t nxt = 1; nxt <= N; ++nxt) {
        for (std::size_t s = 1; s <= K; ++s) {
            if (row >= P.rows()) break;
            const std::size_t col = (nxt - 1) * K + (s - 1);
            if (col >= P.cols()) continue;
            if (num_traits<T>::to_double(P(row, col)) == 0.0) continue;

            std::size_t hop = nxt;
            T entry = T(elapsed + here);
            if (sn.nodes[nxt - 1].nodetype == lang::NodeType::Fork) {
                std::size_t inner_join = 0;
                for (std::size_t k = 0; k < sn.fj.size(); ++k)
                    if (sn.fj[k].first == nxt) inner_join = sn.fj[k].second;
                if (inner_join != 0) {
                    const std::size_t jst = sn.nodes[inner_join - 1].station;
                    if (jst != 0 && num_traits<T>::to_double(RN(jst - 1, s - 1)) == 0.0) {
                        std::vector<std::pair<std::size_t, std::size_t>> inner_path;
                        const PathsResult<T> in =
                            paths_cs(sn, P, nxt, inner_join, s, RN, zero,
                                     std::vector<std::size_t>(), inner_path);
                        RN(jst - 1, s - 1) = exp_max_mean(in.times);
                        for (std::size_t q = 0; q < in.stations.size(); ++q)
                            RN(in.stations[q] - 1, s - 1) = zero;
                    }
                    hop = inner_join;
                }
            }
            const PathsResult<T> sub = paths_cs(sn, P, hop, stop, s, RN, entry, stations, on_path);
            out.times.insert(out.times.end(), sub.times.begin(), sub.times.end());
            out.stations.insert(out.stations.end(), sub.stations.begin(), sub.stations.end());
        }
    }
    on_path.pop_back();
    return out;
}

}  // namespace chain_detail

/** `@@NetworkSolver/getAvgSys`: one response time and one throughput per chain. */
template <class T>
struct SysResult {
    std::vector<T> CN;  ///< (nchains) system response time, i.e. the cycle time
    std::vector<T> XN;  ///< (nchains) system throughput at the reference station
};

/** The station- or node-level table aggregated by chain. */
template <class T>
struct ChainResult {
    Matrix<T> QN, UN, RN, WN, AN, TN;  ///< (rows x nchains), rows = stations or nodes
};

/**
 * Port of `@@NetworkSolver/getAvgSys.m`.
 *
 * FORK-JOIN WITH AN OPEN CHAIN IS SERVED (2026-08-15). The reference fills the
 * join station's response time with the order statistic of the parallel branch
 * times -- `d0`, an inclusion-exclusion sum over every path from the fork to
 * the join, enumerated by `ModelAdapter.pathsCS` -- and both halves of that are
 * ported above as `chain_detail::paths_cs` and `chain_detail::exp_max_mean`.
 * A CLOSED fork-join chain needs neither: its cycle time comes from Little's
 * law, which reads the population and the throughput and never touches the
 * join's response time, and the reference skips the walk there too. The JAR
 * still refuses this case and then fills RN with NaN.
 */
template <class T>
SysResult<T> solver_get_avg_sys(const qn::NetworkStruct<T>& sn, const mva::AvgResult<T>& r) {
    const std::size_t M = sn.nstations, K = sn.nclasses, C = sn.nchains;
    const T zero = num_traits<T>::from_int(0);
    const std::vector<double> njobs = sn.njobs();

    bool has_join = false, has_fork = false;
    for (std::size_t i = 0; i < sn.nodes.size(); ++i) {
        if (sn.nodes[i].nodetype == lang::NodeType::Join) has_join = true;
        if (sn.nodes[i].nodetype == lang::NodeType::Fork) has_fork = true;
    }

    // `RN(join, :) = 0` of the reference, applied before anything reads RN: the
    // join holds no service, and whatever the solver reported there is the
    // synchronization wait, which the cycle time accounts for at the fork.
    Matrix<T> RN = r.RN;
    if (has_join)
        for (std::size_t i = 0; i < M; ++i)
            if (sn.stations[i].nodetype == lang::NodeType::Join)
                for (std::size_t c = 0; c < K; ++c) RN(i, c) = zero;

    // ---- The synchronization delay of each fork-join section, on OPEN chains
    // only. A closed chain gets its cycle time from Little's law, which never
    // reads the join's response time, so the reference skips the whole walk
    // there and so does this.
    //
    // ORDER MATTERS: this rewrites RN -- the join gets the delay, the branch
    // stations get zero -- and CNclass below reads RN. Computing CNclass first
    // would total the branch times ALONG the paths, which double-counts the
    // parallel section instead of taking its maximum.
    if (has_fork && has_join) {
        for (std::size_t f = 0; f < sn.fj.size(); ++f) {
            const std::size_t fork = sn.fj[f].first, join = sn.fj[f].second;
            if (fork == 0 || join == 0) continue;
            const std::size_t jst = sn.nodes[join - 1].station;
            if (jst == 0) continue;
            for (std::size_t c = 0; c < C; ++c) {
                double nJobsChain = 0.0;
                for (std::size_t k = 0; k < K; ++k)
                    if (sn.chains[c][k]) nJobsChain += njobs[k];
                if (!std::isinf(nJobsChain)) continue;
                for (std::size_t q = 0; q < sn.inchain[c].size(); ++q) {
                    const std::size_t rr = sn.inchain[c][q];
                    if (num_traits<T>::to_double(RN(jst - 1, rr - 1)) != 0.0) continue;
                    std::vector<std::pair<std::size_t, std::size_t>> on_path;
                    const chain_detail::PathsResult<T> paths =
                        chain_detail::paths_cs(sn, sn.rtnodes, fork, join, rr, RN, zero,
                                               std::vector<std::size_t>(), on_path);
                    if (paths.times.empty()) continue;
                    // d0 already accounts for the time spent on the branches,
                    // which is why they are zeroed rather than left to be added.
                    RN(jst - 1, rr - 1) = chain_detail::exp_max_mean(paths.times);
                    for (std::size_t j = 0; j < paths.stations.size(); ++j)
                        RN(paths.stations[j] - 1, rr - 1) = zero;
                }
            }
        }
    }

    // ---- CNclass: the per-class system time, visits-weighted to the reference
    // station. Computed for every model, used only by the open branch below,
    // exactly as the reference computes it.
    std::vector<T> CNclass(K, zero);
    for (std::size_t c = 0; c < C; ++c) {
        for (std::size_t j = 0; j < sn.inchain[c].size(); ++j) {
            const std::size_t rr = sn.inchain[c][j];  // 1-based class
            const std::size_t refst = sn.classes[rr - 1].refstat;
            if (refst == 0) continue;
            const std::size_t refsf = sn.stateful_of_station(refst);
            const T vref = sn.visits[c](refsf - 1, rr - 1);
            if (num_traits<T>::to_double(vref) == 0.0) continue;
            for (std::size_t i = 0; i < M; ++i) {
                // The source of an open class is not part of its system time.
                if (std::isinf(njobs[rr - 1]) && i + 1 == refst) continue;
                const std::size_t isf = sn.stateful_of_station(i + 1);
                CNclass[rr - 1] =
                    T(CNclass[rr - 1] + T(T(sn.visits[c](isf - 1, rr - 1) * RN(i, rr - 1)) / vref));
            }
        }
    }

    // ---- alpha: the share of its chain's reference-station completions that a
    // class accounts for, per station. `refclass` is indexed by CHAIN and
    // `inchain` holds CLASS indices; the reference intersects the two anyway,
    // and this port reproduces that intersection rather than repairing it,
    // because the weights it produces are the ones every reference number was
    // computed with.
    Matrix<T> alpha(M, K, zero);
    std::vector<std::size_t> refclass_nz;
    for (std::size_t c = 0; c < C; ++c)
        if (sn.refclass[c] != 0) refclass_nz.push_back(c + 1);
    for (std::size_t c = 0; c < C; ++c) {
        // The classes of the chain that COMPLETE at the reference station: a
        // non-completing class passes through without ending a cycle, so it
        // contributes visits but no completions to divide by.
        std::vector<std::size_t> completing;
        for (std::size_t k = 0; k < K; ++k)
            if (sn.chains[c][k] && sn.classes[k].completes) completing.push_back(k + 1);

        std::vector<std::size_t> ks;
        for (std::size_t j = 0; j < refclass_nz.size(); ++j)
            for (std::size_t q = 0; q < sn.inchain[c].size(); ++q)
                if (refclass_nz[j] == sn.inchain[c][q]) ks.push_back(refclass_nz[j]);
        if (sn.refclass[c] == 0) ks = sn.inchain[c];

        for (std::size_t i = 0; i < M; ++i) {
            for (std::size_t j = 0; j < ks.size(); ++j) {
                const std::size_t k = ks[j];
                const std::size_t refst = sn.classes[k - 1].refstat;
                if (refst == 0) continue;
                const std::size_t refsf = sn.stateful_of_station(refst);
                T denom = zero;
                for (std::size_t q = 0; q < completing.size(); ++q)
                    denom = T(denom + sn.visits[c](refsf - 1, completing[q] - 1));
                if (num_traits<T>::to_double(denom) == 0.0) continue;
                const std::size_t isf = sn.stateful_of_station(i + 1);
                alpha(i, k - 1) = T(alpha(i, k - 1) + T(sn.visits[c](isf - 1, k - 1) / denom));
            }
        }
    }
    // `alpha(~isfinite(alpha)) = 0`: a class that never reaches the reference
    // station has no share, and a non-finite weight would poison the sum.
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k)
            if (!std::isfinite(num_traits<T>::to_double(alpha(i, k)))) alpha(i, k) = zero;

    SysResult<T> out;
    out.CN.assign(C, zero);
    out.XN.assign(C, zero);
    const std::size_t nsf = sn.nof_stateful();

    for (std::size_t c = 0; c < C; ++c) {
        const std::vector<std::size_t>& inchain = sn.inchain[c];
        if (inchain.empty()) continue;
        std::vector<std::size_t> completing;
        for (std::size_t k = 0; k < K; ++k)
            if (sn.chains[c][k] && sn.classes[k].completes) completing.push_back(k + 1);

        // ---- XN: the completing flow INTO the chain's reference station.
        // Read off the routing matrix rather than from any one station's
        // throughput, because a chain completes wherever its routing returns to
        // the reference station and that can be several places at once.
        const std::size_t ref = sn.classes[inchain[0] - 1].refstat;
        if (ref != 0 && r.TN.rows() == M) {
            const std::size_t refsf = sn.stateful_of_station(ref);
            std::vector<std::size_t> ss;
            for (std::size_t j = 0; j < refclass_nz.size(); ++j)
                for (std::size_t q = 0; q < inchain.size(); ++q)
                    if (refclass_nz[j] == inchain[q]) ss.push_back(refclass_nz[j]);
            if (ss.empty()) ss = inchain;
            for (std::size_t i = 0; i < M; ++i) {
                const std::size_t isf = sn.stateful_of_station(i + 1);
                if (isf == 0 || isf > nsf) continue;
                for (std::size_t q = 0; q < completing.size(); ++q) {
                    const std::size_t rr = completing[q];
                    const double tn = num_traits<T>::to_double(r.TN(i, rr - 1));
                    if (std::isnan(tn)) continue;
                    for (std::size_t j = 0; j < ss.size(); ++j) {
                        const std::size_t s = ss[j];
                        out.XN[c] = T(out.XN[c] + T(sn.rt((isf - 1) * K + (rr - 1),
                                                          (refsf - 1) * K + (s - 1)) *
                                                   r.TN(i, rr - 1)));
                    }
                }
            }
        }

        // ---- CN: Little's law on a closed chain, the alpha-weighted sum of the
        // per-class system times on an open one.
        double nJobsChain = 0.0;
        for (std::size_t k = 0; k < K; ++k)
            if (sn.chains[c][k]) nJobsChain += njobs[k];

        if (std::isinf(nJobsChain)) {
            if (inchain.size() != completing.size())
                throw UnsupportedError(
                    "getAvgSys: edge-based chain definition is not supported for open queueing "
                    "networks -- the chain holds a non-completing class, so there is no single "
                    "flow whose reciprocal is the cycle time");
            const std::size_t refst = sn.classes[inchain[0] - 1].refstat;
            T acc = zero;
            for (std::size_t j = 0; j < inchain.size(); ++j) {
                const T v = T(alpha(refst - 1, inchain[j] - 1) * CNclass[inchain[j] - 1]);
                // `sumfinite`: a class the solver reported no finite time for is
                // skipped, not propagated as Inf over the whole chain.
                if (std::isfinite(num_traits<T>::to_double(v))) acc = T(acc + v);
            }
            out.CN[c] = acc;
        } else {
            const double x = num_traits<T>::to_double(out.XN[c]);
            out.CN[c] = x == 0.0 ? zero
                                 : T(num_traits<T>::from_double(nJobsChain) / out.XN[c]);
        }
    }
    return out;
}

/**
 * Port of `@@NetworkSolver/getAvgChain.m`: the station table aggregated by chain.
 *
 * QLen, Util, ArvR, Tput and ResidT are row sums over the chain's classes;
 * RespT is the alpha-weighted average, alpha being `sn_get_demands_chain`'s
 * visit share. See the file header for why the two rules differ.
 */
template <class T>
ChainResult<T> solver_get_avg_chain(const qn::NetworkStruct<T>& sn,
                                    const mva::AvgResult<T>& r) {
    const std::size_t M = sn.nstations, C = sn.nchains;
    const T zero = num_traits<T>::from_int(0);
    const mva::ChainDemands<T> dem = mva::sn_get_demands_chain(sn);

    ChainResult<T> out;
    out.QN = Matrix<T>(M, C, zero);
    out.UN = Matrix<T>(M, C, zero);
    out.RN = Matrix<T>(M, C, zero);
    out.WN = Matrix<T>(M, C, zero);
    out.AN = Matrix<T>(M, C, zero);
    out.TN = Matrix<T>(M, C, zero);
    for (std::size_t c = 0; c < C; ++c) {
        for (std::size_t j = 0; j < sn.inchain[c].size(); ++j) {
            const std::size_t k = sn.inchain[c][j] - 1;
            for (std::size_t i = 0; i < M; ++i) {
                out.QN(i, c) = T(out.QN(i, c) + r.QN(i, k));
                out.UN(i, c) = T(out.UN(i, c) + r.UN(i, k));
                out.WN(i, c) = T(out.WN(i, c) + r.WN(i, k));
                out.AN(i, c) = T(out.AN(i, c) + r.AN(i, k));
                out.TN(i, c) = T(out.TN(i, c) + r.TN(i, k));
                out.RN(i, c) = T(out.RN(i, c) + T(r.RN(i, k) * dem.alpha(i, k)));
            }
        }
    }
    return out;
}

/**
 * Port of `@@NetworkSolver/getAvgNodeChain.m`: the NODE table aggregated by chain.
 *
 * The node-level class matrices are the caller's, because the scatter from
 * stations to nodes and the recomputation of ArvR and Tput per node is
 * `getAvgNodeTable`'s work and is not repeated here. Rows that are not stations
 * carry zero response and residence time, which is the reference's construction
 * and not a gap: a node that is not a station holds no jobs.
 */
template <class T>
ChainResult<T> solver_get_avg_node_chain(const qn::NetworkStruct<T>& sn, const Matrix<T>& QNn,
                                         const Matrix<T>& UNn, const Matrix<T>& RNn,
                                         const Matrix<T>& WNn, const Matrix<T>& ANn,
                                         const Matrix<T>& TNn) {
    const std::size_t I = sn.nodes.size(), C = sn.nchains;
    const T zero = num_traits<T>::from_int(0);
    const mva::ChainDemands<T> dem = mva::sn_get_demands_chain(sn);

    ChainResult<T> out;
    out.QN = Matrix<T>(I, C, zero);
    out.UN = Matrix<T>(I, C, zero);
    out.RN = Matrix<T>(I, C, zero);
    out.WN = Matrix<T>(I, C, zero);
    out.AN = Matrix<T>(I, C, zero);
    out.TN = Matrix<T>(I, C, zero);
    // Node -> station, so the alpha weights (which are indexed by station) can
    // be applied to the two time columns.
    std::vector<std::size_t> node_to_station(I, 0);
    for (std::size_t ist = 0; ist < sn.nstations; ++ist) {
        const std::size_t ind = sn.station_to_node[ist];
        if (ind) node_to_station[ind - 1] = ist + 1;
    }
    for (std::size_t c = 0; c < C; ++c) {
        for (std::size_t j = 0; j < sn.inchain[c].size(); ++j) {
            const std::size_t k = sn.inchain[c][j] - 1;
            for (std::size_t i = 0; i < I; ++i) {
                out.QN(i, c) = T(out.QN(i, c) + QNn(i, k));
                out.UN(i, c) = T(out.UN(i, c) + UNn(i, k));
                out.AN(i, c) = T(out.AN(i, c) + ANn(i, k));
                out.TN(i, c) = T(out.TN(i, c) + TNn(i, k));
                const std::size_t ist = node_to_station[i];
                if (!ist) continue;
                out.RN(i, c) = T(out.RN(i, c) + T(RNn(i, k) * dem.alpha(ist - 1, k)));
                out.WN(i, c) = T(out.WN(i, c) + T(WNn(i, k) * dem.alpha(ist - 1, k)));
            }
        }
    }
    return out;
}

/** `Chain1`, `Chain2`, ... -- the reference's own chain labels. */
inline std::vector<std::string> chain_names(std::size_t nchains) {
    std::vector<std::string> out;
    for (std::size_t c = 0; c < nchains; ++c) out.push_back("Chain" + std::to_string(c + 1));
    return out;
}

/** `(ClassA ClassB)`, the JobClasses column: which classes a chain holds. */
template <class T>
std::vector<std::string> chain_class_labels(const qn::NetworkStruct<T>& sn) {
    std::vector<std::string> out;
    for (std::size_t c = 0; c < sn.nchains; ++c) {
        std::string s = "(";
        for (std::size_t j = 0; j < sn.inchain[c].size(); ++j) {
            if (j) s += " ";
            s += sn.classes[sn.inchain[c][j] - 1].name;
        }
        out.push_back(s + ")");
    }
    return out;
}

}  // namespace solvers
}  // namespace line

#endif  // LINE_SOLVERS_SOLVER_CHAIN_TABLES_H
