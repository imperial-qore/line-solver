/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_BA_SOLVER_BA_SNC_H
#define LINE_SOLVERS_BA_SOLVER_BA_SNC_H

/**
 * Stochastic network calculus UPPER bound on the mean response times and queue
 * lengths of a feed-forward open network, valid for EVERY work-conserving
 * scheduling policy at every station.
 *
 * Port of matlab/src/solvers/BA/solver_ba_snc_analyzer.m, cross-checked against
 * jline.solvers.ba.analyzers.Solver_ba_snc_analyzer and the native-Python
 * solver_ba_snc.py. The api domain `line/api/snc` supplies the envelope
 * algebra; this analyzer maps the LINE model onto it, propagates envelopes hop
 * by hop, and reads the bound back per station and class.
 *
 * UNITS ARE JOBS, NOT WORK. The arrival envelope counts jobs and the service
 * element is `snc_srv_exp`, the counting process of an Exp(mu) server. That is
 * what lets a departure envelope from one station be the arrival envelope of the
 * next: a service-time work unit differs from station to station, a job does
 * not. On a single M/M/1 the resulting backlog bound decays as (lambda/mu)^n and
 * the delay bound as exp(-(mu-lambda)*d), both exact rates.
 *
 * BOUND CONVENTION. R(i,r) is `snc_mean_delay` of the (arrival, service)
 * envelope pair, i.e. the integral of the delay tail bound, so each entry is a
 * valid upper bound on its own. Q follows by Little's law from the bounded R and
 * the EXACT throughput T (an open network's per-class rates are fixed by the
 * traffic equations, not by the policy), and so does C. U is exact for the same
 * reason.
 *
 * ARITHMETIC. The api domain is double-only -- every bound is an exp/log
 * expression minimized numerically over theta -- so this arm converts at the
 * boundary with `num_traits<T>::to_double` and `from_double` rather than
 * pretending to be Rational-clean. It is the one BA family that is not
 * instantiable at exact arithmetic, and `registry.h` registers the domain as
 * Double for that reason.
 *
 * Reference: M. Fidler, A. Rizk (2015). A Guide to the Stochastic Network
 * Calculus. IEEE Communications Surveys and Tutorials 17(1), 92-105.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "line/api/sn/sn_rt_stations.h"
#include "line/api/snc/snc_env_map.h"
#include "line/api/snc/snc_env_poisson.h"
#include "line/api/snc/snc_leftover.h"
#include "line/api/snc/snc_mean_delay.h"
#include "line/api/snc/snc_output.h"
#include "line/api/snc/snc_perc_backlog.h"
#include "line/api/snc/snc_perc_delay.h"
#include "line/api/snc/snc_srv_exp.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace ba {

/** The per-pair envelopes the analyzer built, keyed by (station, class). */
struct SncEnvelopes {
    std::map<std::pair<std::size_t, std::size_t>, snc::Envelope> arv;
    std::map<std::pair<std::size_t, std::size_t>, snc::Envelope> srv;
    std::map<std::pair<std::size_t, std::size_t>, double> lam;
    std::map<std::pair<std::size_t, std::size_t>, double> mu;
    std::size_t M = 0, K = 0;
};

namespace detail {

/** Superposition of independent flows: the exponential forms multiply. */
inline snc::Envelope snc_sum(const std::vector<snc::Envelope>& parts) {
    return [parts](double theta) {
        snc::Env e{0.0, 0.0};
        for (std::size_t k = 0; k < parts.size(); ++k) {
            const snc::Env p = parts[k](theta);
            e.sigma += p.sigma;
            e.rho += p.rho;
        }
        return e;
    };
}

/** Kahn's algorithm; empty result means the graph has a cycle. */
inline std::vector<std::size_t> snc_topo_order(const std::vector<std::vector<bool>>& adj) {
    const std::size_t n = adj.size();
    std::vector<int> indeg(n, 0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            if (adj[i][j]) indeg[j]++;
    std::vector<bool> done(n, false);
    std::vector<std::size_t> order;
    while (true) {
        std::size_t cand = n;
        for (std::size_t i = 0; i < n && cand == n; ++i)
            if (!done[i] && indeg[i] == 0) cand = i;
        if (cand == n) break;
        order.push_back(cand);
        done[cand] = true;
        for (std::size_t j = 0; j < n; ++j)
            if (adj[cand][j]) indeg[j]--;
        indeg[cand] = 1;  // keep it out of the candidate set
    }
    if (order.size() < n) return {};
    return order;
}

}  // namespace detail

/**
 * Builds the per-pair (arrival, service) envelopes of a feed-forward model.
 *
 * Every gate of the family is checked here rather than in the caller, so that
 * the quantile accessors refuse an unsupported model with the same reason as
 * the mean columns.
 *
 * @param L the model
 */
template <class T>
SncEnvelopes solver_ba_snc_envelopes(const qn::NetworkStruct<T>& L) {
    const std::size_t M = L.nstations, K = L.nclasses;
    const double TOL = 1e-10;
    SncEnvelopes env;
    env.M = M;
    env.K = K;

    // ---- model gates ----
    for (std::size_t r = 0; r < L.classes.size(); ++r)
        if (std::isfinite(L.classes[r].population))
            throw UnsupportedError(
                "solver_ba_snc: method 'snc.upper' supports fully open networks only "
                "(no closed classes)");
    std::vector<std::size_t> srcList, qstat;
    for (std::size_t i = 0; i < M; ++i) {
        if (L.stations[i].nodetype == qn::NodeType::Source)
            srcList.push_back(i);
        else
            qstat.push_back(i);
    }
    if (srcList.empty())
        throw UnsupportedError(
            "solver_ba_snc: method 'snc.upper' requires an open network with a Source station");
    for (std::size_t a = 0; a < qstat.size(); ++a) {
        const std::size_t i = qstat[a];
        if (L.stations[i].sched == lang::SchedStrategy::INF)
            throw UnsupportedError(
                "solver_ba_snc: method 'snc.upper' does not support delay (infinite-server) "
                "stations: the service envelope is that of a single busy server");
        const double ns = num_traits<T>::to_double(L.stations[i].nservers);
        if (std::isfinite(ns) && ns > 1)
            throw UnsupportedError(
                "solver_ba_snc: method 'snc.upper' does not support multi-server stations");
    }

    // ---- station-space routing, with the Source absorbed into the injections ----
    const Matrix<T> rtst = api::sn_rt_stations(L).rtst;

    std::vector<std::size_t> pairStation, pairClass, pairFlat;
    for (std::size_t a = 0; a < qstat.size(); ++a)
        for (std::size_t r = 0; r < K; ++r) {
            pairStation.push_back(qstat[a]);
            pairClass.push_back(r);
            pairFlat.push_back(qstat[a] * K + r);
        }
    const std::size_t np = pairFlat.size();

    // Injections are kept per (source, class) so that the exogenous process of
    // each stream is still identifiable once the pairs are known.
    std::vector<std::size_t> srcOfCol, clsOfCol;
    for (std::size_t si = 0; si < srcList.size(); ++si)
        for (std::size_t r0 = 0; r0 < K; ++r0) {
            srcOfCol.push_back(srcList[si]);
            clsOfCol.push_back(r0);
        }
    const std::size_t ncols = srcOfCol.size();
    std::vector<std::vector<double>> inject(np, std::vector<double>(ncols, 0.0));
    std::vector<double> lambda0(np, 0.0);
    for (std::size_t c = 0; c < ncols; ++c) {
        const double arr = num_traits<T>::to_double(L.rates(srcOfCol[c], clsOfCol[c]));
        if (!std::isfinite(arr) || arr <= 0) continue;
        const std::size_t srow = srcOfCol[c] * K + clsOfCol[c];
        for (std::size_t p = 0; p < np; ++p) {
            inject[p][c] = arr * num_traits<T>::to_double(rtst(srow, pairFlat[p]));
            lambda0[p] += inject[p][c];
        }
    }

    Matrix<double> P(np, np, 0.0);
    for (std::size_t p = 0; p < np; ++p)
        for (std::size_t q = 0; q < np; ++q)
            P(p, q) = num_traits<T>::to_double(rtst(pairFlat[p], pairFlat[q]));

    // ---- restrict to the pairs that actually carry traffic ----
    Matrix<double> ImPt(np, np, 0.0);
    for (std::size_t i = 0; i < np; ++i)
        for (std::size_t j = 0; j < np; ++j) ImPt(i, j) = (i == j ? 1.0 : 0.0) - P(j, i);
    Matrix<double> rhs0(np, 1, 0.0);
    for (std::size_t p = 0; p < np; ++p) rhs0(p, 0) = lambda0[p];
    const Matrix<double> lamAll = matmul(inverse(ImPt), rhs0);
    double lamMax = 0.0;
    for (std::size_t p = 0; p < np; ++p) lamMax = std::max(lamMax, lamAll(p, 0));
    const double lamTol = 1e-12 * std::max(1.0, lamMax);
    std::vector<std::size_t> keep;
    for (std::size_t p = 0; p < np; ++p)
        if (lamAll(p, 0) > lamTol) keep.push_back(p);
    if (keep.empty()) throw UnsupportedError("solver_ba_snc: the model carries no open traffic");

    const std::size_t nk = keep.size();
    std::vector<double> lam(nk, 0.0), mu(nk, 0.0);
    std::vector<std::size_t> statk(nk), clsk(nk);
    std::vector<std::vector<double>> injk(nk, std::vector<double>(ncols, 0.0));
    Matrix<double> Pk(nk, nk, 0.0);
    for (std::size_t a = 0; a < nk; ++a) {
        const std::size_t p = keep[a];
        lam[a] = lamAll(p, 0);
        statk[a] = pairStation[p];
        clsk[a] = pairClass[p];
        injk[a] = inject[p];
        for (std::size_t b = 0; b < nk; ++b) Pk(a, b) = P(p, keep[b]);
        mu[a] = num_traits<T>::to_double(L.rates(statk[a], clsk[a]));
        if (!std::isfinite(mu[a]) || mu[a] <= 0)
            throw UnsupportedError("solver_ba_snc: station " + std::to_string(statk[a] + 1) +
                                   " has no service rate for class " + std::to_string(clsk[a] + 1) +
                                   " but carries its traffic");
        if (L.procid(statk[a] + 1, clsk[a] + 1) != lang::ProcessType::EXP)
            throw UnsupportedError(
                "solver_ba_snc: method 'snc.upper' requires exponential service: station " +
                std::to_string(statk[a] + 1) + " class " + std::to_string(clsk[a] + 1) +
                " is not exponential");
    }

    // ---- routing restrictions: no split downstream of the Source ----
    for (std::size_t a = 0; a < nk; ++a) {
        std::size_t nsucc = 0, succ = 0;
        for (std::size_t b = 0; b < nk; ++b)
            if (Pk(a, b) > TOL) {
                nsucc++;
                succ = b;
            }
        if (nsucc > 1)
            throw UnsupportedError(
                "solver_ba_snc: method 'snc.upper' requires deterministic routing downstream of "
                "the Source: station " +
                std::to_string(statk[a] + 1) + " class " + std::to_string(clsk[a] + 1) +
                " splits its flow over " + std::to_string(nsucc) + " destinations");
        if (nsucc == 1 && std::fabs(Pk(a, succ) - 1.0) > 1e-8)
            throw UnsupportedError(
                "solver_ba_snc: method 'snc.upper' requires deterministic routing downstream of "
                "the Source: station " +
                std::to_string(statk[a] + 1) + " class " + std::to_string(clsk[a] + 1) +
                " routes onward with probability " + std::to_string(Pk(a, succ)));
    }
    // A Source that splits is exact only when its process is Poisson, since a
    // Bernoulli thinning of a Poisson stream is again Poisson.
    for (std::size_t c = 0; c < ncols; ++c) {
        std::size_t ndest = 0;
        for (std::size_t a = 0; a < nk; ++a)
            if (injk[a][c] > TOL) ndest++;
        if (ndest <= 1) continue;
        if (L.procid(srcOfCol[c] + 1, clsOfCol[c] + 1) != lang::ProcessType::EXP)
            throw UnsupportedError(
                "solver_ba_snc: method 'snc.upper' can split only a Poisson Source: source " +
                std::to_string(srcOfCol[c] + 1) + " class " + std::to_string(clsOfCol[c] + 1) +
                " is not exponential and feeds " + std::to_string(ndest) + " stations");
    }

    // ---- one service rate per station, and a feed-forward station graph ----
    std::vector<std::size_t> stationsUsed;
    for (std::size_t a = 0; a < nk; ++a)
        if (std::find(stationsUsed.begin(), stationsUsed.end(), statk[a]) == stationsUsed.end())
            stationsUsed.push_back(statk[a]);
    for (std::size_t u = 0; u < stationsUsed.size(); ++u) {
        double lo = std::numeric_limits<double>::infinity(), hi = 0.0;
        for (std::size_t a = 0; a < nk; ++a)
            if (statk[a] == stationsUsed[u]) {
                lo = std::min(lo, mu[a]);
                hi = std::max(hi, mu[a]);
            }
        if (hi - lo > 1e-8 * std::max(1.0, hi))
            throw UnsupportedError(
                "solver_ba_snc: method 'snc.upper' requires the classes sharing a station to have "
                "equal service rates: station " +
                std::to_string(stationsUsed[u] + 1) + " carries rates in [" + std::to_string(lo) +
                ", " + std::to_string(hi) + "]");
    }

    const std::size_t ns = stationsUsed.size();
    std::vector<std::vector<bool>> adj(ns, std::vector<bool>(ns, false));
    const auto stIdx = [&](std::size_t station) {
        return static_cast<std::size_t>(
            std::find(stationsUsed.begin(), stationsUsed.end(), station) - stationsUsed.begin());
    };
    for (std::size_t a = 0; a < nk; ++a)
        for (std::size_t b = 0; b < nk; ++b)
            if (Pk(a, b) > TOL) adj[stIdx(statk[a])][stIdx(statk[b])] = true;
    const std::vector<std::size_t> order = detail::snc_topo_order(adj);
    if (order.empty())
        throw UnsupportedError(
            "solver_ba_snc: method 'snc.upper' requires a feed-forward network: the station graph "
            "has a cycle, so a station's cross traffic is not determined upstream of it");

    // ---- envelope propagation, station by station in feed-forward order ----
    // The departure envelopes live in a shared vector read through a pointer, so
    // that a predecessor filled later in the same station pass is still seen.
    auto outH = std::make_shared<std::vector<snc::Envelope>>(nk);
    std::vector<snc::Envelope> arvH(nk), srvH(nk);
    for (std::size_t oi = 0; oi < ns; ++oi) {
        const std::size_t i = stationsUsed[order[oi]];
        std::vector<std::size_t> here;
        for (std::size_t a = 0; a < nk; ++a)
            if (statk[a] == i) here.push_back(a);

        for (std::size_t h = 0; h < here.size(); ++h) {
            const std::size_t a = here[h];
            std::vector<snc::Envelope> parts;
            for (std::size_t c = 0; c < ncols; ++c) {
                if (injk[a][c] <= TOL) continue;
                if (L.procid(srcOfCol[c] + 1, clsOfCol[c] + 1) == lang::ProcessType::EXP) {
                    parts.push_back(snc::snc_env_poisson_fn(injk[a][c]));
                } else {
                    const mam::Map<T> m =
                        lang::dist_to_map(L.service[srcOfCol[c]][clsOfCol[c]]);
                    Matrix<double> D0(m.D0.rows(), m.D0.cols(), 0.0);
                    Matrix<double> D1(m.D1.rows(), m.D1.cols(), 0.0);
                    for (std::size_t x = 0; x < m.D0.rows(); ++x)
                        for (std::size_t y = 0; y < m.D0.cols(); ++y) {
                            D0(x, y) = num_traits<T>::to_double(m.D0(x, y));
                            D1(x, y) = num_traits<T>::to_double(m.D1(x, y));
                        }
                    parts.push_back(snc::snc_env_map_fn(D0, D1));
                }
            }
            for (std::size_t b = 0; b < nk; ++b)
                if (Pk(b, a) > TOL)
                    parts.push_back([outH, b](double theta) { return (*outH)[b](theta); });
            if (parts.empty())
                throw UnsupportedError("solver_ba_snc: station " + std::to_string(i + 1) +
                                       " class " + std::to_string(clsk[a] + 1) +
                                       " carries traffic with no identifiable source");
            arvH[a] = detail::snc_sum(parts);
        }
        for (std::size_t h = 0; h < here.size(); ++h) {
            const std::size_t a = here[h];
            std::vector<snc::Envelope> cross;
            for (std::size_t g = 0; g < here.size(); ++g)
                if (here[g] != a) cross.push_back(arvH[here[g]]);
            const double mua = mu[a];
            const snc::Envelope crossSum = detail::snc_sum(cross);
            const bool alone = cross.empty();
            srvH[a] = [mua, crossSum, alone](double theta) {
                const snc::Env s = snc::snc_srv_exp(mua, theta);
                if (alone) return s;
                return snc::snc_leftover(s, crossSum(theta));
            };
            const snc::Envelope arva = arvH[a], srva = srvH[a];
            (*outH)[a] = [arva, srva](double theta) {
                return snc::snc_output(arva(theta), srva(theta), theta);
            };
        }
    }

    for (std::size_t a = 0; a < nk; ++a) {
        const std::pair<std::size_t, std::size_t> key(statk[a], clsk[a]);
        env.arv[key] = arvH[a];
        env.srv[key] = srvH[a];
        env.lam[key] = lam[a];
        env.mu[key] = mu[a];
    }
    return env;
}

/**
 * @param L   the model
 * @param out the (Q,U,R,Tp,C,X) block to fill; shapes are set here
 */
template <class T, class Solution>
void solver_ba_snc(const qn::NetworkStruct<T>& L, Solution& out) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = L.nstations, K = L.nclasses;
    const SncEnvelopes env = solver_ba_snc_envelopes(L);

    out.Q = Matrix<T>(M, K, zero);
    out.U = Matrix<T>(M, K, zero);
    out.R = Matrix<T>(M, K, zero);
    out.Tp = Matrix<T>(M, K, zero);
    out.C.assign(K, zero);
    out.X.assign(K, zero);
    out.lG = std::numeric_limits<double>::quiet_NaN();
    out.iter = 1;

    for (std::map<std::pair<std::size_t, std::size_t>, snc::Envelope>::const_iterator it =
             env.arv.begin();
         it != env.arv.end(); ++it) {
        const std::size_t i = it->first.first, r = it->first.second;
        const double ed = snc::snc_mean_delay(it->second, env.srv.at(it->first)).value;
        out.R(i, r) = num_traits<T>::from_double(ed);
        out.Tp(i, r) = num_traits<T>::from_double(env.lam.at(it->first));
        out.U(i, r) = num_traits<T>::from_double(env.lam.at(it->first) / env.mu.at(it->first));
    }

    // ---- exact open-network quantities ----
    for (std::size_t i = 0; i < M; ++i) {
        if (L.stations[i].nodetype != qn::NodeType::Source) continue;
        for (std::size_t r = 0; r < K; ++r) {
            const T arr = L.rates(i, r);
            const double ad = num_traits<T>::to_double(arr);
            if (std::isfinite(ad) && ad > 0) {
                out.Tp(i, r) = T(out.Tp(i, r) + arr);
                out.X[r] = T(out.X[r] + arr);
            }
        }
    }
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r) out.Q(i, r) = T(out.Tp(i, r) * out.R(i, r));
    for (std::size_t r = 0; r < K; ++r) {
        if (out.X[r] > zero) {
            T sum = zero;
            for (std::size_t i = 0; i < M; ++i) sum = T(sum + out.Q(i, r));
            out.C[r] = T(sum / out.X[r]);
        }
    }
}

/** The response-time and queue-length QUANTILES of the 'snc' family. */
struct SncPercentiles {
    /** Response-time quantile per (station, class); NaN where no traffic. */
    Matrix<double> D;
    /** Queue-length quantile, in jobs, per (station, class); NaN where no traffic. */
    Matrix<double> B;
};

/**
 * Both quantile matrices from one envelope propagation.
 *
 * This is the native output of the family and has no counterpart in any other
 * BA family, which bound means only. `SolverBA.getDelayPerc` /
 * `getBacklogPerc` / `getPercTable` in MATLAB, the JAR and Python are the same
 * call.
 *
 * @param L   the model
 * @param eps violation probability, 0 < eps < 1
 */
template <class T>
SncPercentiles solver_ba_snc_perc(const qn::NetworkStruct<T>& L, double eps) {
    const SncEnvelopes env = solver_ba_snc_envelopes(L);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    SncPercentiles out;
    out.D = Matrix<double>(L.nstations, L.nclasses, nan);
    out.B = Matrix<double>(L.nstations, L.nclasses, nan);
    for (std::map<std::pair<std::size_t, std::size_t>, snc::Envelope>::const_iterator it =
             env.arv.begin();
         it != env.arv.end(); ++it) {
        const std::size_t i = it->first.first, r = it->first.second;
        out.D(i, r) = snc::snc_perc_delay(it->second, env.srv.at(it->first), eps).value;
        out.B(i, r) = snc::snc_perc_backlog(it->second, env.srv.at(it->first), eps).value;
    }
    return out;
}

}  // namespace ba
}  // namespace line

#endif  // LINE_SOLVERS_BA_SOLVER_BA_SNC_H
