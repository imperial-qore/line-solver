/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_BA_SOLVER_BA_QRF_ANALYZER_H
#define LINE_SOLVERS_BA_SOLVER_BA_QRF_ANALYZER_H

/**
 * Port of `matlab/src/solvers/BA/solver_ba_qrf_analyzer.m`, the adapter that
 * bridges the `sn` struct to the QRF (Quadratic Reduction Framework) bounds.
 *
 * WHAT THE ADAPTER IS FOR. `qrf_noblo_*` speaks in MAPs, phase counts and a
 * station-to-station routing matrix, not in a NetworkStruct; the whole of this
 * file is that translation plus the reconstruction of [Q,U,R,T,C,X] from the
 * utilizations the optimizer returns. The optimization itself lives in
 * `api/mapqn/mapqn_qrf_noblo.h` and is not repeated here.
 *
 * WHAT IT REFUSES AND WHY. The `qrf_noblo_*` formulation models every station as
 * ONE server and has no infinite-server notion, so a delay station reaches the
 * program as an unbounded direction and a c>1 station solved as c=1 is not a
 * bound in either direction. Both are refused by name, as the reference and the
 * Python twin do, rather than answered wrongly.
 *
 * THE THROUGHPUT COMES FROM A UTILIZATION, NOT FROM THE REFERENCE STATION.
 * `UN_qrf(i)` is P(n_i >= 1) marginalised over phase, so `U_i = X V_i S_i` is
 * exact at a single server and any loaded station determines X. Inverting
 * instead at the reference station assumes R == S there, which holds only for
 * an infinite server and returned an X above the bottleneck capacity.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mapqn/mapqn_qr_bounds_bas.h"
#include "line/api/mapqn/mapqn_qr_bounds_rsrd.h"
#include "line/api/mapqn/mapqn_qrf_bas_nlp.h"
#include "line/api/sn/sn_to_qrf_alpha.h"
#include "line/api/sn/sn_to_qrf_blocking.h"
#include "line/api/mapqn/mapqn_qrf_noblo.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/ba/solver_ba_analyzer.h"

namespace line {
namespace ba {

namespace detail {

/** True for one of the four QRF method names this port serves. */
inline bool is_qrf_noblo_method(const std::string& m) {
    return m == "qrf.mmi" || m == "qrf.mem" || m == "qrf.bethe" || m == "qrf.mmi.ld" ||
           m == "qrf.mmi.linear";
}

/** True for the two LP blocking bounds, which need `qrf_params`. */
inline bool is_qrf_lp_method(const std::string& m) { return m == "qrf.bas" || m == "qrf.rsrd"; }

/** True for one of the three nonlinear BAS-blocking bounds. */
inline bool is_qrf_bas_nlp_method(const std::string& m) {
    return m == "qrf.bas.mem" || m == "qrf.bas.bethe" || m == "qrf.bas.mmi";
}

/** mu/v as the per-station K(i) x K(i) blocks the LP bounds take. */
template <class T>
void qrf_blocks_from_maps(const std::vector<std::pair<Matrix<T>, Matrix<T> > >& MAPs,
                          const std::vector<int>& K, std::vector<Matrix<T> >* mu,
                          std::vector<Matrix<T> >* v) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = MAPs.size();
    mu->assign(M, Matrix<T>());
    v->assign(M, Matrix<T>());
    for (std::size_t i = 0; i < M; ++i) {
        const std::size_t k = static_cast<std::size_t>(K[i]);
        (*mu)[i] = Matrix<T>(k, k, zero);
        (*v)[i] = Matrix<T>(k, k, zero);
        for (std::size_t h = 0; h < k; ++h)
            for (std::size_t j = 0; j < k; ++j) {
                (*mu)[i](h, j) = MAPs[i].second(h, j);
                // (from, to), as in qrf_extract_mu_v; the diagonal of D0 is the
                // exit rate, not a background transition.
                (*v)[i](h, j) = (h == j) ? zero : MAPs[i].first(h, j);
            }
    }
}

/** Build the BAS parameters shared by the LP and entropy objectives. */
template <class T>
mapqn::QrBasParams<T> qrf_bas_params(
    const qn::NetworkStruct<T>& L, const BaOptions& opt,
    const std::vector<std::pair<Matrix<T>, Matrix<T> > >& MAPs, const std::vector<int>& K,
    const Matrix<T>& rt, std::size_t N, bool require_capacity) {
    const std::size_t M = L.nstations;
    (void)require_capacity;
    // The blocking tables are DERIVED from sn rather than demanded from the
    // caller: the model fixes every one of them. The refusal this replaces was
    // right only while the alternative was to INVENT them -- substituting no
    // blocking (MR = 1) measured 4.16667 from exact on
    // sanity_CQN_rm_{fcfs,ps}_1class where real tables sit at 0.133333, i.e.
    // 31x closer. options.config.qrf_params stays an explicit override.
    BaOptions::QrfParams derived;
    if (!opt.qrf_params.supplied) {
        int Ktot = 0;
        for (std::size_t i = 0; i < K.size(); ++i) Ktot += std::max(1, K[i]);
        const sn::QrfBlocking blk = sn::sn_to_qrf_blocking(L, Ktot);
        if (!blk.msg.empty())
            throw UnsupportedError("solver_ba_qrf_analyzer: the '" + opt.method +
                                   "' method cannot be applied to this model: " + blk.msg +
                                   " Supply options.config.qrf_params explicitly to override the "
                                   "derivation");
        derived.supplied = true;
        derived.f = blk.f;
        derived.MR = blk.MR;
        derived.BB = blk.BB;
        derived.MM = blk.MM;
        derived.MM1 = blk.MM1;
        derived.ZZ = blk.ZZ;
        derived.F = blk.F;
    }
    const BaOptions::QrfParams& qp = opt.qrf_params.supplied ? opt.qrf_params : derived;

    std::vector<int> F = qp.F;
    if (F.empty()) {
        // F is an OCCUPANCY BOUND, not a declared capacity: sn_to_qrf_capacity
        // decides binding through sn_get_buffer_size, which folds classcap and
        // the reachable population in as stations[i].cap alone does not.
        const sn::QrfCapacity cap = sn::sn_to_qrf_capacity(L);
        if (!cap.msg.empty())
            throw UnsupportedError("solver_ba_qrf_analyzer: the '" + opt.method +
                                   "' method cannot be applied: " + cap.msg);
        F = cap.F;
    }
    std::vector<Matrix<T> > mu, v;
    qrf_blocks_from_maps(MAPs, K, &mu, &v);

    mapqn::QrBasParams<T> bp;
    bp.M = static_cast<int>(M);
    bp.N = static_cast<int>(N);
    bp.F = F;
    bp.K = K;
    bp.mu = mu;
    bp.v = v;
    bp.r = rt;
    bp.MR = qp.MR;
    bp.BB = qp.BB;
    bp.ZZ = qp.ZZ;
    // qrf_params carries the reference's 1-based queue indices, with 0 for an
    // absent MM1 entry; the port indexes queues from 0 and marks absence -1.
    bp.f = qp.f - 1;
    bp.MM.assign(qp.MM.size(), std::vector<int>());
    for (std::size_t m = 0; m < qp.MM.size(); ++m)
        for (std::size_t c = 0; c < qp.MM[m].size(); ++c) bp.MM[m].push_back(qp.MM[m][c] - 1);
    bp.MM1.assign(qp.MM1.size(), std::vector<int>());
    for (std::size_t m = 0; m < qp.MM1.size(); ++m)
        for (std::size_t j = 0; j < qp.MM1[m].size(); ++j) bp.MM1[m].push_back(qp.MM1[m][j] - 1);
    bp.ZM = 0;
    for (std::size_t m = 0; m < bp.ZZ.size(); ++m) bp.ZM = std::max(bp.ZM, bp.ZZ[m]);
    return bp;
}

/** Run one LP blocking bound and return its utilization vector. */
template <class T>
std::vector<T> qrf_lp_utilizations(const qn::NetworkStruct<T>& L, const BaOptions& opt,
                                   const std::vector<std::pair<Matrix<T>, Matrix<T> > >& MAPs,
                                   const std::vector<int>& K, const Matrix<T>& rt,
                                   const Matrix<T>& alpha, std::size_t N) {
    const std::size_t M = L.nstations;
    if (opt.method == "qrf.rsrd") {
        // RS-RD carries NO blocking tables -- QrRsrdParams has no f/MR/BB/MM/
        // MM1/ZZ member at all -- so it never needed qrf_params, and demanding
        // them refused a well-formed call. Its PBB constraint sums over every
        // queue that can be FULL, so it also admits SEVERAL finite buffers
        // where 'qrf.bas' admits one. What it does need is a truthful F.
        std::vector<int> F = opt.qrf_params.F;
        if (F.empty()) {
            const sn::QrfCapacity cap = sn::sn_to_qrf_capacity(L);
            if (!cap.msg.empty())
                throw UnsupportedError("solver_ba_qrf_analyzer: the 'qrf.rsrd' method cannot be "
                                       "applied: " + cap.msg);
            F = cap.F;
        }
        std::vector<Matrix<T> > mu, v;
        qrf_blocks_from_maps(MAPs, K, &mu, &v);
        mapqn::QrRsrdParams<T> rp;
        rp.M = static_cast<int>(M);
        rp.N = static_cast<int>(N);
        rp.F = F;
        rp.K = K;
        rp.mu = mu;
        rp.v = v;
        rp.r = rt;
        if (alpha.rows() == M && alpha.cols() > 0) {
            rp.alpha.assign(M, std::vector<T>());
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t n = 0; n < alpha.cols(); ++n) rp.alpha[i].push_back(alpha(i, n));
        }
        const mapqn::QrRsrdResult<T> res = mapqn::mapqn_qr_bounds_rsrd(rp, 0, mapqn::MapqnSense::Max);
        if (!res.ok)
            throw NumericError("solver_ba_qrf_analyzer: the 'qrf.rsrd' linear program did not "
                               "solve (" + res.status + ")");
        return res.U;
    }

    const mapqn::QrBasParams<T> bp = qrf_bas_params(L, opt, MAPs, K, rt, N, false);
    const mapqn::QrBasResult<T> res = mapqn::mapqn_qr_bounds_bas(bp, 0, mapqn::MapqnSense::Max);
    if (!res.ok)
        throw NumericError("solver_ba_qrf_analyzer: the 'qrf.bas' linear program did not solve (" +
                           res.status + ")");
    return res.U;
}

/**
 * `derive_qn_from_bounds`: Little's law at the throughput the utilizations imply.
 *
 * The M/G/1-like Q = T S / (1 - U) is the reference's estimate; a station at or
 * above capacity is held at the population instead.
 */
template <class T>
std::vector<T> qrf_qn_from_bounds(const qn::NetworkStruct<T>& L, const std::vector<T>& UN,
                                  const std::vector<T>& V, const std::vector<T>& stimes,
                                  std::size_t N) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations;
    T X = zero;
    for (std::size_t i = 0; i < M; ++i) {
        const T c = num_traits<T>::from_double(L.stations[i].nservers);
        if (stimes[i] > zero && V[i] > zero && UN[i] > zero) {
            X = T(UN[i] * c / T(V[i] * stimes[i]));
            break;
        }
    }
    std::vector<T> QN(M, zero);
    for (std::size_t i = 0; i < M; ++i) {
        if (!(stimes[i] > zero)) continue;
        const T Ti = T(X * V[i]);
        QN[i] = (UN[i] < one) ? T(T(Ti * stimes[i]) / T(one - UN[i]))
                              : num_traits<T>::from_int(static_cast<long>(N));
    }
    return QN;
}

}  // namespace detail

/**
 * Port of `solver_ba_qrf_analyzer`.
 *
 * @param L   the refreshed struct of a single-class closed model
 * @param opt the method and, for the load-dependent arms, `qrf_alpha`
 * @return the [Q,U,R,T,C,X] implied by the QRF utilization bound
 */
template <class T>
BaSolution<T> solver_ba_qrf_analyzer(const qn::NetworkStruct<T>& L, const BaOptions& opt) {
    const T zero = num_traits<T>::from_int(0);
    const std::string& method = opt.method;
    const std::size_t M = L.nstations;

    if (L.nclasses != 1)
        throw UnsupportedError("solver_ba_qrf_analyzer: QRF methods support single-class networks "
                               "only");
    for (double n : L.njobs())
        if (std::isinf(n))
            throw UnsupportedError("solver_ba_qrf_analyzer: QRF methods support closed networks "
                                   "only");
    const std::size_t N = static_cast<std::size_t>(L.nclosedjobs());
    if (N == 0)
        throw UnsupportedError("solver_ba_qrf_analyzer: QRF methods support closed networks only");
    // THE LOAD-DEPENDENT ARMS SERVE DELAY, MULTISERVER AND LOAD-DEPENDENT
    // STATIONS; THE REST STILL CANNOT. alpha(i,n) multiplies every rate out of
    // station i at population n, which IS the rate law of a delay (alpha = n),
    // of a c-server station (alpha = min(n,c)) and of limited load dependence,
    // so `qrf.mmi.ld` and `qrf.mmi.linear` answer the model's OWN chain on all
    // three rather than an approximation of it. `sn_to_qrf_alpha` derives alpha
    // and owns the one restriction that survives: a station serving several
    // jobs at once must be exponential, because the QRF local state carries one
    // phase per station and that describes one job in service and no more.
    //
    // Every other arm builds a population-free q, so it models each station as
    // one server and a c>1 station solved as c=1 is not a bound in either
    // direction. They keep refusing, by naming the two arms that do serve the
    // model. Same gate as the MATLAB, python and JAR twins.
    const sn::QrfAlpha alpha_sn = sn::sn_to_qrf_alpha(L);
    const bool method_is_ld = (method == "qrf.mmi.ld" || method == "qrf.mmi.linear");
    if (alpha_sn.ld && !method_is_ld)
        throw UnsupportedError(
            "solver_ba_qrf_analyzer: the '" + method +
            "' method models every station as a single server: its transition rates carry no "
            "population index, so it has nowhere to put the rate of a delay, a multiserver or a "
            "load-dependent station. Use 'qrf.mmi.ld' or 'qrf.mmi.linear', which do");
    if (!alpha_sn.msg.empty())
        throw UnsupportedError("solver_ba_qrf_analyzer: the '" + method +
                               "' method cannot be applied: " + alpha_sn.msg);

    // MAPs{i} = {D0, D1}; a station with no service process falls back to the
    // unit-rate exponential, as the reference does.
    std::vector<std::pair<Matrix<T>, Matrix<T>>> MAPs(M);
    std::vector<int> K_phases(M, 1);
    for (std::size_t i = 0; i < M; ++i) {
        const mam::Map<T> m = lang::dist_to_map(L.service[i][0]);
        if (m.D0.rows() == 0) {
            Matrix<T> D0(1, 1, num_traits<T>::from_int(-1)), D1(1, 1, num_traits<T>::from_int(1));
            MAPs[i] = std::make_pair(D0, D1);
            K_phases[i] = 1;
        } else {
            MAPs[i] = std::make_pair(m.D0, m.D1);
            K_phases[i] = static_cast<int>(m.D0.rows());
        }
    }

    // Single class, so the (MK x MK) routing matrix is already station indexed.
    Matrix<T> rt(M, M, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < M; ++j) rt(i, j) = L.rt(i, j);

    Matrix<T> mu, v;
    mapqn::qrf_extract_mu_v(MAPs, M, K_phases, &mu, &v);

    // options.config.qrf_alpha overrides the derivation, as options.config
    // .qrf_params does for the blocking tables; absent it, alpha comes from the
    // model, which is all ones on a load-independent one.
    Matrix<T> alpha(M, N, zero);
    if (opt.qrf_alpha.rows() > 0) {
        if (opt.qrf_alpha.rows() != M || opt.qrf_alpha.cols() != N)
            throw InputError("solver_ba_qrf_analyzer: qrf_alpha must be (nstations x N)");
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t n = 0; n < N; ++n)
                alpha(i, n) = num_traits<T>::from_double(opt.qrf_alpha(i, n));
    } else {
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t n = 0; n < N; ++n)
                alpha(i, n) = num_traits<T>::from_double(alpha_sn.alpha[i][n]);
    }

    mapqn::QrfMetrics<T> r;
    if (detail::is_qrf_lp_method(method)) {
        r.UN = detail::qrf_lp_utilizations(L, opt, MAPs, K_phases, rt, alpha, N);
    } else if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "solver_ba_qrf_analyzer: the QRF objective is the MEM entropy, the mutual "
            "information or the tree-reweighted free entropy, all of which take a log, and "
            "exact arithmetic has no representation for it; use the double or real backend");
    } else {
        if (method == "qrf.mmi")
            r = mapqn::qrf_noblo_mmi(M, K_phases, N, mu, v, rt);
        else if (method == "qrf.mem")
            r = mapqn::qrf_noblo_mem(MAPs, N, rt);
        else if (method == "qrf.bethe")
            // Same polytope and same phase-1 start as qrf.mmi; the objective is
            // the tree-reweighted (Bethe) free entropy at the uniform
            // spanning-tree weight lambda = 1/M, the largest uniform weight at
            // which the program is convex.
            r = mapqn::qrf_noblo_bethe(M, K_phases, N, mu, v, rt);
        else if (method == "qrf.mmi.ld")
            r = mapqn::qrf_noblo_mmi_ld(M, K_phases, N, mu, v, rt, alpha);
        else if (method == "qrf.mmi.linear")
            r = mapqn::qrf_noblo_mmi_linear(MAPs, N, rt, alpha);
        else if (method == "qrf.bas.mmi")
            r = mapqn::mapqn_qrf_bas_mmi(
                detail::qrf_bas_params(L, opt, MAPs, K_phases, rt, N, true));
        else if (method == "qrf.bas.bethe")
            // qrf.bethe's objective over the BAS decision vector: same tables
            // and same polytope as qrf.bas and qrf.bas.mem, lambda = 1/M.
            r = mapqn::mapqn_qrf_bas_bethe(
                detail::qrf_bas_params(L, opt, MAPs, K_phases, rt, N, true));
        else if (detail::is_qrf_bas_nlp_method(method))
            r = mapqn::mapqn_qrf_bas_mem(
                detail::qrf_bas_params(L, opt, MAPs, K_phases, rt, N, true));
        else
            throw UnsupportedError("solver_ba_qrf_analyzer: unknown QRF method '" + method + "'");
    }

    // cellsum(sn.visits): the visits summed over the chains, at station level
    std::vector<T> V(M, zero);
    for (std::size_t c = 0; c < L.nchains; ++c)
        for (std::size_t i = 0; i < M; ++i)
            V[i] += L.visits[c](L.stateful_of_station(i + 1) - 1, 0);

    std::vector<T> stimes(M, zero);
    for (std::size_t i = 0; i < M; ++i) {
        mam::Map<T> m;
        m.D0 = MAPs[i].first;
        m.D1 = MAPs[i].second;
        stimes[i] = mam::map_mean(m);
    }

    // `derive_qn_from_bounds`: an LP bound returns utilizations only, so the
    // queue lengths come from Little's law at the throughput those imply.
    if (r.QN.empty()) r.QN = detail::qrf_qn_from_bounds(L, r.UN, V, stimes, N);

    // Normalize the queue lengths to the population constraint.
    T qsum = zero;
    for (const T& q : r.QN) qsum += q;
    std::vector<T> QN = r.QN;
    if (qsum > zero)
        for (T& q : QN) q = T(q / qsum * num_traits<T>::from_int(static_cast<long>(N)));

    // The alpha-free arms and the LP ones leave BN empty, and there
    // BN = P(n >= 1) = UN: a single server's departure rate is proportional to
    // the probability that it is busy. Setting it here rather than in each arm
    // keeps the readout below one formula.
    if (r.BN.empty()) r.BN = r.UN;

    // X from the ALPHA-WEIGHTED marginal mean BN, the mean number of jobs
    // actually in service: E[min(n,c)] at a c-server station, E[n] at a delay,
    // P(n >= 1) at a single server. That is what the departure rate is
    // proportional to, so T_i = BN_i / stime_i holds exactly at the relaxed
    // point and X = T_i / V_i. The single-server case is the former
    // UN[i]*c/(V*stime) unchanged, c being 1 and BN being UN there; a delay no
    // longer needs the refstat fallback, since alpha = n makes BN = E[n] = QN,
    // which is what that fallback computed.
    T X = zero;
    for (std::size_t i = 0; i < M; ++i) {
        if (stimes[i] > zero && V[i] > zero && r.BN[i] > zero) {
            X = T(r.BN[i] / T(V[i] * stimes[i]));
            break;
        }
    }
    if (!(X > zero)) {
        const std::size_t rs = L.classes[0].refstat;
        if (rs >= 1 && rs <= M && stimes[rs - 1] > zero && V[rs - 1] > zero)
            X = T(QN[rs - 1] / T(V[rs - 1] * stimes[rs - 1]));
    }

    BaSolution<T> s;
    s.Q = Matrix<T>(M, 1, zero);
    s.U = Matrix<T>(M, 1, zero);
    s.R = Matrix<T>(M, 1, zero);
    s.Tp = Matrix<T>(M, 1, zero);
    for (std::size_t i = 0; i < M; ++i) {
        s.Q(i, 0) = QN[i];
        s.Tp(i, 0) = T(X * V[i]);
        if (std::isinf(alpha_sn.peak[i]))
            // Delay (infinite server): U = QN by LINE convention.
            s.U(i, 0) = QN[i];
        else
            // The busy fraction of the station's DECLARED peak capacity, BN
            // being the mean number of jobs in service. The normalizer is
            // nservers times the reachable lld peak, LINE's one U = T*S/peak
            // convention (see sn_to_qrf_alpha); at peak 1 it is r.UN[i], what
            // the alpha-free arms report directly.
            s.U(i, 0) = T(r.BN[i] / num_traits<T>::from_double(alpha_sn.peak[i]));
        if (s.Tp(i, 0) > zero) s.R(i, 0) = T(s.Q(i, 0) / s.Tp(i, 0));
    }
    s.X.assign(1, X);
    s.C.assign(1, X > zero ? T(num_traits<T>::from_int(static_cast<long>(N)) / X) : zero);
    s.lG = std::numeric_limits<double>::quiet_NaN();
    return s;
}

}  // namespace ba
}  // namespace line

#endif  // LINE_SOLVERS_BA_SOLVER_BA_QRF_ANALYZER_H
