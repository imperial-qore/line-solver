/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MAM_SOLVER_MAM_PASSAGE_TIME_H
#define LINE_SOLVERS_MAM_SOLVER_MAM_PASSAGE_TIME_H

/**
 * Port of `solver_mam_passage_time.m`: the response-time (sojourn-time)
 * distribution of a single open queue, which is what `getCdfRespT`,
 * `getSjrnT` / `sjrnT` and the CDF path of `getPerctRespT` return.
 *
 * THE REGIME IS NARROW AND THE REFERENCE SAYS SO. The whole analyzer is inside
 * `if M == 2 && all(isinf(N))`: exactly two stations, a Source and one queue,
 * every class open. Anything else makes MATLAB warn and return with NO result
 * at all -- an empty `RD` that then surfaces as a confusing failure further up.
 * The port refuses by name instead, which is the same information delivered
 * where it can be acted on.
 *
 * TWO ENGINES, chosen by the queue's discipline:
 *
 *  - FCFS / HOL: the sojourn time is a phase-type law read straight out of the
 *    age process (`mmapph1fcfs_stdistr_ph`). The evaluation grid is the
 *    reference's: start at mean + 5 sigma and widen until the CDF is within
 *    FineTol of one, then lay down `num_cdf_pts` points from zero.
 *  - PS: the MAP/M/1-PS sojourn law of `map_m1ps_cdfrespt`, which returns the
 *    COMPLEMENTARY CDF, so the port takes 1 - W_bar as the reference does. The
 *    grid there is 10x the M/M/1-PS mean, and the reference requires
 *    exponential service (order-1 subgenerator) and, for several classes,
 *    identical service rates -- both refused by name here.
 *
 * WHAT IS NOT PORTED, and refuses: the priority branch. Distinct priorities
 * under HOL reach BUTools' `MMAPPH1NPPR` (`MMAPPH1PRPR` when preemptive),
 * whose sojourn law MATLAB, the JAR and python tabulate through 'stMoms' +
 * 'stDistr' since 2026-08-27 (no vendored analyzer exports a PH form); neither
 * analyzer is ported to C++, exactly as in `solver_mam_basic`, so this arm is
 * the one member of the family still refusing it.
 */

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mam/map_cdf.h"
#include "line/api/mam/map_m1ps.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mam/mmap_assemble.h"
#include "line/api/mam/mmapph1fcfs.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mam/mam_types.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** One class's response-time CDF, the reference's `RD{station, class} = [F, X]`. */
template <class T>
struct RespTCdf {
    std::vector<T> F;  ///< CDF values
    std::vector<T> X;  ///< the points they are evaluated at
};

/**
 * Port of `solver_mam_passage_time.m`.
 *
 * @return one entry per class, in class order; the Source contributes none
 *         (the reference leaves `RD{idx_arv,k}` empty)
 */
template <class T>
std::vector<RespTCdf<T>> solver_mam_passage_time(const qn::NetworkStruct<T>& L,
                                                 const MamOptions& opt) {
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "solver_mam_passage_time: the sojourn-time law comes from the age process, whose "
            "first-return matrix is a tolerance-terminated Riccati doubling; rerun with "
            "--arith double or --arith real");
    } else {
    using lang::GlobalConstants;
    using lang::SchedStrategy;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations, K = L.nclasses;

    bool allopen = true;
    for (const qn::JobClass& c : L.classes)
        if (!std::isinf(c.population)) allopen = false;
    if (M != 2 || !allopen)
        throw UnsupportedError(
            "solver_mam_passage_time: the MAM response-time distribution covers a single open "
            "queue only (exactly two stations, a Source and one queue, every class open); this "
            "model has " + std::to_string(M) +
            " stations. The reference warns and returns no result at all for it");

    std::size_t src = 0, q = 0;
    for (std::size_t i = 1; i <= M; ++i) {
        if (L.stations[i - 1].sched == SchedStrategy::EXT) src = i;
        else q = i;
    }
    if (src == 0 || q == 0)
        throw UnsupportedError(
            "solver_mam_passage_time: the model must have a Source and one queueing station");

    const SchedStrategy qs = L.stations[q - 1].sched;
    if (!(qs == SchedStrategy::FCFS || qs == SchedStrategy::HOL || qs == SchedStrategy::PS))
        throw UnsupportedError(std::string("solver_mam_passage_time: the ") +
                               lang::sched_to_text(qs) +
                               " discipline is not covered; the reference supports FCFS, HOL and "
                               "PS only");
    // Priorities select the law only under a priority DISCIPLINE (HOL here) --
    // a plain FCFS or PS queue serves in arrival or processor order whatever
    // the prio column says, exactly as the other three codebases gate it
    if (qs == SchedStrategy::HOL) {
        bool distinct = false;
        for (std::size_t r = 1; r < K; ++r)
            if (L.classes[r].prio != L.classes[0].prio) distinct = true;
        if (distinct)
            throw UnsupportedError(
                "solver_mam_passage_time: non-identical class priorities under HOL route to "
                "BUTools' MMAPPH1NPPR (MMAPPH1PRPR when preemptive), whose tabulated sojourn law "
                "MATLAB, the JAR and python now serve; neither priority analyzer is ported to "
                "C++, so this arm still refuses by name");
    }

    const std::size_t npts =
        opt.num_cdf_pts > 0 ? opt.num_cdf_pts : static_cast<std::size_t>(100);

    // The arrival MMAP: each class's source process marked as its own class.
    Mmap<T> A;
    {
        std::vector<Mmap<T>> parts;
        for (std::size_t r = 0; r < K; ++r) {
            Mmap<T> m;
            const Map<T> s = lang::dist_to_map(L.service[src - 1][r]);
            m.D0 = s.D0;
            m.D1 = s.D1;
            m.Dc.assign(1, s.D1);
            parts.push_back(m);
        }
        A = parts[0];
        for (std::size_t p = 1; p < parts.size(); ++p)
            A = mmap_super(A, parts[p]);
    }

    std::vector<RespTCdf<T>> out(K);

    if (qs == SchedStrategy::PS) {
        // MAP/M/1-PS: the reference requires exponential service and, with
        // several classes, one shared rate.
        std::vector<T> mu(K, zero);
        for (std::size_t r = 0; r < K; ++r) {
            const Map<T> sv = lang::dist_to_map(L.service[q - 1][r]);
            if (sv.D0.rows() != 1)
                throw UnsupportedError(
                    "solver_mam_passage_time: a PS queue requires exponential (order-1) service "
                    "times; the MAP/M/1-PS sojourn law has no phase-type generalization here");
            mu[r] = T(-sv.D0(0, 0));
        }
        for (std::size_t r = 1; r < K; ++r)
            if (std::fabs(num_traits<T>::to_double(mu[r]) - num_traits<T>::to_double(mu[0])) >
                GlobalConstants::FineTol)
                throw UnsupportedError(
                    "solver_mam_passage_time: multi-class PS requires identical service rates");

        Matrix<T> Dagg(A.order(), A.order(), zero);
        for (std::size_t c = 0; c < A.classes(); ++c)
            for (std::size_t i = 0; i < Dagg.rows(); ++i)
                for (std::size_t j = 0; j < Dagg.cols(); ++j) Dagg(i, j) += A.Dc[c](i, j);
        const T lambda = map_lambda(Map<T>{A.D0, Dagg});
        const T rho = T(lambda / mu[0]);
        if (!(num_traits<T>::to_double(rho) < 1.0))
            throw NumericError(
                "solver_mam_passage_time: the PS queue is unstable, so it has no sojourn-time "
                "distribution");
        const T mean = T(one / T(mu[0] * T(one - rho)));
        const T xmax = T(mean * num_traits<T>::from_int(10));
        std::vector<T> x(npts, zero);
        for (std::size_t i = 0; i < npts; ++i)
            x[i] = (npts == 1) ? xmax
                               : T(xmax * num_traits<T>::from_double(
                                              static_cast<double>(i) /
                                              static_cast<double>(npts - 1)));
        const MapM1psResult<T> res = map_m1ps_cdfrespt(A.D0, Dagg, mu[0], x);
        for (std::size_t r = 0; r < K; ++r) {
            out[r].X = x;
            out[r].F.resize(npts);
            for (std::size_t i = 0; i < npts; ++i) out[r].F[i] = T(one - res.w_bar[i]);
        }
        return out;
    }

    // FCFS / HOL: the sojourn law is phase-type, straight out of the age process.
    std::vector<PhService<T>> svc(K);
    for (std::size_t r = 0; r < K; ++r) {
        const double ns = L.stations[q - 1].nservers;
        const Map<T> raw = lang::dist_to_map(L.service[q - 1][r]);
        const T target = std::isfinite(ns)
                             ? T(map_mean(raw) / num_traits<T>::from_double(ns))
                             : map_mean(raw);
        const Map<T> sc = map_scale(raw, target);
        svc[r].sigma = map_pie(sc);
        svc[r].S = sc.D0;
    }
    const std::vector<StDistrPh<T>> ph = mmapph1fcfs_stdistr_ph(A, svc);

    for (std::size_t r = 0; r < K; ++r) {
        // Read the PH pair back as the MAP {A, (-A e) alpha} the moment and CDF
        // routines take.
        const std::size_t n = ph[r].A.rows();
        Matrix<T> D1(n, n, zero);
        for (std::size_t i = 0; i < n; ++i) {
            T s = zero;
            for (std::size_t j = 0; j < n; ++j) s += ph[r].A(i, j);
            for (std::size_t j = 0; j < n; ++j) D1(i, j) = T(-s * ph[r].alpha[j]);
        }
        const Map<T> RDph{ph[r].A, D1};
        const T mean = map_mean(RDph);
        const T sigma = num_traits<T>::from_double(
            std::sqrt(num_traits<T>::to_double(map_var(RDph))));
        // Widen until the tail beyond the grid is negligible, as the reference does.
        int nsig = 5;
        for (;;) {
            const T probe = T(mean + num_traits<T>::from_int(nsig) * sigma);
            const std::vector<T> c = map_cdf(RDph, std::vector<T>{probe});
            if (num_traits<T>::to_double(c[0]) >= 1.0 - GlobalConstants::FineTol) break;
            if (++nsig > 1000)
                throw NumericError(
                    "solver_mam_passage_time: the sojourn-time CDF did not reach one within 1000 "
                    "standard deviations");
        }
        const T xmax = T(mean + num_traits<T>::from_int(nsig) * sigma);
        std::vector<T> x(npts, zero);
        for (std::size_t i = 0; i < npts; ++i)
            x[i] = (npts == 1) ? xmax
                               : T(xmax * num_traits<T>::from_double(
                                              static_cast<double>(i) /
                                              static_cast<double>(npts - 1)));
        out[r].X = x;
        out[r].F = map_cdf(RDph, x);
    }
    return out;
    }  // if constexpr has_transcendental
}

/**
 * Port of the CDF path of `@@SolverMAM/getPerctRespT.m`: linear interpolation of
 * the response-time CDF at the requested percentile levels.
 *
 * Duplicate CDF values are collapsed keeping the LAST, as MATLAB's
 * `unique(probs,'last')` does, so a flat tail interpolates from the largest
 * time carrying that probability rather than the smallest.
 *
 * The FJ_codes branch of the reference is not reachable here: it reads
 * percentiles stored by `solver_mam_fj`, which is not ported and whose models
 * the dispatch refuses.
 */
template <class T>
std::vector<T> mam_percentiles_from_cdf(const RespTCdf<T>& cdf, const std::vector<double>& pcts) {
    const T zero = num_traits<T>::from_int(0);
    if (cdf.F.empty()) throw InputError("mam_percentiles_from_cdf: the CDF is empty");
    std::vector<double> p, t;
    for (std::size_t i = 0; i < cdf.F.size(); ++i) {
        const double fi = num_traits<T>::to_double(cdf.F[i]);
        if (!p.empty() && std::fabs(fi - p.back()) < 1e-15) {
            t.back() = num_traits<T>::to_double(cdf.X[i]);  // keep the LAST
            continue;
        }
        p.push_back(fi);
        t.push_back(num_traits<T>::to_double(cdf.X[i]));
    }
    std::vector<T> out(pcts.size(), zero);
    for (std::size_t k = 0; k < pcts.size(); ++k) {
        // Percentiles above 1 are read as percentages, as the reference does.
        const double want = pcts[k] > 1.0 ? pcts[k] / 100.0 : pcts[k];
        if (p.size() == 1) {
            out[k] = num_traits<T>::from_double(t[0]);
            continue;
        }
        std::size_t lo = 0;
        while (lo + 2 < p.size() && p[lo + 1] < want) ++lo;
        const double denom = p[lo + 1] - p[lo];
        const double frac = (std::fabs(denom) < 1e-300) ? 0.0 : (want - p[lo]) / denom;
        out[k] = num_traits<T>::from_double(t[lo] + frac * (t[lo + 1] - t[lo]));
    }
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_SOLVERS_MAM_SOLVER_MAM_PASSAGE_TIME_H
