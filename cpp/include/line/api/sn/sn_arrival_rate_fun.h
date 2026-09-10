/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_ARRIVAL_RATE_FUN_H
#define LINE_API_SN_SN_ARRIVAL_RATE_FUN_H

/**
 * Port of matlab/src/api/sn/sn_arrival_rate_fun.m.
 *
 * The arrival rate of a station-class pair AS A FUNCTION OF TIME. The
 * time-varying analyses (Mt/G/inf, the modified offered load, the Gt/Mt/st+GI
 * fluid queue) consume lambda(t) itself, not a mean rate: their whole content
 * is the LAG between when work arrives and when it is felt, and a time-averaged
 * rate has no lag. LINE carries a time-varying arrival as an NHPP or a MAPt,
 * whose schedule is piecewise constant, so lambda(t) is read off the segment in
 * force at t.
 *
 * For any other process the rate is constant and the handle returns it, which
 * is what lets a caller ask for the time-varying analysis of a stationary model
 * and get the stationary answer rather than an error.
 *
 * ARITHMETIC: field. map_pie is a linear solve; no transcendental appears.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"

namespace line {
namespace api {

/** lambda(t), with whether it actually varies and the cycle length. */
template <class T>
struct ArrivalRateFun {
    std::function<T(const T&)> lambda;  ///< the rate as a function of time
    bool timeVarying = false;           ///< whether it depends on t at all
    double period = std::numeric_limits<double>::infinity();  ///< cycle length when cyclic
};

/**
 * Build lambda(t) for station `ist` (0-based), class `r`.
 */
template <class T>
ArrivalRateFun<T> sn_arrival_rate_fun(const qn::NetworkStruct<T>& sn, std::size_t ist,
                                      std::size_t r) {
    ArrivalRateFun<T> out;
    const T rate = sn.rates(ist, r);
    const lang::Distrib<T>& d = sn.service[ist][r];
    if (!d.has_schedule()) {
        out.lambda = [rate](const T&) { return rate; };
        return out;
    }
    // The schedule is ALREADY in MAP form here: `sched_D0[k]`, `sched_D1[k]` are
    // the pair of segment k for a MAPt and for a PHt alike, the PHt having been
    // converted on the way in. An NHPP is the one-phase case, whose D1 entry IS
    // its lambda.
    const std::vector<T> bp = d.sched_bp;
    const std::size_t n = d.sched_D0.size();
    std::vector<T> seg(n, num_traits<T>::from_int(0));
    for (std::size_t k = 0; k < n; ++k) {
        const Matrix<T>& D0 = d.sched_D0[k];
        const Matrix<T>& D1 = d.sched_D1[k];
        if (D0.rows() == 1) {
            seg[k] = D1(0, 0);
        } else {
            // The arrival rate of a segment is pie_k D1_k e, the stationary
            // throughput of that segment's own MAP.
            mam::Map<T> m;
            m.D0 = D0;
            m.D1 = D1;
            const std::vector<T> pie = mam::map_pie(m);
            T v = num_traits<T>::from_int(0);
            for (std::size_t i = 0; i < D1.rows(); ++i)
                for (std::size_t j = 0; j < D1.cols(); ++j) v = T(v + pie[i] * D1(i, j));
            seg[k] = v;
        }
    }
    for (std::size_t k = 1; k < n; ++k)
        if (std::fabs(num_traits<T>::to_double(seg[k]) - num_traits<T>::to_double(seg[0])) > 1e-12)
            out.timeVarying = true;
    const bool cyclic = d.sched_cyclic;
    if (cyclic && bp.size() >= 2)
        out.period = num_traits<T>::to_double(bp.back()) - num_traits<T>::to_double(bp.front());
    out.lambda = [bp, seg, cyclic](const T& t) {
        T u = t;
        if (cyclic && bp.size() >= 2 && bp.back() > bp.front()) {
            const double span = num_traits<T>::to_double(bp.back() - bp.front());
            double x = num_traits<T>::to_double(t) - num_traits<T>::to_double(bp.front());
            x = std::fmod(std::fmod(x, span) + span, span);
            u = T(bp.front() + num_traits<T>::from_double(x));
        }
        // Segment k is in force on [bp[k], bp[k+1]). Before the first breakpoint
        // the first segment holds and after the last the last one does, so a
        // caller integrating over an infinite past (the Mt/G/inf convolution)
        // gets a defined rate everywhere rather than a NaN.
        std::size_t k = 0;
        for (std::size_t i = 0; i + 1 < bp.size(); ++i)
            if (u >= bp[i]) k = i;
        if (k >= seg.size()) k = seg.size() - 1;
        return seg[k];
    };
    return out;
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_ARRIVAL_RATE_FUN_H
