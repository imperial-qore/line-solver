/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_HAS_BURSTY_ARRIVAL_H
#define LINE_API_SN_SN_HAS_BURSTY_ARRIVAL_H

/**
 * Port of matlab/src/api/sn/sn_has_bursty_arrival.m.
 *
 * True iff some external (Source) arrival process is non-renewal (bursty). A MAP
 * (D0, D1) is renewal iff D1 equals its rank-one renewal form t0*pie, with
 * t0 = -D0 e = D1 e and pie the embedded stationary vector; any departure from
 * that form signals correlation between successive inter-arrival times. A
 * single-phase arrival (Poisson) is renewal by construction. This is the exact
 * test SolverMVA.resolveMethod uses to upgrade method='default' to 'rqna'.
 *
 * ARITHMETIC: field. map_pie is a linear solve; no transcendental appears.
 */

#include <cmath>
#include <cstddef>

#include "line/api/mam/map_moment.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"

namespace line {
namespace api {

template <class T>
bool sn_has_bursty_arrival(const qn::NetworkStruct<T>& sn) {
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        if (sn.stations[i].nodetype != qn::NodeType::Source) continue;
        for (std::size_t r = 0; r < sn.nclasses; ++r) {
            if (!sn.disabled.empty() && sn.disabled[i][r]) continue;
            const lang::Distrib<T>& d = sn.service[i][r];
            if (d.disabled || !d.has_map()) continue;
            // Read the (D0, D1) pair directly rather than via dist_to_map: a
            // non-MAP arrival is renewal (skipped above), and dist_to_map's
            // Replayer branch instantiates the transcendental aph_fit, which
            // would forbid this field-arithmetic predicate under Rational.
            mam::Map<T> m;
            m.D0 = d.D0;
            m.D1 = d.D1;
            const std::size_t n = m.D1.rows();
            if (n <= 1) continue;  // single-phase arrival is Poisson, hence renewal
            const std::vector<T> pie = mam::map_pie(m);
            // t0 = D1 e (row sums of D1), renewal form D1ren = t0 * pie.
            std::vector<T> t0(n, num_traits<T>::from_int(0));
            for (std::size_t a = 0; a < n; ++a) {
                T acc = num_traits<T>::from_int(0);
                for (std::size_t b = 0; b < n; ++b) acc = T(acc + m.D1(a, b));
                t0[a] = acc;
            }
            double diff2 = 0.0, norm2 = 0.0;
            for (std::size_t a = 0; a < n; ++a)
                for (std::size_t b = 0; b < n; ++b) {
                    const double dd = num_traits<T>::to_double(m.D1(a, b));
                    const double rr = num_traits<T>::to_double(T(t0[a] * pie[b]));
                    diff2 += (dd - rr) * (dd - rr);
                    norm2 += dd * dd;
                }
            if (std::sqrt(diff2) > 1e-8 * std::max(1.0, std::sqrt(norm2))) return true;
        }
    }
    return false;
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_HAS_BURSTY_ARRIVAL_H
