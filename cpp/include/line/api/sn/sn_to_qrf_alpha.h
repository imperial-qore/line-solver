/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_TO_QRF_ALPHA_H
#define LINE_API_SN_SN_TO_QRF_ALPHA_H

/**
 * The QRF load-dependent rate scaling alpha(i,n), derived from an sn.
 *
 * Port of matlab/src/api/sn/sn_to_qrf_alpha.m.
 *
 * The load-dependent QRF arms carry a scaling alpha(i,n) that multiplies EVERY
 * rate out of station i while it holds n jobs, completions mu and background
 * phase changes v alike -- see the q construction in `qrf_noblo_mmi_ld`. That
 * is exactly the rate law of
 *
 *   an infinite server        alpha(i,n) = n
 *   a c-server station        alpha(i,n) = min(n, c_i)
 *   limited load dependence   alpha(i,n) = sn.lldscaling(i,n)
 *
 * so the three COMPOSE BY MULTIPLICATION and not one of them is an
 * approximation: the relaxed chain is the model's own, and the QRF answer keeps
 * whatever status it had on a single-server model.
 *
 * WHERE IT STOPS BEING THE MODEL'S OWN IS PHASE-TYPE SERVICE AT A STATION THAT
 * SERVES SEVERAL JOBS AT ONCE. The QRF local state carries ONE phase per
 * station, a faithful description of one job in service and of nothing else:
 * min(n,c) jobs served in parallel each advance through a phase of their own,
 * and no scaling of a single-phase process reproduces that joint motion. A
 * multiserver or delay station must therefore be exponential -- scaling a PH
 * server by min(n,c) would answer a DIFFERENT chain, so the relaxation would
 * stop containing the model's stationary distribution and the number would
 * bound nothing. Limited load dependence at a SINGLE server is exempt and
 * admits PH freely: one job is in service whatever the rate.
 *
 * THE UTILIZATION NORMALIZER IS THE DECLARED PEAK, NOT max(alpha). LINE reports
 * U = T*S/peak at every station whose rate scales with the population, one
 * convention shared by multiserver, lld and class dependence. `peak` is
 * therefore nservers(i) times the largest lld scaling the model can REACH, and
 * not max(alpha(i,:)): at c = 3 with N = 2 the reachable alpha peaks at 2 while
 * the station still has three servers, and normalizing by 2 would report a
 * utilization the model never attains. Infinite at a delay, where LINE reports
 * U = QN instead.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"

namespace line {
namespace sn {

/** The scaling, the utilization normalizer, and why they may not exist. */
struct QrfAlpha {
    std::vector<std::vector<double> > alpha;  ///< (M x N) scaling at population n = 1..N
    std::string msg;                          ///< empty on success
    bool ld = false;                          ///< alpha is not identically 1
    std::vector<double> peak;                 ///< (M) utilization normalizer; inf at a delay
};

/**
 * `ld` stays TRUE through a refusal: the model IS load dependent, and the
 * caller has to tell "no arm serves this" from "the arm you asked for does
 * not". Clearing it would report the latter for both.
 */
template <class T>
QrfAlpha sn_to_qrf_alpha(const qn::NetworkStruct<T>& L) {
    QrfAlpha out;
    const std::size_t M = L.nstations;
    out.peak.assign(M, 1.0);

    double Nd = 0.0;
    const std::vector<double> njobs = L.njobs();
    for (std::size_t r = 0; r < njobs.size(); ++r) Nd += njobs[r];
    if (!(Nd >= 1.0) || std::isinf(Nd)) {
        out.alpha.assign(M, std::vector<double>(1, 1.0));
        out.msg = "the QRF bounds need a closed model with a finite population.";
        return out;
    }
    const std::size_t N = static_cast<std::size_t>(Nd + 0.5);
    out.alpha.assign(M, std::vector<double>(N, 1.0));

    for (std::size_t i = 0; i < M; ++i) {
        const double c = L.stations[i].nservers;
        const bool is_delay = std::isinf(c) || L.stations[i].sched == qn::SchedStrategy::INF;
        const bool serves_many = is_delay || c > 1.0;
        const std::size_t smax = L.stations[i].lldscaling.size();
        // The phase count is read where the adapter reads it, from the service
        // process: that {D0,D1} pair is what sizes the local state, so testing
        // it keeps the refusal and the formulation on one quantity.
        std::size_t ki = 1;
        {
            const mam::Map<T> m = lang::dist_to_map(L.service[i][0]);
            if (m.D0.rows() > 0) ki = static_cast<std::size_t>(m.D0.rows());
        }
        if (serves_many && ki > 1) {
            std::ostringstream os;
            os << "station " << (i + 1) << " serves ";
            if (is_delay)
                os << "unboundedly many";
            else
                os << "up to " << static_cast<int>(c);
            os << " jobs at once with " << ki
               << "-phase service, and the QRF local state carries one phase per station, which "
                  "describes one job in service and no more. Give that station exponential "
                  "service, or use a single-server model.";
            out.msg = os.str();
            out.ld = true;
            return out;
        }
        double lldpeak = 1.0;
        for (std::size_t n = 1; n <= N; ++n) {
            double a = 1.0;
            if (is_delay)
                a = static_cast<double>(n);
            else if (c > 1.0)
                a = std::min(static_cast<double>(n), c);
            if (smax > 0) {
                const double s =
                    num_traits<T>::to_double(L.stations[i].lldscaling[std::min(n, smax) - 1]);
                lldpeak = std::max(lldpeak, s);
                a *= s;
            }
            out.alpha[i][n - 1] = a;
            if (a != 1.0) out.ld = true;
        }
        out.peak[i] = is_delay ? std::numeric_limits<double>::infinity() : c * lldpeak;
    }
    return out;
}

}  // namespace sn
}  // namespace line

#endif  // LINE_API_SN_SN_TO_QRF_ALPHA_H
