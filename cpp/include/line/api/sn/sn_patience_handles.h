/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_PATIENCE_HANDLES_H
#define LINE_API_SN_SN_PATIENCE_HANDLES_H

/**
 * Port of matlab/src/api/sn/sn_patience_handles.m.
 *
 * The abandonment solvers need the patience law as FUNCTIONS -- a complementary
 * cdf and a hazard rate -- not as moments, because that is what the underlying
 * theory consumes: Whitt's engineering solution reads the hazard near the
 * origin, and the fluid models integrate the ccdf. This port reads the law from
 * `Station::patience`, which holds the distribution itself rather than the
 * (D0,D1) pair the MATLAB and Java structs carry.
 *
 * ARITHMETIC: transcendental. The ccdf of a phase-type law is a matrix
 * exponential.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>

#include "line/api/mam/map_cdf.h"
#include "line/api/mam/map_pdf.h"
#include "line/api/qsys/qsys_mgisrgi_whitt.h"
#include "line/api/qsys/qsys_types.h"
#include "line/lang/distribution.h"
#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"

namespace line {
namespace api {

/** The patience law of one station-class pair, in the forms the solvers consume. */
template <class T>
struct PatienceHandles {
    /** Whether the pair declares reneging at all; everything below is unset when false. */
    bool present = false;
    /** F^c(t) = P(patience > t). */
    std::function<T(const T&)> ccdf;
    /** The patience density. */
    std::function<T(const T&)> pdf;
    /** The hazard rate h = f/(1-F). */
    std::function<T(const T&)> hazard;
    /** Mean patience. */
    T mean = num_traits<T>::from_int(0);
    /** Whether the law is exponential, in which case the analysis is exact. */
    bool isExponential = false;
    /** The abandonment rate 1/mean. */
    T rate = num_traits<T>::from_int(0);

    /** This law as the `Patience` argument of the abandonment solvers. */
    qsys::Patience<T> as_patience() const {
        if (isExponential) return qsys::Patience<T>::exponential(rate);
        return qsys::Patience<T>::hazard(hazard);
    }
};

/**
 * Build the patience handles of station `ist` (0-based), class `r`.
 *
 * Returns a handle set with `present == false` when the pair declares no
 * reneging patience.
 */
template <class T>
PatienceHandles<T> sn_patience_handles(const qn::NetworkStruct<T>& sn, std::size_t ist,
                                       std::size_t r) {
    PatienceHandles<T> h;
    if (ist >= sn.nstations || r >= sn.nclasses) return h;
    const qn::Station<T>& st = sn.stations[ist];
    if (r >= st.impatience.size() || st.impatience[r] != lang::ImpatienceType::RENEGING) return h;
    if (r >= st.patience.size() || st.patience[r].disabled) return h;

    const lang::Distrib<T>& d = st.patience[r];
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    (void)one;
    h.present = true;
    h.mean = d.mean;
    h.rate = d.mean > zero ? T(one / d.mean) : zero;
    h.isExponential = (d.type == lang::ProcessType::EXP);

    // WHAT IS AVAILABLE UNDER EXACT ARITHMETIC, and what is not. Whether the
    // pair reneges, at what rate and under which family are all read off the
    // struct, so `present`, `mean`, `rate` and `isExponential` are answered for
    // every T -- which is what lets the method list and the method resolution
    // run under Rational. The ccdf and the hazard are a matrix exponential, so
    // they exist only where the arithmetic has one; the abandonment methods
    // themselves are refused by name before they are ever called.
    if constexpr (num_traits<T>::has_transcendental) {
        if (h.isExponential) {
            const T theta = h.rate;
            h.ccdf = [theta](const T& t) { return qsys::detail::num_exp(T(-theta * t)); };
            h.pdf = [theta](const T& t) { return T(theta * qsys::detail::num_exp(T(-theta * t))); };
            h.hazard = [theta](const T&) { return theta; };
            return h;
        }

        const mam::Map<T> m = lang::dist_to_map(d);
        h.ccdf = [m, one](const T& t) {
            std::vector<T> pts(1, t);
            return T(one - mam::map_cdf(m, pts)[0]);
        };
        h.pdf = [m](const T& t) {
            std::vector<T> pts(1, t);
            return mam::map_pdf(m, pts)[0];
        };
        const T asymptotic = h.rate;
        std::function<T(const T&)> ccdf = h.ccdf;
        std::function<T(const T&)> pdf = h.pdf;
        h.hazard = [ccdf, pdf, asymptotic](const T& t) {
            // h = f/(1-F). Past the point where the ccdf underflows the hazard
            // is the asymptotic decay rate, and returning that is better
            // conditioned than dividing two zeros.
            const T c = ccdf(t);
            if (num_traits<T>::to_double(c) <= 1e-300) return asymptotic;
            return T(pdf(t) / c);
        };
    }
    return h;
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_PATIENCE_HANDLES_H
