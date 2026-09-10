/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_IS_DISCRETE_TIME_H
#define LINE_API_SN_SN_IS_DISCRETE_TIME_H

/**
 * Port of matlab/src/api/sn/sn_is_discrete_time.m.
 *
 * Decides whether every interarrival and service law of a model lives on the
 * slot lattice {slotLength, 2*slotLength, ...}, which is what lets SolverMAM
 * take the discrete-time (Q-MAM) path instead of the continuous one.
 *
 * The test reads `procid`, `rates` and `scv` and NOT the fitted `proc`
 * matrices: by the time a solver sees the struct a Geometric has already been
 * converted to a CONTINUOUS MAP, so the lattice is no longer visible there. A
 * DMAP is the exception, its (D0,D1) having MAP shape it survives verbatim.
 * For the same reason this test must run BEFORE any phase-type conversion.
 *
 * ARITHMETIC: field, with one square root in the DiscreteUniform width. Every
 * other test is a comparison, so the family instantiates under Rational apart
 * from that branch, which needs a real square root.
 */

#include <cmath>
#include <cstddef>
#include <string>

#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace api {

/** Which laws the model carries, and why it was refused when it was. */
struct DiscreteTimeInfo {
    bool has_lattice = false;     ///< some law lives on the slot lattice
    bool has_continuous = false;  ///< some law does not
    bool mixed = false;           ///< both, so no single time scale fits
    bool has_dmap = false;        ///< a DMAP is present, which has no continuous reading
    std::string reason;           ///< empty when the model is discrete-time
};

/** How the caller may override the automatic decision. */
struct DiscreteTimeOptions {
    /** "auto", "discrete" or "continuous". */
    std::string timescale = "auto";
    /** Slot length in model time units. */
    double slotlength = 1.0;
};

namespace detail {

/**
 * Integral bounds of a DiscreteUniform recovered from its mean and SCV, using
 * var = ((hi-lo+1)^2-1)/12. Returns false when the SCV is not usable.
 */
inline bool duniform_bounds(double mean_slots, double scv, double* lo, double* hi) {
    if (!std::isfinite(scv) || scv < 0) return false;
    const double var_slots = scv * mean_slots * mean_slots;
    const double width = std::sqrt(std::max(0.0, 12 * var_slots + 1)) - 1;
    *lo = std::round(mean_slots - width / 2);
    *hi = std::round(mean_slots + width / 2);
    return true;
}

}  // namespace detail

/**
 * True when every law of `sn` is lattice-valued on `slot_length`.
 *
 * `slot_length` receives the slot in model time units and `info` the diagnosis.
 * A DMAP mixed with continuous laws is an ERROR rather than a `false`: its
 * (D0,D1) are probability matrices, so the continuous machinery would form
 * inv(-D0) where the law needs inv(I-D0) and return a wrong number in silence.
 */
template <class T>
bool sn_is_discrete_time(const qn::NetworkStruct<T>& sn, const DiscreteTimeOptions& options,
                         double* slot_length, DiscreteTimeInfo* info) {
    using lang::ProcessType;

    const std::string& timescale = options.timescale;
    if (timescale != "auto" && timescale != "discrete" && timescale != "continuous")
        throw InputError(
            "sn_is_discrete_time: timescale must be 'auto', 'discrete' or 'continuous'");
    if (!std::isfinite(options.slotlength) || options.slotlength <= 0)
        throw InputError("sn_is_discrete_time: slotlength must be a positive finite scalar");

    *slot_length = options.slotlength;
    *info = DiscreteTimeInfo();
    if (timescale == "continuous") return false;

    const double tol = 1e-12;
    for (std::size_t ist = 0; ist < sn.nstations; ++ist)
        for (std::size_t r = 0; r < sn.nclasses; ++r) {
            const ProcessType pt = sn.procid(ist + 1, r + 1);
            if (pt == ProcessType::DISABLED || pt == ProcessType::NONE) continue;
            if (sn.disabled[ist][r]) continue;
            const double rate = num_traits<T>::to_double(sn.rates(ist, r));
            if (!std::isfinite(rate) || rate <= 0) continue;  // no interval law here
            const double mean_slots = 1.0 / (rate * options.slotlength);

            if (pt == ProcessType::GEOMETRIC) {
                // mean 1/p slots; p in (0,1] is recovered exactly from the mean
                info->has_lattice = true;
                if (mean_slots < 1 - tol) {
                    info->has_continuous = true;
                    info->reason = "a Geometric has a mean below the one-slot minimum of its "
                                   "support {1,2,...}";
                }
            } else if (pt == ProcessType::DMAP) {
                info->has_lattice = true;
                info->has_dmap = true;
            } else if (pt == ProcessType::DUNIFORM) {
                info->has_lattice = true;
                double lo = 0, hi = 0;
                if (!detail::duniform_bounds(mean_slots, num_traits<T>::to_double(sn.scv(ist, r)),
                                             &lo, &hi) ||
                    lo < 1 - tol) {
                    info->has_continuous = true;
                    info->reason = "a DiscreteUniform spans a range not contained in {1,2,...}";
                }
            } else if (pt == ProcessType::DET) {
                if (std::abs(mean_slots - std::round(mean_slots)) <= tol * std::max(1.0, mean_slots) &&
                    std::round(mean_slots) >= 1) {
                    info->has_lattice = true;
                } else {
                    // a Det off the lattice is what makes the model continuous
                    info->has_continuous = true;
                }
            } else {
                info->has_continuous = true;
            }
        }

    info->mixed = info->has_lattice && info->has_continuous;

    if (timescale == "discrete") {
        if (info->mixed)
            throw InputError("sn_is_discrete_time: timescale='discrete' was requested but the "
                             "model mixes lattice and non-lattice laws: " +
                             info->reason);
        if (!info->has_lattice)
            throw InputError("sn_is_discrete_time: timescale='discrete' was requested but no "
                             "interarrival or service law is lattice-valued on the given slot");
        return true;
    }

    const bool bool_out = info->has_lattice && !info->has_continuous;
    if (!bool_out && info->has_dmap)
        throw InputError("sn_is_discrete_time: the model mixes a DMAP with continuous-time laws. "
                         "A DMAP is only defined on a slotted time scale, so no solver can "
                         "interpret this model: " +
                         info->reason);
    if (!bool_out && info->reason.empty() && info->mixed)
        info->reason = "the model mixes lattice-valued and continuous laws";
    return bool_out;
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_IS_DISCRETE_TIME_H
