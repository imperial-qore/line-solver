/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_GET_BUFFER_SIZE_H
#define LINE_API_SN_SN_GET_BUFFER_SIZE_H

/**
 * Physical buffer size of a station, in jobs, the one in service included.
 *
 * Port of matlab/src/api/sn/sn_get_buffer_size.m. Kendall's K: the tighter of
 * the station capacity sn.cap and the SUM of the per-class capacities
 * sn.classcap, +inf when the station is unbounded. Both fields already fold
 * setCapacity, setClassCapacity, a finite orbit and the closed-chain
 * population, so this is the single place that decides whether a buffer BINDS.
 *
 * Two details that are easy to get wrong and both matter. The per-class
 * capacities are SUMMED, not minimised: a zero marks a class the station does
 * not serve and is dropped first. And the "unreachable capacity" suppression
 * compares against the population that can actually REACH this station, i.e.
 * the classes it serves, not the model's total, without which a MIXED model
 * reads as finite-buffered at every station an open class never visits.
 *
 * This lived inside solvers/nc/solver_nc_mem.h as `nc::detail::buffer_size`.
 * It is the gate of every finite-buffer branch, not an NC internal, so it now
 * sits on the api surface. The BODY is `NetworkStruct::buffer_size`, which this
 * forwards to: the struct's own `has_blocking` (the product-form conjunct the
 * cpp solvers read through `has_product_form`) needs the same rule, and two
 * copies of a predicate this load-bearing is exactly how the member and the
 * free function drift apart.
 */

#include <algorithm>
#include <cstddef>
#include <limits>

#include "line/lang/qn/network_struct.h"

namespace line {
namespace sn {

/** @param ist 1-based station index. */
template <class T>
double sn_get_buffer_size(const qn::NetworkStruct<T>& sn, std::size_t ist) {
    return sn.buffer_size(ist);
}

}  // namespace sn
}  // namespace line

#endif  // LINE_API_SN_SN_GET_BUFFER_SIZE_H
