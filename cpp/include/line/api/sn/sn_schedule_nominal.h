/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_SCHEDULE_NOMINAL_H
#define LINE_API_SN_SN_SCHEDULE_NOMINAL_H

/**
 * Port of `matlab/src/api/sn/sn_schedule_nominal.m`: unpack the MAPt or PHt slot
 * of `sn.proc` at one (station, class) into the per-segment (D0, D1) pairs and
 * their width-weighted time average.
 *
 * IT IS AN ACCESSOR HERE AND A DECODER THERE, and the difference is where the
 * conversion happens. MATLAB stores the slot as `{breakpoints, A, B, cyclic}`
 * with A and B holding (D0, D1) for a MAPt and (alpha, S) for a PHt, so every
 * reader must know which family it is looking at and rebuild the MAP pair
 * (S, (-S e) alpha) itself. This port converts a PHt ONCE, in
 * `Distrib::pht`, so what is stored is always the MAP schedule and every reader
 * sees one representation. The nominal pair is likewise computed once, in the
 * constructor, and lives in `Distrib::D0` / `D1` where a time-blind consumer --
 * the phase count, `sn.rates`, the fluid layout -- already looks for it.
 *
 * What is left for this function is the CONTRACT: return the same six outputs
 * the reference returns, and refuse a (station, class) that carries no schedule
 * rather than inventing a one-segment one. A caller that wants the nominal of an
 * ordinary process reads `sn.service[i][r].D0` directly, which is what the
 * reference's own `nomD0`/`nomD1` fallback does.
 */

#include <cstddef>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace sn {

/** What `sn_schedule_nominal` returns, in the reference's own output order. */
template <class T>
struct ScheduleNominal {
    Matrix<T> D0bar, D1bar;             ///< the width-weighted time average
    std::vector<T> breakpoints;         ///< the boundary vector, nseg + 1 long
    std::vector<Matrix<T> > segD0;      ///< per-segment D0, already in MAP form
    std::vector<Matrix<T> > segD1;      ///< per-segment D1
    bool cyclic = false;
};

/** True when (ist, r) carries a MAPt / PHt / NHPP schedule. Both 0-based. */
template <class T>
bool sn_has_schedule(const qn::NetworkStruct<T>& sn, std::size_t ist, std::size_t r) {
    if (ist >= sn.service.size() || r >= sn.service[ist].size()) return false;
    return sn.service[ist][r].has_schedule();
}

/** Port of `sn_schedule_nominal(sn, ist, r)`; `ist` and `r` are 0-based here. */
template <class T>
ScheduleNominal<T> sn_schedule_nominal(const qn::NetworkStruct<T>& sn, std::size_t ist,
                                       std::size_t r) {
    if (!sn_has_schedule(sn, ist, r))
        throw InputError("sn_schedule_nominal: station " + std::to_string(ist + 1) + ", class " +
                         std::to_string(r + 1) +
                         " carries no MAPt/PHt schedule; read sn.service[i][r].D0 for an "
                         "ordinary process");
    const lang::Distrib<T>& d = sn.service[ist][r];
    ScheduleNominal<T> out;
    out.D0bar = d.D0;
    out.D1bar = d.D1;
    out.breakpoints = d.sched_bp;
    out.segD0 = d.sched_D0;
    out.segD1 = d.sched_D1;
    out.cyclic = d.sched_cyclic;
    return out;
}

}  // namespace sn
}  // namespace line

#endif  // LINE_API_SN_SN_SCHEDULE_NOMINAL_H
