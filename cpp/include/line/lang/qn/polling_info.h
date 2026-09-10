/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * `State.pollingInfo` and the controller description it returns.
 *
 * IT LIVES IN ITS OWN HEADER because both sides of the State package need it
 * and neither can include the other: `from_marginal` (state.h) has to ENUMERATE
 * the controller configurations a station can occupy, and the event handlers
 * (state_events.h) have to MOVE between them. Leaving it in state_events.h left
 * `from_marginal` unable to emit the controller columns at all, so a polling
 * station's states were built one block too narrow and the enumerated space
 * contained no controller -- the station then behaved as a capacity-one queue.
 */
#ifndef LINE_LANG_QN_POLLING_INFO_H
#define LINE_LANG_QN_POLLING_INFO_H

#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace qn {

using lang::ProcessType;
using lang::SchedStrategy;

/**
 * Port of `State.pollingInfo`: the derived description of a polling controller.
 *
 * The controller lives in the shared node block of the local variables and
 * holds, in order, [pos, swk, ctr]. Each column is materialized ONLY when the
 * discipline needs it, so a polling station never carries state its dynamics
 * cannot distinguish:
 *
 *   pos  the buffer the server is at or heading to. Needed only when some
 *        switchover is non-immediate: while a job is in service pos equals
 *        that job's class, and while the station is empty with immediate
 *        switchovers the position is unobservable.
 *   swk  0 when the server sits at pos, else the phase of the switchover into
 *        it. Same condition as pos.
 *   ctr  the visit budget. Every discipline except EXHAUSTIVE bounds a visit.
 *
 * An Immediate switchover is deliberately NOT a state: its rate is ~1e8, and
 * taking that literally would make the generator stiff and add a spurious
 * state per buffer. `polling_next` folds it into the enclosing transition.
 *
 * EVERY POLLING STATION HAS A CONTROLLER, whichever API declared it and even
 * when none did. The reference returns [] only for a node that is not a polling
 * station and otherwise defaults to EXHAUSTIVE with every switchover immediate,
 * so `NetworkStruct::effective_polling` resolves the two C++ spellings --
 * `set_polling` writing `sn.pollingparam`, and setPollingType / setSwitchover
 * writing the Station fields -- into one description.
 *
 * BUILDING NO CONTROLLER WAS NOT A DEGRADED MODEL BUT A CRASH. Reading only
 * `pollingparam` left every station declared the other way (the JSON reader
 * included) with `valid = false` and, with it, an EMPTY `polled`. `polling_next`
 * indexes `polled[pos-1]` with no size test, and an empty std::vector<bool>
 * holds a NULL word pointer, so that index segfaults rather than returning a
 * garbage bit -- which is why the failure surfaced as a SIGSEGV inside
 * `after_event_station_arv` and not as a wrong number.
 */
template <class T>
struct PollingInfo {
    lang::PollingType ptype = lang::PollingType::EXHAUSTIVE;
    std::size_t pk = 1;
    std::vector<bool> polled, has_sw;
    std::vector<std::size_t> ksw;
    std::vector<Matrix<T>> sw_d0, sw_d1;
    std::vector<std::vector<T>> sw_pie;
    std::size_t off = 0, width = 0;
    // 0-based columns inside the local-variable block, or npos when absent.
    std::size_t ipos = static_cast<std::size_t>(-1);
    std::size_t iswk = static_cast<std::size_t>(-1);
    std::size_t ictr = static_cast<std::size_t>(-1);
    bool valid = false;
};

template <class T>
PollingInfo<T> polling_info(const NetworkStruct<T>& sn, std::size_t ind) {
    const std::size_t R = sn.nclasses;
    PollingInfo<T> pi;
    const std::size_t ist = sn.nodes[ind - 1].station;
    if (ist == 0 || sn.stations[ist - 1].sched != SchedStrategy::POLLING) return pi;
    // Whichever API declared the controller, and a default when neither did.
    const typename NetworkStruct<T>::PollingParam pp = sn.effective_polling(ist);
    pi.valid = true;
    pi.ptype = pp.ptype;
    pi.pk = pp.pk;

    // A class disabled at this station can never hold a job, so it is dropped
    // from the cyclic order rather than polled forever.
    pi.polled.assign(R, true);
    for (std::size_t r = 1; r <= R; ++r)
        if (sn.disabled[ist - 1][r - 1] || sn.phases_of(ist, r) == 0) pi.polled[r - 1] = false;
    bool any = false;
    for (std::size_t r = 0; r < R; ++r) any = any || pi.polled[r];
    if (!any)
        throw InputError("polling_info: the polling station has no class with an enabled service");

    pi.has_sw.assign(R, false);
    pi.ksw.assign(R, 0);
    pi.sw_d0.resize(R);
    pi.sw_d1.resize(R);
    pi.sw_pie.resize(R);
    // The switchover of the leg ENTERING each polled buffer is read off the
    // buffer the server leaves to get there.
    std::vector<std::size_t> polled_list;
    for (std::size_t r = 1; r <= R; ++r)
        if (pi.polled[r - 1]) polled_list.push_back(r);
    for (std::size_t j = 0; j < polled_list.size(); ++j) {
        const std::size_t q = polled_list[j];
        const std::size_t prev = polled_list[(j + polled_list.size() - 1) % polled_list.size()];
        if (pp.switchover.size() < prev) continue;
        const lang::Distrib<T>& d = pp.switchover[prev - 1];
        if (d.disabled || d.D0.rows() == 0) continue;
        if (d.type == ProcessType::IMMEDIATE) continue;  // folded, never a state
        pi.has_sw[q - 1] = true;
        pi.ksw[q - 1] = d.D0.rows();
        pi.sw_d0[q - 1] = d.D0;
        pi.sw_d1[q - 1] = d.D1;
        mam::Map<T> mp;
        mp.D0 = d.D0;
        mp.D1 = d.D1;
        pi.sw_pie[q - 1] = mam::map_pie(mp);
    }

    bool anysw = false;
    for (std::size_t r = 0; r < R; ++r) anysw = anysw || pi.has_sw[r];
    // pos and swk exist only to encode SWITCHING(p); with every switchover
    // immediate the server is either serving (pos = the class in service) or
    // parked (pos unobservable), so neither column carries information.
    std::size_t c = 0;
    if (anysw) {
        pi.ipos = c++;
        pi.iswk = c++;
    }
    if (pi.ptype != lang::PollingType::EXHAUSTIVE) pi.ictr = c++;
    pi.width = c;
    // Guarded for the same reason nvars_of is: a hand-built struct can carry a
    // shorter nvars than the node index, and sn.nvars[ind-1] is evaluated before
    // the inner size() test would ever run.
    std::size_t off = 0;
    if (ind <= sn.nvars.size())
        for (std::size_t j = 0; j < 2 * R && j < sn.nvars[ind - 1].size(); ++j)
            off += sn.nvars[ind - 1][j];
    pi.off = off;
    return pi;
}
}  // namespace qn
}  // namespace line

#endif  // LINE_LANG_QN_POLLING_INFO_H
