/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_SSA_SSA_EVENT_CACHE_H
#define LINE_SOLVERS_SSA_SSA_EVENT_CACHE_H

/**
 * Port of `EventCache.m` and of the lookup that `State.afterEvent` performs
 * against it (afterEvent.m lines 24-44 and its write-back sites).
 *
 * WHAT IS BEING MEMOIZED, AND WHY IT IS SOUND. `after_event` is a pure function
 * of `(sn, ind, inspace, event, cls, no_promote, aux_rate)`: it reads no
 * mutable state, draws no random number and touches nothing outside its
 * arguments. So its value may be stored against those arguments and returned
 * again, and a hit is EXACTLY what recomputation would produce. That is the
 * whole correctness condition, and it is the one the tests assert directly,
 * state by state, rather than inferring it from the metrics.
 *
 * THE KEY IS THE WHOLE ARGUMENT LIST. The reference keys on
 * `mat2str([ind, event, class, noPromote, inspace])`, and each of those five
 * pieces earns its place: `noPromote` distinguishes the departure half of an
 * immediate-feedback self-loop (which leaves the vacated server held) from an
 * ordinary departure at the same station in the same state, and dropping it
 * would return one where the other was asked for. This port keys on the same
 * five plus `aux_rate`, which the C++ handler takes and MATLAB does not: it is
 * the rate carried into the RENEGE, RETRY, FAILURE and REPAIR branches, so two
 * calls that differ only in it have different answers.
 *
 * A HIT AND A MISS ARE INDISTINGUISHABLE HERE, WHICH IS STRONGER THAN THE
 * REFERENCE. MATLAB stores the full enumeration and then, on a hit, SAMPLES one
 * successor row from it, spending a `rand`. This port's `after_event` is the
 * enumeration-mode handler in every case and the sampling happens downstream in
 * the engine, so enabling the cache changes neither the value nor the random
 * stream. The seed-fixed sample path is therefore identical with and without
 * caching, which is what lets a test compare the two runs bit for bit.
 *
 * THE CACHE IS BOUND TO ONE `sn`. The reference's `create` takes `sn` and
 * ignores it (its loops are commented out). Here it is retained and checked:
 * the memoized value is computed from the rates, capacities and routing of one
 * network struct, so serving it to a query about a different struct would
 * return another model's successors. The binding is a pointer identity test
 * because the engines hold `sn` by const reference for their whole lifetime.
 *
 * A DISABLED CACHE IS NOT AN ABSENT ONE. `EventCache.create(false, sn)` returns
 * `[]` and `afterEvent` then takes the uncached path. Here a disabled cache is
 * an object that computes on every call and stores nothing, so the caller has
 * one code path rather than two and cannot accidentally diverge between them.
 */

#include <cstddef>
#include <map>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/state_events.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace ssa {

/**
 * The memoization key: the argument list of `after_event`, in fields.
 *
 * The pieces are kept SEPARATE rather than flattened into one vector because a
 * flatten lets a wide state row of one node impersonate a narrow row of another
 * with a different index prefix; separate fields cannot collide by
 * construction, so no separator sentinel is needed.
 */
struct SsaEventKey {
    std::size_t ind = 0;   ///< 1-based node index
    int event = 0;         ///< `EventType`, as its underlying value
    std::size_t cls = 0;   ///< 1-based class index
    bool no_promote = false;
    double aux_rate = 0.0;
    std::vector<double> inspace;  ///< the node's local state row

    bool operator<(const SsaEventKey& o) const {
        if (ind != o.ind) return ind < o.ind;
        if (event != o.event) return event < o.event;
        if (cls != o.cls) return cls < o.cls;
        if (no_promote != o.no_promote) return !no_promote;
        if (aux_rate != o.aux_rate) return aux_rate < o.aux_rate;
        return inspace < o.inspace;
    }
};

/** `EventCache`: the per-state enabled-event memo of the serial SSA engine. */
template <class T>
class SsaEventCache {
public:
    /**
     * `EventCache.create(enabled, sn)`.
     *
     * The struct is held by pointer, so the cache must not outlive it. Every
     * caller in this port constructs the cache inside the scope that already
     * holds `sn` by const reference for the run.
     */
    static SsaEventCache create(bool enabled, const qn::NetworkStruct<T>& sn) {
        return SsaEventCache(enabled, &sn);
    }

    SsaEventCache() : enabled_(false), sn_(0) {}
    SsaEventCache(bool enabled, const qn::NetworkStruct<T>* sn) : enabled_(enabled), sn_(sn) {}

    bool enabled() const { return enabled_; }

    /**
     * `State.afterEvent` with the lookup in front of it.
     *
     * `sn` is passed rather than taken from the binding so the signature is the
     * free function's and a caller can switch between them by changing one
     * token; the binding is only there to catch the mistake of reusing a cache
     * across structs.
     */
    qn::EventOutcome<T> after_event(const qn::NetworkStruct<T>& sn, std::size_t ind,
                                    const std::vector<T>& inspace, lang::EventType event,
                                    std::size_t cls, bool no_promote = false,
                                    const T& aux_rate = num_traits<T>::from_int(0)) {
        if (!enabled_)
            return qn::after_event(sn, ind, inspace, event, cls, no_promote, aux_rate);
        if (sn_ != 0 && sn_ != &sn)
            throw InputError(
                "SsaEventCache: queried with a different NetworkStruct from the one it was "
                "created against. The memoized successors are a function of that struct's rates, "
                "capacities and routing, so serving them here would answer about another model");

        SsaEventKey key;
        key.ind = ind;
        key.event = static_cast<int>(event);
        key.cls = cls;
        key.no_promote = no_promote;
        key.aux_rate = num_traits<T>::to_double(aux_rate);
        key.inspace.resize(inspace.size());
        for (std::size_t j = 0; j < inspace.size(); ++j)
            key.inspace[j] = num_traits<T>::to_double(inspace[j]);

        const typename std::map<SsaEventKey, qn::EventOutcome<T> >::const_iterator it =
            memo_.find(key);
        if (it != memo_.end()) {
            ++hits_;
            return it->second;
        }
        ++misses_;
        const qn::EventOutcome<T> out =
            qn::after_event(sn, ind, inspace, event, cls, no_promote, aux_rate);
        memo_[key] = out;
        return out;
    }

    /**
     * Hits and misses, so a caller can report the memo's effect.
     *
     * They are counters and not a hit RATIO because the ratio alone hides the
     * denominator, and a 100% hit rate over three lookups says nothing.
     */
    std::size_t hits() const { return hits_; }
    std::size_t misses() const { return misses_; }
    std::size_t size() const { return memo_.size(); }

    void clear() {
        memo_.clear();
        hits_ = 0;
        misses_ = 0;
    }

private:
    bool enabled_;
    const qn::NetworkStruct<T>* sn_;
    std::map<SsaEventKey, qn::EventOutcome<T> > memo_;
    std::size_t hits_ = 0, misses_ = 0;
};

}  // namespace ssa
}  // namespace line

#endif  // LINE_SOLVERS_SSA_SSA_EVENT_CACHE_H
