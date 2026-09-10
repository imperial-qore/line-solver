/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_REMOVE_CLASS_H
#define LINE_API_SN_SN_REMOVE_CLASS_H

/**
 * Drop one job class from a model, port of `ModelAdapter.removeClass`.
 *
 * The reference is `matlab/src/lang/@@MNetwork/removeClass.m` (the mutating
 * form) with `matlab/src/io/@@ModelAdapter/removeClass.m` as its non-mutating
 * wrapper; the JAR twins are `Network.removeClass` / `Network.withoutClass`
 * and the python ones `Network.remove_class` / `Network.without_class`. This
 * port has only the NON-MUTATING form, because a `NetworkStruct` is a value
 * here rather than a handle: `V = sn_remove_class(sn, r)` leaves `sn` intact,
 * which is what `withoutClass` buys the reference, and a caller who wants the
 * mutating form assigns over its own struct.
 *
 * WHAT IT IS FOR. Ablation studies (solve the model without class r and read
 * off what r was contributing) and per-chain decomposition, where one class is
 * peeled at a time. Both need the REST of the model to stay solvable, which is
 * why every per-class table has to be re-indexed rather than merely blanked.
 *
 * WHY THIS IS A SLICE AND NOT A DELETION. Every per-class table in the struct
 * is POSITIONAL: `service[i][r]`, `st.schedparam[r]`, `nd.routing[r]`, the
 * (K x K) matrix of a ClassSwitch, the (r,s) key of a routing block. Blanking
 * entry r and leaving the tables at width K would leave a model whose class
 * list is one shorter than the tables that index it, and the refresh would then
 * read class r+1's service under class r's name. That is the same defect the
 * reference hit from the other side -- MATLAB's `removeClass` used to leave the
 * (K x K) class-switching mask at its old size, and `refreshRoutingMatrix`
 * indexed `sn.refstat` out of range on the next `getStruct` -- so every table
 * is sliced here and the struct is then re-refreshed from the sliced stage-one
 * data.
 *
 * ONLY STAGE-ONE DATA IS TOUCHED. `network_builder.h` states the split: stage
 * one constructs, stage two (`refresh_*`) derives. So this file slices the
 * constructed tables (classes, service, the station and node per-class vectors,
 * the ClassSwitch matrices, the routing blocks) and then calls
 * `refresh_struct()`, which re-derives `rates`, `scv`, `disabled`, `chains`,
 * `inchain`, `visits`, `cap`, `rt`, `rtnodes` and the rest at the new width.
 * Nothing derived is patched by hand, so no derived field can be left stale.
 *
 * WHAT IS REFUSED, and why refusing beats slicing. The reference removes a
 * class from five node kinds (`Node`, `Station`, `ServiceStation`, `Source`,
 * `ClassSwitch`) and REFUSES on a Cache, whose item state is indexed by class
 * and lives in the model state rather than in the node. This port keeps that
 * refusal verbatim and adds one for every other construct that carries a class
 * INDEX in something other than a plain per-class vector -- a fork's per-class
 * matrices, a transition's modes, a region's per-class caps, a retrial or
 * polling or PAS block, a signal's target class, a class-dependent scaling
 * function whose argument is a population vector of width K. Those cannot be
 * re-indexed by slicing, and slicing around them would hand back a model that
 * looks well formed and is not. `tag_chain` refuses the same way and for the
 * same reason.
 *
 * ARITHMETIC: none. Structural.
 */

#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace api {

namespace remove_class_detail {

/**
 * Erase class `r` (1-based) from a table that IS per-class of this model.
 *
 * The `size() == K` guard is the reference's own: `@@MNetwork` tests
 * `numel(self.classCap) == K` before slicing, because these vectors are
 * OPTIONAL -- an empty one means the station declares none of that property at
 * all, and a vector of some other width belongs to something that is not the
 * class axis and must be left alone.
 */
template <class V>
void erase_slot(std::vector<V>& v, std::size_t K, std::size_t r) {
    if (v.size() != K) return;
    v.erase(v.begin() + static_cast<std::ptrdiff_t>(r - 1));
}

/** Drop row and column `r` (1-based) of a (K x K) matrix. */
template <class T>
void erase_row_col(Matrix<T>& C, std::size_t K, std::size_t r) {
    if (C.rows() != K || C.cols() != K) return;
    if (K <= 1) {
        C = Matrix<T>();
        return;
    }
    Matrix<T> R(K - 1, K - 1, num_traits<T>::from_int(0));
    std::size_t ai = 0;
    for (std::size_t a = 0; a < K; ++a) {
        if (a == r - 1) continue;
        std::size_t bi = 0;
        for (std::size_t b = 0; b < K; ++b) {
            if (b == r - 1) continue;
            R(ai, bi) = C(a, b);
            ++bi;
        }
        ++ai;
    }
    C = R;
}

/**
 * Re-index a list of CLASS INDICES: drop `r` and shift what came after it.
 *
 * `Source.markedClasses` is such a list, not a per-class flag vector, so
 * erasing position r would drop the wrong entry.
 */
inline void reindex_class_list(std::vector<std::size_t>& list, std::size_t r) {
    std::vector<std::size_t> keep;
    keep.reserve(list.size());
    for (std::size_t a = 0; a < list.size(); ++a) {
        const std::size_t m = list[a];
        if (m == r) continue;
        keep.push_back(m > r ? m - 1 : m);
    }
    list.swap(keep);
}

/** Everything whose class indexing this transform cannot rewrite, named. */
template <class T>
void remove_class_check(const qn::NetworkStruct<T>& sn) {
    // The reference's own refusal, with its own wording: the cache item state
    // is indexed by class and lives in the model state, not in the node.
    if (!sn.nodeparam.empty())
        throw UnsupportedError(
            "sn_remove_class: cannot dynamically remove classes in models with caches. You need "
            "to re-instantiate the model.");
    if (sn.has_fork() || !sn.fj.empty() || !sn.forkparam.empty() || !sn.joindecl.empty())
        throw UnsupportedError(
            "sn_remove_class: a Fork carries per-class branch probabilities and forking levels, "
            "and a Join a per-class quorum, none of which is a plain per-class vector; rebuild "
            "the model without the class instead");
    if (!sn.transparam.empty() || !sn.initmarking.empty() || !sn.statespace.empty() ||
        !sn.stateprior.empty())
        throw UnsupportedError(
            "sn_remove_class: a Petri-net Place or Transition indexes its modes, its arc "
            "multiplicities and its initial marking by token class; rebuild the model without "
            "the class instead");
    if (!sn.regions.empty())
        throw UnsupportedError(
            "sn_remove_class: a finite-capacity region carries a per-class capacity row and "
            "per-class weights; rebuild the model without the class instead");
    if (!sn.retrialparam.empty())
        throw UnsupportedError("sn_remove_class: a retrial orbit is parameterized per class");
    if (!sn.pollingparam.empty())
        throw UnsupportedError(
            "sn_remove_class: a polling station serves one buffer per class and carries a "
            "per-class switchover schedule");
    if (!sn.pasparam.empty())
        throw UnsupportedError(
            "sn_remove_class: a pass-and-swap station is parameterized by the ordered microstate "
            "of class indices, which no slice can rewrite");
    if (!sn.setupparam.empty())
        throw UnsupportedError("sn_remove_class: a setup / delayoff schedule is per class");
    if (!sn.reward.empty())
        throw UnsupportedError(
            "sn_remove_class: a reward function names its classes and would be silently "
            "re-pointed by the shift");
    for (std::size_t r = 0; r < sn.issignal.size(); ++r)
        if (sn.issignal[r])
            throw UnsupportedError(
                "sn_remove_class: a G-network signal names a TARGET class, which the shift would "
                "re-point at a different class");
    for (std::size_t r = 0; r < sn.syncreply.size(); ++r)
        if (sn.syncreply[r] != 0)
            throw UnsupportedError(
                "sn_remove_class: a synchronous reply names its reply class, which the shift "
                "would re-point at a different class");
    if (sn.gdscaling)
        throw UnsupportedError(
            "sn_remove_class: a global-dependence function takes a population vector of width "
            "nclasses and cannot be narrowed");
    for (std::size_t i = 0; i < sn.stations.size(); ++i) {
        if (sn.stations[i].cdscaling || sn.stations[i].jdscaling)
            throw UnsupportedError("sn_remove_class: station '" + sn.stations[i].name +
                                   "' carries a class- or joint-dependent scaling function, "
                                   "whose argument is a population vector of width nclasses");
        if (sn.stations[i].svc_rate_fun)
            throw UnsupportedError("sn_remove_class: station '" + sn.stations[i].name +
                                   "' carries a service-rate function of the ordered microstate "
                                   "of class indices, which no slice can rewrite");
    }
}

}  // namespace remove_class_detail

/**
 * The model without class `cls` (1-based), leaving `sn` untouched.
 *
 * @param sn  a model; it need not be refreshed, since the result is refreshed here
 * @param cls the 1-based class to remove
 */
template <class T>
qn::NetworkStruct<T> sn_remove_class(const qn::NetworkStruct<T>& sn, std::size_t cls) {
    using namespace remove_class_detail;
    const std::size_t K = sn.classes.size();
    if (cls == 0 || cls > K)
        throw InputError("sn_remove_class: class index " + std::to_string(cls) +
                         " is out of range");
    // The reference errors rather than returning an empty model, and says why.
    if (K <= 1)
        throw InputError(
            "The network has a single class, it cannot be removed from the model.");
    remove_class_check(sn);

    qn::NetworkStruct<T> V = sn;
    const std::size_t r = cls;

    // ---- nodes: the output (routing) strategy, `Node.removeJobClass` --------
    for (std::size_t i = 0; i < V.nodes.size(); ++i) {
        qn::NodeDef& nd = V.nodes[i];
        erase_slot(nd.routing, K, r);
        erase_slot(nd.routing_weights, K, r);
        erase_slot(nd.routing_param, K, r);
    }

    // ---- ClassSwitch: the (K x K) matrix, `ClassSwitch.removeJobClass` ------
    //
    // Sliced, NOT renormalised, which is the reference's `csFun(remaining,
    // remaining)`. Renormalising would invent a switching probability the model
    // never declared, so a row that summed to one only because of the removed
    // class is left SUB-STOCHASTIC. That is a real consequence and not a
    // theoretical one: removing a class that another class switches INTO leaves
    // the switcher with nowhere to go at that node, and what the refresh then
    // derives is a model with a leak rather than an error. Remove the whole
    // chain, or rebuild, when the class is a switch target.
    for (typename std::map<std::size_t, Matrix<T> >::iterator it = V.csmatrix.begin();
         it != V.csmatrix.end(); ++it)
        erase_row_col(it->second, K, r);

    // ---- stations: Station / ServiceStation / Source.removeJobClass ---------
    for (std::size_t i = 0; i < V.stations.size(); ++i) {
        qn::Station<T>& st = V.stations[i];
        // Station: the per-class buffer, blocking rule and patience.
        erase_slot(st.classcap, K, r);
        erase_slot(st.droprule, K, r);
        erase_slot(st.patience, K, r);
        erase_slot(st.impatience, K, r);
        erase_slot(st.orbit_impatience, K, r);
        erase_slot(st.balking, K, r);
        erase_slot(st.batch_reject, K, r);
        erase_slot(st.immfeed, K, r);
        // ServiceStation: the scheduling parameter, and the service processes
        // held inside the server section.
        erase_slot(st.schedparam, K, r);
        erase_slot(st.server_parallelism, K, r);
        erase_slot(st.departure_discipline, K, r);
        erase_slot(st.cdscalingpeak, K, r);
        erase_slot(st.jdscalingpeak, K, r);
        for (std::size_t k = 0; k < st.server_types.size(); ++k) {
            erase_slot(st.server_types[k].compatible, K, r);
            erase_slot(st.server_types[k].service, K, r);
        }
        // Source: the arrival batch and the marked classes. The arrival
        // PROCESS is the Source's row of `service` and is sliced below with
        // every other station's, since this port keeps the two in one table
        // exactly as the reference's `sn.rates` does.
        erase_slot(st.arrival_batch, K, r);
        reindex_class_list(st.marked_classes, r);
    }

    // ---- the service / arrival table ---------------------------------------
    for (std::size_t i = 0; i < V.service.size(); ++i) erase_slot(V.service[i], K, r);

    // ---- the class list, and the class indices classes hold ----------------
    V.classes.erase(V.classes.begin() + static_cast<std::ptrdiff_t>(r - 1));
    for (std::size_t q = 0; q < V.classes.size(); ++q) {
        std::size_t& sp = V.classes[q].spawn;
        if (sp == r) {
            // The spawned class is gone, so the completion spawns nothing. 0 is
            // this port's absent-index sentinel (MATLAB stores -1).
            sp = 0;
        } else if (sp > r) {
            --sp;
        }
    }

    // ---- the routing blocks, keyed by the (r,s) class pair ------------------
    //
    // A block that departs in the removed class or arrives in it is dropped
    // whole; the survivors keep their probabilities and only their KEY shifts.
    // `Peff` is not sliced but cleared: `refresh_routing` rebuilds it from `P`
    // on every call, so slicing it would be work that is immediately discarded.
    std::map<std::pair<std::size_t, std::size_t>, Matrix<T> > P2;
    for (typename std::map<std::pair<std::size_t, std::size_t>, Matrix<T> >::const_iterator it =
             V.P.begin();
         it != V.P.end(); ++it) {
        const std::size_t a = it->first.first, b = it->first.second;
        if (a == r || b == r) continue;
        P2[std::make_pair(a > r ? a - 1 : a, b > r ? b - 1 : b)] = it->second;
    }
    V.P.swap(P2);
    V.Peff.clear();

    V.refresh_struct();
    return V;
}

/**
 * The same, naming the class.
 *
 * THE LOOKUP BY NAME IS NOT A CONVENIENCE. `ModelAdapter.removeClass` works on
 * a COPY of the model, whose class objects are distinct from the caller's, and
 * the JAR's identity-only lookup therefore returned the model unchanged and
 * silently: a defect found the first time the two entry points were exercised
 * together. A caller here holding a class of a struct that has since been
 * copied, tagged or transformed is in exactly that position, and the name is
 * the one handle that survives those.
 */
template <class T>
qn::NetworkStruct<T> sn_remove_class(const qn::NetworkStruct<T>& sn, const std::string& name) {
    for (std::size_t r = 0; r < sn.classes.size(); ++r)
        if (sn.classes[r].name == name) return sn_remove_class(sn, r + 1);
    throw InputError("sn_remove_class: the model has no class named '" + name + "'");
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_REMOVE_CLASS_H
