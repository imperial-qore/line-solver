/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_LANG_QN_TAG_CHAIN_H
#define LINE_LANG_QN_TAG_CHAIN_H

/**
 * Port of matlab/src/io/@@ModelAdapter/tagChain.m: the model-to-model transform
 * that isolates ONE job of a chain so that its passage can be observed.
 *
 * WHAT IT DOES. Every class of the chain gets a TAGGED twin: same service, same
 * routing, same reference station, but its own class index. The chain's
 * population is moved by one job -- the class the caller names loses one, its
 * twin gains it -- so the transformed model has the same total population and
 * the same stationary behaviour, with one job now distinguishable from the
 * rest. A twin exists for EVERY class of the chain and not only for the tagged
 * one because the tagged job class-switches as it circulates; without the twins
 * it would leave the tagged block on its first switch and stop being tagged.
 *
 * WHY IT IS A MODEL TRANSFORM AND NOT A POST-PROCESSING STEP. The quantity the
 * response-time CDF needs is the law of the state SEEN BY ONE JOB between its
 * arrival and its departure. The generator of the untagged model has no event
 * that says "this particular job arrived": arrivals of a class are
 * indistinguishable, and the filtration can only split the generator on class-
 * and node-level events. Minting a class with exactly one job in it turns the
 * question into one the filtration can answer.
 *
 * THE POPULATION BOOKKEEPING IS NOT COSMETIC. The twin carries population 1 at
 * the tagged class and 0 at every other class of the chain, and the original
 * class is decremented by 1. Keeping the original population would enlarge the
 * chain by one job and change the very response time being measured.
 *
 * WHAT THE REFERENCE DOES THAT THIS DOES NOT. tagChain.m works on a Network
 * object and re-links it, so it re-derives the class-switch nodes from the
 * linked routing matrix P{r,s}. Here the transform is applied to the refreshed
 * NetworkStruct directly, exactly as fj_mmt does, and the tagged routing block
 * is written into P and into any explicit ClassSwitch matrix before a single
 * refresh_struct() re-derives the chains, the capacities and rt. The two routes
 * agree because a chain is closed under class switching: the block of P (and of
 * the class-switch matrix) over the chain's own classes is already stochastic,
 * so copying it onto the tagged classes needs no rebalancing.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/util/error.h"

namespace line {
namespace qn {

/**
 * The tagged model and the class bookkeeping a caller needs to read it back.
 *
 * `tagged[a]` is the twin of `orig[a]`, both 1-based and in chain order, and
 * `taggedjob` is the one twin that actually holds a job. A caller filtering the
 * generator wants the whole of `tagged`, not just `taggedjob`: the single job
 * moves between the twins as it switches class.
 */
template <class T>
struct TaggedChain {
    NetworkStruct<T> V;
    std::vector<std::size_t> tagged;
    std::vector<std::size_t> orig;
    std::size_t taggedjob = 0;
};

namespace tag_detail {

/**
 * Refuse, by name, every feature whose per-class parameter tables this
 * transform would have to widen and cannot widen meaningfully.
 *
 * The list is not a portability shortfall of the C++ side: tagChain.m widens
 * only the service processes, the output strategies and the capacities, so on
 * any of these models the reference produces a tagged model whose extra class
 * is missing from the very table that decides its dynamics. Silently doing the
 * same here would return a CDF for a model the caller did not describe.
 */
template <class T>
void tag_chain_check(const NetworkStruct<T>& sn) {
    if (!sn.nodeparam.empty())
        throw UnsupportedError(
            "tag_chain: the model has a Cache node, whose item access costs and hit/miss routing "
            "are indexed by class; a tagged twin of a read class has no entry in them");
    if (!sn.transparam.empty())
        throw UnsupportedError(
            "tag_chain: the model has a Transition (SPN) node, whose enabling and firing arcs are "
            "indexed by class and mode; the transform has no rule for the tagged twin's arcs");
    if (sn.has_fork() || !sn.fj.empty())
        throw UnsupportedError(
            "tag_chain: the model has a Fork/Join, where the tagged job becomes several sibling "
            "tasks and the response time is no longer the passage of one job");
    if (!sn.regions.empty())
        throw UnsupportedError(
            "tag_chain: the model has a finite capacity region, whose per-class caps and drop "
            "rules are declared over the original class set and cannot be widened");
    if (!sn.retrialparam.empty())
        throw UnsupportedError(
            "tag_chain: the model has a retrial station, whose orbit process is declared per "
            "class");
    if (!sn.pollingparam.empty())
        throw UnsupportedError(
            "tag_chain: the model has a polling station, whose switchover walks are declared per "
            "buffer, one buffer per class");
    if (!sn.pasparam.empty())
        throw UnsupportedError(
            "tag_chain: the model has a pass-and-swap / order-independent station, whose service "
            "rate function and swap graph are indexed by class");
    for (std::size_t r = 0; r < sn.issignal.size(); ++r)
        if (sn.issignal[r])
            throw UnsupportedError(
                "tag_chain: the model has a G-network signal class, which can annihilate the "
                "tagged job without a departure and leaves its passage undefined");
    for (std::size_t r = 0; r < sn.syncreply.size(); ++r)
        if (sn.syncreply[r] != 0)
            throw UnsupportedError(
                "tag_chain: the model has synchronous call/reply blocking, whose reply class is "
                "declared per calling class and has no tagged twin");
    for (std::size_t i = 0; i < sn.stations.size(); ++i) {
        if (sn.stations[i].cdscaling || sn.stations[i].jdscaling)
            throw UnsupportedError("tag_chain: station '" + sn.stations[i].name +
                                   "' is class dependent; its scaling is a function of the "
                                   "per-class population vector, whose width the tagging changes");
        if (sn.stations[i].svc_rate_fun)
            throw UnsupportedError("tag_chain: station '" + sn.stations[i].name +
                                   "' carries a service rate function over ordered class lists, "
                                   "which has no value at the tagged twin");
    }
}

}  // namespace tag_detail

/**
 * Tag one class of one chain. `chain` and `jobclass` are 1-based; `chain`
 * indexes `sn.inchain`.
 *
 * The chain must be CLOSED. The reference does not check it, and its tagged
 * twin of an open class ends up with a population of 1 on a class that has no
 * population at all -- a model whose state space is not the one the caller
 * asked about. Here that is refused by name instead.
 */
template <class T>
TaggedChain<T> tag_chain(const NetworkStruct<T>& sn, std::size_t chain, std::size_t jobclass,
                         const std::string& suffix = ".tagged") {
    if (chain == 0 || chain > sn.inchain.size())
        throw InputError("tag_chain: chain index out of range");
    const std::vector<std::size_t>& ic = sn.inchain[chain - 1];
    if (jobclass == 0 || jobclass > sn.classes.size())
        throw InputError("tag_chain: class index out of range");
    bool member = false;
    for (std::size_t a = 0; a < ic.size(); ++a)
        if (ic[a] == jobclass) member = true;
    if (!member) throw InputError("tag_chain: the class to tag does not belong to the given chain");
    for (std::size_t a = 0; a < ic.size(); ++a)
        if (!std::isfinite(sn.classes[ic[a] - 1].population))
            throw UnsupportedError("tag_chain: class '" + sn.classes[ic[a] - 1].name +
                                   "' of the chain to tag is open; tagging moves one job out of a "
                                   "finite population and an open chain has none");
    if (sn.classes[jobclass - 1].population < 1.0)
        throw InputError("tag_chain: class '" + sn.classes[jobclass - 1].name +
                         "' has no job to tag");

    tag_detail::tag_chain_check(sn);

    TaggedChain<T> out;
    out.V = sn;
    NetworkStruct<T>& V = out.V;
    const std::size_t I = sn.nodes.size();
    const std::size_t M = sn.stations.size();
    const double inf = std::numeric_limits<double>::infinity();
    const T zero = num_traits<T>::from_int(0);

    out.orig = ic;
    for (std::size_t a = 0; a < ic.size(); ++a) {
        const std::size_t r = ic[a];
        JobClass tc = sn.classes[r - 1];
        tc.name = sn.classes[r - 1].name + suffix;
        tc.population = (r == jobclass) ? 1.0 : 0.0;
        const std::size_t t = V.add_class(tc);
        const std::size_t K2 = V.classes.size();
        out.tagged.push_back(t);
        if (r == jobclass) out.taggedjob = t;

        for (std::size_t i = 1; i <= M; ++i) {
            const Distrib<T>& d = sn.service[i - 1][r - 1];
            Station<T>& st = V.stations[i - 1];
            // An unset capacity is infinite and an unset rule is derived, so
            // padding with those two sentinels leaves every untouched class
            // exactly as the refresh would have found it.
            if (st.classcap.size() < K2) st.classcap.resize(K2, inf);
            if (st.droprule.size() < K2) st.droprule.resize(K2, 0);
            if (!st.schedparam.empty() && st.schedparam.size() < K2)
                st.schedparam.resize(K2, zero);
            if (!st.cdscalingpeak.empty() && st.cdscalingpeak.size() < K2)
                st.cdscalingpeak.resize(K2, zero);
            st.droprule[t - 1] = static_cast<int>(DropStrategy::WAITQ);

            if (d.disabled) {
                // A class the station never serves keeps a twin so the class
                // indexing of the tagged block mirrors the original one, but it
                // is disabled here and given no buffer at all.
                V.set_service(i, t, Distrib<T>::disabled_dist());
                st.classcap[r - 1] = 0.0;
                st.classcap[t - 1] = 0.0;
                if (!st.schedparam.empty()) st.schedparam[t - 1] = zero;
                continue;
            }
            V.set_service(i, t, d);
            // ONE buffer slot for the twin, one fewer for the original: the
            // station still holds the same number of chain jobs, and the tagged
            // job is guaranteed a slot it never has to compete for.
            st.classcap[t - 1] = 1.0;
            if (std::isfinite(st.classcap[r - 1])) st.classcap[r - 1] -= 1.0;
            if (!st.schedparam.empty()) st.schedparam[t - 1] = st.schedparam[r - 1];
            if (!st.cdscalingpeak.empty()) st.cdscalingpeak[t - 1] = st.cdscalingpeak[r - 1];
        }

        // A node whose routing is state dependent decides per CLASS, so the
        // twin has to inherit the strategy and not fall back on the PROB
        // default that an unsized row would give it.
        for (std::size_t i = 0; i < I; ++i) {
            NodeDef& nd = V.nodes[i];
            if (nd.routing.empty()) continue;
            if (nd.routing.size() < K2) nd.routing.resize(K2, RoutingStrategy::PROB);
            nd.routing[t - 1] = sn.nodes[i].routing.size() >= r ? sn.nodes[i].routing[r - 1]
                                                                : RoutingStrategy::PROB;
        }
    }

    // The tagged block of the routing is the chain's own block, copied. Nothing
    // routes between the tagged and the untagged blocks: a tagged job stays
    // tagged for its whole passage, which is the entire point of the transform.
    for (std::size_t x = 0; x < ic.size(); ++x)
        for (std::size_t y = 0; y < ic.size(); ++y)
            for (std::size_t i = 1; i <= I; ++i)
                for (std::size_t j = 1; j <= I; ++j) {
                    const T p = sn.get_route(ic[x], ic[y], i, j);
                    if (num_traits<T>::to_double(p) != 0)
                        V.set_route(out.tagged[x], out.tagged[y], i, j, p);
                }

    // An explicit ClassSwitch matrix is (nclasses x nclasses) and is applied on
    // the way out of the node, so it needs the same block copy. Its rows stay
    // stochastic without rebalancing because a chain is by definition closed
    // under class switching.
    const std::size_t K2 = V.classes.size();
    for (typename std::map<std::size_t, Matrix<T>>::iterator it = V.csmatrix.begin();
         it != V.csmatrix.end(); ++it) {
        const Matrix<T> old = it->second;
        Matrix<T> C(K2, K2, zero);
        for (std::size_t i = 0; i < old.rows() && i < K2; ++i)
            for (std::size_t j = 0; j < old.cols() && j < K2; ++j) C(i, j) = old(i, j);
        for (std::size_t x = 0; x < ic.size(); ++x)
            for (std::size_t y = 0; y < ic.size(); ++y)
                C(out.tagged[x] - 1, out.tagged[y] - 1) = old(ic[x] - 1, ic[y] - 1);
        it->second = C;
    }

    // Last, so that every table above was still read at its original width.
    V.classes[jobclass - 1].population -= 1.0;
    V.refresh_struct();
    return out;
}

}  // namespace qn
}  // namespace line

#endif  // LINE_LANG_QN_TAG_CHAIN_H
