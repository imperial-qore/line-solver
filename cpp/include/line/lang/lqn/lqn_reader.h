/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_LANG_LQN_LQN_READER_H
#define LINE_LANG_LQN_LQN_READER_H

/**
 * .lqnx -> LqnStruct, a port of matlab/src/lang/layered/@@LayeredNetwork/parseXML.m
 * followed by .../getStruct.m.
 *
 * The two MATLAB stages are fused here because the intermediate object graph
 * (Processor / Task / Entry / Activity handles) exists in MATLAB only to be
 * flattened by getStruct, and nothing in this port holds a model object. The
 * ORDER in which the stages walk the document is load-bearing and is preserved
 * exactly, because it fixes the index assignment that every later array is
 * keyed on:
 *
 *   hosts       document order of `<processor>`
 *   tasks       for each processor, document order of its `<task>` descendants
 *   entries     for each task, document order of its `<entry>` descendants
 *   activities  for each task: the `<entry-phase-activities>` activities of each
 *               of its entries, in entry order, THEN its `<task-activities>`
 *               activities
 *
 * That last ordering is not the document order of `<activity>` elements: MATLAB
 * processes all entries of a task before its task-activities block, so an
 * entry-phase activity declared after a task-activities block still receives
 * the lower index. Reproducing it is what makes an index-by-index comparison
 * against a MATLAB dump meaningful.
 *
 * WHAT IS REFUSED. The reader implements the subset of the .lqnx grammar that
 * a layered model needs to reach SolverLN: processors, tasks, entries with
 * phase activities or an activity graph, synchronous and asynchronous calls,
 * forwarding, sequence / AND / OR / loop precedences, replies and open
 * arrivals. Constructs outside it (fan-in and fan-out replication, cache
 * tasks, setup tasks with setup times, service-time distributions declared
 * by histogram) are rejected by name where they would change the answer, and
 * ignored where MATLAB also ignores them.
 */

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "line/api/lqn/lsn_max_multiplicity.h"
#include "line/lang/dist_fitters.h"
#include "line/lang/lqn/lqn_struct.h"
#include "line/util/decimal.h"
#include "line/util/error.h"
#include "line/util/xml.h"

namespace line {
namespace lqn {

namespace detail {

/** Intermediate objects, the C++ stand-in for the MATLAB handle graph. */
template <class T>
struct RawCall {
    std::string dest;
    T mean;
};

template <class T>
struct RawActivity {
    std::string name;
    Distrib<T> hostdem;
    Distrib<T> thinktime;
    std::string bound_to_entry;
    int phase = 1;
    std::size_t task_slot = 0;  ///< 0-based index into the raw task list
    std::vector<RawCall<T>> sync_calls;
    std::vector<RawCall<T>> async_calls;
    /**
     * Routed call groups declared on this activity: the strategy and the target
     * entry NAMES, in declaration order. The member calls themselves are
     * ordinary rows of `sync_calls`; this only records that they are ONE
     * dispatch. Resolved to entry indices by lqn_finalize.
     */
    std::vector<std::pair<lang::RoutingStrategy, std::vector<std::string>>> call_groups;
};

template <class T>
struct RawPrecedence {
    PrecedenceType pretype = PrecedenceType::PRE_SEQ;
    PrecedenceType posttype = PrecedenceType::POST_SEQ;
    std::vector<std::string> preacts;
    std::vector<std::string> postacts;
    std::vector<T> preparams;   ///< PRE_OR branch probabilities, or a PRE_AND quorum
    std::vector<T> postparams;  ///< POST_OR probabilities or POST_LOOP counts
    bool has_quorum = false;
    std::size_t quorum = 0;
};

template <class T>
struct RawEntry {
    std::string name;
    std::size_t task_slot = 0;
    std::vector<std::string> reply_activities;
    bool has_arrival = false;
    Distrib<T> arrival;
    std::vector<std::string> fwd_dest;
    std::vector<T> fwd_prob;
    /** Item entry: cardinality and the popularity pmf over it; 0 = ordinary entry. */
    std::size_t cardinality = 0;
    std::vector<T> popularity;
};

/**
 * One row of an admission constraint, named rather than positional.
 *
 * The operands are entries (on a task) or tasks (on a host); they cannot be
 * resolved to columns until tasksof/entriesof exist, so they are carried by
 * name and resolved at the end of lqn_finalize, as getStruct.m:221-256 does.
 */
template <class T>
struct RawLinConRow {
    std::vector<std::string> names;
    std::vector<T> coeffs;
    T cap;
};

/**
 * One declared server pool, before its operands are resolved to indices.
 *
 * `compatible` names tasks (on a host) or entries (on a task); the names are
 * looked up against the element's own operand list in lqn_finalize, so a pool
 * may be declared before the operand it names.
 */
template <class T>
struct RawServerPool {
    std::string name;
    double count = 1.0;
    T rate;
    std::vector<std::string> compatible;
};

template <class T>
struct RawTask {
    std::string name;
    SchedStrategy sched = SchedStrategy::FCFS;
    double mult = 1.0;
    double repl = 1.0;
    Distrib<T> thinktime;
    std::size_t proc_slot = 0;
    /** fan-out/fan-in as declared: (peer task NAME, value), resolved once indices exist. */
    std::vector<std::pair<std::string, double>> fanout;
    std::vector<std::pair<std::string, double>> fanin;
    std::vector<RawPrecedence<T>> precedences;
    std::vector<RawLinConRow<T>> linconrows;
    Matrix<T> lincon_A;
    std::vector<T> lincon_b;
    /** Cache task: item population, per-list capacity, replacement rule. */
    std::size_t nitems = 0;
    std::vector<int> itemcap;
    ReplacementStrategy replacestrat = ReplacementStrategy::RR;
    /** Setup task: the server shuts down when idle and pays to restart. */
    Distrib<T> setuptime;
    Distrib<T> delayofftime;
    /**
     * Queue-dependent service rates on this task's layer station, over its
     * entries as operands. Empty where not declared; see LqnStruct.
     */
    std::vector<T> lldscaling;
    CdScaling<T> cdscaling;
    std::vector<T> cdscalingpeak;
    CdScaling<T> jdscaling;
    std::vector<T> jdscalingpeak;
    std::vector<RawServerPool<T>> pools;
};

struct RawProc {
    std::string name;
    SchedStrategy sched = SchedStrategy::FCFS;
    double mult = 1.0;
    double repl = 1.0;
    /**
     * `speed-factor` and `quantum`, carried but never read by a solver here.
     *
     * getStruct.m does not read either, so SolverLN cannot see them and this
     * port's LqnStruct has no slot for them. They still have to survive the
     * intermediate model, because SolverLQNS WRITES a .lqnx back out and lqns
     * does honour both: a document declaring `speed-factor="2"` would come back
     * from an unmindful round trip as a processor twice as slow, and nothing in
     * the answer would say so.
     */
    double speed_factor = 1.0;
    double quantum = 0.0;
};

/**
 * The host-demand distribution of an activity, following the MATLAB mapping.
 *
 * MATLAB uses two slightly different ladders, one in the entry-phase branch and
 * one in the task-activities branch: the phase branch sends every scv != 1 to
 * APH, the activity branch splits scv < 1 to APH, scv == 1 to Exp and scv > 1
 * to HyperExp. The activity ladder is the one taken here, because it is the one
 * the writer's own `host-demand-cvsq` round-trips through and it is what
 * `parseXML` applies to every activity outside the phase-indexed form.
 *
 * THE LAW IS FITTED, NOT RELABELLED. This used to build an Exp and then
 * overwrite `type` with APH or HYPEREXP, leaving the Exp's single-entry
 * `params` and its empty (D0, D1) behind a family that indexes three of them --
 * so any consumer reading the parameters read past the end. `dist_scale_rate`
 * does exactly that on the first fixed-point iterate, which aborted line-cli
 * (`vector::operator[]: __n < size()`) on every layered model carrying a
 * non-unit cvsq: lqn_multi_solvers died before its first table.
 */
template <class T>
Distrib<T> host_demand(const std::string& mean_s, const std::string& scv_s) {
    const double mean_d = dbl_from_decimal(mean_s, 0.0);
    if (!(mean_d > 0.0)) return Distrib<T>::immediate();
    const T mean = num_from_decimal<T>(mean_s);
    const double scv_d = dbl_from_decimal(scv_s, 1.0);
    if (!(scv_d > 0.0)) return Distrib<T>::det(mean);
    if (scv_d == 1.0) return Distrib<T>::exp_mean(mean);
    if constexpr (!num_traits<T>::has_transcendental) {
        // Both fits solve a moment condition through a square root, which the
        // exact-arithmetic types have no representation for. Named rather than
        // silently exponential: an activity whose cvsq the document states is
        // not an Exp, and reporting one would be a different model.
        throw UnsupportedError(
            "lqn reader: an activity declares host-demand-cvsq " + scv_s +
            ", whose APH / HyperExp fit needs a square root that exact arithmetic cannot "
            "represent; solve this model in double precision");
    } else {
        const T scv = num_from_decimal<T>(scv_s);
        if (scv_d < 1.0) return lang::aph_fit_mean_scv(mean, scv);
        return lang::hyperexp_fit_mean_scv(mean, scv);
    }
}

/**
 * A `<setup>` / `<delay-off>` mean and SCV back into a distribution.
 *
 * Twin of parseXML.m's `time_from_element` and the JAR's `timeFromElement`:
 * the family is the one the SETTER itself would build, so an SCV of one is the
 * Exp `setSetupTime(mean)` makes and anything else needs a two-moment APH fit.
 * That ladder differs from `host_demand` above, which splits scv > 1 to
 * HyperExp; do not merge them.
 */
template <class T>
Distrib<T> setup_time(const std::string& mean_s, const std::string& scv_s) {
    const double mean_d = dbl_from_decimal(mean_s, 0.0);
    if (!(mean_d > lang::GlobalConstants::FineTol)) return Distrib<T>::immediate();
    const T mean = num_from_decimal<T>(mean_s);
    const double scv_d = dbl_from_decimal(scv_s, 1.0);
    if (std::abs(scv_d - 1.0) <= lang::GlobalConstants::FineTol) return Distrib<T>::exp_mean(mean);
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "lqn reader: a task declares a setup or delay-off scv " + scv_s +
            ", whose APH fit needs a square root that exact arithmetic cannot represent; "
            "solve this model in double precision");
    } else {
        return lang::aph_fit_mean_scv(mean, num_from_decimal<T>(scv_s));
    }
}

/**
 * `<cache replacement="...">` -> ReplacementStrategy.
 *
 * The attribute carries the enum NAME, as MATLAB's writeXML and the JAR's
 * writeXML both emit it; case is not load-bearing, so both spellings are taken.
 * An unknown rule is refused rather than defaulted: the replacement rule is
 * what the hit probability is a function of, so serving CLIMB as RR would
 * answer a different model without saying so.
 */
inline lang::ReplacementStrategy replacement_from_lqnx(const std::string& s) {
    using R = lang::ReplacementStrategy;
    std::string key;
    for (std::size_t i = 0; i < s.size(); ++i) {
        const char c = s[i];
        if (c == ' ' || c == '\t' || c == '\n' || c == '\r') continue;
        key.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
    }
    if (key.empty() || key == "FIFO") return R::FIFO;
    if (key == "RR" || key == "RANDOM") return R::RR;
    if (key == "SFIFO") return R::SFIFO;
    if (key == "LRU") return R::LRU;
    if (key == "HLRU") return R::HLRU;
    if (key == "CLIMB") return R::CLIMB;
    if (key == "QLRU") return R::QLRU;
    throw UnsupportedError("lqn reader: unsupported cache replacement strategy '" + s + "'");
}

/**
 * `<item-entry><access-popularity>` -> the item pmf `lqn.itemproc` carries.
 *
 * The struct keeps the PMF and not the law, so the two discrete families the
 * dialect writes are both reduced here, as `pmf_from_json` does on the JSON
 * side: `DiscreteSampler` states the pmf outright and `Zipf` states (s, n) with
 * p_i = i^-s / H(s,n). A DiscreteSampler written with its support carries 2*card
 * parameters (p then x) and only the first half is the pmf, which is the split
 * the JAR's `readAccessPopularity` makes on the same document.
 */
template <class T>
std::vector<T> popularity_from_lqnx(const xml::Element* item, std::size_t cardinality) {
    std::vector<T> out;
    const std::vector<const xml::Element*> pops = item->child_tags("access-popularity");
    if (pops.empty()) {
        // An item entry with no popularity is still an item entry: the uniform
        // law over its items is what the reference falls back to, and an empty
        // vector would leave the cache with no read at all.
        for (std::size_t k = 0; k < cardinality; ++k)
            out.push_back(num_traits<T>::from_double(1.0 / double(cardinality ? cardinality : 1)));
        return out;
    }
    const xml::Element* pe = pops[0];
    const std::string family = pe->attr("name");
    std::vector<std::string> raw;
    for (const xml::Element* par : pe->child_tags("parameter")) raw.push_back(par->attr("value"));
    if (family == "Zipf") {
        if (raw.size() != 2)
            throw InputError(
                "lqn reader: an access popularity of class Zipf carries exactly two parameters, "
                "the shape s then the item count n, but this one carries " +
                std::to_string(raw.size()));
        const double s = dbl_from_decimal(raw[0], 1.0);
        const std::size_t n = static_cast<std::size_t>(dbl_from_decimal(raw[1], double(cardinality)));
        double h = 0.0;
        for (std::size_t k = 1; k <= n; ++k) h += std::pow(double(k), -s);
        for (std::size_t k = 1; k <= n; ++k)
            out.push_back(num_traits<T>::from_double(std::pow(double(k), -s) / h));
        return out;
    }
    if (family.empty() || family == "DiscreteSampler") {
        const std::size_t n =
            (cardinality > 0 && raw.size() == 2 * cardinality) ? cardinality : raw.size();
        for (std::size_t k = 0; k < n; ++k) out.push_back(num_from_decimal<T>(raw[k]));
        return out;
    }
    throw UnsupportedError(
        "lqn reader: an access popularity is written as '" + family +
        "', and the discrete families the .lqnx dialect carries are DiscreteSampler and Zipf");
}

/**
 * Wire enum name -> RoutingStrategy, the inverse of the writer's
 * `callgroup_to_lqnx`.
 */
inline lang::RoutingStrategy callgroup_from_lqnx(const std::string& name,
                                                 const std::string& act_name) {
    std::string key;
    for (std::size_t i = 0; i < name.size(); ++i) {
        const char c = name[i];
        if (c == ' ' || c == '\t' || c == '\n' || c == '\r') continue;
        key += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    if (key == "RROBIN") return lang::RoutingStrategy::RROBIN;
    if (key == "JSQ") return lang::RoutingStrategy::JSQ;
    throw InputError("lqn reader: activity '" + act_name +
                     "' declares a call group with an unrecognized strategy '" + name +
                     "'; the dialect spells them RROBIN and JSQ");
}

/**
 * Reads the LINE dialect <call-group> children of an activity element.
 *
 * The member calls are ordinary synch-call elements and have already been read,
 * so only the grouping is recorded; issuing them again would double the call
 * rate.
 */
template <class T>
void read_call_groups(const xml::Element* ae, RawActivity<T>& ac) {
    const std::vector<const xml::Element*> groups = ae->child_tags("call-group");
    for (std::size_t g = 0; g < groups.size(); ++g) {
        std::vector<std::string> dests;
        const std::vector<const xml::Element*> de = groups[g]->child_tags("dest");
        for (std::size_t d = 0; d < de.size(); ++d) dests.push_back(de[d]->attr("name"));
        ac.call_groups.push_back(
            std::make_pair(callgroup_from_lqnx(groups[g]->attr("strategy"), ac.name), dests));
    }
}

}  // namespace detail

/**
 * The intermediate model, and the second stage that flattens it.
 *
 * MATLAB reaches the LayeredNetworkStruct along two routes -- parseXML from a
 * .lqnx file, and the Processor/Task/Entry/Activity constructors used directly
 * from a script -- and both end in the SAME getStruct. The split here mirrors
 * that: `LqnModel` is the flat intermediate both routes fill, and
 * `lqn_finalize` is getStruct. It is not a convenience; the .lqnx interchange
 * is LOSSY for models a script can express (it cannot carry a think time on a
 * non-reference task, because lqns rejects one there), so a port that could
 * only read files could not represent every model the reference can.
 */
template <class T>
struct LqnModel {
    std::vector<detail::RawProc> procs;
    std::vector<detail::RawTask<T>> tasks;
    std::vector<detail::RawEntry<T>> entries;
    std::vector<detail::RawActivity<T>> acts;
    /**
     * Admission constraints declared on a HOST, by 0-based processor slot.
     *
     * Kept beside RawProc rather than inside it because RawProc is not a
     * template and the coefficients are T; tasks carry their own rows.
     * `proc_lincon` is the positional (A,b) form, `proc_linconrows` the named
     * one; a host may use either, exactly as a task may.
     */
    std::map<std::size_t, std::vector<detail::RawLinConRow<T>>> proc_linconrows;
    std::map<std::size_t, std::pair<Matrix<T>, std::vector<T>>> proc_lincon;
    /**
     * Queue-dependent service rates and compatibility pools declared on a HOST,
     * by 0-based processor slot.
     *
     * Kept beside RawProc for the same reason as proc_lincon: RawProc is not a
     * template and the scalings are T. A task carries its own in RawTask. The
     * operands of a host are its TASKS, in declaration order.
     */
    std::map<std::size_t, std::vector<T>> proc_lldscaling;
    std::map<std::size_t, CdScaling<T>> proc_cdscaling;
    std::map<std::size_t, std::vector<T>> proc_cdscalingpeak;
    std::map<std::size_t, CdScaling<T>> proc_jdscaling;
    std::map<std::size_t, std::vector<T>> proc_jdscalingpeak;
    std::map<std::size_t, std::vector<detail::RawServerPool<T>>> proc_pools;
};

/** Port of @@LayeredNetwork/getStruct.m: flatten the model into its struct. */
template <class T>
LqnStruct<T> lqn_finalize(const LqnModel<T>& m) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::vector<detail::RawProc>& procs = m.procs;
    const std::vector<detail::RawTask<T>>& tasks = m.tasks;
    const std::vector<detail::RawEntry<T>>& entries = m.entries;
    const std::vector<detail::RawActivity<T>>& acts = m.acts;

    // Stage 2: getStruct
    LqnStruct<T> l;
    l.nhosts = procs.size();
    l.ntasks = tasks.size();
    l.nentries = entries.size();
    l.nacts = acts.size();
    l.hshift = 0;
    l.tshift = l.nhosts;
    l.eshift = l.nhosts + l.ntasks;
    l.ashift = l.eshift + l.nentries;
    l.nidx = l.ashift + l.nacts;
    l.cshift = l.nidx;

    const std::size_t N = l.nidx;
    const std::size_t NT = l.tshift + l.ntasks;
    l.names.assign(N + 1, {});
    l.hashnames.assign(N + 1, {});
    l.type.assign(N + 1, LqnElement::HOST);
    l.parent.assign(N + 1, 0);
    l.sched.assign(NT + 1, SchedStrategy::NONE);
    l.mult.assign(NT + 1, 0.0);
    l.maxmult.assign(NT + 1, 0.0);
    l.repl.assign(NT + 1, 1.0);
    l.lldscaling.assign(NT + 1, {});
    l.cdscaling.assign(NT + 1, CdScaling<T>());
    l.cdscalingpeak.assign(NT + 1, {});
    l.jdscaling.assign(NT + 1, CdScaling<T>());
    l.jdscalingpeak.assign(NT + 1, {});
    l.pools.assign(NT + 1, ServerPools<T>());
    l.isref.assign(NT + 1, false);
    l.iscache.assign(NT + 1, false);
    l.hassetup.assign(NT + 1, false);
    l.nitems.assign(N + 1, 0);
    l.itemcap.assign(NT + 1, {});
    l.replacestrat.assign(NT + 1, ReplacementStrategy::RR);
    l.itemproc.assign(N + 1, {});
    l.setuptime.assign(NT + 1, Distrib<T>::disabled_dist());
    l.delayofftime.assign(NT + 1, Distrib<T>::disabled_dist());
    l.hostdem.assign(N + 1, Distrib<T>::disabled_dist());
    l.think.assign(N + 1, Distrib<T>::disabled_dist());
    l.actthink.assign(N + 1, Distrib<T>::disabled_dist());
    l.has_arrival.assign(N + 1, false);
    l.arrival.assign(N + 1, Distrib<T>::disabled_dist());
    l.tasksof.assign(l.nhosts + 1, {});
    l.entriesof.assign(NT + 1, {});
    l.actsof.assign(l.ashift + 1, {});
    l.callsof.assign(N + 1, {});
    l.precedences.assign(NT + 1, {});
    l.actpretype.assign(N + 1, PrecedenceType::NONE);
    l.actposttype.assign(N + 1, PrecedenceType::NONE);
    l.actquorum.assign(N + 1, 0);
    l.actphase.assign(l.nacts + 1, 1);
    l.graph.resize(N);
    l.taskgraph.resize(NT);
    l.iscaller.resize(N);
    l.issynccaller.resize(N);
    l.isasynccaller.resize(N);

    std::unordered_map<std::string, std::size_t> byhash;

    for (std::size_t p = 0; p < l.nhosts; ++p) {
        const std::size_t idx = p + 1;
        l.sched[idx] = procs[p].sched;
        l.mult[idx] = procs[p].mult;
        l.repl[idx] = procs[p].repl;
        l.names[idx] = procs[p].name;
        l.hashnames[idx] = "P:" + procs[p].name;
        l.type[idx] = LqnElement::HOST;
        // A host keeps its rate dependence in a side map on the model, since
        // RawProc is not a template; the operands are its tasks.
        if (m.proc_lldscaling.count(p)) l.lldscaling[idx] = m.proc_lldscaling.at(p);
        if (m.proc_cdscaling.count(p)) l.cdscaling[idx] = m.proc_cdscaling.at(p);
        if (m.proc_cdscalingpeak.count(p)) l.cdscalingpeak[idx] = m.proc_cdscalingpeak.at(p);
        if (m.proc_jdscaling.count(p)) l.jdscaling[idx] = m.proc_jdscaling.at(p);
        if (m.proc_jdscalingpeak.count(p)) l.jdscalingpeak[idx] = m.proc_jdscalingpeak.at(p);
        byhash[l.hashnames[idx]] = idx;
    }
    for (std::size_t t = 0; t < l.ntasks; ++t) {
        const std::size_t idx = l.tshift + t + 1;
        l.sched[idx] = tasks[t].sched;
        l.hostdem[idx] = Distrib<T>::immediate();
        l.think[idx] = tasks[t].thinktime;
        l.mult[idx] = tasks[t].mult;
        l.repl[idx] = tasks[t].repl;
        l.lldscaling[idx] = tasks[t].lldscaling;
        l.cdscaling[idx] = tasks[t].cdscaling;
        l.cdscalingpeak[idx] = tasks[t].cdscalingpeak;
        l.jdscaling[idx] = tasks[t].jdscaling;
        l.jdscalingpeak[idx] = tasks[t].jdscalingpeak;
        l.names[idx] = tasks[t].name;
        // A cache task is C:, as getStruct.m prefixes it; the reference gives a
        // reference task R: and every other task T:.
        l.nitems[idx] = tasks[t].nitems;
        l.itemcap[idx] = tasks[t].itemcap;
        l.replacestrat[idx] = tasks[t].replacestrat;
        l.iscache[idx] = tasks[t].nitems > 0;
        l.setuptime[idx] = tasks[t].setuptime;
        l.delayofftime[idx] = tasks[t].delayofftime;
        // LQN2QN.m:1252-1275 gates the feature on a setup that is declared, not
        // Immediate, and above tolerance -- a sub-tolerance setup is no setup.
        l.hassetup[idx] = !tasks[t].setuptime.disabled &&
                            num_traits<T>::to_double(tasks[t].setuptime.mean) >
                                lang::GlobalConstants::FineTol;
        l.hashnames[idx] = (l.iscache[idx]                         ? "C:"
                            : l.hassetup[idx]                      ? "T:"
                            : tasks[t].sched == SchedStrategy::REF ? "R:"
                                                                   : "T:") +
                           tasks[t].name;
        l.parent[idx] = tasks[t].proc_slot + 1;
        l.graph.set(idx, l.parent[idx], one);
        l.type[idx] = LqnElement::TASK;
        byhash[l.hashnames[idx]] = idx;
    }
    // a task inherits its host's replication when it declares none of its own
    for (std::size_t t = 0; t < l.ntasks; ++t) {
        const std::size_t tidx = l.tshift + t + 1;
        l.repl[tidx] = std::max(l.repl[tidx], l.repl[l.parent[tidx]]);
    }
    // fan-out/fan-in resolve now that every task carries an element index. A
    // peer that names no task in the model is dropped, as getStruct.m and the
    // Python and JAR struct builders drop it: the declaration is not a call, so
    // an unresolved name removes an edge that never existed.
    {
        std::unordered_map<std::string, std::size_t> task_by_name;
        for (std::size_t t = 0; t < l.ntasks; ++t) {
            const std::size_t tidx = l.tshift + t + 1;
            task_by_name[l.names[tidx]] = tidx;
        }
        for (std::size_t t = 0; t < l.ntasks; ++t) {
            const std::size_t tidx = l.tshift + t + 1;
            for (std::size_t k = 0; k < tasks[t].fanout.size(); ++k) {
                const std::unordered_map<std::string, std::size_t>::const_iterator it =
                    task_by_name.find(tasks[t].fanout[k].first);
                if (it != task_by_name.end())
                    l.fanout[std::make_pair(tidx, it->second)] = tasks[t].fanout[k].second;
            }
            for (std::size_t k = 0; k < tasks[t].fanin.size(); ++k) {
                const std::unordered_map<std::string, std::size_t>::const_iterator it =
                    task_by_name.find(tasks[t].fanin[k].first);
                if (it != task_by_name.end())
                    l.fanin[std::make_pair(tidx, it->second)] = tasks[t].fanin[k].second;
            }
        }
    }
    for (std::size_t p = 1; p <= l.nhosts; ++p)
        for (std::size_t idx = 1; idx <= NT; ++idx)
            if (l.type[idx] == LqnElement::TASK && l.parent[idx] == p) l.tasksof[p].push_back(idx);

    for (std::size_t e = 0; e < l.nentries; ++e) {
        const std::size_t idx = l.eshift + e + 1;
        l.names[idx] = entries[e].name;
        // An item entry is I:, and carries the item population and its pmf
        l.nitems[idx] = entries[e].cardinality;
        l.itemproc[idx] = entries[e].popularity;
        l.hashnames[idx] = (entries[e].cardinality > 0 ? "I:" : "E:") + entries[e].name;
        l.hostdem[idx] = Distrib<T>::immediate();
        l.has_arrival[idx] = entries[e].has_arrival;
        l.arrival[idx] = entries[e].arrival;
        const std::size_t tidx = l.tshift + entries[e].task_slot + 1;
        l.parent[idx] = tidx;
        l.graph.set(tidx, idx, one);
        l.entriesof[tidx].push_back(idx);
        l.type[idx] = LqnElement::ENTRY;
        byhash[l.hashnames[idx]] = idx;
    }
    for (std::size_t a = 0; a < l.nacts; ++a) {
        const std::size_t idx = l.ashift + a + 1;
        l.names[idx] = acts[a].name;
        l.hashnames[idx] = "A:" + acts[a].name;
        l.hostdem[idx] = acts[a].hostdem;
        l.actthink[idx] = acts[a].thinktime;
        const std::size_t tidx = l.tshift + acts[a].task_slot + 1;
        l.parent[idx] = tidx;
        l.actsof[tidx].push_back(idx);
        l.type[idx] = LqnElement::ACTIVITY;
        l.actphase[a + 1] = acts[a].phase;
        byhash[l.hashnames[idx]] = idx;
    }

    auto find_entry = [&](const std::string& name) -> std::size_t {
        auto it = byhash.find("E:" + name);
        if (it != byhash.end()) return it->second;
        it = byhash.find("I:" + name);
        return it == byhash.end() ? 0 : it->second;
    };
    auto find_act = [&](const std::string& name) -> std::size_t {
        auto it = byhash.find("A:" + name);
        return it == byhash.end() ? 0 : it->second;
    };

    // ---- calls, activity binding and precedences, task by task -------------
    std::vector<std::pair<std::size_t, std::size_t>> loop_back_edges;
    std::unordered_map<std::string, std::string> bound_entry_to_act;
    std::size_t cidx = 0;

    auto add_call = [&](std::size_t src, std::size_t dst_e, CallType ct, const T& mean,
                        const std::string& arrow) {
        ++cidx;
        l.callpair_src.push_back(src);
        l.callpair_dst.push_back(dst_e);
        l.calltype.push_back(ct);
        l.callproc_mean.push_back(mean);
        l.callnames.push_back(l.names[src] + arrow + l.names[dst_e]);
        l.callhashnames.push_back(l.hashnames[src] + arrow + l.hashnames[dst_e]);
    };
    // slot 0 of the call arrays is unused, matching the 1-based element arrays
    l.callpair_src.push_back(0);
    l.callpair_dst.push_back(0);
    l.calltype.push_back(CallType::NONE);
    l.callproc_mean.push_back(zero);
    l.callnames.push_back({});
    l.callhashnames.push_back({});

    for (std::size_t t = 0; t < l.ntasks; ++t) {
        const std::size_t tidx = l.tshift + t + 1;
        for (std::size_t a = 0; a < l.nacts; ++a) {
            if (acts[a].task_slot != t) continue;
            const std::size_t aidx = l.ashift + a + 1;

            if (!acts[a].bound_to_entry.empty()) {
                const std::size_t eidx = find_entry(acts[a].bound_to_entry);
                if (eidx > 0) {
                    l.graph.set(eidx, aidx, one);
                    auto it = bound_entry_to_act.find(acts[a].bound_to_entry);
                    if (it != bound_entry_to_act.end())
                        throw InputError("lqn reader: activities '" + it->second + "' and '" +
                                         acts[a].name + "' are both bound to entry '" +
                                         acts[a].bound_to_entry + "'");
                    bound_entry_to_act[acts[a].bound_to_entry] = acts[a].name;
                }
            }

            for (const auto& c : acts[a].sync_calls) {
                const std::size_t te = find_entry(c.dest);
                if (te == 0)
                    throw InputError("lqn reader: activity '" + acts[a].name +
                                     "' calls unknown entry '" + c.dest + "'");
                const std::size_t tt = l.parent[te];
                if (tidx == tt)
                    throw InputError("lqn reader: an entry on a task cannot call another entry on "
                                     "the same task ('" + acts[a].name + "' -> '" + c.dest + "')");
                add_call(aidx, te, CallType::SYNC, c.mean, "=>");
                l.callsof[aidx].push_back(cidx);
                l.iscaller.set(tidx, tt);
                l.iscaller.set(aidx, tt);
                l.iscaller.set(tidx, te);
                l.iscaller.set(aidx, te);
                l.issynccaller.set(tidx, tt);
                l.issynccaller.set(aidx, tt);
                l.issynccaller.set(tidx, te);
                l.issynccaller.set(aidx, te);
                l.taskgraph.set(tidx, tt, one);
                l.graph.set(aidx, te, one);
            }
            for (const auto& c : acts[a].async_calls) {
                const std::size_t te = find_entry(c.dest);
                if (te == 0)
                    throw InputError("lqn reader: activity '" + acts[a].name +
                                     "' has an async call to unknown entry '" + c.dest + "'");
                const std::size_t tt = l.parent[te];
                if (tidx == tt)
                    throw InputError("lqn reader: async self-call from '" + acts[a].name + "'");
                add_call(aidx, te, CallType::ASYNC, c.mean, "->");
                l.callsof[aidx].push_back(cidx);
                l.iscaller.set(aidx, tt);
                l.iscaller.set(aidx, te);
                l.iscaller.set(tidx, tt);
                l.iscaller.set(tidx, te);
                l.isasynccaller.set(tidx, tt);
                l.isasynccaller.set(tidx, te);
                l.isasynccaller.set(aidx, tt);
                l.isasynccaller.set(aidx, te);
                l.taskgraph.set(tidx, tt, one);
                l.graph.set(aidx, te, one);
            }
            // Routed call groups, resolved from target names to entry indices.
            // A group with fewer than two of its targets resolvable is not a
            // dispatch decision and is dropped, which is what the reference
            // does when it filters the group at layer-build time.
            for (const auto& g : acts[a].call_groups) {
                LqnCallGroup grp;
                grp.caller = aidx;
                grp.strategy = g.first;
                for (const std::string& nm : g.second) {
                    const std::size_t te = find_entry(nm);
                    if (te == 0)
                        throw InputError("lqn reader: activity '" + acts[a].name +
                                         "' dispatches a call group to unknown entry '" + nm +
                                         "'");
                    grp.targets.push_back(te);
                }
                if (grp.targets.size() >= 2) l.callgroups.push_back(grp);
            }
        }

        for (const auto& pr : tasks[t].precedences) {
            // Keep the DECLARED precedence: the arc expansion below turns a loop
            // count into a back-edge probability, which method 'srvn.ph' cannot
            // invert -- see LqnStruct::precedences
            {
                LqnPrecedence<T> kept;
                kept.pretype = pr.pretype;
                kept.posttype = pr.posttype;
                kept.preparams = pr.preparams;
                kept.postparams = pr.postparams;
                // An AND-join quorum is declared as has_quorum/quorum, not as a
                // preparam. The workflow composition reads it from pre_params, and
                // a partial join is the one thing it must refuse, so carry it.
                if (pr.pretype == PrecedenceType::PRE_AND && pr.has_quorum &&
                    kept.preparams.empty())
                    kept.preparams.push_back(num_traits<T>::from_double(double(pr.quorum)));
                bool resolved = true;
                for (const std::string& nm : pr.preacts) {
                    const std::size_t ai = find_act(nm);
                    if (ai == 0) { resolved = false; break; }
                    kept.preacts.push_back(ai);
                }
                for (const std::string& nm : pr.postacts) {
                    const std::size_t ai = find_act(nm);
                    if (ai == 0) { resolved = false; break; }
                    kept.postacts.push_back(ai);
                }
                if (resolved) l.precedences[tidx].push_back(kept);
            }
            std::size_t quorum_count = 0;
            if (pr.pretype == PrecedenceType::PRE_AND) {
                if (pr.preacts.empty())
                    throw InputError("lqn reader: PRE_AND precedence with no pre activities in "
                                     "task '" + tasks[t].name + "'");
                quorum_count = (pr.has_quorum && pr.quorum >= 1 && pr.quorum <= pr.preacts.size())
                                   ? pr.quorum
                                   : pr.preacts.size();
            }
            for (std::size_t pa = 0; pa < pr.preacts.size(); ++pa) {
                const std::size_t preaidx = find_act(pr.preacts[pa]);
                if (preaidx == 0)
                    throw InputError("lqn reader: precedence names unknown activity '" +
                                     pr.preacts[pa] + "'");
                switch (pr.posttype) {
                    case PrecedenceType::POST_OR:
                        for (std::size_t po = 0; po < pr.postacts.size(); ++po) {
                            const std::size_t postaidx = find_act(pr.postacts[po]);
                            l.graph.set(preaidx, postaidx, pr.postparams[po]);
                            l.actpretype[preaidx] = pr.pretype;
                            l.actposttype[postaidx] = pr.posttype;
                        }
                        break;
                    case PrecedenceType::POST_AND:
                        for (std::size_t po = 0; po < pr.postacts.size(); ++po) {
                            const std::size_t postaidx = find_act(pr.postacts[po]);
                            l.graph.set(preaidx, postaidx, one);
                            l.actpretype[preaidx] = pr.pretype;
                            l.actposttype[postaidx] = pr.posttype;
                        }
                        break;
                    case PrecedenceType::POST_LOOP: {
                        // postacts = [body..., end]; postparams[0] is the count
                        const T counts = pr.postparams.empty() ? one : pr.postparams[0];
                        const std::size_t enda = pr.postacts.size() - 1;
                        const std::size_t loopentry = find_act(pr.preacts[0]);
                        const std::size_t loopstart = find_act(pr.postacts[0]);
                        const std::size_t loopend = find_act(pr.postacts[enda]);
                        if (counts < one) {
                            l.graph.set(loopentry, loopstart, counts);
                            l.graph.set(loopentry, loopend, T(one - counts));
                            std::size_t cur = loopstart;
                            for (std::size_t po = 1; po + 1 < pr.postacts.size(); ++po) {
                                const std::size_t pi = find_act(pr.postacts[po]);
                                l.graph.set(cur, pi, one);
                                l.actposttype[pi] = pr.posttype;
                                cur = pi;
                            }
                            l.graph.set(cur, loopend, one);
                            l.actposttype[loopstart] = pr.posttype;
                        } else {
                            std::size_t cur = loopentry;
                            for (std::size_t po = 0; po + 1 < pr.postacts.size(); ++po) {
                                const std::size_t pi = find_act(pr.postacts[po]);
                                l.graph.set(cur, pi, one);
                                l.actposttype[pi] = pr.posttype;
                                cur = pi;
                            }
                            loop_back_edges.emplace_back(cur, loopstart);
                            l.graph.set(cur, loopstart, T(one - one / counts));
                            l.graph.set(cur, loopend, T(one / counts));
                        }
                        l.actposttype[loopend] = pr.posttype;
                        break;
                    }
                    default:
                        for (std::size_t po = 0; po < pr.postacts.size(); ++po) {
                            const std::size_t postaidx = find_act(pr.postacts[po]);
                            if (postaidx == 0)
                                throw InputError("lqn reader: precedence names unknown activity '" +
                                                 pr.postacts[po] + "'");
                            l.graph.set(preaidx, postaidx, one);
                            l.actpretype[preaidx] = pr.pretype;
                            l.actposttype[postaidx] = pr.posttype;
                            if (quorum_count > 0) l.actquorum[postaidx] = quorum_count;
                        }
                        break;
                }
            }
        }
    }

    // ---- forwarding calls, after every ordinary call ------------------------
    for (std::size_t e = 0; e < l.nentries; ++e) {
        const std::size_t eidx = l.eshift + e + 1;
        const std::size_t src_t = l.parent[eidx];
        for (std::size_t f = 0; f < entries[e].fwd_dest.size(); ++f) {
            const std::size_t te = find_entry(entries[e].fwd_dest[f]);
            if (te == 0)
                throw InputError("lqn reader: entry '" + entries[e].name +
                                 "' forwards to unknown entry '" + entries[e].fwd_dest[f] + "'");
            if (l.parent[te] == src_t)
                throw InputError("lqn reader: entry '" + entries[e].name +
                                 "' forwards to an entry on the same task");
            add_call(eidx, te, CallType::FWD, entries[e].fwd_prob[f], "~>");
            l.taskgraph.set(src_t, l.parent[te], one);
            l.graph.set(eidx, te, one);
        }
    }
    l.ncalls = cidx;

    // ---- admission constraints, once tasksof/entriesof exist ---------------
    // getStruct.m:221-256. The columns of a host's constraint are its tasks and
    // of a task's its entries, so neither can be resolved before those lists.
    l.lincon_A.assign(l.tshift + l.ntasks + 1, Matrix<T>());
    l.lincon_b.assign(l.tshift + l.ntasks + 1, std::vector<T>());
    for (std::size_t cidx2 = 1; cidx2 <= l.tshift + l.ntasks; ++cidx2) {
        const bool ishost = cidx2 <= l.nhosts;
        const std::vector<std::size_t>& colIdx =
            ishost ? l.tasksof[cidx2] : l.entriesof[cidx2];
        const char* colwhat = ishost ? "tasks on this host" : "entries of this task";
        const Matrix<T>* rawA = nullptr;
        const std::vector<T>* rawb = nullptr;
        const std::vector<detail::RawLinConRow<T>>* rows = nullptr;
        if (ishost) {
            typename std::map<std::size_t, std::vector<detail::RawLinConRow<T>>>::const_iterator
                it = m.proc_linconrows.find(cidx2 - 1);
            if (it != m.proc_linconrows.end()) rows = &it->second;
            typename std::map<std::size_t,
                              std::pair<Matrix<T>, std::vector<T>>>::const_iterator ip =
                m.proc_lincon.find(cidx2 - 1);
            if (ip != m.proc_lincon.end()) {
                rawA = &ip->second.first;
                rawb = &ip->second.second;
            }
        } else {
            const detail::RawTask<T>& rt = tasks[cidx2 - l.tshift - 1];
            rows = &rt.linconrows;
            rawA = &rt.lincon_A;
            rawb = &rt.lincon_b;
        }
        const std::size_t ncols = colIdx.size();
        std::vector<std::vector<T>> Arows;
        std::vector<T> brows;
        if (rawA != nullptr && rawA->rows() > 0) {
            if (rawA->cols() != ncols)
                throw InputError("lqn reader: admission constraint on '" + l.names[cidx2] +
                                 "' has " + std::to_string(rawA->cols()) + " columns but there are " +
                                 std::to_string(ncols) + " " + colwhat);
            for (std::size_t k = 0; k < rawA->rows(); ++k) {
                std::vector<T> row(ncols, zero);
                for (std::size_t j = 0; j < ncols; ++j) row[j] = (*rawA)(k, j);
                Arows.push_back(row);
                brows.push_back(k < rawb->size() ? (*rawb)[k] : zero);
            }
        }
        if (rows != nullptr) {
            for (std::size_t r = 0; r < rows->size(); ++r) {
                std::vector<T> row(ncols, zero);
                for (std::size_t k = 0; k < (*rows)[r].names.size(); ++k) {
                    std::size_t pos = ncols;
                    for (std::size_t j = 0; j < ncols; ++j)
                        if (l.names[colIdx[j]] == (*rows)[r].names[k]) pos = j;
                    if (pos == ncols)
                        throw InputError("lqn reader: admission constraint on '" + l.names[cidx2] +
                                         "' names '" + (*rows)[r].names[k] +
                                         "', which is not one of the " + colwhat);
                    row[pos] = T(row[pos] + (*rows)[r].coeffs[k]);
                }
                Arows.push_back(row);
                brows.push_back((*rows)[r].cap);
            }
        }
        if (Arows.empty()) continue;
        Matrix<T> A(Arows.size(), ncols, zero);
        for (std::size_t k = 0; k < Arows.size(); ++k)
            for (std::size_t j = 0; j < ncols; ++j) A(k, j) = Arows[k][j];
        l.lincon_A[cidx2] = A;
        l.lincon_b[cidx2] = brows;
    }

    // ---- compatibility pools, once tasksof/entriesof exist ------------------
    // A pool names the operands it may serve, and the operands of an element are
    // its tasks (a host) or its entries (a task), so the names cannot become
    // columns before those lists exist. Same staging as the constraints above.
    for (std::size_t pidx = 1; pidx <= l.tshift + l.ntasks; ++pidx) {
        const bool ishost = pidx <= l.nhosts;
        const std::vector<std::size_t>& colIdx =
            ishost ? l.tasksof[pidx] : l.entriesof[pidx];
        const char* colwhat = ishost ? "tasks on this host" : "entries of this task";
        const std::vector<detail::RawServerPool<T>>* raw = nullptr;
        if (ishost) {
            typename std::map<std::size_t,
                              std::vector<detail::RawServerPool<T>>>::const_iterator it =
                m.proc_pools.find(pidx - 1);
            if (it != m.proc_pools.end()) raw = &it->second;
        } else {
            raw = &tasks[pidx - l.tshift - 1].pools;
        }
        if (raw == nullptr || raw->empty()) continue;
        const std::size_t ncols = colIdx.size();
        if (ncols == 0)
            throw InputError("lqn reader: server pools on '" + l.names[pidx] +
                             "' but the element has no operand to serve");
        ServerPools<T> sp;
        sp.compat = Matrix<T>(raw->size(), ncols, zero);
        for (std::size_t t2 = 0; t2 < raw->size(); ++t2) {
            const detail::RawServerPool<T>& rp = (*raw)[t2];
            sp.names.push_back(rp.name);
            sp.counts.push_back(rp.count);
            sp.rates.push_back(rp.rate);
            for (std::size_t k = 0; k < rp.compatible.size(); ++k) {
                bool found = false;
                for (std::size_t j = 0; j < ncols; ++j) {
                    if (l.names[colIdx[j]] == rp.compatible[k]) {
                        sp.compat(t2, j) = one;
                        found = true;
                        break;
                    }
                }
                if (!found)
                    throw InputError("lqn reader: server pool '" + rp.name + "' on '" +
                                     l.names[pidx] + "' names '" + rp.compatible[k] +
                                     "', which is not one of the " + colwhat);
            }
        }
        // A pool nobody can reach is a declaration error, not a zero column to
        // carry: the operand would be served at rate zero and never complete.
        for (std::size_t j = 0; j < ncols; ++j) {
            bool served = false;
            for (std::size_t t2 = 0; t2 < sp.npools(); ++t2)
                if (sp.compat(t2, j) != zero) served = true;
            if (!served)
                throw InputError("lqn reader: '" + l.names[colIdx[j]] + "' on '" + l.names[pidx] +
                                 "' is compatible with no server pool, so it can never be served");
        }
        l.pools[pidx] = sp;
    }

    // ---- every entry must have a bound activity ---------------------------
    for (std::size_t e = 1; e <= l.nentries; ++e) {
        const std::size_t eidx = l.eshift + e;
        bool bound = false;
        for (std::size_t s : l.graph.succ(eidx))
            if (s > l.ashift) bound = true;
        // the message is getStruct.m's, verbatim: this refusal is pinned to the
        // same wording in MATLAB, the JAR and python, so a harness can compare it
        if (!bound) throw InputError("An entry does not have any boundTo activity.");
    }

    // ---- a replying activity must have no PHASE 1 successor ---------------
    // getStruct.m's guard, absent from this port until 2026-08-15. An activity
    // that replies ends phase 1 of its entry, so a successor still marked phase
    // 1 is a graph the struct cannot represent: the reply would be read as the
    // end of the entry and the tail served as though it did not exist. A phase
    // 2 successor is the legitimate case, post-reply processing.
    for (std::size_t e = 0; e < entries.size(); ++e) {
        for (const std::string& rname : entries[e].reply_activities) {
            const std::size_t aidx = find_act(rname);
            if (aidx == 0) continue;
            for (std::size_t succ : l.graph.succ(aidx)) {
                if (succ <= l.ashift) continue;
                if (l.actphase[succ - l.ashift] == 1)
                    throw InputError("Unsupported replyTo in non-terminal activity.");
            }
        }
    }

    // infinite-server multiplicity correction rationale: see _kb/04-networkstruct.md (cpp port notes)
    for (std::size_t tidx = 1; tidx <= NT; ++tidx) {
        if (l.sched[tidx] != SchedStrategy::INF) continue;
        if (l.type[tidx] != LqnElement::TASK) continue;
        double s = 0.0;
        for (std::size_t c = 1; c <= NT; ++c)
            if (l.taskgraph.get(c, tidx) != zero) s += l.mult[c];
        l.mult[tidx] = s;
    }

    for (std::size_t idx = 1; idx <= NT; ++idx) l.isref[idx] = l.sched[idx] == SchedStrategy::REF;

    // ---- the dag ----------------------------------------------------------
    l.dag = l.graph;
    for (std::size_t i = 1; i <= N; ++i) {
        if (l.type[i] != LqnElement::TASK || l.isref[i]) continue;
        std::vector<std::size_t> to_flip;
        for (const auto& e : l.dag.row[i])
            if (l.type[e.first] == LqnElement::ENTRY && e.second != zero)
                to_flip.push_back(e.first);
        for (std::size_t j : to_flip) {
            l.dag.erase(i, j);
            l.dag.set(j, i, one);
        }
    }
    for (const auto& be : loop_back_edges) l.dag.erase(be.first, be.second);

    // ---- entry-to-activity reachability -----------------------------------
    for (std::size_t e = 1; e <= l.nentries; ++e) {
        const std::size_t eidx = l.eshift + e;
        const std::size_t tidx = l.parent[eidx];
        std::vector<bool> visited(N + 1, false);
        std::vector<std::size_t> stack{eidx};
        visited[eidx] = true;
        while (!stack.empty()) {
            const std::size_t v = stack.back();
            stack.pop_back();
            for (std::size_t w : l.graph.succ(v))
                if (!visited[w]) {
                    visited[w] = true;
                    stack.push_back(w);
                }
        }
        std::vector<std::size_t> found;
        for (std::size_t i = 1; i <= N; ++i)
            if (visited[i] && l.type[i] == LqnElement::ACTIVITY && l.parent[i] == tidx)
                found.push_back(i);
        l.actsof[eidx] = found;
    }

    // ---- sustainable multiplicities ---------------------------------------
    {
        lsn::LsnInput<double> in;
        in.dag = Matrix<double>(N, N, 0.0);
        for (std::size_t i = 1; i <= N; ++i)
            for (const auto& ed : l.dag.row[i])
                if (ed.second != zero) in.dag(i - 1, ed.first - 1) = 1.0;
        in.mult.resize(N);
        in.type.resize(N);
        in.isref.assign(N, false);
        in.entry_has_arrival.assign(N, false);
        // A SETUP TASK KEEPS ITS SPARE CAPACITY. lsn_max_multiplicity.m:69-72
        // exempts it from the min against its inflow, because the servers a
        // caller cannot keep busy are exactly the ones that power down and pay
        // the setup -- trimming them away deletes the effect being modelled.
        // Leaving this vector empty silently built every setup layer with one
        // server, which is the whole layer, not a detail of it.
        in.hassetup.assign(N, false);
        for (std::size_t i = 1; i <= N; ++i) {
            const double m = i <= NT ? l.mult[i] : std::numeric_limits<double>::infinity();
            in.mult[i - 1] = std::isinf(m) ? lsn::Multiplicity<double>::inf()
                                           : lsn::Multiplicity<double>::finite(m);
            switch (l.type[i]) {
                case LqnElement::HOST: in.type[i - 1] = lsn::LsnElementType::HOST; break;
                case LqnElement::TASK: in.type[i - 1] = lsn::LsnElementType::TASK; break;
                case LqnElement::ENTRY: in.type[i - 1] = lsn::LsnElementType::ENTRY; break;
                default: in.type[i - 1] = lsn::LsnElementType::ACTIVITY; break;
            }
            if (i <= NT) in.isref[i - 1] = l.isref[i];
            if (i <= NT) in.hassetup[i - 1] = l.hassetup[i];
            in.entry_has_arrival[i - 1] = l.has_arrival[i];
        }
        const std::vector<lsn::Multiplicity<double>> mm = lsn::lsn_max_multiplicity(in);
        for (std::size_t i = 1; i <= NT; ++i)
            l.maxmult[i] = mm[i - 1].infinite ? std::numeric_limits<double>::infinity()
                                              : mm[i - 1].value;
    }

    // ---- an entry must not be called both synchronously and asynchronously --
    for (std::size_t e = 1; e <= l.nentries; ++e) {
        const std::size_t eidx = l.eshift + e;
        bool sync = false, async = false;
        for (std::size_t c = 1; c <= l.ncalls; ++c) {
            if (l.callpair_dst[c] != eidx) continue;
            if (l.calltype[c] == CallType::SYNC) sync = true;
            if (l.calltype[c] == CallType::ASYNC) async = true;
        }
        if (sync && async)
            throw InputError("lqn reader: entry '" + l.names[eidx] +
                             "' is called both synchronously and asynchronously");
    }

    return l;
}

/**
 * Read a .lqnx model into the INTERMEDIATE form, before getStruct flattens it.
 *
 * `read_lqnx` is this followed by `lqn_finalize`, and is what a solver wants.
 * The intermediate form is what a WRITER wants: `write_lqnx` (lqn_writer.h)
 * emits declarations -- precedence blocks, reply entries, fan-out -- that the
 * struct records only in flattened form, so a round trip through the struct
 * alone could not reproduce the document. SolverLQNS hands the file it writes
 * to an external binary, which will reject a document whose precedence blocks
 * were guessed, so the declarative form is not optional there.
 *
 * @param path  file to read
 * @return      the intermediate model, exactly as the document declares it
 */
namespace detail {

/** Renders a number as the MATLAB, JAR and Python readers do, so messages agree. */
inline std::string fmt_num(double v) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%g", v);
    return std::string(buf);
}

/** Reads a numeric attribute; NaN when the text is present but not a number. */
inline double attr_num(const std::string& s, double dflt) {
    if (s.empty()) return dflt;
    try {
        return std::stod(s);
    } catch (const std::exception&) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

/** Case-insensitive equality, for the scheduling attribute. */
inline bool iequals(const std::string& a, const std::string& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (std::tolower(static_cast<unsigned char>(a[i])) !=
            std::tolower(static_cast<unsigned char>(b[i])))
            return false;
    return true;
}

/**
 * Reject a structurally inconsistent LQN document.
 *
 * Run on the parsed document before any object is built, so that a defective
 * input is named at its source instead of surfacing as a downstream failure.
 * The same checks, in the same order and with the same messages, are applied by
 * the MATLAB, JAR and Python readers.
 *
 * @param doc  root element of the parsed document
 */
inline void validate_input_model(const xml::Element& doc) {
    const double tol = 1e-6;
    std::vector<std::string> proc_names;
    std::vector<std::string> task_names;
    std::vector<std::string> entry_names;
    std::vector<std::string> entry_owner;  // task owning entry_names[k]
    std::vector<bool> is_ref_entry;
    std::vector<std::string> call_dests;
    std::vector<std::string> reply_entries;
    bool has_ref_task = false;
    bool has_open_arrival = false;

    for (const xml::Element* pe : doc.by_tag("processor")) {
        const std::string proc_name = pe->attr("name");
        if (std::find(proc_names.begin(), proc_names.end(), proc_name) != proc_names.end())
            throw InputError("Duplicate processor name \"" + proc_name + "\".");
        proc_names.push_back(proc_name);

        for (const xml::Element* te : pe->by_tag("task")) {
            const std::string task_name = te->attr("name");
            if (std::find(task_names.begin(), task_names.end(), task_name) != task_names.end())
                throw InputError("Duplicate task name \"" + task_name + "\".");
            task_names.push_back(task_name);
            const bool is_ref = iequals(te->attr("scheduling"), "ref");
            has_ref_task = has_ref_task || is_ref;

            const std::vector<const xml::Element*> entry_els = te->by_tag("entry");
            if (entry_els.empty())
                throw InputError("Task \"" + task_name + "\" has no entries.");
            for (const xml::Element* ee : entry_els) {
                const std::string entry_name = ee->attr("name");
                if (std::find(entry_names.begin(), entry_names.end(), entry_name) !=
                    entry_names.end())
                    throw InputError("Duplicate entry name \"" + entry_name + "\".");
                entry_names.push_back(entry_name);
                entry_owner.push_back(task_name);
                is_ref_entry.push_back(is_ref);

                const double arrival_rate =
                    attr_num(ee->attr("open-arrival-rate"), std::numeric_limits<double>::quiet_NaN());
                if (arrival_rate > 0.0) {
                    has_open_arrival = true;
                    if (is_ref)
                        throw InputError("Entry \"" + entry_name + "\" belongs to reference task \"" +
                                         task_name + "\" and cannot have open arrivals.");
                }

                const std::vector<const xml::Element*> fwd_els = ee->by_tag("forwarding");
                if (is_ref && !fwd_els.empty())
                    throw InputError("Entry \"" + entry_name + "\" belongs to reference task \"" +
                                     task_name + "\" and cannot forward requests.");
                double fwd_total = 0.0;
                for (const xml::Element* fe : fwd_els) {
                    const double prob = attr_num(fe->attr("prob"), 1.0);
                    if (std::isnan(prob) || prob < 0.0 || prob > 1.0)
                        throw InputError("Forwarding from entry \"" + entry_name + "\" to entry \"" +
                                         fe->attr("dest") + "\" has an invalid probability of " +
                                         fmt_num(prob) + ".");
                    fwd_total += prob;
                }
                if (fwd_total > 1.0 + tol)
                    throw InputError("Entry \"" + entry_name +
                                     "\" has a total forwarding probability of " +
                                     fmt_num(fwd_total) + ".");
            }

            // activity names are unique within their task; a name under a pre or post list is a reference, not a declaration
            std::vector<std::string> act_names;
            for (const xml::Element* ae : te->by_tag("activity")) {
                if (ae->parent == nullptr) continue;
                if (ae->parent->name != "task-activities" &&
                    ae->parent->name != "entry-phase-activities")
                    continue;
                const std::string act_name = ae->attr("name");
                if (std::find(act_names.begin(), act_names.end(), act_name) != act_names.end())
                    throw InputError("Duplicate activity name \"" + act_name + "\" in task \"" +
                                     task_name + "\".");
                act_names.push_back(act_name);
            }

            for (const xml::Element* ce : te->by_tag("synch-call"))
                call_dests.push_back(ce->attr("dest"));
            for (const xml::Element* ce : te->by_tag("asynch-call"))
                call_dests.push_back(ce->attr("dest"));
            for (const xml::Element* fe : te->by_tag("forwarding"))
                call_dests.push_back(fe->attr("dest"));

            for (const xml::Element* oe : te->by_tag("post-OR")) {
                double branch_total = 0.0;
                for (const xml::Element* be : oe->by_tag("activity")) {
                    const double prob = attr_num(be->attr("prob"), 1.0);
                    if (std::isnan(prob) || prob < 0.0 || prob > 1.0)
                        throw InputError("Activity \"" + be->attr("name") + "\" in task \"" +
                                         task_name + "\" has an invalid branch probability of " +
                                         fmt_num(prob) + ".");
                    branch_total += prob;
                }
                if (std::fabs(branch_total - 1.0) > tol)
                    throw InputError("Branch probabilities of an OR-fork in task \"" + task_name +
                                     "\" sum to " + fmt_num(branch_total) + " instead of 1.");
            }

            for (const xml::Element* re : te->by_tag("reply-entry"))
                reply_entries.push_back(re->attr("name"));
        }
    }

    for (const std::string& dest : call_dests) {
        const std::vector<std::string>::const_iterator it =
            std::find(entry_names.begin(), entry_names.end(), dest);
        if (it == entry_names.end()) continue;
        const std::size_t idx = static_cast<std::size_t>(it - entry_names.begin());
        if (is_ref_entry[idx])
            throw InputError("Entry \"" + entry_names[idx] + "\" belongs to reference task \"" +
                             entry_owner[idx] + "\" and cannot receive requests.");
    }

    for (const std::string& reply_name : reply_entries) {
        const std::vector<std::string>::const_iterator it =
            std::find(entry_names.begin(), entry_names.end(), reply_name);
        if (it == entry_names.end()) continue;
        const std::size_t idx = static_cast<std::size_t>(it - entry_names.begin());
        if (is_ref_entry[idx])
            throw InputError("Entry \"" + entry_names[idx] + "\" belongs to reference task \"" +
                             entry_owner[idx] + "\" and cannot be replied to.");
    }

    if (!has_ref_task && !has_open_arrival)
        throw InputError("The model has no reference task and no open arrivals.");
}

}  // namespace detail

template <class T>
LqnModel<T> read_lqnx_model(const std::string& path) {
    const T one = num_traits<T>::from_int(1);
    LqnModel<T> m;
    std::vector<detail::RawProc>& procs = m.procs;
    std::vector<detail::RawTask<T>>& tasks = m.tasks;
    std::vector<detail::RawEntry<T>>& entries = m.entries;
    std::vector<detail::RawActivity<T>>& acts = m.acts;


    std::unique_ptr<xml::Element> doc = xml::parse_file(path);
    detail::validate_input_model(*doc);

    // Stage 1: parseXML

    const std::vector<const xml::Element*> proc_els = doc->by_tag("processor");
    for (const xml::Element* pe : proc_els) {
        detail::RawProc pr;
        pr.name = pe->attr("name");
        const std::string psched = pe->attr("scheduling");
        pr.sched = lang::sched_from_lqnx(psched.empty() ? std::string("fcfs") : psched);
        pr.repl = dbl_from_decimal(pe->attr("replication"), 1.0);
        pr.speed_factor = dbl_from_decimal(pe->attr("speed-factor"), 1.0);
        pr.quantum = dbl_from_decimal(pe->attr("quantum"), 0.0);
        if (pr.sched == SchedStrategy::INF) {
            // A finite multiplicity on an inf-scheduled processor is discarded,
            // as in MATLAB, which warns and overrides it.
            pr.mult = std::numeric_limits<double>::infinity();
        } else {
            pr.mult = dbl_from_decimal(pe->attr("multiplicity"), 1.0);
        }
        const std::size_t proc_slot = procs.size();
        procs.push_back(pr);

        for (const xml::Element* te : pe->by_tag("task")) {
            detail::RawTask<T> tk;
            tk.name = te->attr("name");
            const std::string tsched = te->attr("scheduling");
            tk.sched = lang::sched_from_lqnx(tsched.empty() ? std::string("fcfs") : tsched);
            tk.repl = dbl_from_decimal(te->attr("replication"), 1.0);
            if (tk.sched == SchedStrategy::INF) {
                tk.mult = std::numeric_limits<double>::infinity();
            } else {
                tk.mult = dbl_from_decimal(te->attr("multiplicity"), 1.0);
            }
            const std::string think_s = te->attr("think-time");
            const double think_d = dbl_from_decimal(think_s, 0.0);
            tk.thinktime = think_d > 0.0 ? Distrib<T>::exp_mean(num_from_decimal<T>(think_s))
                                         : Distrib<T>::immediate();
            tk.proc_slot = proc_slot;
            // fan-out names a callee task, fan-in a caller task; both carry the
            // count of peer replicas one replica of THIS task addresses. Stored
            // by name because the callee's element index does not exist yet.
            for (const xml::Element* fe : te->by_tag("fan-out"))
                tk.fanout.push_back(
                    std::make_pair(fe->attr("dest"), dbl_from_decimal(fe->attr("value"), 1.0)));
            for (const xml::Element* fe : te->by_tag("fan-in"))
                tk.fanin.push_back(
                    std::make_pair(fe->attr("source"), dbl_from_decimal(fe->attr("value"), 1.0)));
            // <setup> and <delay-off> are a LINE extension to the schema that
            // MATLAB, the JAR and Python all write and read. Skipping them here
            // dropped the cold start entirely: lqn_setup came back with the bare
            // host demand at the entry (E2 RespT 0.3333 against 1) because
            // `hassetup` stayed false and `setup_charge` returned 0.
            // child_tags and not by_tag: parseXML.m reads them as DIRECT
            // children, and a descendant search would read a nested task's.
            for (const xml::Element* se : te->child_tags("setup"))
                tk.setuptime = detail::setup_time<T>(se->attr("mean"), se->attr("scv"));
            for (const xml::Element* se : te->child_tags("delay-off"))
                tk.delayofftime = detail::setup_time<T>(se->attr("mean"), se->attr("scv"));
            // <cache> is the LINE extension that makes this a CacheTask, exactly
            // as in parseXML.m and the JAR's readXML: `items` is the item
            // population, each <level> one cache list. IGNORING IT DOES NOT
            // DEGRADE THE MODEL GRACEFULLY -- with `nitems` left at 0 the layer
            // builder never marks the host layer a cache layer, so no Cache node
            // is added and the hit/miss branch keeps the even split `link()`
            // offers, independent of capacity, item count and replacement rule
            // (lcq_threehosts came back with hit 0.5 against 0.48331).
            for (const xml::Element* ce : te->child_tags("cache")) {
                const std::string items_s = ce->attr("items");
                const double items_d = dbl_from_decimal(items_s, 1.0);
                tk.nitems = items_d > 0.0 ? static_cast<std::size_t>(items_d) : 1;
                tk.replacestrat = detail::replacement_from_lqnx(ce->attr("replacement"));
                tk.itemcap.clear();
                for (const xml::Element* le : ce->child_tags("level")) {
                    const double cap_d = dbl_from_decimal(le->attr("capacity"), 1.0);
                    tk.itemcap.push_back(cap_d > 0.0 ? static_cast<int>(cap_d) : 1);
                }
                // A <cache> with no <level> is a single list of one item, the
                // same reading the JAR takes; the capacity is what the cache
                // HOLDS, so an empty vector would make it hold nothing.
                if (tk.itemcap.empty()) tk.itemcap.push_back(1);
            }
            const std::size_t task_slot = tasks.size();
            tasks.push_back(tk);

            for (const xml::Element* ee : te->by_tag("entry")) {
                detail::RawEntry<T> en;
                en.name = ee->attr("name");
                en.task_slot = task_slot;
                const std::string arr_s = ee->attr("open-arrival-rate");
                if (!arr_s.empty()) {
                    const double rate_d = dbl_from_decimal(arr_s, 0.0);
                    if (rate_d > 0.0) {
                        en.has_arrival = true;
                        en.arrival = Distrib<T>::exp_rate(num_from_decimal<T>(arr_s));
                    }
                }
                for (const xml::Element* fe : ee->by_tag("forwarding")) {
                    en.fwd_dest.push_back(fe->attr("dest"));
                    const std::string ps = fe->attr("prob");
                    en.fwd_prob.push_back(ps.empty() ? one : num_from_decimal<T>(ps));
                }
                // <item-entry> is what makes this an ItemEntry, the twin of the
                // <cache> test on the task above; the cache layer reads the item
                // pmf off `lqn.itemproc`, so dropping it leaves the read
                // unpopulated even when the cache task itself was recognised.
                for (const xml::Element* ie2 : ee->child_tags("item-entry")) {
                    const double card_d = dbl_from_decimal(ie2->attr("cardinality"), 1.0);
                    en.cardinality = card_d > 0.0 ? static_cast<std::size_t>(card_d) : 1;
                    en.popularity = detail::popularity_from_lqnx<T>(ie2, en.cardinality);
                }
                const std::size_t entry_slot = entries.size();
                entries.push_back(en);

                const std::vector<const xml::Element*> epa = ee->by_tag("entry-phase-activities");
                if (!epa.empty()) {
                    // phase-indexed activity-name rationale: see _kb/04-networkstruct.md (cpp port notes)
                    std::map<int, std::string> by_phase;
                    for (const xml::Element* ae : epa[0]->by_tag("activity")) {
                        const int phase = static_cast<int>(dbl_from_decimal(ae->attr("phase"), 1.0));
                        detail::RawActivity<T> ac;
                        ac.name = ae->attr("name");
                        ac.hostdem = detail::host_demand<T>(ae->attr("host-demand-mean"),
                                                            ae->attr("host-demand-cvsq"));
                        const std::string att = ae->attr("think-time");
                        const double att_d = dbl_from_decimal(att, 0.0);
                        ac.thinktime = att_d > 0.0
                                           ? Distrib<T>::exp_mean(num_from_decimal<T>(att))
                                           : Distrib<T>::immediate();
                        ac.bound_to_entry = phase == 1 ? entries[entry_slot].name : std::string();
                        ac.phase = phase;
                        ac.task_slot = task_slot;
                        for (const xml::Element* ce : ae->by_tag("synch-call"))
                            ac.sync_calls.push_back(
                                {ce->attr("dest"), num_from_decimal<T>(ce->attr("calls-mean"))});
                        for (const xml::Element* ce : ae->by_tag("asynch-call"))
                            ac.async_calls.push_back(
                                {ce->attr("dest"), num_from_decimal<T>(ce->attr("calls-mean"))});
                        detail::read_call_groups<T>(ae, ac);
                        by_phase[phase] = ac.name;
                        acts.push_back(ac);
                    }
                    // implicit precedence between consecutive phases
                    for (auto it = by_phase.begin(); it != by_phase.end(); ++it) {
                        auto nx = std::next(it);
                        if (nx == by_phase.end()) break;
                        detail::RawPrecedence<T> pr;
                        pr.pretype = PrecedenceType::PRE_SEQ;
                        pr.posttype = PrecedenceType::POST_SEQ;
                        pr.preacts.push_back(it->second);
                        pr.postacts.push_back(nx->second);
                        tasks[task_slot].precedences.push_back(pr);
                    }
                    if (!by_phase.empty() && by_phase.count(1))
                        entries[entry_slot].reply_activities.push_back(by_phase[1]);
                }
            }

            const std::vector<const xml::Element*> tal = te->by_tag("task-activities");
            if (!tal.empty()) {
                const xml::Element* ta = tal[0];
                for (const xml::Element* ae : ta->by_tag("activity")) {
                    // descendant-search scope rationale: see _kb/04-networkstruct.md (cpp port notes)
                    if (ae->parent != ta) continue;
                    detail::RawActivity<T> ac;
                    ac.name = ae->attr("name");
                    ac.hostdem = detail::host_demand<T>(ae->attr("host-demand-mean"),
                                                        ae->attr("host-demand-cvsq"));
                    const std::string att = ae->attr("think-time");
                    const double att_d = dbl_from_decimal(att, 0.0);
                    ac.thinktime = att_d > 0.0 ? Distrib<T>::exp_mean(num_from_decimal<T>(att))
                                               : Distrib<T>::immediate();
                    ac.bound_to_entry = ae->attr("bound-to-entry");
                    ac.phase = 1;
                    ac.task_slot = task_slot;
                    for (const xml::Element* ce : ae->by_tag("synch-call"))
                        ac.sync_calls.push_back(
                            {ce->attr("dest"), num_from_decimal<T>(ce->attr("calls-mean"))});
                    for (const xml::Element* ce : ae->by_tag("asynch-call"))
                        ac.async_calls.push_back(
                            {ce->attr("dest"), num_from_decimal<T>(ce->attr("calls-mean"))});
                    detail::read_call_groups<T>(ae, ac);
                    acts.push_back(ac);
                }

                for (const xml::Element* pe2 : ta->by_tag("precedence")) {
                    detail::RawPrecedence<T> pr;
                    const xml::Element* pre = nullptr;
                    if (!pe2->child_tags("pre").empty()) {
                        pre = pe2->child_tags("pre")[0];
                        pr.pretype = PrecedenceType::PRE_SEQ;
                    } else if (!pe2->child_tags("pre-AND").empty()) {
                        pre = pe2->child_tags("pre-AND")[0];
                        pr.pretype = PrecedenceType::PRE_AND;
                    } else if (!pe2->child_tags("pre-OR").empty()) {
                        pre = pe2->child_tags("pre-OR")[0];
                        pr.pretype = PrecedenceType::PRE_OR;
                    } else {
                        throw InputError("lqn reader: <precedence> without a pre element");
                    }
                    for (const xml::Element* ae : pre->by_tag("activity")) {
                        pr.preacts.push_back(ae->attr("name"));
                        if (pr.pretype == PrecedenceType::PRE_OR)
                            pr.preparams.push_back(num_from_decimal<T>(ae->attr("prob")));
                    }
                    if (pr.pretype == PrecedenceType::PRE_SEQ && pr.preacts.size() > 1)
                        pr.preacts.resize(1);
                    if (pr.pretype == PrecedenceType::PRE_AND) {
                        const std::string q = pre->attr("quorum");
                        if (!q.empty()) {
                            pr.has_quorum = true;
                            pr.quorum = static_cast<std::size_t>(dbl_from_decimal(q, 0.0) + 0.5);
                        }
                    }

                    const xml::Element* post = nullptr;
                    if (!pe2->child_tags("post").empty()) {
                        post = pe2->child_tags("post")[0];
                        pr.posttype = PrecedenceType::POST_SEQ;
                    } else if (!pe2->child_tags("post-AND").empty()) {
                        post = pe2->child_tags("post-AND")[0];
                        pr.posttype = PrecedenceType::POST_AND;
                    } else if (!pe2->child_tags("post-OR").empty()) {
                        post = pe2->child_tags("post-OR")[0];
                        pr.posttype = PrecedenceType::POST_OR;
                    } else if (!pe2->child_tags("post-LOOP").empty()) {
                        post = pe2->child_tags("post-LOOP")[0];
                        pr.posttype = PrecedenceType::POST_LOOP;
                    } else if (!pe2->child_tags("post-CACHE").empty()) {
                        post = pe2->child_tags("post-CACHE")[0];
                        pr.posttype = PrecedenceType::POST_CACHE;
                    } else {
                        // minOccurs="0" on the post choice in lqn-core.xsd: a
                        // precedence carrying only a pre element declares a
                        // TERMINAL activity and no successor, so it contributes
                        // no edge and is dropped rather than refused.
                        continue;
                    }
                    for (const xml::Element* ae : post->by_tag("activity")) {
                        pr.postacts.push_back(ae->attr("name"));
                        if (pr.posttype == PrecedenceType::POST_OR)
                            pr.postparams.push_back(num_from_decimal<T>(ae->attr("prob")));
                        if (pr.posttype == PrecedenceType::POST_LOOP)
                            pr.postparams.push_back(num_from_decimal<T>(ae->attr("count")));
                    }
                    if (pr.posttype == PrecedenceType::POST_CACHE) {
                        // THE BRANCH IS NAMED WHERE THE WRITER NAMED IT. `hit`
                        // then `miss` is what the builder and SolverLN key on, so
                        // a file whose activities are in the other order must be
                        // reordered rather than read positionally. A file written
                        // before `cache-result` existed carries no attribute, and
                        // there document order IS the order, as layered.py:4204
                        // also falls back to.
                        std::vector<std::string> hit, miss, unlabelled;
                        const std::vector<const xml::Element*> aes = post->by_tag("activity");
                        for (std::size_t k = 0; k < aes.size(); ++k) {
                            const std::string res = aes[k]->attr("cache-result");
                            if (res == "hit")
                                hit.push_back(pr.postacts[k]);
                            else if (res == "miss")
                                miss.push_back(pr.postacts[k]);
                            else
                                unlabelled.push_back(pr.postacts[k]);
                        }
                        if (!hit.empty() || !miss.empty()) {
                            pr.postacts = hit;
                            pr.postacts.insert(pr.postacts.end(), miss.begin(), miss.end());
                            pr.postacts.insert(pr.postacts.end(), unlabelled.begin(),
                                               unlabelled.end());
                        }
                        if (pr.postacts.size() < 2)
                            throw InputError(
                                "lqn reader: a <post-CACHE> branches on a hit and a miss, so it "
                                "names two activities");
                    }
                    if (pr.posttype == PrecedenceType::POST_LOOP)
                        pr.postacts.push_back(post->attr("end"));
                    tasks[task_slot].precedences.push_back(pr);
                }

                for (const xml::Element* re : ta->by_tag("reply-entry")) {
                    const std::string ename = re->attr("name");
                    std::size_t slot = entries.size();
                    for (std::size_t s = 0; s < entries.size(); ++s)
                        if (entries[s].name == ename) {
                            slot = s;
                            break;
                        }
                    if (slot == entries.size())
                        throw InputError("lqn reader: <reply-entry> names unknown entry '" + ename +
                                         "'");
                    for (const xml::Element* ra : re->by_tag("reply-activity"))
                        entries[slot].reply_activities.push_back(ra->attr("name"));
                }
            }
        }
    }

    return m;
}

/**
 * Read a .lqnx model.
 *
 * @param path  file to read
 * @return      the flattened struct SolverLN consumes
 */
template <class T>
LqnStruct<T> read_lqnx(const std::string& path) {
    return lqn_finalize(read_lqnx_model<T>(path));
}


}  // namespace lqn
}  // namespace line

#endif  // LINE_LANG_LQN_LQN_READER_H
