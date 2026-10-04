/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_LINEMODEL_WRITER_H
#define LINE_IO_LINEMODEL_WRITER_H

/**
 * @file
 * @ingroup line_io
 * model.json writers for a LayeredNetwork, a Workflow and an Environment: the
 * `layered2json`, `workflow2json` and `environment2json` branches of
 * `linemodel_save.m`, and the inverses of `lqn_json_reader.h`,
 * `workflow_reader.h` and `environment_reader.h`. A Network is written by
 * `network_writer.h`, which the Environment writer reuses for every stage.
 *
 * WHICH SPELLING. The four codebases' writers disagree on a few optional keys,
 * and every reader accepts the union. Where they disagree this writer takes the
 * form every reader decodes to the same value:
 *
 *   host multiplicity/scheduling/quantum/speedFactor   always written (Python, JAR);
 *       MATLAB omits the defaults, and the JAR and native readers default them differently
 *   task multiplicity/scheduling, call `mean`           always written (Python, JAR)
 *   task think time        `thinkTime` plus `thinkTimeMean`/`thinkTimeSCV` (Python, JAR)
 *   activity `hostDemand`  always written, an Immediate as `{"type":"Immediate"}` (Python, JAR)
 *   ItemEntry `accessProb` a DiscreteSampler object over x = 1..n; the JAR reads only
 *       the object form, and Python's bare array is unreadable there
 *   entry `forwarding`, activity `thinkTime`            JAR keys, read by the JAR and C++ readers
 *   precedences            the TYPED form (`type` + `activities`), as MATLAB and Python write it
 *   Loop `loopCount`       written whenever the count is declared, as MATLAB does
 *   AND-join quorum        `preParams: [q]`, only when a quorum is set, as MATLAB does
 *   task `priority`, activity `callOrder`   written when not the default (0, STOCHASTIC), as the JAR does
 *
 * WHAT IS DROPPED. An entry's `type` (PH1PH2/NONE) and a host declared as `Host` rather than
 * `Processor` have no key in any codebase or in the schema; no solver reads either, so the
 * document reloads as the same model and they are left out rather than refused.
 *
 * WHAT IS REFUSED. An `LqnModel` field with no wire key in ANY codebase's
 * reader is refused by name rather than dropped, because a document that
 * reloads as a different model is worse than no document: a phase-2 activity,
 * a call group (RROBIN/JSQ), delayed-hit retrieval, the
 * load/class/joint dependence tables and server pools. An Environment stage
 * holding a LayeredNetwork is refused too, since its `LqnStruct` cannot be
 * turned back into the declared model and no reference writer serializes it.
 *
 * KEY ORDER is alphabetical, since `nlohmann::json` keeps objects sorted, where
 * the MATLAB and Python writers keep insertion order. JSON objects are unordered
 * and every reader looks keys up by name, so the order carries no meaning.
 */

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "line/io/network_writer.h"
#include "line/lang/lang_types.h"
#include "line/lang/lqn/lqn_reader.h"
#include "line/lang/qn/environment.h"
#include "line/lang/workflow/workflow.h"
#include "line/util/error.h"

namespace line {
namespace io {

namespace detail {

template <class T>
double to_d(const T& x) {
    return num_traits<T>::to_double(x);
}

/** A count as a JSON integer when it is one, as the reference writers emit it. */
inline json count_to_json(double v) {
    if (std::isfinite(v) && v == std::floor(v) && std::fabs(v) < 9.0e15)
        return json(static_cast<long long>(v));
    return json(v);
}

/** `inf_multiplicity()`: an infinite multiplicity is Java's Integer.MAX_VALUE on the wire. */
inline json lqn_mult_to_json(double m) {
    if (!std::isfinite(m)) return json(2147483647LL);
    return count_to_json(m);
}

/** A distribution that is really there: neither unset, immediate nor of zero mean. */
template <class T>
bool lqn_dist_declared(const lang::Distrib<T>& d) {
    return d.type != lang::ProcessType::DISABLED && d.type != lang::ProcessType::IMMEDIATE &&
           to_d(d.mean) > lang::GlobalConstants::FineTol;
}

/** `prectype_to_str`: the JAR-compatible precedence spelling of the Workflow wire. */
inline const char* prectype_to_str(lang::PrecedenceType t) {
    switch (t) {
        case lang::PrecedenceType::PRE_SEQ: return "pre";
        case lang::PrecedenceType::PRE_AND: return "pre-AND";
        case lang::PrecedenceType::PRE_OR: return "pre-OR";
        case lang::PrecedenceType::POST_SEQ: return "post";
        case lang::PrecedenceType::POST_AND: return "post-AND";
        case lang::PrecedenceType::POST_OR: return "post-OR";
        case lang::PrecedenceType::POST_LOOP: return "post-LOOP";
        case lang::PrecedenceType::POST_CACHE: return "post-CACHE";
        default: break;
    }
    throw UnsupportedError("linemodel_writer: a precedence carries a type with no wire spelling");
}

/**
 * `lincon2json`: admission-constraint rows in the NAMED form. A positional (A, b)
 * is resolved against `cols`, the element's operands in declaration order.
 */
template <class T>
json lincon_to_json(const std::vector<lqn::detail::RawLinConRow<T>>& named, const Matrix<T>* A,
                    const std::vector<T>* b, const std::vector<std::string>& cols,
                    const std::string& elem) {
    json rows = json::array();
    if (A != NULL && !A->empty()) {
        for (std::size_t r = 0; r < A->rows(); ++r) {
            json ops = json::array(), cf = json::array();
            for (std::size_t c = 0; c < A->cols(); ++c) {
                const double a = to_d((*A)(r, c));
                if (a == 0.0) continue;
                if (c >= cols.size())
                    throw InputError("linemodel_save: admission constraint on " + elem +
                                     " references column " + std::to_string(c + 1) +
                                     " but the element has only " + std::to_string(cols.size()) +
                                     " operands");
                ops.push_back(cols[c]);
                cf.push_back(a);
            }
            if (ops.empty()) continue;
            json row;
            row["operands"] = ops;
            row["coeffs"] = cf;
            row["cap"] = to_d(b->at(r));
            rows.push_back(row);
        }
    }
    for (const lqn::detail::RawLinConRow<T>& nr : named) {
        json row;
        row["operands"] = nr.names;
        json cf = json::array();
        for (const T& c : nr.coeffs) cf.push_back(to_d(c));
        row["coeffs"] = cf;
        row["cap"] = to_d(nr.cap);
        rows.push_back(row);
    }
    return rows;
}

/** Refuse, by name, a declared LqnModel feature no reader on the wire can rebuild. */
inline void lqn_refuse(const std::string& elem, const std::string& what) {
    throw UnsupportedError("linemodel_save: " + elem + " declares " + what +
                           ", which no model.json reader (MATLAB, JAR, Python or C++) rebuilds; "
                           "refusing rather than writing a document that reloads as a "
                           "different model");
}

template <class T>
void lqn_refuse_dependence(const std::string& elem, const std::vector<T>& lld,
                           const lang::CdScaling<T>& cd, const lang::CdScaling<T>& jd) {
    if (!lld.empty()) lqn_refuse(elem, "a load-dependent service scaling");
    if (cd) lqn_refuse(elem, "a class-dependent service scaling");
    if (jd) lqn_refuse(elem, "a joint-dependent service scaling");
}

/** One typed-form precedence, following `layered2json`. */
template <class T>
json lqn_precedence_to_json(const lqn::detail::RawPrecedence<T>& p, const std::string& task) {
    typedef lang::PrecedenceType P;
    json pj;
    pj["task"] = task;
    std::vector<std::string> both(p.preacts);
    both.insert(both.end(), p.postacts.begin(), p.postacts.end());
    if (p.pretype == P::PRE_SEQ && p.posttype == P::POST_SEQ) {
        pj["type"] = "Serial";
        pj["activities"] = both;
    } else if (p.pretype == P::PRE_SEQ && p.posttype == P::POST_AND) {
        pj["type"] = "AndFork";
        pj["activities"] = both;
    } else if (p.pretype == P::PRE_AND && p.posttype == P::POST_SEQ) {
        pj["type"] = "AndJoin";
        pj["activities"] = both;
        // the quorum travels as a one-element preParams, as MATLAB linemodel_save writes it
        if (p.has_quorum) pj["preParams"] = json::array({p.quorum});
    } else if (p.pretype == P::PRE_SEQ && p.posttype == P::POST_OR) {
        pj["type"] = "OrFork";
        pj["activities"] = both;
        if (!p.postparams.empty()) {
            json pr = json::array();
            for (const T& x : p.postparams) pr.push_back(to_d(x));
            pj["probabilities"] = pr;
        }
    } else if (p.pretype == P::PRE_OR && p.posttype == P::POST_SEQ) {
        if (!p.preparams.empty()) lqn_refuse("task '" + task + "'", "OR-join branch probabilities");
        pj["type"] = "OrJoin";
        pj["activities"] = both;
    } else if (p.posttype == P::POST_LOOP) {
        // postacts is the body followed by the exit; the trigger travels as preActivity
        pj["type"] = "Loop";
        pj["activities"] = p.postacts;
        if (!p.preacts.empty()) pj["preActivity"] = p.preacts[0];
        if (!p.postparams.empty()) {
            for (const T& c : p.postparams)
                if (to_d(c) != to_d(p.postparams[0]))
                    lqn_refuse("task '" + task + "'", "a loop whose body activities repeat "
                                                      "different numbers of times");
            pj["loopCount"] = to_d(p.postparams[0]);
        }
    } else if (p.pretype == P::PRE_SEQ && p.posttype == P::POST_CACHE) {
        pj["type"] = "CacheAccess";
        pj["activities"] = both;
    } else {
        lqn_refuse("task '" + task + "'", "a precedence of a pre/post type pair with no typed form");
    }
    return pj;
}

}  // namespace detail

/** `layered2json`: the `model` object of a LayeredNetwork document. */
template <class T>
detail::json lqn_model_to_json(const lqn::LqnModel<T>& m) {
    using detail::json;
    using detail::to_d;
    json out;
    out["type"] = "LayeredNetwork";
    out["name"] = m.name.empty() ? std::string("LQN") : m.name;

    // Operands of a host constraint are its tasks, of a task constraint its entries.
    std::vector<std::vector<std::string>> tasks_of(m.procs.size()), entries_of(m.tasks.size());
    for (const auto& t : m.tasks) tasks_of.at(t.proc_slot).push_back(t.name);
    for (const auto& e : m.entries) entries_of.at(e.task_slot).push_back(e.name);

    json hosts = json::array();
    for (std::size_t i = 0; i < m.procs.size(); ++i) {
        const lqn::detail::RawProc& p = m.procs[i];
        const std::string who = "host '" + p.name + "'";
        if (m.proc_lldscaling.count(i) || m.proc_cdscaling.count(i) || m.proc_jdscaling.count(i))
            detail::lqn_refuse(who, "a queue-dependent service scaling");
        if (m.proc_pools.count(i) && !m.proc_pools.at(i).empty())
            detail::lqn_refuse(who, "server pools");
        json h;
        h["name"] = p.name;
        h["multiplicity"] = detail::lqn_mult_to_json(p.mult);
        h["scheduling"] = detail::sched_to_json(p.sched);
        h["quantum"] = p.quantum;
        h["speedFactor"] = p.speed_factor;
        if (p.repl > 1) h["replication"] = detail::count_to_json(p.repl);
        static const std::vector<lqn::detail::RawLinConRow<T>> kNone;
        const auto named = m.proc_linconrows.find(i);
        const auto pos = m.proc_lincon.find(i);
        const json rows = detail::lincon_to_json<T>(
            named == m.proc_linconrows.end() ? kNone : named->second,
            pos == m.proc_lincon.end() ? NULL : &pos->second.first,
            pos == m.proc_lincon.end() ? NULL : &pos->second.second, tasks_of[i], p.name);
        if (!rows.empty()) h["admissionConstraints"] = rows;
        hosts.push_back(h);
    }
    out["hosts"] = hosts;

    json tasks = json::array();
    for (std::size_t i = 0; i < m.tasks.size(); ++i) {
        const lqn::detail::RawTask<T>& t = m.tasks[i];
        const std::string who = "task '" + t.name + "'";
        detail::lqn_refuse_dependence(who, t.lldscaling, t.cdscaling, t.jdscaling);
        if (!t.pools.empty()) detail::lqn_refuse(who, "server pools");
        if (t.retrieval) detail::lqn_refuse(who, "delayed-hit retrieval (setRetrieval)");
        json tj;
        tj["name"] = t.name;
        tj["host"] = m.procs.at(t.proc_slot).name;
        tj["multiplicity"] = detail::lqn_mult_to_json(t.mult);
        tj["scheduling"] = detail::sched_to_json(t.sched);
        if (t.repl > 1) tj["replication"] = detail::count_to_json(t.repl);
        if (t.priority != 0) tj["priority"] = t.priority;
        if (detail::lqn_dist_declared(t.thinktime)) {
            tj["thinkTime"] = detail::dist_to_json(t.thinktime);
            tj["thinkTimeMean"] = to_d(t.thinktime.mean);
            tj["thinkTimeSCV"] = to_d(t.thinktime.scv);
        }
        if (!t.fanin.empty()) {
            json fi = json::object();
            for (const auto& kv : t.fanin) fi[kv.first] = detail::count_to_json(kv.second);
            tj["fanIn"] = fi;
        }
        if (!t.fanout.empty()) {
            json fo = json::object();
            for (const auto& kv : t.fanout) fo[kv.first] = detail::count_to_json(kv.second);
            tj["fanOut"] = fo;
        }
        const json rows = detail::lincon_to_json<T>(t.linconrows, &t.lincon_A, &t.lincon_b,
                                                    entries_of[i], t.name);
        if (!rows.empty()) tj["admissionConstraints"] = rows;
        const bool has_setup = detail::lqn_dist_declared(t.setuptime);
        if (has_setup) tj["setupTime"] = detail::dist_to_json(t.setuptime);
        if (detail::lqn_dist_declared(t.delayofftime))
            tj["delayOffTime"] = detail::dist_to_json(t.delayofftime);
        if (t.nitems > 0) {
            tj["taskType"] = "CacheTask";
            tj["totalItems"] = static_cast<long long>(t.nitems);
            tj["cacheCapacity"] = t.itemcap;
            tj["replacementStrategy"] = detail::replacement_to_json(t.replacestrat);
        } else if (has_setup) {
            tj["taskType"] = "SetupTask";
        }
        tasks.push_back(tj);
    }
    out["tasks"] = tasks;

    json entries = json::array();
    std::map<std::string, std::string> reply_of;  // activity -> entry it replies to
    for (const lqn::detail::RawEntry<T>& e : m.entries) {
        json ej;
        ej["name"] = e.name;
        ej["task"] = m.tasks.at(e.task_slot).name;
        if (e.has_arrival) ej["arrival"] = detail::dist_to_json(e.arrival);
        if (!e.fwd_dest.empty()) {
            json fw = json::array();
            for (std::size_t k = 0; k < e.fwd_dest.size(); ++k) {
                json f;
                f["dest"] = e.fwd_dest[k];
                f["prob"] = k < e.fwd_prob.size() ? to_d(e.fwd_prob[k]) : 1.0;
                fw.push_back(f);
            }
            ej["forwarding"] = fw;
        }
        if (e.cardinality > 0) {
            ej["entryType"] = "ItemEntry";
            ej["totalItems"] = static_cast<long long>(e.cardinality);
            if (!e.popularity.empty()) {
                std::vector<T> x;
                for (std::size_t k = 1; k <= e.popularity.size(); ++k)
                    x.push_back(num_traits<T>::from_int(static_cast<long>(k)));
                ej["accessProb"] =
                    detail::dist_to_json(lang::Distrib<T>::discrete_sampler(e.popularity, x));
            }
        }
        for (const std::string& a : e.reply_activities) reply_of[a] = e.name;
        entries.push_back(ej);
    }
    out["entries"] = entries;

    json acts = json::array();
    for (const lqn::detail::RawActivity<T>& a : m.acts) {
        const std::string who = "activity '" + a.name + "'";
        if (a.phase != 1) detail::lqn_refuse(who, "phase " + std::to_string(a.phase));
        if (!a.call_groups.empty()) detail::lqn_refuse(who, "a routed call group (RROBIN/JSQ)");
        json aj;
        aj["name"] = a.name;
        aj["task"] = m.tasks.at(a.task_slot).name;
        aj["hostDemand"] = detail::dist_to_json(a.hostdem);
        if (!a.bound_to_entry.empty()) aj["boundToEntry"] = a.bound_to_entry;
        const auto r = reply_of.find(a.name);
        if (r != reply_of.end()) aj["repliesTo"] = r->second;
        if (detail::lqn_dist_declared(a.thinktime)) aj["thinkTime"] = detail::dist_to_json(a.thinktime);
        if (a.call_order != "STOCHASTIC") aj["callOrder"] = a.call_order;
        const char* kKey[2] = {"synchCalls", "asynchCalls"};
        const std::vector<lqn::detail::RawCall<T>>* calls[2] = {&a.sync_calls, &a.async_calls};
        for (int w = 0; w < 2; ++w) {
            if (calls[w]->empty()) continue;
            json arr = json::array();
            for (const lqn::detail::RawCall<T>& c : *calls[w]) {
                json cj;
                cj["dest"] = c.dest;
                cj["mean"] = to_d(c.mean);
                arr.push_back(cj);
            }
            aj[kKey[w]] = arr;
        }
        acts.push_back(aj);
    }
    out["activities"] = acts;

    json precs = json::array();
    for (const lqn::detail::RawTask<T>& t : m.tasks)
        for (const lqn::detail::RawPrecedence<T>& p : t.precedences)
            precs.push_back(detail::lqn_precedence_to_json(p, t.name));
    if (!precs.empty()) out["precedences"] = precs;
    return out;
}

/** `workflow2json`: the `model` object of a Workflow document. */
template <class T>
detail::json workflow_to_json(const workflow::Workflow<T>& wf) {
    using detail::json;
    json out;
    out["type"] = "Workflow";
    out["name"] = wf.name();
    json acts = json::array();
    for (const workflow::WorkflowActivity<T>& a : wf.activities()) {
        json aj;
        aj["name"] = a.name();
        aj["hostDemand"] = detail::dist_to_json(a.host_demand());
        acts.push_back(aj);
    }
    out["activities"] = acts;
    json precs = json::array();
    for (const workflow::Precedence<T>& p : wf.precedences()) {
        json pj;
        pj["preActs"] = p.pre_acts;
        pj["postActs"] = p.post_acts;
        pj["preType"] = detail::prectype_to_str(p.pre_type);
        pj["postType"] = detail::prectype_to_str(p.post_type);
        if (!p.pre_params.empty()) pj["preParams"] = detail::vec_to_json(p.pre_params);
        if (!p.post_params.empty()) pj["postParams"] = detail::vec_to_json(p.post_params);
        precs.push_back(pj);
    }
    out["precedences"] = precs;
    return out;
}

/**
 * `environment2json`: the `model` object of an Environment document.
 *
 * The stages go out EXPANDED, each with its own network, plus the `nodeFailures`
 * descriptors that alone carry the reset policies, as both reference writers do.
 * A reset FUNCTION on an ordinary arc has no wire form in any codebase; it is
 * dropped with a warning, as MATLAB warns for a custom node-failure policy.
 */
template <class T>
detail::json environment_to_json(const env::Environment<T>& e) {
    using detail::json;
    json out;
    out["type"] = "Environment";
    out["name"] = e.name();
    const std::size_t E = e.nstages();
    out["numStages"] = static_cast<long long>(E);

    json stages = json::array();
    for (std::size_t s = 0; s < E; ++s) {
        const env::EnvStage<T>& st = e.stage(s);
        if (st.has_lqn)
            throw UnsupportedError(
                "linemodel_save: stage '" + st.name +
                "' holds a LayeredNetwork, whose finalized LqnStruct cannot be written back as "
                "the declared model; no reference writer serializes a layered stage either");
        json sj;
        sj["name"] = st.name;
        if (!st.type.empty()) sj["type"] = st.type;
        if (st.has_model) sj["model"] = network_to_json(st.model);
        stages.push_back(sj);
    }
    out["stages"] = stages;

    // Arcs into or out of a DOWN_<node> stage are a node failure's, whose reset travels in nodeFailures.
    std::vector<bool> failure_stage(E, false);
    for (const env::NodeFailure<T>& nf : e.node_failures()) {
        const std::size_t d = e.find_stage(env::Environment<T>::down_stage_name(nf.node));
        if (d < E) failure_stage[d] = true;
    }
    json trans = json::array();
    for (std::size_t a = 0; a < E; ++a)
        for (std::size_t b = 0; b < E; ++b) {
            const env::EnvArc<T>& arc = e.arc(a, b);
            if (!arc.enabled || arc.dist.type == lang::ProcessType::DISABLED) continue;
            if ((arc.reset || arc.reset_rates) && !failure_stage[a] && !failure_stage[b])
                std::cerr << "[LINE] Warning: linemodel_save: the transition " << e.stage(a).name
                          << " -> " << e.stage(b).name
                          << " carries a reset function, which cannot be serialized to JSON; the "
                             "saved model reloads without it.\n";
            json tj;
            tj["from"] = static_cast<long long>(a);
            tj["to"] = static_cast<long long>(b);
            tj["distribution"] = detail::dist_to_json(arc.dist);
            trans.push_back(tj);
        }
    out["transitions"] = trans;

    json fails = json::array();
    for (const env::NodeFailure<T>& nf : e.node_failures()) {
        json nj;
        nj["node"] = nf.node;
        nj["breakdownRate"] = detail::dist_to_json(nf.breakdown);
        if (nf.has_repair) nj["repairRate"] = detail::dist_to_json(nf.repair);
        nj["downService"] = detail::dist_to_json(nf.down_service);
        if (nf.breakdown_reset == "custom")
            std::cerr << "[LINE] Warning: linemodel_save: node failure on \"" << nf.node
                      << "\" uses a custom breakdown reset function, which cannot be serialized "
                         "to JSON; the saved model falls back to the 'keep' policy on reload.\n";
        else
            nj["breakdownResetPolicy"] = nf.breakdown_reset;
        if (!nf.repair_reset.empty()) {
            if (nf.repair_reset == "custom")
                std::cerr << "[LINE] Warning: linemodel_save: node failure on \"" << nf.node
                          << "\" uses a custom repair reset function, which cannot be serialized "
                             "to JSON; the saved model falls back to the 'keep' policy on reload.\n";
            else
                nj["repairResetPolicy"] = nf.repair_reset;
        }
        fails.push_back(nj);
    }
    if (!fails.empty()) out["nodeFailures"] = fails;
    return out;
}

/** The `{format, version, model}` envelope around a `model` object, with the wire's non-finites. */
inline detail::json linemodel_envelope(const detail::json& model) {
    detail::json root;
    root["format"] = "line-model";
    root["version"] = "1.0";
    root["model"] = model;
    detail::wire_nonfinite(root);
    return root;
}

/** Write an envelope, indented as `write_network_json` indents it. */
inline void write_linemodel_json(const detail::json& root, const std::string& path) {
    std::ofstream out(path.c_str());
    if (!out) throw InputError("linemodel_save: cannot open " + path + " for writing");
    out << root.dump(2) << "\n";
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_LINEMODEL_WRITER_H
