/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_LQN_JSON_READER_H
#define LINE_IO_LQN_JSON_READER_H

/**
 * `model.json` with `type: "LayeredNetwork"` -> `LqnStruct`, via `LqnBuilder`.
 *
 * WHY THIS EXISTS BESIDE `lqn_reader.h`. The .lqnx reader covers the LQNS file
 * format, which is what an external tool writes; this covers the LINE
 * interchange format, which is what `linemodel_save.m`, `save_model` and the
 * JAR `LineModelIO` write for a `LayeredNetwork`. The two are not
 * interchangeable in either direction: .lqnx cannot carry a non-reference
 * task's think time, a CacheTask, an ItemEntry or a SetupTask, and model.json
 * carries all four. Without this reader a layered model could reach the C++
 * port only by being re-exported through XML, losing exactly those constructs.
 *
 * ORDERING IS PART OF THE FORMAT. `LqnBuilder` assigns indices in call order,
 * and every later array -- the layer decomposition, the visit matrix, the
 * result tables -- is keyed on them. The document's own order is used
 * throughout (hosts, then tasks, then entries, then activities), which is the
 * order the writers emit and therefore the order MATLAB's own `getStruct`
 * produced before serialization.
 *
 * WHAT IS REFUSED. Fan-in and fan-out replication (`fanIn`/`fanOut`) and the
 * per-element admission constraints have no `LqnBuilder` counterpart that the
 * layer decomposition honours, and a call multiplicity silently applied at 1
 * would understate every visit through the replicated task. They are named
 * rather than dropped, as everywhere else on this boundary.
 */

#include <map>
#include <string>
#include <vector>

#include "line/io/network_reader.h"
#include <cctype>
#include <fstream>

#include "line/lang/lqn/lqn_builder.h"
#include "line/lang/lqn/lqn_reader.h"
#include "line/lang/lqn/lqn_struct.h"
#include "line/util/error.h"

namespace line {
namespace io {

namespace detail {

/** The `scheduling` name of a host or task, upper case as the writers emit it. */
inline lang::SchedStrategy lqn_sched_from_json(const std::string& s) {
    typedef lang::SchedStrategy S;
    if (s == "INF" || s == "inf") return S::INF;
    if (s == "FCFS" || s == "fcfs") return S::FCFS;
    if (s == "PS" || s == "ps") return S::PS;
    if (s == "HOL" || s == "hol") return S::HOL;
    if (s == "REF" || s == "ref") return S::REF;
    if (s == "SIRO" || s == "siro") return S::SIRO;
    if (s == "LCFS" || s == "lcfs") return S::LCFS;
    if (s == "LCFSPR" || s == "lcfspr") return S::LCFSPR;
    throw UnsupportedError("lqn_json_reader: unsupported scheduling discipline '" + s + "'");
}

/**
 * The wire's infinite multiplicity.
 *
 * The writers emit `Integer.MAX_VALUE` as a literal rather than a float
 * infinity, because an infinity would not survive an integer-valued export;
 * both spellings are accepted here and both mean the same thing.
 */
inline double lqn_mult_from_json(const json& v) {
    if (v.is_string()) {
        const std::string s = v.get<std::string>();
        if (s == "Infinity" || s == "inf") return std::numeric_limits<double>::infinity();
        return std::atof(s.c_str());
    }
    const double m = v.get<double>();
    return m >= 2147483647.0 ? std::numeric_limits<double>::infinity() : m;
}

}  // namespace detail

/**
 * Build an `LqnStruct<T>` from a parsed model.json envelope carrying a
 * LayeredNetwork.
 */
template <class T>
lqn::LqnStruct<T> build_lqn_from_json(const detail::json& root) {
    using detail::json;
    const json& model = root.contains("model") ? root.at("model") : root;
    const std::string mtype = model.value("type", std::string("LayeredNetwork"));
    if (mtype != "LayeredNetwork")
        throw UnsupportedError("lqn_json_reader: model type '" + mtype +
                               "' is not a LayeredNetwork; a Network model is read by "
                               "network_reader.h");

    lqn::LqnBuilder<T> b;

    // -- hosts (processors) --------------------------------------------------
    // `hosts` is the LINE spelling and `processors` the LQNS one; the schema
    // declares both and a writer may use either.
    const json empty = json::array();
    const json& hosts = model.contains("hosts")
                            ? model.at("hosts")
                            : (model.contains("processors") ? model.at("processors") : empty);
    for (const json& h : hosts) {
        const std::string name = h.at("name").get<std::string>();
        if (h.contains("admissionConstraints"))
            throw UnsupportedError(
                "lqn_json_reader: host '" + name +
                "' carries admission constraints, which this reader does not rebuild");
        const lang::SchedStrategy sched =
            detail::lqn_sched_from_json(h.value("scheduling", std::string("PS")));
        const double mult = h.contains("multiplicity")
                                ? detail::lqn_mult_from_json(h.at("multiplicity"))
                                : 1.0;
        b.processor(name, mult, sched, h.value("replication", 1.0));
    }

    // -- tasks ---------------------------------------------------------------
    for (const json& t : model.at("tasks")) {
        const std::string name = t.at("name").get<std::string>();
        if (t.contains("fanIn") || t.contains("fanOut"))
            throw UnsupportedError(
                "lqn_json_reader: task '" + name +
                "' declares fan-in or fan-out replication, which the layer decomposition in this "
                "port does not honour; a call multiplicity applied at 1 would understate every "
                "visit through the replicated task");
        if (t.contains("admissionConstraints"))
            throw UnsupportedError(
                "lqn_json_reader: task '" + name +
                "' carries admission constraints, which this reader does not rebuild");
        const std::string on = t.at("host").get<std::string>();
        const lang::SchedStrategy sched =
            detail::lqn_sched_from_json(t.value("scheduling", std::string("FCFS")));
        const double mult = t.contains("multiplicity")
                                ? detail::lqn_mult_from_json(t.at("multiplicity"))
                                : 1.0;
        const double repl = t.value("replication", 1.0);
        const std::string kind = t.value("taskType", std::string());
        if (kind == "CacheTask") {
            std::vector<int> cap;
            if (t.at("cacheCapacity").is_array())
                cap = t.at("cacheCapacity").get<std::vector<int> >();
            else
                cap.push_back(t.at("cacheCapacity").get<int>());
            b.cache_task(name, mult, sched, on, t.at("totalItems").get<std::size_t>(), cap,
                         detail::replacement_from_json(
                             t.value("replacementStrategy", std::string("FIFO"))),
                         repl);
        } else {
            b.task(name, mult, sched, on, repl);
        }
        // A think time is accepted on ANY task here, as the reference API
        // allows and as .lqnx cannot express.
        if (t.contains("thinkTime"))
            b.think_time(name, detail::dist_from_json<T>(t.at("thinkTime")));
        else if (t.contains("thinkTimeMean") && t.at("thinkTimeMean").get<double>() > 0)
            b.think_time(name, lang::Distrib<T>::exp_mean(num_traits<T>::from_double(
                                   t.at("thinkTimeMean").get<double>())));
        // BOTH SPELLINGS ARE ACCEPTED, the object form PREFERRED. MATLAB,
        // Python and this port all write `setupTime`/`delayOffTime` as full
        // distributions, so that is the canonical wire form; the JAR writes
        // `setupTimeMean`+`setupTimeSCV` instead. Documents already written in
        // the field carry only the scalars, and without this fallback both
        // fields vanish on read and a SetupTask arrives with no setup time at
        // all -- silently, which is the failure mode `.lqnx` already has. An
        // SCV cannot be recovered from a mean, so the scalar path reconstructs
        // an exponential and claims no more, exactly as `thinkTimeMean` does.
        const bool has_setup_obj = t.contains("setupTime");
        const bool has_setup_mean =
            t.contains("setupTimeMean") && t.at("setupTimeMean").get<double>() > 0;
        if (has_setup_obj || has_setup_mean) {
            const bool has_off_obj = t.contains("delayOffTime");
            const bool has_off_mean =
                t.contains("delayOffTimeMean") && t.at("delayOffTimeMean").get<double>() > 0;
            if (!has_off_obj && !has_off_mean)
                throw InputError("lqn_json_reader: task '" + name +
                                 "' declares a setup time with no delay-off time; a server that "
                                 "never shuts down pays the setup at most once");
            const lang::Distrib<T> su =
                has_setup_obj ? detail::dist_from_json<T>(t.at("setupTime"))
                              : lang::Distrib<T>::exp_mean(num_traits<T>::from_double(
                                    t.at("setupTimeMean").get<double>()));
            const lang::Distrib<T> doff =
                has_off_obj ? detail::dist_from_json<T>(t.at("delayOffTime"))
                            : lang::Distrib<T>::exp_mean(num_traits<T>::from_double(
                                  t.at("delayOffTimeMean").get<double>()));
            b.setup_time(name, su, doff);
        }
    }

    // -- entries -------------------------------------------------------------
    for (const json& e : model.at("entries")) {
        const std::string name = e.at("name").get<std::string>();
        const std::string on = e.at("task").get<std::string>();
        if (e.value("entryType", std::string()) == "ItemEntry") {
            const std::size_t card = e.at("totalItems").get<std::size_t>();
            // The popularity crosses as a DISTRIBUTION (a Zipf, typically) and
            // the builder wants the pmf over the items, which is the same
            // object the cache node's `pread` is.
            std::vector<T> pop;
            if (e.contains("accessProb"))
                pop = detail::pmf_from_json<T>(e.at("accessProb"), card);
            if (pop.empty())
                pop.assign(card, num_traits<T>::from_double(1.0 / double(card)));
            b.item_entry(name, on, card, pop);
        } else {
            b.entry(name, on);
        }
        if (e.contains("arrival"))
            b.open_arrival(name, detail::dist_from_json<T>(e.at("arrival")));
    }

    // -- activities ----------------------------------------------------------
    for (const json& a : model.at("activities")) {
        const std::string name = a.at("name").get<std::string>();
        const std::string on = a.at("task").get<std::string>();
        // An activity with no host demand is an IMMEDIATE one; the writers omit
        // the key in exactly that case, so the absence is the value.
        const lang::Distrib<T> dem = a.contains("hostDemand")
                                         ? detail::dist_from_json<T>(a.at("hostDemand"))
                                         : lang::Distrib<T>::immediate();
        b.activity(name, dem, on);
        if (a.contains("boundToEntry"))
            b.bound_to(name, a.at("boundToEntry").get<std::string>());
        else if (a.contains("boundTo"))
            b.bound_to(name, a.at("boundTo").get<std::string>());
        if (a.contains("repliesTo")) b.replies_to(name, a.at("repliesTo").get<std::string>());
        const char* kCallKey[2] = {"synchCalls", "asynchCalls"};
        for (int which = 0; which < 2; ++which) {
            if (!a.contains(kCallKey[which])) continue;
            for (const json& c : a.at(kCallKey[which])) {
                // `dest` is the LINE spelling and `entry` the schema's.
                const std::string dest = c.contains("dest")
                                             ? c.at("dest").get<std::string>()
                                             : c.at("entry").get<std::string>();
                const T mean = num_traits<T>::from_double(c.value("mean", 1.0));
                if (which == 0) b.sync_call(name, dest, mean);
                else b.async_call(name, dest, mean);
            }
        }
    }

    // Entry forwarding is declared on the ENTRY and resolved after every entry
    // exists, since it names one.
    for (const json& e : model.at("entries")) {
        if (!e.contains("forwarding")) continue;
        const std::string name = e.at("name").get<std::string>();
        for (const json& f : e.at("forwarding"))
            b.forward(name, f.at("dest").get<std::string>(),
                      num_traits<T>::from_double(f.value("prob", 1.0)));
    }

    // -- precedences ---------------------------------------------------------
    //
    // The TYPED form the writers emit: one `type` naming the pre/post pair,
    // plus the activity list with the pre-activities first. The EXPLICIT form
    // (`preType`/`postType` with separate lists) is the schema's other spelling
    // and is mapped onto the same builder calls.
    if (model.contains("precedences")) {
        for (const json& p : model.at("precedences")) {
            const std::string kind = p.value("type", std::string());
            if (kind.empty()) {
                // Explicit form.
                const std::string pre = p.at("preType").get<std::string>();
                const std::string post = p.at("postType").get<std::string>();
                std::vector<std::string> pres, posts;
                for (const json& x : p.at("preActs")) pres.push_back(x.get<std::string>());
                for (const json& x : p.at("postActs")) posts.push_back(x.get<std::string>());
                if (pre == "PRE_SEQ" && post == "POST_SEQ") b.serial(pres.at(0), posts.at(0));
                else if (pre == "PRE_SEQ" && post == "POST_AND") b.and_fork(pres.at(0), posts);
                else if (pre == "PRE_AND" && post == "POST_SEQ") b.and_join(pres, posts.at(0));
                else if (pre == "PRE_OR" && post == "POST_SEQ") b.or_join(pres, posts.at(0));
                else if (pre == "PRE_SEQ" && post == "POST_OR") {
                    std::vector<T> probs;
                    if (p.contains("postParams"))
                        probs = detail::num_vec_from_json<T>(p.at("postParams"));
                    b.or_fork(pres.at(0), posts, probs);
                } else if (post == "POST_LOOP") {
                    // `postActs` is the loop BODY followed by the activity the
                    // loop exits to, which is how `ActivityPrecedence.Loop`
                    // stores it and how the builder takes it back apart.
                    if (posts.size() < 2)
                        throw InputError(
                            "lqn_json_reader: a loop's postActs is the body followed by the "
                            "activity it exits to, so it holds at least two entries");
                    const T count = p.contains("postParams")
                                        ? detail::num_vec_from_json<T>(p.at("postParams")).at(0)
                                        : num_traits<T>::from_int(1);
                    b.loop(pres.at(0), std::vector<std::string>(posts.begin(), posts.end() - 1),
                           posts.back(), count);
                } else if (post == "POST_CACHE") {
                    b.cache_access(pres.at(0), posts.at(0), posts.at(1));
                } else {
                    throw UnsupportedError("lqn_json_reader: precedence " + pre + " -> " + post +
                                           " has no builder counterpart");
                }
                continue;
            }
            std::vector<std::string> acts;
            for (const json& x : p.at("activities")) acts.push_back(x.get<std::string>());
            if (kind == "Serial") {
                // A serial chain is written as one list, and the builder takes
                // it a pair at a time.
                for (std::size_t i = 0; i + 1 < acts.size(); ++i) b.serial(acts[i], acts[i + 1]);
            } else if (kind == "AndFork") {
                b.and_fork(acts.at(0), std::vector<std::string>(acts.begin() + 1, acts.end()));
            } else if (kind == "AndJoin") {
                b.and_join(std::vector<std::string>(acts.begin(), acts.end() - 1), acts.back());
            } else if (kind == "OrFork") {
                std::vector<T> probs;
                if (p.contains("probabilities"))
                    probs = detail::num_vec_from_json<T>(p.at("probabilities"));
                b.or_fork(acts.at(0), std::vector<std::string>(acts.begin() + 1, acts.end()),
                          probs);
            } else if (kind == "OrJoin") {
                b.or_join(std::vector<std::string>(acts.begin(), acts.end() - 1), acts.back());
            } else if (kind == "Loop") {
                // The trigger is written separately from the body, so `acts` is
                // the body alone.
                if (!p.contains("preActivity"))
                    throw InputError(
                        "lqn_json_reader: a Loop precedence carries its trigger as 'preActivity'");
                if (acts.size() < 2)
                    throw InputError(
                        "lqn_json_reader: a Loop precedence lists the body followed by the "
                        "activity it exits to, so it holds at least two entries");
                b.loop(p.at("preActivity").get<std::string>(),
                       std::vector<std::string>(acts.begin(), acts.end() - 1), acts.back(),
                       num_traits<T>::from_double(p.value("loopCount", 1.0)));
            } else if (kind == "CacheAccess") {
                if (acts.size() < 3)
                    throw InputError(
                        "lqn_json_reader: a CacheAccess precedence names the read activity plus "
                        "its hit and miss branches");
                b.cache_access(acts.at(0), acts.at(1), acts.at(2));
            } else {
                throw UnsupportedError("lqn_json_reader: precedence type '" + kind +
                                       "' has no builder counterpart");
            }
        }
    }
    return b.build();
}

/**
 * True when a file is a LayeredNetwork model.json rather than an .lqnx.
 *
 * Decided on the CONTENT, not the extension: `linemodel_save` writes `.json`
 * for both model kinds and the two readers cannot be told apart by name. The
 * first non-space character settles it -- `{` is JSON, `<` is XML -- and the
 * `type` field settles which JSON.
 */
inline bool is_layered_json(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in) return false;
    char c = 0;
    while (in.get(c) && std::isspace(static_cast<unsigned char>(c))) {}
    if (c != '{') return false;
    in.seekg(0);
    detail::json root;
    try {
        in >> root;
    } catch (const detail::json::parse_error&) {
        return false;
    }
    const detail::json& model = root.contains("model") ? root.at("model") : root;
    return model.value("type", std::string()) == "LayeredNetwork";
}

/** Parse a LayeredNetwork model.json file into an `LqnStruct<T>`. */
template <class T>
lqn::LqnStruct<T> read_lqn_json(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in) throw InputError("lqn_json_reader: cannot open " + path);
    detail::json root;
    try {
        in >> root;
    } catch (const detail::json::parse_error& e) {
        throw InputError("lqn_json_reader: malformed JSON in " + path + ": " + e.what());
    }
    return build_lqn_from_json<T>(root);
}

/**
 * Read a layered model from either interchange: the LINE `model.json` or the
 * LQNS `.lqnx`.
 *
 * Every layered entry point goes through here, so a caller never has to know
 * which of the two it was handed -- and, more to the point, a `model.json`
 * written by `linemodel_save` for a LayeredNetwork stops being unreadable
 * merely because the XML reader was the only one wired up.
 */
template <class T>
lqn::LqnStruct<T> read_layered_model(const std::string& path) {
    return is_layered_json(path) ? read_lqn_json<T>(path) : lqn::read_lqnx<T>(path);
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_LQN_JSON_READER_H
