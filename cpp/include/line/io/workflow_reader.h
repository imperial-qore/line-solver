/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_WORKFLOW_READER_H
#define LINE_IO_WORKFLOW_READER_H

/**
 * Reader for the LINE `model.json` interchange (a Workflow model) into a
 * `workflow::Workflow<T>` built through the programmatic builder.
 *
 * The wire format is the one `linemodel_save.m::workflow2json` and the Python
 * `_workflow_to_json` emit: a `{format, version, model}` envelope whose
 * `model.type == "Workflow"` carries `activities` (name plus `hostDemand`) and
 * `precedences` (`preActs`, `postActs`, `preType`, `postType` and the optional
 * `preParams` / `postParams`). Activities and precedences are fed to the SAME
 * `add_activity` / `add_precedence` the programmatic API uses, so a workflow
 * that reaches C++ this way is indistinguishable from one authored in code.
 *
 * WHY THIS EXISTS SEPARATELY from `network_reader.h`: a Workflow is not a
 * queueing network and shares none of its node, class or routing structure. It
 * reduces to ONE phase-type law, so it has no stations to solve. The two
 * readers share only `detail::dist_from_json`, which is why this header
 * includes `network_reader.h` rather than duplicating the distribution table.
 *
 * WHAT IT REFUSES, and by name rather than by silent degradation: a precedence
 * type the wire spells but this port has no composition rule for, an activity
 * whose `hostDemand` names a family `dist_from_json` cannot reconstruct, and
 * any unknown key at model / activity / precedence level. The last of these is
 * the same rule `network_reader.h` states at length: a key carrying model
 * semantics that no branch consumes would otherwise be dropped in silence, and
 * `to_ph()` would then return a confident law for a different workflow.
 *
 * The SEMANTIC checks (a loop count that is not a positive scalar, a quorum
 * join, a graph that is not series-parallel) are NOT repeated here. They belong
 * to `Workflow::validate()` and `to_ph()`, which every construction path runs,
 * so duplicating them would let the two drift apart.
 */

#include <string>
#include <vector>

#include "json.hpp"
#include "line/io/network_reader.h"
#include "line/lang/lang_types.h"
#include "line/lang/workflow/workflow.h"
#include "line/util/error.h"

namespace line {
namespace io {

namespace detail {

/**
 * Map the wire's precedence spelling to `lang::PrecedenceType`.
 *
 * The strings are the JAR-compatible ones `linemodel_save.m::prectype_to_str`
 * writes; the Python writer emits the same set. `POST_CACHE` is spelled by the
 * wire and accepted here, because refusing it at the READER would misreport a
 * layered cache-queueing workflow as unrepresentable when the truth is that
 * `to_ph` has no rule for it -- which `Workflow::validate()` says in its own
 * words, at the point where it is actually true.
 */
inline lang::PrecedenceType precedence_type_from_str(const std::string& s) {
    if (s == "pre") return lang::PrecedenceType::PRE_SEQ;
    if (s == "pre-AND") return lang::PrecedenceType::PRE_AND;
    if (s == "pre-OR") return lang::PrecedenceType::PRE_OR;
    if (s == "post") return lang::PrecedenceType::POST_SEQ;
    if (s == "post-AND") return lang::PrecedenceType::POST_AND;
    if (s == "post-OR") return lang::PrecedenceType::POST_OR;
    if (s == "post-LOOP") return lang::PrecedenceType::POST_LOOP;
    if (s == "post-CACHE") return lang::PrecedenceType::POST_CACHE;
    throw UnsupportedError("workflow_reader: unsupported precedence type '" + s +
                           "'; the wire spells one of pre, pre-AND, pre-OR, post, post-AND, "
                           "post-OR, post-LOOP, post-CACHE");
}

/** Read a `preActs` / `postActs` array as a list of activity names. */
inline std::vector<std::string> activity_names_from_json(const json& arr, const char* key) {
    if (!arr.is_array())
        throw InputError(std::string("workflow_reader: '") + key + "' must be an array of "
                         "activity names");
    std::vector<std::string> out;
    out.reserve(arr.size());
    for (const json& v : arr) {
        if (!v.is_string())
            throw InputError(std::string("workflow_reader: '") + key +
                             "' must hold activity NAMES, not objects; the wire references an "
                             "activity by the name its 'activities' entry declares");
        out.push_back(v.get<std::string>());
    }
    return out;
}

/** Read an optional `preParams` / `postParams` array of scalars. */
template <class T>
std::vector<T> params_from_json(const json& obj, const char* key) {
    std::vector<T> out;
    if (!obj.contains(key)) return out;
    const json& arr = obj.at(key);
    // A scalar is accepted as the one-element form: `Loop` writes its count as
    // a bare number in some writers and as a singleton array in others.
    if (arr.is_number()) {
        out.push_back(num_traits<T>::from_double(num_from_json(arr)));
        return out;
    }
    if (!arr.is_array())
        throw InputError(std::string("workflow_reader: '") + key +
                         "' must be a number or an array of numbers");
    for (const json& v : arr) out.push_back(num_traits<T>::from_double(num_from_json(v)));
    return out;
}

/**
 * Refuse any key this reader does not consume, at every level.
 *
 * Same rule as `network_reader.h::reject_unknown_keys`, and for the same
 * reason: a dropped key is indistinguishable from a key that was never there.
 */
inline void reject_unknown_workflow_keys(const json& model) {
    static const char* kModelKeys[] = {"type", "name", "activities", "precedences"};
    static const char* kActKeys[] = {"name", "hostDemand"};
    static const char* kPrecKeys[] = {"preActs", "postActs", "preType", "postType",
                                      "preParams", "postParams"};
    auto known = [](const char* const* tab, std::size_t n, const std::string& k) {
        for (std::size_t i = 0; i < n; ++i)
            if (k == tab[i]) return true;
        return false;
    };
    const std::string why =
        "', which this reader does not implement. Refusing rather than dropping it: a "
        "constraint silently discarded here would make to_ph() return a confident law "
        "for a different workflow";
    for (auto it = model.begin(); it != model.end(); ++it)
        if (!known(kModelKeys, sizeof(kModelKeys) / sizeof(*kModelKeys), it.key()))
            throw UnsupportedError("workflow_reader: the model carries '" + it.key() + why);
    if (model.contains("activities"))
        for (const json& act : model.at("activities"))
            for (auto it = act.begin(); it != act.end(); ++it)
                if (!known(kActKeys, sizeof(kActKeys) / sizeof(*kActKeys), it.key()))
                    throw UnsupportedError("workflow_reader: activity '" +
                                           act.value("name", std::string("?")) + "' carries '" +
                                           it.key() + why);
    if (model.contains("precedences"))
        for (const json& prec : model.at("precedences"))
            for (auto it = prec.begin(); it != prec.end(); ++it)
                if (!known(kPrecKeys, sizeof(kPrecKeys) / sizeof(*kPrecKeys), it.key()))
                    throw UnsupportedError("workflow_reader: a precedence carries '" + it.key() +
                                           why);
}

}  // namespace detail

/**
 * Build a `workflow::Workflow<T>` from a parsed model.json envelope.
 *
 * One pass: an activity is declared before any precedence may reference it,
 * which the wire guarantees by carrying `activities` as its own block. The
 * result is NOT validated here -- `validate()` and `to_ph()` do that on the
 * composed model, where the series-parallel and quorum rules actually apply.
 */
template <class T>
workflow::Workflow<T> build_workflow_from_json(const detail::json& root) {
    using detail::json;
    const json& model = root.contains("model") ? root.at("model") : root;
    const std::string mtype = model.value("type", std::string("Workflow"));
    if (mtype != "Workflow")
        throw UnsupportedError("workflow_reader: model type '" + mtype +
                               "' is not a Workflow; a Network is read by network_reader.h and an "
                               "Environment by environment_reader.h");

    detail::reject_unknown_workflow_keys(model);

    workflow::Workflow<T> wf(model.value("name", std::string("workflow")));

    if (model.contains("activities")) {
        for (const json& act : model.at("activities")) {
            if (!act.contains("name"))
                throw InputError("workflow_reader: every activity needs a 'name'");
            const std::string name = act.at("name").get<std::string>();
            if (act.contains("hostDemand"))
                wf.add_activity(name, detail::dist_from_json<T>(act.at("hostDemand")));
            else
                // The writers ALWAYS emit a hostDemand (a float mean is written
                // as the Exp of that mean), so an absent one is a hand-authored
                // document. Unit mean matches what the Python reader does.
                wf.add_activity(name, lang::Distrib<T>::exp_mean(num_traits<T>::from_double(1.0)));
        }
    }

    if (model.contains("precedences")) {
        for (const json& prec : model.at("precedences")) {
            if (!prec.contains("preActs") || !prec.contains("postActs"))
                throw InputError("workflow_reader: every precedence needs 'preActs' and "
                                 "'postActs'");
            workflow::Precedence<T> p;
            p.pre_acts = detail::activity_names_from_json(prec.at("preActs"), "preActs");
            p.post_acts = detail::activity_names_from_json(prec.at("postActs"), "postActs");
            p.pre_type =
                detail::precedence_type_from_str(prec.value("preType", std::string("pre")));
            p.post_type =
                detail::precedence_type_from_str(prec.value("postType", std::string("post")));
            p.pre_params = detail::params_from_json<T>(prec, "preParams");
            p.post_params = detail::params_from_json<T>(prec, "postParams");
            wf.add_precedence(p);
        }
    }

    return wf;
}

/** Read a Workflow model.json off disk. */
template <class T>
workflow::Workflow<T> read_workflow_json(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in)
        throw InputError("workflow_reader: cannot open '" + path + "'");
    detail::json root;
    in >> root;
    return build_workflow_from_json<T>(root);
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_WORKFLOW_READER_H
