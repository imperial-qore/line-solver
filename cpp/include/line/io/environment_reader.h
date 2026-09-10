/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_ENVIRONMENT_READER_H
#define LINE_IO_ENVIRONMENT_READER_H

/**
 * Reader for the LINE `model.json` interchange of an ENVIRONMENT model into an
 * `env::Environment<T>`.
 *
 * The wire format is the one `environment2json` in `linemodel_save.m` and
 * `_environment_to_json` in `python/line_solver/io/linemodel_io.py` emit:
 * `{type: "Environment", name, numStages, stages: [{name, type, model}],
 *   transitions: [{from, to, distribution}]}`, with `from`/`to` ZERO-BASED on
 * the wire in both writers (MATLAB converts from its 1-based stage index when
 * writing and back when reading, so a C++ reader that treated them as 1-based
 * would silently rotate the environment process).
 *
 * Each stage's `model` is an ordinary Network envelope, so it is handed to
 * `build_network_from_json` unchanged: a stage network that MVA or the fluid
 * solver can read on its own is exactly the stage network ENV runs, and no
 * second, weaker parser exists for it.
 *
 * `nodeFailures` HAS TWO ROLES, and which one applies is decided exactly as
 * both readers decide it: by whether a `DOWN_<node>` stage is already declared.
 *
 *   MACRO FORM (no DOWN stage declared): the block IS the environment. The one
 *   declared stage holds the base model, and each entry expands into a
 *   `DOWN_<node>` stage carrying the degraded service plus the breakdown and
 *   repair arcs. `transitions` must be absent, because the block implies them.
 *
 *   EXPANDED FORM (the writers' own output): the stages and arcs came off the
 *   wire and the block adds only what the wire cannot express, the QUEUE-LENGTH
 *   RESET POLICIES `breakdownResetPolicy` and `repairResetPolicy`. They decide
 *   what the next stage starts from at a switch, so dropping them would run the
 *   identity reset and report a confident answer for a different model -- which
 *   is why this reader refused the whole block before the policies were ported.
 *
 * A `custom` policy is refused by `env_reset_policy`, and correctly: both
 * writers warn and OMIT the key when the policy is a function handle, so the
 * only way `custom` reaches here is a hand-written file claiming a policy that
 * no file can carry.
 */

#include <fstream>
#include <string>
#include <vector>

#include "json.hpp"
#include "line/io/network_reader.h"
#include "line/lang/qn/environment.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace io {

namespace detail {

/** The keys an Environment envelope may carry; anything else is refused. */
inline void reject_unconsumed_env_keys(const json& model) {
    static const char* kEnvKeys[] = {"name",        "type",   "numStages", "stages",
                                     "transitions", "format", "version",   "nodeFailures"};
    static const char* kStageKeys[] = {"name", "type", "model"};
    static const char* kArcKeys[] = {"from", "to", "distribution"};
    static const char* kFailKeys[] = {"node",        "breakdownRate",        "repairRate",
                                      "downService", "breakdownResetPolicy", "repairResetPolicy"};
    auto known = [](const char* const* tab, std::size_t n, const std::string& k) {
        for (std::size_t i = 0; i < n; ++i)
            if (k == tab[i]) return true;
        return false;
    };
    const std::string why =
        "', which this reader does not implement. Refusing rather than dropping it: a "
        "constraint silently discarded here would make the solver return a confident answer "
        "for a different environment";
    for (json::const_iterator it = model.begin(); it != model.end(); ++it)
        if (!known(kEnvKeys, sizeof(kEnvKeys) / sizeof(*kEnvKeys), it.key()))
            throw UnsupportedError("environment_reader: the model carries '" + it.key() + why);
    if (model.contains("nodeFailures"))
        for (const json& nf : model.at("nodeFailures"))
            for (json::const_iterator it = nf.begin(); it != nf.end(); ++it)
                if (!known(kFailKeys, sizeof(kFailKeys) / sizeof(*kFailKeys), it.key()))
                    throw UnsupportedError("environment_reader: the node failure on '" +
                                           nf.value("node", std::string("?")) + "' carries '" +
                                           it.key() + why);
    if (model.contains("stages"))
        for (const json& st : model.at("stages"))
            for (json::const_iterator it = st.begin(); it != st.end(); ++it)
                if (!known(kStageKeys, sizeof(kStageKeys) / sizeof(*kStageKeys), it.key()))
                    throw UnsupportedError("environment_reader: stage '" +
                                           st.value("name", std::string("?")) + "' carries '" +
                                           it.key() + why);
    if (model.contains("transitions"))
        for (const json& tr : model.at("transitions"))
            for (json::const_iterator it = tr.begin(); it != tr.end(); ++it)
                if (!known(kArcKeys, sizeof(kArcKeys) / sizeof(*kArcKeys), it.key()))
                    throw UnsupportedError("environment_reader: an environment transition carries '" +
                                           it.key() + why);
}

/** One decoded `nodeFailures` entry: `nodefailure_fields` in both readers. */
template <class T>
struct NodeFailureSpec {
    std::string node;
    lang::Distrib<T> breakdown, down_service, repair;
    bool has_repair = false;
    std::string breakdown_reset = "keep";
    std::string repair_reset = "keep";
};

/**
 * Decode one entry.
 *
 * `breakdownRate` and `repairRate` carry FULL DISTRIBUTIONS and not scalar
 * rates, despite the names; both writers note it and both readers do the same.
 */
template <class T>
NodeFailureSpec<T> node_failure_fields(const json& nf) {
    NodeFailureSpec<T> s;
    if (!nf.contains("node"))
        throw InputError(
            "environment_reader: a 'nodeFailures' entry is missing the required 'node' field, "
            "which names the node that breaks down");
    s.node = nf.at("node").get<std::string>();
    if (!nf.contains("breakdownRate") || nf.at("breakdownRate").is_null())
        throw InputError("environment_reader: the node failure on '" + s.node +
                         "' is missing the required 'breakdownRate' field, the time to failure");
    if (!nf.contains("downService") || nf.at("downService").is_null())
        throw InputError("environment_reader: the node failure on '" + s.node +
                         "' is missing the required 'downService' field, the service the node "
                         "gives while it is down");
    s.breakdown = dist_from_json<T>(nf.at("breakdownRate"));
    s.down_service = dist_from_json<T>(nf.at("downService"));
    if (nf.contains("repairRate") && !nf.at("repairRate").is_null()) {
        s.repair = dist_from_json<T>(nf.at("repairRate"));
        s.has_repair = true;
    }
    if (nf.contains("breakdownResetPolicy") && !nf.at("breakdownResetPolicy").is_null()) {
        const std::string p = nf.at("breakdownResetPolicy").get<std::string>();
        if (!p.empty()) s.breakdown_reset = p;
    }
    if (nf.contains("repairResetPolicy") && !nf.at("repairResetPolicy").is_null()) {
        const std::string p = nf.at("repairResetPolicy").get<std::string>();
        if (!p.empty()) s.repair_reset = p;
    }
    return s;
}

/**
 * The MACRO form: one base stage plus a `nodeFailures` block that expands into
 * the UP and DOWN_<node> stages and the arcs between them.
 *
 * THE DECLARED STAGE NAME IS DROPPED, and deliberately: the reference's
 * `addNodeBreakdown` names the stage it creates `UP` whatever the base model was
 * called, and the repair arcs are then found by that name. Keeping the declared
 * name would leave `addNodeRepair`'s lookup with nothing to find, on both sides.
 */
template <class T>
env::Environment<T> expand_node_failures(const json& model, const json& stage0,
                                         const std::vector<NodeFailureSpec<T> >& fails) {
    if (model.at("stages").size() != 1)
        throw InputError(
            "environment_reader: 'nodeFailures' expands the base model into the UP and "
            "DOWN_<node> stages, so 'stages' must declare exactly one stage, holding the base "
            "(UP) model");
    if (model.contains("transitions") && !model.at("transitions").empty())
        throw InputError(
            "environment_reader: 'nodeFailures' implies the breakdown and repair transitions; "
            "'transitions' must not be declared alongside it");
    if (!stage0.contains("model") || stage0.at("model").is_null())
        throw InputError(
            "environment_reader: 'nodeFailures' requires the base stage to carry a 'model'");

    qn::Network<T> base_net = build_network_from_json<T>(stage0.at("model"));
    const qn::NetworkStruct<T> base = base_net.get_struct();
    const std::size_t E = 1 + fails.size();
    // Either count is a consistent statement of the same file: `numStages` may
    // count the stages as DECLARED (one) or as EXPANDED, and both writers emit
    // the expanded form, where the question does not arise. Anything else means
    // the file counts stages this expansion would not produce.
    if (model.contains("numStages")) {
        const std::size_t declared = model.at("numStages").get<std::size_t>();
        if (declared != E && declared != 1)
            throw InputError("environment_reader: 'numStages' is " + std::to_string(declared) +
                             " but the 'nodeFailures' block expands one base stage into " +
                             std::to_string(E) + " stages");
    }

    env::Environment<T> e(model.value("name", std::string("env")), E);
    for (std::size_t k = 0; k < fails.size(); ++k) {
        const NodeFailureSpec<T>& nf = fails[k];
        e.add_node_breakdown(0, k + 1, base, nf.node, nf.breakdown, nf.down_service,
                             nf.breakdown_reset);
        if (nf.has_repair) e.add_node_repair(nf.node, nf.repair, nf.repair_reset);
    }
    return e;
}

}  // namespace detail

/**
 * Build an `env::Environment<T>` from a parsed model.json envelope.
 *
 * `init()` is NOT called here. It superposes the arcs into the marked processes
 * the analyzer integrates against, and it throws when a stage has no finite
 * holding time; that diagnostic belongs to the solve, next to the options that
 * chose it, and not to the parse. Every caller in this port calls `init()`
 * itself, exactly as the reference's `env.init()` is a separate step from
 * building the object.
 */
template <class T>
env::Environment<T> build_environment_from_json(const detail::json& root) {
    using detail::json;
    const json& model = root.contains("model") ? root.at("model") : root;
    const std::string mtype = model.value("type", std::string(""));
    if (mtype != "Environment")
        throw UnsupportedError("environment_reader: model type '" + mtype +
                               "' is not an Environment; -s env solves a random-environment "
                               "model and the Network path solves the rest");
    detail::reject_unconsumed_env_keys(model);

    if (!model.contains("stages"))
        throw InputError("environment_reader: the environment declares no 'stages'");
    const json& stages = model.at("stages");
    const std::size_t D = stages.size();
    if (D == 0) throw InputError("environment_reader: the environment declares no stage");
    std::vector<std::string> declared_names(D);
    for (std::size_t s = 0; s < D; ++s)
        declared_names[s] = stages[s].value("name", std::string("Stage") + std::to_string(s + 1));

    // The node-failure block, and which of its two roles applies. The rule is
    // both readers': the block is a MACRO to expand only while no DOWN_<node>
    // stage has been declared for it; once one has, the stages already carry
    // the structure and the block carries only the reset policies.
    std::vector<detail::NodeFailureSpec<T> > fails;
    if (model.contains("nodeFailures"))
        for (const json& nf : model.at("nodeFailures"))
            fails.push_back(detail::node_failure_fields<T>(nf));
    bool macro = !fails.empty();
    for (std::size_t k = 0; k < fails.size() && macro; ++k) {
        const std::string down = env::Environment<T>::down_stage_name(fails[k].node);
        for (std::size_t s = 0; s < D; ++s)
            if (declared_names[s] == down) macro = false;
    }
    for (std::size_t k = 0; k < fails.size(); ++k)
        for (std::size_t j = k + 1; j < fails.size(); ++j)
            if (fails[k].node == fails[j].node)
                throw InputError("environment_reader: 'nodeFailures' declares node '" +
                                 fails[k].node +
                                 "' twice, and a node has ONE down stage; merge the two entries");

    if (macro) return detail::expand_node_failures<T>(model, stages[0], fails);

    const std::size_t E = D;
    // `numStages` is written by both writers and is the count the transitions
    // are indexed against, so a disagreement with the array length means the
    // file is inconsistent and the arcs cannot be trusted to point where they
    // say. Checking it costs nothing and turns a wrong answer into a message.
    if (model.contains("numStages")) {
        const std::size_t declared = model.at("numStages").get<std::size_t>();
        if (declared != E)
            throw InputError("environment_reader: 'numStages' is " + std::to_string(declared) +
                             " but 'stages' holds " + std::to_string(E) + " entries");
    }

    env::Environment<T> e(model.value("name", std::string("env")), E);
    for (std::size_t s = 0; s < E; ++s) {
        const json& st = stages[s];
        const std::string nm = st.value("name", std::string("Stage") + std::to_string(s + 1));
        // A stage with no network is not a stage ENV can solve: the analyzer
        // runs a transient solve per stage per iteration, so the missing model
        // would surface as an empty drift rather than as a missing input.
        if (!st.contains("model"))
            throw InputError("environment_reader: stage '" + nm +
                             "' carries no 'model'; every stage of a random environment holds "
                             "the network in force while it lasts");
        qn::Network<T> net = build_network_from_json<T>(st.at("model"));
        e.set_stage(s, nm, st.value("type", std::string("")), net.get_struct());
    }

    if (model.contains("transitions")) {
        for (const json& tr : model.at("transitions")) {
            if (!tr.contains("from") || !tr.contains("to"))
                throw InputError(
                    "environment_reader: an environment transition is missing 'from' or 'to'");
            const long long from = tr.at("from").get<long long>();
            const long long to = tr.at("to").get<long long>();
            if (from < 0 || to < 0 || static_cast<std::size_t>(from) >= E ||
                static_cast<std::size_t>(to) >= E)
                throw InputError("environment_reader: transition " + std::to_string(from) + " -> " +
                                 std::to_string(to) +
                                 " names a stage outside 0.." + std::to_string(E - 1) +
                                 " (the wire index is zero-based in both writers)");
            if (!tr.contains("distribution"))
                throw InputError("environment_reader: the transition " + std::to_string(from) +
                                 " -> " + std::to_string(to) +
                                 " carries no 'distribution'; the arc holding time is what the "
                                 "environment process is made of");
            e.add_transition(static_cast<std::size_t>(from), static_cast<std::size_t>(to),
                             detail::dist_from_json<T>(tr.at("distribution")));
        }
    }

    // Re-attach the descriptors and their reset policies to the stages just
    // built, so the environment is the one it was written from: the arcs came
    // off the wire, the policies could not, and only these entries carry them.
    for (std::size_t k = 0; k < fails.size(); ++k) {
        const detail::NodeFailureSpec<T>& nf = fails[k];
        e.register_node_failure(nf.node, nf.breakdown, nf.down_service, nf.has_repair, nf.repair,
                                nf.breakdown_reset, nf.repair_reset);
    }
    return e;
}

/** Parse a model.json file into an `env::Environment<T>`. */
template <class T>
env::Environment<T> read_environment_json(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in) throw InputError("environment_reader: cannot open " + path);
    detail::json root;
    try {
        in >> root;
    } catch (const detail::json::parse_error& err) {
        throw InputError("environment_reader: malformed JSON in " + path + ": " + err.what());
    }
    return build_environment_from_json<T>(root);
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_ENVIRONMENT_READER_H
