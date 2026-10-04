/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_LINEMODEL_H
#define LINE_IO_LINEMODEL_H

/**
 * @file
 * @ingroup line_io
 * Port of `linemodel_load.m` / `linemodel_save.m`: ONE entry point for the LINE
 * `model.json` interchange, whatever model kind the document carries.
 *
 * The readers stay split by kind (`network_reader.h`, `lqn_json_reader.h`,
 * `workflow_reader.h`, `environment_reader.h`); `linemodel_load` parses the file
 * once and hands the tree to the builder the envelope's `model.type` names, as
 * the reference's `switch mtype` does. The result is a `std::variant` over the
 * four MODEL-LEVEL kinds. A LayeredNetwork comes back as the `lqn::LqnModel`
 * rather than the finalized `LqnStruct`, because the model is what the reference
 * returns and what `write_lqnx` accepts; `lqn::lqn_finalize` turns it into the
 * struct SolverLN solves.
 *
 * SAVING covers the same four kinds: a Network through `network_writer.h`, the
 * other three through `linemodel_writer.h`, whose header says which optional
 * spellings it takes and which undeclarable LqnModel fields it refuses.
 */

#include <fstream>
#include <string>
#include <utility>
#include <variant>

#include "json.hpp"
#include "line/io/environment_reader.h"
#include "line/io/linemodel_writer.h"
#include "line/io/lqn_json_reader.h"
#include "line/io/network_reader.h"
#include "line/io/network_writer.h"
#include "line/io/workflow_reader.h"
#include "line/lang/lqn/lqn_reader.h"
#include "line/lang/qn/environment.h"
#include "line/lang/qn/network_builder.h"
#include "line/lang/workflow/workflow.h"
#include "line/util/error.h"

namespace line {
namespace io {

/** What `linemodel_load` returns: MATLAB's Network, LayeredNetwork, Workflow or Environment. */
template <class T>
using LineModel =
    std::variant<qn::Network<T>, lqn::LqnModel<T>, workflow::Workflow<T>, env::Environment<T>>;

/** The envelope's `model.type` of a parsed document; empty when it declares none. */
inline std::string linemodel_type(const nlohmann::json& root) {
    const nlohmann::json& model = root.contains("model") ? root.at("model") : root;
    if (!model.is_object() || !model.contains("type") || !model.at("type").is_string())
        return std::string();
    return model.at("type").get<std::string>();
}

/** Port of `linemodel_load(filename)` on an already parsed document. */
template <class T>
LineModel<T> linemodel_from_json(const nlohmann::json& root) {
    const std::string mtype = linemodel_type(root);
    if (mtype == "Network") return LineModel<T>(build_network_from_json<T>(root));
    if (mtype == "LayeredNetwork") return LineModel<T>(build_lqn_model_from_json<T>(root));
    if (mtype == "Workflow") return LineModel<T>(build_workflow_from_json<T>(root));
    if (mtype == "Environment") return LineModel<T>(build_environment_from_json<T>(root));
    if (mtype.empty())
        throw InputError("linemodel_load: the document's model declares no 'type'; expected "
                         "Network, LayeredNetwork, Workflow or Environment");
    throw UnsupportedError("linemodel_load: Unsupported model type: " + mtype);
}

/**
 * Port of `linemodel_load(filename)`.
 *
 * @return the model, as whichever alternative of `LineModel<T>` its `model.type` names
 */
template <class T>
LineModel<T> linemodel_load(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in) throw InputError("linemodel_load: cannot open " + path);
    nlohmann::json root;
    try {
        in >> root;
    } catch (const nlohmann::json::parse_error& e) {
        throw InputError("linemodel_load: malformed JSON in " + path + ": " + e.what());
    }
    return linemodel_from_json<T>(root);
}

/** Port of `linemodel_save(model, filename)` for a refreshed Network struct. */
template <class T>
void linemodel_save(const qn::NetworkStruct<T>& sn, const std::string& path) {
    write_network_json(sn, path);
}

/** Port of `linemodel_save(model, filename)` for a Network; refreshes it first. */
template <class T>
void linemodel_save(qn::Network<T>& model, const std::string& path) {
    write_network_json(model.get_struct(), path);
}

/** Port of `linemodel_save(model, filename)` for a LayeredNetwork (`layered2json`). */
template <class T>
void linemodel_save(const lqn::LqnModel<T>& model, const std::string& path) {
    write_linemodel_json(linemodel_envelope(lqn_model_to_json(model)), path);
}

/** Port of `linemodel_save(model, filename)` for a Workflow (`workflow2json`). */
template <class T>
void linemodel_save(const workflow::Workflow<T>& model, const std::string& path) {
    write_linemodel_json(linemodel_envelope(workflow_to_json(model)), path);
}

/** Port of `linemodel_save(model, filename)` for an Environment (`environment2json`). */
template <class T>
void linemodel_save(const env::Environment<T>& model, const std::string& path) {
    write_linemodel_json(linemodel_envelope(environment_to_json(model)), path);
}

/** `linemodel_save` on whatever `linemodel_load` returned, dispatched on the alternative held. */
template <class T>
void linemodel_save(LineModel<T>& model, const std::string& path) {
    std::visit([&path](auto& m) { linemodel_save(m, path); }, model);
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_LINEMODEL_H
