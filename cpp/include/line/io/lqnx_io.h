/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_LQNX_IO_H
#define LINE_IO_LQNX_IO_H

/**
 * @file
 * @ingroup line_io
 * The `.lqnx` (LQNS XML) interchange as a symmetric MODEL-LEVEL pair in `io::`.
 *
 * The implementations live in `lang/lqn/`: `lqn::read_lqnx_model` returns the
 * `LqnModel` a builder produces, `lqn::write_lqnx` takes one, and
 * `lqn::read_lqnx` returns the FINALIZED `LqnStruct` SolverLN consumes. The
 * struct reader is the odd one out of that trio: its output cannot be handed
 * back to the writer. So the io pair is `read_lqnx_model` / `write_lqnx`,
 * both on `LqnModel`, keeping the `lqn::` names verbatim so a caller moving
 * between the namespaces meets no second spelling. MATLAB spells the same pair
 * `LayeredNetwork.parseXML` / `writeXML`, which name methods of a class this
 * port does not have.
 */

#include <string>

#include "line/lang/lqn/lqn_reader.h"
#include "line/lang/lqn/lqn_writer.h"

namespace line {
namespace io {

/** Read a `.lqnx` file into the model-level `LqnModel`, which `write_lqnx` accepts back. */
template <class T>
lqn::LqnModel<T> read_lqnx_model(const std::string& path) {
    return lqn::read_lqnx_model<T>(path);
}

/**
 * Write an `LqnModel` as `.lqnx`; forwards to `lqn::write_lqnx`.
 *
 * @return what the schema could not carry (see `lqn::LqnWriteReport`)
 */
template <class T>
lqn::LqnWriteReport write_lqnx(const lqn::LqnModel<T>& m, const std::string& path,
                               const std::string& model_name = std::string("LQN"),
                               bool use_abstract_names = false, std::string* out_text = NULL) {
    return lqn::write_lqnx(m, path, model_name, use_abstract_names, out_text);
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_LQNX_IO_H
