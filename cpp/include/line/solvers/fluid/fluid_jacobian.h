/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_FLUID_JACOBIAN_H
#define LINE_SOLVERS_FLUID_FLUID_JACOBIAN_H

/**
 * Port of `@@SolverFLD/getJacobian`, all four of its outputs.
 *
 * The reference builds the drift locally with `getSymbolicDrift` and then posts
 * it to the line-sage-rest backend, asking for the Jacobian and, when a fourth
 * output is requested, the solutions of f(x) = 0. This header does the same:
 * `fluid_symbolic_drift` supplies rhs and vars, `sym::sym_resolve` supplies the
 * backend, and `SymEngine::fluidODEs` answers.
 *
 * WHY THE BACKEND IS NOT OPTIONAL FOR THE EQUILIBRIA. Differentiating the drift
 * is structural -- `fluid_symbolic_jacobian` does it exactly by the chain rule
 * over the typed factors of `FluidSymSystem`, with no algebra system in sight --
 * but SOLVING f(x) = 0 in closed form is not. There is no local answer to fall
 * back to, so a request for equilibria without a backend is refused by name,
 * with the same guidance the reference's `SAGE.require()` prints. Returning an
 * empty list instead would read as "this system has no equilibria", which is a
 * different and false statement.
 *
 * DIVERGENCE, DELIBERATE: when equilibria are NOT asked for and no backend
 * resolves, this returns the locally differentiated Jacobian rather than
 * erroring as the reference does. The entries are the same derivative, obtained
 * without a 3 GB container; refusing to answer a question we can answer exactly
 * would be a worse port than answering it. `FluidJacobian::engine` names which
 * path produced them, so a caller comparing text against MATLAB knows whether it
 * is comparing against SAGE's normal form or this port's.
 *
 * The smoothness gate is upstream, in `fluid_symbolic_drift`: a drift that
 * scales rates by min(n_i, S_i) has no Jacobian at n_i = S_i and is refused
 * before any backend is contacted.
 */

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "line/api/sym/sage_rest_engine.h"
#include "line/api/sym/sym_engine.h"
#include "line/api/sym/sym_engines.h"
#include "line/solvers/fluid/fluid_symodes.h"
#include "line/util/error.h"

namespace line {
namespace fluid {

/** Backend selection, mirroring `options.config.symbolic` and its timeout. */
struct FluidSymbolicOptions {
    /** `auto` to search, a URL, an image name, or `none` to stay local. */
    std::string backend = "auto";
    /** `options.config.symbolic_timeout`, seconds; the reference defaults to 300. */
    int timeout_s = 300;
    /** The reference's `nargout >= 4`: ask the backend to solve f(x) = 0. */
    bool equilibria = false;
};

/** The four outputs of `@@SolverFLD/getJacobian`. */
struct FluidJacobian {
    std::vector<std::string> vars;             ///< state variable names
    std::vector<std::string> rhs;              ///< the drift, one expression per variable
    std::vector<std::vector<std::string> > J;  ///< J[i][j] = d f_i / d x_j
    /**
     * Solutions of f(x) = 0, each a variable -> expression map.
     *
     * EMPTY IS NOT AN ASSERTION THAT NONE EXIST. The backend returns what it
     * solves in closed form, and a system beyond it answers with nothing; that
     * is a limit of the solve. `has_equilibria` separates "asked and answered
     * with none" from "never asked".
     */
    std::vector<std::map<std::string, std::string> > equilibria;
    bool has_equilibria = false;  ///< the backend answered the equilibria request
    std::string engine;           ///< `sage` or `local`, whichever produced J
};

namespace detail {

/** The reference's `SAGE.require()` message, verbatim in substance. */
inline std::string fluid_sym_backend_missing() {
    return std::string(
        "fluid_jacobian: solving f(x) = 0 in closed form needs a symbolic backend, and none is "
        "available. Start one with\n  docker run -d -p 8080:8080 ") +
        sym::SYM_DOCKER_IMAGE +
        "\npoint the " + sym::SYM_URL_ENV +
        " environment variable at a running service, or set the backend to its URL. The Jacobian "
        "alone needs no backend: ask for it without the equilibria";
}

}  // namespace detail

/**
 * Jacobian, drift and equilibria of the mean-field vector field.
 *
 * @param sys the symbolic system, from `fluid_symodes`
 * @param opt backend selection and whether equilibria are wanted
 * @return the four outputs; `engine` names what produced the Jacobian
 */
inline FluidJacobian fluid_jacobian(const FluidSymSystem& sys,
                                    const FluidSymbolicOptions& opt = FluidSymbolicOptions()) {
    // The drift comes first and locally, exactly as the reference's getJacobian
    // calls getSymbolicDrift before touching the backend: it is what carries the
    // smoothness gate, so a non-smooth method is refused without a round trip.
    const FluidSymbolicDrift d = fluid_symbolic_drift(sys);

    FluidJacobian out;
    out.vars = d.vars;
    out.rhs = d.rhs;

    std::shared_ptr<sym::SymEngine> engine = sym::sym_resolve(opt.backend);
    if (!engine) {
        if (opt.equilibria) throw sym::SymEngineError(detail::fluid_sym_backend_missing());
        const FluidSymbolicJacobian local = fluid_symbolic_jacobian(sys);
        out.J = local.J;
        out.engine = "local";
        return out;
    }

    // The timeout travels in the request body, so it must be set on the engine
    // before the call rather than passed to it; only the REST engine has one.
    if (opt.timeout_s > 0) {
        sym::SageRestEngine* rest = dynamic_cast<sym::SageRestEngine*>(engine.get());
        if (rest != nullptr) rest->setTimeoutSeconds(opt.timeout_s);
    }

    std::vector<std::string> want;
    want.push_back("jacobian");
    if (opt.equilibria) want.push_back("equilibria");
    const sym::FluidODEs odes = engine->fluidODEs(d.rhs, d.vars, want);

    if (!odes.hasJacobian)
        throw sym::SymEngineError("fluid_jacobian: backend '" + engine->name() +
                                  "' answered without the jacobian it was asked for");
    if (odes.jacobian.size() != sys.nstates)
        throw sym::SymEngineError("fluid_jacobian: backend returned " +
                                  detail::sym_state_index(odes.jacobian.size()) + " rows for a " +
                                  detail::sym_state_index(sys.nstates) + " state system");
    for (std::size_t i = 0; i < odes.jacobian.size(); ++i)
        if (odes.jacobian[i].size() != sys.nstates)
            throw sym::SymEngineError("fluid_jacobian: backend row " +
                                      detail::sym_state_index(i + 1) + " has " +
                                      detail::sym_state_index(odes.jacobian[i].size()) +
                                      " columns, expected " + detail::sym_state_index(sys.nstates));
    out.J = odes.jacobian;
    out.engine = engine->name();
    if (opt.equilibria) {
        out.equilibria = odes.equilibria;
        out.has_equilibria = odes.hasEquilibria;
    }
    return out;
}

}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_FLUID_JACOBIAN_H
