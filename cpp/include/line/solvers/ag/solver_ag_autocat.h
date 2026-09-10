/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_AG_SOLVER_AG_AUTOCAT_H
#define LINE_SOLVERS_AG_SOLVER_AG_AUTOCAT_H

/**
 * `solver_ag_autocat.m`, the optimization-based search for RCAT product-form
 * solutions. Refused by name, on two independent grounds either of which is
 * sufficient.
 *
 * IT IS ALREADY DEAD IN THE REFERENCE. `SolverMAM.listValidMethods` no longer
 * lists 'exact', and its comment says why: "'exact' method removed - autocat
 * moved to line-legacy.git". `solver_mam_analyzer.m` never calls autocat on any
 * path, and the one remaining reference to it, the 'exact' case of
 * `solver_mam_ag.m`, does not call it either -- it emits a line_warning and
 * falls back to INAP. The .m file still sitting in matlab/src/solvers/MAM is a
 * leftover of that move, not live code. Porting it would put a method into the
 * C++ solver that the reference cannot reach, so the two would disagree by
 * construction on every model.
 *
 * IT NEEDS A MATHEMATICAL PROGRAMMING STACK THIS PORT DOES NOT HAVE. Eleven of
 * its fourteen relaxations bottom out in `linprog` or `fmincon` from the MATLAB
 * Optimization Toolbox: LP relaxation with McCormick envelopes, the tightened
 * LP and zero-potential families with cutting planes, and the interior-point
 * nonlinear programs behind 'ens' and 'qcp'. There is no LP or NLP solver
 * anywhere under cpp/include/line, and writing one to serve a dead code path
 * would be a large piece of numerics whose correctness nothing here could
 * check.
 *
 * The live RCAT path is `solver_mam_ag.h`, which is ported: INAP, INAP+ and the
 * matrix-geometric INAPINF. Callers wanting a product form should use those.
 */

#include <string>

#include "line/lang/qn/network_struct.h"
#include "line/solvers/mva/mva_types.h"
#include "line/solvers/ag/ag_types.h"
#include "line/util/error.h"

namespace line {
namespace ag {

/**
 * The reference signature takes RCAT rate matrices and an action-process map
 * rather than a NetworkStruct, because autocat sits below the network layer.
 * The port keeps the struct-level signature the dispatch would need, since
 * there is no body for the other one to feed.
 *
 * It returns the plain `mva::MvaSolution<T>` the other analyzers return, not a
 * wrapper carrying `actualmethod`: nothing can observe the return of a function
 * that always throws, and the wrapper would leave the dispatch with a third
 * result shape to unpack for no gain. `solver_mam_ag.h`'s `AgResult` stays a
 * wrapper because it carries a payload of its own.
 */
template <class T>
mva::MvaSolution<T> solver_ag_autocat(const qn::NetworkStruct<T>&, const AgOptions&) {
    throw UnsupportedError(
        "SolverAG: solver_ag_autocat is not ported. It is dead in the reference -- "
        "SolverMAM.listValidMethods dropped 'exact' when autocat moved to line-legacy.git, "
        "solver_mam_analyzer never calls it, and solver_mam_ag's 'exact' case warns and falls "
        "back to INAP -- and it needs MATLAB's linprog and fmincon (McCormick-envelope LP "
        "relaxations, cutting planes, interior-point nonlinear programs), for which this port has "
        "no LP or NLP solver. Use method 'inap', 'inapplus' or 'inapinf', which are ported in "
        "solver_mam_ag.h");
}

}  // namespace ag
}  // namespace line

#endif  // LINE_SOLVERS_AG_SOLVER_AG_AUTOCAT_H
