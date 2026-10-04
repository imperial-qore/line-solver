"""Seed the layered fixed point with the Method of Layers.

`config['warmstart']='mol'` runs :func:`lqn_mol`, the compact Rolia-Sevcik
Method of Layers on the SRVN decomposition, and writes its converged state onto
`SolverLN`'s iterate. It is the second starting point this solver takes, beside
`config['warmstart']='nlp'`, and the two differ in HOW they reach the same law
rather than in which law they reach: both price a station by QD-AMVA, the
program by rooting its stationarity conditions in one shot and the Method of
Layers by sweeping software then hardware submodels until the state settles.

WHY IT TRANSFERS CLEANLY. `lqn_mol` carries exactly the four vectors the iterate
is made of, under the same names: `servt`, `residt`, `callservt` and `thinkt`.
Nothing has to be re-derived from a different parameterization, as the program's
answer does; the only translation is SolverLN's own surrogate-delay law, which
is applied here to `lqn_mol`'s busy-thread counts so the delays land in the
convention `update_think_times` maintains rather than in the one a submodel's
`Z` vector uses.

SCOPE IS NARROWER THAN THE PROGRAM'S. `lqn_mol` refuses activity graphs, second
phases, asynchronous and forwarding calls, caches, setup times, replication,
admission constraints and open arrivals, naming the element it refuses. On the
`doc/latex/exdata` suite that is 9 of the 20 models that parse, against 17 for
the program. Where both serve, this one is the cheaper seed to build.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import math
from typing import Any, Dict

import numpy as np

from ...api.io.logging import line_debug
from ...api.lqn import lqn_mol
from ...constants import GlobalConstants, SchedStrategy
from ...distributions import Exp

__all__ = ['ln_mol_warm_start']


def _mol_options(solver) -> Dict[str, Any]:
    """The Method of Layers settings, read off the LN options.

    The seed is built to ITS OWN tolerance, not to the layered run's. A warm
    start is worth having only if it lands near the fixed point, and `iter_tol`
    on `SolverLN` is a stopping rule for a different iteration whose default of
    5e-3 would hand over a state still two digits out. `lqn_mol` is cheap enough
    that its own 1e-6 costs a fraction of one layered sweep.
    """
    cfg = getattr(getattr(solver, 'options', None), 'config', None)
    cfg = cfg if isinstance(cfg, dict) else {}
    return {
        'iter_max': int(cfg.get('mol_iter_max', 200)),
        'iter_tol': float(cfg.get('mol_iter_tol', 1e-6)),
        'relax_factor': float(cfg.get('mol_relax_factor', 0.5)),
    }


def ln_mol_warm_start(solver) -> float:
    """Seed the layered fixed point with the Method of Layers' solution.

    This is not a method but a STARTING POINT. The layers, the layer solvers and
    the iteration stay the layered method's own; only the state they begin from
    changes, from the model's bare host demands to `lqn_mol`'s answer. The
    iteration then runs normally and converges to ITS fixed point, so the answer
    is the layered one and only the path to it is short.

    `converged()` drops `iter_min` for a run seeded this way, the same as for the
    program's seed: that floor keeps a cold iterate from reading an early plateau
    as convergence, and a warm one starts past the plateau.

    Args:
        solver: the SolverLN whose iterate is to be seeded, already constructed

    Returns:
        `lqn_mol`'s final residual, for the caller to log

    Raises:
        ValueError: the model uses a feature `lqn_mol` does not model, named
    """
    lqn = solver.lqn
    QN, UN, _RN, TN, info = lqn_mol(lqn, _mol_options(solver))
    _seed_times(solver, info, UN, TN)
    _seed_think_times(solver, QN, TN)
    # the same laws the tail of post() pushes, from the Method of Layers instead
    # of from a layered sweep. update_routing_probabilities is NOT among them,
    # for the reason the program's seed skips it too: the call frequencies it
    # writes are already right at construction and it reads a layer result that
    # does not exist yet. Iteration 1 runs it normally.
    solver.update_layers(0)
    solver._refresh_ensemble()
    solver.warmstarted = True
    return float(info['resid'])


def _seed_times(solver, info, UN, TN):
    """Write the service, residence, throughput and utilization onto the iterate.

    ENTRY-ONLY IS WHAT MAKES THIS A COPY. Each entry binds one activity executed
    once per invocation, so the activity's service time is the entry's whole
    execution, its residence is the entry's host residence with a visit count of
    one, and its utilization is the one `lqn_mol` already reports for the entry,
    busy SERVERS at a delay processor and a busy FRACTION at a finite one. On a
    model with an activity graph none of those equalities hold, which is why
    `lqn_mol` refuses one rather than approximating it.
    """
    lqn = solver.lqn
    servt, residt, callservt = info['servt'], info['residt'], info['callservt']

    for eidx in range(int(lqn.eshift), int(lqn.eshift) + int(lqn.nentries)):
        aidx = int(lqn.actsof[eidx][0])
        x = float(TN[eidx])

        solver.servt[eidx] = solver.residt[eidx] = float(servt[eidx])
        solver.tput[eidx] = x

        solver.servt[aidx] = float(servt[eidx])
        solver.residt[aidx] = float(residt[eidx])
        solver.tput[aidx] = x
        solver.util[aidx] = float(UN[eidx])

    for tidx in range(int(lqn.tshift), int(lqn.tshift) + int(lqn.ntasks)):
        solver.tput[tidx] = float(TN[tidx])

    for cidx in range(int(lqn.ncalls)):
        y = float(lqn.callpair[cidx, 2]) if lqn.callpair.shape[1] > 2 else 1.0
        y = 0.0 if math.isnan(y) else y
        solver.callservt[cidx] = float(callservt[cidx])
        solver.callresidt[cidx] = y * float(callservt[cidx])

    for i in range(int(lqn.nidx)):
        if solver.servt[i] > 0:
            solver.servtproc[i] = Exp.fit_mean(float(solver.servt[i]))
    for cidx in range(int(lqn.ncalls)):
        if solver.callservt[cidx] > 0:
            solver.callservtproc[cidx] = Exp.fit_mean(float(solver.callservt[cidx]))


def _seed_think_times(solver, QN, TN):
    """Write the surrogate client delays onto the iterate.

    These are what a starting point is MADE of: a seed of the service times
    alone moves the iterate hardly at all, because the delays are the rest of a
    caller's cycle and run orders of magnitude larger than a service time.
    `update_think_times` will not supply them here, its assignment sitting behind
    a `len(self.results) > 0` guard with no layer yet solved.

    `lqn_mol` KEEPS ITS OWN `thinkt` IN A DIFFERENT CONVENTION: it is the `Z` a
    submodel hands `pfqn_qdamva`, which for a reference task is the declared
    think time rather than a surrogate. So the law applied here is SolverLN's
    own, idle threads over throughput less the declared think time, evaluated at
    `lqn_mol`'s busy-thread counts. `QN` at a task index is exactly that count,
    `sum_e X(e) * S(e)` by Little's law over the entries the task serves.
    """
    lqn = solver.lqn
    for tidx in range(int(lqn.tshift), int(lqn.tshift) + int(lqn.ntasks)):
        x = float(TN[tidx])
        if x <= GlobalConstants.Zero:
            continue
        busy = float(QN[tidx])
        z = _think_of(lqn, tidx)
        mult = _mult_of(lqn, tidx)
        # the POPULATION IS THE LAYER'S, not the task's multiplicity: an
        # infinite-server task holds as many threads as can reach it, and reading
        # `inf` off the task here would leave the delay undefined
        njobs = float(np.max(solver.njobs[tidx, :]))
        if solver._get_sched(tidx) == SchedStrategy.INF or not math.isfinite(mult):
            val = (njobs - busy) / x - z
        else:
            val = njobs * abs(1.0 - busy / max(1.0, mult)) / x - z
        solver.thinkt[tidx] = max(GlobalConstants.Zero,
                                  max(0.0, val) + solver._setup_charge(tidx))
        if solver.thinkt[tidx] + z > 0:
            solver.thinktproc[tidx] = Exp.fit_mean(solver.thinkt[tidx] + z)


def _mult_of(lqn, idx: int) -> float:
    """The multiplicity of element IDX, with a non-finite value meaning INF."""
    mult = np.asarray(getattr(lqn, 'mult', None)).ravel()
    if idx >= mult.size:
        return 1.0
    m = float(mult[idx])
    return 1.0 if m == 0.0 else m


def _think_of(lqn, tidx: int) -> float:
    """The DECLARED think time of task TIDX, zero for a non-reference task."""
    d = (getattr(lqn, 'think', None) or {}).get(tidx)
    if d is None:
        return 0.0
    if isinstance(d, (int, float, np.floating)):
        z = float(d)
    else:
        z = float(d.getMean()) if hasattr(d, 'getMean') else 0.0
    return 0.0 if (not math.isfinite(z) or z < 0.0) else z
