"""
The QD-AMVA analyzer of SolverLN: the whole layered model as ONE program.

Every other SolverLN method decomposes the model into layers, hands each one to
a queueing-network solver and drives the fixed point between them. This one does
not decompose at all. It writes the model's flow, composition and closure
relations as linear constraints, its congestion and Little laws as residuals,
and hands the lot to a root-find. There are no submodels, no layer solvers and
no outer iteration, so none of the machinery `_build_layers` sets up is built.

The program is the one `line_solver.api.lqn.export_nlp` emits, executed rather
than written to a file. That is deliberate: the exported script IS the
implementation, so the method a user reads in a generated file and the method
this analyzer runs cannot drift apart. `config['nlp_source']` writes the same
text out for inspection.

The laws are QD-AMVA (Casale, Perez and Wang, PERFORMANCE 2015), its
arrival-instant estimate included: a station is served at rate r_i(n) with n
jobs present, and a multiserver is the case r_i(n) = softmin(n, m) rather than a
clamp of its own. `config['lld']` declares a further smooth rate multiplier per
station. See the module docstring of `lqn_export_nlp` for the full form.

What this trades away against the layered fixed point: interlocking, phase 2,
forwarding beyond what the flow relations carry, replication, caches, open
arrivals and anything else the algebraic form cannot state. Those are refused BY
NAME when the program is built, not silently dropped.

AND-fork/join is refused by default and carried by `config['nlp_fork']='exp'`,
which prices a forking entry at E[max] of its concurrent branches taken as
exponential. That is an assumption about branch shape rather than a fact the
program carries, and it leaves the task congestion law approximate for a
forking task; see THE AND-FORK CLOSURE in `lqn_export_nlp`.
"""

from typing import Any, Dict, Optional, Tuple

import numpy as np

from ...api.io.logging import line_debug, line_warning
from ...api.lqn import LqnNlpExportError, export_nlp
from ...constants import GlobalConstants, SchedStrategy, VerboseLevel
from ...distributions import Exp

__all__ = ['LqnNlpExportError', 'ln_nlp_program', 'ln_nlp_warm_start',
           'solver_ln_nlp_analyzer']


def ln_nlp_program(lqn, name: Optional[str] = None,
                   lld: Optional[Dict[str, str]] = None,
                   alpha: float = 20.0,
                   fork: str = 'refuse',
                   filename: Optional[str] = None) -> Tuple[str, Dict[str, Any]]:
    """
    Build the nonlinear program of a layered model and make it callable.

    Args:
        lqn: a LayeredNetwork or its LayeredNetworkStruct
        name: model name recorded in the program
        lld: queue-dependent rate multipliers, ``{station name: expression}``
        alpha: sharpness of the softmin standing in for ``min(n, m)``
        fork: ``'refuse'`` or ``'exp'``, how to treat an AND-fork block
        filename: also write the program out, for inspection

    Returns:
        ``(source, namespace)``; the namespace holds ``solve``, ``report``,
        ``objective`` and the model's matrices under the names the emitted
        script gives them.

    Raises:
        LqnNlpExportError: the model uses a feature the form cannot express
    """
    src = export_nlp(lqn, filename=filename, name=name, lld=lld, alpha=alpha,
                     fork=fork)
    # `__name__` is set to anything but '__main__' so the script's own CLI entry
    # point does not run and call sys.exit on the caller.
    ns: Dict[str, Any] = {'__name__': 'line_solver.solver_ln.nlp'}
    exec(compile(src, filename or '<lqn-nlp>', 'exec'), ns)
    return src, ns


def _block(ns, x, name):
    """One named variable block of the program's iterate, as an array.

    The layout is read off the program's own `_BLOCKS`, never restated here: a
    second copy of it would drift the moment the exporter gains a variable.
    """
    off = ns['OFFSET'][name]
    return x[off:off + dict(ns['_BLOCKS'])[name]]


def _nlp_options(solver):
    """The analyzer's settings, read off the LN options."""
    opts = getattr(solver, 'options', None)
    cfg = getattr(opts, 'config', None)
    cfg = cfg if isinstance(cfg, dict) else {}
    lld = cfg.get('lld', None)
    if lld is not None and not isinstance(lld, dict):
        raise LqnNlpExportError("config['lld'] must be a dict of {station name: "
                                "expression in n}, got %r" % type(lld).__name__)
    return {
        'lld': lld,
        'alpha': float(cfg.get('nlp_alpha', 20.0)),
        # Concurrency. 'refuse' declines an AND-fork by name; 'exp' carries it
        # at E[max] over exponential branches. It is off by default because the
        # closure is an assumption about branch shape, not a reading of the
        # model, and the user should be the one making it.
        'fork': str(cfg.get('nlp_fork', 'refuse')),
        # Warm sweeps of the QD-AMVA fixed point taken before the root-find.
        # Zero is a floor, not a cap: solve() escalates on its own when the
        # cold start lands in the wrong basin, so this only sets where the
        # ladder begins.
        'warm': int(cfg.get('nlp_warm', 0)),
        'source': cfg.get('nlp_source', None),
        # `iter_max` bounds the layered FIXED POINT, which this method does not
        # run; the budget that means anything here is the fallback descent's,
        # and LN's own default of 200 is far too small for it. It gets a knob of
        # its own rather than a reading of one whose meaning does not carry.
        'maxiter': int(cfg.get('nlp_maxiter', 3000)),
        # The trust-region trace is a per-iteration table, so it belongs to the
        # DEBUG level and not to the STD one every solver banner prints at.
        'verbose': 2 if (bool(getattr(opts, 'verbose', False))
                         and GlobalConstants.getVerbose() == VerboseLevel.DEBUG) else 0,
    }


def solver_ln_nlp_analyzer(solver) -> Tuple[np.ndarray, ...]:
    """
    Solve a layered model as one nonlinear program and report its averages.

    Args:
        solver: the SolverLN whose ``lqn`` struct is to be analyzed

    Returns:
        ``(QN, UN, RN, TN, AN, WN)``, each indexed by ABSOLUTE LQN index over
        ``0 .. lqn.nidx``, as ``SolverLN.get_ensemble_avg`` returns them

    Raises:
        LqnNlpExportError: the model uses a feature the form cannot express
    """
    lqn = solver.lqn
    cfg = _nlp_options(solver)
    name = str(getattr(getattr(solver, 'model', None), 'name', '') or '')

    src, ns = ln_nlp_program(lqn, name=name, lld=cfg['lld'], alpha=cfg['alpha'],
                             fork=cfg['fork'],
                             filename=cfg['source'])
    line_debug("LN nlp: %d variables, %d equalities, %d inequalities, %d residuals",
               ns['NVAR'], ns['AEQ'].shape[0], ns['AUB'].shape[0], ns['NRES'])

    res = ns['solve'](warm=cfg['warm'], maxiter=cfg['maxiter'], verbose=cfg['verbose'])
    x = np.maximum(np.asarray(res.x, dtype=float), 0.0)
    obj = float(ns['objective'](x))
    # A root of the square system is this method's counterpart of the layered
    # fixed point converging, and is reported through the same flag.
    solver.hasconverged = bool(ns['solved'](x))
    if not solver.hasconverged:
        line_warning("SolverLN", "method='nlp' did not reach a root of the QD-AMVA "
                     "laws (objective %.3e over %d residuals); the reported means "
                     "satisfy the linear relations but not every congestion law. "
                     "solve() already escalated the warm start to 100 sweeps, so "
                     "raise config['nlp_warm'] past that to start it further along."
                     % (obj, ns['NRES']))

    # kept for inspection, the way the layered path keeps its per-layer results
    solver.nlp_source = src
    solver.nlp_solution = x
    solver.nlp_objective = obj

    return _map_to_lqn(solver, ns, x)


def _map_to_lqn(solver, ns, x) -> Tuple[np.ndarray, ...]:
    """
    Place the program's blocks on the LQN index space.

    The program numbers hosts, tasks, entries, activities and calls from zero in
    declaration order, which is the order `_collect` reads them off the struct,
    so element i of a block is absolute index ``shift + i``. The metrics
    themselves follow SolverLN's own table conventions:

        processor   Util only, the sum of the shares its activities hold
        task        Tput, QLen and Util summed over its entries and activities,
                    ResidT the host residence its activities accumulate
        entry       Tput X(e), RespT S(e), QLen X(e)*S(e)
        activity    Tput X(a), RespT W(a) + sum_c y(c)*R(c), ResidT v(a)*W(a)

    A processor has no response time, a task no response time of its own and an
    entry no residence time, and each stays NaN here exactly as the layered path
    leaves it.
    """
    lqn = solver.lqn
    n = int(lqn.nidx)
    QN = np.full(n, np.nan)
    UN = np.full(n, np.nan)
    RN = np.full(n, np.nan)
    TN = np.full(n, np.nan)
    AN = np.full(n, np.nan)     # no arrival-rate metric in the algebraic form
    WN = np.full(n, np.nan)

    hshift = int(getattr(lqn, 'hshift', 0) or 0)
    tshift, eshift, ashift = int(lqn.tshift), int(lqn.eshift), int(lqn.ashift)

    Xe = _block(ns, x, 'Xe')
    Xa = _block(ns, x, 'Xa')
    W = _block(ns, x, 'W')
    S = _block(ns, x, 'S')
    R = _block(ns, x, 'R')

    NH, NT, NE, NA = ns['NH'], ns['NT'], ns['NE'], ns['NA']
    ACT_ENTRY, ACT_TASK = ns['ACT_ENTRY'], ns['ACT_TASK']
    ACT_DEMAND, ACT_VISITS = ns['ACT_DEMAND'], ns['ACT_VISITS']
    TASK_HOST, HOST_DELAY = ns['TASK_HOST'], ns['HOST_DELAY']
    HOST_MAXMULT = ns['HOST_MAXMULT']
    CALLS_OF_ACT, CALL_MEAN = ns['CALLS_OF_ACT'], ns['CALL_MEAN']

    # activities: the leaves everything else is accumulated from
    util_a = np.zeros(NA)
    for a in range(NA):
        h = TASK_HOST[ACT_TASK[a]]
        busy = ACT_DEMAND[a] * Xa[a]
        # an infinite server reports the mean number busy, a finite pool the
        # fraction of its REACHABLE servers that are: HOST_MAXMULT, not the
        # declared count, is what the layered methods divide by, so a pool
        # wider than the work reaching it does not read as idle capacity
        util_a[a] = busy if HOST_DELAY[h] else busy / HOST_MAXMULT[h]
        # the response time of ONE execution: its own host residence plus the
        # response of every call it dispatches
        resp = W[a] + sum(CALL_MEAN[c] * R[c] for c in CALLS_OF_ACT[a])
        aidx = ashift + a
        TN[aidx] = Xa[a]
        RN[aidx] = resp
        QN[aidx] = Xa[a] * resp
        UN[aidx] = util_a[a]
        # residence is per invocation of the entry, so it carries the visits
        WN[aidx] = ACT_VISITS[a] * W[a]

    # entries
    for e in range(NE):
        eidx = eshift + e
        TN[eidx] = Xe[e]
        RN[eidx] = S[e]
        QN[eidx] = Xe[e] * S[e]
        UN[eidx] = sum(util_a[a] for a in range(NA) if ACT_ENTRY[a] == e)

    # tasks
    for t in range(NT):
        tidx = tshift + t
        ents = [e for e in range(NE) if ns['ENTRY_TASK'][e] == t]
        acts = [a for a in range(NA) if ACT_TASK[a] == t]
        TN[tidx] = float(sum(Xe[e] for e in ents))
        QN[tidx] = float(sum(Xe[e] * S[e] for e in ents))
        UN[tidx] = float(sum(util_a[a] for a in acts))
        WN[tidx] = float(sum(ACT_VISITS[a] * W[a] for a in acts))

    # processors
    for h in range(NH):
        hidx = hshift + h
        UN[hidx] = float(sum(util_a[a] for a in ns['ACTS_OF_HOST'][h]))

    # elements outside the workload-anchored component report zero, not NaN,
    # exactly as the layered path zeroes them
    ignore = getattr(solver, 'ignore', None)
    if ignore is not None:
        for idx in range(n):
            if ignore[idx]:
                QN[idx] = UN[idx] = RN[idx] = TN[idx] = AN[idx] = WN[idx] = 0.0

    return QN, UN, RN, TN, AN, WN


def ln_nlp_warm_start(solver) -> float:
    """
    Seed the layered fixed point with the QD-AMVA program's solution.

    This is not a method but a STARTING POINT. The layers, the layer solvers and
    the iteration are the layered method's own; only the laws they begin from
    change, from the model's bare host demands to the program's answer. The
    iteration then runs normally and converges to ITS fixed point, not to the
    program's, so the answer is the layered one and only the path to it is short.

    `converged()` drops `iter_min` for a run seeded this way. That floor exists to
    keep a cold iterate from reading an early plateau as convergence, and a warm
    one starts past the plateau.

    Args:
        solver: the SolverLN whose iterate is to be seeded, already constructed

    Returns:
        the program's objective at its solution, for the caller to log

    Raises:
        LqnNlpExportError: the model uses a feature the form cannot express
    """
    lqn = solver.lqn
    cfg = _nlp_options(solver)
    name = str(getattr(getattr(solver, 'model', None), 'name', '') or '')
    _, ns = ln_nlp_program(lqn, name=name, lld=cfg['lld'], alpha=cfg['alpha'],
                           fork=cfg['fork'])
    res = ns['solve'](warm=cfg['warm'], maxiter=cfg['maxiter'], verbose=cfg['verbose'])
    x = np.maximum(np.asarray(res.x, dtype=float), 0.0)

    _seed_times(solver, ns, x)
    _seed_think_times(solver, ns, x)
    # the same laws the tail of post() pushes, from the program instead of from
    # an iteration: update_routing_probabilities is NOT among them, because the
    # call frequencies it writes are already right at construction and it reads a
    # layer result that does not exist yet. Iteration 1 runs it normally.
    solver.update_layers(0)
    solver._refresh_ensemble()
    solver.warmstarted = True
    return float(ns['objective'](x))


def _seed_times(solver, ns, x):
    """Write the program's service, residence and throughput onto the iterate."""
    lqn = solver.lqn
    Xe, Xa, W, S, R = (_block(ns, x, n) for n in ('Xe', 'Xa', 'W', 'S', 'R'))
    NA, NE, NC, NT = ns['NA'], ns['NE'], ns['NC'], ns['NT']
    y, visits = ns['CALL_MEAN'], ns['ACT_VISITS']

    # An activity's SERVICE time in the iterate is its whole execution, host
    # residence plus the calls it dispatches: a zero-demand call activity still
    # carries the callee's service time, which is what the layered state holds
    # for it. Its RESIDENCE is the host part alone, per invocation of the entry,
    # so it carries the visits.
    for a in range(NA):
        idx = lqn.ashift + a
        solver.servt[idx] = W[a] + sum(y[c] * R[c] for c in ns['CALLS_OF_ACT'][a])
        solver.residt[idx] = visits[a] * W[a]
        solver.tput[idx] = Xa[a]
        h = ns['TASK_HOST'][ns['ACT_TASK'][a]]
        busy = ns['ACT_DEMAND'][a] * Xa[a]
        solver.util[idx] = busy if ns['HOST_DELAY'][h] else busy / ns['HOST_MAXMULT'][h]
    for e in range(NE):
        idx = lqn.eshift + e
        solver.servt[idx] = solver.residt[idx] = S[e]
        solver.tput[idx] = Xe[e]
    for t in range(NT):
        idx = lqn.tshift + t
        solver.tput[idx] = float(sum(Xe[e] for e in range(NE)
                                     if ns['ENTRY_TASK'][e] == t))
    for c in range(NC):
        solver.callservt[c] = R[c]
        solver.callresidt[c] = y[c] * R[c]

    for i in range(lqn.nidx):
        if solver.servt[i] > 0:
            solver.servtproc[i] = Exp.fit_mean(float(solver.servt[i]))
    for c in range(lqn.ncalls):
        if solver.callservt[c] > 0:
            solver.callservtproc[c] = Exp.fit_mean(float(solver.callservt[c]))


def _seed_think_times(solver, ns, x):
    """Write the surrogate client delays onto the iterate.

    These are what a starting point is MADE of, and the reason a seed of the
    service times alone moves nothing: they are the rest of a caller's cycle, so
    on a model like ofbiz they run from 1 to 14 against service times of 0.01.
    `update_think_times` will not supply them here, since its assignment sits
    behind a `len(self.results) > 0` guard and no layer has been solved yet.

    The law is that function's own, idle threads over throughput less the
    declared think time, and the program carries the busy threads directly:
    `BC(c) = X(c)*S(e)` is call c's hold on its CALLEE and `BR(r)` the threads of
    reference task r held.
    """
    lqn = solver.lqn
    BC, BR = _block(ns, x, 'BC'), _block(ns, x, 'BR')
    refs = list(ns['REF_TASKS'])
    INF = ns['INF']

    for t in range(ns['NT']):
        idx = lqn.tshift + t
        X = float(solver.tput[idx])
        if X <= GlobalConstants.Zero:
            continue
        busy = float(BR[refs.index(t)] if t in refs
                     else sum(BC[c] for c in ns['CALLS_INTO_TASK'][t]))
        z = float(ns['TASK_THINK'][t])
        mult = ns['TASK_MULT'][t]
        # the POPULATION IS THE LAYER'S, not the task's multiplicity: an
        # infinite-server task holds as many threads as can reach it, and reading
        # `inf` off the task here would leave the delay undefined
        njobs = float(np.max(solver.njobs[idx, :]))
        if solver._get_sched(idx) == SchedStrategy.INF or mult == INF:
            val = (njobs - busy) / X - z
        else:
            val = njobs * abs(1.0 - busy / float(mult)) / X - z
        solver.thinkt[idx] = max(GlobalConstants.Zero,
                                 max(0.0, val) + solver._setup_charge(idx))
        if solver.thinkt[idx] + z > 0:
            solver.thinktproc[idx] = Exp.fit_mean(solver.thinkt[idx] + z)
