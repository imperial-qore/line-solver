"""Static applicability test for the 'minnormal' fluid method.

Python twin of the MATLAB ``fluid_minnormal_applicable`` and the JAR
``FluidMinNormalApplicable``. Used by the 'default' method of SolverFLD to
prefer 'minnormal' over 'matrix' wherever the second-order closure applies.

The test is STATIC: it inspects the model and the options, never the solution.
A method chosen by trial and rollback would make the reported method depend on
a failed run, and the fallback would silently absorb real defects in the
closure; anything the check cannot decide in advance is left to fail loudly
under an explicit method='minnormal'.

Every condition below mirrors a guard that MinNormalSolver or ClosingMethod
would otherwise raise, so this module and those guards must move together.
"""

from typing import Tuple

from line_solver.api.sn import NodeType, SchedStrategy

from .methods.closing import ClosingMethod


def fluid_minnormal_applicable(sn, options) -> Tuple[bool, str]:
    """Whether the moment-closure method can answer this model.

    Parameters
    ----------
    sn : NetworkStruct
        Network structure, after sn_nonmarkov_toph.
    options : SolverFLDOptions
        Solver options.

    Returns
    -------
    (ok, reason)
        ok is True when 'minnormal' can be selected; reason names the blocking
        feature otherwise and is empty when ok.
    """
    probe = ClosingMethod(sn, options)
    _, _, phases = probe._extract_service_params()
    sched = probe._get_sched()
    M = sn.nstations
    K = sn.nclasses

    # Open and mixed models are supported: the moment terms project the EXT
    # source pool out of the covariance. What they cannot take is a NON-POISSON
    # arrival stream, whose source coordinates track the phase of a single
    # arrival process rather than a population.
    for i in range(M):
        if sched[i] != SchedStrategy.EXT:
            continue
        for r in range(K):
            if phases[i, r] > 1:
                return False, ('class %d has a %d-phase (non-Poisson) arrival process'
                               % (r + 1, phases[i, r]))

    # A cache model is answered through the decomposition path, with the closure
    # in its network step, so cache nodes no longer decline the method. What
    # still declines is a replacement strategy with no drift-based fluid model,
    # mirroring the runtime guard in _solve_rmf: LRU, HLRU, CLIMB and QLRU are
    # answered by a characteristic-time fixed point, which is not a fluid method
    # and has no covariance.
    from ...lang.base import ReplacementStrategy
    # `x or []` is ambiguous when x is an ndarray of more than one element, and
    # nodetype is a list on a Network struct but an array on structs built by
    # other routes, so test for absence explicitly.
    nodetype = getattr(sn, 'nodetype', None)
    if nodetype is None:
        nodetype = []
    _drift_strats = (ReplacementStrategy.RR, ReplacementStrategy.FIFO,
                     ReplacementStrategy.SFIFO)
    for ind, nt in enumerate(nodetype):
        if nt != NodeType.CACHE:
            continue
        ch = sn.nodeparam.get(ind) if getattr(sn, 'nodeparam', None) else None
        strat = getattr(ch, 'replacestrat', None) if ch is not None else None
        if strat is None and isinstance(ch, dict):
            strat = ch.get('replacestrat')
        if strat not in _drift_strats:
            return False, 'a cache uses a replacement strategy with no drift-based fluid model'

    for i in range(M):
        if sched[i] not in ClosingMethod._DRIFT_SCHEDS:
            return False, ('station %d uses %s, which has no fluid drift branch'
                           % (i + 1, sched[i]))

    # the Lyapunov solve is cubic in the phase-resolved state, so the same cap
    # MinNormalSolver enforces decides selection rather than being hit later
    maxstate = getattr(options, 'moment_maxstate', None) or 200
    nstate = int(phases.sum())
    if nstate > maxstate:
        return False, ('the phase-resolved state has %d coordinates, above the '
                       'moment_maxstate limit of %d' % (nstate, maxstate))

    # the moment methods need an autonomous drift
    cfg = getattr(options, 'config', None)
    if isinstance(cfg, dict):
        if cfg.get('rate_traj') is not None:
            return False, 'options.config.rate_traj makes the drift time-varying'
        if cfg.get('nhpp_sched') or cfg.get('rate_sched'):
            return False, 'options.config.nhpp_sched makes the drift time-varying'

    return True, ''
