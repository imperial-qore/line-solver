"""
SSA Solver handler.

Native Python SSA (Stochastic Simulation Algorithm) entry point. The simulation
core lives in :mod:`line_solver.api.solvers.ssa.serial`, an afterEvent/sync-action
Gillespie engine that mirrors MATLAB ``solver_ssa.m`` and drives the same
``State.afterEvent`` machinery as the CTMC solver (full parity across
disciplines, including order-/phase-dependent state, PAS, and SPN). This module
keeps the public result/options dataclasses and the method dispatcher; the
``solver_ssa_basic``/``solver_ssa_with_cache``/``solver_ssa_parallel`` names are
retained as thin delegators to the unified engine.
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional, Dict, List, Tuple

from ...sn import NetworkStruct, NodeType


@dataclass
class SolverSSAOptions:
    """Options for SSA solver."""
    method: str = 'default'
    tol: float = 1e-6
    verbose: bool = False
    samples: int = 10000          # Number of simulated events (fired transitions), one state sample per event
    warmupfrac: float = 0.0       # Warmup discard: drop the first floor(warmupfrac*samples) events from the steady-state tallies (0 = disabled)
    timespan: Tuple[float, float] = (0.0, float('inf'))
    seed: int = 0                 # Random seed for reproducibility
    cutoff: float = np.inf        # Open-class truncation (inf = no truncation)
    confidence_level: float = 0.95
    record_events: bool = False
    timeout: float = float('inf')  # Wall-clock time budget in seconds (inf = no budget)


@dataclass
class SolverSSAReturn:
    """Result of the SSA solver handler (M x K metric matrices)."""
    Q: Optional[np.ndarray] = None
    U: Optional[np.ndarray] = None
    R: Optional[np.ndarray] = None
    T: Optional[np.ndarray] = None
    A: Optional[np.ndarray] = None
    C: Optional[np.ndarray] = None
    X: Optional[np.ndarray] = None
    Q_ci: Optional[np.ndarray] = None
    U_ci: Optional[np.ndarray] = None
    R_ci: Optional[np.ndarray] = None
    T_ci: Optional[np.ndarray] = None
    total_time: float = 0.0
    runtime: float = 0.0
    method: str = "default"
    samples: int = 0
    event_log: Optional[List] = None
    state_log: Optional[List] = None
    timedOut: bool = False


def _has_cache_nodes(sn: NetworkStruct) -> bool:
    """Check if the network has any cache nodes."""
    if not hasattr(sn, 'nodetype') or sn.nodetype is None:
        return False
    return any(nt == NodeType.CACHE for nt in sn.nodetype)


# ---------------------------------------------------------------------------
# Public engine entry points (thin delegators to the afterEvent serial engine).
# Cache and SPN models are handled natively by the unified engine, so no
# separate cache/aggregate code path is required.
# ---------------------------------------------------------------------------

def solver_ssa_basic(sn, options=None, model=None) -> SolverSSAReturn:
    """Serial afterEvent Gillespie simulation."""
    from .serial import solver_ssa_run
    if options is None:
        options = SolverSSAOptions()
    return solver_ssa_run(sn, options, method='serial')


def solver_ssa_with_cache(sn, options=None, model=None) -> SolverSSAReturn:
    """Cache networks are simulated by the same engine (afterEventCache)."""
    return solver_ssa_basic(sn, options, model)


def solver_ssa_parallel(sn, options=None) -> SolverSSAReturn:
    """Average several seeded serial replicas."""
    from .serial import solver_ssa_parallel as _parallel
    if options is None:
        options = SolverSSAOptions()
    return _parallel(sn, options)


def solver_ssa(
    sn: NetworkStruct,
    options: Optional[SolverSSAOptions] = None,
    model=None,
) -> SolverSSAReturn:
    """Main SSA dispatcher. Routes by ``options.method``.

    - ``default`` / ``serial`` / ``ssa`` -> afterEvent serial engine
    - ``parallel`` / ``para``            -> seeded serial replicas
    - ``nrm``                            -> dedicated next-reaction method
                                            (falls back to serial if unsupported)
    """
    if options is None:
        options = SolverSSAOptions()
    method = (options.method or 'default').lower()

    # see _kb/06-solver-catalog.md (SSA main section) -- record_events is serial-only
    if getattr(options, 'record_events', False):
        return solver_ssa_basic(sn, options, model)

    if method in ('parallel', 'para', 'ssa.parallel'):
        # NRM-eligible model runs on the fast single-run NRM rather than replicated serial.
        from .nrm import solver_ssa_nrm, _nrm_eligible
        if _nrm_eligible(sn):
            return solver_ssa_nrm(sn, options)
        return solver_ssa_parallel(sn, options)
    if method == 'nrm':
        from .nrm import (solver_ssa_nrm, _ALLOWED_SCHED, _fcr_nrm_ok,
                          _routing_nrm_ok, _impatience_nrm_ok, _phase_nrm_ok,
                          _cache_nrm_ok)
        try:
            supported = all(
                (sn.sched.get(ist) if isinstance(sn.sched, dict) else sn.sched[ist]) in _ALLOWED_SCHED
                for ist in range(sn.nstations))
        except Exception:
            supported = False
        # see _kb/06-solver-catalog.md (SSA: "NRM engine now supports FCR directly")
        if supported and not _fcr_nrm_ok(sn):
            supported = False
        # NRM resolves JSQ/memoryless SQ at firing time; others need serial.
        if supported and not _routing_nrm_ok(sn):
            supported = False
        # QUEUE_LENGTH balking and memoryless reneging only
        if supported and not _impatience_nrm_ok(sn):
            supported = False
        # see _kb/06-solver-catalog.md (SSA main section) -- NRM phase expansion scope
        if supported and not _phase_nrm_ok(sn):
            supported = False
        # Cache nodes are simulated natively (immediate read -> hit/miss class
        # switch); only the retrieval (delayed-hit) system still needs serial.
        if supported and not _cache_nrm_ok(sn):
            supported = False
        if supported:
            return solver_ssa_nrm(sn, options)
        return solver_ssa_basic(sn, options, model)
    if method == 'default':
        # see _kb/06-solver-catalog.md (SSA main section) -- NRM eligibility scope
        from .nrm import solver_ssa_nrm, _nrm_eligible
        if _nrm_eligible(sn):
            return solver_ssa_nrm(sn, options)
    if method not in ('default', 'serial', 'ssa') and getattr(options, 'verbose', False):
        print(f"Warning: Unknown SSA method '{method}'. Using serial.")
    return solver_ssa_basic(sn, options, model)


__all__ = [
    'solver_ssa',
    'solver_ssa_basic',
    'solver_ssa_parallel',
    'solver_ssa_with_cache',
    'SolverSSAReturn',
    'SolverSSAOptions',
    '_has_cache_nodes',
]
