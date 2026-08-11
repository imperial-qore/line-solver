"""
CTMC Solver handler.

Native Python implementation of CTMC (Continuous-Time Markov Chain) solver
handler that analyzes queueing networks through exact state-space enumeration.

The CTMC solver builds the complete state space and infinitesimal generator
matrix, then solves for steady-state probabilities to compute performance metrics.

Port from:

"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Dict, List, Tuple, Any
from itertools import product
import time

import warnings
from ...sn import (
    NetworkStruct,
    SchedStrategy,
    NodeType,
    RoutingStrategy,
    sn_is_open_model,
    sn_is_closed_model,
    sn_has_open_classes,
    sn_refresh_visits,
)
from .signal_util import signal_lossy_classes, busy_fraction
from ....constants import ProcessType, GlobalConstants, VerboseLevel
from ...mc import ctmc_solve, ctmc_makeinfgen
from ...mc.dtmc import dtmc_stochcomp
from ...cache.miss import cache_xi_fp
from ...mam import map_mean


@dataclass
class SolverCTMCOptions:
    """Options for CTMC solver."""
    method: str = 'default'
    tol: float = 1e-6
    verbose: bool = False
    cutoff: int = 10  # Cutoff for open class populations
    hide_immediate: bool = True  # Hide immediate transitions
    state_space_gen: str = 'default'  # 'default', 'full', 'reachable'
    force: bool = False  # Force solver to run even if state space may be too large
    gen_method: str = 'default'  # 'default' = monolithic builder, 'sync' = sync-action-based builder


@dataclass
class SolverCTMCReturn:
    """
    Result of CTMC solver handler.

    Attributes:
        Q: Mean queue lengths (M x K)
        U: Utilizations (M x K)
        R: Response times (M x K)
        T: Throughputs (M x K)
        C: Cycle times (1 x K)
        X: System throughputs (1 x K)
        pi: Steady-state distribution
        infgen: Infinitesimal generator matrix
        space: State space matrix
        station_col_ranges: List of (start, end) tuples for each station's columns in state space
        rrobin_info: Round-robin and cache state information
        runtime: Runtime in seconds
        method: Method used
    """
    Q: Optional[np.ndarray] = None
    U: Optional[np.ndarray] = None
    R: Optional[np.ndarray] = None
    T: Optional[np.ndarray] = None
    C: Optional[np.ndarray] = None
    X: Optional[np.ndarray] = None
    pi: Optional[np.ndarray] = None
    infgen: Optional[np.ndarray] = None
    space: Optional[np.ndarray] = None
    space_aggr: Optional[np.ndarray] = None
    # per-(state,stateful,class) rates and hashed state space, exposed for the SolverENV state-vector analyzer reusing _compute_metrics_sync.
    arvRates: Optional[np.ndarray] = None
    depRates: Optional[np.ndarray] = None
    space_hashed: Optional[np.ndarray] = None
    sn: Optional[object] = None
    station_col_ranges: Optional[List[Tuple[int, int]]] = None
    rrobin_info: Optional[dict] = None
    eventFilt: Optional[List[np.ndarray]] = None
    runtime: float = 0.0
    method: str = "default"


def _refresh_phase_fields(sn, immediate_as_rate=False):
    """
    Pre-compute phasessz, phaseshift, mu, phi, pie from sn.proc.

    The afterEvent handlers require these fields but getStruct() doesn't
    populate them in the Python native implementation. This function
    extracts them from sn.proc (MAP/PH representation).

    immediate_as_rate: when True (SSA and CTMC paths), an immediate-service
    class (procid == IMMEDIATE, which has no (D0,D1) proc) is given a
    fast-exponential mu equal to its rate in sn.rates (~1e8), so the afterEvent
    DEP handler produces a non-zero departure rate. Without this, mu stays 0,
    the immediate departure never fires, and the job in that class is trapped:
    a Gillespie simulation deadlocks (empty result), and a CTMC accumulates the
    job with wrong QN/throughput. Mirrors MATLAB's mu~1e7 for immediate. Note
    the CTMC stochastic complementation folds only immediate pass-through NODES
    (Routers), not immediate station service, so it cannot substitute for this.
    The SSA caller (solver_ssa_run) snapshots/restores mu/phi/pie so this does
    not leak; the CTMC caller works on a private deepcopy of sn.
    """
    M = sn.nstations
    R = sn.nclasses

    if sn.phasessz is None:
        sn.phasessz = np.ones((M, R), dtype=int)
    if sn.phaseshift is None:
        sn.phaseshift = np.zeros((M, R), dtype=int)
    # Use nested dicts: mu[ist][r] = vec, matching MATLAB cell mu{ist,r}
    if sn.mu is None or not isinstance(sn.mu, dict):
        sn.mu = {}
    if sn.phi is None or not isinstance(sn.phi, dict):
        sn.phi = {}
    if sn.pie is None or not isinstance(sn.pie, dict):
        sn.pie = {}
    if sn.phases is None:
        sn.phases = np.ones((M, R), dtype=int)

    for ist in range(M):
        if ist not in sn.mu:
            sn.mu[ist] = {}
        if ist not in sn.phi:
            sn.phi[ist] = {}
        if ist not in sn.pie:
            sn.pie[ist] = {}

        for r in range(R):
            n_phases = 1
            D0 = None

            # marked (MMAP) source: proc holds the M3A cell [D0,D1_agg,D11..D1K]; never rewrite or misdetect it as [alpha,T].
            is_marked_ir = (getattr(sn, 'markidx', None) is not None
                            and ist < sn.markidx.shape[0]
                            and sn.markidx[ist, r] > 0)

            if sn.proc is not None and ist < len(sn.proc) and r < len(sn.proc[ist]):
                proc_ir = sn.proc[ist][r]
                if proc_ir is not None:
                    if is_marked_ir and isinstance(proc_ir, (list, tuple)) and len(proc_ir) >= 2:
                        D0 = np.atleast_2d(np.array(proc_ir[0], dtype=float))
                        D1 = np.atleast_2d(np.array(proc_ir[1], dtype=float))
                        n_phases = D0.shape[0]
                    elif isinstance(proc_ir, (list, tuple)) and len(proc_ir) >= 2:
                        first_elem = np.atleast_2d(np.array(proc_ir[0], dtype=float))
                        second_elem = np.atleast_2d(np.array(proc_ir[1], dtype=float))
                        # Check [alpha, T] format vs [D0, D1] format
                        if first_elem.shape[0] == 1 and second_elem.shape[0] == second_elem.shape[1]:
                            # [alpha, T] format: convert to D0/D1
                            alpha = first_elem.flatten()
                            T = second_elem
                            D0 = T
                            exit_rates = -np.sum(T, axis=1)
                            D1 = np.outer(exit_rates, alpha)
                        else:
                            D0 = first_elem
                            D1 = second_elem
                        n_phases = D0.shape[0]
                    elif isinstance(proc_ir, dict):
                        if 'k' in proc_ir and 'mu' in proc_ir:
                            # Erlang distribution
                            n_phases = int(proc_ir['k'])
                            mu_val = float(proc_ir['mu'])
                            if n_phases > 1:
                                D0 = np.zeros((n_phases, n_phases))
                                D1 = np.zeros((n_phases, n_phases))
                                for p in range(n_phases):
                                    D0[p, p] = -mu_val
                                    if p < n_phases - 1:
                                        D0[p, p + 1] = mu_val
                                D1[n_phases - 1, 0] = mu_val
                            else:
                                D0 = np.array([[-mu_val]])
                                D1 = np.array([[mu_val]])
                        elif 'probs' in proc_ir and 'rates' in proc_ir:
                            # HyperExponential distribution
                            probs = np.array(proc_ir['probs'])
                            rates = np.array(proc_ir['rates'])
                            n_phases = len(rates)
                            D0 = np.diag(-rates)
                            D1 = np.outer(rates, probs)
                        elif 'rate' in proc_ir:
                            # Exponential distribution
                            rate = proc_ir.get('rate', 1.0)
                            if rate is None or np.isnan(rate):
                                D0 = np.array([[np.nan]])
                                D1 = np.array([[np.nan]])
                            else:
                                D0 = np.array([[-rate]])
                                D1 = np.array([[rate]])
                        else:
                            n_phases = 1

            # converted (D0,D1) stored back into sn.proc for afterEvent handlers, except marked classes (M3A cell keeps per-mark matrices at 2..K+1).
            if D0 is not None and not is_marked_ir:
                proc_ir = sn.proc[ist][r] if (ist < len(sn.proc) and r < len(sn.proc[ist])) else None
                if isinstance(proc_ir, dict):
                    sn.proc[ist][r] = [D0, D1]
                elif isinstance(proc_ir, (list, tuple)) and len(proc_ir) >= 2:
                    # Ensure [alpha, T] is stored as [D0, D1]
                    sn.proc[ist][r] = [D0, D1]

            if is_marked_ir and sn.markidx[ist, r] > 1:
                # marked non-carrier class holds a single always-zero state column; the modulating chain lives in the carrier's phase block; mirrors MATLAB refreshProcessRepresentations.
                sn.phasessz[ist, r] = 1
            else:
                sn.phasessz[ist, r] = max(n_phases, 1)
            sn.phases[ist, r] = n_phases

            # Preserve authoritative mu/phi/pie already populated by refreshStruct,
            # per field. sn.proc[ist][r]'s second element may be the exit-rate
            # VECTOR (n x 1) rather than the full D1 matrix, and D0 sign/format is
            # not always the canonical generator form — recomputing the entry
            # distribution from that collapses multi-phase pie (e.g. HyperExp's
            # [0.5,0.5]) to [1.0]. Only (re)derive a field when it is missing.
            def _present(field):
                return (r in field[ist] and field[ist][r] is not None
                        and len(np.atleast_1d(field[ist][r])) == max(n_phases, 1))
            _has_pie = _present(sn.pie)
            _has_mu = _present(sn.mu)
            _has_phi = _present(sn.phi)

            if D0 is not None and not np.any(np.isnan(D0)):
                D1m = np.atleast_2d(D1)
                # exit-rate per phase: full D1 -> row sums; exit vector (n x 1) -> its entries
                if D1m.shape[1] == 1 and D0.shape[0] > 1:
                    exit_rate = D1m[:, 0]
                else:
                    exit_rate = np.sum(D1m, axis=1)
                mu_vec = np.abs(np.diag(D0))
                if not _has_mu:
                    sn.mu[ist][r] = mu_vec

                if not _has_phi:
                    phi_vec = np.zeros(n_phases)
                    for k in range(n_phases):
                        if mu_vec[k] != 0:
                            phi_vec[k] = abs(exit_rate[k]) / mu_vec[k]
                    sn.phi[ist][r] = phi_vec

                if not _has_pie:
                    # initial phase distribution needs a full D1 (column sums); an exit-vector-only representation falls back to entering phase 0.
                    if D1m.shape[1] == n_phases and n_phases > 1:
                        d1_col_sums = np.sum(D1m, axis=0)
                        total = np.sum(d1_col_sums)
                        pie_vec = d1_col_sums / total if total > 0 else np.eye(1, n_phases)[0]
                    else:
                        pie_vec = np.zeros(n_phases)
                        pie_vec[0] = 1.0
                    sn.pie[ist][r] = pie_vec
            else:
                imm_rate = None
                if immediate_as_rate and sn.procid is not None:
                    pid = np.asarray(sn.procid)
                    pid_ir = (pid[ist, r] if pid.ndim == 2 and ist < pid.shape[0]
                              and r < pid.shape[1] else None)
                    is_imm = (pid_ir == ProcessType.IMMEDIATE
                              or getattr(pid_ir, 'value', None) == ProcessType.IMMEDIATE.value)
                    if is_imm and sn.rates is not None:
                        rate = float(np.asarray(sn.rates)[ist, r])
                        imm_rate = rate if (np.isfinite(rate) and rate > 0) else 1e7
                if imm_rate is not None:
                    # Immediate service forced to a fast exponential (mu>0) so the afterEvent DEP handler fires; proc stays None since INF/PS departure uses mu directly.
                    sn.mu[ist][r] = np.array([imm_rate])
                    sn.phi[ist][r] = np.array([1.0])
                    sn.pie[ist][r] = np.array([1.0])
                else:
                    if not _has_mu:
                        sn.mu[ist][r] = np.array([0.0])
                    if not _has_phi:
                        sn.phi[ist][r] = np.array([1.0])
                    if not _has_pie:
                        sn.pie[ist][r] = np.array([1.0])

        # Compute phaseshift: cumulative sum of phases for earlier classes
        cum = 0
        for r in range(R):
            sn.phaseshift[ist, r] = cum
            cum += int(sn.phasessz[ist, r])


def _build_generator_sync(
    sn: NetworkStruct,
    state_space: np.ndarray,
    state_space_hashed: np.ndarray,
    options: SolverCTMCOptions,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List]:
    """
    Build the infinitesimal generator matrix using sync-action-based architecture.

    This uses the afterEvent dispatch mechanism (matching MATLAB/JAR). Each sync
    action pairs an active event (DEP/PHASE/READ) with a passive event
    (ARV/LOCAL) and iterates over all states.

    Port from MATLAB solver_ctmc.m:52-236 / JAR Solver_ctmc.kt:64-290.

    Args:
        sn: Network structure
        state_space: Full state space matrix (one row per global state)
        state_space_hashed: Hashed state space (one col per stateful node, values = row indices in sn.space)
        options: Solver options

    Returns:
        Tuple of (Q, arvRates, depRates, Dfilt):
        - Q: Infinitesimal generator matrix
        - arvRates: Arrival rates per (state, stateful_node, class)
        - depRates: Departure rates per (state, stateful_node, class)
        - Dfilt: List of per-action rate matrices
    """
    from ...state.after_event import after_event_hashed, build_space_hash
    from ...state.ctmc_ssg import _is_bas_station
    from ....lang.sync import refresh_sync
    from ...mc import ctmc_makeinfgen
    from ....constants import EventType

    nstateful = sn.nstateful
    nclasses = sn.nclasses
    nnodes = sn.nnodes
    local = nnodes  # Sentinel for "no passive node"

    # Build sync actions
    sync = refresh_sync(sn)
    A = len(sync)

    # Build hash maps for O(1) state lookup
    hash_maps = build_space_hash(sn)

    n_states = state_space_hashed.shape[0]

    # Build global state hash for matchrow
    global_hash = {}
    for s in range(n_states):
        key = tuple(state_space_hashed[s].astype(int))
        global_hash[key] = s

    # Initialize Q and Dfilt
    Q = np.eye(n_states)  # Will be corrected later
    Dfilt = [np.zeros((n_states, n_states)) for _ in range(A)]
    # true-BAS become-blocked transitions change the chain but are not station departures; kept out of Dfilt so throughput is not double-counted; see _kb/06-solver-catalog.md True BAS blocking section.
    basBlockQ = np.zeros((n_states, n_states))

    arvRates = np.zeros((n_states, nstateful, nclasses))
    depRates = np.zeros((n_states, nstateful, nclasses))

    # immediate-provenance accumulator Qimm; see _kb/06-solver-catalog.md Vanishing states: the row purge before stochastic complementation.
    Qimm = np.zeros((n_states, n_states))
    immAction = [False] * A

    # Check class-switching mask
    csmask = sn.csmask if hasattr(sn, 'csmask') and sn.csmask is not None else np.ones((nclasses, nclasses))

    for a in range(A):
        act = sync[a]
        node_a = act.active.node
        class_a = act.active.job_class
        event_a = act.active.event

        node_p = act.passive.node
        class_p = act.passive.job_class
        event_p = act.passive.event

        if not sn.isstateful[node_a]:
            continue
        isf_a = int(sn.nodeToStateful[node_a])

        # Immediate-scale actions are Router/Join pass-through, matching the vanishing states marked by _find_immediate_states_sync.
        immAction[a] = _is_immediate_passthrough_node(sn, int(node_a))

        is_local = (node_p >= nnodes)  # Passive is LOCAL sentinel
        isf_p = -1
        if not is_local and sn.isstateful[node_p]:
            isf_p = int(sn.nodeToStateful[node_p])

        for s in range(n_states):
            state = state_space_hashed[s].copy().astype(int)
            state_a = state[isf_a]

            # Fire active event
            new_state_a, rate_a, outprob_a = after_event_hashed(
                sn, node_a, state_a, event_a, class_a, hash_maps)

            if np.all(new_state_a < 0):
                continue

            for ia in range(len(new_state_a)):
                if new_state_a[ia] < 0 or rate_a[ia] == 0:
                    continue

                if is_local:
                    # Local action: only active node changes
                    new_state = state.copy()
                    new_state[isf_a] = new_state_a[ia]
                    ns_key = tuple(new_state)
                    ns = global_hash.get(ns_key, -1)
                    if ns >= 0:
                        prob_sync = 1.0
                        Dfilt[a][s, ns] += rate_a[ia] * prob_sync

                else:
                    # Non-local: fire passive event
                    if isf_p < 0:
                        continue

                    if node_p == node_a:
                        # Same node: passive sees post-active state
                        state_p = int(new_state_a[ia])
                    else:
                        state_p = state[isf_p]

                    new_state_p, rate_p, outprob_p = after_event_hashed(
                        sn, node_p, state_p, event_p, class_p, hash_maps)

                    for ip in range(len(new_state_p)):
                        if new_state_p[ip] < 0:
                            continue

                        # Compute sync probability
                        if callable(act.passive.prob):
                            # State-dependent routing
                            # Build state cell arrays for rtfun evaluation
                            state_before = [None] * nstateful
                            state_after = [None] * nstateful
                            for isf in range(nstateful):
                                state_before[isf] = sn.space[isf][state[isf]:state[isf] + 1] if isf in sn.space else np.array([[]])
                                state_after[isf] = state_before[isf].copy()
                            state_after[isf_a] = sn.space[isf_a][new_state_a[ia]:new_state_a[ia] + 1]
                            if isf_p >= 0:
                                state_after[isf_p] = sn.space[isf_p][new_state_p[ip]:new_state_p[ip] + 1]
                            prob_sync = act.passive.prob(state_before, state_after)
                            if isinstance(prob_sync, np.ndarray):
                                prob_sync = float(prob_sync.ravel()[0])
                        else:
                            prob_sync = float(act.passive.prob) if act.passive.prob is not None else 1.0

                        prob_sync *= outprob_p[ip] if ip < len(outprob_p) else 1.0

                        # Build new global state
                        new_state = state.copy()
                        new_state[isf_a] = new_state_a[ia]
                        new_state[isf_p] = new_state_p[ip]

                        ns_key = tuple(new_state)
                        ns = global_hash.get(ns_key, -1)
                        if ns >= 0:
                            Dfilt[a][s, ns] += rate_a[ia] * prob_sync

                    # true-BAS become-blocked arc: hold the completed job at the server (marker 0->1) at rate mu; recorded in basBlockQ only, not Dfilt; see _kb/06-solver-catalog.md True BAS blocking section.
                    if (event_a == EventType.DEP
                            and _is_bas_station(sn, int(sn.nodeToStation[node_a]))
                            and all(new_state_p[ip2] < 0 for ip2 in range(len(new_state_p)))):
                        curVec = sn.space[isf_a][state[isf_a]]
                        if curVec[-1] == 0:
                            blockedVec = curVec.copy()
                            blockedVec[-1] = 1
                            matches = np.where((sn.space[isf_a] == blockedVec).all(axis=1))[0]
                            if len(matches) > 0:
                                nsb_state = state.copy()
                                nsb_state[isf_a] = matches[0]
                                nsb = global_hash.get(tuple(nsb_state), -1)
                                if nsb >= 0:
                                    basBlockQ[s, nsb] += rate_a[ia]

    # Dfilt summed into Q; DEP/ARV accounting deferred until after the vanishing-row purge below.
    for a in range(A):
        Q = Q + Dfilt[a]
        if immAction[a]:
            Qimm = Qimm + Dfilt[a]
    # Fold in true-BAS become-blocked transitions (not counted as departures).
    Q = Q + basBlockQ

    # SPN gsync (ENABLE/FIRE) events; see _kb/06-solver-catalog.md Vanishing states section. Mirrors jar Solver_ctmc.kt:298-395 / matlab solver_ctmc.m:242-325.
    from ....lang.sync import refresh_global_sync as _refresh_gsync
    from ...state.after_global_event import after_global_event as _after_gsync
    if not hasattr(sn, 'gsync') or sn.gsync is None or len(sn.gsync) == 0:
        sn.gsync = _refresh_gsync(sn)
    G = len(sn.gsync) if sn.gsync is not None else 0
    DfiltGsyncComp = [np.zeros((n_states, n_states)) for _ in range(G)]
    immGsync = [False] * G

    for g in range(G):
        glevent = sn.gsync[g]
        if not glevent.active:
            continue
        gind = int(glevent.active[0].node)
        isf_transition = int(sn.nodeToStateful[gind])
        # ENABLE phase moves and firings of a TimingStrategy.IMMEDIATE mode are
        # the two gsync sources emitted at the Immediate scale.
        immGsync[g] = (glevent.active[0].event == EventType.ENABLE) or (
            glevent.active[0].event == EventType.FIRE
            and _is_immediate_mode(sn, gind, int(glevent.active[0].mode)))

        for s in range(n_states):
            state = state_space_hashed[s].astype(int)
            glspace = []
            for isf in range(nstateful):
                state_idx = int(state[isf])
                space_isf = np.atleast_2d(sn.space[isf]) if isf in sn.space else np.zeros((1, 0))
                if state_idx < 0 or state_idx >= space_isf.shape[0]:
                    glspace.append(np.zeros(space_isf.shape[1]))
                else:
                    glspace.append(space_isf[state_idx])

            result = _after_gsync(sn, gind, glspace, glevent, False)
            if result.outrate.size == 0:
                continue
            if not np.any(result.outrate > 0):
                continue

            for io in range(result.outrate.size):
                rate_io = float(result.outrate[io])
                if rate_io <= 0:
                    continue
                prob_io = float(result.outprob[io]) if io < result.outprob.size else 1.0
                eff = rate_io * prob_io
                if eff == 0:
                    continue

                # Build the new hashed state by re-hashing every node row
                # whose contents changed.
                new_state = state.copy()
                gl_io = result.outglspace[io]
                ok = True
                for isf in range(nstateful):
                    if isf not in sn.space or sn.space[isf] is None:
                        continue
                    space_isf = np.atleast_2d(sn.space[isf])
                    new_row = np.asarray(gl_io[isf], dtype=float).ravel()
                    if np.array_equal(new_row, glspace[isf]):
                        continue
                    matches = np.where(np.all(np.isclose(space_isf, new_row), axis=1))[0]
                    if matches.size == 0:
                        ok = False
                        break
                    new_state[isf] = int(matches[0])
                if not ok:
                    continue

                ns_key = tuple(int(x) for x in new_state)
                ns = global_hash.get(ns_key, -1)
                if ns < 0:
                    continue

                Q[s, ns] += eff
                if immGsync[g]:
                    Qimm[s, ns] += eff
                if bool(result.is_completion[io]):
                    DfiltGsyncComp[g][s, ns] += eff

    # fork firing synchronizations: each firing consumes the Fork's held parent and emits one sibling per branch atomically; mirrors matlab solver_ctmc.m:324-374.
    from ...state.after_fj_event import after_fj_event as _after_fj
    fjsync = sn.fjsync if getattr(sn, 'fjsync', None) else []
    FJ = len(fjsync)
    Dfilt_fjsync = [np.zeros((n_states, n_states)) for _ in range(FJ)]
    for k in range(FJ):
        for s in range(n_states):
            state = state_space_hashed[s].astype(int)
            glspace = []
            for isf in range(nstateful):
                space_isf = np.atleast_2d(sn.space[isf]) if isf in sn.space else np.zeros((1, 0))
                idx = int(state[isf])
                glspace.append(space_isf[idx] if 0 <= idx < space_isf.shape[0] else np.zeros(space_isf.shape[1]))
            fj_states, fj_rate, fj_prob = _after_fj(sn, fjsync[k], glspace, False)
            for io in range(len(fj_states)):
                if io < len(fj_prob) and fj_prob[io] <= 0:
                    continue
                new_state = state.copy()
                ok = True
                gl_io = fj_states[io]
                for isf in range(nstateful):
                    space_isf = np.atleast_2d(sn.space[isf]) if isf in sn.space else np.zeros((1, 0))
                    new_row = np.asarray(gl_io[isf], dtype=float).ravel()
                    if np.array_equal(new_row, np.asarray(glspace[isf], dtype=float).ravel()):
                        continue
                    matches = np.where(np.all(np.isclose(space_isf, new_row), axis=1))[0]
                    if matches.size == 0:
                        ok = False
                        break
                    new_state[isf] = int(matches[0])
                if not ok:
                    continue
                ns = global_hash.get(tuple(int(x) for x in new_state), -1)
                if ns < 0:
                    continue
                eff = float(fj_rate[io]) * (float(fj_prob[io]) if io < len(fj_prob) else 1.0)
                Q[s, ns] += eff
                Qimm[s, ns] += eff
                Dfilt_fjsync[k][s, ns] += eff

    # vanishing-row purge; see _kb/06-solver-catalog.md Vanishing states: the row purge before stochastic complementation.
    _, imm_rows = _find_immediate_states_sync(sn, state_space_hashed)
    if imm_rows:
        _immidx = np.asarray(imm_rows, dtype=int)
        # a vanishing state with an empty Qimm row means the predicate and provenance tagging drifted apart; left unpurged (with a reported gap) rather than manufacturing an absorbing state.
        _gap = _immidx[np.asarray(Qimm[_immidx, :].sum(axis=1)).ravel() <= 0]
        if _gap.size:
            from ...io.logging import line_warning_always
            line_warning_always(
                'solver_ctmc',
                "CTMC: %d vanishing state(s) have no immediate outgoing arc; "
                "the vanishing predicate and the immediate-arc tagging disagree, "
                "so those rows keep their timed arcs." % _gap.size)
            _immidx = np.setdiff1d(_immidx, _gap)
        Q[_immidx, :] = Qimm[_immidx, :]
        for a in range(A):
            if not immAction[a]:
                Dfilt[a][_immidx, :] = 0.0
        for g in range(G):
            if not immGsync[g]:
                DfiltGsyncComp[g][_immidx, :] = 0.0

    # Deferred DEP/ARV accounting, computed from the purged filters.
    for a in range(A):
        act = sync[a]
        node_a = act.active.node
        if not sn.isstateful[node_a]:
            continue
        isf_a = int(sn.nodeToStateful[node_a])
        node_p = act.passive.node
        is_local = (node_p >= nnodes)
        isf_p = -1
        if not is_local and sn.isstateful[node_p]:
            isf_p = int(sn.nodeToStateful[node_p])
        row_sums = np.sum(Dfilt[a], axis=1)
        if act.active.event == EventType.DEP:
            depRates[:, isf_a, act.active.job_class] += row_sums
        if not is_local and act.passive.event == EventType.ARV and isf_p >= 0:
            arvRates[:, isf_p, act.passive.job_class] += row_sums

    for g in range(G):
        glevent = sn.gsync[g]
        if not glevent.active or glevent.active[0].event != EventType.FIRE:
            continue
        row_sums = np.sum(DfiltGsyncComp[g], axis=1)
        for pev in glevent.passive:
            pev_node = int(pev.node)
            if pev_node >= sn.nnodes or not sn.isstateful[pev_node]:
                continue
            pev_isf = int(sn.nodeToStateful[pev_node])
            pev_class = int(pev.job_class)
            if pev_class < 0 or pev_class >= nclasses:
                continue
            if pev.event == EventType.PRE:
                depRates[:, pev_isf, pev_class] += row_sums
            elif pev.event == EventType.POST:
                arvRates[:, pev_isf, pev_class] += row_sums

    # Remove initial identity diagonal
    Q = Q - np.diag(np.diag(Q))

    # Make proper infinitesimal generator (row sums = 0)
    Q = ctmc_makeinfgen(Q)

    # DfiltGsyncComp lets the rate complement recover firings that occur only in vanishing markings.
    return Q, arvRates, depRates, Dfilt, Dfilt_fjsync, sync, DfiltGsyncComp


def _is_immediate_mode(sn, ind, mode):
    """True when mode of Transition node ind is declared TimingStrategy.IMMEDIATE.

    timingstrategies holds a TimingStrategy enum once set, but defaults to the
    string 'TIMED', so the comparison is by name.
    """
    if sn.nodeparam is None or ind not in sn.nodeparam or sn.nodeparam[ind] is None:
        return False
    timing = getattr(sn.nodeparam[ind], 'timingstrategies', None)
    if timing is None or mode >= len(timing):
        return False
    return str(getattr(timing[mode], 'name', timing[mode])).upper() == 'IMMEDIATE'


def _is_immediate_passthrough_node(sn, ind):
    """True when node `ind` emits its outgoing arcs at the GlobalConstants.
    Immediate scale, i.e. it is one of the node types whose occupied states
    _find_immediate_states_sync marks as vanishing. Keep the two in step: a node
    marked vanishing but not recognised here would have its whole row zeroed.
    """
    nt = sn.nodetype[ind]
    nt_val = int(nt.value) if hasattr(nt, 'value') else int(nt)
    router_val = int(NodeType.ROUTER.value) if hasattr(NodeType.ROUTER, 'value') else int(NodeType.ROUTER)
    fork_val = int(NodeType.FORK.value) if hasattr(NodeType.FORK, 'value') else int(NodeType.FORK)
    join_val = int(NodeType.JOIN.value) if hasattr(NodeType.JOIN, 'value') else int(NodeType.JOIN)
    if nt_val == router_val:
        return True
    if getattr(sn, 'isfjaugmented', False) and nt_val in (fork_val, join_val):
        return True
    return False


def _find_immediate_states_sync(sn, state_space_hashed):
    """
    Find immediate states for stochastic complementation (sync-based path).

    Port from MATLAB solver_ctmc.m lines 286-298.
    Immediate states are those where a non-station, non-Cache stateful node has jobs.

    Args:
        sn: NetworkStruct with sn.space populated
        state_space_hashed: Hashed state space (n_states x nstateful)

    Returns:
        Tuple of (nonimm_indices, imm_indices) - lists of global state indices
    """
    nclasses = sn.nclasses
    n_states = state_space_hashed.shape[0]

    # Design Y: positive list of immediate pass-through node types; see _kb/06-solver-catalog.md Design Y section.
    router_val = int(NodeType.ROUTER.value) if hasattr(NodeType.ROUTER, 'value') else int(NodeType.ROUTER)
    fork_val = int(NodeType.FORK.value) if hasattr(NodeType.FORK, 'value') else int(NodeType.FORK)
    join_val = int(NodeType.JOIN.value) if hasattr(NodeType.JOIN, 'value') else int(NodeType.JOIN)
    _isFJ = getattr(sn, 'isfjaugmented', False)

    imm_set = set()

    for ind in range(sn.nnodes):
        if not sn.isstateful[ind] or sn.isstation[ind]:
            continue
        nt_val = int(sn.nodetype[ind].value) if hasattr(sn.nodetype[ind], 'value') else int(sn.nodetype[ind])
        # a Fork holding the FJ parent job is immediate (one vanishing state before the fork firing, sn.fjsync).
        is_immediate_pass_through = (nt_val == router_val) or (_isFJ and nt_val == fork_val)
        if not is_immediate_pass_through:
            continue

        isf = int(sn.nodeToStateful[ind])
        if isf not in sn.space or sn.space[isf] is None:
            continue

        space_isf = np.atleast_2d(sn.space[isf])
        # Find per-node states where node has jobs (sum of first nclasses cols > 0)
        n_cols = min(nclasses, space_isf.shape[1])
        imm_st = set()
        for row_idx in range(space_isf.shape[0]):
            if np.sum(space_isf[row_idx, :n_cols]) > 0:
                imm_st.add(row_idx)

        if not imm_st:
            continue

        # Find global states where this node's hash is in imm_st
        for s in range(n_states):
            h = int(state_space_hashed[s, isf])
            if h in imm_st:
                imm_set.add(s)

    # a Join holding a complete sibling set fires immediately (vanishing); incomplete sibling sets are genuine synchronization-delay states. Mirrors solver_ctmc.m:493-518.
    if _isFJ:
        from ...state.after_event_join import after_event_join as _aej
        from ....constants import EventType
        for ind in range(sn.nnodes):
            nt_val = int(sn.nodetype[ind].value) if hasattr(sn.nodetype[ind], 'value') else int(sn.nodetype[ind])
            if nt_val != join_val:
                continue
            isf = int(sn.nodeToStateful[ind])
            if isf not in sn.space or sn.space[isf] is None:
                continue
            fjp = sn.nodeparam[ind].get('fj') if (sn.nodeparam is not None and ind in sn.nodeparam
                                                  and isinstance(sn.nodeparam[ind], dict)) else None
            if fjp is None:
                continue
            origcl = np.asarray(fjp['origclasses']).ravel()
            space_isf = np.atleast_2d(sn.space[isf])
            firable_rows = set()
            for row_idx in range(space_isf.shape[0]):
                for r in origcl:
                    ospace, _, _ = _aej(sn, ind, space_isf[row_idx:row_idx + 1], EventType.DEP, int(r), False)
                    if ospace is not None and np.atleast_2d(ospace).size > 0:
                        firable_rows.add(row_idx)
                        break
            if firable_rows:
                for s in range(n_states):
                    if int(state_space_hashed[s, isf]) in firable_rows:
                        imm_set.add(s)

    # a Transition state is immediate whenever an ENABLE event would change its own row. Mirrors jar Solver_ctmc.kt:556-595.
    if hasattr(sn, 'gsync') and sn.gsync:
        from ...state.after_global_event import after_global_event as _after_gsync
        from ....constants import EventType
        for s in range(n_states):
            if s in imm_set:
                continue
            state = state_space_hashed[s].astype(int)
            glspace = []
            for isf in range(int(sn.nstateful)):
                if isf in sn.space and sn.space[isf] is not None:
                    space_isf = np.atleast_2d(sn.space[isf])
                    idx = int(state[isf])
                    glspace.append(space_isf[idx] if 0 <= idx < space_isf.shape[0] else np.zeros(space_isf.shape[1]))
                else:
                    glspace.append(np.zeros(0))
            for g_idx, glevent in sn.gsync.items():
                if not glevent.active:
                    continue
                if glevent.active[0].event != EventType.ENABLE:
                    continue
                gind = int(glevent.active[0].node)
                isf_t = int(sn.nodeToStateful[gind])
                orig_row = np.asarray(glspace[isf_t], dtype=float).ravel().copy()
                result = _after_gsync(sn, gind, glspace, glevent, False)
                if result.outrate.size == 0:
                    continue
                if not np.any(result.outrate > 0):
                    continue
                changed = False
                for io in range(result.outrate.size):
                    if float(result.outrate[io]) <= 0:
                        continue
                    new_row = np.asarray(result.outglspace[io][isf_t], dtype=float).ravel()
                    if not np.array_equal(new_row, orig_row):
                        changed = True
                        break
                if changed:
                    imm_set.add(s)
                    break

        # SPN IMMEDIATE-timed markings are vanishing; see _kb/06-solver-catalog.md Vanishing states section.
        for s in range(n_states):
            if s in imm_set:
                continue
            state = state_space_hashed[s].astype(int)
            glspace = []
            for isf in range(int(sn.nstateful)):
                if isf in sn.space and sn.space[isf] is not None:
                    space_isf = np.atleast_2d(sn.space[isf])
                    idx = int(state[isf])
                    glspace.append(space_isf[idx] if 0 <= idx < space_isf.shape[0] else np.zeros(space_isf.shape[1]))
                else:
                    glspace.append(np.zeros(0))
            for g_idx, glevent in sn.gsync.items():
                if not glevent.active:
                    continue
                if glevent.active[0].event != EventType.FIRE:
                    continue
                gind = int(glevent.active[0].node)
                if not _is_immediate_mode(sn, gind, int(glevent.active[0].mode)):
                    continue
                result = _after_gsync(sn, gind, glspace, glevent, False)
                if result.outrate.size and np.any(result.outrate > 0):
                    imm_set.add(s)
                    break

    imm_indices = sorted(imm_set)
    nonimm_indices = sorted(set(range(n_states)) - imm_set)
    return nonimm_indices, imm_indices



def _sn_all_phasetype(sn) -> bool:
    """True when every station-class process admits a phase-type reading.

    False as soon as one service or arrival process is a matrix exponential or a
    rational arrival process, in which case the stationary vector of the
    generator is a signed measure. See sn_is_phasetype and _kb/04-networkstruct.md.
    """
    isph = getattr(sn, 'isph', None)
    if isph is None:
        return True
    return bool(np.all(np.asarray(isph, dtype=bool)))


def _compute_metrics_sync(sn, pi, arvRates, depRates, state_space_aggr,
                          state_space, state_space_hashed, options=None):
    """
    Compute performance metrics from CTMC steady-state distribution (sync-based path).

    Port from MATLAB solver_ctmc_analyzer.m.

    Args:
        sn: NetworkStruct
        pi: Steady-state distribution (1D array, length = n_states)
        arvRates: Arrival rates (n_states x nstateful x nclasses)
        depRates: Departure rates (n_states x nstateful x nclasses)
        state_space_aggr: Aggregated state space (n_states x M*R)
        state_space: Full state space
        state_space_hashed: Hashed state space
            the optional trailing options argument of solver_ctmc_avg_from_pi.m)

    Returns:
        Dict with keys 'Q', 'U', 'R', 'T'
    """
    M = sn.nstations
    K = sn.nclasses
    n_states = len(pi)

    QN = np.zeros((M, K))
    UN = np.zeros((M, K))
    RN = np.zeros((M, K))
    TN = np.zeros((M, K))

    # Ensure pi is 1D
    pi = pi.ravel()
    if len(pi) != n_states:
        return {'Q': QN, 'U': UN, 'R': RN, 'T': TN}

    # Clean up pi. The clamp removes the tiny negative residues a genuine CTMC
    # solve leaves behind, but a station whose service is a matrix exponential
    # makes the stationary vector a genuinely SIGNED measure: the generator has
    # negative off-diagonal entries by construction, and only the aggregates
    # over each phase block are probabilities. Clamping there deletes real mass
    # (an M/CME/1 came out with utilization 3.19), so the sign is kept and every
    # metric below, being linear in pi, stays exact.
    signed = not _sn_all_phasetype(sn)
    if not signed:
        pi[pi < 1e-14] = 0
    pi_sum = np.sum(pi)
    if pi_sum > 0:
        pi = pi / pi_sum

    # System throughput: XN(k) = pi * arvRates[:, refsf, k]
    XN = np.zeros(K)
    for k in range(K):
        ref_stat = int(sn.refstat[k]) if hasattr(sn, 'refstat') and k < len(sn.refstat) else 0
        if ref_stat < M:
            refsf = int(sn.stationToStateful[ref_stat]) if hasattr(sn, 'stationToStateful') else ref_stat
            if refsf < arvRates.shape[1]:
                XN[k] = np.dot(pi, arvRates[:, refsf, k])

    # Per-station metrics
    for ist in range(M):
        isf = int(sn.stationToStateful[ist]) if hasattr(sn, 'stationToStateful') else ist
        ind = int(sn.stationToNode[ist]) if hasattr(sn, 'stationToNode') else ist

        # Check if Source
        nt_val = int(sn.nodetype[ind].value) if hasattr(sn.nodetype[ind], 'value') else int(sn.nodetype[ind])
        source_val = int(NodeType.SOURCE.value) if hasattr(NodeType.SOURCE, 'value') else int(NodeType.SOURCE)
        is_source = (nt_val == source_val)

        for k in range(K):
            # Throughput: TN = pi * depRates
            if isf < depRates.shape[1]:
                TN[ist, k] = np.dot(pi, depRates[:, isf, k])

            # Queue length: QN = pi * stateSpaceAggr
            if not is_source:
                col = ist * K + k
                if col < state_space_aggr.shape[1]:
                    QN[ist, k] = np.dot(pi, state_space_aggr[:, col])

        if is_source:
            continue

        # Utilization
        S = sn.nservers[ist] if hasattr(sn, 'nservers') else 1
        sched = sn.sched[ist] if hasattr(sn, 'sched') else SchedStrategy.FCFS

        if sched == SchedStrategy.PAS:
            # PAS utilization = time-average in-service jobs per class / servers; sir already counts positive-rate-increment positions so multi-server-type jobs count once.
            from ...state.marginal import toMarginal
            space_isf = (np.atleast_2d(sn.space[isf])
                         if (sn.space is not None and isf in sn.space and sn.space[isf] is not None)
                         else None)
            if space_isf is not None:
                for s_idx in range(n_states):
                    h = int(state_space_hashed[s_idx, isf])
                    _, _, sir, _ = toMarginal(sn, ind, space_isf[h:h + 1])
                    sir = np.atleast_1d(sir).flatten()
                    for k in range(K):
                        if k < len(sir):
                            UN[ist, k] += pi[s_idx] * sir[k] / S
            continue

        # load/class/joint-dependent stations use the PS-share formula UN[k]=E[(n_k*w_k)/sum_j(n_j*w_j)], not the plain rho-style T*mean/S.
        _has_sd = False
        if hasattr(sn, 'cdscaling') and sn.cdscaling is not None:
            _cd = sn.cdscaling
            _cdf = (_cd.get(ist) if isinstance(_cd, dict) else (_cd[ist] if ist < len(_cd) else None))
            if callable(_cdf):
                _has_sd = True
        if not _has_sd and getattr(sn, 'jdscaling', None) is not None:
            _jd = sn.jdscaling
            _jdf = (_jd.get(ist) if isinstance(_jd, dict) else (_jd[ist] if ist < len(_jd) else None))
            if callable(_jdf):
                _has_sd = True
        if not _has_sd and hasattr(sn, 'lldscaling') and sn.lldscaling is not None:
            _lld = np.asarray(sn.lldscaling)
            if _lld.ndim >= 2 and ist < _lld.shape[0] and not np.allclose(_lld[ist, :], 1.0):
                _has_sd = True

        _ps_like = sched in (SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS)

        if _has_sd and _ps_like:
            # busy-server fraction under load-dependent scaling weights the PS share by the current scaling and normalizes by the peak (NC convention).
            weights = np.ones(K)
            if hasattr(sn, 'schedparam') and sn.schedparam is not None:
                for r in range(K):
                    sp = sn.schedparam[ist, r]
                    if sp is not None and float(sp) > 0:
                        weights[r] = float(sp)
            _lldrow = None
            ceff = float(S)
            if hasattr(sn, 'lldscaling') and sn.lldscaling is not None:
                _lld = np.asarray(sn.lldscaling)
                if _lld.ndim >= 2 and ist < _lld.shape[0]:
                    _lldrow = _lld[ist, :]
                    ceff = max(ceff, float(np.max(_lldrow)))
            for s_idx in range(n_states):
                nrow = np.array([state_space_aggr[s_idx, ist * K + j]
                                 if ist * K + j < state_space_aggr.shape[1] else 0.0
                                 for j in range(K)])
                wtot = float(np.sum(nrow * weights))
                if wtot > 0:
                    lldnow = 1.0
                    if _lldrow is not None:
                        lldnow = float(_lldrow[min(max(int(np.sum(nrow)), 1) - 1, len(_lldrow) - 1)])
                    for k in range(K):
                        UN[ist, k] += pi[s_idx] * nrow[k] * weights[k] / wtot * lldnow / ceff
        else:
            # effective server count for a load-level-dependent-scaling station (nservers stays 1) is the max scaling, matching the flat analyzer.
            S_eff = S
            if hasattr(sn, 'lldscaling') and sn.lldscaling is not None:
                _lld = np.asarray(sn.lldscaling)
                if _lld.ndim >= 2 and ist < _lld.shape[0]:
                    _mx = float(np.max(_lld[ist, :]))
                    if _mx > S_eff:
                        S_eff = _mx
            # class-dependent scaling's effective server count is the cap over reachable states.
            _sdf = None
            if hasattr(sn, 'cdscaling') and sn.cdscaling is not None:
                _cd = sn.cdscaling
                _sdf = (_cd.get(ist) if isinstance(_cd, dict) else (_cd[ist] if ist < len(_cd) else None))
            if not callable(_sdf) and getattr(sn, 'jdscaling', None) is not None:
                _jd = sn.jdscaling
                _sdf = (_jd.get(ist) if isinstance(_jd, dict) else (_jd[ist] if ist < len(_jd) else None))
            if callable(_sdf):
                _cdmax = 0.0
                for s_idx in range(n_states):
                    nrow = np.array([state_space_aggr[s_idx, ist * K + j]
                                     if ist * K + j < state_space_aggr.shape[1] else 0.0
                                     for j in range(K)])
                    try:
                        _val = np.atleast_1d(np.asarray(_sdf(nrow), dtype=float))
                        _cdmax = max(_cdmax, float(np.max(_val)))
                    except Exception:
                        pass
                if _cdmax > S_eff:
                    S_eff = _cdmax
            for k in range(K):
                if sched == SchedStrategy.INF:
                    UN[ist, k] = QN[ist, k]
                else:
                    # MAP/MMPP service mean comes from (D0,D1) since sn.rates is not its reciprocal; matches the flat analyzer's UN=TN*map_mean/c.
                    mean_s = 0.0
                    is_map_k = False
                    if hasattr(sn, 'procid') and sn.procid is not None:
                        try:
                            _pid = sn.procid[ist, k]
                            _pv = int(_pid.value) if hasattr(_pid, 'value') else int(_pid)
                            is_map_k = _pv in (
                                int(ProcessType.MAP.value) if hasattr(ProcessType.MAP, 'value') else int(ProcessType.MAP),
                                int(ProcessType.MMPP2.value) if hasattr(ProcessType.MMPP2, 'value') else int(ProcessType.MMPP2))
                        except (TypeError, KeyError, IndexError):
                            is_map_k = False
                    if is_map_k:
                        try:
                            _pr = sn.proc[ist][k]
                            _D0 = np.atleast_2d(np.asarray(_pr[0], dtype=float))
                            _D1 = np.atleast_2d(np.asarray(_pr[1], dtype=float))
                            mean_s = map_mean(_D0, _D1)
                        except (TypeError, KeyError, IndexError, AttributeError):
                            mean_s = 0.0
                    elif hasattr(sn, 'rates') and sn.rates is not None:
                        rate = sn.rates[ist, k]
                        if rate > 0 and not np.isnan(rate) and not np.isinf(rate):
                            mean_s = 1.0 / rate

                    if mean_s > 0:
                        UN[ist, k] = TN[ist, k] * mean_s / S_eff

        # T*E[S]/c is exact only for exponential service under a destructive signal; busy-server occupancy is read off the state space instead (the lld/cd branch already accumulates it; INF reports queue length).
        if not _has_sd and sched != SchedStrategy.INF:
            _arv_at_station = np.zeros(K)
            if isf < arvRates.shape[1]:
                for r in range(K):
                    _arv_at_station[r] = float(np.dot(pi, arvRates[:, isf, r]))
            _lossy = signal_lossy_classes(sn, _arv_at_station)
            if _lossy.any():
                _space_isf = (np.atleast_2d(sn.space[isf])
                              if (sn.space is not None and isf in sn.space
                                  and sn.space[isf] is not None)
                              else None)
                _UNb = busy_fraction(sn, ind, ist, sched, S, _space_isf,
                                     state_space_hashed, isf, pi, K)
                for k in range(K):
                    if _lossy[k]:
                        UN[ist, k] = _UNb[k]

    # limited class-dependent stations report Util=T*S/peak (sn.cdscalingpeak), matching ordinary multiserver convention; mirrors MATLAB solver_ctmc_analyzer.m/solver_ctmc_avg_from_pi.m.
    _cd = getattr(sn, 'cdscaling', None)
    _njobs = np.asarray(sn.njobs, dtype=float).flatten() if sn.njobs is not None else np.array([])
    if _cd is not None and len(_cd) > 0 and _njobs.size > 0 and not np.any(np.isinf(_njobs)):
        _peak = sn.cdscalingpeak
        for ist in range(M):
            _cdf = (_cd.get(ist) if isinstance(_cd, dict)
                    else (_cd[ist] if ist < len(_cd) else None))
            if _cdf is None or not callable(_cdf):
                continue
            for k in range(K):
                bmax = _peak[ist, k]
                rate = sn.rates[ist, k]
                if np.isfinite(rate) and rate > 0 and bmax > 0:
                    UN[ist, k] = TN[ist, k] / rate / bmax
                else:
                    UN[ist, k] = 0.0

    # Joint-dependent stations report Util=T*S/peak using sn.jdscalingpeak,
    # mirroring the class-dependence block above (non-product-form eta_i).
    _jd = getattr(sn, 'jdscaling', None)
    if _jd is not None and len(_jd) > 0 and _njobs.size > 0 and not np.any(np.isinf(_njobs)):
        _peak = sn.jdscalingpeak
        for ist in range(M):
            _jdf = (_jd.get(ist) if isinstance(_jd, dict)
                    else (_jd[ist] if ist < len(_jd) else None))
            if _jdf is None or not callable(_jdf):
                continue
            for k in range(K):
                bmax = _peak[ist, k]
                rate = sn.rates[ist, k]
                if np.isfinite(rate) and rate > 0 and bmax > 0:
                    UN[ist, k] = TN[ist, k] / rate / bmax
                else:
                    UN[ist, k] = 0.0

    # BUG-82 BAS destination-convention QLen attribution (reporting only); see _kb/06-solver-catalog.md True BAS blocking and _kb/04-networkstruct.md isbasblocking field.
    isb = getattr(sn, 'isbasblocking', None)
    if (isb is not None and np.asarray(isb).size > 0
            and getattr(sn, 'connmatrix', None) is not None
            and sn.space is not None):
        isb = np.asarray(isb).flatten()
        for ist in range(M):
            ind = int(sn.stationToNode[ist])
            if ind >= len(isb) or int(isb[ind]) != 1:
                continue  # no blocked marker at this station
            # destinations: stations j with connmatrix(ind, j) == 1
            dests = [j for j in range(sn.nnodes)
                     if sn.connmatrix[ind, j] == 1 and int(sn.isstation[j]) == 1]
            if len(dests) != 1:
                continue  # ambiguous destination: leave the job where it sits
            jst = int(sn.nodeToStation[dests[0]])
            if jst < 0:
                continue
            isf = int(sn.stationToStateful[ist])
            if isf not in sn.space or sn.space[isf] is None:
                continue
            space_isf = np.atleast_2d(sn.space[isf])
            if space_isf.shape[1] == 0:
                continue
            # marker is the trailing column of the station's local state block
            blocked = np.zeros(n_states, dtype=bool)
            for s_idx in range(n_states):
                row = int(state_space_hashed[s_idx, isf])
                if 0 <= row < space_isf.shape[0]:
                    blocked[s_idx] = space_isf[row, -1] == 1
            if not np.any(blocked):
                continue
            for k in range(K):
                col = ist * K + k
                if col >= state_space_aggr.shape[1]:
                    continue
                # a blocked state holds exactly ONE completed job, so cap the
                # per-state count at 1: only the held job moves, not the queue.
                counts = np.minimum(state_space_aggr[blocked, col], 1.0)
                shift = float(np.dot(pi[blocked], counts))
                if shift > 0:
                    QN[ist, k] -= shift
                    QN[jst, k] += shift

    # Synchronous calls (REPLY signals): a caller that is blocked awaiting a
    # reply still HOLDS its server and, by the LDES/LQN convention, still owns
    # the job it sent to the callee -- that simultaneous resource possession is
    # the point of the feature. The held servers live in the reply block at the
    # tail of the caller's local state, so add their time average to the
    # caller's utilization and queue length; without it CTMC reports only the
    # carried load (0.32915 against LDES 0.54997 on the closed client/server
    # model, with QLen 0.46395 against 0.68499).
    QNblocked = np.zeros((M, K))
    if getattr(sn, 'replyblock', None) is not None and np.any(np.asarray(sn.replyblock) > 0):
        from ...state.reply_block import reply_block_info
        for ist in range(M):
            ind = int(sn.stationToNode[ist])
            rinfo = reply_block_info(sn, ind)
            if rinfo.width == 0:
                continue
            isf = int(sn.stationToStateful[ist])
            space_isf = np.atleast_2d(sn.space[isf])
            rows = state_space_hashed[:, isf].astype(int)
            S_ist = float(sn.nservers[ist])
            pos = 0
            for r in rinfo.classes:
                col = space_isf.shape[1] - rinfo.width + pos
                pos += 1
                bmean = float(np.dot(pi, space_isf[rows, col]))
                QN[ist, r] += bmean
                UN[ist, r] += bmean / S_ist
                QNblocked[ist, r] = bmean

    # Response time via Little's law
    for ist in range(M):
        for k in range(K):
            if TN[ist, k] > 1e-14:
                # Response time is time spent AT the station, so the job blocked
                # out at the callee is excluded even though QLen/Util count it
                # (LDES measures the sojourn directly and reports the same).
                RN[ist, k] = (QN[ist, k] - QNblocked[ist, k]) / TN[ist, k]

    return {'Q': QN, 'U': UN, 'R': RN, 'T': TN}


def _update_cache_routing_probabilities(sn: NetworkStruct) -> None:
    """
    Update routing matrix with actual cache hit/miss probabilities.

    The default routing matrix uses 0.5/0.5 split for cache hit/miss.
    This function computes the actual probabilities using cache_xi_fp
    and updates sn.rtnodes accordingly.

    This matches the logic in MVA solver's _update_cache_routing method.

    Args:
        sn: Network structure to update
    """
    if sn.rtnodes is None:
        return

    I = sn.nnodes
    K = sn.nclasses

    # Find cache nodes
    cache_indices = []
    for ind in range(I):
        if ind < len(sn.nodetype) and sn.nodetype[ind] == NodeType.CACHE:
            cache_indices.append(ind)

    if not cache_indices:
        return

    # Update routing for each cache
    for ind in cache_indices:
        ch = sn.nodeparam.get(ind) if sn.nodeparam else None
        if ch is None:
            continue

        hitclass = getattr(ch, 'hitclass', None)
        missclass = getattr(ch, 'missclass', None)
        if hitclass is None or missclass is None:
            continue

        hitclass = np.atleast_1d(hitclass).flatten()
        missclass = np.atleast_1d(missclass).flatten()

        # Find input classes (classes that have hit/miss mappings)
        input_classes = []
        for r in range(K):
            if r < len(hitclass) and r < len(missclass):
                hc = int(hitclass[r]) if hitclass[r] >= 0 else -1
                mc = int(missclass[r]) if missclass[r] >= 0 else -1
                if hc >= 0 and mc >= 0:
                    input_classes.append(r)

        if not input_classes:
            continue

        # Get cache parameters
        n_items = getattr(ch, 'nitems', 5)
        m_cap = getattr(ch, 'cap', None)
        if m_cap is None:
            m_cap = np.array([1])
        m_cap = np.atleast_1d(m_cap)

        # Build gamma matrix from pread
        pread = getattr(ch, 'pread', None)
        if pread is None:
            # Default uniform access
            pread_arr = np.full(n_items, 1.0 / n_items)
        elif isinstance(pread, dict):
            # pread[class_idx] = probability array
            pread_arr = pread.get(input_classes[0], np.full(n_items, 1.0 / n_items))
            pread_arr = np.atleast_1d(pread_arr)
        elif isinstance(pread, (list, tuple)):
            pread_arr = np.atleast_1d(pread[input_classes[0]] if input_classes[0] < len(pread) else pread[0])
        else:
            pread_arr = np.atleast_1d(pread)

        # Ensure pread has correct length
        if len(pread_arr) != n_items:
            pread_arr = np.full(n_items, 1.0 / n_items)

        # Build gamma (n_items x h) where h = len(m_cap)
        h = len(m_cap)
        gamma = np.zeros((n_items, h))
        for i in range(n_items):
            for l in range(h):
                gamma[i, l] = pread_arr[i]  # Same popularity at each level for simple models

        # Compute hit/miss probabilities using cache_xi_fp
        try:
            xi, pi0, pij, it_fp = cache_xi_fp(gamma, m_cap)

            # Overall miss rate = sum(pread * pi0) where pi0 is miss prob per item
            overall_miss_rate = np.sum(pread_arr * pi0)
            overall_hit_rate = 1 - overall_miss_rate
        except Exception:
            # Fall back to default
            overall_hit_rate = 0.5
            overall_miss_rate = 0.5

        # Update routing matrix
        for r in input_classes:
            hc = int(hitclass[r])
            mc = int(missclass[r])

            # Zero out the row for input class at cache
            sn.rtnodes[ind * K + r, :] = 0

            # Find connected nodes
            for jnd in range(I):
                if sn.connmatrix is not None and ind < sn.connmatrix.shape[0] and jnd < sn.connmatrix.shape[1]:
                    if sn.connmatrix[ind, jnd]:
                        # Route to hit class with hit probability
                        if hc >= 0 and hc < K:
                            sn.rtnodes[ind * K + r, jnd * K + hc] = overall_hit_rate
                        # Route to miss class with miss probability
                        if mc >= 0 and mc < K:
                            sn.rtnodes[ind * K + r, jnd * K + mc] = overall_miss_rate

        # Store hit/miss probabilities for result reporting
        ch.actualhitprob = np.zeros(K)
        ch.actualmissprob = np.zeros(K)
        for r in input_classes:
            ch.actualhitprob[r] = overall_hit_rate
            ch.actualmissprob[r] = overall_miss_rate

    # Cache routing updated directly on sn.rt (not recomputed from rtnodes, whose Sink->Source visit-ratio edge would misdirect traffic back to Source).
    if sn.rt is None:
        return

    # Find next hop stateful nodes from Cache for hit/miss classes
    # This is done by looking at existing routing in rt for hit/miss classes
    for ind in cache_indices:
        # Get Cache's stateful index
        cache_sf = int(sn.nodeToStateful[ind]) if sn.nodeToStateful is not None and ind < len(sn.nodeToStateful) else -1
        if cache_sf < 0:
            continue

        ch = sn.nodeparam.get(ind) if sn.nodeparam else None
        if ch is None:
            continue

        hitclass = getattr(ch, 'hitclass', None)
        missclass = getattr(ch, 'missclass', None)
        if hitclass is None or missclass is None:
            continue

        hitclass = np.atleast_1d(hitclass).flatten()
        missclass = np.atleast_1d(missclass).flatten()

        # Get the actual hit/miss probabilities that were computed
        actualhitprob = getattr(ch, 'actualhitprob', None)
        actualmissprob = getattr(ch, 'actualmissprob', None)
        if actualhitprob is None or actualmissprob is None:
            continue

        # retrieval-system cache routing rows are left as link() built them (collapsing to hit/miss-only would erase the Cache->Queue edges and orphan the retrieval stations).
        rc_mat = getattr(ch, 'retrieval_classes', None)
        retrieval_input_classes = set()
        if rc_mat is not None:
            rc_arr = np.atleast_2d(np.asarray(rc_mat))
            for r in range(rc_arr.shape[1]):
                if np.any(rc_arr[:, r] >= 0):
                    retrieval_input_classes.add(r)

        # For each input class, update rt routing from Cache
        for r in range(K):
            if r >= len(hitclass) or r >= len(missclass):
                continue
            if r in retrieval_input_classes:
                continue

            hc = int(hitclass[r]) if hitclass[r] >= 0 else -1
            mc = int(missclass[r]) if missclass[r] >= 0 else -1
            if hc < 0 or mc < 0:
                continue

            hit_prob = actualhitprob[r] if r < len(actualhitprob) else 0.5
            miss_prob = actualmissprob[r] if r < len(actualmissprob) else 0.5

            # Source indices in rt (stateful-indexed)
            input_src_idx = cache_sf * K + r
            hit_src_idx = cache_sf * K + hc
            miss_src_idx = cache_sf * K + mc

            if input_src_idx >= sn.rt.shape[0]:
                continue

            # Zero out input class routing first
            sn.rt[input_src_idx, :] = 0

            # Combine hit and miss class routing with their probabilities
            for dst_idx in range(sn.rt.shape[1]):
                hit_route = sn.rt[hit_src_idx, dst_idx] if hit_src_idx < sn.rt.shape[0] else 0
                miss_route = sn.rt[miss_src_idx, dst_idx] if miss_src_idx < sn.rt.shape[0] else 0

                combined_prob = hit_prob * hit_route + miss_prob * miss_route
                if combined_prob > 1e-10:
                    sn.rt[input_src_idx, dst_idx] = combined_prob


def _forward_reachable(Q, start):
    """State indices forward-reachable from ``start`` over the nonzero off-diagonal
    entries of the generator ``Q`` (an ME embeds with negative ones, which are arcs
    all the same). For a closed pass-and-swap network whose initial placement lies in
    a recurrent class, this is that recurrent class."""
    Qd = np.asarray(Q)
    n = Qd.shape[0]
    visited = np.zeros(n, dtype=bool)
    visited[start] = True
    stack = [int(start)]
    while stack:
        u = stack.pop()
        row = np.abs(Qd[u])
        for v in np.nonzero(row > 1e-12)[0]:
            v = int(v)
            if v != u and not visited[v]:
                visited[v] = True
                stack.append(v)
    return np.nonzero(visited)[0]


def _initial_state_index(sn, state_space_hashed):
    """Global index of the initial state (sn.state) within the hashed state
    space, by matching each node's initial state into its per-node space.

    Used wherever a reducible generator has to be restricted to the recurrent
    class the model actually starts in (closed pass-and-swap placements,
    synchronous-call blocked-server counters)."""
    if not hasattr(sn, 'state') or sn.state is None or not hasattr(sn, 'space') or sn.space is None:
        return None
    nstateful = state_space_hashed.shape[1]
    init_hashed = np.zeros(nstateful, dtype=int)
    for isf in range(nstateful):
        node_state = sn.state.get(isf) if hasattr(sn.state, 'get') else sn.state[isf]
        node_space = sn.space.get(isf) if hasattr(sn.space, 'get') else sn.space[isf]
        if node_state is None or node_space is None:
            return None
        ns = np.atleast_2d(np.asarray(node_space, dtype=float))
        target = np.atleast_1d(np.asarray(node_state, dtype=float)).ravel()
        # The per-station initial row carries only as many buffer slots as the
        # initial population needs, while the enumerated local space is sized
        # for the full capacity. Left-pad to the space width (empty buffer slots
        # pad the left, so the server-phase and local-variable tail stays
        # aligned); without this the lookup fails and any pruning keyed on the
        # initial state is silently skipped.
        if ns.shape[1] > target.shape[0]:
            target = np.concatenate([np.zeros(ns.shape[1] - target.shape[0]), target])
        idx = -1
        for r in range(ns.shape[0]):
            if ns.shape[1] == target.shape[0] and np.allclose(ns[r], target):
                idx = r
                break
        if idx < 0:
            return None
        init_hashed[isf] = idx
    for s in range(state_space_hashed.shape[0]):
        if np.array_equal(state_space_hashed[s].astype(int), init_hashed):
            return s
    return None


def _build_fcr_waitq_ssg(sn, options):
    """
    Reachability-based augmented state space and per-action rate filters for
    finite capacity regions with the WAITQ (waiting queue) drop rule.

    Port of MATLAB solver_ctmc_fcr_waitq.m. JMT semantics: a job refused
    entry to a full region leaves the upstream station and waits in a
    per-region FIFO of (class, destination) tokens; releases are strictly
    FIFO head-of-line after every transition that frees capacity; a fresh
    cap-admissible arrival overtakes a stuck head; blocked jobs are counted
    at no station. Class-switching hops between members of the same region
    are an exit (with release) followed by a gated re-entry. DROP classes
    keep transition censoring.

    The augmented CTMC state is [h(0:nstateful-1), buf_0, ..., buf_{F-1}]
    with buf_f the token FIFO of region f padded with -1 to its max length.

    Returns:
        (state_space, state_space_aggr, state_space_hashed_aug, Q,
         arvRates, depRates, Dfilt)
    """
    from ...state.after_event import after_event_hashed, build_space_hash
    from ...state.marginal import toMarginal, toMarginalAggr
    from ....lang.sync import refresh_sync
    from ....constants import EventType, NodeType
    from ....lang.base import DropStrategy as BaseDropStrategy

    nstateful = int(sn.nstateful)
    K = int(sn.nclasses)
    nnodes = int(sn.nnodes)
    F = int(sn.nregions)

    # feature gates: combinations whose semantics are not defined here
    if getattr(sn, 'fjsync', None) is not None and len(sn.fjsync) > 0:
        raise RuntimeError('WAITQ finite capacity regions are not supported together with fork-join in SolverCTMC.')
    if hasattr(sn, 'nodetype') and sn.nodetype is not None:
        for ind in range(nnodes):
            if sn.nodetype[ind] == NodeType.Transition:
                raise RuntimeError('WAITQ finite capacity regions are not supported together with stochastic Petri net transitions in SolverCTMC.')

    sync = refresh_sync(sn)
    A = len(sync)
    for act in sync:
        if callable(act.passive.prob):
            raise RuntimeError('WAITQ finite capacity regions are not supported together with state-dependent routing in SolverCTMC.')
    hash_maps = build_space_hash(sn)
    # true-BAS blocking is orthogonal to region rules and handled via the become-blocked arcs under sentinel action -1, as in the default generator.
    isb_nodes = getattr(sn, 'isbasblocking', None)
    if isb_nodes is not None:
        isb_nodes = np.asarray(isb_nodes).ravel()
        if isb_nodes.size == 0:
            isb_nodes = None

    # region data
    M = int(sn.nstations)
    member_mask = np.zeros((F, M), dtype=bool)
    ccap = np.full((F, K), np.inf)
    gcap = np.full(F, np.inf)
    memcap = np.full(F, np.inf)
    szrow = np.ones((F, K))
    linA = [None] * F
    linb = [None] * F
    iswaitq = np.zeros((F, K), dtype=bool)
    drop_id = int(BaseDropStrategy.DROP)
    rr = np.asarray(sn.regionrule, dtype=float)
    for f in range(F):
        Rmat = np.asarray(sn.region[f], dtype=float)  # M x (K+1)
        memvec = -np.ones(M)
        if getattr(sn, 'regionmaxmem', None) and len(sn.regionmaxmem) > f and sn.regionmaxmem[f] is not None:
            memvec = np.asarray(sn.regionmaxmem[f], dtype=float).ravel()[:M]
            if memvec.size < M:
                memvec = np.concatenate([memvec, -np.ones(M - memvec.size)])
        members = [i for i in range(M) if np.any(Rmat[i, :] != -1) or memvec[i] != -1]
        member_mask[f, members] = True
        for r in range(K):
            cv = [Rmat[i, r] for i in members if Rmat[i, r] != -1]
            if cv:
                ccap[f, r] = min(cv)
            iswaitq[f, r] = (rr[f, r] != drop_id)
        gv = [Rmat[i, K] for i in members if Rmat[i, K] != -1]
        if gv:
            gcap[f] = min(gv)
        mv = [memvec[i] for i in members if memvec[i] != -1]
        if mv:
            memcap[f] = min(mv)
        if getattr(sn, 'regionsz', None) is not None and np.size(sn.regionsz) > 0:
            szrow[f, :] = np.asarray(sn.regionsz, dtype=float)[f].ravel()[:K]
        if (getattr(sn, 'regionlincon', None) and len(sn.regionlincon) > f
                and sn.regionlincon[f] is not None):
            linA[f] = np.atleast_2d(np.asarray(sn.regionlincon[f][0], dtype=float))
            linb[f] = np.asarray(sn.regionlincon[f][1], dtype=float).ravel()

    def violates(f, x):
        if np.any(x > ccap[f, :]) or x.sum() > gcap[f] or float(x @ szrow[f, :]) > memcap[f]:
            return True
        if linA[f] is not None:
            return bool(np.any(linA[f] @ x.reshape(-1, 1) > linb[f].reshape(-1, 1)))
        return False

    # token FIFO length bound: closed chains bound by the chain population
    # (jobs may switch into a class with njobs=0); open classes by the cutoff
    cutoff = options.cutoff if getattr(options, 'cutoff', None) is not None else 0
    if np.isscalar(cutoff):
        cutoff_mat = np.full((M, K), float(cutoff))
    else:
        cutoff_mat = np.atleast_2d(np.asarray(cutoff, dtype=float))
    chains = np.atleast_2d(np.asarray(sn.chains, dtype=float))
    njobs = np.asarray(sn.njobs, dtype=float).ravel()
    tokbound = np.zeros(K)
    for r in range(K):
        if np.any(iswaitq[:, r]):
            crow = np.flatnonzero(chains[:, r])
            chainpop = njobs[chains[crow[0], :] > 0].sum() if crow.size else njobs[r]
            if np.isfinite(chainpop):
                tokbound[r] = chainpop
            else:
                tokbound[r] = cutoff_mat[:, r].max() if cutoff_mat.size else 0
    Lmax = np.zeros(F, dtype=int)
    for f in range(F):
        Lmax[f] = int(tokbound[iswaitq[f, :]].sum())
    bufoff = nstateful + np.concatenate([[0], np.cumsum(Lmax[:-1])]).astype(int)
    width = nstateful + int(Lmax.sum())

    # initial augmented state (buffers empty, tokens padded with -1)
    h0 = np.zeros(nstateful, dtype=int)
    for ind in range(nnodes):
        if not sn.isstateful[ind]:
            continue
        isf = int(sn.nodeToStateful[ind])
        spc = np.atleast_2d(np.asarray(sn.space[isf], dtype=float))
        if spc.shape[0] == 1:
            # single local state (e.g. Source): nothing to match
            h0[isf] = 0
            continue
        node_state = sn.state.get(isf) if hasattr(sn.state, 'get') else sn.state[isf]
        if node_state is None:
            raise RuntimeError('WAITQ finite capacity regions need the initial state (sn.state) to be set.')
        row = np.atleast_2d(np.asarray(node_state, dtype=float))[0].ravel()
        w = spc.shape[1]
        if row.size < w:
            row = np.concatenate([np.zeros(w - row.size), row])
        key = tuple(row.astype(float))
        idx = hash_maps[isf].get(key, -1) if isf in hash_maps else -1
        if idx < 0:
            # fall back to a row match
            for rix in range(spc.shape[0]):
                if np.allclose(spc[rix], row):
                    idx = rix
                    break
        if idx < 0:
            raise RuntimeError('Initial state of a stateful node not found in its local state space (WAITQ FCR).')
        h0[isf] = idx
    row0 = -np.ones(width, dtype=int)
    row0[:nstateful] = h0

    def region_aggr(hvec, f):
        x = np.zeros(K)
        for ist in np.flatnonzero(member_mask[f, :]):
            ind = int(sn.stationToNode[ist])
            isf = int(sn.nodeToStateful[ind])
            srow = np.atleast_2d(sn.space[isf])[int(hvec[isf]):int(hvec[isf]) + 1]
            _, nir = toMarginalAggr(sn, ind, srow)
            x += np.asarray(nir, dtype=float).ravel()[:K]
        return x

    # BFS structures
    keymap = {tuple(row0): 0}
    SSH = [row0.copy()]
    frontier = [0]
    trips = []  # (a, src, dst, w)

    def emit(a, src, newh, newbufs, w, pend=None):
        """Release cascade + optional pending gated re-entry, then register."""
        if w <= 0:
            return
        work = [(newh.copy(), [list(b) for b in newbufs], 1.0, pend)]
        while work:
            hh, bb, pw, pd = work.pop(0)
            progressed = False
            for f in range(F):
                if not bb[f]:
                    continue
                x = region_aggr(hh, f)
                tok = bb[f][0]
                dest = tok // K
                r = tok % K
                xn = x.copy()
                xn[r] += 1
                if violates(f, xn):
                    continue  # head-of-line: this FIFO stays blocked
                isf_d = int(sn.nodeToStateful[dest])
                hd, _, opd = after_event_hashed(sn, dest, int(hh[isf_d]), EventType.ARV, r, hash_maps)
                if hd is None or np.all(np.asarray(hd) < 0):
                    continue
                hd = np.asarray(hd).ravel()
                opd = np.asarray(opd, dtype=float).ravel()
                for idd in range(hd.size):
                    if hd[idd] < 0 or (idd < opd.size and opd[idd] <= 0):
                        continue
                    hh2 = hh.copy()
                    hh2[isf_d] = int(hd[idd])
                    bb2 = [list(b) for b in bb]
                    bb2[f].pop(0)
                    work.append((hh2, bb2, pw * (opd[idd] if idd < opd.size else 1.0), pd))
                progressed = True
                break
            if progressed:
                continue
            if pd is not None:
                # pending gated re-entry (class-switching intra-region hop)
                f, cls, dest = pd[0], pd[1], pd[2]
                pd_waitq = pd[3] if len(pd) > 3 else True
                x = region_aggr(hh, f)
                xn = x.copy()
                xn[cls] += 1
                admitted = False
                if not violates(f, xn):
                    isf_d = int(sn.nodeToStateful[dest])
                    hd, _, opd = after_event_hashed(sn, dest, int(hh[isf_d]), EventType.ARV, cls, hash_maps)
                    if hd is not None and not np.all(np.asarray(hd) < 0):
                        hd = np.asarray(hd).ravel()
                        opd = np.asarray(opd, dtype=float).ravel()
                        for idd in range(hd.size):
                            if hd[idd] < 0 or (idd < opd.size and opd[idd] <= 0):
                                continue
                            hh2 = hh.copy()
                            hh2[isf_d] = int(hd[idd])
                            work.append((hh2, [list(b) for b in bb], pw * (opd[idd] if idd < opd.size else 1.0), None))
                            admitted = True
                if not admitted:
                    if pd_waitq:
                        bb2 = [list(b) for b in bb]
                        bb2[f].append(dest * K + cls)
                        work.append((hh.copy(), bb2, pw, None))
                    else:
                        # DROP: the switching job is destroyed
                        work.append((hh.copy(), [list(b) for b in bb], pw, None))
                continue
            # settled: register augmented state and transition
            rr_ = -np.ones(width, dtype=int)
            rr_[:nstateful] = hh
            for f in range(F):
                for j, tok in enumerate(bb[f]):
                    rr_[bufoff[f] + j] = tok
            kk = tuple(rr_)
            if kk in keymap:
                dst = keymap[kk]
            else:
                dst = len(SSH)
                keymap[kk] = dst
                SSH.append(rr_.copy())
                frontier.append(dst)
            trips.append((a, src, dst, w * pw))

    if any(violates(f, region_aggr(h0, f)) for f in range(F)):
        raise RuntimeError('The initial state violates the finite capacity region constraints.')

    while frontier:
        s = frontier.pop(0)
        row = SSH[s]
        h = row[:nstateful].copy()
        bufs = []
        for f in range(F):
            bf = row[bufoff[f]:bufoff[f] + Lmax[f]]
            bufs.append([int(t) for t in bf if t >= 0])
        xf = np.array([region_aggr(h, f) for f in range(F)]) if F > 0 else np.zeros((0, K))
        for a in range(A):
            act = sync[a]
            node_a = int(act.active.node)
            if not sn.isstateful[node_a]:
                continue
            isf_a = int(sn.nodeToStateful[node_a])
            class_a = int(act.active.job_class)
            event_a = act.active.event
            new_state_a, rate_a, _ = after_event_hashed(sn, node_a, int(h[isf_a]), event_a, class_a, hash_maps)
            if new_state_a is None or np.all(np.asarray(new_state_a) < 0):
                continue
            new_state_a = np.asarray(new_state_a).ravel()
            rate_a = np.asarray(rate_a, dtype=float).ravel()
            node_p = int(act.passive.node)
            is_local = node_p >= nnodes
            for ia in range(new_state_a.size):
                if new_state_a[ia] < 0 or not np.isfinite(rate_a[ia]) or rate_a[ia] <= 0:
                    continue
                if is_local:
                    newh = h.copy()
                    newh[isf_a] = int(new_state_a[ia])
                    emit(a, s, newh, bufs, float(rate_a[ia]))
                    continue
                class_p = int(act.passive.job_class)
                event_p = act.passive.event
                isf_p = int(sn.nodeToStateful[node_p]) if sn.isstateful[node_p] else -1
                if isf_p < 0:
                    continue
                stat_a = int(sn.nodeToStation[node_a]) if sn.isstation[node_a] else -1
                stat_p = int(sn.nodeToStation[node_p]) if sn.isstation[node_p] else -1
                blockedf = -1
                droppedf = -1
                switchf = -1
                if event_p == EventType.ARV and stat_p >= 0:
                    for f in range(F):
                        if member_mask[f, stat_p] and (stat_a < 0 or not member_mask[f, stat_a]):
                            xn = xf[f].copy()
                            xn[class_p] += 1
                            if violates(f, xn):
                                if not iswaitq[f, class_p]:
                                    droppedf = f  # DROP: the job is destroyed
                                else:
                                    blockedf = f
                                break
                        elif (member_mask[f, stat_p] and stat_a >= 0 and member_mask[f, stat_a]
                              and class_p != class_a):
                            switchf = f
                            break
                prob_p = float(act.passive.prob) if act.passive.prob is not None else 1.0
                if droppedf >= 0:
                    # DROP rule (JMT): only the active part applies, job vanishes
                    newh = h.copy()
                    newh[isf_a] = int(new_state_a[ia])
                    emit(a, s, newh, bufs, float(rate_a[ia]) * prob_p)
                    continue
                if switchf >= 0:
                    newh = h.copy()
                    newh[isf_a] = int(new_state_a[ia])
                    emit(a, s, newh, bufs, float(rate_a[ia]) * prob_p,
                         (switchf, class_p, node_p, bool(iswaitq[switchf, class_p])))
                    continue
                if blockedf >= 0:
                    if len(bufs[blockedf]) >= Lmax[blockedf]:
                        continue  # FIFO truncation boundary (open-class cutoff)
                    newh = h.copy()
                    newh[isf_a] = int(new_state_a[ia])
                    newbufs = [list(b) for b in bufs]
                    newbufs[blockedf].append(node_p * K + class_p)
                    emit(a, s, newh, newbufs, float(rate_a[ia]) * prob_p)
                    continue
                # normal passive application
                if node_p == node_a:
                    state_p = int(new_state_a[ia])
                else:
                    state_p = int(h[isf_p])
                new_state_p, _, outprob_p = after_event_hashed(sn, node_p, state_p, event_p, class_p, hash_maps)
                refused = new_state_p is None or np.all(np.asarray(new_state_p).ravel() < 0)
                if refused:
                    # become-blocked arc on refusal: hold the completed job at the server rather than voiding the departure (voiding would let the server re-serve immediately, understating throughput vs a fresh completion draw). See _kb/06-solver-catalog.md True BAS blocking section.
                    if (event_a == EventType.DEP and isb_nodes is not None
                            and node_a < isb_nodes.size and int(isb_nodes[node_a]) == 1):
                        space_a = np.atleast_2d(sn.space[isf_a])
                        cur_vec_a = space_a[int(h[isf_a])]
                        if cur_vec_a.size > 0 and cur_vec_a[-1] == 0:
                            blocked_vec = cur_vec_a.copy()
                            blocked_vec[-1] = 1
                            blocked_idx = hash_maps.get(isf_a, {}).get(
                                tuple(np.asarray(blocked_vec, dtype=float)), -1)
                            if blocked_idx >= 0:
                                newh = h.copy()
                                newh[isf_a] = int(blocked_idx)
                                emit(-1, s, newh, bufs, float(rate_a[ia]) * prob_p)
                    continue
                new_state_p = np.asarray(new_state_p).ravel()
                outprob_p = np.asarray(outprob_p, dtype=float).ravel()
                for ip in range(new_state_p.size):
                    if new_state_p[ip] < 0:
                        continue
                    psync = prob_p * (outprob_p[ip] if ip < outprob_p.size else 1.0)
                    if psync <= 0:
                        continue
                    newh = h.copy()
                    newh[isf_a] = int(new_state_a[ia])
                    newh[isf_p] = int(new_state_p[ip])
                    emit(a, s, newh, bufs, float(rate_a[ia]) * psync)

    n = len(SSH)
    SSH_arr = np.asarray(SSH, dtype=int)
    Dfilt = [np.zeros((n, n)) for _ in range(A)]
    # Sentinel action -1 collects the true-BAS become-blocked arcs: part of the
    # generator, but not a departure of any action, so kept out of Dfilt.
    basBlockQ = np.zeros((n, n))
    for (a, i, j, v) in trips:
        if a < 0:
            basBlockQ[i, j] += v
        else:
            Dfilt[a][i, j] += v
    Q = np.eye(n)
    for a in range(A):
        Q = Q + Dfilt[a]
    Q = Q + basBlockQ

    arvRates = np.zeros((n, nstateful, K))
    depRates = np.zeros((n, nstateful, K))
    for a in range(A):
        act = sync[a]
        node_a = int(act.active.node)
        if act.active.event == EventType.DEP and sn.isstateful[node_a]:
            isf_a = int(sn.nodeToStateful[node_a])
            depRates[:, isf_a, int(act.active.job_class)] += Dfilt[a].sum(axis=1)
            node_p = int(act.passive.node)
            if node_p < nnodes and sn.isstateful[node_p]:
                isf_p = int(sn.nodeToStateful[node_p])
                arvRates[:, isf_p, int(act.passive.job_class)] += Dfilt[a].sum(axis=1)

    # full state matrix (concatenated per-node states + buffer columns) and aggr
    cols = 0
    widths = []
    for ind in range(nnodes):
        if sn.isstateful[ind]:
            isf = int(sn.nodeToStateful[ind])
            w = np.atleast_2d(sn.space[isf]).shape[1]
            widths.append((ind, isf, w))
            cols += w
    state_space = np.zeros((n, cols + int(Lmax.sum())))
    state_space_aggr = np.zeros((n, M * K))
    for s in range(n):
        pos = 0
        for (ind, isf, w) in widths:
            srow = np.atleast_2d(sn.space[isf])[SSH_arr[s, isf]:SSH_arr[s, isf] + 1]
            state_space[s, pos:pos + w] = srow.ravel()
            pos += w
            if sn.isstation[ind]:
                ist = int(sn.nodeToStation[ind])
                try:
                    result = toMarginal(sn, ind, srow)
                    nir_vec = np.asarray(result[1], dtype=float).ravel()[:K]
                    state_space_aggr[s, ist * K:(ist + 1) * K] = nir_vec
                except Exception:
                    pass
        state_space[s, cols:] = SSH_arr[s, nstateful:]

    return state_space, state_space_aggr, SSH_arr, Q, arvRates, depRates, Dfilt


def solver_ctmc_basic(
    sn: NetworkStruct,
    options: Optional[SolverCTMCOptions] = None
) -> SolverCTMCReturn:
    """
    Basic CTMC solver using state-space enumeration.

    Enumerates all valid states, builds the infinitesimal generator,
    and solves for steady-state distribution.

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        SolverCTMCReturn with performance metrics
    """
    start_time = time.time()

    if options is None:
        options = SolverCTMCOptions()

    M = sn.nstations
    K = sn.nclasses

    # Check state space size estimate (matches MATLAB's runAnalyzer.m)
    # Uses gamma function to estimate worst-case state space size
    NK = sn.njobs if sn.njobs is not None else np.ones(K)
    size_estimator = 0.0
    from scipy.special import gammaln
    for k in range(K):
        if np.isfinite(NK[k]):
            # log(C(NK[k]+M-1, M-1)) = gammaln(NK[k]+M) - gammaln(M) - gammaln(NK[k]+1)
            size_estimator += gammaln(1 + NK[k] + M - 1) - gammaln(1 + M - 1) - gammaln(1 + NK[k])

    if size_estimator > 6 and not options.force:
        raise RuntimeError(
            f"CTMC state space may be too large to solve (size estimate: {size_estimator:.1f} > 6). "
            f"Use force=True option to bypass this check, e.g., CTMC(model, force=True). "
            f"Alternative solvers: MVA, NC, or SSA."
        )

    # LCFS+LCFSPR networks use the standard state-space/generator path (buffer ordering already handles LCFS correctly).

    # Detect Cache nodes early (needed for visits refresh and rt recomputation)
    has_cache = False
    if hasattr(sn, 'nodetype') and sn.nodetype is not None:
        for ind in range(int(sn.nnodes)):
            if ind < len(sn.nodetype):
                nt = sn.nodetype[ind]
                nt_val = int(nt.value) if hasattr(nt, 'value') else int(nt)
                if nt_val == 6:  # NodeType.CACHE
                    has_cache = True
                    break

    # Update cache routing probabilities (replaces default 0.5/0.5 split)
    # This must be done before building the generator matrix
    _update_cache_routing_probabilities(sn)

    # Refresh visit ratios with updated cache routing probabilities
    # (matches MATLAB refreshChains which recomputes visits after cache prob update)
    if has_cache:
        sn_refresh_visits(sn)

    if GlobalConstants.getVerbose() == VerboseLevel.DEBUG and sn.rt is not None:
        rt = np.asarray(sn.rt)
        nstateful = sn.nstateful if hasattr(sn, 'nstateful') else M
        print("  sn.rt after cache routing update (non-zero entries):")
        for i in range(rt.shape[0]):
            sf_i = i // K
            k_i = i % K
            for j in range(rt.shape[1]):
                if abs(rt[i, j]) > 1e-10:
                    sf_j = j // K
                    k_j = j % K
                    print(f"    rt[sf{sf_i},c{k_i} -> sf{sf_j},c{k_j}] = {rt[i,j]:.6f}")

    # sn.rt recomputed from rtnodes via stochastic complementation, except Cache networks (already updated directly, since rtnodes' Sink->Source edge would misdirect traffic).
    if sn.rtnodes is not None and not has_cache:
        I = sn.nnodes
        K = sn.nclasses
        nstateful = sn.nstateful if hasattr(sn, 'nstateful') else M

        # Get stateful node indices from stationToStateful and Router nodes
        stateful_nodes = set()
        if hasattr(sn, 'stationToStateful') and sn.stationToStateful is not None:
            for sf in sn.stationToStateful:
                stateful_nodes.add(int(sf))
        # Add any Router nodes (they are stateful but not stations)
        if hasattr(sn, 'nodetype') and sn.nodetype is not None:
            node_to_stateful = {}
            sf_idx = 0
            for nidx in range(I):
                # Check if this node is stateful (station or stateful non-station)
                if hasattr(sn, 'nodeToStateful') and sn.nodeToStateful is not None:
                    if nidx < len(sn.nodeToStateful):
                        sf = sn.nodeToStateful[nidx]
                        if sf >= 0:
                            stateful_nodes.add(int(sf))
                elif hasattr(sn, 'isstateful') and sn.isstateful is not None:
                    if nidx < len(sn.isstateful) and sn.isstateful[nidx]:
                        stateful_nodes.add(sf_idx)
                        sf_idx += 1

        # Build stateful_nodes_classes for dtmc_stochcomp
        stateful_node_list = sorted(stateful_nodes)
        if not stateful_node_list:
            # Fallback: use station indices
            stateful_node_list = list(range(M))

        # (stateful_node_idx, class_idx) pairs built for dtmc_stochcomp, which indexes into rtnodes by node (not stateful) index.
        stateful_to_node = {}
        if hasattr(sn, 'statefulToNode') and sn.statefulToNode is not None:
            for sf_idx, n_idx in enumerate(sn.statefulToNode):
                stateful_to_node[sf_idx] = int(n_idx)

        # Build the list of node*K + class indices for stateful nodes
        stateful_nodes_classes = []
        for sf_idx in range(nstateful):
            node_idx = stateful_to_node.get(sf_idx, sf_idx)
            for k in range(K):
                stateful_nodes_classes.append(node_idx * K + k)
        stateful_nodes_classes = np.array(stateful_nodes_classes, dtype=int)

        try:
            new_rt = dtmc_stochcomp(sn.rtnodes, stateful_nodes_classes)
            sn.rt = new_rt
            # CRITICAL: Also update rt_visits since sn_refresh_visits uses it
            if hasattr(sn, 'rt_visits') and sn.rt_visits is not None:
                sn.rt_visits = new_rt.copy()
        except Exception:
            pass  # Keep existing rt if stochcomp fails

    # Detect Router-like nodes (non-station stateful, non-Cache)
    router_nodes = []
    if sn_has_open_classes(sn) and has_cache:
        if hasattr(sn, 'isstateful') and hasattr(sn, 'isstation'):
            for ind in range(int(sn.nnodes)):
                if ind < len(sn.isstateful) and ind < len(sn.isstation):
                    is_stateful = sn.isstateful[ind]
                    is_station = sn.isstation[ind]
                    if is_stateful and not is_station:
                        if hasattr(sn, 'nodetype') and ind < len(sn.nodetype):
                            nt = sn.nodetype[ind]
                            nt_val = int(nt.value) if hasattr(nt, 'value') else int(nt)
                            if nt_val != 6:  # Not CACHE
                                router_nodes.append(ind)

    M = int(sn.nstations)
    K = int(sn.nclasses)

    # The sync (sync-action) builder is the CTMC state-space generator for every
    # model; it has been validated against MATLAB across the full example suite.
    options.gen_method = 'sync'

    # ---- Sync-action-based builder path ----
    if options.gen_method == 'sync':
        from ...state.ctmc_ssg import ctmc_ssg as ctmc_ssg_fn
        # generator built on a private sn copy: the sync builder mutates phase fields and per-node state spaces in place, which must not corrupt the struct shared with downstream solvers.
        import copy as _copy
        sn = _copy.deepcopy(sn)

        # Mandatory truncation warning for open/mixed models
        if sn_has_open_classes(sn):
            print(f"CTMC solver using state space cutoff = {options.cutoff} for open/mixed model.")
            warnings.warn(
                "State space truncation may cause inaccurate results. "
                "Consider varying cutoff to assess sensitivity.",
                UserWarning
            )

        # phase fields precomputed with immediate_as_rate=True so an Immediate-service class gets a fast mu (~1e8) and its DEP can fire; stochastic complementation still only folds immediate pass-through NODES, not station service; mutation stays on the private deepcopy.
        _refresh_phase_fields(sn, immediate_as_rate=True)

        # native fork-join breaks per-chain population conservation, so it uses the reachability-based generator instead of the population-lattice ctmc_ssg.
        _isFJ = getattr(sn, 'fjsync', None) is not None and len(sn.fjsync) > 0
        if _isFJ:
            from ...state.reachable_ssg import reachable_ssg
            state_space, state_space_aggr, state_space_hashed, sn = reachable_ssg(sn, options)
        else:
            state_space, state_space_aggr, state_space_hashed, sn = ctmc_ssg_fn(sn, options.cutoff)

        if state_space_hashed.size == 0:
            raise RuntimeError("CTMC sync builder: empty state space generated.")

        # FCR: DROP removes over-cap states (boundary transitions into them are dropped at lookup, exact for memoryless sources); WAITQ is inexact by censoring alone and uses the augmented _build_fcr_waitq_ssg generator instead; mirrors solver_ctmc_fcr_waitq.m.
        _fcr_waitq = bool(getattr(sn, 'nregions', 0) and int(sn.nregions) > 0)
        if _fcr_waitq:
            (state_space, state_space_aggr, state_space_hashed, Q, arvRates,
             depRates, Dfilt) = _build_fcr_waitq_ssg(sn, options)
            Dfilt_fjsync, _fj_sync = [], []
        if (not _fcr_waitq) and getattr(sn, 'nregions', 0) and int(sn.nregions) > 0 and getattr(sn, 'region', None):
            Kf = int(sn.nclasses)
            nS = state_space_hashed.shape[0]
            feas = np.ones(nS, dtype=bool)
            for f in range(int(sn.nregions)):
                Rmat = np.asarray(sn.region[f], dtype=float)  # M x (K+1)
                M = Rmat.shape[0]
                # membership: any job-count cap OR the region memory budget set
                # on the station row (a memory-only region has all caps at -1)
                memvec = -np.ones(M)
                if (getattr(sn, 'regionmaxmem', None) and len(sn.regionmaxmem) > f
                        and sn.regionmaxmem[f] is not None and np.size(sn.regionmaxmem[f]) > 0):
                    mv_ = np.asarray(sn.regionmaxmem[f], dtype=float).ravel()[:M]
                    memvec[:mv_.size] = mv_
                members = [i for i in range(M) if np.any(Rmat[i, :] != -1) or memvec[i] != -1]
                if not members:
                    continue
                xagg = np.zeros((nS, Kf))
                for i in members:
                    for r in range(Kf):
                        col = i * Kf + r
                        if col < state_space_aggr.shape[1]:
                            xagg[:, r] += state_space_aggr[:, col]
                for r in range(Kf):
                    capvals = [Rmat[i, r] for i in members if Rmat[i, r] != -1]
                    if capvals:
                        feas &= (xagg[:, r] <= min(capvals))
                gvals = [Rmat[i, Kf] for i in members if Rmat[i, Kf] != -1]
                if gvals:
                    feas &= (xagg.sum(axis=1) <= min(gvals))
                if (getattr(sn, 'regionmaxmem', None) and len(sn.regionmaxmem) > f
                        and sn.regionmaxmem[f] is not None and np.size(sn.regionmaxmem[f]) > 0):
                    memmat = np.asarray(sn.regionmaxmem[f], dtype=float).ravel()
                    mvals = [memmat[i] for i in members if i < len(memmat) and memmat[i] != -1]
                    if mvals:
                        if getattr(sn, 'regionsz', None) is not None and np.size(sn.regionsz) > 0:
                            szrow = np.asarray(sn.regionsz, dtype=float)[f].ravel()
                        else:
                            szrow = np.ones(Kf)
                        feas &= (xagg.dot(szrow[:Kf]) <= min(mvals))
            if not feas.all():
                if not feas.any():
                    raise RuntimeError("Finite Capacity Region constraints leave no feasible "
                                       "state; check the region caps and the initial state.")
                state_space = state_space[feas]
                state_space_aggr = state_space_aggr[feas]
                state_space_hashed = state_space_hashed[feas]

        if GlobalConstants.getVerbose() == VerboseLevel.DEBUG:
            print(f"CTMC sync state space size: {state_space_hashed.shape[0]} states, "
                  f"{state_space_hashed.shape[1]} stateful nodes")

        # Build generator matrix using sync actions
        if not _fcr_waitq:
            Q, arvRates, depRates, Dfilt, Dfilt_fjsync, _fj_sync, _DfiltGsyncComp = _build_generator_sync(
                sn, state_space, state_space_hashed, options)

        if GlobalConstants.getVerbose() == VerboseLevel.DEBUG:
            print(f"  Q matrix shape: {Q.shape}, nnz: {np.count_nonzero(Q)}")

        # Apply stochastic complementation to eliminate immediate states
        nonimm_indices, imm_indices = _find_immediate_states_sync(sn, state_space_hashed)
        _isFJ = getattr(sn, 'fjsync', None) is not None and len(sn.fjsync) > 0

        if imm_indices:
            Q_dense = np.asarray(Q)
            imm = np.array(imm_indices, dtype=int)
            nonimm = np.array(nonimm_indices, dtype=int)

            Q11 = Q_dense[np.ix_(nonimm, nonimm)]
            Q12 = Q_dense[np.ix_(nonimm, imm)]
            Q21 = Q_dense[np.ix_(imm, nonimm)]
            Q22 = Q_dense[np.ix_(imm, imm)]

            # -Q22 (the immediate-state sub-generator) is identical across every Dfilt action; factorized ONCE (sparse LU) and reused, avoiding an O(imm^3) dense solve per action that turns large immediate blocks into an effective hang.
            from scipy.sparse import csc_matrix as _csc
            from scipy.sparse.linalg import splu as _splu

            _lu = None
            try:
                _lu = _splu(_csc(-Q22))
            except (RuntimeError, ValueError):
                _lu = None
            _pinv_cache = {}

            def _stochcomp_solve(rhs):
                rhs = np.asarray(rhs, dtype=float)
                if _lu is not None:
                    out = np.zeros_like(rhs)
                    nz = np.flatnonzero(np.any(rhs != 0.0, axis=0))
                    if nz.size:
                        out[:, nz] = _lu.solve(rhs[:, nz])
                    return out
                if 'p' not in _pinv_cache:
                    _pinv_cache['p'] = np.linalg.pinv(-Q22)
                return _pinv_cache['p'] @ rhs

            Q_reduced = Q11 + Q12 @ _stochcomp_solve(Q21)

            # native fork-join rate recomputation: r_a = tangible exit rate via a + Q12*(-Q22)^-1*vanishing exit rate via a; must run BEFORE the Dfilt-overwrite loop below. Mirrors solver_ctmc.m:549-620.
            if _isFJ or imm.size:
                from ....constants import EventType
                nstateful_ = int(sn.nstateful)
                nclasses_ = int(sn.nclasses)
                arvRates = np.zeros((len(nonimm), nstateful_, nclasses_))
                depRates = np.zeros((len(nonimm), nstateful_, nclasses_))

                def _fj_rate(Dfull):
                    r_direct = np.asarray(Dfull[nonimm, :]).sum(axis=1).ravel()
                    r_chain_imm = np.asarray(Dfull[imm, :]).sum(axis=1).ravel()
                    if imm.size:
                        r_chain = Q12 @ _stochcomp_solve(r_chain_imm.reshape(-1, 1))
                        r_chain = np.asarray(r_chain).ravel()
                    else:
                        r_chain = np.zeros(len(nonimm))
                    return r_direct + r_chain

                # SPN IMMEDIATE-mode firings survive only through the chain term above (they fire only in vanishing markings).
                if getattr(sn, 'gsync', None):
                    for g_idx, glevent_rc in sn.gsync.items():
                        if not glevent_rc.active or glevent_rc.active[0].event != EventType.FIRE:
                            continue
                        if g_idx >= len(_DfiltGsyncComp):
                            continue
                        r_g = _fj_rate(_DfiltGsyncComp[g_idx])
                        for pev in glevent_rc.passive:
                            pev_node = int(pev.node)
                            if pev_node >= sn.nnodes or not sn.isstateful[pev_node]:
                                continue
                            pev_isf = int(sn.nodeToStateful[pev_node])
                            pev_class = int(pev.job_class)
                            if pev_class < 0 or pev_class >= nclasses_:
                                continue
                            if pev.event == EventType.PRE:
                                depRates[:, pev_isf, pev_class] += r_g
                            elif pev.event == EventType.POST:
                                arvRates[:, pev_isf, pev_class] += r_g

                for a in range(len(_fj_sync)):
                    act = _fj_sync[a]
                    if act.active.event != EventType.DEP:
                        continue
                    node_a = int(act.active.node)
                    if not sn.isstateful[node_a]:
                        continue
                    class_a = int(act.active.job_class)
                    node_p = int(act.passive.node)
                    node_a_sf = int(sn.nodeToStateful[node_a])
                    # a DEP is a tangible-only timed completion and must not be complemented over vanishing rows, EXCEPT an FJ Join's DEP (which fires only from the complete-sibling-set vanishing marking); see _kb/06-solver-catalog.md Vanishing states section (Join DEP paragraph).
                    _join_val = (int(NodeType.JOIN.value) if hasattr(NodeType.JOIN, 'value')
                                 else int(NodeType.JOIN))
                    if _isFJ and int(sn.nodetype[node_a]) == _join_val:
                        r_a = _fj_rate(Dfilt[a])
                    else:
                        r_a = np.asarray(Dfilt[a][nonimm, :]).sum(axis=1).ravel()
                    depRates[:, node_a_sf, class_a] += r_a
                    if node_p < sn.nnodes and sn.isstateful[node_p]:
                        node_p_sf = int(sn.nodeToStateful[node_p])
                        class_p = int(act.passive.job_class)
                        arvRates[:, node_p_sf, class_p] += r_a

                for k in range(len(Dfilt_fjsync)):
                    fjentry = sn.fjsync[k]
                    isf_fork = int(sn.nodeToStateful[int(fjentry['fork'])])
                    r_k = _fj_rate(Dfilt_fjsync[k])
                    depRates[:, isf_fork, int(fjentry['class'])] += r_k
                    branchheads = np.asarray(fjentry['branchheads']).ravel()
                    auxclasses = np.asarray(fjentry['auxclasses']).ravel()
                    for b in range(len(branchheads)):
                        isf_bh = int(sn.nodeToStateful[int(branchheads[b])])
                        arvRates[:, isf_bh, int(auxclasses[b])] += r_k

            # Also apply stochcomp to Dfilt for accurate rates
            for a in range(len(Dfilt)):
                Q21a = Dfilt[a][np.ix_(imm, nonimm)]
                Ta = Q12 @ _stochcomp_solve(Q21a)
                Dfilt[a] = Dfilt[a][np.ix_(nonimm, nonimm)] + Ta

            if not (_isFJ or imm.size):
                depRates = depRates[nonimm_indices, :, :]
                arvRates = arvRates[nonimm_indices, :, :]
            state_space = state_space[nonimm_indices]
            state_space_aggr = state_space_aggr[nonimm_indices]
            state_space_hashed = state_space_hashed[nonimm_indices]

            Q = Q_reduced

            if options.verbose:
                print(f"  Stochcomp: {len(imm_indices)} immediate states eliminated, "
                      f"{len(nonimm_indices)} non-immediate states remain.")

        # REPLY blocked-server counters are enumerated independently of the marginals, so unreachable (absorbing) configurations are pruned from the initial state; mirrors MATLAB solver_ctmc.m.
        if (getattr(sn, 'replyblock', None) is not None
                and np.any(np.asarray(sn.replyblock) > 0) and Q.shape[0] > 1):
            _init_idx = _initial_state_index(sn, state_space_hashed)
            if _init_idx is None or _init_idx < 0:
                raise RuntimeError(
                    "Synchronous calls (REPLY signals): the initial state was not found in the "
                    "enumerated state space, so the unreachable blocked-server configurations "
                    "cannot be pruned and the generator would be reducible.")
            _reach = _forward_reachable(Q, _init_idx)
            if 0 < len(_reach) < Q.shape[0]:
                _reach = np.asarray(sorted(int(x) for x in _reach), dtype=int)
                Q = Q[np.ix_(_reach, _reach)]
                state_space = state_space[_reach]
                state_space_aggr = state_space_aggr[_reach]
                state_space_hashed = state_space_hashed[_reach]
                arvRates = arvRates[_reach]
                depRates = depRates[_reach]
                Dfilt = [d[np.ix_(_reach, _reach)] for d in Dfilt]

        # closed PAS networks are reducible (placement order is conserved); the generator is restricted to the recurrent class reachable from the initial placement before solving.
        _pas_present = False
        if hasattr(sn, 'sched') and sn.sched is not None:
            for _ist in range(int(sn.nstations)):
                if sn.sched[_ist] == SchedStrategy.PAS:
                    _pas_present = True
                    break
        if _pas_present:
            _init_idx = _initial_state_index(sn, state_space_hashed)
            if _init_idx is not None and _init_idx >= 0:
                _reach = _forward_reachable(Q, _init_idx)
                if 0 < len(_reach) < Q.shape[0]:
                    _reach = np.asarray(sorted(int(x) for x in _reach), dtype=int)
                    Q = Q[np.ix_(_reach, _reach)]
                    state_space = state_space[_reach]
                    state_space_aggr = state_space_aggr[_reach]
                    state_space_hashed = state_space_hashed[_reach]
                    arvRates = arvRates[_reach]
                    depRates = depRates[_reach]
                    Dfilt = [d[np.ix_(_reach, _reach)] for d in Dfilt]
            elif Q.shape[0] > 1 and len(_forward_reachable(Q, 0)) < Q.shape[0]:
                # a reducible closed-PAS generator with no valid initial placement is a required-input error, not something to average over the mirror recurrent components.
                raise RuntimeError(
                    "A closed pass-and-swap network with a non-empty swapping graph requires "
                    "an explicit initial job placement. Call setState on the PAS station with "
                    "the ordered class list (oldest first) before solving.")

        # Solve for steady-state distribution
        pi = ctmc_solve(Q)

        # Compute metrics from sync-based data
        metrics = _compute_metrics_sync(
            sn, pi, arvRates, depRates, state_space_aggr,
            state_space, state_space_hashed, options)

        QN = metrics['Q']
        UN = metrics['U']
        RN = metrics['R']
        TN = metrics['T']

        CN = np.sum(RN, axis=0).reshape(1, -1)
        XN = np.zeros((1, K))
        for k in range(K):
            ref_stat = int(sn.refstat[k]) if hasattr(sn, 'refstat') and k < len(sn.refstat) else 0
            if ref_stat < M:
                XN[0, k] = TN[ref_stat, k]

        QN = np.nan_to_num(QN, nan=0.0)
        UN = np.nan_to_num(UN, nan=0.0)
        RN = np.nan_to_num(RN, nan=0.0)
        TN = np.nan_to_num(TN, nan=0.0)
        CN = np.nan_to_num(CN, nan=0.0)
        XN = np.nan_to_num(XN, nan=0.0)

        result = SolverCTMCReturn()
        result.Q = QN
        result.U = UN
        result.R = RN
        result.T = TN
        result.C = CN
        result.X = XN
        result.pi = pi
        result.infgen = Q
        result.space = state_space
        result.space_aggr = state_space_aggr
        result.arvRates = arvRates
        result.depRates = depRates
        result.space_hashed = state_space_hashed
        result.sn = sn
        # per-station column ranges in the concatenated sync state vector, needed by getProb*/getProbAggr.
        _scr = [(0, 0)] * M
        _col = 0
        for _isf in range(int(sn.nstateful)):
            _w = 0
            if sn.space is not None and _isf in sn.space and sn.space[_isf] is not None:
                _w = int(np.atleast_2d(sn.space[_isf]).shape[1])
            _ind = -1
            for _n in range(int(sn.nnodes)):
                if int(sn.nodeToStateful[_n]) == _isf:
                    _ind = _n
                    break
            _ist = int(sn.nodeToStation[_ind]) if _ind >= 0 else -1
            if 0 <= _ist < M:
                _scr[_ist] = (_col, _col + _w)
            _col += _w
        result.station_col_ranges = _scr
        result.depRates = depRates
        result.rrobin_info = {}
        result.eventFilt = Dfilt
        result.runtime = time.time() - start_time
        result.method = "sync"

        return result



def solver_ctmc(
    sn: NetworkStruct,
    options: Optional[SolverCTMCOptions] = None
) -> SolverCTMCReturn:
    """
    Main CTMC solver handler.

    Routes to appropriate method based on options and network characteristics.

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        SolverCTMCReturn with performance metrics
    """
    if options is None:
        options = SolverCTMCOptions()

    method = options.method.lower()

    if method in ['default', 'basic']:
        return solver_ctmc_basic(sn, options)
    else:
        # Unknown method - use basic
        if options.verbose:
            print(f"Warning: Unknown CTMC method '{method}'. Using basic.")
        return solver_ctmc_basic(sn, options)


__all__ = [
    'solver_ctmc',
    'solver_ctmc_basic',
    'SolverCTMCReturn',
    'SolverCTMCOptions',
]
