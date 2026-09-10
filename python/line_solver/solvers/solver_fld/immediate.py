"""
Elimination of immediate transitions from the fluid ODE system.

An immediate transition carries the rate GlobalConstants.Immediate, which is
large but finite, so leaving it in the drift makes the system stiff without
adding any dynamics: the mass it moves is gone within O(1/Immediate) of
arriving. Folding it away by STOCHASTIC COMPLEMENTATION removes the stiff
coordinates while preserving the stationary law of the ones that remain.

Two representations need it, and they are two functions rather than one because
the fluid solver builds the drift in two different forms:

- `ode_eliminate_immediate` works on the EVENT representation (all_jumps,
  rateBase, eventIdx) that the closing family integrates. It is handed rates it
  knows are `rateBase`, so it recognises an immediate transition by comparing
  against the constant itself.
- `eliminate_immediate_matrix` works on the GENERATOR W that the matrix method
  builds, where the same information is a row maximum.

Port from:
    - matlab/src/solvers/FLD/ode_eliminate_immediate.m
    - matlab/src/solvers/FLD/eliminate_immediate_matrix.m
    - matlab/src/solvers/FLD/generator_to_jumps.m
The JAR twin is jline.solvers.fluid.handlers.ImmediateElimination; the C++ twin
is fluid_eliminate_immediate / fluid_eliminate_immediate_matrix in
cpp/include/line/solvers/fluid/fluid_stiff.h.
"""

import warnings
from typing import Any, Optional, Tuple

import numpy as np

from ...constants import GlobalConstants

__all__ = [
    'fluid_hide_immediate',
    'generator_to_jumps',
    'ode_eliminate_immediate',
    'eliminate_immediate_matrix',
    'ode_expand_state',
]


def fluid_hide_immediate(sn: Any, options: Any) -> bool:
    """Whether this model's fluid drift is built on the stochastic complement.

    Every fluid route that builds its drift from the station/class/phase event set
    or from the linear generator asks here rather than reading the flag directly,
    so the answer is the same across matrix, closing, statedep, tbi, minnormal,
    refined and dae. The flag defaults to True for SolverFLD: a coordinate whose
    exit rate is ``GlobalConstants.Immediate`` is LINE's stand-in for infinity, and
    integrating it is meaningless work no integrator does well.

    THE STOCHASTIC PETRI NET ROUTE IS THE ONE EXCEPTION, and it is not a refusal.
    It carries immediate firings as ALGEBRAIC unknowns of an index-1 DAE, a
    stronger treatment than absorbing them, and never builds the event set this
    reduction acts on, so the answer here is simply False.
    """
    cfg = getattr(options, 'config', None) or {}
    flag = bool(getattr(options, 'hide_immediate', False))
    if isinstance(cfg, dict) and cfg.get('hide_immediate') is not None:
        flag = bool(cfg.get('hide_immediate'))
    if not flag:
        return False
    nodetype = getattr(sn, 'nodetype', None)
    if nodetype is not None:
        try:
            from ...constants import NodeType
            if any(int(nt) == int(NodeType.Transition) for nt in np.asarray(nodetype).ravel()):
                return False
        except Exception:  # noqa: BLE001 - a model without node types is not an SPN
            pass
    return True


def _immediate_tolerance_rates() -> float:
    """A transition is immediate when its rate is within 1% of the constant."""
    return GlobalConstants.Immediate * (1.0 - 0.01)


def _immediate_tolerance_matrix(options: Any) -> float:
    """Same rule as the event form: only the InfRate sentinel qualifies."""
    cfg = getattr(options, 'config', None) or {}
    if isinstance(cfg, dict) and cfg.get('immediate_tol') is not None:
        return float(cfg['immediate_tol'])
    tol = getattr(options, 'immediate_tol', None)
    if tol is not None:
        return float(tol)
    return GlobalConstants.Immediate * (1.0 - 0.01)


def generator_to_jumps(W: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Infinitesimal generator -> (all_jumps, rateBase, eventIdx).

    A jump vector holds -1 at the source coordinate and +1 at the destination,
    which is the form the closing drift multiplies its rates by.
    """
    W = np.asarray(W, dtype=float)
    n_states = W.shape[0]
    src, dst = np.nonzero(W)
    rates = W[src, dst]
    keep = (src != dst) & (rates > 0)
    src, dst, rates = src[keep], dst[keep], rates[keep]

    all_jumps = np.zeros((n_states, len(rates)))
    for t in range(len(rates)):
        all_jumps[src[t], t] = -1.0
        all_jumps[dst[t], t] = 1.0
    return all_jumps, rates.copy(), src.copy()


def _stochcomp(W: np.ndarray, timed_states: np.ndarray):
    """ctmc_stochcomp on the kept states, returning (S, Q22)."""
    from ...api.mc.ctmc import ctmc_stochcomp
    out = ctmc_stochcomp(np.asarray(W, dtype=float), np.asarray(timed_states, dtype=int))
    return out['S'], out['Q22']


def ode_eliminate_immediate(all_jumps, rateBase, eventIdx, sn, options):
    """Stochastic-complement the INSTANTANEOUS coordinates out of the event set.

    Returns ``(all_jumps_red, rateBase_red, eventIdx_red, state_map, Emap,
    absorb)``. The reduced jump matrix keeps the ORIGINAL coordinate layout, with
    the eliminated rows identically zero, so a caller that indexes the drift by
    the original coordinates keeps working; ``state_map`` names the coordinates
    that survived.

    THE REDUCTION IS EXACT, not an approximation. A coordinate whose exit rate is
    ``GlobalConstants.Immediate`` is LINE's stand-in for infinity, written by
    SolverLN for the branch of an activity that takes no time, and the flow that
    would enter it is routed straight to where it would have sent it.

    WHY THIS IS A STRUCTURAL COMPOSITION AND NOT A GENERATOR ROUND TRIP. Every
    event is a single -1 at ``eventIdx`` and a single +1 at its destination, so a
    path through the immediate block composes to one event that keeps the original
    source's GATING. Rebuilding the events from a reduced generator instead loses
    that identity, and with it the event ORDER the moment-closure methods read
    throughputs off -- which is why they used to refuse the reduction outright.
    ``Emap[e, o]`` is the expected number of times the ORIGINAL event ``o`` fires
    per firing of the reduced event ``e``, so a caller maps any per-event quantity
    with ``new = Emap @ old``. It is the identity when nothing is eliminated.

    A COMPOSED EVENT CAN BE A DEPARTURE AT TWO STATIONS AT ONCE: a job that leaves
    the delay, passes through a queue's immediate phase and returns has completed
    at both, and both throughputs must count it.

    ``absorb`` projects an initial condition onto the surviving coordinates. Mass
    parked on an eliminated one would otherwise be frozen there for the whole
    integration, because nothing moves it any more.

    Every failure path returns the system UNCHANGED rather than a partially
    reduced one: an elimination that cannot be carried out is not a reason to
    integrate a different model than the caller built.
    """
    all_jumps = np.asarray(all_jumps, dtype=float)
    rateBase = np.asarray(rateBase, dtype=float).ravel()
    eventIdx = np.asarray(eventIdx, dtype=int).ravel()
    n_states = all_jumps.shape[0]
    n_events = rateBase.size

    def _identity():
        return (all_jumps, rateBase, eventIdx, np.arange(n_states, dtype=int),
                np.eye(n_events), np.eye(n_states))

    imm_tol = _immediate_tolerance_rates()
    cfg = getattr(options, 'config', None) or {}
    if isinstance(cfg, dict) and cfg.get('immediate_tol') is not None:
        imm_tol = float(cfg['immediate_tol'])

    imm_idx = np.nonzero(rateBase >= imm_tol)[0]
    if imm_idx.size == 0:
        return _identity()

    try:
        return _eliminate_structural(all_jumps, rateBase, eventIdx, imm_idx)
    except Exception as exc:  # noqa: BLE001 - the reference falls back on any failure
        warnings.warn(
            "Immediate coordinate elimination failed: %s. Using original system." % exc)
        return _identity()


def _eliminate_structural(all_jumps, rateBase, eventIdx, imm_idx):
    n_states = all_jumps.shape[0]
    n_events = rateBase.size

    def _identity():
        return (all_jumps, rateBase, eventIdx, np.arange(n_states, dtype=int),
                np.eye(n_events), np.eye(n_states))

    # The immediate coordinates are the SOURCES of the immediate events: it is the
    # coordinate that empties instantaneously, not the event.
    is_imm = np.zeros(n_states, dtype=bool)
    is_imm[eventIdx[imm_idx]] = True

    # Destination of each event. A departure that re-enters its own coordinate
    # cancels to an all-zero column, whose destination is that same coordinate.
    src = eventIdx.astype(int)
    dst = src.copy()
    for e in range(n_events):
        pos = np.nonzero(all_jumps[:, e] > 0)[0]
        if pos.size:
            dst[e] = int(pos[0])

    # A coordinate with no outflow cannot be complemented away, and one whose
    # outflow is entirely a self-loop would make the fundamental matrix singular.
    for f in np.nonzero(is_imm)[0]:
        out_f = np.nonzero(src == f)[0]
        if rateBase[out_f].sum() <= 0 or np.all(dst[out_f] == f):
            is_imm[f] = False

    Fidx = np.nonzero(is_imm)[0]
    Sidx = np.nonzero(~is_imm)[0]
    if Fidx.size == 0 or Sidx.size == 0:
        return _identity()
    nF, nS = Fidx.size, Sidx.size
    posF = np.zeros(n_states, dtype=int)
    posS = np.zeros(n_states, dtype=int)
    posF[Fidx] = np.arange(nF)
    posS[Sidx] = np.arange(nS)

    # Branching of the embedded jump chain out of each immediate coordinate. The
    # probabilities are the rate shares, so a coordinate carrying both an immediate
    # and an ordinary exit gives the ordinary one its (vanishing) share.
    PFF = np.zeros((nF, nF))
    PFS = np.zeros((nF, nS))
    cnt = np.zeros((nF, n_events))
    for a in range(nF):
        f = Fidx[a]
        out_f = np.nonzero(src == f)[0]
        tot = rateBase[out_f].sum()
        for e in out_f:
            p = rateBase[e] / tot
            cnt[a, e] += p
            if is_imm[dst[e]]:
                PFF[a, posF[dst[e]]] += p
            else:
                PFS[a, posS[dst[e]]] += p

    # Fundamental matrix of the instantaneous chain. A closed cycle of immediate
    # coordinates has no absorption distribution and is left unreduced.
    Nfm = np.linalg.solve(np.eye(nF) - PFF, np.eye(nF))
    Aabs = Nfm @ PFS
    exp_cnt = Nfm @ cnt
    if not np.all(np.isfinite(Aabs)) or np.any(Aabs.sum(axis=1) < 0.5):
        warnings.warn("the immediate coordinates form a closed cycle, so they have no "
                      "absorption distribution; integrating the unreduced system instead")
        return _identity()

    jumps_new = []
    rate_new = []
    evidx_new = []
    emap_rows = []
    for e in range(n_events):
        if is_imm[src[e]]:
            continue  # its flow is already carried by whichever event feeds it
        if not is_imm[dst[e]]:
            jumps_new.append(all_jumps[:, e].copy())
            rate_new.append(rateBase[e])
            evidx_new.append(src[e])
            row = np.zeros(n_events)
            row[e] = 1.0
            emap_rows.append(row)
            continue
        # The event feeds an immediate coordinate: one event per absorbing
        # destination, keeping the original source and so the original gating,
        # since the rate of the composed flow IS the rate of the inflow.
        a = posF[dst[e]]
        for b in np.nonzero(Aabs[a] > 0)[0]:
            jump = np.zeros(n_states)
            jump[src[e]] -= 1.0
            jump[Sidx[b]] += 1.0
            jumps_new.append(jump)
            rate_new.append(rateBase[e] * Aabs[a, b])
            evidx_new.append(src[e])
            # Weighting every absorbing branch by the SAME unconditional expected
            # counts is what makes the rate accounting exact: the branch rates sum
            # back to rateBase[e], so the mapped total is rateBase[e] * counts.
            row = exp_cnt[a].copy()
            row[e] += 1.0
            emap_rows.append(row)

    all_jumps_red = np.column_stack(jumps_new) if jumps_new else np.zeros((n_states, 0))
    rateBase_red = np.asarray(rate_new, dtype=float)
    eventIdx_red = np.asarray(evidx_new, dtype=int)
    Emap = np.vstack(emap_rows) if emap_rows else np.zeros((0, n_events))

    absorb = np.eye(n_states)
    absorb[Fidx, :] = 0.0
    for a in range(nF):
        absorb[Fidx[a], Sidx] = Aabs[a]

    return all_jumps_red, rateBase_red, eventIdx_red, Sidx.astype(int), Emap, absorb


def eliminate_immediate_matrix(W, sn, options) -> Tuple[np.ndarray, np.ndarray]:
    """The same elimination on the generator the matrix method builds.

    Returns (W_red, state_map). A state is immediate when its largest outgoing
    rate reaches the threshold, which is the generator-form reading of the same
    test the event form makes against rateBase.
    """
    W = np.asarray(W, dtype=float)
    n = W.shape[0]
    identity_map = np.arange(n, dtype=int)

    imm_tol = _immediate_tolerance_matrix(options)
    imm_states = np.nonzero(np.abs(W).max(axis=1) >= imm_tol)[0]
    if imm_states.size == 0:
        return W, identity_map

    timed_states = np.setdiff1d(np.arange(n), imm_states)
    if timed_states.size <= 1:
        return W, identity_map

    try:
        W_red, _ = _stochcomp(W, timed_states)
        if not np.all(np.isfinite(W_red)):
            warnings.warn(
                "Stochastic complementation produced non-finite values. "
                "Using original system.")
            return W, identity_map
        return W_red, timed_states.astype(int)
    except Exception as exc:  # noqa: BLE001 - the reference falls back on any failure
        warnings.warn(
            "Immediate transition elimination failed: %s. Using original system." % exc)
        return W, identity_map


def ode_expand_state(x_reduced, state_map, n_original) -> np.ndarray:
    """Reduced state vector back to the original dimension, eliminated slots at 0."""
    x_full = np.zeros(int(n_original))
    x_full[np.asarray(state_map, dtype=int)] = np.asarray(x_reduced, dtype=float).ravel()
    return x_full
