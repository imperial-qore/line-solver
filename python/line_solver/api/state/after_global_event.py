"""
SPN global event handler for CTMC state-space generation.

Mirrors MATLAB matlab/src/lang/+State/afterGlobalEvent.m and JAR
jline/lang/state/AfterGlobalEvent.{kt,java}. Used by the CTMC sync builder to
process ENABLE and FIRE events at Transition nodes.

State vector at a Transition (per row of sn.space[isf_transition]):

    [ idle_m1, phase_counts_m1(fK_1),
      idle_m2, phase_counts_m2(fK_2),
      ...
      idle_mM, phase_counts_mM(fK_M),
      fired_m1, fired_m2, ..., fired_mM,
      <local_vars> ]

For CTMC the fired counts are not used (initialized to zero); they are kept
for parity with the simulation path. ``nmodeservers[m]`` may be infinite
(modeled as MaxInt in the state vector).

State vectors at Places hold per-class job counts in the first ``nclasses``
columns (count-based layout matches the Place's row in ``sn.space``).
"""

from __future__ import annotations

import numpy as np
from typing import List, Optional, Sequence, Tuple

from ...constants import EventType, GlobalConstants


def _split_transition_state(sn, ind: int, state_row: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int, int, np.ndarray]:
    """
    Slice a Transition's row into its components.

    Returns:
        (idle, phase, fired, var, nmodes, V, fK)
    """
    nparam = sn.nodeparam[ind]
    nmodes = int(getattr(nparam, 'nmodes', 0) if not isinstance(nparam, dict)
                 else nparam.get('nmodes', 0))
    fphases = np.asarray(getattr(nparam, 'firingphases', np.array([])), dtype=float)

    fK = np.zeros(nmodes, dtype=int)
    fp = getattr(nparam, 'firingproc', None)
    for m in range(nmodes):
        if m < fphases.size and not np.isnan(fphases[m]):
            fK[m] = int(fphases[m])
        elif fp is not None and m < len(fp) and fp[m] is not None:
            D0 = np.atleast_2d(np.asarray(fp[m][0]))
            fK[m] = int(D0.shape[0])
        else:
            fK[m] = 1
    fK = np.maximum(fK, 1)

    sumK = int(np.sum(fK))
    V = int(np.sum(sn.nvars[ind])) if (sn.nvars is not None and ind < sn.nvars.shape[0]) else 0

    state_row = np.asarray(state_row, dtype=float).ravel()
    expected = nmodes + sumK + nmodes + V
    if state_row.size < expected:
        # Allow legacy layout without fired block (older MATLAB exports).
        legacy = nmodes + sumK + V
        if state_row.size == legacy:
            idle = state_row[:nmodes].astype(int)
            phase = state_row[nmodes:nmodes + sumK].astype(int)
            fired = np.zeros(nmodes, dtype=int)
            var = state_row[nmodes + sumK:].astype(int)
            return idle, phase, fired, var, nmodes, V, fK
        raise ValueError(
            f"Transition state row of length {state_row.size} too short "
            f"(expected {expected} for nmodes={nmodes}, sum(fK)={sumK}, V={V})"
        )

    idle = state_row[:nmodes].astype(int)
    phase = state_row[nmodes:nmodes + sumK].astype(int)
    fired = state_row[nmodes + sumK:nmodes + sumK + nmodes].astype(int)
    var = state_row[nmodes + sumK + nmodes:].astype(int)
    return idle, phase, fired, var, nmodes, V, fK


def _join_transition_state(idle: np.ndarray, phase: np.ndarray, fired: np.ndarray, var: np.ndarray) -> np.ndarray:
    """Recompose a Transition state row from components."""
    return np.concatenate([
        idle.astype(float),
        phase.astype(float),
        fired.astype(float),
        var.astype(float),
    ])


def _phase_slice(fK: np.ndarray, mode: int) -> Tuple[int, int]:
    """Return (start, end) index pair into the per-mode phase block."""
    start = int(np.sum(fK[:mode]))
    return start, start + int(fK[mode])


def _enabling_count_at_place(glspace: Sequence[np.ndarray], isf: int, job_class: int) -> int:
    """Read the number of class `job_class` tokens at the place row glspace[isf]."""
    if isf is None or isf < 0 or isf >= len(glspace):
        return 0
    row = np.asarray(glspace[isf], dtype=float).ravel()
    if job_class < 0 or job_class >= row.size:
        return 0
    return int(row[job_class])


def _multichoose_capped(n_bins: int, total: int, caps: Optional[np.ndarray] = None) -> np.ndarray:
    """Enumerate length-n vectors of non-negatives summing to ``total`` with
    optional per-bin upper bounds ``caps`` (used only by the running>en case)."""
    if n_bins <= 0:
        return np.zeros((1 if total == 0 else 0, 0), dtype=int)
    if n_bins == 1:
        if caps is not None and caps[0] < total:
            return np.zeros((0, 1), dtype=int)
        return np.array([[total]], dtype=int)
    rows = []
    upper = total if caps is None else min(total, int(caps[0]))
    for i in range(upper + 1):
        sub_caps = caps[1:] if caps is not None else None
        sub = _multichoose_capped(n_bins - 1, total - i, sub_caps)
        if sub.shape[0] == 0:
            continue
        first = np.full((sub.shape[0], 1), i, dtype=int)
        rows.append(np.hstack([first, sub]))
    if not rows:
        return np.zeros((0, n_bins), dtype=int)
    return np.vstack(rows)


def _factln(n: int) -> float:
    if n <= 1:
        return 0.0
    return float(np.sum(np.log(np.arange(2, n + 1))))


def _multinomial_log_prob(comb: np.ndarray, pentry: np.ndarray) -> float:
    """logP for the multinomial distribution P(comb | sum(comb), pentry).

    Returns -inf if any pentry[k]==0 yet comb[k]>0.
    """
    n = int(np.sum(comb))
    logp = _factln(n)
    for k in range(comb.size):
        if pentry[k] > 0:
            logp += comb[k] * np.log(pentry[k]) - _factln(int(comb[k]))
        elif pentry[k] == 0 and comb[k] == 0:
            continue
        else:
            return -np.inf
    return logp


def after_global_event(sn, ind: int, glspace: List[np.ndarray], glevent,
                       is_simulation: bool = False):
    """
    Process an SPN global synchronization event.

    Args:
        sn: NetworkStruct (with nodeparam[ind] populated incl. firingproc /
            firingphases / firingpie).
        ind: Node index of the active Transition.
        glspace: list of state-row vectors, indexed by isf, one row per
            stateful node.
        glevent: GlobalSync entry from sn.gsync (active=[ENABLE|FIRE], passive=[...]).
        is_simulation: when True, returns a single sampled outcome (for SSA).
            CTMC always sets this False.

    Returns:
        Object with attributes:
            outglspace : List[List[np.ndarray]]
                Per-outcome list of glspace lists. outglspace[i][isf] is the
                row for stateful node isf in outcome i.
            outrate    : np.ndarray, shape (n_outcomes,)
            outprob    : np.ndarray, shape (n_outcomes,)
            is_completion : np.ndarray of bool (True iff the outcome triggered
                D1 firing — informs depRates/arvRates accounting on the caller).
    """
    R = int(sn.nclasses)

    nparam = sn.nodeparam[ind]
    nmodes = int(getattr(nparam, 'nmodes', 0) if not isinstance(nparam, dict)
                 else nparam.get('nmodes', 0))
    if nmodes <= 0:
        return _empty_result()

    isf_transition = int(sn.nodeToStateful[ind])
    trans_row = np.asarray(glspace[isf_transition], dtype=float).ravel()
    idle, phase, fired, var, nmodes, V, fK = _split_transition_state(sn, ind, trans_row)

    active = glevent.active[0]
    mode = int(active.mode)
    if mode < 0 or mode >= nmodes:
        return _empty_result()

    fp = getattr(nparam, 'firingproc', None)
    fpie_list = getattr(nparam, 'firingpie', None)
    nmodeservers = np.asarray(getattr(nparam, 'nmodeservers', np.ones(nmodes)), dtype=float)

    if active.event == EventType.ENABLE:
        return _handle_enable(sn, ind, isf_transition, glspace, glevent, mode,
                              idle, phase, fired, var, nmodes, fK,
                              nmodeservers, fp, fpie_list, R, is_simulation)
    elif active.event == EventType.FIRE:
        return _handle_fire(sn, ind, isf_transition, glspace, glevent, mode,
                            idle, phase, fired, var, nmodes, fK,
                            nmodeservers, fp, R, is_simulation)
    else:
        return _empty_result()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

class _GlobalEventResult:
    __slots__ = ('outglspace', 'outrate', 'outprob', 'is_completion')

    def __init__(self, outglspace, outrate, outprob, is_completion):
        self.outglspace = outglspace
        self.outrate = np.asarray(outrate, dtype=float)
        self.outprob = np.asarray(outprob, dtype=float)
        self.is_completion = np.asarray(is_completion, dtype=bool)


def _empty_result() -> _GlobalEventResult:
    return _GlobalEventResult([], np.zeros(0), np.zeros(0), np.zeros(0, dtype=bool))


def _enabling_degree_for_mode(sn, ind: int, mode: int, glspace: Sequence[np.ndarray],
                              glevent, R: int) -> int:
    """Compute the maximum number of times mode m is concurrently enabled by
    the current token marking at the input places (passive entries)."""
    nparam = sn.nodeparam[ind]
    enabling_list = getattr(nparam, 'enabling', None) if not isinstance(nparam, dict) else nparam.get('enabling', None)
    if enabling_list is None or mode >= len(enabling_list):
        return 0
    en = np.atleast_2d(np.asarray(enabling_list[mode], dtype=float))
    inhibiting_list = getattr(nparam, 'inhibiting', None) if not isinstance(nparam, dict) else nparam.get('inhibiting', None)

    # ep_space[ind, r] = current token count at (place=ind, class=r).
    ep_space = np.zeros((sn.nnodes, R), dtype=float)
    for pev in glevent.passive:
        ep_ind = pev.node
        if ep_ind >= sn.nnodes:
            continue
        if not sn.isstateful[ep_ind]:
            continue
        isf = int(sn.nodeToStateful[ep_ind])
        row = np.asarray(glspace[isf], dtype=float).ravel()
        # A Place holds its tokens in the [buffer(R), server(R)] slots (arrivals
        # deposit into the server slot, a firing relocates the remainder into the
        # buffer slot). Count both, matching toMarginalAggr; reading only the
        # first R columns saw the (empty) buffer slot and left every mode disabled.
        if row.size >= 2 * R:
            ep_space[ep_ind, :R] = row[:R] + row[R:2 * R]
        elif row.size >= R:
            ep_space[ep_ind, :R] = row[:R]
        else:
            ep_space[ep_ind, :R] = np.pad(row, (0, R - row.size))

    # Inhibitor arcs: mode cannot fire while any inhibited input place has
    # reached its threshold (Inf default => never true).
    if inhibiting_list is not None and mode < len(inhibiting_list):
        inh = np.atleast_2d(np.asarray(inhibiting_list[mode], dtype=float))
        if np.any((inh < np.inf) & (ep_space >= inh)):
            return 0

    en_degree = 1
    while np.all(ep_space >= en_degree * en):
        en_degree += 1
    return en_degree - 1


def _marking_matrix(sn, glspace, glevent, R):
    """Node-indexed input-place marking (sn.nnodes x R) over the mode's passive
    places, the argument passed to a firing-rate dependence handle. Mirrors the
    ep_space construction in _enabling_degree_for_mode."""
    m = np.zeros((sn.nnodes, R), dtype=float)
    for pev in glevent.passive:
        ep_ind = pev.node
        if ep_ind >= sn.nnodes or not sn.isstateful[ep_ind]:
            continue
        isf = int(sn.nodeToStateful[ep_ind])
        row = np.asarray(glspace[isf], dtype=float).ravel()
        if row.size >= 2 * R:
            m[ep_ind, :R] = row[:R] + row[R:2 * R]
        elif row.size >= R:
            m[ep_ind, :R] = row[:R]
        else:
            m[ep_ind, :R] = np.pad(row, (0, R - row.size))
    return m


def _handle_enable(sn, ind, isf_transition, glspace, glevent, mode,
                   idle, phase, fired, var, nmodes, fK, nmodeservers,
                   fp, fpie_list, R, is_simulation) -> _GlobalEventResult:
    """ENABLE event: adjust idle/phase counts so that running servers in mode
    `mode` matches its enabling degree."""
    en_degree = _enabling_degree_for_mode(sn, ind, mode, glspace, glevent, R)

    # Cap by per-mode server count (treating Inf as MaxInt).
    nmsv = nmodeservers[mode]
    if not np.isfinite(nmsv):
        nmsv = int(GlobalConstants.MaxInt)
    en_degree = min(en_degree, int(nmsv))

    p_start, p_end = _phase_slice(fK, mode)
    running_m = int(np.sum(phase[p_start:p_end]))

    # Disable: tokens at PRE places are below the enabling threshold.
    nparam = sn.nodeparam[ind]
    enabling_list = getattr(nparam, 'enabling', None) if not isinstance(nparam, dict) else nparam.get('enabling', None)
    inhibiting_list = getattr(nparam, 'inhibiting', None) if not isinstance(nparam, dict) else nparam.get('inhibiting', None)
    en = np.atleast_2d(np.asarray(enabling_list[mode], dtype=float))
    ep_space = np.zeros((sn.nnodes, R), dtype=float)
    for pev in glevent.passive:
        ep_ind = pev.node
        if ep_ind >= sn.nnodes or not sn.isstateful[ep_ind]:
            continue
        isf = int(sn.nodeToStateful[ep_ind])
        row = np.asarray(glspace[isf], dtype=float).ravel()
        # A Place holds its tokens in the [buffer(R), server(R)] slots (arrivals
        # deposit into the server slot, a firing relocates the remainder into the
        # buffer slot). Count both, matching toMarginalAggr; reading only the
        # first R columns saw the (empty) buffer slot and left every mode disabled.
        if row.size >= 2 * R:
            ep_space[ep_ind, :R] = row[:R] + row[R:2 * R]
        elif row.size >= R:
            ep_space[ep_ind, :R] = row[:R]
        else:
            ep_space[ep_ind, :R] = np.pad(row, (0, R - row.size))

    # Inhibitor arcs: mode is disabled while any inhibited input place has
    # reached its threshold (Inf default => never true).
    inhibited = False
    if inhibiting_list is not None and mode < len(inhibiting_list):
        inh = np.atleast_2d(np.asarray(inhibiting_list[mode], dtype=float))
        inhibited = bool(np.any((inh < np.inf) & (ep_space >= inh)))

    if np.any(ep_space < en) or inhibited:
        # Disable: park all mode-m servers in idle pool, zero phase counters.
        new_idle = idle.copy()
        new_phase = phase.copy()
        new_idle[mode] = int(nmsv)
        new_phase[p_start:p_end] = 0
        new_state = _join_transition_state(new_idle, new_phase, fired, var)
        if np.array_equal(new_state, _join_transition_state(idle, phase, fired, var)):
            # Already disabled — emit no-op outcome (rate 0, prob 1).
            return _emit_single(glspace, isf_transition, _join_transition_state(idle, phase, fired, var),
                                rate=0.0, prob=1.0, is_completion=False)
        return _emit_single(glspace, isf_transition, new_state,
                            rate=GlobalConstants.Immediate if hasattr(GlobalConstants, 'Immediate') else 1.0 / GlobalConstants.FineTol,
                            prob=1.0, is_completion=False)

    if running_m == en_degree:
        return _emit_single(glspace, isf_transition,
                            _join_transition_state(idle, phase, fired, var),
                            rate=0.0, prob=1.0, is_completion=False)

    # Recover entry probabilities.
    fK_m = int(fK[mode])
    if fpie_list is not None and mode < len(fpie_list) and fpie_list[mode] is not None:
        pentry = np.atleast_1d(np.asarray(fpie_list[mode], dtype=float)).ravel()
        if pentry.size < fK_m:
            pentry = np.pad(pentry, (0, fK_m - pentry.size))
    else:
        pentry = np.zeros(fK_m, dtype=float)
        if fK_m > 0:
            pentry[0] = 1.0

    immediate_rate = GlobalConstants.Immediate if hasattr(GlobalConstants, 'Immediate') else 1.0 / GlobalConstants.FineTol

    if running_m < en_degree:
        # Activate (en_degree - running_m) extra servers, distributing phase entries
        # multinomially with weights ``pentry``.
        nadd = en_degree - running_m
        combs = _multichoose_capped(fK_m, nadd, None)
        # Sort descending to match MATLAB's sortrows(descend) ordering.
        combs = combs[np.lexsort(combs.T[::-1])][::-1]

        outcomes = []
        for c in range(combs.shape[0]):
            comb = combs[c]
            new_phase = phase.copy()
            for k in range(fK_m):
                new_phase[p_start + k] += int(comb[k])
            new_idle = idle.copy()
            new_idle[mode] = max(0, int(idle[mode]) - nadd)
            new_state = _join_transition_state(new_idle, new_phase, fired, var)
            logp = _multinomial_log_prob(comb, pentry)
            if not np.isfinite(logp):
                continue
            outcomes.append((new_state, immediate_rate, float(np.exp(logp)), False, None))
        return _emit_many(glspace, isf_transition, outcomes)

    # running_m > en_degree: stop ndiff servers chosen via hypergeometric weights.
    ndiff = running_m - en_degree
    srv_vec = phase[p_start:p_end].astype(int)
    n = int(srv_vec.size)
    all_combs = _multichoose_capped(n, ndiff, srv_vec.copy())
    if all_combs.shape[0] == 0:
        return _emit_single(glspace, isf_transition, _join_transition_state(idle, phase, fired, var),
                            rate=0.0, prob=1.0, is_completion=False)
    logW = np.zeros(all_combs.shape[0])
    for i in range(all_combs.shape[0]):
        for k in range(n):
            logW[i] += _factln(int(srv_vec[k])) - _factln(int(all_combs[i, k])) - _factln(int(srv_vec[k] - all_combs[i, k]))
    W = np.exp(logW - np.max(logW))
    W = W / np.sum(W)

    outcomes = []
    for i in range(all_combs.shape[0]):
        comb = all_combs[i]
        new_phase = phase.copy()
        for k in range(n):
            new_phase[p_start + k] = int(srv_vec[k] - comb[k])
        new_idle = idle.copy()
        new_idle[mode] = int(idle[mode]) + ndiff
        new_state = _join_transition_state(new_idle, new_phase, fired, var)
        outcomes.append((new_state, immediate_rate, float(W[i]), False, None))
    return _emit_many(glspace, isf_transition, outcomes)


def _handle_fire(sn, ind, isf_transition, glspace, glevent, mode,
                 idle, phase, fired, var, nmodes, fK, nmodeservers, fp,
                 R, is_simulation) -> _GlobalEventResult:
    """FIRE event: enumerate D0 phase moves and D1 completions for the active mode."""
    # An IMMEDIATE mode fires in zero time and carries no firing process: add_mode
    # leaves the distribution None unless setDistribution is called, and Immediate
    # derives from Det, so nothing populates a (D0,D1) proc for it. Returning empty
    # here would make the mode unable to ever fire, leaving the marking that enables
    # it absorbing and the whole generator degenerate. Its rate does not come from a
    # firing process anyway (see below), so a single-phase placeholder is enough to
    # drive the enumeration; _split_transition_state already reserves one phase for
    # such a mode.
    nparam_f = sn.nodeparam[ind]
    is_immediate_mode = False
    timing_m = getattr(nparam_f, 'timingstrategies', None)
    if timing_m is not None and mode < len(timing_m):
        # timingstrategies holds a TimingStrategy enum once set, but defaults to
        # the string 'TIMED', so compare by name as the JMT handler does.
        is_immediate_mode = str(getattr(timing_m[mode], 'name', timing_m[mode])).upper() == 'IMMEDIATE'

    has_proc = not (fp is None or mode >= len(fp) or fp[mode] is None)
    if not has_proc and not is_immediate_mode:
        return _empty_result()
    if has_proc:
        D0 = np.atleast_2d(np.asarray(fp[mode][0], dtype=float))
        D1 = np.atleast_2d(np.asarray(fp[mode][1], dtype=float))
    else:
        _imm = GlobalConstants.Immediate if hasattr(GlobalConstants, 'Immediate') else 1.0 / GlobalConstants.FineTol
        D0 = np.array([[-_imm]], dtype=float)
        D1 = np.array([[_imm]], dtype=float)
    fK_m = int(fK[mode])
    p_start, p_end = _phase_slice(fK, mode)

    outcomes = []  # list of (new_state, rate, prob, is_completion, glspace_overrides)

    # NOTE: D0 off-diagonal phase moves are emitted by the PHASE sync action
    # (refresh_sync emits one PHASE event per Transition mode, which dispatches
    # to after_event_transition and enumerates the same D0 moves). Adding them
    # here too would double-count — e.g., for Erlang(k) closed loop with N=1,
    # the observed Tput would be 2k/(k+1) instead of the correct 1.0.

    # D1 completions: source phase k -> firing, plus PRE/POST place updates.
    immediate_rate = GlobalConstants.Immediate if hasattr(GlobalConstants, 'Immediate') else 1.0 / GlobalConstants.FineTol

    # An IMMEDIATE mode fires in zero time, so its firing process is not the
    # rate: the mode's distribution is ignored and the firing is emitted at
    # GlobalConstants.Immediate scaled by the mode's firing weight. Only the
    # ratio of these rates matters, since it sets the branching probabilities
    # among the immediate modes enabled together, and the states in which they
    # fire are vanishing and removed by stochastic complementation in the CTMC
    # handler. The Immediate scale makes any timed mode enabled in the same
    # marking lose the race, which is the priority of immediate over timed
    # transitions required by GSPN semantics.
    immediate_weight = 1.0
    if is_immediate_mode:
        fw = getattr(nparam_f, 'fireweight', None)
        if fw is not None and mode < len(fw):
            immediate_weight = float(fw[mode])

    # Latching costs no time for a timed mode, whose ENABLE runs at
    # GlobalConstants.Immediate while its firing runs at an ordinary rate. For an
    # immediate mode the two are the same scale, so a marking that enables two
    # immediate modes can fire the one that happens to be latched before the other
    # is latched at all, and the branching then follows the latching order instead
    # of the firing weights: with weights 1 and 3 the split came out at 0.344
    # rather than 0.25. Gating on the marking makes every latch variant of a
    # marking offer the same firing alternatives, which is the weight ratio GSPN
    # semantics require. mark_degree_m is how many concurrent firings the marking
    # alone supports, before any cap by the servers actually latched. Mirrors
    # MATLAB afterGlobalEvent.m:315-330.
    mark_degree_m = _enabling_degree_for_mode(sn, ind, mode, glspace, glevent, R)
    if is_immediate_mode:
        nmsv_m = float(nmodeservers[mode])
        if not np.isfinite(nmsv_m):
            nmsv_m = float(GlobalConstants.MaxInt)
        imm_servers_m = min(mark_degree_m, int(nmsv_m))
    else:
        imm_servers_m = 0

    # Marking-dependent firing-rate multiplier g_mode(marking): evaluate once at
    # the current node-indexed input-place marking. Exact per enumerated CTMC
    # state; the SSA serial handler samples from the same dependent rates.
    firing_mult = 1.0
    firingdep_list = getattr(nparam_f, 'firingdep', None)
    if (not is_immediate_mode and firingdep_list is not None
            and mode < len(firingdep_list) and firingdep_list[mode] is not None):
        firing_mult = float(firingdep_list[mode](_marking_matrix(sn, glspace, glevent, R)))

    for k in range(fK_m):
        cnt_k = int(phase[p_start + k])
        # A completion consumes the enabling weight from the input places (the PRE
        # below), so it can only fire while those places still hold it. A server can
        # be left latched in a marking that no longer enables it: two immediate
        # modes contend for one token, the first consumes it, and the second stays
        # in execution with cnt_k > 0 though its input place is now empty. Gating on
        # cnt_k alone then fired it and drove the place negative. Given cnt_k > 0 a
        # server is running in this mode, so mark_degree_m >= 1 is exactly MATLAB's
        # en_degree_m >= 1 = min(mark_degree_m, running_m) >= 1.
        if is_immediate_mode:
            if k != 0 or imm_servers_m < 1:
                continue
            rate_kd = immediate_rate * immediate_weight * imm_servers_m
        else:
            if cnt_k <= 0 or mark_degree_m < 1:
                continue
            rate_kd = float(np.sum(D1[k, :])) * cnt_k * firing_mult
        if rate_kd <= 0:
            continue

        # Update Transition row. Only a latched server can be retired: an
        # immediate mode gated on the marking may fire from a state where none
        # is latched yet. Mirrors MATLAB afterGlobalEvent.m:339.
        new_phase = phase.copy()
        new_idle = idle.copy()
        if cnt_k > 0:
            new_phase[p_start + k] -= 1
            new_idle[mode] += 1
        new_fired = fired.copy()
        if is_simulation:
            new_fired[mode] += 1
        trans_state = _join_transition_state(new_idle, new_phase, new_fired, var)

        # Apply PRE/POST place edits to a fresh glspace copy. One D1 completion
        # is ONE firing: PRE consumes ``weight`` and POST produces ``weight``.
        place_overrides = {isf_transition: trans_state}
        if not _apply_pre_post(sn, glevent, glspace, place_overrides):
            # PRE consume failed (insufficient tokens) — skip.
            continue

        outcomes.append((trans_state, rate_kd, 1.0, True, place_overrides))

    return _emit_many(glspace, isf_transition, outcomes)


def _apply_pre_post(sn, glevent, glspace: Sequence[np.ndarray],
                    overrides: dict) -> bool:
    """Apply PRE consume and POST produce edits to overrides[isf] = new_row.
    Returns False if any PRE place lacks the required tokens."""
    R = int(sn.nclasses)

    # PRE: consume `weight` tokens from the (place, class) entry. The concurrency
    # of the other running servers is already carried by the completion rate, so
    # scaling by the enabling degree here would let a single firing consume every
    # enabled token and would not conserve tokens against POST.
    for pev in glevent.passive:
        if pev.event != EventType.PRE:
            continue
        ep_ind = pev.node
        if ep_ind >= sn.nnodes or not sn.isstateful[ep_ind]:
            continue
        isf = int(sn.nodeToStateful[ep_ind])
        row = np.array(overrides.get(isf, glspace[isf]), dtype=float).ravel().copy()
        ep_class = int(pev.job_class)
        consume = int(pev.weight)
        if row.size >= 2 * R:
            # [buffer(R), server(R)] Place: consume from the total and keep the
            # surviving tokens in the buffer slot (mirrors MATLAB
            # afterGlobalEvent.m), so the folded marginal (buf+srv) drops by the
            # arc weight without producing negative counts.
            total = row[ep_class] + row[R + ep_class]
            if total < consume:
                return False
            row[ep_class] = total - consume
            row[R + ep_class] = 0
        else:
            if row.size <= ep_class or row[ep_class] < consume:
                return False
            row[ep_class] -= consume
        overrides[isf] = row

    # POST: produce `weight` tokens (per firing).
    for pev in glevent.passive:
        if pev.event != EventType.POST:
            continue
        fp_ind = pev.node
        if fp_ind >= sn.nnodes or not sn.isstateful[fp_ind]:
            continue
        isf = int(sn.nodeToStateful[fp_ind])
        row = np.array(overrides.get(isf, glspace[isf]), dtype=float).ravel().copy()
        fp_class = int(pev.job_class)
        produce = int(pev.weight)
        if row.size <= fp_class:
            return False
        row[fp_class] += produce
        overrides[isf] = row

    return True


def _emit_single(glspace, isf_transition, trans_state, rate, prob, is_completion) -> _GlobalEventResult:
    out = [list(glspace)]
    out[0][isf_transition] = trans_state
    return _GlobalEventResult(out, [rate], [prob], [is_completion])


def _emit_many(glspace, isf_transition, outcomes: list) -> _GlobalEventResult:
    if not outcomes:
        return _empty_result()
    out = []
    rates = []
    probs = []
    is_comp = []
    base = list(glspace)
    for trans_state, rate, prob, completion, overrides in outcomes:
        gl = list(base)
        if overrides is not None:
            for isf, row in overrides.items():
                gl[isf] = row
        else:
            gl[isf_transition] = trans_state
        out.append(gl)
        rates.append(rate)
        probs.append(prob)
        is_comp.append(completion)
    return _GlobalEventResult(out, rates, probs, is_comp)
