"""QRF blocking tables derived from a NetworkStruct.

Twin of ``matlab/src/api/sn/sn_to_qrf_capacity.m`` and
``matlab/src/api/sn/sn_to_qrf_blocking.m``.

``qrf_bas`` describes a Blocking-After-Service network by a finite-capacity
queue ``f`` and an enumeration of the BLOCKING CONFIGURATIONS reachable behind
it. Everything in that enumeration is implied by the model, so it is derived
here rather than demanded from the caller; ``options.config['qrf_params']``
remains an explicit override.

The tables, and the constraint that reads each one in ``qrf_bas``:

======  =====================================================================
f       the ONE finite-capacity queue. The formulation carries a scalar f
        (ZERO4/ZERO7/ZERO8, THM30, THM3I, THM3L all index it), so a model with
        two binding buffers is refused here.
F(i)    min(buffer size, N) for every queue; N where the buffer is unbounded,
        since no queue can hold more than the population.
BB(m,i) 1 iff queue i is blocked in configuration m.
ZZ(m)   nnz(BB(m,:)), the blocking depth of configuration m.
MM(m,0) head of the FIFO blocking order: the queue that takes the slot when f
        completes. ``qrf_bas`` reads ONLY the first column of MM.
MM1(m,j) index of the configuration reached from m when j becomes blocked.
        Read by THM3L alone, at depth ZM-1.
======  =====================================================================

THREE INVARIANTS, each a correctness condition rather than a convention:

1. Configuration 1 (index 0 here, 1 in the emitted tables) MUST be the empty
   one. ZERO4 iterates ``m = 2:MR`` and ZERO5/ZERO7/ZERO8 test ``m >= 2`` to
   mean "some queue is blocked".
2. ``ZM = max(ZZ)`` MUST be the reachable maximum. ``qrf_bas`` recomputes ZM
   from ZZ and closes the depth ladder there, so a truncated enumeration
   excises states the real chain visits and the polytope stops containing the
   true distribution -- the bound stops bounding. The size guard therefore
   REFUSES; it never truncates.
3. Blocking APPENDS at the tail: a queue that becomes blocked joins behind
   those already waiting, so MM1's successor is the configuration with j
   appended, and the head -- the queue MM(m,0) names -- never moves.

THE ENUMERATION IS THE FULL ORDERED ONE, and it has to be. A (set, head)
collapse looks sound -- the LP reads configurations only through BB, ZZ, MM(:,0)
and MM1, and both objective and readout sum over m -- and it would shrink MR
from ``sum_z P(B,z)`` to ``1 + sum_z C(B,z)*z``. It was tried and it is WRONG.
Merging the depth-ZM configurations that share a set and a head makes several
THM3L rows, one per depth-(ZM-1) predecessor, reference the SAME merged
successor block. That is extra coupling the fine system does not have, so the
collapsed polytope is strictly SMALLER, not a projection of the fine one, and it
can cut off the true distribution. Measured on a 4-station model with three
feeders (B=3, ZM=3, MR 13 collapsed vs 16 full), the collapse reported upper
bounds of 0.681/0.979/0.768 where the full enumeration gives 0.709/0.982/0.800:
tighter, from a coarser state space, which is the signature of a cut that is not
valid.

So MR is factorial in the number of feeders B, and the size guard is what keeps
that honest: it REFUSES an oversized instance rather than trimming the
enumeration, because trimming is the same unsound cut by another name.

WHO CAN BE BLOCKED is read from ``sn.isbasblocking`` / ``sn.isbasdestination``,
not from ``sn.droprule``: LINE accepts the BAS declaration on the upstream
station or on the full destination, and reading droprule at the capped station
sees only the second (BUG-83).
"""

from itertools import combinations, permutations

import numpy as np

__all__ = ['sn_to_qrf_capacity', 'sn_to_qrf_blocking']

# Variable-count ceiling of the derived LP. A guard, not a tuning knob: the
# enumeration cannot be truncated (invariant 2), so an oversized model is
# refused rather than approximated.
QRF_DEFAULT_MAXVARS = 5e5


def sn_to_qrf_capacity(sn):
    """Return ``(F, binding, msg)``: the QRF occupancy bound of every station.

    F is an occupancy bound rather than a declared capacity: the station's
    buffer where that buffer BINDS, and the population N everywhere else, since
    no queue of a closed model can hold more than N jobs.

    Binding is decided by ``sn_get_buffer_size``, the single place in LINE that
    makes that call: refreshCapacity derives a finite classcap (the chain
    population) at every station of every closed model, so a plain finiteness
    test on ``sn.cap`` would report a buffer at every station.

    Both QRF blocking bounds need this. ``qrf.bas`` needs it beside the
    blocking tables; ``qrf.rsrd`` needs it ALONE, since its PBB constraint
    reads only which queues can be full and it carries no blocking tables.
    """
    from ..me.solver_nc_mem import sn_get_buffer_size

    M = int(sn.nstations)
    F = np.zeros(M, dtype=int)
    binding = np.zeros(M, dtype=bool)

    N = float(np.sum(np.asarray(sn.njobs, dtype=float)))
    if not np.isfinite(N) or N < 1:
        return F, binding, 'the QRF bounds need a closed model with a finite population.'
    N = int(round(N))

    for i in range(M):
        b = sn_get_buffer_size(sn, i)
        binding[i] = bool(np.isfinite(b))
        F[i] = N if (not np.isfinite(b) or b > N) else int(round(b))
        if F[i] < 1:
            return F, binding, (
                'station %d has capacity %d: the QRF bounds need every queue to be able to '
                'hold at least one job.' % (i + 1, F[i]))
    return F, binding, ''


def sn_to_qrf_blocking(sn, options=None):
    """Return ``(blk, msg)``: the derived BAS blocking tables, or why not.

    ``blk`` is a dict with keys f (1-based, as the LP indexes it), F, MR, BB,
    MM, ZZ, MM1, ZM and blockers (1-based station indices). ``msg`` is empty on
    success and otherwise says why the tables cannot be derived.
    """
    config = {}
    if options is not None:
        config = getattr(options, 'config', None) or {}

    M = int(sn.nstations)
    N = int(round(float(np.sum(np.asarray(sn.njobs, dtype=float)))))

    F, binding, msg = sn_to_qrf_capacity(sn)
    if msg:
        return None, msg

    fcand = [i for i in range(M) if binding[i]]
    if not fcand:
        # No binding buffer: callers gate on sn_has_blocking first, so this is
        # a defensive branch rather than a normal path.
        return _empty_blocking(F, 1, M), ''
    if len(fcand) > 1:
        names = ', '.join(_station_name(sn, i) for i in fcand)
        return None, (
            "'qrf.bas' models a single finite-capacity queue (its f is a scalar), but %d "
            "stations have a binding buffer: %s. Use 'qrf.rsrd', whose PBB constraint sums "
            "over every full queue and therefore admits several, or cap only one station."
            % (len(fcand), names))
    f = fcand[0]

    blockers = _qrf_blockers(sn, f, M)

    # Blocking needs f at capacity plus one held job per blocked queue, so the
    # population caps the depth as tightly as the feeder count does.
    ZM = min(len(blockers), N - int(F[f]))
    if ZM < 0:
        ZM = 0
    if ZM == 0:
        blk = _empty_blocking(F, f + 1, M)
        blk['blockers'] = [b + 1 for b in blockers]
        return blk, ''

    cfg = _enumerate_permutations(blockers, ZM)

    MR = len(cfg)

    # Size guard: refuse, never truncate (invariant 2).
    Ktot = _total_phases(sn, M)
    n_vars = MR * (N + 1) ** 2 * Ktot ** 2 + Ktot
    max_vars = config.get('qrf_maxvars', QRF_DEFAULT_MAXVARS) or QRF_DEFAULT_MAXVARS
    if n_vars > max_vars:
        return None, (
            'the QRF BAS linear program for this model would carry %.3g variables (MR=%d '
            'blocking configurations, N=%d, %d service phases in total), above the '
            "config['qrf_maxvars'] limit of %.3g. The enumeration cannot be truncated -- a "
            'depth below the reachable maximum ZM=%d excises states the chain visits, and the '
            'result would no longer bound. Reduce the population, the number of stations '
            "feeding %s, or the phase counts; or raise config['qrf_maxvars'] deliberately."
            % (n_vars, MR, N, Ktot, max_vars, ZM, _station_name(sn, f)))

    BB = np.zeros((MR, M), dtype=int)
    ZZ = np.zeros(MR, dtype=int)
    MM = np.zeros((MR, max(2, len(blockers))), dtype=int)
    MM1 = np.zeros((MR, M), dtype=int)

    for m, seq in enumerate(cfg):
        ZZ[m] = len(seq)
        for z, j in enumerate(seq):
            BB[m, j] = 1
            MM[m, z] = j + 1  # only column 0 is read; the rest records the full order

    # MM1: successor under "j becomes blocked". Blocking appends at the tail,
    # so the head is preserved (invariant 3) and the successor is unique.
    index = {_cfg_key(seq): m for m, seq in enumerate(cfg)}
    for m, seq in enumerate(cfg):
        if ZZ[m] >= ZM:
            continue  # THM3L reads MM1 only below ZM
        for j in blockers:
            if BB[m, j]:
                continue
            mp = index.get(_cfg_key(tuple(seq) + (j,)))
            if mp is not None:
                MM1[m, j] = mp + 1

    return {
        'f': f + 1,
        'F': F,
        'MR': MR,
        'BB': BB,
        'MM': MM,
        'ZZ': ZZ,
        'MM1': MM1,
        'ZM': ZM,
        'blockers': [b + 1 for b in blockers],
    }, ''


def _empty_blocking(F, f_one_based, M):
    """The one-configuration table for a model in which no blocking state is
    reachable. ``qrf_bas`` reads it as a plain finite-buffer network."""
    return {
        'f': f_one_based,
        'F': F,
        'MR': 1,
        'BB': np.zeros((1, M), dtype=int),
        'MM': np.zeros((1, 2), dtype=int),
        'ZZ': np.zeros(1, dtype=int),
        'MM1': np.zeros((1, M), dtype=int),
        'ZM': 0,
        'blockers': [],
    }


def _qrf_blockers(sn, f, M):
    """0-based station indices that hold a completed job when f is full.

    A blocker must route into f, must not be f, and must not be an infinite
    server (which has a server per job and cannot be held). BAS itself is read
    from the BUG-83 fields, with a structural fallback for an sn built without
    the local-variable refresh.
    """
    from .network_struct import DropStrategy, SchedStrategy

    declared = np.zeros(M, dtype=bool)
    marker = getattr(sn, 'isbasblocking', None)
    station_to_node = getattr(sn, 'stationToNode', None)
    if marker is not None and station_to_node is not None:
        marker = np.asarray(marker).ravel()
        station_to_node = np.asarray(station_to_node).ravel()
        for i in range(M):
            if i >= station_to_node.size:
                continue
            ind = int(station_to_node[i])
            if 0 <= ind < marker.size and marker[ind] == 1:
                declared[i] = True

    if not declared.any():
        # Fallback: BAS declared on the upstream station or on the full
        # destination, read straight off sn.droprule -- the same two forms
        # declaresBlockedMarker resolves.
        droprule = getattr(sn, 'droprule', None)
        if droprule is not None:
            droprule = np.asarray(droprule)
            if droprule.ndim == 2 and f < droprule.shape[0]:
                dest_bas = bool(np.any(droprule[f, :] == int(DropStrategy.BAS)))
                for i in range(M):
                    if i == f or i >= droprule.shape[0]:
                        continue
                    if dest_bas or np.any(droprule[i, :] == int(DropStrategy.BAS)):
                        declared[i] = True

    R = int(sn.nclasses)
    rt = np.asarray(sn.rt, dtype=float)
    sched = getattr(sn, 'sched', None)

    blockers = []
    for i in range(M):
        if i == f or not declared[i]:
            continue
        if sched is not None and int(sched[i]) == int(SchedStrategy.INF):
            continue
        # Summed over class pairs so the test survives a multiclass sn, even
        # though the QRF gate upstream admits one class only.
        block = rt[i * R:(i + 1) * R, f * R:(f + 1) * R]
        if np.any(block > 0):
            blockers.append(i)
    return sorted(blockers)


def _enumerate_permutations(blockers, ZM):
    """Every ordered sequence of distinct blockers up to length ZM, first entry
    the head.

    Deterministic order -- depth ascending, then subsets lexicographic by
    ascending station index, then the orders of each subset sorted -- so every
    codebase emits identical tables. The empty configuration is first
    (invariant 1).
    """
    cfg = [()]
    for z in range(1, ZM + 1):
        for members in combinations(blockers, z):
            for p in sorted(permutations(members)):
                cfg.append(tuple(p))
    return cfg


def _cfg_key(seq):
    """Identity of a configuration: the whole blocking order, since that is what
    distinguishes configurations in the enumeration ``qrf_bas`` is entitled to."""
    return tuple(seq)


def _total_phases(sn, M):
    """Total service phases across stations, which sizes the QRF variable space
    together with MR and the population."""
    total = 0
    proc = getattr(sn, 'proc', None)
    for i in range(M):
        ki = 1
        try:
            entry = proc[i][0]
            if entry is not None and len(entry) and entry[0] is not None:
                ki = int(np.asarray(entry[0]).shape[0])
        except (TypeError, KeyError, IndexError):
            ki = 1
        total += max(1, ki)
    return total


def _station_name(sn, ist):
    """Printable station name, falling back to the index."""
    names = getattr(sn, 'nodenames', None)
    station_to_node = getattr(sn, 'stationToNode', None)
    if names is not None and station_to_node is not None:
        station_to_node = np.asarray(station_to_node).ravel()
        if ist < station_to_node.size:
            ind = int(station_to_node[ist])
            if 0 <= ind < len(names):
                return str(names[ind])
    return 'station %d' % (ist + 1)
