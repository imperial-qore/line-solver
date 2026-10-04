"""
Conservation laws of a layered queueing network, enumerated from its structure.

A layered model is not free to report any tuple of throughputs, think times and
utilizations: five families of relations tie them together, and every one of them
is fixed by the STRUCTURE of the model alone. This module walks a
LayeredNetworkStruct and emits them, one record per relation, with the index sets
to aggregate over, the constant coefficients and a printable form. Nothing is
solved here.

The families, with ``kind`` as emitted:

``little``
    Little's law on a task's THREAD POOL. The threads of task t form a closed
    cycle of one delay stage (the surrogate think time SolverLN imputes to the
    task, plus the declared think time of a reference task) and one service stage
    (holding a request from above). With B(t,k) the mean number of threads of t
    busy serving caller class k -- the per-class utilization in JOB units --

        X(t)*(Z(t) + z(t)) + sum_k B(t,k) = N(t)

    which is the update ``SolverLN.update_think_times`` iterates on. In the
    [0,1]-normalized utilization LINE reports for a queueing station,
    B(t,k) = N(t)*U(t,k), giving X(t)*(Z(t)+z(t)) = N(t)*(1 - sum_k U(t,k)); at an
    infinite server the utilization is already a job count, so B = U. The caller
    classes k are the CALLS targeting an entry of t -- the in-edges of t in the
    call graph -- plus each entry of t carrying an OPEN ARRIVAL, a stream that
    holds a thread exactly as a call does and that a task can have alongside its
    callers. Both are structural neighbours of t, which makes the relation
    node-local; ``termisentry`` says which of the two a term is.

``callflow``
    Throughput conservation across one call, X(c) = X(src(c))*y(c), with y(c) the
    mean number of calls and src(c) the dispatching activity (the dispatching
    ENTRY for a forwarding call).

``entryflow``
    The requests an entry serves are the calls reaching it plus its open-arrival
    stream, X(e) = sum_c X(c) + lambda(e).

``actflow``
    An activity executes v(a) times per invocation of its entry,
    X(a) = X(e)*v(a), with v from the activity precedence graph. An AND-JOIN is
    the one place where flow does not add up -- its target executes once per fork,
    not once per branch -- so the arcs into a join are scaled by 1/(number of
    joined branches).

``hostutil``
    The utilization law at a processor, sum_a X(a)*D(a) = m(h)*U(h); m(h)*U(h) is
    a job count and the factor m(h) drops at an infinite server.

Together these close the system: ``little`` alone is one equation per task and
admits the all-zero solution, so a physics-informed loss built on it should carry
the flow and utilization families as well.

Twin of the MATLAB ``lqn_balance_equations.m``, the JAR
``jline.api.lqn.LqnBalanceEquations`` and the C++
``line/api/lqn/lqn_balance_equations.h``.
"""

from typing import Any, Dict, List, Optional, Sequence

import numpy as np

from ...constants import GlobalConstants, SchedStrategy
from .lqn_ph import call_hashname

__all__ = ['BalanceRelation', 'BalanceEquations', 'balance_equations']

# CallType ids on the wire, matching MATLAB CallType and the JAR enum ordinals.
_SYNC = 1
_ASYNC = 2
_FWD = 3

# ActivityPrecedenceType.PRE_AND, carried on the JOINED predecessors of a join.
_PRE_AND = 2


class BalanceRelation:
    """
    One conservation law, as an aggregation over a node's structural neighbours.

    Attributes:
        kind: 'little' | 'callflow' | 'entryflow' | 'actflow' | 'hostutil'
        branch: which case produced a ``little`` record ('ref', 'inf',
            'queueing', 'fwd', 'arrival'), the call type of a ``callflow``, or
            the server kind of a ``hostutil``
        target: absolute index the relation is anchored on
        terms: absolute indices to aggregate over (calls and arrival-carrying
            entries for ``little``, activities for ``hostutil``)
        termisentry: per term, True when it is an ENTRY class (an open-arrival
            stream, or the self-driven cycle of a reference task) rather than a
            CALL class; the two read their rate from different vectors
        coeff: constant coefficient of each term
        const: right-hand-side constant
        scaled: the per-class utilization must be multiplied by ``mult`` to reach
            job units (a queueing server); False at an infinite server
        degenerate: the relation cannot serve as a residual (infinite
            multiplicity, empty term set, zero call count)
        clamped: set once instantiated; the equality is unattainable because the
            task is saturated and the think time is pinned at zero
    """

    __slots__ = ('kind', 'branch', 'target', 'targetname', 'terms', 'termisentry',
                 'coeff', 'const', 'mult', 'maxmult', 'repl', 'scaled', 'phase2',
                 'setup', 'degenerate', 'text', 'lhs', 'rhs', 'residual',
                 'relresidual', 'clamped', 'perclassutil')

    def __init__(self, kind: str = '', target: int = 0, targetname: str = ''):
        self.kind = kind
        self.branch = ''
        self.target = target
        self.targetname = targetname
        self.terms: List[int] = []
        self.termisentry: List[bool] = []
        self.coeff: List[float] = []
        self.const = 0.0
        self.mult = float('nan')
        self.maxmult = float('nan')
        self.repl = 1.0
        self.scaled = False
        self.phase2 = False
        self.setup = False
        self.degenerate = False
        self.text = ''
        self.lhs = float('nan')
        self.rhs = float('nan')
        self.residual = float('nan')
        self.relresidual = float('nan')
        self.clamped = False
        self.perclassutil: List[float] = []

    def __repr__(self) -> str:
        return 'BalanceRelation(%s, %s, %s)' % (self.kind, self.targetname, self.text)


class BalanceEquations:
    """
    The relation set of one layered model.

    Attributes:
        eqs: the relations, grouped by family in emission order
        visits: (nidx,) expected executions of each activity per invocation of
            its entry
        A_little: (ntasks, ncalls) 1 where call c is a caller class of task t
        A_flow: (ncalls, nidx) call-to-source incidence weighted by y(c)
        A_host: (nhosts, nidx) host-to-activity incidence weighted by D(a)
        maxresidual: largest absolute residual over the non-degenerate relations,
            NaN when no solution was supplied
    """

    def __init__(self):
        self.eqs: List[BalanceRelation] = []
        self.visits: Optional[np.ndarray] = None
        self.convention: List[str] = []
        self.A_little: Optional[np.ndarray] = None
        self.A_flow: Optional[np.ndarray] = None
        self.A_host: Optional[np.ndarray] = None
        self.maxresidual = float('nan')
        self.text: List[str] = []

    def __str__(self) -> str:
        return '\n'.join(self.text)

    def print(self) -> None:
        print(str(self))


def balance_equations(lqn, sol: Any = None, un: Any = None) -> BalanceEquations:
    """
    Enumerate the conservation laws of the layered model LQN.

    Args:
        lqn: LayeredNetworkStruct, from ``LayeredNetwork.get_struct()``
        sol: optional; a solved SolverLN, or any object exposing the same
            iterates (``tput``, ``util``, ``thinkt``, ``servt``, ``residt``),
            used to instantiate every relation and report its residual
        un: optional; the REPORTED utilization per element, i.e. the second
            return of ``get_ensemble_avg()``. Needed by the ``hostutil`` family,
            whose right-hand side is the reported processor utilization and not
            the ``util`` iterate (which a host leaves at zero). It is a separate
            argument because ``get_ensemble_avg`` re-enters ``iterate()``, and a
            diagnostic must not re-run a fixed point on its own.

    Returns:
        BalanceEquations

    Conventions:
        Rates and populations in a ``little`` record are PER REPLICA, matching
        ``update_think_times``: X is tput/repl and N is the multiplicity of one
        copy. Elsewhere throughputs and utilizations are as the solver reports
        them, totalled over replicas. N(t) is ``lqn.mult``; SolverLN iterates on
        ``njobs``, which carries the interlocking corrections and may be
        ``maxmult`` under replication, so both are returned per record. S(k) is
        the entry SERVICE time (phase 1 plus phase 2), the time a thread is held,
        not the residence time the caller waits for; the difference is the
        phase-2 tail, flagged by ``phase2``.
    """
    nidx = int(lqn.nidx)
    visits = _act_visits(lqn)
    calls_into = _calls_into(lqn)
    s = _read_solution(lqn, sol, un) if sol is not None else None

    eqs: List[BalanceRelation] = []

    for t in range(int(lqn.ntasks)):
        tidx = int(lqn.tshift) + t
        branch, terms = _thread_pool_branch(lqn, tidx, calls_into.get(tidx, []))
        if branch is None:
            continue
        eqs.append(_little_record(lqn, tidx, branch, terms, s))

    for cidx in range(int(lqn.ncalls)):
        eqs.append(_callflow_record(lqn, cidx, s))
    for e in range(int(lqn.nentries)):
        eidx = int(lqn.eshift) + e
        inc = _incoming_calls(lqn, eidx)
        lam = _arrival_rate(lqn, eidx)
        if not inc and lam == 0.0:
            continue    # a reference entry is driven by its own cycle, not by flow
        eqs.append(_entryflow_record(lqn, eidx, inc, lam, s))
    for a in range(int(lqn.nacts)):
        aidx = int(lqn.ashift) + a
        eidx = _entry_of_activity(lqn, aidx)
        if eidx < 0:
            continue
        eqs.append(_actflow_record(lqn, aidx, eidx, float(visits[aidx]), s))

    for h in range(int(lqn.nhosts)):
        eqs.append(_hostutil_record(lqn, h, s))

    out = BalanceEquations()
    out.eqs = eqs
    out.visits = visits
    out.convention = _convention_text()
    out.A_little, out.A_flow, out.A_host = _incidences(lqn, eqs, nidx)
    if s is not None:
        res = [abs(r.residual) for r in eqs
               if not r.degenerate and not np.isnan(r.residual)]
        if res:
            out.maxresidual = max(res)
    out.text = _report(lqn, out, s is not None)
    return out


# ----------------------------------------------------------------------
# record construction
# ----------------------------------------------------------------------

def _thread_pool_branch(lqn, tidx: int, cin: Sequence[int]):
    """
    Which thread-pool case task TIDX falls in, and the caller classes of its pool.

    The term set is uniform across the cases: every request stream that can hold
    a thread of the task contributes one class. That is each CALL targeting one of
    its entries, plus each entry carrying an OPEN ARRIVAL -- a stream that holds
    threads exactly as a call does, and that a task can have alongside its callers
    -- plus, for a reference task, its own entries, since the cycle of a reference
    task closes on itself and no layer above drives it.

    The ``branch`` label follows the case analysis of update_think_times, which is
    what decides whether the utilization is a job count (infinite server) or is
    normalized to [0,1] (every other discipline).
    """
    ents = [int(e) for e in lqn.entriesof.get(tidx, [])]
    arv = [e for e in ents if _arrival_rate(lqn, e) > 0]
    terms = [('call', c) for c in cin] + [('entry', e) for e in arv]
    if _isref(lqn, tidx):
        return 'ref', terms + [('entry', e) for e in ents if e not in arv]
    blocking = [c for c in cin if _calltype(lqn, c) != _FWD]
    if blocking:
        branch = 'inf' if _sched_of(lqn, tidx) == SchedStrategy.INF else 'queueing'
        return branch, terms     # a forwarded request holds a thread too
    if cin:
        return 'fwd', terms
    if arv:
        return 'arrival', terms
    return None, []   # no caller, no arrival: no cycle to close


def _little_record(lqn, tidx, branch, terms, s):
    r = BalanceRelation('little', tidx, _elem_name(lqn, tidx))
    r.branch = branch
    r.termisentry = [k == 'entry' for k, _ in terms]
    r.terms = [i for _, i in terms]
    r.coeff = [1.0] * len(terms)
    r.mult = _mult(lqn, tidx)
    r.maxmult = _maxmult(lqn, tidx)
    r.repl = max(1.0, _repl(lqn, tidx))
    r.const = r.mult
    r.scaled = branch != 'inf'
    r.phase2 = _has_phase2(lqn, tidx)
    r.setup = _has_setup(lqn, tidx)
    r.degenerate = (not np.isfinite(r.const)) or not terms

    z = _ref_thinktime(lqn, tidx)
    names = [_term_name(lqn, k, i) for k, i in terms]
    r.text = _little_text(r, z, names)

    if s is not None:
        X, B, U = _little_values(lqn, tidx, terms, s)
        zt = s['thinkt'][tidx]
        if np.isnan(zt):
            zt = 0.0
        r.lhs = X * (zt + z) + float(sum(B))
        r.rhs = r.const
        r.residual = r.lhs - r.rhs
        r.relresidual = r.residual / max(abs(r.rhs), GlobalConstants.FineTol)
        r.clamped = np.isfinite(r.const) and (r.const - sum(B) - X * z) < 0
        r.perclassutil = list(U)
    return r


def _little_text(r, z, names) -> str:
    X = 'X(%s)' % r.targetname
    Z = 'Z(%s)' % r.targetname
    zs = ' + %g' % z if z > 0 else ''
    lhs = '%s*(%s%s)' % (X, Z, zs)
    for nm in names:
        lhs += ' + B(%s)' % nm
    txt = '%s = %s' % (lhs, _num(r.const))
    if r.scaled and names:
        us = ''.join(' - U(%s,%s)' % (r.targetname, nm) for nm in names)
        txt += '\n           equivalently  %s*%s = %s*(1%s)' % (X, Z, _num(r.const), us)
    return txt


def _little_values(lqn, tidx, terms, s):
    """
    Per-caller-class busy threads and normalized utilization of task TIDX.

    A class holds a thread for the SERVICE time of the entry it targets, at the
    rate it drives that entry: the call rate for a call class, the entry rate for
    an entry class (an open-arrival stream, or the self-driven cycle of a
    reference task).
    """
    X = s['tput'][tidx] / max(1.0, _repl(lqn, tidx))
    if np.isnan(X):
        X = 0.0
    B = []
    for kind, i in terms:
        if kind == 'entry':
            b = s['tput'][i] * s['servt'][i]
        else:
            b = s['calltput'][i] * s['servt'][int(lqn.callpair[i, 1])]
        B.append(0.0 if np.isnan(b) else float(b))
    m = _mult(lqn, tidx)
    if _sched_of(lqn, tidx) == SchedStrategy.INF or not np.isfinite(m) or m <= 0:
        U = list(B)
    else:
        U = [b / m for b in B]
    return float(X), B, U


def _callflow_record(lqn, cidx, s):
    r = BalanceRelation('callflow', _cshift(lqn) + cidx, _call_name(lqn, cidx))
    r.branch = _calltype_name(_calltype(lqn, cidx))
    src = int(lqn.callpair[cidx, 0])
    y = float(lqn.callpair[cidx, 2])
    r.terms = [src]
    r.coeff = [y]
    r.degenerate = (y == 0.0)
    r.text = 'X(%s) = X(%s) * %g' % (r.targetname, _elem_name(lqn, src), y)
    if s is not None:
        r.lhs = _z(s['calltput'][cidx])
        r.rhs = _z(s['tput'][src]) * y
        r.residual = r.lhs - r.rhs
        r.relresidual = r.residual / max(abs(r.rhs), GlobalConstants.FineTol)
    return r


def _entryflow_record(lqn, eidx, inc, lam, s):
    r = BalanceRelation('entryflow', eidx, _elem_name(lqn, eidx))
    r.terms = list(inc)
    r.coeff = [1.0] * len(inc)
    r.const = lam
    txt = 'X(%s) =' % r.targetname
    for i, c in enumerate(inc):
        txt += (' X(%s)' if i == 0 else ' + X(%s)') % _call_name(lqn, c)
    if lam > 0:
        txt += (' %g' if not inc else ' + %g') % lam
        txt += '   (open arrival)'
    r.text = txt
    if s is not None:
        r.lhs = _z(s['tput'][eidx])
        r.rhs = lam + sum(_z(s['calltput'][c]) for c in inc)
        r.residual = r.lhs - r.rhs
        r.relresidual = r.residual / max(abs(r.rhs), GlobalConstants.FineTol)
    return r


def _actflow_record(lqn, aidx, eidx, v, s):
    r = BalanceRelation('actflow', aidx, _elem_name(lqn, aidx))
    r.terms = [eidx]
    r.coeff = [v]
    r.degenerate = not np.isfinite(v)
    r.text = 'X(%s) = X(%s) * %g' % (r.targetname, _elem_name(lqn, eidx), v)
    if s is not None:
        r.lhs = _z(s['tput'][aidx])
        r.rhs = _z(s['tput'][eidx]) * v
        r.residual = r.lhs - r.rhs
        r.relresidual = r.residual / max(abs(r.rhs), GlobalConstants.FineTol)
    return r


def _hostutil_record(lqn, h, s):
    hidx = int(getattr(lqn, 'hshift', 0)) + h
    r = BalanceRelation('hostutil', hidx, _elem_name(lqn, hidx))
    r.mult = _mult(lqn, hidx)
    r.repl = max(1.0, _repl(lqn, hidx))
    r.scaled = _sched_of(lqn, hidx) != SchedStrategy.INF
    acts: List[int] = []
    dem: List[float] = []
    for tidx in lqn.tasksof.get(hidx, []):
        for aidx in lqn.actsof.get(int(tidx), []):
            d = _hostdem(lqn, int(aidx))
            if d == 0.0 or np.isnan(d):
                continue
            acts.append(int(aidx))
            dem.append(d)
    r.terms = acts
    r.coeff = dem
    # The server count is the declared multiplicity ALONE, not mult*repl: a
    # replicated host reports its throughputs and its utilization as TOTALS over
    # the copies, so the extra factor would double-count the replication. Checked
    # on the two-replica processor of lqn_sockshop, where mult*repl overshoots the
    # reported utilization by exactly the replication factor.
    m = r.mult
    if not r.scaled or not np.isfinite(m):
        m = 1.0
    r.const = m
    r.branch = 'queueing' if r.scaled else 'inf'
    r.degenerate = not acts
    lhs = ' + '.join('X(%s)*%g' % (_elem_name(lqn, a), d) for a, d in zip(acts, dem))
    r.text = '%s = %s*U(%s)' % (lhs if lhs else '0', _num(m), r.targetname)
    if s is not None and not np.isnan(s['un'][hidx]):
        r.lhs = sum(_z(s['tput'][a]) * d for a, d in zip(acts, dem))
        r.rhs = m * float(s['un'][hidx])
        r.residual = r.lhs - r.rhs
        r.relresidual = r.residual / max(abs(r.rhs), GlobalConstants.FineTol)
    return r


# ----------------------------------------------------------------------
# structural helpers
# ----------------------------------------------------------------------

def _act_visits(lqn) -> np.ndarray:
    """
    Expected executions of every activity per invocation of its entry.

    The activity precedence arcs of a task are a transient Markov chain whose
    absorbing state is the reply, so the visit counts of the block solve
    v = e0*(I-P)^-1: a loop back-edge of weight 1-1/count returns count, and an
    AND-fork row summing above one returns the branching expectation. Call arcs
    leave the block and drop out.

    An AND-JOIN is the one place where flow does not add up. Its target executes
    ONCE per fork, not once per branch, so summing the inbound arcs would count
    it as many times as there are branches; the arcs into a join are therefore
    scaled by 1/(number of joined branches), which recovers the rate of one
    branch exactly when the branches carry equal rate -- the case for a
    well-formed fork/join block.
    """
    nidx = int(lqn.nidx)
    v = np.zeros(nidx)
    G = _join_scaled_graph(lqn)
    for t in range(int(lqn.ntasks)):
        tidx = int(lqn.tshift) + t
        A = sorted(set(int(a) for a in lqn.actsof.get(tidx, [])))
        if not A:
            continue
        pos = dict((a, i) for i, a in enumerate(A))
        P = G[np.ix_(A, A)]
        M = np.eye(len(A)) - P
        for eidx in lqn.entriesof.get(tidx, []):
            bound = [a for a in A if G[int(eidx), a] != 0]
            if not bound:
                continue
            e0 = np.zeros(len(A))
            for a in bound:
                e0[pos[a]] = 1.0
            v[A] += np.linalg.solve(M.T, e0)
    return v


def _join_scaled_graph(lqn) -> np.ndarray:
    """
    The precedence graph with the arcs into every AND-join target divided by the
    number of branches the join waits for. ``actpretype`` marks the joined
    PREDECESSORS (value PRE_AND = 2), so a join target is any successor of one.
    """
    G = _graph(lqn).copy()
    pre = getattr(lqn, 'actpretype', None)
    if pre is None:
        return G
    pre = np.asarray(pre).ravel()
    andpre = [i for i in range(len(pre)) if int(pre[i]) == _PRE_AND]
    if not andpre:
        return G
    n = G.shape[0]
    for j in range(n):
        joined = [i for i in andpre if i < n and G[i, j] != 0]
        if len(joined) > 1:
            for i in joined:
                G[i, j] = G[i, j] / len(joined)
    return G


def _calls_into(lqn) -> Dict[int, List[int]]:
    """Calls targeting each task, i.e. the caller classes of its thread pool."""
    out: Dict[int, List[int]] = {}
    for cidx in range(int(lqn.ncalls)):
        dst = int(lqn.callpair[cidx, 1])
        if dst < 0 or dst >= int(lqn.nidx):
            continue
        tidx = _parent(lqn, dst)
        out.setdefault(tidx, []).append(cidx)
    return out


def _incoming_calls(lqn, eidx: int) -> List[int]:
    return [c for c in range(int(lqn.ncalls))
            if int(lqn.callpair[c, 1]) == eidx]


def _arrival_rate(lqn, eidx: int) -> float:
    arv = getattr(lqn, 'arrival', None)
    if not arv:
        return 0.0
    d = arv.get(eidx, None) if isinstance(arv, dict) else None
    if d is None:
        return 0.0
    m = float(d.getMean()) if hasattr(d, 'getMean') else float(d)
    if np.isfinite(m) and m > GlobalConstants.FineTol:
        return 1.0 / m
    return 0.0


def _entry_of_activity(lqn, aidx: int) -> int:
    """
    Entry an activity belongs to: the one whose activity block contains it. An
    activity shared by two entries is attributed to the first, matching how
    _act_visits accumulates its visit counts.
    """
    tidx = _parent(lqn, aidx)
    for e in lqn.entriesof.get(tidx, []):
        if aidx in [int(a) for a in lqn.actsof.get(int(e), [])]:
            return int(e)
    return -1      # 0 is the first host in the 0-based band, so it cannot mean "none"


def _has_phase2(lqn, tidx: int) -> bool:
    ph = getattr(lqn, 'actphase', None)
    if ph is None:
        return False
    ph = np.asarray(ph).ravel()
    for eidx in lqn.entriesof.get(tidx, []):
        for aidx in lqn.actsof.get(int(eidx), []):
            # actphase is 0-based over activities here, 1-based in MATLAB
            a = int(aidx) - int(lqn.ashift) - 1
            if 0 <= a < len(ph) and ph[a] > 1:
                return True
    return False


def _has_setup(lqn, tidx: int) -> bool:
    hs = getattr(lqn, 'hassetup', None)
    if hs is None:
        return False
    hs = np.asarray(hs).ravel()
    return bool(tidx < len(hs) and hs[tidx])


def _incidences(lqn, eqs, nidx):
    ntasks, ncalls, nhosts = int(lqn.ntasks), int(lqn.ncalls), int(lqn.nhosts)
    Al = np.zeros((ntasks, ncalls))
    Af = np.zeros((ncalls, nidx))
    Ah = np.zeros((nhosts, nidx))
    for r in eqs:
        if r.kind == 'little':
            # the call classes only; an entry class (open arrival, or the
            # self-driven cycle of a reference task) is not a call
            cols = [i for e, i in zip(r.termisentry, r.terms) if not e]
            if cols:
                Al[r.target - int(lqn.tshift), cols] = 1.0
        elif r.kind == 'callflow':
            Af[r.target - _cshift(lqn), r.terms] = r.coeff
        elif r.kind == 'hostutil' and r.terms:
            Ah[r.target - int(getattr(lqn, 'hshift', 0)), r.terms] = r.coeff
    return Al, Af, Ah


# ----------------------------------------------------------------------
# solution adapter
# ----------------------------------------------------------------------

def _read_solution(lqn, sol, un) -> Dict[str, np.ndarray]:
    """
    Read the iterate vectors a solved SolverLN carries. Any object exposing
    tput/util/thinkt/servt/residt works, which is what makes this usable on an
    external prediction as well as on LINE's own fixed point.
    """
    n = int(lqn.nidx)
    s = {}
    for name in ('tput', 'util', 'thinkt', 'servt', 'residt'):
        s[name] = _solvec(sol, name, n)
    s['un'] = _reported_util(lqn, sol, un)
    # A call throughput is not an iterate of SolverLN: a call inherits the rate
    # of its dispatching element scaled by the mean call count.
    ct = np.full(int(lqn.ncalls), np.nan)
    for cidx in range(int(lqn.ncalls)):
        src = int(lqn.callpair[cidx, 0])
        y = float(lqn.callpair[cidx, 2])
        if 0 <= src < n:
            ct[cidx] = s['tput'][src] * y
    s['calltput'] = ct
    return s


def _reported_util(lqn, sol, un) -> np.ndarray:
    """
    The REPORTED utilization of every element, which is not the same quantity as
    the ``util`` iterate: the iterate holds a task's utilization as a SERVER in
    its own task layer -- the U that closes the thread-pool cycle -- and it is
    left at zero on a host. The utilization law at a processor is about the
    reported UN, which the caller supplies.

    It is NOT read off the solver here. ``get_ensemble_avg`` re-enters
    ``iterate()`` in every codebase, and a diagnostic must not re-run a fixed
    point as a side effect of being asked a question. Pass ``un`` explicitly:

        QN, UN, RN, TN, AN, WN = solver.get_ensemble_avg()
        out = balance_equations(lqn, solver, UN)
    """
    n = int(lqn.nidx)
    if un is not None:
        return _solvec({'UN': un}, 'UN', n)
    if isinstance(sol, dict):
        for key in ('un', 'UN'):
            if key in sol:
                return _solvec(sol, key, n)
    return np.full(n, np.nan)


def _solvec(sol, name: str, n: int) -> np.ndarray:
    v = np.full(n, np.nan)
    w = sol.get(name, None) if isinstance(sol, dict) else getattr(sol, name, None)
    if w is None:
        return v
    w = np.asarray(w, dtype=float).ravel()
    m = min(len(w), n)
    v[:m] = w[:m]
    return v


# ----------------------------------------------------------------------
# reporting
# ----------------------------------------------------------------------

def _convention_text() -> List[str]:
    return [
        'Conventions:',
        '  little   : rates and populations PER REPLICA (X = tput/repl, N = mult of one copy).',
        '             B(t,k) is the per-class utilization in JOB units; for a queueing task',
        '             B = mult*U with U in [0,1], at an infinite server B = U directly.',
        '             S(k) is the entry SERVICE time (phase 1 + phase 2), the thread hold time.',
        '  other    : throughputs and utilizations as the solver reports them, totalled over replicas.',
        '  N(t)     : lqn.mult. SolverLN iterates on njobs (interlocking corrections, maxmult under',
        '             replication), reported per record as mult/maxmult.',
    ]


_TITLES = [
    ('little', "thread-pool Little's law"),
    ('callflow', 'call-flow balance'),
    ('entryflow', 'entry-flow balance'),
    ('actflow', 'activity-flow balance'),
    ('hostutil', 'host utilization law'),
]


def _report(lqn, out, has_sol: bool) -> List[str]:
    txt = ['LQN balance equations: %d hosts, %d tasks, %d entries, %d activities, %d calls'
           % (lqn.nhosts, lqn.ntasks, lqn.nentries, lqn.nacts, lqn.ncalls)]
    txt += out.convention
    for kind, title in _TITLES:
        sel = [(i, r) for i, r in enumerate(out.eqs) if r.kind == kind]
        if not sel:
            continue
        txt.append('')
        txt.append("--- %s  (kind='%s', %d relations) ---" % (title, kind, len(sel)))
        for i, r in sel:
            head = '[%3d] %-14s' % (i + 1, r.targetname)
            if r.kind in ('little', 'hostutil'):
                head += ' %-9s mult=%s repl=%g' % (r.branch, _num(r.mult), r.repl)
            elif r.branch:
                head += ' %-9s' % r.branch
            txt.append(head)
            for line in r.text.split('\n'):
                txt.append('       %s' % line)
            notes = []
            if r.phase2:
                notes.append('phase-2 tail on Z')
            if r.setup:
                notes.append('setup charge on Z')
            if r.degenerate:
                notes.append('DEGENERATE (not usable as a residual)')
            if r.clamped:
                notes.append('SATURATED (equality unattainable, use a hinge)')
            if has_sol and not r.degenerate and np.isnan(r.residual):
                notes.append('not instantiated (pass un for the host utilization law)')
            if notes:
                txt.append('       note: %s' % ', '.join(notes))
            if has_sol and not r.degenerate and not np.isnan(r.residual):
                txt.append('       lhs=%.6g  rhs=%.6g  residual=%.3e  rel=%.3e'
                           % (r.lhs, r.rhs, r.residual, r.relresidual))
    if has_sol:
        txt.append('')
        txt.append('max |residual| over %d non-degenerate relations: %.3e'
                   % (sum(1 for r in out.eqs if not r.degenerate), out.maxresidual))
    return txt


# ----------------------------------------------------------------------
# small utilities
# ----------------------------------------------------------------------

def _num(x) -> str:
    if np.isinf(x):
        return 'Inf'
    if float(x) == int(x):
        return '%d' % int(x)
    return '%g' % x


def _z(x) -> float:
    return 0.0 if (x is None or np.isnan(x)) else float(x)


def _graph(lqn) -> np.ndarray:
    G = np.asarray(lqn.graph, dtype=float)
    n = int(lqn.nidx)
    if G.shape[0] < n or G.shape[1] < n:
        H = np.zeros((n, n))
        H[:G.shape[0], :G.shape[1]] = G
        return H
    return G


def _cshift(lqn) -> int:
    """Absolute offset of the call index space; MATLAB places it after nidx."""
    return int(getattr(lqn, 'cshift', 0) or int(lqn.nidx))


def _vecfield(lqn, name: str, idx: int, default=float('nan')) -> float:
    f = getattr(lqn, name, None)
    if f is None:
        return default
    if isinstance(f, dict):
        v = f.get(idx, default)
        return float(v) if v is not None else default
    f = np.asarray(f).ravel()
    return float(f[idx]) if idx < len(f) else default


def _mult(lqn, idx: int) -> float:
    return _vecfield(lqn, 'mult', idx)


def _maxmult(lqn, idx: int) -> float:
    return _vecfield(lqn, 'maxmult', idx)


def _repl(lqn, idx: int) -> float:
    r = _vecfield(lqn, 'repl', idx, 1.0)
    return 1.0 if np.isnan(r) else r


def _hostdem(lqn, aidx: int) -> float:
    """Mean host demand of an activity. Stored as a mean, or as a Distribution."""
    hd = getattr(lqn, 'hostdem', None)
    if hd is None:
        return 0.0
    d = hd.get(aidx, 0.0) if isinstance(hd, dict) else _vecfield(lqn, 'hostdem', aidx, 0.0)
    if d is None:
        return 0.0
    if hasattr(d, 'getMean'):
        d = d.getMean()
    d = float(d)
    return 0.0 if np.isnan(d) else d


def _parent(lqn, idx: int) -> int:
    p = np.asarray(lqn.parent).ravel()
    return int(p[idx]) if idx < len(p) else 0


def _isref(lqn, tidx: int) -> bool:
    f = getattr(lqn, 'isref', None)
    if f is None:
        return False
    f = np.asarray(f).ravel()
    return bool(tidx < len(f) and f[tidx])


def _calltype(lqn, cidx: int) -> int:
    ct = getattr(lqn, 'calltype', None)
    if ct is None:
        return _SYNC
    ct = np.asarray(ct).ravel()
    return int(ct[cidx]) if cidx < len(ct) else _SYNC


def _calltype_name(ct: int) -> str:
    return {_SYNC: 'sync', _ASYNC: 'async', _FWD: 'fwd'}.get(ct, '')


def _sched_of(lqn, idx: int):
    sched = getattr(lqn, 'sched', None)
    if sched is None:
        return SchedStrategy.PS
    v = sched.get(idx, None) if isinstance(sched, dict) else None
    if v is None and not isinstance(sched, dict):
        arr = np.asarray(sched).ravel()
        v = int(arr[idx]) if idx < len(arr) else None
    if v is None:
        return SchedStrategy.PS
    if isinstance(v, SchedStrategy):
        return v
    for st in SchedStrategy:
        if st.value == v:
            return st
    return SchedStrategy.PS


def _ref_thinktime(lqn, tidx: int) -> float:
    """
    Declared think time of task TIDX as it enters the thread cycle: the value for
    a REFERENCE task, zero for any other. Twin of MATLAB lqn_ref_thinktime.
    """
    if not _isref(lqn, tidx):
        return 0.0
    think = getattr(lqn, 'think', None)
    if not think:
        return 0.0
    d = think.get(tidx, None) if isinstance(think, dict) else None
    if d is None:
        return 0.0
    z = float(d.getMean()) if hasattr(d, 'getMean') else float(d)
    return z if np.isfinite(z) and z >= 0 else 0.0


def _elem_name(lqn, idx: int) -> str:
    hn = getattr(lqn, 'hashnames', None)
    if hn is not None:
        hn = np.asarray(hn).ravel()
        if 0 <= idx < len(hn) and str(hn[idx]):
            return str(hn[idx])
    return '#%d' % idx


def _call_name(lqn, cidx: int) -> str:
    chn = getattr(lqn, 'callhashnames', None)
    if chn is not None:
        chn = np.asarray(chn).ravel()
        if cidx < len(chn) and str(chn[cidx]):
            return str(chn[cidx])
    return call_hashname(lqn, cidx)


def _term_name(lqn, kind: str, idx: int) -> str:
    return _elem_name(lqn, idx) if kind == 'entry' else _call_name(lqn, idx)
