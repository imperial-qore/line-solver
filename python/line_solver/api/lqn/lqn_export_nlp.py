"""
Export a layered queueing network as a standalone nonlinear program.

The balance equations of ``lqn_balance_equations`` say what a layered model's
throughputs, utilizations and think times must satisfy; they do not say what the
model's performance IS, because the service times they are stated over are
themselves outputs. This module closes that gap ALGEBRAICALLY, with no solver
and no layer oracle, by adding the two congestion laws that a mean-value
argument supplies, and emits the result as one self-contained Python script.

The split the emitted program is built on:

    LINEAR CONSTRAINTS
        call flow            X(c) = y(c)*X(a)
        entry flow           X(e) = sum_c X(c)
        activity flow        X(a) = v(a)*X(e)
        entry composition    S(e) = sum_a v(a)*[W(a) + sum_c y(c)*R(c)]
        delay residence      W(a) = D(a) and R(c) = S(e) at an infinite server
        reference closure    z(t)*X(e_t) + B(t) = N(t)
        capacity             sum_a D(a)*X(a) <= m(h),  sum_c B(c) <= N(t)

    NONLINEAR OBJECTIVE (sum of squares, zero at the fixed point)
        Little at a host     Q(a)  = X(a)*W(a)
        Little at a task     Qt(c) = X(c)*R(c)
        busy threads         B(c)  = X(c)*S(e)
        host congestion      W(a)*G(h) = D(a)*(1 + Qbar(h,a))
        task congestion      R(c)*G(t) = S(e)*(1 + Qbar(t,c))
        queue-dependent rate G(i) = r_i(1 + delta_i*sum_j Q_j)

The congestion laws are QD-AMVA (Casale, Perez and Wang, PERFORMANCE 2015),
including the arrival-instant estimate they are stated over,

    Qbar(h,a) = sum_a' Q(a') - (1/N) sum_{a' in task(a)} Q(a'),

linear in the queue lengths. The class an arrival belongs to is its TASK, not
its activity: a thread is one customer whichever of its task's activities it is
running, so the share the estimate removes is the task's whole share of that
station. What the queue-dependence adds is that a station is
not served at a fixed rate but at r_i(n), its own rate when n jobs are present.
Then a MULTISERVER is not a special case and needs no clamp: r_i(n) =
softmin(n, m_i) is one job's worth of service while the pool has spare servers
and m_i once it is full, so a lightly loaded pool charges the bare demand, a
saturated one charges m_i-way sharing, and nothing in between has to be decided
by a branch. Any further load-dependent multiplier declared for the station
multiplies into r_i next to it; see the ``lld`` argument of :func:`export_nlp`.
Delay stations keep the exact law W = D and stay linear.

The price is that the host law is no longer linear. A rate r_i(n) evaluated at a
queue length is not an affine function of the iterate whatever its shape, so the
congestion laws at BOTH levels are now residuals, carrying one variable G(i) per
queueing station for the rate at the arrival instant. Each such law is a plain
bilinear residual W*G - D*(1+Qbar), and only the definition of G is a general
smooth term.

Counting: after the linear block the program has exactly as many free
dimensions as the objective has residuals, so its zeros are isolated rather
than a manifold, and a zero objective is an exact solution of the approximate
laws above. This is an APPROXIMATION by construction: QD-AMVA, not exact MVA and
not Linearizer, and furthest off
where a station sits at saturation. What it is not is a black box -- the
emitted script carries the model's parameters as named literals, writes each
station's rate function out in full, and builds every matrix from them in the
open.

THE AND-FORK CLOSURE (``fork='exp'``, off by default)
-----------------------------------------------------
The entry composition law above is a SUM, and a sum is exact however the graph
branches or loops as long as its activities run one at a time: expectation is
linear along a path, and ``v(a)`` absorbs the branching. An AND-fork breaks
that, because its branches run CONCURRENTLY and the entry is released by the
LAST of them, not by all of them in turn.

``max`` is not linear, so the law cannot stay in the equality block. Worse,
``E[max]`` IS NOT A FUNCTION OF THE BRANCH MEANS AT ALL, and means are all this
program carries: two branches of mean D give 1.5*D if they are exponential and
D if they are deterministic, against the 2*D a sum would charge. So no value
put in that row is right for more than one assumption about branch shape, and
``fork='exp'`` states the assumption instead of hiding it: each branch
completion is taken as exponential, and

    E[max] = sum over nonempty S of (-1)^(|S|+1) / sum_{i in S} 1/m(i)

is a rational, analytic, degree-one-homogeneous function of the branch means,
carried as ONE RESIDUAL per forking entry. The linear block loses that entry's
equality and the objective gains a residual, so the counting above is
preserved. This is the closure ``fj_quorum_moments`` already applies inside
SolverLN; the difference is that here it sits INSIDE the fixed point rather
than being applied to ``entry_servt`` after the sweep.

WHAT IT DOES NOT FIX. A forking entry holds SEVERAL THREADS of its task at
once, and the thread-seconds it occupies stay the additive sum even as the
elapsed time collapses to the max. The two are therefore split: ``S(e)`` is the
elapsed hold, which is what a caller waits for and what the reference closure
counts, and ``SA(e)`` is the additive thread-seconds, which is what ``B(c)``
and the thread-pool capacity bound are built from. That makes occupancy right.
What stays approximate is the task congestion law ``R(c)*G(t) = S(e)*(1 +
Qbar(t,c))``: it is a multiserver law over the pool, and it prices a request at
one server for its elapsed time, where a forking entry in fact ties up several
at once. Simultaneous resource possession is outside QD-AMVA, so a task whose
entries fork is served right on time and approximately on thread contention.

A QUORUM join, which fires on k of n branches, is refused: the k-th order
statistic is a different closure and is not this one. So is a branch that
carries neither host demand nor a call, whose completion time is identically
zero, and so is a fork nested inside another fork's branch.

The emitted ``solve()`` runs Newton on the square system first, since the
minimum wanted is a root. Its starting point is the empty network, which is the
wrong basin for a model that runs congested, so it escalates the number of
fixed-point sweeps taken first (0, 5, 20, 100) before falling back to a
constrained descent (``trust-constr``): a sweep costs one pass over the call
graph and the descent costs orders of magnitude more. The Jacobian stays
analytic across the rate functions: they are differentiated by complex step,
which is exact to machine precision for anything smooth.

Features the algebraic form cannot express are refused BY NAME rather than
silently dropped; see :class:`LqnNlpExportError`.
"""

from typing import Dict, List, Optional

import numpy as np

from ...constants import SchedStrategy
from .lqn_balance_equations import (_act_visits, _call_name, _calltype, _elem_name,
                                    _entry_of_activity, _graph, _hostdem, _isref,
                                    _maxmult, _mult, _parent, _ref_thinktime, _repl,
                                    _sched_of, _SYNC)

__all__ = ['LqnNlpExportError', 'export_nlp']


#: precedence ids the struct builder writes; they are shared with MATLAB and
#: the JAR. PRE_AND marks the JOINED PREDECESSORS of a join, POST_AND the
#: BRANCH HEADS of a fork, and ``actquorum`` sits on the join TARGET.
_PRE_AND = 2
_POST_AND = 12

#: what ``export_nlp(fork=...)`` accepts
FORK_MODES = ('refuse', 'exp')


class LqnNlpExportError(ValueError):
    """
    A model feature the algebraic NLP form cannot express.

    Raised with every offending element named. The alternative to refusing is
    emitting a script that runs, reports numbers and is quietly wrong about the
    feature it dropped, which is worse than no script at all.
    """


# ----------------------------------------------------------------------
# public API
# ----------------------------------------------------------------------

def export_nlp(lqn, filename: Optional[str] = None, name: Optional[str] = None,
               lld: Optional[Dict[str, str]] = None, alpha: float = 20.0,
               fork: str = 'refuse') -> str:
    """
    Emit a standalone Python script that approximates LQN by a nonlinear program.

    Args:
        lqn: a LayeredNetwork, or its LayeredNetworkStruct
        filename: where to write the script; the text is returned either way
        name: model name to record in the script, defaults to the model's own
        lld: queue-dependent rate multipliers, ``{station name: expression}``.
            The station is a processor or a task, and the expression is Python
            source in the free variable ``n``, the number of jobs present. It
            multiplies the multiserver rate ``softmin(n, m)``, so ``'1.0'`` is
            the default, ``'1 + 0.1*n'`` a station that speeds up under load and
            ``'2 - n/10'`` one that slows down. It MUST BE SMOOTH: the Jacobian
            differentiates it by complex step, so it has to be written in
            operations that accept a complex argument. A tabulated
            ``lldscaling`` row is refused rather than interpolated, because the
            piecewise-linear interpolant of one is not differentiable at its
            knots and would cost the Newton step its quadratic convergence.
        alpha: sharpness of the softmin standing in for ``min(n, m)``. LINE's
            own QD-AMVA uses 20, which places the multiserver knee within about
            1e-9 of the corner one job's worth away from it.
        fork: how to treat an AND-fork/join block. ``'refuse'``, the default,
            declines the model by name. ``'exp'`` carries it under the
            EXPONENTIAL-BRANCH closure: the entry is held for E[max] of its
            concurrent branches, taken as if each branch completion were
            exponential with the mean the program solves for. That is an
            assumption about branch SHAPE, which the program does not otherwise
            carry, and it is the same one ``fj_quorum_moments`` makes inside
            SolverLN; read ``THE AND-FORK CLOSURE`` below before using it.

    Returns:
        the script source

    Raises:
        LqnNlpExportError: the model uses a feature the form cannot express
    """
    if hasattr(lqn, 'getStruct'):
        # A struct carries no name of its own, so take the model's while it is here.
        name = name or str(getattr(lqn, 'name', '') or '')
        lqn = lqn.getStruct()
    if fork not in FORK_MODES:
        raise LqnNlpExportError('fork takes %s, got %r'
                                % (' or '.join(repr(m) for m in FORK_MODES), fork))
    p = _collect(lqn, name, fork)
    _attach_lld(p, lld, alpha)
    text = _emit(p)
    if filename:
        with open(filename, 'w') as fh:
            fh.write(text)
    return text


# ----------------------------------------------------------------------
# scope
# ----------------------------------------------------------------------

def _scope_check(lqn, fork: str) -> List[Dict]:
    bad: List[str] = []
    nidx = int(lqn.nidx)
    tshift, ashift = int(lqn.tshift), int(lqn.ashift)

    for cidx in range(int(lqn.ncalls)):
        ct = _calltype(lqn, cidx)
        if ct != _SYNC:
            kind = 'asynchronous' if ct == 2 else 'forwarding'
            bad.append('%s call: %s' % (kind, _call_name(lqn, cidx)))

    arv = getattr(lqn, 'arrival', None) or {}
    for eidx in sorted(arv):
        bad.append('open arrival at entry %s' % _elem_name(lqn, int(eidx)))

    for idx in range(nidx):
        if _repl(lqn, idx) > 1:
            bad.append('replication of %s (repl=%g)' % (_elem_name(lqn, idx), _repl(lqn, idx)))

    ph = getattr(lqn, 'actphase', None)
    if ph is not None:
        ph = np.asarray(ph).ravel()
        for a in range(int(lqn.nacts)):
            # actphase is 0-based over activities, as every band is now
            if a < len(ph) and float(ph[a]) > 1:
                bad.append('phase-2 activity %s' % _elem_name(lqn, ashift + a))

    hs = getattr(lqn, 'hassetup', None)
    if hs is not None:
        hs = np.asarray(hs).ravel()
        for t in range(int(lqn.ntasks)):
            if tshift + t < len(hs) and hs[tshift + t]:
                bad.append('setup time on task %s' % _elem_name(lqn, tshift + t))

    ic = getattr(lqn, 'iscache', None)
    if ic is not None:
        ic = np.asarray(ic).ravel()
        for t in range(int(lqn.ntasks)):
            if tshift + t < len(ic) and ic[tshift + t]:
                bad.append('cache task %s' % _elem_name(lqn, tshift + t))

    # Concurrency. Under fork='refuse' every AND marking is named and declined;
    # under fork='exp' a well-formed block is carried and only what cannot be
    # read as one is named. See THE AND-FORK CLOSURE in the module docstring.
    blocks, forkbad = _fork_blocks(lqn, fork)
    bad += forkbad

    lc = getattr(lqn, 'lincon', None) or {}
    for idx in sorted(lc):
        bad.append('admission region on %s' % _elem_name(lqn, int(idx)))

    # An activity in two entry blocks has no single entry throughput to be a
    # multiple of, so its flow relation is not expressible here.
    for a in range(int(lqn.nacts)):
        aidx = ashift + a
        owners = [int(e) for e in lqn.entriesof.get(_parent(lqn, aidx), [])
                  if aidx in [int(x) for x in lqn.actsof.get(int(e), [])]]
        if len(owners) > 1:
            bad.append('activity %s shared by entries %s'
                       % (_elem_name(lqn, aidx),
                          ', '.join(_elem_name(lqn, e) for e in owners)))
        elif not owners:
            bad.append('activity %s bound to no entry' % _elem_name(lqn, aidx))

    refs = [tshift + t for t in range(int(lqn.ntasks)) if _isref(lqn, tshift + t)]
    if not refs:
        bad.append('no reference task: nothing drives the model')
    for tidx in refs:
        ents = [int(e) for e in lqn.entriesof.get(tidx, [])]
        if len(ents) != 1:
            bad.append('reference task %s has %d entries, expected 1'
                       % (_elem_name(lqn, tidx), len(ents)))

    if bad:
        raise LqnNlpExportError(
            'cannot express this model as linear constraints plus a nonlinear '
            'objective:\n  ' + '\n  '.join(bad))
    return blocks


# ----------------------------------------------------------------------
# AND-fork blocks
# ----------------------------------------------------------------------

def _vecf(lqn, name: str) -> np.ndarray:
    v = getattr(lqn, name, None)
    return np.zeros(0) if v is None else np.asarray(v).ravel()


def _atv(v: np.ndarray, i: int) -> float:
    return float(v[i]) if 0 <= i < len(v) else 0.0


def _fork_blocks(lqn, fork: str):
    """
    The AND-fork/join blocks of the model, and every AND marking outside one.

    Returns ``(blocks, bad)``. A block is a dict carrying the entry it runs in,
    the join target, and the branches as lists of ABSOLUTE activity indices;
    ``bad`` names what could not be read as a well-formed block, for the scope
    check to refuse. Under ``fork='refuse'`` every marking lands in ``bad``.

    The struct marks a fork on its BRANCH HEADS (``actposttype`` = POST_AND) and
    a join on its JOINED PREDECESSORS (``actpretype`` = PRE_AND), with the
    quorum on the join TARGET, so the block is recovered by grouping heads that
    share a predecessor and following each to where the branches reconverge.
    """
    ashift, nacts = int(lqn.ashift), int(lqn.nacts)
    post, pre, quo = (_vecf(lqn, 'actposttype'), _vecf(lqn, 'actpretype'),
                      _vecf(lqn, 'actquorum'))
    allacts = set(range(ashift, ashift + nacts))
    marked = sorted(a for a in allacts
                    if int(_atv(post, a)) == _POST_AND or int(_atv(pre, a)) == _PRE_AND)
    if not marked:
        return [], []
    if fork != 'exp':
        return [], ['AND-fork or AND-join at activity %s (pass fork=\'exp\' to '
                    'carry it under the exponential-branch closure)' % _elem_name(lqn, a)
                    for a in marked]

    G = _graph(lqn)
    blocks: List[Dict] = []
    bad: List[str] = []
    placed = set()

    for tidx in [int(lqn.tshift) + t for t in range(int(lqn.ntasks))]:
        for eidx in sorted(int(e) for e in lqn.entriesof.get(tidx, [])):
            A = sorted(set(int(a) for a in lqn.actsof.get(eidx, [])) & allacts)
            if not A:
                continue
            succ = dict((a, [b for b in A if G[a, b] != 0]) for a in A)
            prd = dict((a, [b for b in A if G[b, a] != 0]) for a in A)

            def reach(start):
                seen, stack = set(), [start]
                while stack:
                    x = stack.pop()
                    if x in seen:
                        continue
                    seen.add(x)
                    stack += succ[x]
                return seen

            groups: Dict = {}
            for h in [a for a in A if int(_atv(post, a)) == _POST_AND]:
                groups.setdefault(tuple(sorted(prd[h])), []).append(h)
            for key, H in sorted(groups.items()):
                blk, why = _read_block(lqn, A, H, key, reach, prd, post, pre, quo)
                if why:
                    bad.append(why)
                    continue
                blk['entry'] = eidx
                blk['label'] = '%s -> %s' % (
                    '+'.join(_elem_name(lqn, a) for a in key) or _elem_name(lqn, eidx),
                    _elem_name(lqn, blk['join']))
                blocks.append(blk)
                placed |= set(H) | set(blk['tails'])

    for a in marked:
        if a not in placed:
            bad.append('AND marking at activity %s outside any well-formed '
                       'fork/join block' % _elem_name(lqn, a))
    return blocks, bad


def _read_block(lqn, A, H, srcs, reach, prd, post, pre, quo):
    """One candidate fork/join block, or the reason it cannot be read as one."""

    def nm(a):
        return _elem_name(lqn, a)

    lbl = '+'.join(nm(h) for h in sorted(H))
    if len(H) < 2:
        return None, 'AND-fork at %s has one branch' % lbl

    Rs = [reach(h) for h in H]
    common = set(Rs[0])
    for r in Rs[1:]:
        common &= r
    if not common:
        return None, ('the AND-fork branches at %s never reconverge, so the entry '
                      'has no single completion' % lbl)
    cands = [x for x in sorted(common)
             if int(_atv(quo, x)) > 0 and common <= reach(x)]
    if len(cands) != 1:
        return None, ('the AND-fork branches at %s do not meet at one join '
                      'target' % lbl)
    j = cands[0]

    tails = sorted(a for a in A if int(_atv(pre, a)) == _PRE_AND and a in prd[j])
    if len(tails) != len(H):
        return None, ('the join %s takes %d AND predecessors where the fork at %s '
                      'opens %d branches' % (nm(j), len(tails), lbl, len(H)))
    q = int(_atv(quo, j))
    if q != len(H):
        return None, ('quorum join at %s fires on %d of %d branches; the k-th '
                      'order statistic is a different closure' % (nm(j), q, len(H)))

    branches = [sorted(r - common) for r in Rs]
    seen = set()
    for i, br in enumerate(branches):
        if not br:
            return None, 'branch %s of the AND-fork at %s is empty' % (nm(H[i]), lbl)
        if seen & set(br):
            return None, ('the AND-fork branches at %s overlap, so they are not '
                          'independent' % lbl)
        seen |= set(br)
    for i, br in enumerate(branches):
        for a in br:
            for b in prd[a]:
                if b not in br and b not in srcs:
                    return None, ('activity %s is inside a branch of the AND-fork at '
                                  '%s but is also entered from %s, so the branch is '
                                  'not a block' % (nm(a), lbl, nm(b)))
            if int(_atv(post, a)) == _POST_AND and a != H[i]:
                return None, ('AND-fork at %s nested inside a branch of the one at '
                              '%s' % (nm(a), lbl))
            if int(_atv(pre, a)) == _PRE_AND and a not in tails:
                return None, ('AND-join at %s nested inside a branch of the fork at '
                              '%s' % (nm(a), lbl))
        if not any(_hostdem(lqn, a) > 0 or any(int(lqn.callpair[c, 0]) == a
                                               for c in range(int(lqn.ncalls)))
                   for a in br):
            return None, ('branch %s of the AND-fork at %s has no host demand and '
                          'no call, so its completion time is identically zero'
                          % (nm(H[i]), lbl))
    return {'join': j, 'tails': tails, 'branches': branches}, ''


# ----------------------------------------------------------------------
# parameter collection
# ----------------------------------------------------------------------

def _collect(lqn, name: Optional[str], fork: str = 'refuse') -> Dict:
    blocks = _scope_check(lqn, fork)

    nhosts, ntasks = int(lqn.nhosts), int(lqn.ntasks)
    nentries, nacts, ncalls = int(lqn.nentries), int(lqn.nacts), int(lqn.ncalls)
    hshift = int(getattr(lqn, 'hshift', 0) or 0)
    tshift, eshift, ashift = int(lqn.tshift), int(lqn.eshift), int(lqn.ashift)

    hidx = [hshift + h for h in range(nhosts)]
    tidx = [tshift + t for t in range(ntasks)]
    eidx = [eshift + e for e in range(nentries)]
    aidx = [ashift + a for a in range(nacts)]

    tpos = dict((x, i) for i, x in enumerate(tidx))
    epos = dict((x, i) for i, x in enumerate(eidx))
    apos = dict((x, i) for i, x in enumerate(aidx))
    hpos = dict((x, i) for i, x in enumerate(hidx))

    visits = _act_visits(lqn)

    p: Dict = {}
    p['name'] = name or str(getattr(lqn, 'name', '') or 'lqn')
    p['hosts'] = [_elem_name(lqn, i) for i in hidx]
    p['host_mult'] = [_mult(lqn, i) for i in hidx]
    p['host_maxmult'] = [_eff_mult(lqn, i) for i in hidx]
    p['host_inf'] = [_is_inf(lqn, i) for i in hidx]

    p['tasks'] = [_elem_name(lqn, i) for i in tidx]
    p['task_host'] = [hpos[_parent(lqn, i)] for i in tidx]
    p['task_mult'] = [_mult(lqn, i) for i in tidx]
    p['task_inf'] = [_is_inf(lqn, i) for i in tidx]
    p['task_isref'] = [bool(_isref(lqn, i)) for i in tidx]
    p['task_think'] = [_ref_thinktime(lqn, i) for i in tidx]

    p['entries'] = [_elem_name(lqn, i) for i in eidx]
    p['entry_task'] = [tpos[_parent(lqn, i)] for i in eidx]

    p['acts'] = [_elem_name(lqn, i) for i in aidx]
    p['act_task'] = [tpos[_parent(lqn, i)] for i in aidx]
    p['act_entry'] = [epos[_entry_of_activity(lqn, i)] for i in aidx]
    p['act_demand'] = [_hostdem(lqn, i) for i in aidx]
    p['act_visits'] = [float(visits[i]) for i in aidx]

    p['calls'] = [_call_name(lqn, c) for c in range(ncalls)]
    p['call_act'] = [apos[int(lqn.callpair[c, 0])] for c in range(ncalls)]
    p['call_entry'] = [epos[int(lqn.callpair[c, 1])] for c in range(ncalls)]
    p['call_mean'] = [_call_mean(lqn, c) for c in range(ncalls)]
    p['call_caller_task'] = [p['act_task'][a] for a in p['call_act']]

    # AND-fork blocks, in the program's own 0-based positions. An entry named
    # here composes by E[max] over its branches instead of by a plain sum.
    p['fork_entry'] = [epos[b['entry']] for b in blocks]
    p['fork_branches'] = [[[apos[a] for a in br] for br in b['branches']]
                          for b in blocks]
    p['fork_label'] = [b['label'] for b in blocks]

    p['entry_order'] = _entry_order(p)
    return p


def _is_inf(lqn, idx: int) -> bool:
    """A delay station: an infinite server, or an unbounded multiplicity."""
    m = _mult(lqn, idx)
    return bool(np.isinf(m)) or _sched_of(lqn, idx) == SchedStrategy.INF


def _eff_mult(lqn, idx: int) -> float:
    """Servers that can ACTUALLY be busy, which is what a utilization divides by.

    `maxmult` is `lsn_max_multiplicity`: the reference-task populations pushed
    along the call DAG, each node keeping `min(inflow, mult)`, so a pool wider
    than the work reaching it is cut down to the work. SolverLN hands exactly
    this to the layer's station (`_get_nservers`), so dividing by anything else
    reports a different quantity from every other layered method. It is absent
    on a struct built without it and 0 at a delay, where all capacity is usable;
    both fall back to the declared count, as the layer construction does.
    """
    mm = _maxmult(lqn, idx)
    return _mult(lqn, idx) if (np.isnan(mm) or mm <= 0) else mm


def _call_mean(lqn, cidx: int) -> float:
    """Mean number of calls y(c). A missing or non-finite count is not a rate."""
    y = float(lqn.callpair[cidx, 2]) if lqn.callpair.shape[1] > 2 else 1.0
    if not np.isfinite(y) or y < 0:
        raise LqnNlpExportError('call %s has a non-finite mean count %r'
                                % (_call_name(lqn, cidx), y))
    return y


def _entry_order(p: Dict) -> List[int]:
    """
    Entries with every callee before its caller.

    The entry composition law reads the response times of the calls an entry
    dispatches, so a topological order over the call graph is what lets the
    emitted script build a starting point in one pass instead of iterating. A
    cycle here is a call-path deadlock in the model, not a numerical issue.
    """
    ne = len(p['entries'])
    dep = [set() for _ in range(ne)]     # dep[e] = entries e must follow
    for c in range(len(p['calls'])):
        caller = p['act_entry'][p['call_act'][c]]
        dep[caller].add(p['call_entry'][c])
    order, done = [], [False] * ne
    while len(order) < ne:
        progressed = False
        for e in range(ne):
            if not done[e] and all(done[d] for d in dep[e]):
                done[e] = True
                order.append(e)
                progressed = True
        if not progressed:
            stuck = [p['entries'][e] for e in range(ne) if not done[e]]
            raise LqnNlpExportError('cyclic call graph over entries: %s' % ', '.join(stuck))
    return order


# ----------------------------------------------------------------------
# queue-dependent rates
# ----------------------------------------------------------------------

# Shared verbatim between the validation below and the emitted script, so the
# expression an LLD is checked against is the one it will be evaluated by.
_SOFTMIN_SRC = r'''
def softmin(a, b, alpha=SOFTMIN_ALPHA):
    # Smooth stand-in for min(a,b): (a*e^-alpha*a + b*e^-alpha*b) over the same
    # two exponentials, factored through the smaller argument so nothing
    # overflows. Written to accept a complex argument, since `drate`
    # differentiates the rate functions through it by complex step.
    lo, hi = (a, b) if np.real(a) <= np.real(b) else (b, a)
    gap = hi - lo
    if np.real(gap) * alpha > 700.0:
        return lo                       # the exponential has underflowed: min
    w = np.exp(-alpha * gap)
    return lo + gap * w / (1.0 + w)
'''


def _lld_env(alpha: float) -> Dict:
    env = {'np': np, 'SOFTMIN_ALPHA': float(alpha)}
    exec(compile(_SOFTMIN_SRC, '<softmin>', 'exec'), env)
    return env


def _attach_lld(p: Dict, lld: Optional[Dict[str, str]], alpha: float) -> None:
    """Resolve the declared rate multipliers onto the stations they name."""
    alpha = float(alpha)
    if not np.isfinite(alpha) or alpha <= 0:
        raise LqnNlpExportError('the softmin sharpness alpha must be finite and '
                                'positive, got %r' % alpha)
    p['alpha'] = alpha
    p['host_lld'] = ['1.0'] * len(p['hosts'])
    p['task_lld'] = ['1.0'] * len(p['tasks'])
    if not lld:
        return

    env = _lld_env(alpha)
    bad: List[str] = []
    for key in sorted(lld):
        src = lld[key]
        if isinstance(src, (list, tuple, np.ndarray)):
            bad.append('%s: a tabulated lldscaling row is not smooth -- its interpolant '
                       'has a corner at every knot. Pass an expression in n instead.' % key)
            continue
        if not isinstance(src, str):
            bad.append('%s: a queue-dependent rate is an expression in n, not %s'
                       % (key, type(src).__name__))
            continue
        where = [('host', i) for i, nm in enumerate(p['hosts']) if nm == key]
        where += [('task', i) for i, nm in enumerate(p['tasks']) if nm == key]
        if not where:
            bad.append('%s: no processor or task by that name' % key)
            continue
        if len(where) > 1:
            bad.append('%s: names both a processor and a task' % key)
            continue
        why = _check_lld(src, env)
        if why:
            bad.append('%s: %s' % (key, why))
            continue
        kind, i = where[0]
        p['host_lld' if kind == 'host' else 'task_lld'][i] = src

    if bad:
        raise LqnNlpExportError('cannot use these queue-dependent rates:\n  '
                                + '\n  '.join(bad))


def _check_lld(src: str, env: Dict) -> str:
    """
    What the emitted program needs of an LLD expression, checked here.

    Positive and real on the range of populations a station can hold, and
    ANALYTIC, since the Jacobian differentiates it by complex step. The last is
    the one that is easy to get wrong by accident: abs, min, max, a comparison
    or an interpolation all evaluate perfectly well and all silently return a
    complex step of zero or nonsense, which would leave the Newton step wrong
    rather than merely slow. Comparing the complex step against a central
    difference catches every one of them.
    """
    try:
        code = compile(src, '<lld>', 'eval')
    except SyntaxError as exc:
        return 'not a Python expression (%s)' % exc.msg

    def at(n):
        return eval(code, dict(env), {'n': n})

    for n in (1.0, 2.5, 10.0):
        try:
            v = at(n)
        except Exception as exc:
            return 'fails at n=%g (%s: %s)' % (n, type(exc).__name__, exc)
        if np.imag(v) != 0 or not np.isfinite(np.real(v)) or float(np.real(v)) <= 0:
            return 'must be real and positive, is %r at n=%g' % (v, n)

    n0, step = 2.5, 1e-6
    try:
        dcs = float(np.imag(at(complex(n0, 1e-20)))) / 1e-20
    except Exception as exc:
        return 'rejects a complex argument, so it cannot be differentiated by ' \
               'complex step (%s: %s)' % (type(exc).__name__, exc)
    dfd = (float(np.real(at(n0 + step))) - float(np.real(at(n0 - step)))) / (2 * step)
    if abs(dcs - dfd) > 1e-4 * (1.0 + abs(dfd)):
        return ('is not analytic at n=%g: its complex step gives %g where a finite '
                'difference gives %g. abs, min, max, a comparison or an interpolation '
                'will do this.' % (n0, dcs, dfd))
    return ''


# ----------------------------------------------------------------------
# emission
# ----------------------------------------------------------------------

_Q3 = '"' * 3


def _lit(x) -> str:
    if isinstance(x, bool):
        return 'True' if x else 'False'
    if isinstance(x, float):
        if np.isinf(x):
            return "float('inf')" if x > 0 else "-float('inf')"
        return repr(float(x))
    if isinstance(x, str):
        return repr(x)
    return repr(x)


def _arr(nm: str, vals, comment: str = '') -> str:
    body = '[' + ', '.join(_lit(v) for v in vals) + ']'
    line = '%s = %s' % (nm, body)
    if comment:
        line += '  # ' + comment
    return line


def _fork_header(p: Dict) -> List[str]:
    """The concurrency paragraph, present only when the model actually forks."""
    if not p['fork_entry']:
        return []
    return [
        '',
        'THIS MODEL FORKS. %d AND-fork block%s run branches CONCURRENTLY, so the'
        % (len(p['fork_entry']), '' if len(p['fork_entry']) == 1 else 's'),
        'entries holding them are released by the LAST branch and not by all of them',
        'in turn. Their composition law is therefore not the sum above but',
        '',
        '  entry service   S(e) = sum_a v(a)*[...] + sum_blocks E[max of branches]',
        '  thread-seconds  SA(e) = sum_a v(a)*[...]   over every branch, additively',
        '  busy threads    B(c) = X(c)*SA(e)',
        '',
        'E[MAX] IS NOT A FUNCTION OF THE BRANCH MEANS, and means are all this program',
        'carries: two branches of mean D give 1.5*D if exponential and D if',
        'deterministic, against the 2*D a sum would charge. fj_max_exp below takes',
        'the EXPONENTIAL case, which is the closure LINE applies elsewhere, and says',
        'so rather than hiding it. Edit that one function to assume something else.',
        '',
        'S(e) and SA(e) are split because a forking entry HOLDS SEVERAL THREADS AT',
        'ONCE. S is the elapsed hold a caller waits for and the reference closure',
        'counts; SA is the thread-seconds its pool is occupied for, which stay',
        'additive. What is still approximate is the task congestion law, which',
        'prices a request at one thread for its elapsed time where a fork ties up',
        'several: simultaneous resource possession is outside QD-AMVA. A task whose',
        'entries fork is served right on time and approximately on contention.',
    ]


def _header(p: Dict) -> List[str]:
    return [
        _Q3,
        'Nonlinear program approximating the layered queueing network %r.' % p['name'],
        '',
        'Generated by line_solver.api.lqn.export_nlp. Self-contained: numpy and',
        'scipy only, with the model parameters below as literals.',
        '',
        'The program is a set of LINEAR EQUALITY AND INEQUALITY CONSTRAINTS over the',
        'flows, the service-time composition and the closures, plus a NONLINEAR',
        'OBJECTIVE that is the sum of squares of the relations the linear block cannot',
        'carry. Its global minimum is zero, attained where',
        '',
        '  flow            X(c) = y(c)*X(a),  X(e) = sum_c X(c),  X(a) = v(a)*X(e)',
        '  entry service   S(e) = sum_a v(a)*[W(a) + sum_c y(c)*R(c)]',
        '  host residence  W(a)*G(h) = D(a)*(1 + Qbar(h,a))',
        '  task residence  R(c)*G(t) = S(e)*(1 + Qbar(t,c))',
        '  station rate    G(i) = r_i(1 + delta_i*sum_j Q_j)',
        '  queue lengths   Q(a) = X(a)*W(a),  Qt(c) = X(c)*R(c)',
        '  busy threads    B(c) = X(c)*S(e),  z(t)*X(e_t) + B(t) = N(t)',
        '',
        'This is QD-AMVA. Qbar is its arrival-instant queue-length estimate,',
        'Qbar(h,a) = sum_a. Q(a.) - (1/N) sum over the activities of the arriving',
        "job's own task, since a thread is one customer whichever of its task's",
        'activities it is running. r_i(n) is the rate station i serves',
        'at while n jobs are present: HOST_RATE and TASK_RATE below write one out per',
        'station, in full, for you to read or edit.',
        '',
        'A MULTISERVER needs no special case and no clamp. Its rate is',
        'softmin(n, m), which is n while the pool has spare servers and m once it is',
        'full, so W = D*(1+Qbar)/softmin(1+Qbar, m) charges the bare demand at light',
        'load and m-way sharing at saturation, with nothing in between decided by a',
        'branch. Any load-dependent multiplier declared for the station multiplies',
        'into its rate alongside the softmin. Delay stations keep the exact W = D.',
        '',
        'This is an APPROXIMATION: QD-AMVA. It is',
        'not exact MVA and does not reproduce SolverLN to the last digit; it',
        'reproduces its structure, and is furthest off where a station is near',
        'saturation.',
    ] + _fork_header(p) + [
        '',
        'Run it directly to solve and print the metrics:',
        '',
        '    python %s.py              solve and report' % _slug(p['name']),
        '    python %s.py --warm 50    start from 50 fixed-point sweeps' % _slug(p['name']),
        '    python %s.py --verbose    show the optimizer trace' % _slug(p['name']),
        _Q3,
        '',
        'import sys',
        '',
        'import numpy as np',
        'from scipy.optimize import Bounds, LinearConstraint, minimize',
        '',
    ]


def _slug(name: str) -> str:
    return ''.join(ch if (ch.isalnum() or ch == '_') else '_' for ch in str(name)) or 'lqn'


def _params(p: Dict) -> List[str]:
    L = ['# ' + '-' * 70,
         '# model parameters, as declared in the LQN',
         '# ' + '-' * 70,
         '',
         'MODEL = %s' % _lit(p['name']),
         '',
         '# processors']
    L.append(_arr('HOSTS', p['hosts']))
    L.append(_arr('HOST_MULT', p['host_mult'], 'servers at the processor'))
    L.append(_arr('HOST_MAXMULT', p['host_maxmult'],
                  'servers reachable by the offered work; the utilization divisor'))
    L.append(_arr('HOST_INF', p['host_inf'], 'infinite server, no queueing'))
    L += ['', '# tasks']
    L.append(_arr('TASKS', p['tasks']))
    L.append(_arr('TASK_HOST', p['task_host'], 'index into HOSTS'))
    L.append(_arr('TASK_MULT', p['task_mult'], 'threads, or population of a reference task'))
    L.append(_arr('TASK_INF', p['task_inf'], 'unbounded thread pool'))
    L.append(_arr('TASK_ISREF', p['task_isref']))
    L.append(_arr('TASK_THINK', p['task_think'], 'declared think time, reference tasks only'))
    L += ['', '# entries']
    L.append(_arr('ENTRIES', p['entries']))
    L.append(_arr('ENTRY_TASK', p['entry_task'], 'index into TASKS'))
    L.append(_arr('ENTRY_ORDER', p['entry_order'], 'callees before callers'))
    L += ['', '# activities']
    L.append(_arr('ACTS', p['acts']))
    L.append(_arr('ACT_TASK', p['act_task'], 'index into TASKS'))
    L.append(_arr('ACT_ENTRY', p['act_entry'], 'index into ENTRIES'))
    L.append(_arr('ACT_DEMAND', p['act_demand'], 'D(a), mean host demand per execution'))
    L.append(_arr('ACT_VISITS', p['act_visits'], 'v(a), executions per invocation of the entry'))
    L += ['', '# calls']
    L.append(_arr('CALLS', p['calls']))
    L.append(_arr('CALL_ACT', p['call_act'], 'dispatching activity, index into ACTS'))
    L.append(_arr('CALL_ENTRY', p['call_entry'], 'target entry, index into ENTRIES'))
    L.append(_arr('CALL_MEAN', p['call_mean'], 'y(c), mean calls per execution of the activity'))
    L.append(_arr('CALL_CALLER_TASK', p['call_caller_task'], 'index into TASKS'))
    L += ['', '# AND-fork blocks: the branches of each run CONCURRENTLY, so the entry',
          '# is held for the longest of them rather than for their sum']
    L.append(_arr('FORK_ENTRY', p['fork_entry'], 'index into ENTRIES'))
    L.append(_arr('FORK_BRANCHES', p['fork_branches'], 'per block, per branch, into ACTS'))
    L.append(_arr('FORK_LABEL', p['fork_label'], 'fork source -> join target'))
    L += _rates(p)
    return L


def _rate_lambda(mult: float, inf: bool, lld: str) -> str:
    if inf or np.isinf(mult):
        return 'lambda n: n'         # a delay serves every job at once
    body = 'softmin(n, %s)' % _lit(float(mult))
    if lld.strip() != '1.0':
        body += ' * (%s)' % lld.strip()
    return 'lambda n: ' + body


def _rate_note(nm: str, mult: float, inf: bool, lld: str) -> str:
    if inf or np.isinf(mult):
        return '%s: delay, never queues' % nm
    what = '%s: %g server%s' % (nm, mult, '' if mult == 1 else 's')
    if lld.strip() != '1.0':
        what += ', scaled by %s' % lld.strip()
    return what


def _rates(p: Dict) -> List[str]:
    L = ['',
         '# ' + '-' * 70,
         '# queue-dependent service rates (QD-AMVA)',
         '#',
         '# r_i(n) is the rate station i serves at while n jobs are present. The',
         '# softmin is the multiserver term: smooth everywhere, and min(n, m) to',
         '# within 1e-9 a job either side of the knee. Any factor after it is the',
         '# load-dependent multiplier declared for the station. Edit either and the',
         '# program follows, since the Jacobian differentiates these by complex step',
         '# rather than assuming a shape; anything smooth in n will do.',
         '# ' + '-' * 70,
         '',
         'SOFTMIN_ALPHA = %s' % _lit(float(p['alpha'])),
         '',
         _SOFTMIN_SRC.strip('\n'),
         '',
         _arr('HOST_RATE_CAPPED', [s.strip() == '1.0' for s in p['host_lld']],
              'rate never exceeds HOST_MULT, so the utilization bound holds'),
         '',
         'HOST_RATE = [']
    for i, nm in enumerate(p['hosts']):
        L.append('    %s,  # %s' % (_rate_lambda(p['host_mult'][i], p['host_inf'][i],
                                                 p['host_lld'][i]),
                                    _rate_note(nm, p['host_mult'][i], p['host_inf'][i],
                                               p['host_lld'][i])))
    L += [']', '', 'TASK_RATE = [']
    for i, nm in enumerate(p['tasks']):
        L.append('    %s,  # %s' % (_rate_lambda(p['task_mult'][i], p['task_inf'][i],
                                                 p['task_lld'][i]),
                                    _rate_note(nm, p['task_mult'][i], p['task_inf'][i],
                                               p['task_lld'][i])))
    L += [']', '']
    return L


_BODY = r'''
# ----------------------------------------------------------------------
# derived index sets
# ----------------------------------------------------------------------

INF = float('inf')
NH, NT, NE, NA, NC = len(HOSTS), len(TASKS), len(ENTRIES), len(ACTS), len(CALLS)

ACTS_OF_HOST = [[a for a in range(NA) if TASK_HOST[ACT_TASK[a]] == h] for h in range(NH)]
ACTS_OF_ENTRY = [[a for a in range(NA) if ACT_ENTRY[a] == e] for e in range(NE)]
CALLS_OF_ACT = [[c for c in range(NC) if CALL_ACT[c] == a] for a in range(NA)]
CALLS_INTO_ENTRY = [[c for c in range(NC) if CALL_ENTRY[c] == e] for e in range(NE)]
CALLS_INTO_TASK = [[c for c in range(NC) if ENTRY_TASK[CALL_ENTRY[c]] == t] for t in range(NT)]
ENTRIES_OF_TASK = [[e for e in range(NE) if ENTRY_TASK[e] == t] for t in range(NT)]

# AND-fork blocks. The branches of block b run CONCURRENTLY inside the entry
# FORK_ENTRY[b], so its activities are taken out of the entry's additive sum and
# priced by the closure below instead. Everything before the fork and from the
# join on stays sequential and stays in SEQ_ACTS_OF_ENTRY.
NFORK = len(FORK_ENTRY)
FORK_ACTS = set(a for blk in FORK_BRANCHES for br in blk for a in br)
BLOCKS_OF_ENTRY = [[b for b in range(NFORK) if FORK_ENTRY[b] == e] for e in range(NE)]
SEQ_ACTS_OF_ENTRY = [[a for a in ACTS_OF_ENTRY[e] if a not in FORK_ACTS]
                     for e in range(NE)]


def fj_max_exp(m):
    # E[max] of K independent EXPONENTIAL branch completion times of means m.
    # The max of exponentials is not exponential, but its mean is a rational
    # function of theirs:
    #
    #   E[max] = sum over nonempty S of (-1)^(|S|+1) / sum_{i in S} 1/m(i)
    #
    # written here with the reciprocals cleared, so a branch of zero mean
    # contributes zero instead of dividing by it. Two equal branches of mean D
    # give 1.5*D, against the 2*D a sequential sum charges and the D a max of
    # the means would. Homogeneous of degree one, so it does not matter whether
    # the branch times are read per fork execution or per entry invocation.
    #
    # THIS IS AN ASSUMPTION ABOUT BRANCH SHAPE, and the program carries no shape
    # otherwise: deterministic branches of the same means would give D. It is
    # the same closure fj_quorum_moments applies inside SolverLN.
    k = len(m)
    tot = 0.0
    for mask in range(1, 1 << k):
        S = [i for i in range(k) if (mask >> i) & 1]
        num, den = 1.0, 0.0
        for i in S:
            num = num * m[i]
            q = 1.0
            for j in S:
                if j != i:
                    q = q * m[j]
            den = den + q
        tot = tot + (1.0 if len(S) % 2 else -1.0) * num / den
    return tot


REF_TASKS = [t for t in range(NT) if TASK_ISREF[t]]
NR = len(REF_TASKS)
REF_ENTRY = [ENTRIES_OF_TASK[t][0] for t in REF_TASKS]

# A task is a delay if it is declared one or has no thread limit; a call into it
# never queues, which turns its congestion law from a residual into a constraint.
TASK_DELAY = [TASK_INF[t] or TASK_MULT[t] == INF for t in range(NT)]
HOST_DELAY = [HOST_INF[h] or HOST_MULT[h] == INF for h in range(NH)]

# Population of the chain a class belongs to, used only by the arrival-instant
# correction. A class of task t carries t's own multiplicity, which is what the
# layer decomposition gives it; an unbounded pool leaves the correction off.
TASK_POP = list(TASK_MULT)

# THE CLASS AT A STATION IS THE CALLER, NOT THE REQUEST. The arrival estimate
# removes the arriving job's own share of ITS OWN CLASS's queue, and at a host
# the class is the TASK whose threads are the customers: a thread running B0 is
# the same customer as the thread running B1, so the queue to correct by is the
# task's whole share of the station, spread over however many activities it runs
# there. Correcting by the single activity's queue instead charges a one-thread
# task for queueing behind itself, which is how a model whose reference task has
# population 1 came out at half the throughput it must have.
SIBLING_ACTS = [[ap for ap in ACTS_OF_HOST[TASK_HOST[ACT_TASK[a]]]
                 if ACT_TASK[ap] == ACT_TASK[a]] for a in range(NA)]
SIBLING_CALLS = [[cp for cp in CALLS_INTO_TASK[ENTRY_TASK[CALL_ENTRY[c]]]
                  if CALL_CALLER_TASK[cp] == CALL_CALLER_TASK[c]] for c in range(NC)]

# Every queueing station carries one variable for the rate it serves at when a
# request arrives, G(i) = r_i(1 + delta_i*sum_j Q_j). A delay station carries
# none: it never queues, so its residence law is a linear constraint instead.
QD_HOSTS = [h for h in range(NH) if not HOST_DELAY[h] and ACTS_OF_HOST[h]]
QD_TASKS = [t for t in range(NT) if not TASK_DELAY[t] and CALLS_INTO_TASK[t]]
GH_OF = dict((h, i) for i, h in enumerate(QD_HOSTS))
GT_OF = dict((t, i) for i, t in enumerate(QD_TASKS))
NGH, NGT = len(QD_HOSTS), len(QD_TASKS)


def _tot(pops):
    # The jobs that can ever be at a station at once. The classes of a HOST
    # layer are the tasks running on it; those of a TASK layer are the tasks
    # that call into it.
    tot = 0.0
    for n in pops:
        if n == INF:
            return INF
        tot += n
    return tot


def _delta(tot):
    # The arrival estimate removes the arriving job's own share of its class,
    # 1/N of it.
    # The rate is read at the station's TOTAL population, which belongs to no
    # one class, so the same correction is taken over the layer as a whole:
    # delta = (Nt-1)/Nt. An unbounded population leaves the correction off.
    if tot == INF:
        return 1.0
    return (tot - 1.0) / tot if tot > 0 else 0.0


HOST_NTOT = [_tot([TASK_POP[t] for t in sorted(set(ACT_TASK[a]
                                                   for a in ACTS_OF_HOST[h]))])
             for h in range(NH)]
TASK_NTOT = [_tot([TASK_POP[t] for t in sorted(set(CALL_CALLER_TASK[c]
                                                   for c in CALLS_INTO_TASK[t]))])
             for t in range(NT)]
HOST_DELTA = [_delta(n) for n in HOST_NTOT]
TASK_DELTA = [_delta(n) for n in TASK_NTOT]

# ----------------------------------------------------------------------
# variable layout
# ----------------------------------------------------------------------

# Xe(e)  entry throughput                 Q(a)   jobs of activity a at its host
# Xa(a)  activity throughput              Qt(c)  requests of call c at its task
# Xc(c)  call throughput                  BC(c)  threads of the callee held by c
# W(a)   host residence per execution     BR(r)  threads of reference task r held
# S(e)   entry service (thread hold) time GH(h)  arrival-instant rate of host h
# R(c)   call response time               GT(t)  arrival-instant rate of task t
#
# SA(e) appears ONLY when the model forks. S(e) is then the ELAPSED hold, the
# time a caller waits and the time the reference closure counts, while SA(e) is
# the additive THREAD-SECONDS the entry occupies of its pool: a forking entry
# holds several threads at once, so the two part company and each law takes the
# one it is about. Without a fork they are equal and S(e) carries both.
_BLOCKS = (('Xe', NE), ('Xa', NA), ('Xc', NC), ('W', NA), ('S', NE),
           ('R', NC), ('Q', NA), ('Qt', NC), ('BC', NC), ('BR', NR),
           ('GH', NGH), ('GT', NGT), ('SA', NE if NFORK else 0))

OFFSET = {}
_off = 0
for _nm, _n in _BLOCKS:
    OFFSET[_nm] = _off
    _off += _n
NVAR = _off


def col(block, i):
    # absolute column of variable BLOCK at index I
    return OFFSET[block] + i


def sent(e):
    # Column of the entry's THREAD-SECONDS per invocation, which is what a
    # thread pool is occupied for. Equal to the entry service time unless the
    # entry forks, so a model without concurrency carries no extra variable.
    return col('SA', e) if NFORK else col('S', e)


def var_names():
    out = []
    for nm, n in _BLOCKS:
        labels = {'Xe': ENTRIES, 'S': ENTRIES, 'SA': ENTRIES, 'Xa': ACTS,
                  'W': ACTS, 'Q': ACTS,
                  'Xc': CALLS, 'R': CALLS, 'Qt': CALLS, 'BC': CALLS,
                  'BR': [TASKS[t] for t in REF_TASKS],
                  'GH': [HOSTS[h] for h in QD_HOSTS],
                  'GT': [TASKS[t] for t in QD_TASKS]}[nm]
        out += ['%s[%s]' % (nm, labels[i]) for i in range(n)]
    return out


def qbar_host_terms(h, a, g):
    # g times the QD-AMVA arrival-instant estimate at host H as seen by an
    # arrival of activity A, sum_a' Q(a') - Qclass(a)/N, as (column,
    # coefficient) pairs: linear in Q, and used as one factor of a bilinear
    # residual term. Qclass is the queue of the arriving job's own CLASS, which
    # is its task -- see SIBLING_ACTS.
    out = [(col('Q', ap), g) for ap in ACTS_OF_HOST[h]]
    pop = TASK_POP[ACT_TASK[a]]
    if pop != INF and pop > 0:
        out += [(col('Q', ap), -g / pop) for ap in SIBLING_ACTS[a]]
    return out


def qbar_task_terms(t, c, g):
    out = [(col('Qt', cp), g) for cp in CALLS_INTO_TASK[t]]
    pop = TASK_POP[CALL_CALLER_TASK[c]]
    if pop != INF and pop > 0:
        out += [(col('Qt', cp), -g / pop) for cp in SIBLING_CALLS[c]]
    return out


def qd_arg_host(h):
    # Where the rate of host H is read: 1 + delta*sum_a Q(a), the station's own
    # population at the arrival instant. Linear in Q; returned as (const, row).
    r = np.zeros(NVAR)
    for a in ACTS_OF_HOST[h]:
        r[col('Q', a)] += HOST_DELTA[h]
    return 1.0, r


def qd_arg_task(t):
    r = np.zeros(NVAR)
    for c in CALLS_INTO_TASK[t]:
        r[col('Qt', c)] += TASK_DELTA[t]
    return 1.0, r


def drate(f, n):
    # Derivative of a rate function by COMPLEX STEP: f(n + ih) has imaginary
    # part h*f'(n) + O(h^3), so dividing by h is exact to machine precision for
    # any analytic f, with no subtraction and so no cancellation to trade
    # against the step size. This is what keeps the Newton step below quadratic
    # across an arbitrary load-dependent rate.
    return float(np.imag(f(complex(float(n), 1e-20)))) / 1e-20


def mgrad(f, u):
    # Gradient of a MULTIVARIATE smooth term by complex step, one coordinate at
    # a time: exact to machine precision for any analytic f, as in `drate`, at
    # one evaluation per argument. Only the AND-fork closure needs this.
    g = np.zeros(u.size)
    for j in range(u.size):
        z = u.astype(complex)
        z[j] = z[j] + 1e-20j
        g[j] = float(np.imag(f(z))) / 1e-20
    return g


def mhess(f, u, h=1e-5):
    # Second derivatives of a multivariate smooth term, for the exact Hessian
    # only: a central difference of the complex-step gradient, symmetrized.
    H = np.zeros((u.size, u.size))
    for j in range(u.size):
        up = u.copy(); up[j] += h
        um = u.copy(); um[j] -= h
        H[:, j] = (mgrad(f, up) - mgrad(f, um)) / (2.0 * h)
    return 0.5 * (H + H.T)


def d2rate(f, n, h=1e-5):
    # Second derivative, for the exact Hessian only. Central difference of the
    # complex-step first derivative: O(h^2) on a quantity that is already exact,
    # which is far more accurate than differencing f twice and is all the trust
    # region needs.
    return (drate(f, n + h) - drate(f, n - h)) / (2.0 * h)


# ----------------------------------------------------------------------
# linear constraints
# ----------------------------------------------------------------------

def build_linear():
    # Returns (Aeq, beq, Aub, bub): the relations whose every term is a variable
    # times a model parameter, so they hold exactly at the optimum.
    Aeq, beq, Aub, bub = [], [], [], []

    def row():
        return np.zeros(NVAR)

    # call flow: a call inherits the rate of its dispatching activity
    for c in range(NC):
        r = row()
        r[col('Xc', c)] = 1.0
        r[col('Xa', CALL_ACT[c])] -= CALL_MEAN[c]
        Aeq.append(r); beq.append(0.0)

    # entry flow: the requests an entry serves are the calls reaching it. A
    # reference task's entry is driven by its own population instead.
    for e in range(NE):
        if TASK_ISREF[ENTRY_TASK[e]]:
            continue
        r = row()
        r[col('Xe', e)] = 1.0
        for c in CALLS_INTO_ENTRY[e]:
            r[col('Xc', c)] -= 1.0
        Aeq.append(r); beq.append(0.0)

    # activity flow: v(a) executions per invocation of the entry
    for a in range(NA):
        r = row()
        r[col('Xa', a)] = 1.0
        r[col('Xe', ACT_ENTRY[a])] -= ACT_VISITS[a]
        Aeq.append(r); beq.append(0.0)

    # an activity on a delay host never queues: its residence is the bare demand
    for h in range(NH):
        if not HOST_DELAY[h]:
            continue                      # queueing: a residual, see build_residuals
        for a in ACTS_OF_HOST[h]:
            r = row()
            r[col('W', a)] = 1.0
            Aeq.append(r); beq.append(ACT_DEMAND[a])

    # entry composition: a thread is held for the host time of every execution
    # plus the response time of every call those executions dispatch. An entry
    # whose graph FORKS is held for the longest of its concurrent branches
    # instead, which is not linear and moves to build_residuals.
    for e in range(NE):
        if BLOCKS_OF_ENTRY[e]:
            continue                      # concurrent: a residual, see build_residuals
        r = row()
        r[col('S', e)] = 1.0
        for a in ACTS_OF_ENTRY[e]:
            r[col('W', a)] -= ACT_VISITS[a]
            for c in CALLS_OF_ACT[a]:
                r[col('R', c)] -= ACT_VISITS[a] * CALL_MEAN[c]
        Aeq.append(r); beq.append(0.0)

    # thread-seconds: the SAME composition, always additive, because a forking
    # entry occupies its pool for every branch at once even while it is held for
    # only the longest of them. This is the quantity B(c) is built from.
    for e in (range(NE) if NFORK else ()):
        r = row()
        r[col('SA', e)] = 1.0
        for a in ACTS_OF_ENTRY[e]:
            r[col('W', a)] -= ACT_VISITS[a]
            for c in CALLS_OF_ACT[a]:
                r[col('R', c)] -= ACT_VISITS[a] * CALL_MEAN[c]
        Aeq.append(r); beq.append(0.0)

    # a call into a delay task never queues: the congestion law is linear there
    for c in range(NC):
        if TASK_DELAY[ENTRY_TASK[CALL_ENTRY[c]]]:
            r = row()
            r[col('R', c)] = 1.0
            r[col('S', CALL_ENTRY[c])] -= 1.0
            Aeq.append(r); beq.append(0.0)

    # reference closure: the population of a reference task is split between
    # thinking and holding one of its own threads
    for ri, t in enumerate(REF_TASKS):
        r = row()
        r[col('Xe', REF_ENTRY[ri])] = TASK_THINK[t]
        r[col('BR', ri)] = 1.0
        Aeq.append(r); beq.append(TASK_MULT[t])

    # capacity: a processor cannot be busier than it has servers. Only sound
    # while its rate is capped by its server count, which a declared multiplier
    # need not respect, so a station carrying one is left to the bound below.
    for h in range(NH):
        if HOST_DELAY[h] or not HOST_RATE_CAPPED[h]:
            continue
        r = row()
        for a in ACTS_OF_HOST[h]:
            r[col('Xa', a)] += ACT_DEMAND[a]
        Aub.append(r); bub.append(HOST_MULT[h])

    # capacity: a task cannot hold more threads than its pool has
    for t in range(NT):
        if TASK_ISREF[t] or TASK_DELAY[t] or not CALLS_INTO_TASK[t]:
            continue
        r = row()
        for c in CALLS_INTO_TASK[t]:
            r[col('BC', c)] = 1.0
        Aub.append(r); bub.append(TASK_MULT[t])

    # population: a station holds no more jobs than can reach it. Independent of
    # how fast it serves, so this one holds whatever rate was declared.
    for h in range(NH):
        if HOST_DELAY[h] or HOST_NTOT[h] == INF or not ACTS_OF_HOST[h]:
            continue
        r = row()
        for a in ACTS_OF_HOST[h]:
            r[col('Q', a)] = 1.0
        Aub.append(r); bub.append(HOST_NTOT[h])

    for t in range(NT):
        if TASK_DELAY[t] or TASK_NTOT[t] == INF or not CALLS_INTO_TASK[t]:
            continue
        r = row()
        for c in CALLS_INTO_TASK[t]:
            r[col('Qt', c)] = 1.0
        Aub.append(r); bub.append(TASK_NTOT[t])

    Aeq = np.array(Aeq) if Aeq else np.zeros((0, NVAR))
    Aub = np.array(Aub) if Aub else np.zeros((0, NVAR))
    return Aeq, np.array(beq), Aub, np.array(bub)


# ----------------------------------------------------------------------
# nonlinear residuals
# ----------------------------------------------------------------------

def build_residuals():
    # Each residual is  const + A@x + sum coef*x[i]*x[j] + s*f(c + a@x). i == j
    # is allowed in the bilinear part and means a square; the derivative code
    # below handles it without a special case. The last term is the general
    # smooth one, and only the queue-dependent rate laws carry it: exactly one
    # each, a rate function read at a linear form in the queue lengths.
    lin, const, bil, nlt, mnl = [], [], [], [], []

    def row():
        return np.zeros(NVAR)

    k = 0
    # Little at a host: the jobs of an activity are its rate times its residence
    for a in range(NA):
        r = row(); r[col('Q', a)] = 1.0
        lin.append(r); const.append(0.0)
        bil.append((k, col('Xa', a), col('W', a), -1.0)); k += 1

    # Little at a task: the requests of a call are its rate times its response
    for c in range(NC):
        r = row(); r[col('Qt', c)] = 1.0
        lin.append(r); const.append(0.0)
        bil.append((k, col('Xc', c), col('R', c), -1.0)); k += 1

    # busy threads: a call holds a thread of the callee for its service time,
    # and every thread of it a fork inside that entry runs concurrently, so the
    # quantity here is thread-seconds and not elapsed time
    for c in range(NC):
        r = row(); r[col('BC', c)] = 1.0
        lin.append(r); const.append(0.0)
        bil.append((k, col('Xc', c), sent(CALL_ENTRY[c]), -1.0)); k += 1

    # busy threads of a reference task, holding its own entry
    for ri, t in enumerate(REF_TASKS):
        e = REF_ENTRY[ri]
        r = row(); r[col('BR', ri)] = 1.0
        lin.append(r); const.append(0.0)
        bil.append((k, col('Xe', e), col('S', e), -1.0)); k += 1

    # entry composition where the graph FORKS. The branches run concurrently, so
    # the entry is held for the longest of them; E[max] is not linear and is not
    # a function of the branch means alone, so it is taken under the closure
    # fj_max_exp states and sits here rather than in the linear block. One
    # equality is given up and one residual gained, so the counting holds.
    for e in range(NE):
        if not BLOCKS_OF_ENTRY[e]:
            continue
        r = row(); r[col('S', e)] = 1.0
        for a in SEQ_ACTS_OF_ENTRY[e]:
            r[col('W', a)] -= ACT_VISITS[a]
            for c in CALLS_OF_ACT[a]:
                r[col('R', c)] -= ACT_VISITS[a] * CALL_MEAN[c]
        lin.append(r); const.append(0.0)
        for b in BLOCKS_OF_ENTRY[e]:
            A = np.zeros((len(FORK_BRANCHES[b]), NVAR))
            for i, br in enumerate(FORK_BRANCHES[b]):
                for a in br:
                    A[i, col('W', a)] += ACT_VISITS[a]
                    for c in CALLS_OF_ACT[a]:
                        A[i, col('R', c)] += ACT_VISITS[a] * CALL_MEAN[c]
            mnl.append((k, fj_max_exp, -1.0, np.zeros(A.shape[0]), A))
        k += 1

    # host congestion, QD-AMVA. W(a)*G(h) = D(a)*(1 + Qbar(h,a)): the work an
    # execution brings, inflated by what it finds queued, divided by the rate the
    # host runs at then. Bilinear on the left, linear on the right, because D(a)
    # is a model parameter.
    for h in range(NH):
        if HOST_DELAY[h]:
            continue                      # already a linear constraint
        g = col('GH', GH_OF[h])
        for a in ACTS_OF_HOST[h]:
            r = row()
            for j, w in qbar_host_terms(h, a, -ACT_DEMAND[a]):
                r[j] += w
            lin.append(r); const.append(-ACT_DEMAND[a])
            bil.append((k, col('W', a), g, 1.0))
            k += 1

    # task congestion, the same law one layer up. Bilinear on both sides now,
    # because the work a call brings is the callee's entry service time, which
    # is a variable rather than a parameter.
    for c in range(NC):
        e = CALL_ENTRY[c]
        t = ENTRY_TASK[e]
        if TASK_DELAY[t]:
            continue                      # already a linear constraint
        r = row(); r[col('S', e)] = -1.0
        lin.append(r); const.append(0.0)
        bil.append((k, col('R', c), col('GT', GT_OF[t]), 1.0))
        for j, w in qbar_task_terms(t, c, -1.0):
            bil.append((k, col('S', e), j, w))
        k += 1

    # the rate laws themselves: G(i) = r_i(1 + delta_i*sum_j Q_j). The only
    # place the program leaves the quadratic, and the only place the shape of a
    # load-dependent station enters.
    for h in QD_HOSTS:
        r = row(); r[col('GH', GH_OF[h])] = 1.0
        lin.append(r); const.append(0.0)
        c0, arow = qd_arg_host(h)
        nlt.append((k, HOST_RATE[h], -1.0, c0, arow))
        k += 1

    for t in QD_TASKS:
        r = row(); r[col('GT', GT_OF[t])] = 1.0
        lin.append(r); const.append(0.0)
        c0, arow = qd_arg_task(t)
        nlt.append((k, TASK_RATE[t], -1.0, c0, arow))
        k += 1

    L = np.array(lin) if lin else np.zeros((0, NVAR))
    return L, np.array(const), bil, nlt, mnl


AEQ, BEQ, AUB, BUB = build_linear()
RES_LIN, RES_CONST, _BIL, _NLT, _MNL = build_residuals()
NRES = RES_LIN.shape[0]
BIL_K = np.array([t[0] for t in _BIL], dtype=int)
BIL_I = np.array([t[1] for t in _BIL], dtype=int)
BIL_J = np.array([t[2] for t in _BIL], dtype=int)
BIL_C = np.array([t[3] for t in _BIL], dtype=float)
NL_K = [t[0] for t in _NLT]
NL_F = [t[1] for t in _NLT]
NL_S = np.array([t[2] for t in _NLT], dtype=float)
NL_C = np.array([t[3] for t in _NLT], dtype=float)
NL_A = np.array([t[4] for t in _NLT]) if _NLT else np.zeros((0, NVAR))
# The multivariate smooth terms, one per fork block: f is read at a VECTOR of
# linear forms, one per branch, so each carries its own (branches x NVAR) map.
MNL_K = [t[0] for t in _MNL]
MNL_F = [t[1] for t in _MNL]
MNL_S = np.array([t[2] for t in _MNL], dtype=float)
MNL_C = [np.asarray(t[3], dtype=float) for t in _MNL]
MNL_A = [np.asarray(t[4], dtype=float) for t in _MNL]


def nl_args(x):
    # Where each rate function is read at the point X.
    return NL_A @ x + NL_C


def residuals(x):
    r = RES_CONST + RES_LIN @ x
    if BIL_K.size:
        np.add.at(r, BIL_K, BIL_C * x[BIL_I] * x[BIL_J])
    u = nl_args(x)
    for i, k in enumerate(NL_K):
        r[k] += NL_S[i] * float(np.real(NL_F[i](float(u[i]))))
    for i, k in enumerate(MNL_K):
        r[k] += MNL_S[i] * float(np.real(MNL_F[i](MNL_C[i] + MNL_A[i] @ x)))
    return r


def jacobian(x):
    J = RES_LIN.copy()
    if BIL_K.size:
        np.add.at(J, (BIL_K, BIL_I), BIL_C * x[BIL_J])
        np.add.at(J, (BIL_K, BIL_J), BIL_C * x[BIL_I])
    u = nl_args(x)
    for i, k in enumerate(NL_K):
        J[k] += NL_S[i] * drate(NL_F[i], float(u[i])) * NL_A[i]
    for i, k in enumerate(MNL_K):
        J[k] += MNL_S[i] * (mgrad(MNL_F[i], MNL_C[i] + MNL_A[i] @ x) @ MNL_A[i])
    return J


def objective(x):
    r = residuals(x)
    return float(r @ r)


def gradient(x):
    return 2.0 * (jacobian(x).T @ residuals(x))


def hessian(x):
    # 2*J'J from the squares, plus the curvature each bilinear term and each
    # rate function contributes through its own residual.
    J = jacobian(x)
    H = J.T @ J
    r = residuals(x) if (BIL_K.size or NL_K or MNL_K) else None
    if BIL_K.size:
        w = BIL_C * r[BIL_K]
        np.add.at(H, (BIL_I, BIL_J), w)
        np.add.at(H, (BIL_J, BIL_I), w)
    if NL_K:
        u = nl_args(x)
        for i, k in enumerate(NL_K):
            H += (r[k] * NL_S[i] * d2rate(NL_F[i], float(u[i]))) * np.outer(NL_A[i], NL_A[i])
    for i, k in enumerate(MNL_K):
        A = MNL_A[i]
        Hf = mhess(MNL_F[i], MNL_C[i] + A @ x)
        H += (r[k] * MNL_S[i]) * (A.T @ Hf @ A)
    return 2.0 * H


# ----------------------------------------------------------------------
# starting point
# ----------------------------------------------------------------------

def initial_point(warm=0):
    # A point that satisfies every LINEAR relation exactly, built in one pass
    # over the call graph with the queues at zero. WARM fixed-point sweeps of
    # the QD-AMVA laws can be run first for a congested model; the linear
    # relations are then re-derived from the queues so the start stays feasible.
    Q = np.zeros(NA)
    Qt = np.zeros(NC)
    for _ in range(max(0, int(warm)) + 1):
        W, S, SA, R, Xe, Xa, Xc = _sweep(Q, Qt)
        Q = Xa * W
        Qt = Xc * R
    W, S, SA, R, Xe, Xa, Xc = _sweep(Q, Qt)

    x = np.zeros(NVAR)
    x[OFFSET['Xe']:OFFSET['Xe'] + NE] = Xe
    x[OFFSET['Xa']:OFFSET['Xa'] + NA] = Xa
    x[OFFSET['Xc']:OFFSET['Xc'] + NC] = Xc
    x[OFFSET['W']:OFFSET['W'] + NA] = W
    x[OFFSET['S']:OFFSET['S'] + NE] = S
    if NFORK:
        x[OFFSET['SA']:OFFSET['SA'] + NE] = SA
    x[OFFSET['R']:OFFSET['R'] + NC] = R
    x[OFFSET['Q']:OFFSET['Q'] + NA] = Q if warm else np.zeros(NA)
    x[OFFSET['Qt']:OFFSET['Qt'] + NC] = Qt if warm else np.zeros(NC)
    for c in range(NC):
        x[col('BC', c)] = Xc[c] * SA[CALL_ENTRY[c]]
    for ri in range(NR):
        e = REF_ENTRY[ri]
        x[col('BR', ri)] = Xe[e] * S[e]
    Qs = x[OFFSET['Q']:OFFSET['Q'] + NA]
    Qts = x[OFFSET['Qt']:OFFSET['Qt'] + NC]
    for h in QD_HOSTS:
        x[col('GH', GH_OF[h])] = HOST_RATE[h](qd_pop_host(Qs, h))
    for t in QD_TASKS:
        x[col('GT', GT_OF[t])] = TASK_RATE[t](qd_pop_task(Qts, t))
    return x


def qd_pop_host(Q, h):
    # The population the rate of host H is read at, as a number.
    return 1.0 + HOST_DELTA[h] * sum(Q[a] for a in ACTS_OF_HOST[h])


def qd_pop_task(Qt, t):
    return 1.0 + TASK_DELTA[t] * sum(Qt[c] for c in CALLS_INTO_TASK[t])


def qbar_host(Q, h, a):
    # Arrival-instant estimate at host H seen by an arrival of activity A.
    pop = TASK_POP[ACT_TASK[a]]
    tot = sum(Q[ap] for ap in ACTS_OF_HOST[h])
    own = sum(Q[ap] for ap in SIBLING_ACTS[a])
    return tot - (own / pop if pop != INF and pop > 0 else 0.0)


def qbar_task(Qt, t, c):
    pop = TASK_POP[CALL_CALLER_TASK[c]]
    tot = sum(Qt[cp] for cp in CALLS_INTO_TASK[t])
    own = sum(Qt[cp] for cp in SIBLING_CALLS[c])
    return tot - (own / pop if pop != INF and pop > 0 else 0.0)


def _sweep(Q, Qt):
    # One pass of the relations at the given queue lengths: host residence, then
    # entry service and call response callees-first, then the flows
    # callers-first off the reference closure. This is the QD-AMVA fixed-point
    # iteration written out, and it is what the program's residuals say.
    W = np.zeros(NA)
    for h in range(NH):
        for a in ACTS_OF_HOST[h]:
            if HOST_DELAY[h]:
                W[a] = ACT_DEMAND[a]
            else:
                W[a] = (ACT_DEMAND[a] * (1.0 + qbar_host(Q, h, a))
                        / HOST_RATE[h](qd_pop_host(Q, h)))

    S = np.zeros(NE)
    SA = np.zeros(NE)
    R = np.zeros(NC)
    for e in ENTRY_ORDER:
        SA[e] = _compose(ACTS_OF_ENTRY[e], W, R)
        if BLOCKS_OF_ENTRY[e]:
            # concurrent branches: held for the longest, not for their sum
            s = _compose(SEQ_ACTS_OF_ENTRY[e], W, R)
            for b in BLOCKS_OF_ENTRY[e]:
                s += fj_max_exp(np.array([_compose(br, W, R)
                                          for br in FORK_BRANCHES[b]]))
            S[e] = s
        else:
            S[e] = SA[e]
        t = ENTRY_TASK[e]
        for c in CALLS_INTO_ENTRY[e]:
            if TASK_DELAY[t]:
                R[c] = S[e]
            else:
                R[c] = S[e] * (1.0 + qbar_task(Qt, t, c)) / TASK_RATE[t](qd_pop_task(Qt, t))

    Xe = np.zeros(NE)
    Xa = np.zeros(NA)
    Xc = np.zeros(NC)
    for ri, t in enumerate(REF_TASKS):
        e = REF_ENTRY[ri]
        cyc = TASK_THINK[t] + S[e]
        Xe[e] = TASK_MULT[t] / cyc if cyc > 0 else 0.0
    for e in reversed(ENTRY_ORDER):
        if not TASK_ISREF[ENTRY_TASK[e]]:
            Xe[e] = sum(Xc[c] for c in CALLS_INTO_ENTRY[e])
        for a in ACTS_OF_ENTRY[e]:
            Xa[a] = ACT_VISITS[a] * Xe[e]
            for c in CALLS_OF_ACT[a]:
                Xc[c] = CALL_MEAN[c] * Xa[a]
    return W, S, SA, R, Xe, Xa, Xc


def _compose(acts, W, R):
    # The additive composition over a set of activities: host time per
    # execution plus the response of every call it dispatches, at v(a)
    # executions. This is one branch of a fork, or a whole sequential entry.
    s = 0.0
    for a in acts:
        s += ACT_VISITS[a] * W[a]
        for c in CALLS_OF_ACT[a]:
            s += ACT_VISITS[a] * CALL_MEAN[c] * R[c]
    return s


# ----------------------------------------------------------------------
# solve
# ----------------------------------------------------------------------

def polish(x, iters=50):
    # The minimum wanted here is a ROOT, not merely a stationary point: the
    # residuals and the linear equalities together are a SQUARE system, so a few
    # Newton steps finish to machine precision what the trust region leaves at
    # its own tolerance. The step is taken in least-squares form so a rank
    # deficiency degrades rather than raises. Any step that leaves the feasible
    # set ends the polish and the best feasible point so far is kept.
    best, bobj = x.copy(), objective(x)
    scale = 1.0 + float(np.max(np.abs(x))) if x.size else 1.0
    for _ in range(iters):
        F = np.concatenate([residuals(x), AEQ @ x - BEQ])
        JF = np.vstack([jacobian(x), AEQ])
        dx = np.linalg.lstsq(JF, -F, rcond=None)[0]
        x = x + dx
        if np.min(x) < -1e-6 * scale:
            break
        if AUB.shape[0] and np.max(AUB @ x - BUB) > 1e-6 * scale:
            break
        obj = objective(x)
        if obj < bobj:
            best, bobj = x.copy(), obj
        if np.linalg.norm(dx) <= 1e-14 * (1.0 + np.linalg.norm(x)):
            break
    return np.maximum(best, 0.0)


class Result(object):
    # What the Newton fast path returns, shaped like a scipy OptimizeResult.
    def __init__(self, x, message, niter):
        self.x, self.message, self.niter = x, message, niter


def solved(x):
    # A root of the square system, feasible to the bounds and the inequalities.
    scale = 1.0 + (float(np.max(np.abs(x))) if x.size else 0.0)
    if np.min(x) < -1e-9 * scale:
        return False
    if AUB.shape[0] and np.max(AUB @ x - BUB) > 1e-7 * scale:
        return False
    if AEQ.shape[0] and np.max(np.abs(AEQ @ x - BEQ)) > 1e-9 * scale:
        return False
    return bool(np.max(np.abs(residuals(x))) <= 1e-10 * scale) if NRES else True


def solve(warm=0, maxiter=3000, verbose=0):
    # The global minimum is ZERO, so the first thing to try is the root-find:
    # Newton on the square system reaches it directly whenever the starting
    # point is already in the right basin, and is far cheaper than a constrained
    # descent. It can land outside the feasible set, and `solved` is what
    # catches that.
    #
    # WHAT PUTS IT IN THE WRONG BASIN IS A COLD START ON A CONGESTED MODEL: the
    # starting point is the empty network, whose residence times are the bare
    # demands, and a station that actually runs near saturation is nowhere near
    # them. A sweep of the QD-AMVA fixed point costs one pass over the call
    # graph, so ESCALATING THE WARM START IS TRIED FIRST and trust-constr stays
    # the last resort: on a two-layer model at 97% processor utilization the
    # ladder below finds the root in milliseconds where the descent spends
    # seconds and still stops on its evaluation limit.
    ladder = [max(0, int(warm))] + [w for w in (5, 20, 100) if w > warm]
    best_x0, best_obj = None, np.inf
    for w in ladder:
        x0 = initial_point(w)
        x = polish(x0)
        if solved(x):
            return Result(x, 'Newton polish reached a feasible root from %d '
                             'warm sweep%s.' % (w, '' if w == 1 else 's'), 0)
        obj = objective(x0)
        if obj < best_obj:
            best_x0, best_obj = x0, obj

    cons = [LinearConstraint(AEQ, BEQ, BEQ)]
    if AUB.shape[0]:
        cons.append(LinearConstraint(AUB, -np.inf, BUB))
    res = minimize(objective, best_x0, jac=gradient, hess=hessian,
                   method='trust-constr', bounds=Bounds(np.zeros(NVAR), np.inf),
                   constraints=cons,
                   options={'maxiter': maxiter, 'verbose': verbose,
                            'gtol': 1e-12, 'xtol': 1e-14})
    res.x = polish(res.x)
    return res


def report(x):
    out = []
    out.append('LQN %s: nonlinear-program approximation' % MODEL)
    out.append('')
    Xe = x[OFFSET['Xe']:OFFSET['Xe'] + NE]
    Xa = x[OFFSET['Xa']:OFFSET['Xa'] + NA]
    Xc = x[OFFSET['Xc']:OFFSET['Xc'] + NC]
    W = x[OFFSET['W']:OFFSET['W'] + NA]
    S = x[OFFSET['S']:OFFSET['S'] + NE]
    R = x[OFFSET['R']:OFFSET['R'] + NC]
    Q = x[OFFSET['Q']:OFFSET['Q'] + NA]
    Qt = x[OFFSET['Qt']:OFFSET['Qt'] + NC]
    BC = x[OFFSET['BC']:OFFSET['BC'] + NC]
    BR = x[OFFSET['BR']:OFFSET['BR'] + NR]

    GH = x[OFFSET['GH']:OFFSET['GH'] + NGH]
    GT = x[OFFSET['GT']:OFFSET['GT'] + NGT]

    # `rate` is the queue-dependent rate the station is serving at, r_i read at
    # the arrival instant. Below its server count the pool has spare capacity
    # and charges no queueing; at it the station is sharing every server out.
    # `servers` is what the model declares, `usable` what the offered work can
    # reach and what `util` therefore divides by; they differ wherever a pool is
    # wider than its callers.
    out.append('%-16s %12s %12s %12s %12s %12s'
               % ('processor', 'util', 'queue', 'servers', 'usable', 'rate'))
    for h in range(NH):
        busy = sum(ACT_DEMAND[a] * Xa[a] for a in ACTS_OF_HOST[h])
        u = busy if HOST_DELAY[h] else busy / HOST_MAXMULT[h]
        q = sum(Q[a] for a in ACTS_OF_HOST[h])
        out.append('%-16s %12.6g %12.6g %12s %12s %12s'
                   % (HOSTS[h], u, q,
                      'inf' if HOST_DELAY[h] else '%g' % HOST_MULT[h],
                      'inf' if HOST_DELAY[h] else '%g' % HOST_MAXMULT[h],
                      '-' if h not in GH_OF else '%.6g' % GH[GH_OF[h]]))

    out.append('')
    out.append('%-16s %12s %12s %12s %12s %12s'
               % ('task', 'tput', 'util', 'busy', 'threads', 'rate'))
    for t in range(NT):
        xt = sum(Xe[e] for e in ENTRIES_OF_TASK[t])
        if TASK_ISREF[t]:
            busy = sum(BR[ri] for ri, tt in enumerate(REF_TASKS) if tt == t)
        else:
            busy = sum(BC[c] for c in CALLS_INTO_TASK[t])
        u = busy if TASK_DELAY[t] else busy / TASK_MULT[t]
        out.append('%-16s %12.6g %12.6g %12.6g %12s %12s'
                   % (TASKS[t], xt, u, busy,
                      'inf' if TASK_DELAY[t] else '%g' % TASK_MULT[t],
                      '-' if t not in GT_OF else '%.6g' % GT[GT_OF[t]]))

    out.append('')
    if NFORK:
        # `service` is the ELAPSED hold a caller waits for; `threads.s` the
        # thread-seconds the entry occupies of its pool, which stay additive
        # across a fork and are what the pool's occupancy is built from.
        SA = x[OFFSET['SA']:OFFSET['SA'] + NE]
        out.append('%-16s %12s %12s %12s' % ('entry', 'tput', 'service', 'threads.s'))
        for e in range(NE):
            out.append('%-16s %12.6g %12.6g %12.6g' % (ENTRIES[e], Xe[e], S[e], SA[e]))
    else:
        out.append('%-16s %12s %12s' % ('entry', 'tput', 'service'))
        for e in range(NE):
            out.append('%-16s %12.6g %12.6g' % (ENTRIES[e], Xe[e], S[e]))

    out.append('')
    out.append('%-16s %12s %12s %12s' % ('activity', 'tput', 'residence', 'queue'))
    for a in range(NA):
        out.append('%-16s %12.6g %12.6g %12.6g' % (ACTS[a], Xa[a], W[a], Q[a]))

    if NFORK:
        out.append('')
        out.append('%-24s %8s %12s %12s' % ('fork block', 'branches', 'sum',
                                            'concurrent'))
        for b in range(NFORK):
            m = np.array([_compose(br, W, R) for br in FORK_BRANCHES[b]])
            out.append('%-24s %8d %12.6g %12.6g'
                       % (FORK_LABEL[b], m.size, float(m.sum()),
                          float(fj_max_exp(m))))
        out.append('the concurrent time is E[max] taken over EXPONENTIAL branches; '
                   'see the header')

    if NC:
        out.append('')
        out.append('%-16s %12s %12s %12s' % ('call', 'tput', 'response', 'queue'))
        for c in range(NC):
            out.append('%-16s %12.6g %12.6g %12.6g' % (CALLS[c], Xc[c], R[c], Qt[c]))
    return '\n'.join(out)


def main(argv):
    warm, verbose = 0, 0
    for i, a in enumerate(argv):
        if a == '--warm' and i + 1 < len(argv):
            warm = int(argv[i + 1])
        elif a == '--verbose':
            verbose = 2
    res = solve(warm=warm, verbose=verbose)
    x = np.maximum(res.x, 0.0)
    print(report(x))
    print('')
    obj = objective(x)
    eqv = float(np.max(np.abs(AEQ @ x - BEQ))) if AEQ.shape[0] else 0.0
    ubv = float(np.max(AUB @ x - BUB)) if AUB.shape[0] else -np.inf
    print('objective %.6e over %d residuals, %d variables, %d equalities, %d inequalities'
          % (obj, NRES, NVAR, AEQ.shape[0], AUB.shape[0]))
    print('max equality violation %.3e, max inequality margin %.3e (negative = slack)'
          % (eqv, ubv))
    print(res.message if not res.niter else '%s (%d iterations)' % (res.message, res.niter))
    if obj > 1e-8 * max(1.0, float(x @ x)):
        print('WARNING: the objective did not reach zero, so the laws above are not')
        print('         all satisfied. solve() already escalates the warm start up to')
        print('         100 sweeps on its own; retry with a higher --warm, or raise')
        print('         --verbose to inspect the descent it fell back to.')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
'''


def _emit(p: Dict) -> str:
    lines = _header(p) + _params(p)
    return '\n'.join(lines) + _BODY
