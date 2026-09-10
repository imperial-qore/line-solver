"""
Phase-type composition of an LQN activity graph, the machinery behind SolverLN
method 'srvn.ph'.

An entry becomes a Workflow whose leaves are its activities and, when requested,
its synchronous calls; the series-parallel reduction of that workflow is then
the exact law of the entry service time. Twin of the MATLAB
lqn_entry_workflow.m, lqn_ph_serial_law.m and lqn_ph_moments.m, and of the JAR
jline.api.lqn.LqnPh.
"""

from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from ...constants import GlobalConstants
from ...lang.workflow import ActivityPrecedence as WfPrecedence
from ...lang.workflow import Workflow

__all__ = ['EntryWorkflow', 'entry_workflow', 'serial_law', 'ph_moments']


class EntryWorkflow:
    """Activity graph of one entry, as a workflow plus its execution counts."""

    def __init__(self):
        self.wf: Optional[Workflow] = None
        self.act_idx_of: Dict[int, int] = {}
        self.call_idx_of: Dict[int, int] = {}
        self.execs: Dict[int, float] = {}
        self.callexecs: Dict[int, float] = {}


def entry_workflow(model, lqn, eidx: int, with_calls: bool = True) -> EntryWorkflow:
    """
    Activity graph of LQN entry EIDX as a Workflow.

    The activity precedences are read from the task object rather than
    reconstructed from lqn.graph, whose loop back-edges carry probabilities and
    not counts.

    Args:
        model: the LayeredNetwork the struct LQN was obtained from
        lqn: the LayeredNetworkStruct
        eidx: absolute index of the entry
        with_calls: True expands every synchronous call of an activity into a
            leaf of its own, placed in series after the activity, so that the
            call response law and the host demand law stay separable across
            iterations. False keeps only the host demands, which is the
            processor-demand law of the entry: the host is released while a call
            is outstanding.

    Returns:
        an EntryWorkflow

    Raises:
        ValueError: when the entry binds no activity, or its precedence graph is
            not series-parallel
    """
    tidx = int(lqn.parent[eidx, 0])
    task = _task_named(model, _name_of(lqn, tidx))
    if task is None:
        raise ValueError("Entry %s has no task in the layered model." % _hash_of(lqn, eidx))
    acts = list(lqn.actsof.get(eidx, []))
    if not acts:
        raise ValueError("Entry %s binds no activity." % _hash_of(lqn, eidx))

    out = EntryWorkflow()
    out.wf = Workflow(_name_of(lqn, eidx) + '.Workflow')
    head_name: Dict[int, str] = {}
    tail_name: Dict[int, str] = {}
    act_of_name: Dict[str, int] = {}

    for aidx in acts:
        nm = _name_of(lqn, aidx)
        a = out.wf.addActivity(nm, _host_demand(lqn, aidx))
        out.act_idx_of[aidx] = a.index
        head_name[aidx] = nm
        tail_name[aidx] = nm
        act_of_name[nm] = aidx

    if with_calls:
        for aidx in acts:
            chain = [head_name[aidx]]
            for cidx in lqn.callsof.get(aidx, []):
                if _call_type(lqn, cidx) != 1:
                    continue  # an asynchronous call blocks the caller for no time
                cnm = call_hashname(lqn, cidx)
                c = out.wf.addActivity(cnm, _immediate())
                out.call_idx_of[cidx] = c.index
                chain.append(cnm)
            if len(chain) > 1:
                for k in range(1, len(chain)):
                    out.wf.addPrecedence(WfPrecedence.Serial(chain[k - 1], chain[k]))
                tail_name[aidx] = chain[-1]

    # Precedences of the task, restricted to the activities of this entry and
    # rewritten so that a predecessor is entered at its head and left at its tail
    for prec in task.precedences:
        for wprec in _translate(prec, act_of_name, head_name, tail_name):
            out.wf.addPrecedence(wprec)

    ok, msg = out.wf.validate()
    if not ok:
        raise ValueError("Entry %s cannot be composed into a phase-type law: %s"
                         % (_hash_of(lqn, eidx), msg))
    tree = out.wf.getSPTree()
    if tree is None:
        raise ValueError("Entry %s has a precedence graph that is not series-parallel, so its "
                         "activity graph has no exact phase-type reduction. Use method='default'."
                         % _hash_of(lqn, eidx))

    for aidx in acts:
        out.execs[aidx] = float(tree['execs'][tree['leaf_of'][out.act_idx_of[aidx]]])
    for cidx, widx in out.call_idx_of.items():
        out.callexecs[cidx] = float(tree['execs'][tree['leaf_of'][widx]])
    return out


def serial_law(wf: Workflow) -> Tuple[np.ndarray, np.ndarray]:
    """
    Composed law of a workflow in which the branches of an AND fork are SERIAL
    rather than concurrent, that is, the total work the branches request rather
    than the elapsed time until the last of them finishes.

    This is the law of the PROCESSOR demand of an LQN entry. Two branches of an
    AND fork are two activity threads of the same task instance: they overlap in
    time, so the entry response time is the maximum of the branches, but they run
    on ONE processor, so the demand they place on it is the sum. Composing the
    host law with Workflow.toPH would charge the processor the maximum and let
    the layer report a utilization below the true one, which no amount of
    iterating recovers.

    Every other node keeps its own composition rule: an OR fork is a mixture, a
    loop is a geometric compound, so the correlation within a branch survives.
    """
    tree = wf.getSPTree()
    if tree is None:
        raise ValueError("Workflow %s is not series-parallel, so it has no exact phase-type "
                         "reduction." % wf.name)
    return _compose_serialized(wf, tree, tree['root'])


def _compose_serialized(wf: Workflow, tree: Dict[str, Any], k: int) -> Tuple[np.ndarray, np.ndarray]:
    kids = tree['kids'][k]
    ntype = tree['type'][k]
    if ntype == 'leaf':
        return wf.getActivities()[tree['act'][k]].getPHRepresentation()
    if ntype in ('serial', 'par'):
        alpha, T = _compose_serialized(wf, tree, kids[0])
        for i in range(1, len(kids)):
            a2, T2 = _compose_serialized(wf, tree, kids[i])
            alpha, T = Workflow._composeSerial(alpha, T, a2, T2)
        return alpha, T
    if ntype == 'or':
        alphas, Ts = [], []
        for k2 in kids:
            a2, T2 = _compose_serialized(wf, tree, k2)
            alphas.append(a2)
            Ts.append(T2)
        return Workflow._composeMixture(alphas, Ts, tree['probs'][k])
    if ntype == 'loop':
        a1, T1 = _compose_serialized(wf, tree, kids[0])
        return Workflow._composeLoopGeometric(a1, T1, tree['count'][k])
    raise ValueError('Unknown series-parallel node type "%s".' % ntype)


def ph_moments(alpha: np.ndarray, T: np.ndarray) -> Tuple[float, float]:
    """
    First two moments of the phase-type law (ALPHA, T) without building a
    Distribution object, which is what the layered fixed point needs at every
    iteration for every composed entry law. A defective ALPHA carries an atom at
    zero and contributes nothing to either moment.

    Returns:
        (mean, squared coefficient of variation)
    """
    a = np.asarray(alpha, dtype=float).reshape(1, -1)
    Tm = np.asarray(T, dtype=float)
    e = np.ones((Tm.shape[0], 1))
    try:
        x1 = -np.linalg.solve(Tm, e)
        m1 = float((a @ x1)[0, 0])
        x2 = -np.linalg.solve(Tm, x1)
        m2 = 2.0 * float((a @ x2)[0, 0])
    except np.linalg.LinAlgError:
        return GlobalConstants.FineTol, 1.0
    if not np.isfinite(m1) or m1 <= GlobalConstants.FineTol:
        return GlobalConstants.FineTol, 1.0
    scv = m2 / (m1 * m1) - 1.0
    if not np.isfinite(scv) or scv <= GlobalConstants.FineTol:
        scv = GlobalConstants.FineTol
    return m1, scv


def call_hashname(lqn, cidx: int) -> str:
    """Name of call CIDX, in the caller=>callee form MATLAB uses."""
    src = int(lqn.callpair[cidx, 0])
    dst = int(lqn.callpair[cidx, 1])
    return "%s=>%s" % (_hash_of(lqn, src), _hash_of(lqn, dst))


# ----------------------------------------------------------------------
# translation of an LQN precedence onto the expanded workflow chains
# ----------------------------------------------------------------------

def _translate(prec, act_of_name: Dict[str, int], head_name: Dict[int, str],
               tail_name: Dict[int, str]) -> List[WfPrecedence]:
    """
    Rewrite one LQN activity precedence as workflow precedences, or return an
    empty list when it names an activity outside this entry.

    The two representations differ: an LQN precedence carries a PrecedenceType
    plus activity lists, a workflow one carries the pre/post type pair of an
    activity graph. A predecessor is left at its TAIL (the last call leaf spliced
    after it) and a successor entered at its HEAD.
    """
    from ...layered import PrecedenceType

    def tail(a):
        idx = act_of_name.get(_act_name(a))
        return None if idx is None else tail_name[idx]

    def head(a):
        idx = act_of_name.get(_act_name(a))
        return None if idx is None else head_name[idx]

    ptype = getattr(prec, 'prec_type', None)
    pre_acts = list(getattr(prec, 'pre_activities', []) or [])
    post_acts = list(getattr(prec, 'post_activities', []) or [])
    acts = list(getattr(prec, 'activities', []) or [])

    if ptype == PrecedenceType.SERIAL:
        names = [(_act_name(a)) for a in acts]
        if any(nm not in act_of_name for nm in names):
            return []
        out = []
        for k in range(1, len(names)):
            out.append(WfPrecedence.Serial(tail(acts[k - 1]), head(acts[k])))
        return out

    if ptype == PrecedenceType.LOOP:
        if not pre_acts or not acts:
            return []
        pre = tail(pre_acts[0])
        body = [head(a) for a in acts]
        if pre is None or any(b is None for b in body):
            return []
        return [WfPrecedence.Loop(pre, body, float(getattr(prec, 'count', 1.0)))]

    if ptype == PrecedenceType.PARALLEL:
        if len(pre_acts) == 1 and len(post_acts) >= 1 and len(post_acts) > 1:
            pre = tail(pre_acts[0])
            posts = [head(a) for a in post_acts]
            if pre is None or any(p is None for p in posts):
                return []
            return [WfPrecedence.AndFork(pre, posts)]
        if len(post_acts) == 1 and len(pre_acts) >= 1:
            pres = [tail(a) for a in pre_acts]
            post = head(post_acts[0])
            if post is None or any(p is None for p in pres):
                return []
            quorum = getattr(prec, 'pre_params', None)
            return [WfPrecedence.AndJoin(pres, post, quorum)]
        return []

    if ptype == PrecedenceType.CHOICE:
        probs = list(getattr(prec, 'probabilities', []) or [])
        if len(pre_acts) == 1 and len(post_acts) > 1:
            pre = tail(pre_acts[0])
            posts = [head(a) for a in post_acts]
            if pre is None or any(p is None for p in posts):
                return []
            if not probs:
                probs = [1.0 / len(posts)] * len(posts)
            return [WfPrecedence.OrFork(pre, posts, np.asarray(probs, dtype=float))]
        if len(post_acts) == 1 and len(pre_acts) >= 1:
            pres = [tail(a) for a in pre_acts]
            post = head(post_acts[0])
            if post is None or any(p is None for p in pres):
                return []
            return [WfPrecedence.OrJoin(pres, post)]
        return []

    return []


# ----------------------------------------------------------------------
# small struct accessors, tolerant of the dict/array duality of the struct
# ----------------------------------------------------------------------

def _act_name(a) -> str:
    return a if isinstance(a, str) else a.name


def _name_of(lqn, idx: int) -> str:
    v = lqn.names
    if isinstance(v, dict):
        return v.get(idx, 'Node_%d' % idx)
    return str(v[idx])


def _hash_of(lqn, idx: int) -> str:
    v = getattr(lqn, 'hashnames', None)
    if v is None:
        return _name_of(lqn, idx)
    if isinstance(v, dict):
        return v.get(idx, _name_of(lqn, idx))
    return str(v[idx])


def _call_type(lqn, cidx: int) -> int:
    ct = getattr(lqn, 'calltype', None)
    if ct is None:
        return 1
    arr = np.asarray(ct).ravel()
    return int(arr[cidx]) if cidx < len(arr) else 1


def _host_demand(lqn, aidx: int):
    proc = getattr(lqn, 'hostdem_proc', None)
    if isinstance(proc, dict) and proc.get(aidx) is not None:
        return proc[aidx]
    hd = lqn.hostdem
    m = hd.get(aidx, 0.0) if isinstance(hd, dict) else float(np.asarray(hd).ravel()[aidx])
    return max(float(m), GlobalConstants.FineTol)


def _immediate():
    from ...distributions import Immediate
    return Immediate()


def _task_named(model, name: str):
    for t in getattr(model, 'tasks', []):
        if str(t.name) == name:
            return t
    return None
