"""
Validate a fork-join model against the native CTMC/SSA supported feature set.

Port of matlab/src/api/fj/sn_fj_validate.m.
"""

import numpy as np
from ...lang.base import NodeType, JoinStrategy


def _nt_val(nt):
    return int(nt.value) if hasattr(nt, 'value') else int(nt)


def sn_fj_validate(sn):
    """
    Validate that a fork-join model is within the supported feature set of
    the native CTMC/SSA fork-join implementation (v1): closed classes only,
    non-nested fork-join pairs, standard join strategy (wait for all
    siblings), one task per output link. Unsupported models raise a
    RuntimeError with a specific message; further structural checks (class
    switching on a branch, nesting, sibling traps) are performed during the
    branch discovery in ModelAdapter.fjtag.
    """
    from ...io.model_adapter import ModelAdapter

    K = sn.nclasses
    fork_val = _nt_val(NodeType.FORK)
    join_val = _nt_val(NodeType.JOIN)
    nodetype_vals = np.array([_nt_val(nt) for nt in sn.nodetype])
    forkIndexes = np.where(nodetype_vals == fork_val)[0]
    joinIndexes = np.where(nodetype_vals == join_val)[0]
    Vnodes = ModelAdapter._compute_node_visits(sn)

    njobs = np.asarray(sn.njobs).ravel()
    chains = np.asarray(sn.chains) if sn.chains is not None else np.eye(K, dtype=bool)

    for f in forkIndexes:
        j = np.where(sn.fj[f, :])[0]
        if j.size == 0:
            raise RuntimeError('Fork nodes without a matched Join are not supported by the native CTMC/SSA fork-join implementation.')
        if j.size > 1:
            raise RuntimeError('Multiple Join nodes per Fork are not supported by the native CTMC/SSA fork-join implementation.')
        # TASKS PER LINK. The tag-augmented construction carries an integer
        # weight per branch, so a link whose count is an integer is served
        # exactly, whether or not it differs from its siblings'. What it cannot
        # carry is a count that is not fixed at build time: a DISTRIBUTION has
        # to be drawn per firing.
        fp = sn.nodeparam.get(f) if sn.nodeparam is not None else None
        fp = fp if isinstance(fp, dict) else None
        if fp is not None and fp.get('fanOutDist', None) is not None:
            if any(d is not None for row in fp['fanOutDist'] for d in row):
                raise RuntimeError(
                    'A random tasks-per-link distribution '
                    '(Fork.setTasksPerLinkDistribution) is not supported by the native '
                    'CTMC/SSA fork-join implementation; use SolverJMT or SolverLDES, '
                    'which draw the degree at the fork epoch.')
        if fp is not None and fp.get('fanOutLink', None) is not None:
            taken = fp['fanOutProb'] > 0
            vals = fp['fanOutLink'][taken]
            if np.any(vals != np.round(vals)):
                raise RuntimeError('Non-integer tasksPerLink is not supported by the native CTMC/SSA fork-join implementation.')
        else:
            fanOut = _get_fanout(sn, f)
            if fanOut is not None and np.any(np.asarray(fanOut) != np.round(np.asarray(fanOut))):
                raise RuntimeError('Non-integer tasksPerLink is not supported by the native CTMC/SSA fork-join implementation.')

        # BRANCH PROBABILITIES. A branch that may decline makes the SET of
        # siblings random, so the firing has 2^B outcomes and the Join's
        # required count is a function of which subset fired. The tag
        # construction records no such per-firing state, so this path would emit
        # every branch anyway and answer with the certain-fork number. Refused
        # by name instead: silently returning the answer for a DIFFERENT model
        # is the one outcome worth avoiding.
        if fp is not None and fp.get('fanOutProb', None) is not None:
            pb = fp['fanOutProb']
            taken = pb > 0
            if np.any(pb[taken] < 1.0):
                raise RuntimeError(
                    'A branch activation probability below one '
                    '(Fork.setBranchProbability) is not supported by the native CTMC/SSA '
                    'fork-join implementation: the sibling SET would be random and the tag '
                    'construction fixes it when the state space is built. Use SolverJMT or '
                    'SolverLDES, which draw the activation at the fork epoch, or SolverMVA, '
                    'whose MMT transform sees the expected degree.')
        for r in np.where(Vnodes[f, :] > 0)[0]:
            c = _find_chain(chains, r)
            classes_in_chain = np.where(chains[c, :])[0] if c is not None else np.array([r])
            if not np.isfinite(njobs[r]) or np.any(~np.isfinite(njobs[classes_in_chain])):
                raise RuntimeError('Open classes routed through a Fork are not supported by the native CTMC/SSA fork-join implementation.')

    for j in joinIndexes:
        if not np.any(sn.fj[:, j]):
            raise RuntimeError('Join nodes without a matched Fork are not supported by the native CTMC/SSA fork-join implementation.')
        js = _get_join_strategy(sn, j)
        if js is not None:
            for r in range(K):
                # PARTIAL is served: the per-branch required count is that
                # count lowered, not a different mechanism.
                if r in js and js[r] is not None and _js_val(js[r]) not in (
                        int(JoinStrategy.STD), int(JoinStrategy.PARTIAL)):
                    raise RuntimeError('Only JoinStrategy.STD and JoinStrategy.PARTIAL are supported by the native CTMC/SSA fork-join implementation.')


def _get_fanout(sn, f):
    if sn.nodeparam is not None and f in sn.nodeparam:
        p = sn.nodeparam[f]
        if isinstance(p, dict):
            return p.get('fanOut', None)
        return getattr(p, 'fanOut', None)
    return None


def _get_join_strategy(sn, j):
    if sn.nodeparam is not None and j in sn.nodeparam:
        p = sn.nodeparam[j]
        if isinstance(p, dict):
            return p.get('joinStrategy', None)
        return getattr(p, 'joinStrategy', None)
    return None


def _js_val(js):
    return int(js.value) if hasattr(js, 'value') else int(js)


def _find_chain(chains, r):
    idx = np.where(chains[:, r])[0]
    return int(idx[0]) if idx.size > 0 else None


def sn_fj_supports(sn):
    """Can the exact fork-join construction be asked for this model?

    The fork-join model class :func:`sn_fj_validate` admits, asked as a
    predicate rather than raised. ``SolverCTMC.supportsModelMethod`` and
    ``SolverSSA.supportsModelMethod`` call it so that a caller (model.help,
    findSolver, SolverAUTO) sees the verdict before paying for a run, and BOTH
    analyzers reach the SAME rules through ``ModelAdapter.fjtag``. The sentence
    the validator raises names 'the native CTMC/SSA fork-join implementation',
    which is why this predicate lives beside it rather than inside either solver.

    IT WRAPS THE VALIDATOR RATHER THAN RESTATING IT, and that is the point: the
    rules are eight and they move (pairing, join strategy, tasks-per-link,
    branch probability, open classes through a fork), so a second copy would be
    a second thing to keep in step. There is exactly one body of rules and two
    ways in -- one that raises, for the run, and this one, which answers.

    WHAT IT REFUSES AND WHY THE ANALYZER IS RIGHT TO. The fork-join PAIRING is a
    declaration carried by the Join (``Join(model, name, fork)`` in all four
    codebases), not a derivation from the routing: a nested model such as
    examples/basic/forkJoin/fj_basic_nesting has two forks and two joins whose
    pairing the routing alone does not determine. So a Join built without naming
    its fork leaves ``sn.fj`` empty, and 'Fork nodes without a matched Join' is
    the honest answer to a model that declares none -- not a topology test that
    failed to see one.

    Args:
        sn: NetworkStruct of the model.

    Returns:
        (ok, reason); reason is '' when ok is True.
    """
    nodetype = np.ravel(np.asarray(sn.nodetype, dtype=int))
    if not (np.any(nodetype == _nt_val(NodeType.FORK))
            or np.any(nodetype == _nt_val(NodeType.JOIN))):
        return True, ''
    try:
        sn_fj_validate(sn)
    except Exception as err:      # the validator's own refusal, turned into an answer
        return False, str(err)
    return True, ''
