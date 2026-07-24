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
        fanOut = _get_fanout(sn, f)
        if fanOut is not None and np.any(np.asarray(fanOut) != np.round(np.asarray(fanOut))):
            raise RuntimeError('Non-integer tasksPerLink is not supported by the native CTMC/SSA fork-join implementation.')
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
                if r in js and js[r] is not None and _js_val(js[r]) != int(JoinStrategy.STD):
                    raise RuntimeError('Only JoinStrategy.STD is supported by the native CTMC/SSA fork-join implementation.')


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
