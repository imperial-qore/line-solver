"""Shared helpers for cache-aware average tables (native Python).

Provides:
  * retrieval_hidden_classes(sn): indices of auxiliary retrieval classes that
    Cache.set_retrieval_system creates internally and that must be hidden from
    the node table (plumbing, not user-facing classes).
  * build_cache_avg_table(solver): the getAvgCacheTable DataFrame, with a total
    row (List=0) per cache+read class plus per-list rows where available.
"""

import numpy as np
import pandas as pd

from ..api.sn.network_struct import NodeType


def retrieval_hidden_classes(sn):
    """Return a set of 0-based class indices that are auxiliary retrieval classes."""
    hidden = set()
    if sn is None or sn.nodeparam is None:
        return hidden
    for ind in range(sn.nnodes):
        if sn.nodetype is not None and ind < len(sn.nodetype) and sn.nodetype[ind] == NodeType.CACHE:
            np_ = sn.nodeparam.get(ind)
            if np_ is None:
                continue
            rci = getattr(np_, 'retrieval_class_indices', None)
            if rci:
                for idx in rci:
                    if idx is not None and 0 <= int(idx) < sn.nclasses:
                        hidden.add(int(idx))
    return hidden


def _nan_at(vec, r):
    if vec is None:
        return np.nan
    v = np.atleast_1d(np.asarray(vec)).flatten()
    if r >= len(v):
        return np.nan
    return float(v[r])


def build_cache_avg_table(solver):
    """Build the detailed per-class cache performance table for a solver.

    Columns: Node, JobClass, List, ListCap, Items, HitProb, DelayedHitProb,
    MissProb, HitRate, DelayedHitRate, MissRate, ArvR, ResidT. One total row
    (List=0) per cache+read class plus, where per-list hit probabilities are
    available and the cache has more than one list, one row per cache list.
    """
    sn = solver._sn
    cols = ['Node', 'JobClass', 'List', 'ListCap', 'Items', 'HitProb',
            'DelayedHitProb', 'MissProb', 'HitRate', 'DelayedHitRate',
            'MissRate', 'ArvR', 'ResidT']
    rows = []
    any_cache = sn is not None and any(
        sn.nodetype is not None and i < len(sn.nodetype) and sn.nodetype[i] == NodeType.CACHE
        for i in range(sn.nnodes))
    if any_cache:
        # Node-level throughput; the cache read-class arrival equals that class's
        # source throughput (every read request enters the cache).
        _, _, _, _, _, TNn = solver.getAvgNode()
        src_node = None
        for i in range(sn.nnodes):
            if sn.nodetype is not None and i < len(sn.nodetype) and sn.nodetype[i] == NodeType.SOURCE:
                src_node = i
                break
        nodenames = list(sn.nodenames) if hasattr(sn, 'nodenames') and sn.nodenames else []
        model_nodes = getattr(getattr(solver, 'model', None), '_nodes', None)
        for ind in range(sn.nnodes):
            if sn.nodetype[ind] != NodeType.CACHE:
                continue
            np_ = sn.nodeparam.get(ind)
            if np_ is None:
                continue
            hitclass = np.atleast_1d(np.asarray(getattr(np_, 'hitclass', []))).flatten()
            missclass = np.atleast_1d(np.asarray(getattr(np_, 'missclass', []))).flatten()
            itemcap = np.atleast_1d(np.asarray(getattr(np_, 'itemcap', []))).flatten()
            h = len(itemcap)
            totcap = float(np.sum(itemcap)) if h > 0 else 0.0
            nitems = int(getattr(np_, 'nitems', 0))
            node = model_nodes[ind] if model_nodes is not None and ind < len(model_nodes) else None
            hitp = node.get_hit_ratio() if node is not None else getattr(np_, 'actualhitprob', None)
            missp = node.get_miss_ratio() if node is not None else getattr(np_, 'actualmissprob', None)
            dhitp = node.get_delayed_hit_ratio() if (node is not None and hasattr(node, 'get_delayed_hit_ratio')) else getattr(np_, 'actualdelayedhitprob', None)
            lat = node.get_residt() if node is not None else getattr(np_, 'actualresidt', None)
            hplist = node.get_hit_ratio_by_list() if (node is not None and hasattr(node, 'get_hit_ratio_by_list')) else getattr(np_, 'actualhitproblist', None)
            node_name = nodenames[ind] if ind < len(nodenames) else f'Node{ind}'
            for r in range(sn.nclasses):
                if r >= len(hitclass) or hitclass[r] <= 0:
                    continue
                ph = _nan_at(hitp, r); pm = _nan_at(missp, r); pd_ = _nan_at(dhitp, r)
                if np.isnan(ph) and np.isnan(pm) and np.isnan(pd_):
                    continue
                if np.isnan(ph): ph = 0.0
                if np.isnan(pm): pm = 0.0
                if np.isnan(pd_): pd_ = 0.0
                # Read-class arrival = that class's source throughput (every read
                # request enters the cache); robust across solvers, including
                # simulators where delayed hits are not folded into hit/miss tput.
                arvr = 0.0
                if TNn is not None and src_node is not None and src_node < TNn.shape[0] and r < TNn.shape[1]:
                    arvr = TNn[src_node, r]
                cnames = getattr(solver, 'class_names', None) or getattr(solver, '_class_names', None) or []
                class_name = cnames[r] if r < len(cnames) else f'Class{r}'
                latr = _nan_at(lat, r)
                # ArvR is the retrieval-system throughput lambda*(missprob+delayedprob)
                # = MissRate + DelayedHitRate, i.e. the rate of requests that enter the
                # retrieval system (a miss or a delayed hit). This is the arrival rate that
                # is Little-consistent with ResidT = Z (the delayed-hit expected latency,
                # eq:latency tot / retrieval_fpi_latency): ArvR*ResidT = sum_i(phi_i+d_i),
                # the mean number of requests in the retrieval system (fetch job included).
                arvr_retr = arvr * (pm + pd_)
                rows.append(dict(zip(cols, [node_name, class_name, 0, totcap, nitems,
                                            ph, pd_, pm, arvr * ph, arvr * pd_, arvr * pm, arvr_retr, latr])))
                # per-list rows
                if h > 1 and hplist is not None:
                    hpl = np.atleast_2d(np.asarray(hplist))
                    if r < hpl.shape[0] and not np.all(np.isnan(hpl[r, :])):
                        for l in range(h):
                            phl = hpl[r, l] if l < hpl.shape[1] else np.nan
                            if np.isnan(phl): phl = 0.0
                            capl = float(itemcap[l]) if l < len(itemcap) else np.nan
                            rows.append(dict(zip(cols, [node_name, class_name, l + 1, capl, nitems,
                                                        phl, np.nan, np.nan, arvr * phl, np.nan, np.nan,
                                                        arvr, np.nan])))
    df = pd.DataFrame(rows, columns=cols)
    if not getattr(solver, '_table_silent', False):
        print(df.to_string(index=False))
    return df


def build_item_avg_table(solver):
    """Build the item-level cache occupancy table for a solver.

    Columns: Node, Item, List, ListCap, Prob. One row per cache node, item and
    cache list, where Prob is the steady-state probability the item resides in
    that list. Populated only where the solver computes a per-item distribution.
    """
    sn = solver._sn
    cols = ['Node', 'Item', 'List', 'ListCap', 'Prob']
    rows = []
    any_cache = sn is not None and any(
        sn.nodetype is not None and i < len(sn.nodetype) and sn.nodetype[i] == NodeType.CACHE
        for i in range(sn.nnodes))
    if any_cache:
        solver.getAvgNode()  # ensure solved so cache item probabilities are filled
        nodenames = list(sn.nodenames) if hasattr(sn, 'nodenames') and sn.nodenames else []
        model_nodes = getattr(getattr(solver, 'model', None), '_nodes', None)
        for ind in range(sn.nnodes):
            if sn.nodetype[ind] != NodeType.CACHE:
                continue
            np_ = sn.nodeparam.get(ind)
            if np_ is None:
                continue
            itemcap = np.atleast_1d(np.asarray(getattr(np_, 'itemcap', []))).flatten()
            h = len(itemcap)
            node = model_nodes[ind] if model_nodes is not None and ind < len(model_nodes) else None
            itemprob = node.get_item_prob() if (node is not None and hasattr(node, 'get_item_prob')) else getattr(np_, 'actualitemprob', None)
            if itemprob is None or h == 0:
                continue
            itemprob = np.atleast_2d(np.asarray(itemprob))
            n = itemprob.shape[0]
            node_name = nodenames[ind] if ind < len(nodenames) else f'Node{ind}'
            for i in range(n):
                for l in range(h):
                    p = itemprob[i, l + 1] if (l + 1) < itemprob.shape[1] else np.nan
                    rows.append(dict(zip(cols, [node_name, i + 1, l + 1, float(itemcap[l]), p])))
    df = pd.DataFrame(rows, columns=cols)
    if not getattr(solver, '_table_silent', False):
        print(df.to_string(index=False))
    return df
