"""
Quasi-reduced ordered Multi-valued Decision Diagram.

Compact symbolic store for a set of discrete states, after
A.S. Miner, G. Ciardo, "Efficient Reachability Set Generation and Storage
Using Decision Diagrams", ICATPN 1999, LNCS 1639, pp.6-25.

A global state is a K-tuple of *local* state values (one per level/submodel),
state[k] in {0,...,domain[k]-1}. The set is stored as a directed acyclic graph
with K variable levels plus a terminal level: level 1 is the top (root), a node
at level k has domain[k] outgoing arcs to level k+1 nodes, and a state belongs
to the set iff its path of arcs reaches the TRUE terminal. Canonicity is
enforced by a per-level unique table (no duplicate nodes) and by collapsing the
all-FALSE node to the FALSE terminal. Storage is O(#nodes), typically
O(K * #local-states), instead of O(|S|) as in an explicit state list.

Two constant terminals encode the boolean value of a completed path:
    TERM_FALSE = 0  (empty subgraph / state not in set)
    TERM_TRUE  = -1 (state in set)
Arcs of a level-k node hold ids of level-(k+1) nodes when k<K, or a terminal
when k==K.

Node ids are 1-based positive integers, matching the MATLAB reference
implementation, so that 0 can serve as TERM_FALSE.

See also: mdd_reachset, mdd_mcd.
"""

# Copyright (c) 2012-2026, Imperial College London
# All rights reserved.

from typing import Dict, List, Sequence, Tuple

import numpy as np

TERM_TRUE = -1
TERM_FALSE = 0


class MDD:
    """Quasi-reduced ordered multi-valued decision diagram over K levels."""

    TERM_TRUE = TERM_TRUE
    TERM_FALSE = TERM_FALSE

    def __init__(self, domain: Sequence[int]):
        """Create an empty set over the given per-level domains.

        domain is a length-K sequence; values at level k are 0..domain[k]-1.
        """
        self.domain: List[int] = [int(d) for d in np.ravel(np.asarray(domain))]
        self.K: int = len(self.domain)
        # node[k] is a list of arc rows; row r-1 holds the arcs of node id r
        self.node: List[List[List[int]]] = [[] for _ in range(self.K)]
        self.uniq: List[Dict[Tuple[int, ...], int]] = [dict() for _ in range(self.K)]
        self.root: int = TERM_FALSE
        self._cnt: List[List[int]] = [[] for _ in range(self.K)]
        self._dirty: bool = True

    # -- construction ----------------------------------------------------

    def insert(self, state: Sequence[int]) -> None:
        """Add a K-tuple of 0-based local values to the set."""
        self.root = self._add_state(0, self.root, state)
        self._dirty = True

    @staticmethod
    def from_states(domain: Sequence[int], states) -> "MDD":
        """Build an MDD from a matrix whose rows are 0-based state tuples."""
        obj = MDD(domain)
        for row in np.atleast_2d(np.asarray(states)):
            obj.insert([int(v) for v in row])
        return obj

    # -- queries ---------------------------------------------------------

    def member(self, state: Sequence[int]) -> bool:
        """Return True iff state is in the set (O(K))."""
        node_id = self.root
        for k in range(self.K):
            if node_id == TERM_FALSE:
                return False
            node_id = self.node[k][node_id - 1][int(state[k])]
        return node_id == TERM_TRUE

    def cardinality(self) -> int:
        """Return |S|, the number of stored states."""
        self._ensure_counts()
        return self._child_count(0, self.root)

    def index(self, state: Sequence[int]) -> int:
        """0-based lexicographic rank of state (level 1 most significant).

        Returns -1 when state is not stored. This is a bijection
        S <-> {0,...,|S|-1}, so a generator can be assembled without an
        explicit state list.
        """
        self._ensure_counts()
        idx = 0
        node_id = self.root
        for k in range(self.K):
            if node_id == TERM_FALSE:
                return -1
            arcs = self.node[k][node_id - 1]
            v = int(state[k])
            for vv in range(v):
                idx += self._child_count(k + 1, arcs[vv])
            node_id = arcs[v]
        if node_id != TERM_TRUE:
            return -1
        return idx

    def enumerate(self) -> np.ndarray:
        """All stored states as rows, in index() order."""
        if self.root == TERM_FALSE:
            return np.zeros((0, self.K), dtype=int)
        return np.asarray(self._enum_below(0, self.root), dtype=int).reshape(-1, self.K)

    # -- export ----------------------------------------------------------

    def to_struct(self) -> "MDDStruct":
        """Export the diagram as plain arrays for downstream algorithms."""
        return MDDStruct(
            K=self.K,
            domain=list(self.domain),
            root=self.root,
            nnodes=[len(self.node[k]) for k in range(self.K)],
            node=[np.asarray(self.node[k], dtype=np.int64).reshape(len(self.node[k]),
                                                                   self.domain[k])
                  for k in range(self.K)],
        )

    def stats(self) -> Dict[str, object]:
        """Storage description of the current set; only reachable nodes count."""
        vis = self._reachable_ids()
        nodes_per_level = [int(np.count_nonzero(v)) for v in vis]
        num_nodes = int(sum(nodes_per_level))
        num_states = self.cardinality()
        mdd_ints = int(sum(nodes_per_level[k] * self.domain[k] for k in range(self.K)))
        explicit_ints = num_states * self.K
        return {
            'levels': self.K,
            'nodesPerLevel': nodes_per_level,
            'numNodes': num_nodes,
            'liveNodes': num_nodes,
            'tableNodes': int(sum(len(self.node[k]) for k in range(self.K))),
            'numStates': num_states,
            'mddInts': mdd_ints,
            'explicitInts': explicit_ints,
            'compression': explicit_ints / max(mdd_ints, 1),
        }

    def compact(self) -> None:
        """Reclaim dead nodes left by the append-only build.

        Membership/index/enumerate are unchanged. A production MDD would
        reference-count instead and never accumulate dead nodes; this is the
        basic sweep.
        """
        vis = self._reachable_ids()
        newnode: List[List[List[int]]] = [[] for _ in range(self.K)]
        remap: List[Dict[int, int]] = [dict() for _ in range(self.K)]
        for k in range(self.K):
            ids = [r + 1 for r in range(len(self.node[k])) if vis[k][r]]
            for new_pos, old_id in enumerate(ids, start=1):
                remap[k][old_id] = new_pos
            newnode[k] = [list(self.node[k][old_id - 1]) for old_id in ids]
        for k in range(self.K - 1):
            for row in newnode[k]:
                for v in range(self.domain[k]):
                    if row[v] > 0:
                        row[v] = remap[k + 1][row[v]]
        self.node = newnode
        if self.root != TERM_FALSE:
            self.root = remap[0][self.root]
        self.uniq = [dict() for _ in range(self.K)]
        for k in range(self.K):
            for p, row in enumerate(self.node[k], start=1):
                self.uniq[k][tuple(row)] = p
        self._dirty = True

    def __repr__(self) -> str:
        s = self.stats()
        lines = [
            "  MDD  %d levels, domains [%s]" % (self.K, ' '.join(str(d) for d in self.domain)),
            "       %d states stored in %d nodes (%s per level)"
            % (s['numStates'], s['numNodes'], ' '.join(str(n) for n in s['nodesPerLevel'])),
            "       footprint %d ints vs %d explicit (%.1fx compression)"
            % (s['mddInts'], s['explicitInts'], s['compression']),
        ]
        if s['tableNodes'] > s['liveNodes']:
            lines.append("       (%d dead nodes in tables; call compact() to reclaim)"
                         % (s['tableNodes'] - s['liveNodes']))
        return '\n'.join(lines)

    # -- internals -------------------------------------------------------

    def _make_node(self, k: int, arcs: List[int]) -> int:
        """Canonical node creation through the per-level unique table."""
        if all(a == TERM_FALSE for a in arcs):
            return TERM_FALSE
        key = tuple(arcs)
        found = self.uniq[k].get(key)
        if found is not None:
            return found
        self.node[k].append(list(arcs))
        node_id = len(self.node[k])
        self.uniq[k][key] = node_id
        return node_id

    def _add_state(self, k: int, node_id: int, state: Sequence[int]) -> int:
        """Recursively add one state below node_id at level k.

        Nodes are immutable and shared, so this rebuilds the path bottom-up
        rather than mutating in place.
        """
        if k >= self.K:
            return TERM_TRUE
        if node_id == TERM_FALSE:
            arcs = [TERM_FALSE] * self.domain[k]
        else:
            arcs = list(self.node[k][node_id - 1])
        v = int(state[k])
        arcs[v] = self._add_state(k + 1, arcs[v], state)
        return self._make_node(k, arcs)

    def _child_count(self, k: int, child_id: int) -> int:
        """Accepted states below a child reference at level k."""
        if k >= self.K:
            return 1 if child_id == TERM_TRUE else 0
        if child_id == TERM_FALSE:
            return 0
        return self._count_node(k, child_id)

    def _count_node(self, k: int, node_id: int) -> int:
        c = self._cnt[k][node_id - 1]
        if c >= 0:
            return c
        arcs = self.node[k][node_id - 1]
        c = 0
        for v in range(self.domain[k]):
            c += self._child_count(k + 1, arcs[v])
        self._cnt[k][node_id - 1] = c
        return c

    def _ensure_counts(self) -> None:
        if not self._dirty and all(len(self._cnt[k]) == len(self.node[k])
                                   for k in range(self.K)):
            return
        self._cnt = [[-1] * len(self.node[k]) for k in range(self.K)]
        self._dirty = False

    def _reachable_ids(self) -> List[List[bool]]:
        """Per-level masks of nodes reachable from the root."""
        vis = [[False] * len(self.node[k]) for k in range(self.K)]
        if self.root == TERM_FALSE:
            return vis
        vis[0][self.root - 1] = True
        stack = [(0, self.root)]
        while stack:
            k, node_id = stack.pop()
            if k == self.K - 1:
                continue
            arcs = self.node[k][node_id - 1]
            for v in range(self.domain[k]):
                ch = arcs[v]
                if ch > 0 and not vis[k + 1][ch - 1]:
                    vis[k + 1][ch - 1] = True
                    stack.append((k + 1, ch))
        return vis

    def _enum_below(self, k: int, node_id: int) -> List[List[int]]:
        """All states over levels k..K-1 reaching TRUE below node_id."""
        arcs = self.node[k][node_id - 1]
        if k == self.K - 1:
            return [[v] for v in range(self.domain[k]) if arcs[v] == TERM_TRUE]
        out: List[List[int]] = []
        for v in range(self.domain[k]):
            child = arcs[v]
            if child != TERM_FALSE:
                for sub in self._enum_below(k + 1, child):
                    out.append([v] + sub)
        return out


class MDDStruct:
    """Plain-array export of an MDD, the input contract of mdd_mcd.

    Fields mirror MDD.toStruct in MATLAB:
        K, domain, root - as the MDD properties
        nnodes          - per-level live node count
        node            - per-level (nnodes[k] x domain[k]) arrays of child ids
                          (level k+1 ids, or terminals at k = K-1)
    """

    __slots__ = ('K', 'domain', 'root', 'nnodes', 'node')

    def __init__(self, K, domain, root, nnodes, node):
        self.K = K
        self.domain = domain
        self.root = root
        self.nnodes = nnodes
        self.node = node
