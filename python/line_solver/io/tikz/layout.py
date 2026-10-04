"""
Automatic layered layout for queueing network diagrams (Sugiyama style).

Port of jar/src/main/java/jline/io/tikz/TikZLayoutEngine.java. Nodes are
identified by their 0-based index in ``model.get_nodes()``, which is also the
row/column order of ``sn.connmatrix``.
"""

from ._jcompat import JAVA_DOUBLE_MAX_VALUE, node_kind


class TikZLayoutEngine:
    """Assigns layers, orders nodes within layers and computes (x, y) in cm."""

    def __init__(self, model, options):
        self.model = model
        self.options = options
        self.positions = {}
        self.layers = []

    def compute_layout(self):
        nodes = self.model.get_nodes()
        if len(nodes) == 0:
            return
        successors = self._build_successor_map(nodes)
        predecessors = self._build_predecessor_map(nodes, successors)
        self._assign_layers(nodes, predecessors)
        self._minimize_crossings(successors, predecessors)
        self._assign_coordinates()

    computeLayout = compute_layout

    def get_position(self, idx):
        return self.positions.get(idx)

    def get_all_positions(self):
        return dict(self.positions)

    def get_layers(self):
        return [list(layer) for layer in self.layers]

    getPosition = get_position
    getAllPositions = get_all_positions
    getLayers = get_layers

    def _build_successor_map(self, nodes):
        n = len(nodes)
        successors = {i: [] for i in range(n)}
        conn = _connmatrix(self.model)
        if conn is None:
            return successors
        for i in range(n):
            for j in range(n):
                if conn[i][j] > 0:
                    successors[i].append(j)
        return successors

    @staticmethod
    def _build_predecessor_map(nodes, successors):
        predecessors = {i: [] for i in range(len(nodes))}
        for frm, tos in successors.items():
            for to in tos:
                predecessors[to].append(frm)
        return predecessors

    def _assign_layers(self, nodes, predecessors):
        n = len(nodes)
        kinds = [node_kind(node) for node in nodes]
        self.layers = []
        assigned = set()

        layer0 = []
        for i in range(n):
            if kinds[i] == 'source' or len(predecessors[i]) == 0:
                layer0.append(i)
                assigned.add(i)
        if not layer0 and n > 0:
            layer0.append(0)
            assigned.add(0)
        self.layers.append(layer0)

        current = 0
        while len(assigned) < n:
            next_layer = [i for i in range(n)
                          if i not in assigned and all(p in assigned for p in predecessors[i])]
            if not next_layer:
                # cycle: take the first unassigned node
                for i in range(n):
                    if i not in assigned:
                        next_layer.append(i)
                        break
            if next_layer:
                assigned.update(next_layer)
                self.layers.append(next_layer)
            current += 1
            if current > n:
                break

        # sinks go to the last layer, then empty layers are dropped
        if len(self.layers) > 1:
            last = self.layers[-1]
            for layer in self.layers[:-1]:
                sinks = [i for i in layer if kinds[i] == 'sink']
                for i in sinks:
                    layer.remove(i)
                    if i not in last:
                        last.append(i)
            self.layers = [layer for layer in self.layers if layer]

    def _minimize_crossings(self, successors, predecessors):
        for _ in range(4):
            for i in range(1, len(self.layers)):
                self._reorder_by_barycenter(self.layers[i], predecessors, self.layers[i - 1])
            for i in range(len(self.layers) - 2, -1, -1):
                self._reorder_by_barycenter(self.layers[i], successors, self.layers[i + 1])

    @staticmethod
    def _reorder_by_barycenter(layer, connections, reference):
        ref_pos = {node: k for k, node in enumerate(reference)}
        bary = {}
        for node in layer:
            connected = connections.get(node)
            if not connected:
                bary[node] = JAVA_DOUBLE_MAX_VALUE
                continue
            total = 0.0
            count = 0
            for c in connected:
                pos = ref_pos.get(c)
                if pos is not None:
                    total += pos
                    count += 1
            bary[node] = total / count if count > 0 else JAVA_DOUBLE_MAX_VALUE
        layer.sort(key=lambda node: bary[node])  # stable, like Collections.sort

    def _assign_coordinates(self):
        self.positions = {}
        spacing = float(self.options.node_spacing)
        for layer_idx, layer in enumerate(self.layers):
            x = layer_idx * float(self.options.layer_spacing)
            start_y = (len(layer) - 1) * spacing / 2.0
            for k, node in enumerate(layer):
                self.positions[node] = (x, start_y - k * spacing)


def _connmatrix(model):
    """Return sn.connmatrix as a nested list, or None when the model has no struct/connections."""
    import numpy as np
    sn = model.get_struct()
    if sn is None or getattr(sn, 'connmatrix', None) is None:
        return None
    return np.asarray(sn.connmatrix, dtype=float).tolist()
