"""
Network visualization for the native implementation.

Port of the wrapper Network.getGraph/plot methods (python-wrapper
line_solver/_lang.py): builds station- and node-level directed-graph
dictionaries from the compiled NetworkStruct and renders them with
networkx/matplotlib. The TikZ export pipeline (toTikZ/plotTikZ/exportTikZ*)
remains JAR-backed and is intentionally not available natively.
"""

import numpy as np


def _nodetype_name(sn, ind):
    """Best-effort readable node type name for node ind."""
    try:
        nt = np.asarray(sn.nodetype).flatten()[ind]
    except (AttributeError, IndexError, TypeError):
        return 'Unknown'
    if hasattr(nt, 'name'):
        return nt.name
    try:
        from ..api.sn.network_struct import NodeType
        return NodeType(int(nt)).name
    except (ValueError, TypeError, ImportError):
        return str(nt)


def _nservers_at(sn, ind):
    ns = np.asarray(sn.nservers).flatten()
    if ind >= len(ns):
        return 0
    v = ns[ind]
    if np.isinf(v):
        return -1
    return int(v)


def network_get_graph(self):
    """Build (H, G) graph dictionaries for the model.

    H is the station-level graph (id/name/type/jobs/servers per node, edges
    with weight/rate/class), G is the node-level graph. Mirrors the wrapper
    Network.getGraph return format.
    """
    sn = self.get_struct()
    P = np.asarray(sn.rt) if sn.rt is not None else np.zeros((0, 0))
    Pnodes = np.asarray(sn.rtnodes) if sn.rtnodes is not None else np.zeros((0, 0))
    nodenames = list(sn.nodenames)
    classnames = list(sn.classnames)
    njobs = np.asarray(sn.njobs).flatten()
    refstat = np.asarray(sn.refstat).flatten()
    rates = np.asarray(sn.rates) if sn.rates is not None else np.zeros((0, 0))

    G = {'nodes': [], 'edges': []}
    for ind in range(sn.nnodes):
        G['nodes'].append({
            'id': ind,
            'name': nodenames[ind] if ind < len(nodenames) else 'Node%d' % ind,
            'type': _nodetype_name(sn, ind),
            'servers': _nservers_at(sn, ind),
        })
    for ind in range(sn.nnodes):
        for jnd in range(sn.nnodes):
            for k in range(sn.nclasses):
                i1 = ind * sn.nclasses + k
                i2 = jnd * sn.nclasses + k
                if i1 < Pnodes.shape[0] and i2 < Pnodes.shape[1] and Pnodes[i1, i2] > 0:
                    G['edges'].append({
                        'source': ind, 'target': jnd,
                        'weight': float(Pnodes[i1, i2]), 'class': k,
                    })

    H = {'nodes': [], 'edges': []}
    stationToNode = np.asarray(sn.stationToNode).flatten()
    for ist in range(sn.nstations):
        ind = int(stationToNode[ist]) if ist < len(stationToNode) else ist
        jobs = 0.0
        for k in range(sn.nclasses):
            # refstat holds the reference STATION of each class (0-based)
            if k < len(refstat) and int(refstat[k]) == ist and k < len(njobs) \
                    and np.isfinite(njobs[k]):
                jobs += njobs[k]
        H['nodes'].append({
            'id': ist,
            'name': nodenames[ind] if ind < len(nodenames) else 'Station%d' % ist,
            'type': _nodetype_name(sn, ind),
            'jobs': int(jobs),
            'servers': _nservers_at(sn, ist),
        })
    statefulToStation = np.asarray(sn.statefulToStation).flatten() \
        if getattr(sn, 'statefulToStation', None) is not None else None
    stationToStateful = np.asarray(sn.stationToStateful).flatten() \
        if getattr(sn, 'stationToStateful', None) is not None else None
    for ist in range(sn.nstations):
        for jst in range(sn.nstations):
            isf = int(stationToStateful[ist]) if stationToStateful is not None else ist
            jsf = int(stationToStateful[jst]) if stationToStateful is not None else jst
            for k in range(sn.nclasses):
                i1 = isf * sn.nclasses + k
                i2 = jsf * sn.nclasses + k
                if i1 < P.shape[0] and i2 < P.shape[1] and P[i1, i2] > 0:
                    H['edges'].append({
                        'source': ist, 'target': jst,
                        'weight': float(P[i1, i2]),
                        'rate': float(rates[ist, k])
                        if ist < rates.shape[0] and k < rates.shape[1] else 0.0,
                        'class': classnames[k] if k < len(classnames) else 'Class%d' % k,
                    })
    return H, G


def network_plot(self, graph_type='station', method='names', **kwargs):
    """Plot the network as a directed graph (matplotlib + networkx).

    Args:
        graph_type: 'station' or 'node' level graph (default 'station').
        method: node labeling, one of 'names', 'types', 'ids'.
        kwargs: figsize, node_color, node_size, font_size, font_weight,
                edge_color, arrowsize, title_fontsize, show.
    """
    try:
        import matplotlib.pyplot as plt
        import networkx as nx
    except ImportError:
        raise ImportError(
            "Matplotlib and NetworkX are required for plotting. "
            "Install with: pip install matplotlib networkx")

    H, G = self.get_graph()
    graph_data = H if graph_type == 'station' else G

    nx_graph = nx.DiGraph()
    node_labels = {}
    for node in graph_data['nodes']:
        nx_graph.add_node(node['id'])
        if method == 'names':
            node_labels[node['id']] = node['name']
        elif method == 'types':
            node_labels[node['id']] = node['type']
        else:
            node_labels[node['id']] = str(node['id'])
    for edge in graph_data['edges']:
        nx_graph.add_edge(edge['source'], edge['target'], weight=edge['weight'])

    plt.figure(figsize=kwargs.get('figsize', (12, 8)))
    try:
        pos = nx.nx_agraph.graphviz_layout(nx_graph, prog='dot')
    except Exception:
        pos = nx.spring_layout(nx_graph, k=5, iterations=100, seed=23000, scale=2)

    nx.draw(nx_graph, pos,
            with_labels=True,
            labels=node_labels,
            node_color=kwargs.get('node_color', 'lightblue'),
            node_size=kwargs.get('node_size', 1000),
            font_size=kwargs.get('font_size', 8),
            font_weight=kwargs.get('font_weight', 'bold'),
            arrows=True,
            edge_color=kwargs.get('edge_color', 'gray'),
            arrowsize=kwargs.get('arrowsize', 20))

    plt.title('%s Graph - %s' % (graph_type.capitalize(), self.get_name()),
              fontsize=kwargs.get('title_fontsize', 14))
    plt.axis('off')
    if kwargs.get('show', True):
        plt.show()
    return plt.gcf()
