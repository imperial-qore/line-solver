"""
Ensemble model utilities for native Python LINE.

Provides Ensemble.merge, which combines several independent Network models into
a single Network containing disconnected subnetworks. Node and class names are
prefixed with their model name to avoid conflicts, and all Source/Sink nodes are
merged into single MergedSource/MergedSink nodes.

Mirrors the JAR jline.lang.Ensemble.merge implementation.
"""

import warnings

import numpy as np


class Ensemble:
    """Static helpers for combining independent network models."""

    @staticmethod
    def merge(models):
        """Combine a list of independent Network models into a single Network.

        Node and class names are prefixed with their model name to prevent
        conflicts; all Source and Sink nodes are merged into single
        MergedSource and MergedSink nodes.

        Args:
            models: list of Network models to merge.

        Returns:
            A new Network containing the disconnected subnetworks.

        References:
            MATLAB: matlab/src/lang/Ensemble.m (merge)
            Java:   java/src/main/java/jline/lang/Ensemble.java (merge)
        """
        from .network import Network
        from .nodes import (Source, Sink, Delay, Queue, Router, ClassSwitch,
                            Fork, Join, Cache, Logger, Place, Transition)
        from .classes import OpenClass, ClosedClass, SelfLoopingClass
        from ..distributions import Disabled

        if not models:
            warnings.warn("Empty ensemble provided, returning empty Network")
            return Network("EmptyMergedNetwork")

        union_name = "_".join(m.get_name() for m in models)
        if len(union_name) > 50:
            union_name = "MergedNetwork_%d" % len(models)
        union = Network(union_name)

        has_open = any(m.has_open_classes() for m in models)
        union_source = Source(union, "MergedSource") if has_open else None
        union_sink = Sink(union, "MergedSink") if has_open else None

        node_maps = [dict() for _ in models]
        class_maps = [dict() for _ in models]
        join_pending = []  # (model_idx, old_join, new_name)

        # Phase 1: create nodes
        for mi, model in enumerate(models):
            prefix = model.get_name() + "_"
            nmap = node_maps[mi]
            for old in model.get_nodes():
                new_name = prefix + old.name
                if isinstance(old, Source):
                    nmap[old.name] = union_source
                    continue
                if isinstance(old, Sink):
                    nmap[old.name] = union_sink
                    continue
                if isinstance(old, Delay):
                    new = Delay(union, new_name)
                elif isinstance(old, Cache):
                    new = Cache(union, new_name, old._num_items,
                                old._item_level_cap, old._replacement_strategy)
                elif isinstance(old, Queue):
                    new = Queue(union, new_name, old.get_sched_strategy())
                    ns = getattr(old, 'number_of_servers', 1)
                    if ns is not None and ns != float('inf') and ns > 1:
                        new.set_number_of_servers(int(ns))
                elif isinstance(old, Router):
                    new = Router(union, new_name)
                elif isinstance(old, ClassSwitch):
                    if getattr(old, '_auto_added', False) or getattr(old, 'auto_added', False):
                        continue
                    new = ClassSwitch(union, new_name)
                elif isinstance(old, Fork):
                    new = Fork(union, new_name)
                    tpl = old.get_tasks_per_link()
                    if tpl is not None:
                        new.set_tasks_per_link(tpl)
                elif isinstance(old, Join):
                    join_pending.append((mi, old, new_name))
                    continue
                elif isinstance(old, Logger):
                    new = Logger(union, new_name)
                elif isinstance(old, Place):
                    new = Place(union, new_name)
                elif isinstance(old, Transition):
                    new = Transition(union, new_name)
                else:
                    warnings.warn("Unsupported node type: %s, skipping" % type(old).__name__)
                    continue
                nmap[old.name] = new

        # Phase 2: pending Join nodes (after Forks exist)
        for mi, old_join, new_name in join_pending:
            nmap = node_maps[mi]
            fork = getattr(old_join, '_fork', None)
            new_fork = nmap.get(fork.name) if fork is not None else None
            new = Join(union, new_name, new_fork) if new_fork is not None else Join(union, new_name)
            nmap[old_join.name] = new

        # Phase 3: create job classes
        for mi, model in enumerate(models):
            prefix = model.get_name() + "_"
            cmap = class_maps[mi]
            nmap = node_maps[mi]
            for old in model.get_classes():
                new_name = prefix + old.name
                if isinstance(old, SelfLoopingClass):
                    refname = old.get_reference_station().name
                    newref = nmap.get(refname)
                    if newref is None:
                        warnings.warn("Reference station %s not found for class %s" % (refname, old.name))
                        continue
                    new = SelfLoopingClass(union, new_name, old.getPopulation(), newref, old.priority)
                elif isinstance(old, ClosedClass):
                    refname = old.get_reference_station().name
                    newref = nmap.get(refname)
                    if newref is None:
                        warnings.warn("Reference station %s not found for class %s" % (refname, old.name))
                        continue
                    new = ClosedClass(union, new_name, old.getPopulation(), newref, old.priority)
                elif isinstance(old, OpenClass):
                    new = OpenClass(union, new_name, old.priority)
                else:
                    warnings.warn("Unknown class type: %s, skipping" % type(old).__name__)
                    continue
                cmap[old.name] = new

        # Phase 4: service and arrival distributions
        for mi, model in enumerate(models):
            nmap = node_maps[mi]
            cmap = class_maps[mi]
            for old in model.get_nodes():
                if isinstance(old, Source):
                    for oc in model.get_classes():
                        if isinstance(oc, OpenClass) and oc.name in cmap:
                            dist = old.get_arrival(oc)
                            if dist is not None and not isinstance(dist, Disabled):
                                union_source.set_arrival(cmap[oc.name], dist)
                    continue
                if isinstance(old, Sink):
                    continue
                if isinstance(old, ClassSwitch) and (getattr(old, '_auto_added', False) or getattr(old, 'auto_added', False)):
                    continue
                if old.name not in nmap:
                    continue
                new = nmap[old.name]
                if hasattr(old, 'get_service') and hasattr(new, 'set_service'):
                    for oc in model.get_classes():
                        if oc.name not in cmap:
                            continue
                        dist = old.get_service(oc)
                        if dist is not None and not isinstance(dist, Disabled):
                            new.set_service(cmap[oc.name], dist)

        # Phase 5: routing matrix, reconstructed from each model's rtnodes
        P = union.init_routing_matrix()
        for mi, model in enumerate(models):
            sn = model.get_struct()
            rtnodes = sn.rtnodes
            if rtnodes is None or getattr(rtnodes, 'size', 0) == 0:
                continue
            rtnodes = np.asarray(rtnodes)
            nodes = model.get_nodes()
            classes = model.get_classes()
            K = len(classes)
            nmap = node_maps[mi]
            cmap = class_maps[mi]
            for i, ni in enumerate(nodes):
                if ni.name not in nmap:
                    continue
                for r, cr in enumerate(classes):
                    if cr.name not in cmap:
                        continue
                    for j, nj in enumerate(nodes):
                        if nj.name not in nmap:
                            continue
                        for s, cs in enumerate(classes):
                            if cs.name not in cmap:
                                continue
                            prob = rtnodes[i * K + r, j * K + s]
                            if prob > 0:
                                P.set(cmap[cr.name], cmap[cs.name],
                                      nmap[ni.name], nmap[nj.name], float(prob))
        union.link(P)
        return union
