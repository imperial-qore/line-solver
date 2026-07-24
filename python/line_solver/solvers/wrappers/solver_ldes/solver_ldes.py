"""
Native Python LDES solver using subprocess to call ldes.jar CLI.

This module provides the SolverLDES class that runs discrete event
simulation by invoking `java -jar ldes.jar solve model.json` and
parsing the JSON result.
"""

import json
import math
import os
import platform
import re
import shutil
import subprocess
import tempfile
import time
import urllib.request
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
from typing import Optional, Any, List, Tuple

from .ldes_options import LDESOptions, LDESResult
from ...base import NetworkSolver
from ....api.sn.transforms import sn_get_residt_from_respt
from ....api.sn.network_struct import NodeType
from ....constants import GlobalConstants
from ....io.linemodel_io import save_model


# ---------------------------------------------------------------------------
# JSIMG/JSIM/JSIMW → Python native Network conversion
# ---------------------------------------------------------------------------

def _jsim_file_to_network(filename: str):
    """
    Parse a JMT jsimg/jsim/jsimw file and construct a Python native Network.

    Mirrors the logic of M2M.JSIM2LINE() in the JAR implementation, using
    Python-native Network/node/class/distribution APIs.

    Args:
        filename: Path to a JMT file (.jsimg, .jsim, .jsimw, or .jmva).

    Returns:
        A fully constructed Network object ready for simulation.
    """
    from ....lang.network import Network
    from ....lang.nodes import Source, Queue, Delay, Sink, Router, Fork, Join
    from ....lang.classes import OpenClass, ClosedClass
    from ....lang.base import SchedStrategy, RoutingStrategy
    from ....distributions.continuous import (
        Exp, Erlang, HyperExp, Det, Gamma, Pareto, Lognormal, Uniform, Weibull, Disabled
    )
    from ....distributions.markovian import Coxian, APH, MAP, MMPP2

    ext = os.path.splitext(filename)[1].lower()
    if ext == '.jmva':
        raise ValueError(
            f"JMVA file format is not yet supported for direct LDES loading: {filename}. "
            "Convert to JSIMG first using JMT."
        )

    tree = ET.parse(filename)
    root = tree.getroot()

    # Get <sim> element (may be direct root or child of <archive>)
    sim = root.find('sim')
    if sim is None:
        sim = root

    # Derive model name from archive name attribute
    arch_name = root.get('name', os.path.basename(filename))
    model_name = os.path.splitext(os.path.basename(arch_name))[0]
    model = Network(model_name)

    # -------------------------------------------------------------------------
    # Step 1: Collect node elements and their sections
    # -------------------------------------------------------------------------
    node_elems = sim.findall('node')
    orig_names: List[str] = [ne.get('name', f'Node{i}') for i, ne in enumerate(node_elems)]
    san_names: List[str] = [n.replace('/', '_').replace('\\', '_') for n in orig_names]
    sections_list: List[List[ET.Element]] = [ne.findall('section') for ne in node_elems]

    # -------------------------------------------------------------------------
    # Step 2: Create node objects (first pass)
    # -------------------------------------------------------------------------
    node_objs: List[Optional[object]] = [None] * len(orig_names)
    fork_stack: List[Fork] = []

    for i, (name, sects) in enumerate(zip(san_names, sections_list)):
        if not sects:
            continue
        first_cls = sects[0].get('className', '')

        if first_cls == 'JobSink':
            node_objs[i] = Sink(model, name)

        elif first_cls == 'RandomSource':
            node_objs[i] = Source(model, name)

        elif first_cls == 'Join':
            fork_obj = fork_stack[-1] if fork_stack else None
            node_objs[i] = Join(model, name, fork_obj)

        elif first_cls == 'Storage':
            # Place nodes (Petri nets) — not parsed here
            pass

        elif first_cls == 'Enabling':
            # Transition nodes (Petri nets) — not parsed here
            pass

        elif first_cls == 'Queue':
            second_cls = sects[1].get('className', '') if len(sects) > 1 else ''
            third_cls = sects[2].get('className', '') if len(sects) > 2 else ''

            if third_cls == 'Fork':
                f = Fork(model, name)
                fork_stack.append(f)
                node_objs[i] = f

            elif second_cls == 'ServiceTunnel':
                node_objs[i] = Router(model, name)

            elif second_cls == 'Delay':
                d = Delay(model, name)
                cap = _jsim_parse_queue_size(sects[0])
                if cap is not None and cap > 0:
                    d.set_capacity(cap)
                node_objs[i] = d

            elif second_cls in ('Server', 'PSServer'):
                strat, nservers = _jsim_parse_server_info(sects, second_cls)
                q = Queue(model, name, strat)
                cap = _jsim_parse_queue_size(sects[0])
                if cap is not None and cap > 0:
                    q.set_capacity(cap)
                if nservers > 1:
                    q.set_number_of_servers(nservers)
                node_objs[i] = q

            else:
                # Fallback: generic FCFS queue
                node_objs[i] = Queue(model, name, SchedStrategy.FCFS)

    # -------------------------------------------------------------------------
    # Step 3: Create job classes
    # -------------------------------------------------------------------------
    class_elems = sim.findall('userClass')
    jmt_prios = [int(ce.get('priority', '0')) for ce in class_elems]
    max_prio = max(jmt_prios) if jmt_prios else 0

    job_classes: List[object] = []
    for ce, jmt_prio in zip(class_elems, jmt_prios):
        cname = ce.get('name', '')
        ctype = ce.get('type', 'open').lower()
        line_prio = max_prio - jmt_prio  # invert JMT convention

        if ctype == 'open':
            c = OpenClass(model, cname, line_prio)
        else:
            ref_src = ce.get('referenceSource', '')
            customers_str = ce.get('customers', ce.get('population', '1'))
            customers = int(customers_str) if customers_str else 1
            ref_idx = orig_names.index(ref_src) if ref_src in orig_names else 0
            ref_node = node_objs[ref_idx] if 0 <= ref_idx < len(node_objs) else None
            c = ClosedClass(model, cname, customers, ref_node, line_prio)
        job_classes.append(c)

    # -------------------------------------------------------------------------
    # Step 4: Set service / arrival distributions
    # -------------------------------------------------------------------------
    for i, (node, sects) in enumerate(zip(node_objs, sections_list)):
        if node is None or not sects:
            continue

        first_cls = sects[0].get('className', '')

        if isinstance(node, Source):
            # RandomSource section: ServiceStrategy parameter contains arrivals
            _jsim_set_source_arrivals(node, sects[0], job_classes, Disabled)

        elif isinstance(node, (Queue, Delay)):
            # Server/PSServer/Delay section at index 1
            if len(sects) > 1:
                _jsim_set_service_distributions(node, sects[1], job_classes, Disabled)

    # -------------------------------------------------------------------------
    # Step 5: Create connections (topology)
    # -------------------------------------------------------------------------
    for ce in sim.findall('connection'):
        src_name = ce.get('source', '')
        tgt_name = ce.get('target', '')
        si = orig_names.index(src_name) if src_name in orig_names else -1
        ti = orig_names.index(tgt_name) if tgt_name in orig_names else -1
        if si >= 0 and ti >= 0 and node_objs[si] is not None and node_objs[ti] is not None:
            model.add_link(node_objs[si], node_objs[ti])

    # -------------------------------------------------------------------------
    # Step 6: Set routing strategies
    # -------------------------------------------------------------------------
    for i, (node, sects) in enumerate(zip(node_objs, sections_list)):
        if node is None or isinstance(node, Sink) or not sects:
            continue

        first_cls = sects[0].get('className', '')

        # Routing section is the last section for Source/Queue nodes
        if first_cls == 'RandomSource' and len(sects) > 2:
            routing_sect = sects[2]
        elif first_cls == 'Queue' and len(sects) > 2:
            routing_sect = sects[2]
        else:
            continue

        _jsim_set_node_routing(node, routing_sect, job_classes, orig_names, node_objs)

    return model


def _jsim_parse_queue_size(queue_sect: ET.Element) -> Optional[int]:
    """Parse the queue buffer capacity from a Queue section element."""
    for param in queue_sect.findall('parameter'):
        if param.get('name', '') == 'size':
            val_elem = param.find('value')
            if val_elem is not None and val_elem.text:
                try:
                    v = int(val_elem.text)
                    return v if v > 0 else None
                except ValueError:
                    pass
    return None


def _jsim_parse_server_info(sects: List[ET.Element], second_cls: str) -> Tuple[Any, int]:
    """
    Determine scheduling strategy and number of servers for a Queue node.

    Returns:
        (SchedStrategy, nservers)
    """
    from ....lang.base import SchedStrategy

    strategy = SchedStrategy.FCFS
    nservers = 1

    if second_cls == 'Server':
        # Scheduling from QueuePutStrategy in section[0]
        for param in sects[0].findall('parameter'):
            if param.get('name', '') == 'QueuePutStrategy':
                subs = param.findall('subParameter')
                if subs:
                    strategy = _jsim_put_strategy_to_sched(subs[0].get('name', 'TailStrategy'))
                break
        # Number of servers from maxJobs in section[1] (Server)
        for param in sects[1].findall('parameter'):
            if param.get('name', '') == 'maxJobs':
                val = param.find('value')
                if val is not None and val.text:
                    try:
                        nservers = int(val.text)
                    except ValueError:
                        pass
                break

    elif second_cls == 'PSServer':
        strategy = _jsim_parse_ps_strategy(sects[1])
        for param in sects[1].findall('parameter'):
            if param.get('name', '') == 'maxJobs':
                val = param.find('value')
                if val is not None and val.text:
                    try:
                        nservers = int(val.text)
                    except ValueError:
                        pass
                break

    return strategy, nservers


def _jsim_put_strategy_to_sched(name: str) -> Any:
    """Map a JMT QueuePutStrategy name to a LINE SchedStrategy."""
    from ....lang.base import SchedStrategy
    mapping = {
        'TailStrategy': SchedStrategy.FCFS,
        'TailStrategyPriority': SchedStrategy.HOL,
        'HeadStrategy': SchedStrategy.LCFS,
        'RandStrategy': SchedStrategy.SIRO,
        'SJFStrategy': SchedStrategy.SJF,
        'SEPTStrategy': SchedStrategy.SEPT,
        'LJFStrategy': SchedStrategy.LJF,
        'LEPTStrategy': SchedStrategy.LEPT,
    }
    return mapping.get(name, SchedStrategy.FCFS)


def _jsim_parse_ps_strategy(server_sect: ET.Element) -> Any:
    """Parse PS/DPS/GPS strategy from a PSServer section."""
    from ....lang.base import SchedStrategy
    for sp in server_sect.iter('subParameter'):
        name = sp.get('name', '')
        if name == 'EPSStrategy':
            return SchedStrategy.PS
        elif name == 'DPSStrategy':
            return SchedStrategy.DPS
        elif name == 'GPSStrategy':
            return SchedStrategy.GPS
    return SchedStrategy.PS


def _jsim_parse_distribution(subs: List[ET.Element]) -> Any:
    """
    Parse a service/arrival distribution from the two subParameter children of a
    ServiceTimeStrategy element.

    Args:
        subs: List of subParameter elements — [distribution_type_elem, distrPar_elem].

    Returns:
        A distribution object or None if not parsable (caller should use Disabled).
    """
    from ....distributions.continuous import (
        Exp, Erlang, HyperExp, Det, Gamma, Pareto, Lognormal, Uniform, Weibull, Disabled
    )
    from ....distributions.markovian import Coxian

    if not subs:
        return None

    distr_name = subs[0].get('name', '')

    if distr_name in ('ZeroServiceTimeStrategy', 'Disabled', ''):
        return None  # caller sets Disabled

    if len(subs) < 2:
        return None

    distr_params = subs[1].findall('subParameter')

    def _val(idx: int) -> Optional[float]:
        if idx < len(distr_params):
            v = distr_params[idx].find('value')
            if v is not None and v.text:
                try:
                    return float(v.text)
                except ValueError:
                    pass
        return None

    if distr_name == 'Exponential':
        lam = _val(0)
        return Exp(lam) if lam is not None else None

    elif distr_name == 'Erlang':
        lam = _val(0)
        phases = _val(1)
        return Erlang(lam, int(phases)) if lam is not None and phases is not None else None

    elif distr_name == 'Hyperexponential':
        p = _val(0)
        r1 = _val(1)
        r2 = _val(2)
        return HyperExp(p, r1, r2) if None not in (p, r1, r2) else None

    elif distr_name == 'Coxian':
        lam = _val(0)
        lam1 = _val(1)
        phi = _val(2)
        return Coxian([lam, lam1], [phi, 1.0]) if None not in (lam, lam1, phi) else None

    elif distr_name == 'Deterministic':
        val = _val(0)
        return Det(val) if val is not None else None

    elif distr_name == 'Pareto':
        alpha = _val(0)
        k = _val(1)
        return Pareto(alpha, k) if None not in (alpha, k) else None

    elif distr_name == 'Gamma':
        shape = _val(0)
        scale = _val(1)
        return Gamma(shape, scale) if None not in (shape, scale) else None

    elif distr_name == 'Weibull':
        # Java M2M note: "scale and shape are inverted in the constructor"
        # JMT par[0]=shape(lambda), par[1]=scale(lambda1) → Weibull(scale, shape)
        shape = _val(0)
        scale = _val(1)
        return Weibull(scale, shape) if None not in (shape, scale) else None

    elif distr_name == 'Lognormal':
        mu = _val(0)
        sigma = _val(1)
        return Lognormal(mu, sigma) if None not in (mu, sigma) else None

    elif distr_name == 'Uniform':
        a = _val(0)
        b = _val(1)
        return Uniform(a, b) if None not in (a, b) else None

    return None


def _jsim_set_source_arrivals(source_node, source_sect: ET.Element,
                               job_classes: list, Disabled) -> None:
    """Set arrival distributions for a Source node from its RandomSource section."""
    for param in source_sect.findall('parameter'):
        if param.get('name', '') == 'ServiceStrategy':
            subparams = param.findall('subParameter')  # one ServiceTimeStrategy per class
            for r, cls in enumerate(job_classes):
                if r < len(subparams):
                    sp = subparams[r]
                    subs = sp.findall('subParameter')  # [dist_type, distrPar]
                    dist = _jsim_parse_distribution(subs)
                    source_node.set_arrival(cls, dist if dist is not None else Disabled())
            break


def _jsim_set_service_distributions(queue_node, srv_sect: ET.Element,
                                     job_classes: list, Disabled) -> None:
    """Set service distributions for a Queue or Delay node from its Server/Delay section."""
    for param in srv_sect.findall('parameter'):
        if param.get('name', '') == 'ServiceStrategy':
            subparams = param.findall('subParameter')  # one ServiceTimeStrategy per class
            for r, cls in enumerate(job_classes):
                if r < len(subparams):
                    sp = subparams[r]
                    subs = sp.findall('subParameter')  # [dist_type, distrPar]
                    dist = _jsim_parse_distribution(subs)
                    queue_node.set_service(cls, dist if dist is not None else Disabled())
            break


def _jsim_set_node_routing(node, routing_sect: ET.Element, job_classes: list,
                            orig_names: List[str], node_objs: list) -> None:
    """Set routing strategies for all classes on a node from its Router section."""
    from ....lang.base import RoutingStrategy

    for param in routing_sect.findall('parameter'):
        if param.get('name', '') == 'RoutingStrategy':
            subparams = param.findall('subParameter')  # one per class
            for r, cls in enumerate(job_classes):
                if r < len(subparams):
                    sp = subparams[r]
                    strat_name = sp.get('name', 'Random')
                    if strat_name == 'Random':
                        node.set_routing(cls, RoutingStrategy.RAND)
                    elif strat_name == 'Probabilities':
                        node.set_routing(cls, RoutingStrategy.PROB)
                        _jsim_set_prob_routing(node, cls, sp, orig_names, node_objs)
                    elif strat_name == 'Round Robin':
                        node.set_routing(cls, RoutingStrategy.RROBIN)
                    elif strat_name == 'Weighted Round Robin':
                        node.set_routing(cls, RoutingStrategy.WRROBIN)
                        _jsim_set_wrrobin_routing(node, cls, sp, orig_names, node_objs)
                    elif strat_name in ('Join the Shortest Queue (JSQ)',
                                        'Shortest Response Time', 'Fastest Service'):
                        node.set_routing(cls, RoutingStrategy.JSQ)
                    elif strat_name == 'Disabled':
                        node.set_routing(cls, RoutingStrategy.DISABLED)
                    else:
                        node.set_routing(cls, RoutingStrategy.RAND)
            break


def _jsim_set_prob_routing(node, cls, prob_param: ET.Element,
                            orig_names: List[str], node_objs: list) -> None:
    """Parse probabilistic routing destinations and set them on the node."""
    # XML structure under the Probabilities subParameter:
    # <subParameter name="Probabilities">
    #   <subParameter classPath="Probabilities_class">
    #     <subParameter name="dest1">
    #       <value>dest1_name</value>
    #       <value>probability</value>
    #     </subParameter>
    #     ...
    #   </subParameter>
    # </subParameter>
    outer = prob_param.find('subParameter')
    if outer is None:
        return
    for dest_sp in outer.findall('subParameter'):
        values = dest_sp.findall('value')
        if len(values) >= 2:
            dest_name = values[0].text or ''
            try:
                prob = float(values[1].text or '0')
            except ValueError:
                continue
            # Sanitize name
            dest_san = dest_name.replace('/', '_').replace('\\', '_')
            dest_idx = -1
            if dest_name in orig_names:
                dest_idx = orig_names.index(dest_name)
            elif dest_san in [n.replace('/', '_').replace('\\', '_') for n in orig_names]:
                for j, n in enumerate(orig_names):
                    if n.replace('/', '_').replace('\\', '_') == dest_san:
                        dest_idx = j
                        break
            if dest_idx >= 0 and node_objs[dest_idx] is not None:
                node.set_prob_routing(cls, node_objs[dest_idx], prob)


def _jsim_set_wrrobin_routing(node, cls, wrr_param: ET.Element,
                               orig_names: List[str], node_objs: list) -> None:
    """Parse weighted round-robin weights and set them on the node."""
    from ....lang.base import RoutingStrategy
    outer = wrr_param.find('subParameter')
    if outer is None:
        return
    for dest_sp in outer.findall('subParameter'):
        values = dest_sp.findall('value')
        if len(values) >= 2:
            dest_name = values[0].text or ''
            try:
                weight = float(values[1].text or '1')
            except ValueError:
                continue
            dest_idx = orig_names.index(dest_name) if dest_name in orig_names else -1
            if dest_idx >= 0 and node_objs[dest_idx] is not None:
                node.set_routing(cls, RoutingStrategy.WRROBIN,
                                 node_objs[dest_idx], weight)


# ---------------------------------------------------------------------------

class SolverLDES(NetworkSolver):
    """
    Native Python LINE Discrete Event Simulator (LDES) solver.

    Runs simulation by calling ldes.jar via subprocess and parsing
    the JSON result. This ensures Python native uses the same
    simulation engine as MATLAB and python-wrapper.

    Args:
        model: Network model (Python native)
        options: LDESOptions configuration (optional)
        **kwargs: Additional solver options (seed, samples, verbose, etc.)

    Example:
        >>> from line_solver.solvers.solver_ldes import SolverLDES, LDESOptions
        >>> options = LDESOptions(seed=23000, samples=200000)
        >>> solver = SolverLDES(model, options)
        >>> solver.runAnalyzer()
        >>> table = solver.getAvgTable()

    Passing an auxiliary solver in place of the options warm-starts the
    simulation from that solver's steady-state solution (see initFromSolver):
        >>> solver = SolverLDES(model, SolverMVA(model), samples=50000)
    """

    def __init__(
        self,
        model: Any,
        options: Optional[LDESOptions] = None,
        **kwargs
    ):
        # Accept a JMT file path (.jsimg/.jsim/.jsimw) in place of a Network object
        if isinstance(model, str):
            model = _jsim_file_to_network(model)

        # An auxiliary solver passed in place of the options requests a warm
        # start: its steady-state distribution decides the initial simulation
        # state (see initFromSolver).
        init_solver = None
        if options is not None and not isinstance(options, LDESOptions) \
                and hasattr(options, 'getAvgQLen'):
            init_solver = options
            options = None

        self.model = model
        self.options = options or LDESOptions()

        # Apply kwargs to options
        for key, value in kwargs.items():
            if hasattr(self.options, key):
                setattr(self.options, key, value)

        # options.events (DES event budget) overrides options.samples when set;
        # samples remains accepted as a deprecated alias for the event budget.
        if getattr(self.options, 'events', None):
            self.options.samples = int(self.options.events)

        self._result: Optional[LDESResult] = None

        # Extract network structure
        self._sn = self._get_network_struct()

        # Station and class names for table output
        self._station_names: list = []
        self._class_names: list = []
        self._extract_names()

        if init_solver is not None:
            self.initFromSolver(init_solver)

    def _get_network_struct(self) -> Any:
        """Get NetworkStruct from model."""
        model = self.model

        if hasattr(model, '_sn') and model._sn is not None:
            return model._sn

        if hasattr(model, 'refresh_struct'):
            model.refresh_struct()
            if hasattr(model, '_sn') and model._sn is not None:
                return model._sn

        if hasattr(model, 'getStruct'):
            try:
                return model.getStruct()
            except Exception:
                pass

        if hasattr(model, 'obj'):
            try:
                return model.getStruct()
            except Exception:
                pass

        if hasattr(model, 'nclasses') and hasattr(model, 'nstations'):
            return model

        raise ValueError("Cannot extract network structure from model")

    def _extract_names(self) -> None:
        """Extract station and class names from network structure."""
        sn = self._sn
        num_stations = int(sn.nstations) if hasattr(sn, 'nstations') else 1

        nodenames = list(sn.nodenames) if hasattr(sn, 'nodenames') and sn.nodenames is not None else []
        stationToNode = sn.stationToNode if hasattr(sn, 'stationToNode') and sn.stationToNode is not None else None

        if stationToNode is not None and nodenames:
            stationToNode = np.asarray(stationToNode).flatten()
            self._station_names = []
            for i in range(num_stations):
                if i < len(stationToNode):
                    node_idx = int(stationToNode[i])
                    if node_idx < len(nodenames):
                        self._station_names.append(nodenames[node_idx])
                    else:
                        self._station_names.append(f'Station{i}')
                else:
                    self._station_names.append(f'Station{i}')
        elif nodenames:
            self._station_names = nodenames[:num_stations]
            while len(self._station_names) < num_stations:
                self._station_names.append(f'Station{len(self._station_names)}')
        else:
            self._station_names = [f'Station{i}' for i in range(num_stations)]

        if hasattr(sn, 'classnames') and sn.classnames is not None:
            try:
                self._class_names = list(sn.classnames)
            except Exception:
                self._class_names = [f'Class{i}' for i in range(sn.nclasses)]
        else:
            num_classes = int(sn.nclasses) if hasattr(sn, 'nclasses') else 1
            self._class_names = [f'Class{i}' for i in range(num_classes)]

    def initFromSolver(self, init_solver: Any) -> 'SolverLDES':
        """Warm-start the simulation from the steady-state solution of an
        auxiliary solver.

        If the auxiliary solver is a SolverCTMC, the exact stationary
        distribution over the aggregate state space is computed and the
        initial state is set to the mode of that distribution (the most
        probable aggregate state). For any other network solver, the
        steady-state mean queue lengths are used instead and rounded to an
        integer placement that conserves each closed-class population.

        Since the simulation starts (approximately) in steady state, the
        transient removal filter is disabled (tranfilter='fixed' with
        warmupfrac=0), so every simulated sample contributes to the
        estimators.

        Args:
            init_solver: auxiliary solver used to compute the steady-state
                distribution (e.g. SolverCTMC or SolverMVA on the same model)

        Returns:
            self, for chaining
        """
        from ...warmstart import warm_start_placement

        placement = warm_start_placement(init_solver, self._sn)

        # LDES consumes the placement through options.init_sol (station-major
        # vector) rather than the model initial state.
        self.options.init_sol = placement.flatten()
        self.options.tranfilter = 'fixed'
        self.options.warmupfrac = 0.0
        return self

    @staticmethod
    def _get_package_bin_dir() -> str:
        """Get path to bin/ directory bundled inside the installed package."""
        package_dir = os.path.dirname(os.path.abspath(__file__))
        # Navigate: solver_ldes -> wrappers -> solvers -> line_solver
        line_solver_dir = os.path.dirname(os.path.dirname(
            os.path.dirname(package_dir)))
        return os.path.join(line_solver_dir, 'bin')

    @staticmethod
    def _get_common_dir() -> str:
        """Get path to common/ directory in the repository root."""
        package_dir = os.path.dirname(os.path.abspath(__file__))
        # Navigate: solver_ldes -> wrappers -> solvers -> line_solver -> python
        python_dir = os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.dirname(package_dir))))
        root_dir = os.path.dirname(python_dir)
        return os.path.join(root_dir, 'common')

    @staticmethod
    def _native_candidate_paths() -> List[str]:
        """Candidate locations for the native LDES binary, in preference order.

        1. Bundled in package (line_solver/bin/ldes) for pip-installed wheels
        2. Repository common/ directory for development checkouts
        """
        return [
            os.path.join(SolverLDES._get_package_bin_dir(), 'ldes'),
            os.path.join(SolverLDES._get_common_dir(), 'ldes'),
        ]

    # ELF e_machine identifiers (ELF header offset 0x12, 2 bytes)
    _EM_386 = 0x03
    _EM_ARM = 0x28
    _EM_X86_64 = 0x3E
    _EM_AARCH64 = 0xB7

    @staticmethod
    def _elf_machine(path: str) -> Optional[int]:
        """Return the ELF e_machine value of a binary, or None if not an ELF."""
        try:
            with open(path, 'rb') as f:
                hdr = f.read(20)
        except OSError:
            return None
        if len(hdr) < 20 or hdr[:4] != b'\x7fELF':
            return None
        # EI_DATA (index 5): 1 = little-endian, 2 = big-endian.
        endian = 'little' if hdr[5] == 1 else 'big'
        return int.from_bytes(hdr[18:20], endian)

    @staticmethod
    def _host_elf_machine() -> Optional[int]:
        """Expected ELF e_machine for the current host CPU, or None if unknown."""
        machine = platform.machine().lower()
        if machine in ('x86_64', 'amd64'):
            return SolverLDES._EM_X86_64
        if machine in ('aarch64', 'arm64'):
            return SolverLDES._EM_AARCH64
        if machine in ('i386', 'i486', 'i586', 'i686', 'x86'):
            return SolverLDES._EM_386
        if machine.startswith('arm'):
            return SolverLDES._EM_ARM
        return None

    @staticmethod
    def _get_ldes_native_path() -> Optional[str]:
        """Get path to a runnable native LDES binary, or None if none is usable.

        The native binary is only used on Linux (macOS/Windows run via
        java -jar ldes.jar). A binary whose ELF architecture does not match the
        host CPU (e.g. an x86-64 build on an aarch64 host) is ignored so the
        caller can fall back to ldes.jar. If the host or binary architecture
        cannot be determined, the binary is used as a best effort.
        """
        if platform.system() != 'Linux':
            return None

        # Escape hatch: force the ldes.jar path (e.g. when the prebuilt native
        # binary is stale relative to the current jline sources).
        if os.environ.get('LINE_LDES_FORCE_JAR'):
            return None

        host_machine = SolverLDES._host_elf_machine()
        for native_path in SolverLDES._native_candidate_paths():
            if not (os.path.isfile(native_path) and os.access(native_path, os.X_OK)):
                continue
            bin_machine = SolverLDES._elf_machine(native_path)
            if (host_machine is not None and bin_machine is not None
                    and bin_machine != host_machine):
                # Present but built for a different CPU architecture: skip it.
                continue
            return native_path

        return None

    @staticmethod
    def _incompatible_native_present() -> bool:
        """True if a native ldes binary exists but targets a different CPU."""
        if platform.system() != 'Linux':
            return False
        host_machine = SolverLDES._host_elf_machine()
        if host_machine is None:
            return False
        for native_path in SolverLDES._native_candidate_paths():
            if not os.path.isfile(native_path):
                continue
            bin_machine = SolverLDES._elf_machine(native_path)
            if bin_machine is not None and bin_machine != host_machine:
                return True
        return False

    @staticmethod
    def _find_java() -> Optional[str]:
        """Return a runnable java executable, or None if no JVM is available.

        Resolution order: $LINE_JAVA, then $JAVA_HOME/bin/java, then `java` on
        PATH. Mirrors the JAR dispatch convention so a single JVM setting works
        across the native-Python solvers.
        """
        cand = os.environ.get('LINE_JAVA')
        if cand:
            if os.path.isfile(cand) and os.access(cand, os.X_OK):
                return cand
            return None
        java_home = os.environ.get('JAVA_HOME')
        if java_home:
            cand = os.path.join(java_home, 'bin', 'java')
            if os.path.isfile(cand) and os.access(cand, os.X_OK):
                return cand
        return shutil.which('java')

    @staticmethod
    def _no_backend_message() -> str:
        """Human-readable error when neither a native binary nor a JVM is usable."""
        host = "%s/%s" % (platform.system(), platform.machine())
        if SolverLDES._incompatible_native_present():
            reason = (
                "A native LDES binary is present but was built for a different "
                "CPU architecture than this host (%s), and no Java runtime was "
                "found to run ldes.jar instead." % host
            )
        else:
            reason = (
                "No native LDES binary is available for this platform (%s), and "
                "no Java runtime was found to run ldes.jar." % host
            )
        return (
            "Cannot run the LDES solver. %s\n"
            "Please install a Java runtime (JRE/JDK 8 or newer), for example:\n"
            "  Debian/Ubuntu: sudo apt-get install default-jre\n"
            "  Fedora/RHEL:   sudo dnf install java-latest-openjdk\n"
            "  macOS:         brew install openjdk\n"
            "  Windows:       install from https://adoptium.net/\n"
            "Then ensure `java` is on your PATH, or set $LINE_JAVA or $JAVA_HOME."
            % reason
        )

    @staticmethod
    def _get_ldes_jar_path() -> str:
        """Get path to ldes.jar, downloading if necessary.

        Lookup order:
        1. Bundled in package (line_solver/bin/ldes.jar)
        2. Repository common/ directory
        3. Auto-download from SourceForge
        """
        # 1. Check package bin/ directory (pip-installed wheel)
        bin_dir = SolverLDES._get_package_bin_dir()
        jar_path = os.path.join(bin_dir, 'ldes.jar')
        if os.path.isfile(jar_path):
            return jar_path

        # 2. Check common/ directory (development/repo checkout)
        common_dir = SolverLDES._get_common_dir()
        ldes_path = os.path.join(common_dir, 'ldes.jar')
        if os.path.isfile(ldes_path):
            return ldes_path

        # 3. Try to download to common/
        os.makedirs(common_dir, exist_ok=True)
        ldes_url = 'https://line-solver.sourceforge.net/latest/ldes.jar'
        try:
            print(f"ldes.jar not found in {common_dir}")
            print("Attempting to download ldes.jar...")
            urllib.request.urlretrieve(ldes_url, ldes_path)
            print(f"Successfully downloaded ldes.jar to {common_dir}")
            return ldes_path
        except Exception as e:
            raise RuntimeError(
                f"ldes.jar not found and download failed: {e}\n"
                f"Please manually download from {ldes_url} and place in {common_dir}"
            )

    def _build_cli_args(self, model_path: str, result_path: str,
                         trajectory: bool = False,
                         export_histogram: bool = False) -> list:
        """Build CLI command arguments for LDES solver.

        Prefers native binary (common/ldes) for faster startup;
        falls back to java -jar ldes.jar if native binary not found.
        """
        return self._build_cli_runners(model_path, result_path,
                                       trajectory=trajectory,
                                       export_histogram=export_histogram)[0]

    def _build_cli_runners(self, model_path: str, result_path: str,
                           trajectory: bool = False,
                           export_histogram: bool = False) -> list:
        """Ordered CLI command candidates: native binary first, then ldes.jar.

        The native GraalVM binary may lack reflective features the JVM jar
        has (e.g. the fork-join MMT transformation deep-copies the model via
        Java serialization, unsupported in the AOT image), so the caller runs
        each candidate in order and falls through on failure, mirroring the
        MATLAB solveCli runner list.
        """
        opts = self.options

        native_path = self._get_ldes_native_path()
        # The --initsol flag was added after the prebuilt GraalVM binary was
        # cut, so a warm-start placement forces the ldes.jar path until the
        # native binary is rebuilt.
        if getattr(opts, 'init_sol', None) is not None:
            native_path = None
        flags = self._build_flag_args(trajectory=trajectory,
                                      export_histogram=export_histogram)
        runners = []
        if native_path is not None:
            runners.append([native_path, 'solve', model_path, '-o', result_path] + flags)
        java_exe = self._find_java()
        if java_exe is not None:
            ldes_jar = self._get_ldes_jar_path()
            runners.append([java_exe, '-jar', ldes_jar, 'solve', model_path,
                            '-o', result_path] + flags)
        if not runners:
            # Neither a usable native binary (absent, wrong CPU architecture,
            # or a non-Linux host) nor a JVM for ldes.jar: tell the user to
            # install a runtime rather than fail with a cryptic error.
            raise RuntimeError(self._no_backend_message())
        return runners

    def _build_flag_args(self, trajectory: bool = False,
                         export_histogram: bool = False) -> list:
        """Build the LDES option flags that follow "solve <model> -o <result>".

        Shared by the subprocess runner and the REST client: the server accepts
        the same long-form flags verbatim, so both paths derive their options
        from this single mapping.
        """
        opts = self.options
        cmd = []

        if opts.samples != 200_000:
            cmd.extend(['-s', str(opts.samples)])
        # Always pass the seed. The CLI default is -1 (auto/random seed, see
        # LdesCLI), NOT 23000, so omitting --seed at the nominal seed=23000 makes
        # the native binary pick a random seed and produces non-reproducible
        # results that diverge from the in-JVM wrapper. Emitting it unconditionally
        # keeps native LDES bit-for-bit with the wrapper for any seed (>=0), and
        # still forwards seed=-1 to request a random seed when the user wants one.
        cmd.extend(['--seed', str(opts.seed)])
        if opts.method != 'default':
            cmd.extend(['--method', opts.method])
        if opts.cnvgon:
            cmd.append('--cnvgon')
        if opts.cnvgtol != 0.05:
            cmd.extend(['--cnvgtol', str(opts.cnvgtol)])
        if opts.tranfilter != 'mser5':
            cmd.extend(['--tranfilter', opts.tranfilter])
        if opts.warmupfrac != 0.2:
            cmd.extend(['--warmupfrac', str(opts.warmupfrac)])
        if getattr(opts, 'mserbatch', 5) != 5:
            cmd.extend(['--mserbatch', str(opts.mserbatch)])
        if opts.cimethod != 'obm':
            cmd.extend(['--cimethod', opts.cimethod])
        if getattr(opts, 'obmoverlap', 0.5) != 0.5:
            cmd.extend(['--obmoverlap', str(opts.obmoverlap)])
        if getattr(opts, 'ciminbatch', 10) != 10:
            cmd.extend(['--ciminbatch', str(opts.ciminbatch)])
        if getattr(opts, 'ciminobs', 100) != 100:
            cmd.extend(['--ciminobs', str(opts.ciminobs)])
        if opts.spectral_low_freq_frac != 0.25:
            cmd.extend(['--spectrallowfreqfrac', str(opts.spectral_low_freq_frac)])
        if getattr(opts, 'cnvgbatch', 20) != 20:
            cmd.extend(['--cnvgbatch', str(opts.cnvgbatch)])
        if getattr(opts, 'cnvgchk', 0) != 0:
            cmd.extend(['--cnvgchk', str(opts.cnvgchk)])
        # Discrete-time (slotted) mode. --slotlength implies --slotted on the CLI
        # side, but both are emitted when the length is non-default so the
        # command line states the intent explicitly.
        if getattr(opts, 'slotted', False):
            cmd.append('--slotted')
            slot_len = float(getattr(opts, 'slot_length', 1.0))
            if slot_len != 1.0:
                cmd.extend(['--slotlength', repr(slot_len)])
        if hasattr(opts, 'replications') and opts.replications is not None:
            cmd.extend(['--replications', str(opts.replications)])
        if hasattr(opts, 'numthreads') and opts.numthreads is not None:
            cmd.extend(['--numthreads', str(opts.numthreads)])
        if opts.timespan is not None:
            cmd.extend(['--timespan', f'{opts.timespan[0]},{opts.timespan[1]}'])
        # Cooperative wall-clock budget: the SSJ event loop checks this and stops
        # early with stopping_reason='max_time' (see LdesCLI --maxtime). The outer
        # subprocess timeout is the hard bound if this is not honored in time.
        _tmo = float(getattr(opts, 'timeout', float('inf')))
        if _tmo == _tmo and _tmo not in (float('inf'),) and _tmo > 0:
            cmd.extend(['--maxtime', str(_tmo)])
        if getattr(opts, 'init_sol', None) is not None:
            init_vals = np.asarray(opts.init_sol).flatten()
            cmd.extend(['--initsol', ','.join(repr(float(v)) for v in init_vals)])
        if trajectory:
            cmd.append('--trajectory')
        if export_histogram:
            cmd.append('--export-histogram')

        return cmd

    def _fetch_rest_result(self, rest_url: str, model_path: str, result_path: str,
                           trajectory: bool = False,
                           export_histogram: bool = False) -> None:
        """Solve through an LDES REST server and write its result to result_path.

        The wire format is the model.json the CLI reads and the ldes-result
        document it writes, so a fixed seed and sample count give bit-for-bit
        the same metrics as the subprocess runner.
        """
        import urllib.error

        url = rest_url.rstrip('/')
        if not re.search(r'/api/v\d+/solve$', url):
            url = url + '/api/v1/solve'

        with open(model_path, 'r') as f:
            model_text = f.read()
        payload = {
            'model': {'content': model_text, 'base64': False},
            'flags': self._build_flag_args(trajectory=trajectory,
                                           export_histogram=export_histogram),
        }
        body = json.dumps(payload).encode('utf-8')

        if self.options.verbose not in ('silent',):
            print(f"SolverLDES REST: POST {url} ({' '.join(payload['flags'])})")

        _tmo = float(getattr(self.options, 'timeout', float('inf')))
        req_timeout = (_tmo + 30.0) if math.isfinite(_tmo) and _tmo > 0 else 600

        req = urllib.request.Request(
            url, data=body, method='POST',
            headers={'Content-Type': 'application/json'})
        try:
            with urllib.request.urlopen(req, timeout=req_timeout) as resp:
                payload_out = json.loads(resp.read().decode('utf-8'))
        except urllib.error.HTTPError as e:
            detail = e.read().decode('utf-8', errors='ignore')
            raise RuntimeError(
                f"LDES REST solve failed (HTTP {e.code}): {detail}") from None
        except urllib.error.URLError as e:
            raise RuntimeError(
                f"LDES REST request to {url} failed: {e.reason}") from None

        if payload_out.get('status') != 'ok':
            raise RuntimeError(
                "LDES REST solve failed: %s %s" % (payload_out.get('message', 'unspecified error'),
                                                   payload_out.get('stderr', '')))
        with open(result_path, 'w') as f:
            json.dump(payload_out['result'], f)

    def _parse_result_json(self, path: str) -> LDESResult:
        """Parse LDES result JSON into LDESResult dataclass."""
        with open(path, 'r') as f:
            data = json.load(f)

        # Check for error
        if 'error' in data:
            raise RuntimeError(f"LDES solver error: {data['error']}")

        result = LDESResult()

        # Parse dimensions
        dims = data.get('dimensions', {})
        nstations = dims.get('nstations', 0)
        nclasses = dims.get('nclasses', 0)

        # Update names from result if available
        if 'stationNames' in dims:
            self._station_names = dims['stationNames']
        if 'classNames' in dims:
            self._class_names = dims['classNames']

        # Parse metrics
        metrics = data.get('metrics', {})
        for key in ('QN', 'UN', 'RN', 'TN', 'AN', 'WN', 'CN', 'XN', 'DropRateJoin'):
            val = metrics.get(key)
            if val is not None:
                arr = self._json_array_to_numpy(val)
                setattr(result, key, arr)

        # Parse the exact joint-state residence-time histogram (present when the
        # simulation was run with --export-histogram). Used by get_avg_reward to
        # evaluate arbitrary Markov rewards on the empirical state distribution.
        hist = data.get('stateHistogram')
        if hist is not None:
            sp = hist.get('space')
            tm = hist.get('time')
            if sp is not None:
                result.state_histogram_space = self._json_array_to_numpy(sp)
            if tm is not None:
                result.state_histogram_time = self._json_array_to_numpy(tm)
            tsp = hist.get('trajSpace')
            ttm = hist.get('trajTime')
            if tsp is not None:
                result.state_trajectory_space = self._json_array_to_numpy(tsp)
            if ttm is not None:
                result.state_trajectory_time = self._json_array_to_numpy(ttm)

        # Parse per-cache hit/miss/latency metrics (keyed by cache node name)
        cm = data.get('cacheMetrics', {})
        if cm:
            result.cache_metrics = {}
            for cname, cdata in cm.items():
                entry = {}
                for k in ('hit', 'delayed', 'miss', 'latency', 'hitList', 'itemProb'):
                    v = cdata.get(k)
                    if v is not None:
                        entry[k] = self._json_array_to_numpy(v)
                result.cache_metrics[cname] = entry

        # Parse finite capacity region (FCR) metrics [regions x classes]
        fcr = data.get('fcr', {})
        if fcr:
            for key in ('QNfcr', 'UNfcr', 'RNfcr', 'TNfcr', 'ANfcr', 'WNfcr',
                        'WeightNfcr', 'MemOccNfcr', 'DropRateNfcr'):
                val = fcr.get(key)
                if val is not None:
                    setattr(result, key, self._json_array_to_numpy(val))

        # Parse confidence intervals
        ci = data.get('confidenceIntervals', {})
        for key in ('QNCI', 'UNCI', 'RNCI', 'TNCI', 'ANCI', 'WNCI'):
            val = ci.get(key)
            if val is not None:
                arr = self._json_array_to_numpy(val)
                setattr(result, key, arr)

        # Parse relative precision
        rp = data.get('relativePrecision', {})
        for key in ('QNRelPrec', 'UNRelPrec', 'RNRelPrec', 'TNRelPrec'):
            val = rp.get(key)
            if val is not None:
                arr = self._json_array_to_numpy(val)
                setattr(result, key, arr)

        # Parse metadata
        result.method = data.get('method', 'default')
        result.runtime = data.get('runtime', 0.0)
        result.converged = data.get('converged', False)
        result.stopping_reason = data.get('stoppingReason', '')
        result.convergence_batches = data.get('convergenceBatches', 0)

        # Parse transient trajectory data
        tran = data.get('transient', {})
        if tran:
            t_arr = tran.get('t')
            if t_arr is not None:
                result.t = np.array(t_arr).reshape(-1, 1)

            for key in ('QNt', 'UNt', 'TNt'):
                val = tran.get(key)
                if val is not None:
                    parsed = []
                    for station_arr in val:
                        station_list = []
                        for class_data in station_arr:
                            if class_data is None:
                                station_list.append(None)
                            else:
                                station_list.append(self._json_array_to_numpy(class_data))
                        parsed.append(station_list)
                    setattr(result, key, parsed)

            rts = tran.get('respTimeSamples')
            if rts is not None:
                parsed_rts = []
                for station_arr in rts:
                    station_list = []
                    for class_samples in station_arr:
                        if class_samples is None:
                            station_list.append(None)
                        else:
                            station_list.append([float(v) for v in class_samples])
                    parsed_rts.append(station_list)
                result.respTimeSamples = parsed_rts

        return result

    @staticmethod
    def _json_array_to_numpy(val) -> np.ndarray:
        """Convert a JSON 2D array (with possible null) to numpy array."""
        if isinstance(val, list) and len(val) > 0 and isinstance(val[0], list):
            # 2D array
            rows = []
            for row in val:
                rows.append([float('nan') if v is None else float(v) for v in row])
            return np.array(rows)
        elif isinstance(val, list):
            # 1D array
            return np.array([float('nan') if v is None else float(v) for v in val])
        else:
            return np.array([[float('nan') if val is None else float(val)]])

    def runAnalyzer(self, trajectory: bool = False, export_histogram: bool = False) -> LDESResult:
        """
        Run the LDES simulation via subprocess.

        Args:
            trajectory: If True, request trajectory data (QNt, UNt, TNt, t) from ldes.jar.
            export_histogram: If True, request the exact joint-state residence-time
                histogram (used by get_avg_reward to evaluate arbitrary rewards).

        Returns:
            LDESResult containing performance metrics
        """
        start_time = time.time()

        # Server breakdown (set_breakdown): the joint (queue, server status)
        # chain is expanded only by SolverCTMC, and the LDES JSON wire has no
        # field for a failure or repair process, so the engine would simply
        # simulate an always-up server and return a result indistinguishable
        # from a correct one. 'Breakdown' is absent from getFeatureSet(), but
        # supports() is not consulted here, so reject by name.
        _sn_bd = self.model.getStruct() if hasattr(self.model, 'getStruct') else None
        _hasbd = getattr(_sn_bd, 'hasbreakdown', None) if _sn_bd is not None else None
        if _hasbd is not None and np.any(np.asarray(_hasbd).ravel() == 1):
            _bd_nodes = [str(_sn_bd.nodenames[_i])
                         for _i in np.nonzero(np.asarray(_hasbd).ravel() == 1)[0]]
            raise RuntimeError(
                "Station(s) %s declare server breakdowns, which SolverLDES does not "
                "simulate: the LDES engine has no failure or repair process on its JSON "
                "interface. Use SolverCTMC for models with set_breakdown."
                % ', '.join(_bd_nodes))

        with tempfile.TemporaryDirectory(prefix='line_ldes_') as temp_dir:
            model_path = os.path.join(temp_dir, 'model.json')
            result_path = os.path.join(temp_dir, 'result.json')

            # Save model to JSON
            save_model(self.model, model_path)

            # Remote engine (options.rest_url): an LDES REST server, typically
            # the imperialqore/ldes container, runs the simulation and returns
            # the same ldes-result document. It is written to result_path so the
            # parsing path below is shared with the subprocess runner.
            rest_url = getattr(self.options, 'rest_url', None)
            if rest_url:
                self._fetch_rest_result(rest_url, model_path, result_path,
                                        trajectory=trajectory,
                                        export_histogram=export_histogram)
            else:
                # Build the ordered CLI runner candidates (native, then jar)
                runners = self._build_cli_runners(model_path, result_path,
                                                  trajectory=trajectory,
                                                  export_histogram=export_histogram)
                cmd = runners[0]

                if self.options.verbose not in ('silent',):
                    print(f"SolverLDES command: {' '.join(cmd)}")

                # Wall-clock time budget (options.timeout, seconds). The CLI also gets
                # a cooperative --maxtime flag (see _build_cli_args); this subprocess
                # timeout is the hard outer bound. On expiry the process is killed and
                # an empty result flagged as timed out is returned. Infinite budget
                # falls back to the 600s safety cap.
                import math as _math
                _tmo = float(getattr(self.options, 'timeout', float('inf')))
                # The cooperative --maxtime flag does the real early stop; give the
                # process grace to finish writing results and shut down the JVM before
                # the hard subprocess bound kills it (otherwise no result file is
                # produced). Infinite budget falls back to the 600s safety cap.
                _sub_timeout = (_tmo + 30.0) if _math.isfinite(_tmo) and _tmo > 0 else 600
                # Try each runner in order; the native binary may lack some
                # reflective features (e.g. fork-join MMT serialization), so
                # fall through to the JVM jar on failure.
                last_err = ''
                for ri, cmd in enumerate(runners):
                    if os.path.isfile(result_path):
                        os.remove(result_path)
                    try:
                        proc = subprocess.run(
                            cmd,
                            capture_output=True,
                            cwd=temp_dir,
                            timeout=_sub_timeout,
                        )
                    except subprocess.TimeoutExpired:
                        import warnings as _warnings
                        _warnings.warn("SolverLDES exceeded the wall-clock time budget "
                                       "(timeout=%gs) and was terminated; returning an empty result." % _tmo)
                        empty = LDESResult()
                        empty.timedOut = True
                        empty.stopping_reason = 'max_time'
                        return empty

                    if proc.returncode == 0 and os.path.isfile(result_path):
                        break
                    stdout = proc.stdout.decode('utf-8', errors='ignore')
                    stderr = proc.stderr.decode('utf-8', errors='ignore')
                    last_err = (f"exit code {proc.returncode}\nstdout: {stdout}\n"
                                f"stderr: {stderr}")
                else:
                    raise RuntimeError(
                        f"LDES simulation failed on all {len(runners)} runner(s): {last_err}")

            # Parse result JSON
            self._result = self._parse_result_json(result_path)

        # Restore per-cache hit/miss ratios onto the Cache nodes / nodeparam so
        # cache.get_hit_ratio()/get_miss_ratio() and the node-level hit/miss
        # throughput table work for cache (incl. retrieval) models.
        self._apply_cache_metrics()

        self._result.runtime = time.time() - start_time
        return self._result

    def _apply_cache_metrics(self):
        """Propagate parsed per-cache hit/miss/latency onto the Cache nodes and
        the network struct's nodeparam (read by the node-throughput getters)."""
        result = self._result
        if result is None or not getattr(result, 'cache_metrics', None):
            return
        sn = self._sn
        if sn is None or getattr(sn, 'nodetype', None) is None:
            return
        nodenames = list(sn.nodenames) if getattr(sn, 'nodenames', None) else []

        def _flat(a):
            return np.asarray(a).reshape(-1) if a is not None else None

        for ind in range(sn.nnodes):
            if ind >= len(sn.nodetype) or sn.nodetype[ind] != NodeType.CACHE:
                continue
            name = nodenames[ind] if ind < len(nodenames) else None
            cm = result.cache_metrics.get(name)
            if cm is None:
                continue
            hit = _flat(cm.get('hit'))
            delayed = _flat(cm.get('delayed'))
            miss = _flat(cm.get('miss'))
            lat = _flat(cm.get('latency'))
            if getattr(sn, 'nodeparam', None) is not None and ind in sn.nodeparam \
                    and sn.nodeparam[ind] is not None:
                if hit is not None:
                    sn.nodeparam[ind].actualhitprob = hit
                if delayed is not None:
                    sn.nodeparam[ind].actualdelayedhitprob = delayed
                if miss is not None:
                    sn.nodeparam[ind].actualmissprob = miss
            hit_list = cm.get('hitList')     # [classes x lists], keep 2D
            item_prob = cm.get('itemProb')   # [items x (lists+1)], keep 2D
            if hasattr(self.model, '_nodes') and ind < len(self.model._nodes):
                cnode = self.model._nodes[ind]
                if hit is not None and hasattr(cnode, 'set_result_hit_prob'):
                    cnode.set_result_hit_prob(hit)
                if delayed is not None and hasattr(cnode, 'set_result_delayed_hit_prob'):
                    cnode.set_result_delayed_hit_prob(delayed)
                if miss is not None and hasattr(cnode, 'set_result_miss_prob'):
                    cnode.set_result_miss_prob(miss)
                if lat is not None and hasattr(cnode, 'set_result_residt'):
                    cnode.set_result_residt(lat)
                if hit_list is not None and hasattr(cnode, 'set_result_hit_prob_list'):
                    cnode.set_result_hit_prob_list(np.atleast_2d(hit_list))
                if item_prob is not None and hasattr(cnode, 'set_result_item_prob'):
                    cnode.set_result_item_prob(np.atleast_2d(item_prob))

    def run_analyzer(self) -> LDESResult:
        """Alias for runAnalyzer (Python convention)."""
        return self.runAnalyzer()

    def getAvgReward(self) -> Tuple[np.ndarray, List[str]]:
        """Steady-state expected Markov reward via LDES simulation.

        The simulator exports the exact joint-state residence-time histogram; the
        reward functions defined via model.set_reward are evaluated here on each
        visited state, so E[r] = sum_s (t_s / sum_t) r(state_s) is correct also for
        nonlinear rewards (e.g. E[n^2]).

        Returns:
            Tuple of (R, names): expected reward values and reward names.
        """
        if hasattr(self.model, 'get_rewards'):
            rewards_dict = self.model.get_rewards()
        elif hasattr(self.model, '_rewards'):
            rewards_dict = self.model._rewards
        else:
            rewards_dict = {}
        if not rewards_dict:
            return np.array([]), []

        # Run the simulation requesting the state-residence histogram.
        self._result = None
        result = self.runAnalyzer(export_histogram=True)

        names = list(rewards_dict.keys())
        space = getattr(result, 'state_histogram_space', None)
        time_arr = getattr(result, 'state_histogram_time', None)
        if space is None or time_arr is None or np.sum(time_arr) <= 0:
            return np.array([0.0] * len(names)), names

        space = np.atleast_2d(np.asarray(space, dtype=float))
        time_arr = np.asarray(time_arr, dtype=float).flatten()
        w = time_arr / np.sum(time_arr)

        from ....lang.reward_state import RewardState
        import inspect
        sn = self._sn

        nodes_to_station = {}
        if hasattr(self.model, 'get_nodes') and hasattr(sn, 'nodeToStation'):
            for node in self.model.get_nodes():
                node_idx = node.get_index() if hasattr(node, 'get_index') else getattr(node, 'index', None)
                if node_idx is None:
                    continue
                node_idx0 = node_idx - 1
                if node_idx0 < len(sn.nodeToStation):
                    station_idx = int(sn.nodeToStation[node_idx0])
                    if station_idx >= 0:
                        nodes_to_station[node_idx] = station_idx + 1
        classes_to_idx = {}
        if hasattr(self.model, 'get_classes'):
            for i, jobclass in enumerate(self.model.get_classes()):
                class_idx = jobclass.get_index() if hasattr(jobclass, 'get_index') else getattr(jobclass, 'index', i + 1)
                classes_to_idx[class_idx] = i + 1

        R = []
        for name in names:
            reward_fn = rewards_dict[name]
            try:
                n_params = len(inspect.signature(reward_fn).parameters)
            except (ValueError, TypeError):
                n_params = 1
            acc = 0.0
            for s in range(space.shape[0]):
                reward_state = RewardState(space[s, :], sn, nodes_to_station, classes_to_idx)
                try:
                    val = reward_fn(reward_state, sn) if n_params >= 2 else reward_fn(reward_state)
                    acc += w[s] * float(val)
                except Exception:
                    pass
            R.append(acc)
        return np.array(R), names

    get_avg_reward = getAvgReward

    def _reward_state_maps(self):
        """Build (nodes_to_station, classes_to_idx) 1-based maps for RewardState."""
        sn = self._sn
        nodes_to_station = {}
        if hasattr(self.model, 'get_nodes') and hasattr(sn, 'nodeToStation'):
            for node in self.model.get_nodes():
                node_idx = node.get_index() if hasattr(node, 'get_index') else getattr(node, 'index', None)
                if node_idx is None:
                    continue
                node_idx0 = node_idx - 1
                if node_idx0 < len(sn.nodeToStation):
                    station_idx = int(sn.nodeToStation[node_idx0])
                    if station_idx >= 0:
                        nodes_to_station[node_idx] = station_idx + 1
        classes_to_idx = {}
        if hasattr(self.model, 'get_classes'):
            for i, jobclass in enumerate(self.model.get_classes()):
                class_idx = jobclass.get_index() if hasattr(jobclass, 'get_index') else getattr(jobclass, 'index', i + 1)
                classes_to_idx[class_idx] = i + 1
        return nodes_to_station, classes_to_idx

    def getTranReward(self, reward_name: str = None):
        """Transient reward r(X(t)) via LDES simulation.

        Returns the reward trajectory along the simulated path for each reward
        function defined via model.set_reward, evaluated on the exact integer joint
        state (correct also for nonlinear rewards).

        Returns:
            Tuple (Rt, t, names): Rt is a dict name -> numpy reward series, t is the
            time vector, names is the list of reward names.
        """
        if hasattr(self.model, 'get_rewards'):
            rewards_dict = self.model.get_rewards()
        elif hasattr(self.model, '_rewards'):
            rewards_dict = self.model._rewards
        else:
            rewards_dict = {}
        if not rewards_dict:
            return {}, np.array([]), []

        self._result = None
        result = self.runAnalyzer(export_histogram=True)

        names = [reward_name] if reward_name is not None else list(rewards_dict.keys())
        space = getattr(result, 'state_trajectory_space', None)
        t = getattr(result, 'state_trajectory_time', None)
        if space is None or t is None:
            return {n: np.array([]) for n in names}, np.array([]), names

        space = np.atleast_2d(np.asarray(space, dtype=float))
        t = np.asarray(t, dtype=float).flatten()

        from ....lang.reward_state import RewardState
        import inspect
        sn = self._sn
        nodes_to_station, classes_to_idx = self._reward_state_maps()

        Rt = {n: np.zeros(space.shape[0]) for n in names}
        n_params = {}
        for n in names:
            try:
                n_params[n] = len(inspect.signature(rewards_dict[n]).parameters)
            except (ValueError, TypeError):
                n_params[n] = 1
        for s in range(space.shape[0]):
            reward_state = RewardState(space[s, :], sn, nodes_to_station, classes_to_idx)
            for n in names:
                fn = rewards_dict[n]
                try:
                    Rt[n][s] = float(fn(reward_state, sn) if n_params[n] >= 2 else fn(reward_state))
                except Exception:
                    Rt[n][s] = 0.0
        return Rt, t, names

    get_tran_reward = getTranReward

    def getAvgQLen(self) -> np.ndarray:
        """Average queue lengths [stations x classes], running the simulation
        if needed (JAR/MATLAB SolverLDES API parity)."""
        if self._result is None:
            self.runAnalyzer()
        M = len(self._station_names)
        K = len(self._class_names)
        return self._result.QN.copy() if self._result.QN is not None else np.zeros((M, K))

    def getAvg(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Average station metrics (QN, UN, RN, TN, AN, WN) as station x class
        matrices, applying the same visit-based masking used by getAvgTable.
        """
        return self._computeAvgMetrics()

    get_avg = getAvg

    def _computeAvgMetrics(self):
        """Assemble the station x class average matrices from the last
        simulation result (running it if needed)."""
        if self._result is None:
            self.runAnalyzer()

        result = self._result
        M = len(self._station_names)
        K = len(self._class_names)
        sn = self._sn

        # Make copies of station-level metrics (to avoid modifying originals)
        QN = result.QN.copy() if result.QN is not None else np.zeros((M, K))
        UN = result.UN.copy() if result.UN is not None else np.zeros((M, K))
        RN = result.RN.copy() if result.RN is not None else np.zeros((M, K))
        TN = result.TN.copy() if result.TN is not None else np.zeros((M, K))
        AN = result.AN.copy() if result.AN is not None else np.zeros((M, K))

        # Zero out metrics for classes that don't visit stations based on visit ratios
        hasForkJoin = False
        hasSPN = False
        hasCache = False
        if hasattr(sn, 'nodetype') and sn.nodetype is not None:
            hasForkJoin = np.any(sn.nodetype == NodeType.FORK) and np.any(sn.nodetype == NodeType.JOIN)
            hasSPN = np.any(sn.nodetype == NodeType.PLACE) or np.any(sn.nodetype == NodeType.TRANSITION)
            # A Cache node switches classes by item state, which the
            # routing-based visit equations cannot represent: visits
            # downstream of the cache solve to garbage, so trust the
            # simulation there, like fork-join.
            hasCache = np.any(sn.nodetype == NodeType.CACHE)

        if sn is not None and hasattr(sn, 'nchains') and sn.nchains > 0 and not hasSPN:
            if hasattr(sn, 'chains') and sn.chains is not None and hasattr(sn, 'visits') and sn.visits:
                chains_arr = np.asarray(sn.chains)
                # A chain fed by spawn-on-completion (LQN phase-2) gets its jobs
                # by injection at a station rather than by routing, so the
                # routing-based visit equations solve to zero for its classes
                # even though they are served. Trust the simulation there, like
                # fork-join.
                spawn_fed = np.zeros(chains_arr.shape[0], dtype=bool)
                try:
                    classes = self.model.get_classes()
                    idx_of = {id(jc): j for j, jc in enumerate(classes)}
                    for jc in classes:
                        tgt = jc.get_spawn_class() if hasattr(jc, 'get_spawn_class') else None
                        if tgt is not None:
                            tj = idx_of.get(id(tgt))
                            if tj is not None and tj < chains_arr.shape[1]:
                                spawn_fed |= chains_arr[:, tj] > 0
                except Exception:
                    pass
                for k in range(K):
                    chains_with_class = np.where(chains_arr[:, k] > 0)[0] if k < chains_arr.shape[1] else []
                    if len(chains_with_class) > 0:
                        c = chains_with_class[0]
                        if c in sn.visits and sn.visits[c] is not None:
                            visits_c = np.asarray(sn.visits[c])
                            for i in range(M):
                                if i < visits_c.shape[0] and k < visits_c.shape[1]:
                                    if visits_c[i, k] == 0:
                                        if (hasForkJoin or hasCache or spawn_fed[c]) and (
                                                QN[i, k] > GlobalConstants.FineTol or
                                                UN[i, k] > GlobalConstants.FineTol or
                                                TN[i, k] > GlobalConstants.FineTol):
                                            continue
                                        QN[i, k] = 0
                                        UN[i, k] = 0
                                        RN[i, k] = 0
                                        TN[i, k] = 0
                                        AN[i, k] = 0

        # Compute residence times from response times using visit ratios
        WN = None
        if RN is not None and sn is not None and hasattr(sn, 'visits') and sn.visits:
            try:
                WN = sn_get_residt_from_respt(sn, RN, None)
            except Exception:
                WN = RN.copy() if RN is not None else None

        return QN, UN, RN, TN, AN, WN

    def getAvgTable(self) -> pd.DataFrame:
        """
        Get average performance metrics as a DataFrame.

        Returns:
            DataFrame with columns: Station, JobClass, QLen, Util, RespT, ResidT, ArvR, Tput
        """
        QN, UN, RN, TN, AN, WN = self._computeAvgMetrics()
        M = len(self._station_names)
        K = len(self._class_names)

        # Identify source stations
        source_stations = set()
        nodetype = self._sn.nodetype if hasattr(self._sn, 'nodetype') else None
        if nodetype is not None:
            stationToNode = self._sn.stationToNode
            if stationToNode is not None:
                stationToNode = np.asarray(stationToNode).flatten()
                nodetype = np.asarray(nodetype).flatten()
                for i in range(M):
                    if i < len(stationToNode):
                        node_idx = int(stationToNode[i])
                        if node_idx < len(nodetype):
                            if int(nodetype[node_idx]) == 0:  # SOURCE = 0
                                source_stations.add(i)

        rows = []
        for i in range(M):
            for r in range(K):
                qlen = QN[i, r] if i < QN.shape[0] and r < QN.shape[1] else 0
                util = UN[i, r] if i < UN.shape[0] and r < UN.shape[1] else 0
                respt = RN[i, r] if i < RN.shape[0] and r < RN.shape[1] else 0
                tput = TN[i, r] if i < TN.shape[0] and r < TN.shape[1] else 0
                arvr = AN[i, r] if i < AN.shape[0] and r < AN.shape[1] else 0

                is_source = i in source_stations

                if is_source:
                    if hasattr(self._sn, 'rates') and self._sn.rates is not None:
                        rates = np.asarray(self._sn.rates)
                        stationToNode = np.asarray(self._sn.stationToNode).flatten()
                        node_idx = int(stationToNode[i])
                        if node_idx < rates.shape[0] and r < rates.shape[1]:
                            tput = rates[node_idx, r]
                    arvr = 0.0

                if abs(qlen) < 1e-12 and abs(util) < 1e-12 and abs(tput) < 1e-12:
                    continue

                residt = respt
                if WN is not None and i < WN.shape[0] and r < WN.shape[1]:
                    residt_val = WN[i, r]
                    if not np.isnan(residt_val) and residt_val >= 0:
                        residt = residt_val

                rows.append({
                    'Station': self._station_names[i],
                    'JobClass': self._class_names[r],
                    'QLen': qlen,
                    'Util': util,
                    'RespT': respt,
                    'ResidT': residt,
                    'ArvR': arvr,
                    'Tput': tput,
                })

        df = pd.DataFrame(rows)

        if not self._table_silent and len(df) > 0:
            print(df.to_string(index=False))

        return df

    def getAvgNode(self):
        """Average metrics per node (incl. non-station nodes such as Cache).

        Cache hit/miss class throughputs come from the actual hit/miss
        probabilities restored from the simulation result (see
        _apply_cache_metrics); the JAR LDES simulator already folds delayed
        hits into the hit ratio, so hit + miss sum to one at the cache.
        """
        from ....api.sn.getters import (sn_get_node_arvr_from_tput,
                                       sn_get_node_tput_from_tput)

        if self._result is None:
            self.runAnalyzer()

        sn = self._sn
        I = sn.nnodes
        M = sn.nstations
        R = sn.nclasses
        result = self._result

        QN = result.QN.copy() if result.QN is not None else np.zeros((M, R))
        UN = result.UN.copy() if result.UN is not None else np.zeros((M, R))
        RN = result.RN.copy() if result.RN is not None else np.zeros((M, R))
        TN = result.TN.copy() if result.TN is not None else np.zeros((M, R))
        AN = None
        if result.AN is not None and np.any(result.AN > 0):
            AN = result.AN

        WN = sn_get_residt_from_respt(sn, RN, None)
        TH = np.zeros_like(TN)
        TH[TN > 0] = 1.0

        ANn = sn_get_node_arvr_from_tput(sn, TN, TH, AN)
        TNn = sn_get_node_tput_from_tput(sn, TN, TH, ANn)

        QNn = np.zeros((I, R))
        UNn = np.zeros((I, R))
        RNn = np.zeros((I, R))
        WNn = np.zeros((I, R))
        for ist in range(M):
            ind = sn.stationToNode[ist]
            if 0 <= ind < I:
                QNn[ind, :] = QN[ist, :]
                UNn[ind, :] = UN[ist, :]
                RNn[ind, :] = RN[ist, :]
                WNn[ind, :] = WN[ist, :]

        return QNn, UNn, RNn, WNn, ANn, TNn

    def getAvgNodeTable(self) -> pd.DataFrame:
        """Average metrics by node as a DataFrame (one row per node per class).

        Includes non-station nodes such as Cache, whose hit/miss class
        throughputs are reconstructed from the simulation hit/miss ratios.
        """
        QNn, UNn, RNn, WNn, ANn, TNn = self.getAvgNode()
        sn = self._sn
        nodenames = list(sn.nodenames) if getattr(sn, 'nodenames', None) else []

        from ...cache_table import retrieval_hidden_classes
        hidden = retrieval_hidden_classes(sn)

        rows = []
        for ind in range(sn.nnodes):
            node_name = nodenames[ind] if ind < len(nodenames) else f'Node{ind}'
            for r in range(sn.nclasses):
                if r in hidden:
                    continue  # auxiliary retrieval class - omit from node table
                class_name = self._class_names[r] if r < len(self._class_names) else f'Class{r}'
                if abs(QNn[ind, r]) < 1e-10 and abs(UNn[ind, r]) < 1e-10 and \
                   abs(RNn[ind, r]) < 1e-10 and abs(ANn[ind, r]) < 1e-10 and abs(TNn[ind, r]) < 1e-10:
                    continue
                rows.append({
                    'Node': node_name,
                    'JobClass': class_name,
                    'QLen': QNn[ind, r],
                    'Util': UNn[ind, r],
                    'RespT': RNn[ind, r],
                    'ResidT': WNn[ind, r],
                    'ArvR': ANn[ind, r],
                    'Tput': TNn[ind, r],
                })

        df = pd.DataFrame(rows)
        if not self._table_silent and len(df) > 0:
            print(df.to_string(index=False))
        return df

    def getAvgRegionTable(self) -> pd.DataFrame:
        """Per-region finite capacity region (FCR) metrics as a DataFrame.

        One row per region per class. In addition to the queue-length, response
        time and throughput of the region, the table reports the two FCR-specific
        metrics: Weight (time-average weighted occupation, i.e. FCR Total Weight)
        and MemOcc (time-average memory occupation). Row sums over the classes of
        a region give the region aggregates.
        """
        result = self._result
        sn = self._sn
        nregions = int(getattr(sn, 'nregions', 0) or 0)
        rows = []
        if result is not None and nregions > 0 and getattr(result, 'QNfcr', None) is not None:
            QN = np.asarray(result.QNfcr, dtype=float)
            RN = np.asarray(result.RNfcr, dtype=float) if result.RNfcr is not None else None
            WN = np.asarray(result.WNfcr, dtype=float) if result.WNfcr is not None else None
            TN = np.asarray(result.TNfcr, dtype=float) if result.TNfcr is not None else None
            AN = np.asarray(result.ANfcr, dtype=float) if result.ANfcr is not None else None
            WG = np.asarray(result.WeightNfcr, dtype=float) if result.WeightNfcr is not None else None
            MO = np.asarray(result.MemOccNfcr, dtype=float) if result.MemOccNfcr is not None else None
            nclasses = QN.shape[1]
            for f in range(nregions):
                region_name = f'FCRegion{f + 1}'
                for r in range(nclasses):
                    class_name = self._class_names[r] if r < len(self._class_names) else f'Class{r}'
                    rows.append({
                        'Region': region_name,
                        'JobClass': class_name,
                        'QLen': QN[f, r],
                        'RespT': RN[f, r] if RN is not None else float('nan'),
                        'ResidT': WN[f, r] if WN is not None else float('nan'),
                        'ArvR': AN[f, r] if AN is not None else float('nan'),
                        'Tput': TN[f, r] if TN is not None else float('nan'),
                        'Weight': WG[f, r] if WG is not None else float('nan'),
                        'MemOcc': MO[f, r] if MO is not None else float('nan'),
                    })
        df = pd.DataFrame(rows)
        if not self._table_silent and len(df) > 0:
            print(df.to_string(index=False))
        return df

    get_avg_region_table = getAvgRegionTable
    avg_region_table = getAvgRegionTable

    def getAvgCacheTable(self) -> pd.DataFrame:
        """Detailed per-class cache performance metrics (see cache_table)."""
        from ...cache_table import build_cache_avg_table
        return build_cache_avg_table(self)

    get_avg_cache_table = getAvgCacheTable
    avg_cache_table = getAvgCacheTable

    def getAvgItemTable(self) -> pd.DataFrame:
        """Item-level cache occupancy table (see cache_table)."""
        from ...cache_table import build_item_avg_table
        return build_item_avg_table(self)

    get_avg_item_table = getAvgItemTable
    avg_item_table = getAvgItemTable

    def getAvgChainTable(self) -> pd.DataFrame:
        """
        Get average performance metrics aggregated by chain.

        Returns:
            DataFrame with columns: Chain, QLen, Util, RespT, Tput
        """
        if self._result is None:
            self.runAnalyzer()

        result = self._result
        nchains = self._sn.nchains if hasattr(self._sn, 'nchains') else self._sn.nclasses
        inchain = self._sn.inchain if hasattr(self._sn, 'inchain') else None

        rows = []
        for c in range(nchains):
            chain_name = f'Chain{c+1}'

            if inchain is not None and c in inchain:
                chain_classes = np.asarray(inchain[c]).flatten().astype(int)
            else:
                chain_classes = [c]

            total_qlen = 0.0
            total_util = 0.0
            total_respt = 0.0
            total_tput = 0.0

            M = self._sn.nstations
            for i in range(M):
                for k in chain_classes:
                    if k < result.QN.shape[1] if result.QN is not None else 0:
                        total_qlen += result.QN[i, k] if result.QN is not None and not np.isnan(result.QN[i, k]) else 0.0
                        total_util += result.UN[i, k] if result.UN is not None and not np.isnan(result.UN[i, k]) else 0.0
                        total_respt += result.RN[i, k] if result.RN is not None and not np.isnan(result.RN[i, k]) else 0.0
                        total_tput = max(total_tput, result.TN[i, k] if result.TN is not None and not np.isnan(result.TN[i, k]) else 0.0)

            rows.append({
                'Chain': chain_name,
                'QLen': total_qlen,
                'Util': total_util,
                'RespT': total_respt,
                'Tput': total_tput,
            })

        return pd.DataFrame(rows)

    def getAvgSysTable(self) -> pd.DataFrame:
        """
        Get system-level average performance metrics.

        Returns:
            DataFrame with columns: Chain, SysRespT, SysTput
        """
        if self._result is None:
            self.runAnalyzer()

        chain_table = self.getAvgChainTable()
        CN = []
        XN = []
        for _, row in chain_table.iterrows():
            CN.append(row['RespT'])
            XN.append(row['Tput'])
        return self._make_sys_table(CN, XN)

    # =========================================================================
    # Transient analysis
    # =========================================================================

    def getTranAvg(self) -> LDESResult:
        """
        Run transient analysis and return result with time-series metrics.

        If options.timespan is not set, defaults to [0, 30/min_rate].
        The result will contain QNt, UNt, TNt, t fields with trajectory data.

        Returns:
            LDESResult with transient time-series data populated.
        """
        if self.options.timespan is None:
            sn = self._sn
            min_rate = float('inf')
            if hasattr(sn, 'rates') and sn.rates is not None:
                rates = np.asarray(sn.rates).flatten()
                positive_rates = rates[rates > 0]
                if len(positive_rates) > 0:
                    min_rate = float(np.min(positive_rates))
            if min_rate == float('inf'):
                min_rate = 1.0
            self.options.timespan = [0.0, 30.0 / min_rate]

        return self.runAnalyzer(trajectory=True)

    # =========================================================================
    # Sampling methods
    # =========================================================================

    def _run_transient(self, num_events: int) -> LDESResult:
        """Run transient simulation with the given number of events as time horizon.

        The horizon is passed explicitly via timespan; options.samples (the
        steady-state event budget) is left untouched, since the engine ignores
        it in transient mode.
        """
        original_timespan = self.options.timespan
        horizon = float(num_events) if num_events > 0 else float(self.options.samples)
        self.options.timespan = [0.0, horizon]

        try:
            result = self.runAnalyzer(trajectory=True)
        finally:
            self.options.timespan = original_timespan

        return result

    def sample(self, node, num_events: int = 0) -> Optional[dict]:
        """
        Generate a sample path (state trajectory) for a specific node.

        Args:
            node: The stateful node to sample (or node index).
            num_events: Number of events for the simulation time horizon.
                        If 0, uses options.samples.

        Returns:
            Dict with keys 't' (time vector), 'state' (numpy array [timepoints x classes]),
            or None if simulation produced no trajectory data.
        """
        result = self._run_transient(num_events)

        if result.t is None or result.QNt is None:
            return None

        sn = self._sn
        node_idx = node if isinstance(node, int) else getattr(node, '_node_idx', 0)

        isf = int(np.asarray(sn.nodeToStateful).flatten()[node_idx])
        num_time_points = result.t.shape[0]
        num_classes = int(sn.nclasses)
        state = np.zeros((num_time_points, num_classes))

        if isf < len(result.QNt):
            for k in range(min(num_classes, len(result.QNt[isf]))):
                class_data = result.QNt[isf][k]
                if class_data is not None:
                    n = min(num_time_points, class_data.shape[0])
                    state[:n, k] = class_data[:n, 0]

        return {'t': result.t.flatten(), 'state': state}

    def sampleAggr(self, node, num_events: int = 0) -> Optional[dict]:
        """
        Generate an aggregated sample path for a specific node.

        For LDES, sample() already returns per-class queue lengths, so this
        is the same as sample() but marked as aggregated.

        Returns:
            Dict with keys 't', 'state', 'aggregated'=True.
        """
        result = self.sample(node, num_events)
        if result is not None:
            result['aggregated'] = True
        return result

    def sampleSys(self, num_events: int = 0) -> Optional[dict]:
        """
        Generate a system-wide sample path for all stateful nodes.

        Returns:
            Dict with keys 't' (time vector), 'states' (list of numpy arrays,
            one per stateful node), or None if no trajectory data.
        """
        result = self._run_transient(num_events)

        if result.t is None or result.QNt is None:
            return None

        sn = self._sn
        num_time_points = result.t.shape[0]
        num_classes = int(sn.nclasses)

        states = []
        for isf in range(int(sn.nstateful)):
            node_state = np.zeros((num_time_points, num_classes))
            if isf < len(result.QNt):
                for k in range(min(num_classes, len(result.QNt[isf]))):
                    class_data = result.QNt[isf][k]
                    if class_data is not None:
                        n = min(num_time_points, class_data.shape[0])
                        node_state[:n, k] = class_data[:n, 0]
            states.append(node_state)

        return {'t': result.t.flatten(), 'states': states}

    def sampleSysAggr(self, num_events: int = 0) -> Optional[dict]:
        """
        Generate an aggregated system-wide sample path.

        Returns:
            Dict with keys 't', 'states', 'aggregated'=True.
        """
        result = self.sampleSys(num_events)
        if result is not None:
            result['aggregated'] = True
        return result

    # =========================================================================
    # Probability estimation methods
    # =========================================================================

    def getProb(self, node, state: Optional[np.ndarray] = None) -> float:
        """
        Estimate steady-state probability of a specific state at a node.

        Uses time-weighted fraction from simulation trajectory.

        Args:
            node: The stateful node (or node index).
            state: Target state vector (per-class job counts). If None, uses
                   current model state.

        Returns:
            Estimated probability (0.0 if state not observed).
        """
        sample_result = self.sample(node, 0)
        if sample_result is None:
            return 0.0

        t = sample_result['t']
        state_matrix = sample_result['state']
        num_time_points = len(t)

        if num_time_points < 2:
            return 0.0

        if state is None:
            sn = self._sn
            node_idx = node if isinstance(node, int) else getattr(node, '_node_idx', 0)
            isf = int(np.asarray(sn.nodeToStateful).flatten()[node_idx])
            state = np.asarray(sn.state[isf]).flatten() if hasattr(sn, 'state') and sn.state else np.zeros(state_matrix.shape[1])

        target = np.asarray(state).flatten()
        total_time = t[-1] - t[0]
        if total_time <= 0:
            return 0.0

        time_in_state = 0.0
        for ti in range(num_time_points - 1):
            dt = t[ti + 1] - t[ti]
            if np.allclose(state_matrix[ti, :len(target)], target, atol=1e-10):
                time_in_state += dt

        return time_in_state / total_time

    def getProbAggr(self, node, state_aggr: Optional[np.ndarray] = None) -> float:
        """
        Estimate aggregated state probability at a node.

        For LDES, same as getProb since sample paths are already per-class.

        Returns:
            Estimated probability.
        """
        return self.getProb(node, state_aggr)

    def getProbSys(self) -> float:
        """
        Estimate joint steady-state probability of the current system state.

        Uses system-wide trajectory to compute time-weighted fraction.

        Returns:
            Estimated joint probability.
        """
        sys_result = self.sampleSys(0)
        if sys_result is None:
            return 0.0

        t = sys_result['t']
        states = sys_result['states']
        num_time_points = len(t)

        if num_time_points < 2 or not states:
            return 0.0

        sn = self._sn
        target_states = []
        for isf in range(len(states)):
            if hasattr(sn, 'state') and sn.state and isf in sn.state:
                target_states.append(np.asarray(sn.state[isf]).flatten())
            else:
                target_states.append(np.zeros(states[isf].shape[1]))

        total_time = t[-1] - t[0]
        if total_time <= 0:
            return 0.0

        time_in_state = 0.0
        for ti in range(num_time_points - 1):
            dt = t[ti + 1] - t[ti]
            all_match = True
            for isf in range(len(states)):
                if not np.allclose(states[isf][ti, :len(target_states[isf])],
                                   target_states[isf], atol=1e-10):
                    all_match = False
                    break
            if all_match:
                time_in_state += dt

        return time_in_state / total_time

    def getProbSysAggr(self) -> float:
        """
        Estimate aggregated joint system probability.

        For LDES, same as getProbSys since trajectories are already per-class.

        Returns:
            Estimated joint probability.
        """
        return self.getProbSys()

    # =========================================================================
    # Transient CDF methods
    # =========================================================================

    def getTranCdfRespT(self) -> Optional[dict]:
        """
        Get empirical CDF of response times from simulation samples.

        Returns:
            Dict with keys 'station_names', 'class_names', 'cdfs'
            where cdfs[i][k] is a numpy array (n x 2) with [time, CDF_value]
            columns, or None if no response time samples available.
        """
        if self._result is None or self._result.respTimeSamples is None:
            result = self._run_transient(0)
        else:
            result = self._result

        if result.respTimeSamples is None:
            return None

        cdfs = []
        for i, station_samples in enumerate(result.respTimeSamples):
            station_cdfs = []
            for k, samples in enumerate(station_samples):
                if samples and len(samples) > 0:
                    sorted_samples = np.sort(samples)
                    n = len(sorted_samples)
                    cdf_values = np.arange(1, n + 1) / n
                    cdf_matrix = np.column_stack([sorted_samples, cdf_values])
                    station_cdfs.append(cdf_matrix)
                else:
                    station_cdfs.append(None)
            cdfs.append(station_cdfs)

        return {
            'station_names': self._station_names,
            'class_names': self._class_names,
            'cdfs': cdfs,
        }

    def getTranCdfPassT(self) -> Optional[dict]:
        """
        Get empirical CDF of passage times from simulation samples.

        For LDES, passage times are approximated by response times.

        Returns:
            Same format as getTranCdfRespT().
        """
        return self.getTranCdfRespT()

    # =========================================================================
    # Transient probability methods
    # =========================================================================

    def getTranProb(self, node) -> Optional[dict]:
        """
        Estimate transient state probabilities at a node over time.

        Args:
            node: The stateful node (or node index).

        Returns:
            Dict with 't' (time points), 'states' (unique state vectors),
            'probabilities' (time-windowed probabilities), or None.
        """
        sample_result = self.sample(node, 0)
        if sample_result is None:
            return None

        t = sample_result['t']
        state_matrix = sample_result['state']

        unique_states = np.unique(np.round(state_matrix, 10), axis=0)
        num_windows = min(100, len(t) // 2)
        if num_windows < 1:
            return None

        window_size = len(t) // num_windows
        t_windows = np.zeros(num_windows)
        probs = np.zeros((num_windows, len(unique_states)))

        for w in range(num_windows):
            start = w * window_size
            end = min((w + 1) * window_size, len(t) - 1)
            t_windows[w] = t[start]
            window_total = t[end] - t[start]
            if window_total <= 0:
                continue
            for si, s in enumerate(unique_states):
                time_in = 0.0
                for ti in range(start, end):
                    dt = t[ti + 1] - t[ti]
                    if np.allclose(state_matrix[ti, :len(s)], s, atol=1e-10):
                        time_in += dt
                probs[w, si] = time_in / window_total

        return {'t': t_windows, 'states': unique_states, 'probabilities': probs}

    def getTranProbAggr(self, node) -> Optional[dict]:
        """Transient aggregated state probabilities. For LDES, same as getTranProb."""
        return self.getTranProb(node)

    def getTranProbSys(self) -> Optional[dict]:
        """
        Estimate transient joint system probabilities over time.

        Returns:
            Dict with 't', 'probabilities' for system state, or None.
        """
        sys_result = self.sampleSys(0)
        if sys_result is None:
            return None

        t = sys_result['t']
        states = sys_result['states']
        if len(t) < 2 or not states:
            return None

        # Build combined system state at each time point
        combined = np.hstack(states)
        unique_sys_states = np.unique(np.round(combined, 10), axis=0)

        num_windows = min(100, len(t) // 2)
        if num_windows < 1:
            return None

        window_size = len(t) // num_windows
        t_windows = np.zeros(num_windows)
        probs = np.zeros((num_windows, len(unique_sys_states)))

        for w in range(num_windows):
            start = w * window_size
            end = min((w + 1) * window_size, len(t) - 1)
            t_windows[w] = t[start]
            window_total = t[end] - t[start]
            if window_total <= 0:
                continue
            for si, s in enumerate(unique_sys_states):
                time_in = 0.0
                for ti in range(start, end):
                    dt = t[ti + 1] - t[ti]
                    if np.allclose(combined[ti], s, atol=1e-10):
                        time_in += dt
                probs[w, si] = time_in / window_total

        return {'t': t_windows, 'states': unique_sys_states, 'probabilities': probs}

    def getTranProbSysAggr(self) -> Optional[dict]:
        """Transient aggregated system probabilities. For LDES, same as getTranProbSys."""
        return self.getTranProbSys()

    # Method aliases (consistent with other solvers)
    avgT = getAvgTable
    aT = getAvgTable
    get_avg_node_table = getAvgNodeTable
    avg_node_table = getAvgNodeTable
    get_avg_node = getAvgNode
    avg_node = getAvgNode
    get_avg_chain_table = getAvgChainTable
    avg_chain_table = getAvgChainTable
    aCT = getAvgChainTable
    chainAvgT = getAvgChainTable
    get_avg_sys_table = getAvgSysTable
    avg_sys_table = getAvgSysTable

    @property
    def result(self) -> Optional[LDESResult]:
        """Get the LDES result (after runAnalyzer is called)."""
        return self._result

    def getName(self) -> str:
        """Get solver name."""
        return "LDES"

    def get_name(self) -> str:
        """Get solver name (Python convention)."""
        return self.getName()

    def isStochasticMethod(self, method):
        """LDES is a discrete-event simulator; all methods are stochastic."""
        return True

    is_stochastic_method = isStochasticMethod

    @staticmethod
    def defaultOptions() -> LDESOptions:
        """Get default LDES solver options."""
        return LDESOptions()

    @staticmethod
    def default_options() -> LDESOptions:
        """Get default options (Python convention)."""
        return LDESOptions()

    @staticmethod
    def getFeatureSet() -> set:
        """Get set of features supported by the LDES solver.

        Returns the canonical feature names (mirrors MATLAB
        SolverLDES.getFeatureSet and the JAR SolverLDES; every row calls the
        same JAR simulation engine, so the sets must stay identical).
        """
        return {
            'Sink', 'Source',
            'Queue', 'Delay',
            'Fork', 'Join', 'Forker', 'Joiner',
            'Place', 'Transition',
            'Linkage', 'Enabling', 'Inhibiting', 'Timing', 'Firing', 'Storage',
            'Logger', 'LogTunnel',
            'Cache', 'CacheClassSwitcher',
            'ReplacementStrategy_LRU', 'ReplacementStrategy_FIFO', 'ReplacementStrategy_RR', 'ReplacementStrategy_SFIFO',
            'ReplacementStrategy_HLRU', 'ReplacementStrategy_CLIMB', 'ReplacementStrategy_QLRU',
            'Buffer', 'Region',
            'Exp', 'Erlang', 'HyperExp', 'PH', 'APH', 'Coxian', 'Cox2', 'MAP', 'DMAP', 'MMAP', 'BMAP', 'MMPP2', 'ME', 'RAP', 'Immediate', 'Disabled', 'Replayer', 'Trace',  # Trace is an alias of Replayer
            'Det', 'Uniform', 'Gamma', 'Pareto', 'Weibull', 'Lognormal',
            'Geometric',  # Lattice-valued interarrival/service time on {1,2,...} (Geo/Geo/1 and slotted models)
            'Bernoulli', 'Binomial', 'Poisson',  # Counting distributions; zero atom becomes an immediate interval (continuous mode only)
            'NHPP',
            'Server', 'JobSink', 'RandomSource',
            'InfiniteServer', 'SharedServer', 'ServiceTunnel', 'DelayStation',  # internal station-section markers
            'SchedStrategy_FCFS', 'SchedStrategy_INF',
            'SchedStrategy_HOL', 'SchedStrategy_FCFSPRIO',
            'SchedStrategy_PS', 'SchedStrategy_DPS', 'SchedStrategy_GPS',
            'SchedStrategy_PSPRIO', 'SchedStrategy_DPSPRIO', 'SchedStrategy_GPSPRIO',
            'SchedStrategy_LCFS', 'SchedStrategy_LCFSPR', 'SchedStrategy_LCFSPI',
            'SchedStrategy_LCFSPRIO', 'SchedStrategy_LCFSPRPRIO', 'SchedStrategy_LCFSPIPRIO',
            'SchedStrategy_SIRO',
            'SchedStrategy_SJF', 'SchedStrategy_LJF',
            'SchedStrategy_LEPT', 'SchedStrategy_SEPT',
            'SchedStrategy_SRPT', 'SchedStrategy_SRPTPRIO',
            'SchedStrategy_PSJF', 'SchedStrategy_FB', 'SchedStrategy_LRPT',
            'SchedStrategy_FSP',
            'SchedStrategy_EDD', 'SchedStrategy_EDF', 'SchedStrategy_SETF',
            'SchedStrategy_FCFSPR', 'SchedStrategy_FCFSPI',
            'SchedStrategy_FCFSPRPRIO', 'SchedStrategy_FCFSPIPRIO',
            'SchedStrategy_LPS', 'SchedStrategy_EXT', 'SchedStrategy_POLLING',
            'SchedStrategy_PAS',
            'SchedStrategy_OI',
            'Router', 'Dispatcher',
            'ClassSwitch', 'StatelessClassSwitcher',
            'RoutingStrategy_PROB', 'RoutingStrategy_RAND',
            'RoutingStrategy_RROBIN', 'RoutingStrategy_WRROBIN',
            'RoutingStrategy_JSQ', 'RoutingStrategy_SQ',
            'OpenClass', 'ClosedClass', 'SelfLoopingClass',
            'OpenSignal', 'ClosedSignal',
            'SignalType_NEGATIVE', 'SignalType_REPLY', 'SignalType_CATASTROPHE',
            'SignalBatchRemoval', 'SignalRemovalPolicy',
            'LoadDependence', 'ClassDependence', 'JointDependence', 'Balking', 'Reneging', 'Retrial',
            # Engine simulates the SETUP/DELAYOFF server states.
            'SetupDelayOff',
        }

    @staticmethod
    def supports(model) -> bool:
        """Check if model is supported.

        Mirrors MATLAB SolverLDES.supports. No supports() existed anywhere in
        the MRO, so calling it raised AttributeError and the solver had no gate.
        """
        from ...base import supports_via_featureset
        return supports_via_featureset(SolverLDES, model)
