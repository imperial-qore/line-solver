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

from .ldes_options import LDESOptions, LDESResult, LNLDESResult
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
    Native Python LDES solver.

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
            # The launcher is java.exe on Windows, where a bare "java" is also
            # not executable by os.access and the JAVA_HOME branch would be
            # skipped even with a perfectly good JDK installed.
            exe_name = 'java.exe' if platform.system() == 'Windows' else 'java'
            cand = os.path.join(java_home, 'bin', exe_name)
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
                         export_histogram: bool = False,
                         respt_samples: bool = False) -> list:
        """Build CLI command arguments for LDES solver.

        Prefers native binary (common/ldes) for faster startup;
        falls back to java -jar ldes.jar if native binary not found.
        """
        return self._build_cli_runners(model_path, result_path,
                                       respt_samples=respt_samples,
                                       trajectory=trajectory,
                                       export_histogram=export_histogram)[0]

    def _build_cli_runners(self, model_path: str, result_path: str,
                           trajectory: bool = False,
                           export_histogram: bool = False,
                           respt_samples: bool = False) -> list:
        """Ordered CLI command candidates: native binary first, then ldes.jar.

        The native C++ binary may lack features the JVM jar
        has (e.g. the fork-join MMT transformation deep-copies the model via
        Java serialization, unsupported in the AOT image), so the caller runs
        each candidate in order and falls through on failure, mirroring the
        MATLAB solveCli runner list.
        """
        opts = self.options

        native_path = self._get_ldes_native_path()
        # --respt-samples rides on a result block the C++ engine does not fill
        # and which it refuses by name, so a warm-start of THAT kind still goes
        # to the jar; a runner list that tried the native binary first would
        # fall through only after paying for a failed process. Mirrors
        # solveCli.m and ldes_flags in the C++ client.
        #
        # --initsol and --busyperiod were on this list too, because the retired
        # GraalVM image predated them. The C++ engine honours both (warm-start
        # placement 2026-08-02, busy periods once the CLI stopped dropping the
        # block), verified against ldes.jar to the last digit, so they no longer
        # force the jar.
        if respt_samples:
            native_path = None
        # A LayeredNetwork is served by the JVM engine only: the C++ LDES CLI
        # reads the Network document alone, so offering the native binary first
        # would buy a failed process before the fall-through.
        if self._is_layered():
            native_path = None
        # Storage cost caps postdate the prebuilt binary too; running them on
        # it would silently simulate the UNCAPPED cache.
        for nd in getattr(self.model, '_nodes', []) or []:
            if getattr(nd, '_cost_cap', None) is not None:
                native_path = None
                break
        flags = self._build_flag_args(respt_samples=respt_samples,
                                      trajectory=trajectory,
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
                         export_histogram: bool = False,
                         respt_samples: bool = False) -> list:
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
        # METHOD='PARALLEL' IS A REPLICATION COUNT, not a second engine. The
        # engine's parallel analyzer is selected by --replications > 1 and by
        # nothing else, so the method name has to resolve to one: it takes the
        # count the caller supplied, and 8 when there is none -- the default the
        # SSA parallel analyzer uses for its replica count.
        _reps = getattr(opts, 'replications', None)
        if str(getattr(opts, 'method', '') or '').lower() == 'parallel' \
                and not (isinstance(_reps, (int, float)) and _reps > 1):
            _reps = 8
        if _reps is not None:
            cmd.extend(['--replications', str(int(_reps))])
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
        if respt_samples:
            cmd.append('--respt-samples')
        if trajectory:
            cmd.append('--trajectory')
        if export_histogram:
            cmd.append('--export-histogram')
        if getattr(opts, 'busy_period_orders', 0) > 0:
            cmd.extend(['--busyperiod', str(int(opts.busy_period_orders))])
            for subnet in getattr(opts, 'busy_period_subnets', []) or []:
                # station indexes cross the wire zero-based
                cmd.extend(['--busyperiod-subnet',
                            ','.join(str(int(v)) for v in subnet)])

        return cmd

    def _fetch_rest_result(self, rest_url: str, model_path: str, result_path: str,
                           trajectory: bool = False,
                           export_histogram: bool = False,
                           respt_samples: bool = False) -> None:
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
            'flags': self._build_flag_args(respt_samples=respt_samples,
                                           trajectory=trajectory,
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

    def _is_layered(self) -> bool:
        """True when this solver holds a LayeredNetwork rather than a Network."""
        from ....layered import LayeredNetwork
        return isinstance(self.model, LayeredNetwork)

    def _parse_ln_result_json(self, data: dict) -> 'LNLDESResult':
        """Parse the LAYERED `ldes-result` document.

        Its metrics are vectors over the LQN element index space (hosts, tasks,
        entries, activities, with the shifts carried in `dimensions`), not the
        (station, class) matrices of the Network document, so it has its own
        container and its own getters. Written by
        `jline.io.LDESResultIO.saveLN`.
        """
        result = LNLDESResult()
        dims = data.get('dimensions', {})
        for key in ('nidx', 'nhosts', 'ntasks', 'nentries', 'nacts', 'ncalls',
                    'tshift', 'eshift', 'ashift'):
            setattr(result, key, int(dims.get(key, 0)))
        result.names = list(dims.get('names', []))

        metrics = data.get('metrics', {})
        for key in ('QLN', 'ULN', 'RLN', 'WLN', 'TLN', 'ALN', 'ZLN',
                    'UCallLN', 'TCallLN', 'UEntryClassLN'):
            val = metrics.get(key)
            if val is not None:
                setattr(result, key, np.asarray(
                    self._json_array_to_numpy(val), dtype=float).reshape(-1))

        ci = data.get('confidenceIntervals', {})
        for key in ('QLNCI', 'ULNCI', 'RLNCI', 'TLNCI'):
            val = ci.get(key)
            if val is not None:
                setattr(result, key, np.asarray(
                    self._json_array_to_numpy(val), dtype=float).reshape(-1))

        result.cache_metrics = data.get('cacheMetrics', [])

        samples = data.get('entryRespTimeSamples')
        if samples is not None:
            result.entryRespTimeSamples = [
                np.asarray(row, dtype=float) if row else np.empty(0)
                for row in samples]
        return result

    def _parse_result_json(self, path: str) -> LDESResult:
        """Parse LDES result JSON into LDESResult dataclass."""
        with open(path, 'r') as f:
            data = json.load(f)

        # Check for error
        if 'error' in data:
            raise RuntimeError(f"LDES solver error: {data['error']}")

        # A layered run emits a DIFFERENT document, keyed by LQN element rather
        # than by (station, class); it is declared by modelType.
        if data.get('modelType') == 'LayeredNetwork':
            return self._parse_ln_result_json(data)

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

        # Parse the busy period measurement (present when the simulation was run
        # with --busyperiod). One entry per target: a station, a station-class
        # pair, or a declared subnetwork.
        bp = data.get('busyPeriods')
        if bp is not None:
            result.busy_periods = []
            for tgt in bp.get('targets', []):
                result.busy_periods.append({
                    'name': tgt.get('name'),
                    'stations': [int(v) for v in tgt.get('stations', [])],
                    'class': int(tgt.get('class', -1)),
                    'mean': np.asarray(tgt.get('mean', []), dtype=float),
                    'count': np.asarray(tgt.get('count', []), dtype=float),
                })

        # Parse per-cache hit/miss/latency metrics (keyed by cache node name)
        cm = data.get('cacheMetrics', {})
        if cm:
            result.cache_metrics = {}
            for cname, cdata in cm.items():
                entry = {}
                for k in ('hit', 'delayed', 'miss', 'latency', 'hitList', 'itemProb', 'listCost'):
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

        # HOW LONG THE RUN SHOULD HAVE BEEN, when the caller asked for it. The
        # engine's half-width at the configured confidence over the events it
        # ran pins the ASYMPTOTIC variance, which is the quantity a run length
        # is planned from -- not the stationary variance, which on M/M/1 differs
        # from it by a factor blowing up like (1-rho)^-2.
        plan_spec = None
        cfg = getattr(self.options, 'config', None)
        if cfg is not None:
            plan_spec = cfg.get('runLengthPlan') if isinstance(cfg, dict) \
                else getattr(cfg, 'runLengthPlan', None)
        if plan_spec is not None and getattr(result, 'QNCI', None) is not None:
            from ....api.sim import sim_runlength_plan
            rel = 0.05
            conf = 0.95
            if isinstance(plan_spec, dict):
                rel = float(plan_spec.get('relprecision', rel))
                conf = float(plan_spec.get('confidence', conf))
            elif isinstance(plan_spec, (int, float)) and plan_spec > 0:
                rel = float(plan_spec)
            # The ACTUAL number of events simulated where the engine reports
            # it, not the budget: LDES stops early on convergence, and planning
            # from a budget it never spent would overstate N and so overstate
            # the asymptotic variance.
            used = int(data.get('totalSimulatedEvents', 0) or data.get('events', 0) or 0)
            if used <= 0:
                used = int(getattr(self.options, 'samples', 0) or 0)
            if used > 0:
                result.runLengthPlan = sim_runlength_plan(
                    np.asarray(result.QN, dtype=float), np.asarray(result.QNCI, dtype=float),
                    used, relPrecision=rel, confidence=conf)

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


        # `respTimeSamples` appears at TOP LEVEL under --respt-samples and
        # inside `transient` under --trajectory; the two are the same
        # measurement and whichever is present is read. Looking only in the
        # transient block, as this reader did, makes every --respt-samples run
        # report no samples at all -- the document has them and the parser
        # never sees them.
        rts = tran.get('respTimeSamples') if tran else None
        if rts is None:
            rts = data.get('respTimeSamples')
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

    def supportsTransientAnalysis(self):
        """Transient averages are available (simulation restricted to options.timespan)."""
        return True

    supports_transient_analysis = supportsTransientAnalysis

    def runAnalyzer(self, trajectory: bool = False, export_histogram: bool = False,
                    respt_samples: bool = False) -> LDESResult:
        """
        Run the LDES simulation via subprocess.

        Args:
            trajectory: If True, request trajectory data (QNt, UNt, TNt, t) from ldes.jar.
            respt_samples: If True, request the per-job response time samples the
                empirical response-time CDF is built from.
            export_histogram: If True, request the exact joint-state residence-time
                histogram (used by get_avg_reward to evaluate arbitrary rewards).

        Returns:
            LDESResult containing performance metrics
        """
        start_time = time.time()

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
                                        respt_samples=respt_samples,
                                        trajectory=trajectory,
                                        export_histogram=export_histogram)
            else:
                # Build the ordered CLI runner candidates (native, then jar)
                runners = self._build_cli_runners(model_path, result_path,
                                                  respt_samples=respt_samples,
                                                  trajectory=trajectory,
                                                  export_histogram=export_histogram)
                cmd = runners[0]

                if self.options.verbose not in ('silent',):
                    print(f"SolverLDES command: {' '.join(cmd)}")

                # Wall-clock time budget (options.timeout, seconds). The CLI also gets
                # a cooperative --maxtime flag (see _build_cli_args); this subprocess
                # timeout is the hard outer bound. On expiry the process is killed and
                # an empty result flagged as timed out is returned. An infinite
                # budget, the default, imposes no wall-clock bound at all.
                import math as _math
                _tmo = float(getattr(self.options, 'timeout', float('inf')))
                # The cooperative --maxtime flag does the real early stop; give the
                # process grace to finish writing results and shut down the JVM before
                # the hard subprocess bound kills it (otherwise no result file is
                # produced). An INFINITE budget is no budget: the old 600 s fallback
                # killed a merely SLOW simulation and handed back an empty result, so
                # a loaded host read as a solver that produced nothing.
                _sub_timeout = (_tmo + 30.0) if _math.isfinite(_tmo) and _tmo > 0 else None
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
                        empty.runtime = time.time() - start_time
                        # Cache it like the success path: _computeAvgMetrics reads
                        # self._result after calling runAnalyzer, so returning without
                        # assigning left it None and raised AttributeError on .QN.
                        self._result = empty
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
            list_cost = _flat(cm.get('listCost'))  # [lists] mean storage cost
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
                if list_cost is not None and hasattr(cnode, 'set_result_list_cost'):
                    cnode.set_result_list_cost(np.atleast_1d(list_cost))

    def getAvgBusyPeriod(self, stations=None, jobclass=-1, n=1):
        """Mean busy period of order n for a set of stations, measured by LDES.

        A busy period of order n runs from the instant an arrival raises the jobs
        held by the set to n up to the instant the set falls back below n
        (H. Daduna, "Busy Periods for Subnetworks in Stochastic Networks: Mean
        Value Analysis", J. ACM 35(3), 1988). Periods in progress when the warmup
        ends are discarded, since their start is not observable.

        Args:
            stations: stations forming the subnetwork, as objects, names or station
                indexes; None returns the whole measured table.
            jobclass: job class object, name or index; -1 counts every class.
            n: busy period order or sequence of orders.

        Returns:
            (b, count) with the mean duration(s) and the number of completed
            periods behind each mean, or (table, targets) when stations is None.
            A mean with no completed period is NaN.
        """
        orders = np.atleast_1d(np.asarray(n, dtype=int))
        subnet = self._busy_period_stations(stations)
        cls = self._busy_period_class(jobclass)

        opts = self.options
        if getattr(opts, 'busy_period_orders', 0) < int(orders.max()):
            opts.busy_period_orders = int(orders.max())
            self._result = None
        if len(subnet) > 1:
            wanted = sorted(subnet)
            if not any(sorted(sub) == wanted for sub in opts.busy_period_subnets):
                opts.busy_period_subnets.append(wanted)
                self._result = None
        if self._result is None or getattr(self._result, 'busy_periods', None) is None:
            self.runAnalyzer()
        targets = getattr(self._result, 'busy_periods', None)
        if not targets:
            raise RuntimeError('The LDES engine returned no busy period measurement.')

        if stations is None:
            table = np.full((len(targets), int(orders.max())), np.nan)
            for t, tgt in enumerate(targets):
                mean = np.asarray(tgt['mean'], dtype=float)
                count = np.asarray(tgt['count'], dtype=float)
                row = mean[:table.shape[1]].copy()
                row[count[:table.shape[1]] == 0] = np.nan
                table[t, :] = row
            return table, [tgt['name'] for tgt in targets]

        wanted = sorted(subnet)
        for tgt in targets:
            if sorted(tgt['stations']) == wanted and tgt['class'] == cls:
                mean = np.asarray(tgt['mean'], dtype=float)[orders - 1]
                count = np.asarray(tgt['count'], dtype=float)[orders - 1]
                mean = np.where(count > 0, mean, np.nan)
                if mean.size == 1:
                    return float(mean[0]), float(count[0])
                return mean, count
        raise ValueError('No busy period target matches the requested stations and class.')

    def _busy_period_stations(self, stations):
        """Station indexes of a busy period target, as declared by the caller."""
        if stations is None:
            return []
        if not isinstance(stations, (list, tuple, np.ndarray)):
            stations = [stations]
        out = []
        for st in stations:
            if isinstance(st, (int, np.integer)):
                out.append(int(st))
            else:
                # get_node_index and get_station_index are 1-based, the wire is not
                out.append(int(self.model.get_station_index(st)) - 1)
        return out

    def _busy_period_class(self, jobclass):
        """Zero-based class index of a busy period target, -1 for the aggregate."""
        if jobclass is None:
            return -1
        if isinstance(jobclass, (int, np.integer)):
            return int(jobclass)
        if isinstance(jobclass, str):
            return list(self._class_names).index(jobclass)
        return int(jobclass.get_index0())

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
            self._ensureAvgResults()
        M = len(self._station_names)
        K = len(self._class_names)
        return self._result.QN.copy() if self._result.QN is not None else np.zeros((M, K))

    def getAvgUtil(self) -> np.ndarray:
        """Average utilizations [stations x classes]."""
        return self._computeAvgMetrics()[1]

    def getAvgRespT(self) -> np.ndarray:
        """Average response times [stations x classes]."""
        return self._computeAvgMetrics()[2]

    def getAvgTput(self) -> np.ndarray:
        """Average throughputs [stations x classes]."""
        return self._computeAvgMetrics()[3]

    def getAvgArvR(self) -> np.ndarray:
        """Average arrival rates [stations x classes]."""
        return self._computeAvgMetrics()[4]

    def getAvgWaitT(self) -> np.ndarray:
        """Average waiting times [stations x classes]."""
        return self._computeAvgMetrics()[5]

    get_avg_util = getAvgUtil
    get_avg_respt = getAvgRespT
    get_avg_tput = getAvgTput
    get_avg_arvr = getAvgArvR
    get_avg_waitt = getAvgWaitT

    def getAvg(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Average station metrics (QN, UN, RN, TN, AN, WN) as station x class
        matrices, applying the same visit-based masking used by getAvgTable.
        """
        return self._computeAvgMetrics()

    get_avg = getAvg

    def _computeAvgMetrics(self):
        """Assemble the station x class average matrices from the last
        simulation result (running it if needed)."""
        # A LayeredNetwork is measured over the LQN ELEMENT index space -- hosts,
        # tasks, entries, activities -- and the run parks an LNLDESResult here,
        # which carries QLN..ZLN and no QN..WN at all. Reading it as a (station,
        # class) grid used to die on the missing attribute; refuse by name and
        # send the caller to the layered table instead.
        if self._is_layered():
            raise RuntimeError(
                "the average metrics of a LayeredNetwork are indexed by LQN element, "
                "not by (station, class); use getLNAvgTable() for the table, "
                "getEnsembleAvg() for the raw vectors, and getCdfRespTLN() for the "
                "per-entry response time distribution.")
        if self._result is None:
            self._ensureAvgResults()

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
            # sn.nodetype is a list, and `list == member` is a scalar False, so the
            # three flags below stay off unless the comparison is vectorised first
            nodetype = np.asarray(sn.nodetype)
            hasForkJoin = np.any(nodetype == NodeType.FORK) and np.any(nodetype == NodeType.JOIN)
            hasSPN = np.any(nodetype == NodeType.PLACE) or np.any(nodetype == NodeType.TRANSITION)
            # A Cache node switches classes by item state, which the
            # routing-based visit equations cannot represent: visits
            # downstream of the cache solve to garbage, so trust the
            # simulation there, like fork-join.
            hasCache = np.any(nodetype == NodeType.CACHE)

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
            self._ensureAvgResults()

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
            self._ensureAvgResults()

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

        # five SIGNIFICANT digits like MATLAB's table, not pandas' five decimals
        from line_solver.indexed_table import IndexedTable
        return IndexedTable(pd.DataFrame(rows))

    def getAvgSysTable(self) -> pd.DataFrame:
        """
        Get system-level average performance metrics.

        Returns:
            DataFrame with columns: Chain, SysRespT, SysTput
        """
        if self._result is None:
            self._ensureAvgResults()

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
        # `_node_idx` is NOT an attribute a node has -- the 0-based index is
        # `get_index0()` / `_node_index` -- so the getattr default silently made
        # EVERY call sample node 0 whatever node was asked for.
        node_idx = node if isinstance(node, int) else node.get_index0()

        # The engine writes QNt as [STATION][class] (`Solver_ssj`:
        # `result.QNt = new Matrix[numStations][numClasses]`, indexed by
        # `serviceStation`), so the row is found by station index. Reading it
        # with the STATEFUL index returned another station's trajectory, or all
        # zeros past the end, on any model holding a stateful node that is not a
        # station -- a Router, a Cache, a stateful Fork, a Place -- because
        # sn.isstateful admits those and sn.isstation does not.
        ist = int(np.asarray(sn.nodeToStation).flatten()[node_idx])
        num_time_points = result.t.shape[0]
        num_classes = int(sn.nclasses)
        state = np.zeros((num_time_points, num_classes))

        if 0 <= ist < len(result.QNt):
            for k in range(min(num_classes, len(result.QNt[ist]))):
                class_data = result.QNt[ist][k]
                # An EMPTY series parses to a 1-D array of size zero (the class
                # does not visit the station), and indexing column 0 of it raises
                # rather than yielding nothing, so the shape is checked and not
                # only the None.
                if class_data is not None and class_data.ndim == 2 and class_data.shape[0] > 0:
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

        # One entry per STATEFUL node, as the contract says, but each filled from
        # that node's STATION row: the two index spaces coincide only when every
        # stateful node is a station (see sample() above). A stateful node that is
        # not a station -- a Router, a Cache -- holds no queue-length series in
        # the engine's output and keeps its zero block.
        node_to_stateful = np.asarray(sn.nodeToStateful).flatten()
        node_to_station = np.asarray(sn.nodeToStation).flatten()
        stateful_to_station = {}
        for nd in range(len(node_to_stateful)):
            isf_nd = int(node_to_stateful[nd])
            if isf_nd >= 0:
                stateful_to_station[isf_nd] = int(node_to_station[nd])

        states = []
        for isf in range(int(sn.nstateful)):
            node_state = np.zeros((num_time_points, num_classes))
            ist = stateful_to_station.get(isf, -1)
            if 0 <= ist < len(result.QNt):
                for k in range(min(num_classes, len(result.QNt[ist]))):
                    class_data = result.QNt[ist][k]
                    if class_data is not None and class_data.ndim == 2 \
                            and class_data.shape[0] > 0:
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

    def _state_histogram(self):
        """
        Exact joint-state residence-time histogram of one LDES run.

        Returns ``(space, time)``: ``space`` is (nstates x nstations*nclasses)
        in the station-major, class-minor layout ``ctmc_state_space_aggr``
        builds, and ``time`` the residence time of each row, so
        ``P(state) = t(state)/sum(t)`` is exact on the sampled path.

        WHY THIS AND NOT THE TRAJECTORY. The engine's transient ``QNt`` series
        holds INTERVAL TIME-AVERAGES of the queue length, not the integer states
        the path visits, so comparing it against a state matches only where a
        bucket mean happens to land on an integer. Every probability this
        wrapper reports goes through the histogram, as ``get_avg_reward``
        already does.
        """
        self._result = None
        result = self.runAnalyzer(export_histogram=True)
        space = getattr(result, 'state_histogram_space', None)
        time_arr = getattr(result, 'state_histogram_time', None)
        if space is None or time_arr is None:
            raise RuntimeError(
                'SolverLDES: the engine returned no state histogram, so no state '
                'probability can be read from this run.')
        space = np.atleast_2d(np.asarray(space, dtype=float))
        time_arr = np.asarray(time_arr, dtype=float).flatten()
        return space, time_arr

    def _hist_prob(self, space, time_arr, stations, targets) -> float:
        """
        Residence-time probability of an aggregate joint state.

        ``stations`` is a sequence of 0-based station indices to constrain and
        ``targets`` the matching per-class job-count vectors; a target shorter
        than nclasses constrains only the classes it names.
        """
        total = float(np.sum(time_arr))
        if total <= 0:
            return 0.0
        R = int(self._sn.nclasses)
        match = np.ones(space.shape[0], dtype=bool)
        for ist, tgt in zip(stations, targets):
            tgt = np.asarray(tgt, dtype=float).flatten()
            L = min(len(tgt), R)
            lo = ist * R
            if lo + L > space.shape[1]:
                raise ValueError(
                    'SolverLDES: station %d is past the end of the state histogram, '
                    'which holds %d columns for %d classes.'
                    % (ist + 1, space.shape[1], R))
            match &= np.all(np.abs(space[:, lo:lo + L] - tgt[:L]) < 1e-9, axis=1)
        return float(np.sum(time_arr[match]) / total)

    def _station_of(self, node) -> int:
        """0-based station index of a node, or a refusal naming it."""
        sn = self._sn
        node_idx = node if isinstance(node, int) else node.get_index0()
        ist = int(np.asarray(sn.nodeToStation).flatten()[node_idx])
        if ist < 0:
            raise ValueError(
                'SolverLDES: node %d is not a station; the LDES state histogram '
                'records station queue lengths only.' % (node_idx + 1))
        return ist

    def getProb(self, node, state: Optional[np.ndarray] = None) -> float:
        """
        Steady-state probability of a state at a node.

        The residence-time fraction the exact joint-state histogram of the run
        assigns to the state. ``state`` is a DETAILED node state and is
        aggregated here to per-class job counts, which is the resolution the
        engine's histogram carries (integer queue lengths, no service phases).
        If None, the model's current state for the node is used.

        This is NOT computed from sample(): the transient QNt series holds
        interval time-averages of the queue length, so comparing it against an
        integer state matched almost nowhere and reported a near-zero
        probability for a state the chain mostly occupies.

        Returns:
            Estimated probability (0.0 if the state is never visited).
        """
        from ....lang.state import State

        sn = self._sn
        node_idx = node if isinstance(node, int) else node.get_index0()
        ist = self._station_of(node)

        if state is None:
            # sn.state IS stateful-indexed, unlike the engine's QNt.
            isf = int(np.asarray(sn.nodeToStateful).flatten()[node_idx])
            sn_state = getattr(sn, 'state', None)
            if sn_state is None or isf >= len(sn_state) or sn_state[isf] is None:
                raise ValueError(
                    'SolverLDES.getProb: no state was given and the model carries '
                    'none for this node.')
            state = sn_state[isf]

        _, nir, _, _ = State.toMarginal(sn, node_idx, np.atleast_2d(np.asarray(state)))
        nir = np.atleast_2d(np.asarray(nir, dtype=float))

        # A multi-row state names a SET of states and its probability is the sum
        # over the set. The run is made ONCE and every row weighed against the
        # same histogram; re-solving per row would draw a fresh sample path for
        # each and the sum would not be a probability of anything.
        space, time_arr = self._state_histogram()
        return float(sum(self._hist_prob(space, time_arr, [ist], [nir[r, :]])
                         for r in range(nir.shape[0])))

    def getProbAggr(self, node, state_aggr: Optional[np.ndarray] = None) -> float:
        """
        Aggregated state probability at a node.

        ``state_aggr`` is a per-class job-count vector, already aggregated over
        service phases, which is exactly the resolution the engine's histogram
        carries, so unlike getProb no conversion is applied. It is therefore NOT
        a synonym for getProb: delegating to it would put the counts through
        toMarginal a second time.

        Returns:
            Estimated probability.
        """
        from ....lang.state import State

        sn = self._sn
        node_idx = node if isinstance(node, int) else node.get_index0()
        ist = self._station_of(node)

        if state_aggr is None:
            isf = int(np.asarray(sn.nodeToStateful).flatten()[node_idx])
            sn_state = getattr(sn, 'state', None)
            if sn_state is None or isf >= len(sn_state) or sn_state[isf] is None:
                raise ValueError(
                    'SolverLDES.getProbAggr: no state was given and the model '
                    'carries none for this node.')
            _, state_aggr, _, _ = State.toMarginal(
                sn, node_idx, np.atleast_2d(np.asarray(sn_state[isf])))

        state_aggr = np.atleast_2d(np.asarray(state_aggr, dtype=float))
        space, time_arr = self._state_histogram()
        return float(sum(self._hist_prob(space, time_arr, [ist], [state_aggr[r, :]])
                         for r in range(state_aggr.shape[0])))

    def getProbSys(self) -> float:
        """
        Joint steady-state probability of the current system state.

        Every STATION is constrained to the per-class job counts its current
        state aggregates to. A stateful node that is not a station -- a Router,
        a Cache, a Place -- holds no queue length in the engine's histogram and
        cannot be constrained, so it is excluded and the answer is the joint law
        of the station queue lengths alone.

        Returns:
            Estimated joint probability.
        """
        from ....lang.state import State

        sn = self._sn
        node_to_station = np.asarray(sn.nodeToStation).flatten()
        node_to_stateful = np.asarray(sn.nodeToStateful).flatten()
        sn_state = getattr(sn, 'state', None)

        stations = []
        targets = []
        for nd in range(len(node_to_station)):
            ist = int(node_to_station[nd])
            if ist < 0:
                continue
            isf = int(node_to_stateful[nd])
            if sn_state is None or isf < 0 or isf >= len(sn_state) or sn_state[isf] is None:
                raise ValueError(
                    'SolverLDES.getProbSys: the model carries no current state for '
                    'station %d.' % (ist + 1))
            _, nir, _, _ = State.toMarginal(sn, nd, np.atleast_2d(np.asarray(sn_state[isf])))
            stations.append(ist)
            targets.append(np.atleast_2d(np.asarray(nir, dtype=float))[0, :])

        if not stations:
            return 0.0

        space, time_arr = self._state_histogram()
        return self._hist_prob(space, time_arr, stations, targets)

    def getProbSysAggr(self) -> float:
        """
        Aggregated joint system probability.

        This equals getProbSys() because the engine's state histogram is
        aggregated already: it records integer queue lengths per station and
        class, with no phase resolution, so the detailed and the aggregate joint
        question have the same answer here. getProbSys aggregates the model's
        current state before matching, so no second aggregation is needed.

        Returns:
            Estimated joint probability.
        """
        return self.getProbSys()

    # =========================================================================
    # Transient CDF methods
    # =========================================================================

    def getTranCdfRespT(self):
        """
        Get empirical CDF of response times from simulation samples.

        The engine keeps one set of per-job response time samples, so this is
        the same measured ecdf getCdfRespT reports, under the transient
        getter's name -- the LDES arms serve one curve under every CDF name.
        A LayeredNetwork is refused: the layered engine reports steady state
        only and its response times belong to entries.

        Returns:
            RD[station][class], an (n x 2) array of [F(t), t], or None for a
            (station, class) pair the run observed nothing at.
        """
        if self._is_layered():
            raise RuntimeError(
                "getTranCdfRespT is indexed by (station, class) and this solver holds a "
                "LayeredNetwork, whose response times belong to entries and whose engine "
                "reports steady state only. Use getCdfRespTLN().")
        return self.getCdfRespT()

    def getCdfRespTLN(self):
        """Empirical response time distribution of every ENTRY of a LayeredNetwork.

        The simulated counterpart of ``SolverLN.getCdfRespT``: where the moment3
        pass fits an APH to three moments and convolves, this is the ecdf of the
        response times the run observed, so its tail is measured rather than
        extrapolated. The engine times each request from the instant the entry
        acquires a thread to the instant it replies -- the interval ``RLN``
        averages -- so the mean of this law reproduces that row.

        Returns:
            A list of ``nentries`` items, one per entry in the entry-local index
            space (``lsn.eshift + i``). Each is an ``(n, 2)`` array whose columns
            are ``[F(t), t]``, or None for an entry the run observed nothing at.

        Raises:
            RuntimeError: if the model is not a LayeredNetwork, or the run kept
                no samples.
        """
        if not self._is_layered():
            raise RuntimeError(
                "getCdfRespTLN requires a LayeredNetwork; this solver holds a Network, "
                "whose response times are indexed by (station, class). Use getCdfRespT().")

        result = self.runAnalyzer(respt_samples=True)
        samples_all = getattr(result, 'entryRespTimeSamples', None)
        if not samples_all:
            raise RuntimeError(
                "The LDES run returned no entry response time samples, so an empirical "
                "CDF cannot be built. Increase options.samples, or shorten the warmup.")

        RD = []
        for samples in samples_all:
            if samples is None or len(samples) == 0:
                RD.append(None)
                continue
            x = np.sort(np.asarray(samples, dtype=float).ravel())
            F = np.arange(1, x.size + 1, dtype=float) / x.size
            # Collapse repeated observations, keeping the LARGEST CDF value at
            # each distinct time: a tie left expanded makes the ecdf multivalued.
            xu, last_idx = np.unique(x[::-1], return_index=True)
            last_idx = x.size - 1 - last_idx
            RD.append(np.column_stack([F[last_idx], xu]))
        return RD

    get_cdf_resp_t_ln = getCdfRespTLN

    def _lnResultOrRun(self, caller: str) -> 'LNLDESResult':
        """The LNLDESResult of this model, running the simulation if needed."""
        if not self._is_layered():
            raise RuntimeError(
                "%s requires a LayeredNetwork; this solver holds a Network." % caller)
        if not isinstance(self._result, LNLDESResult):
            self.runAnalyzer()
        if not isinstance(self._result, LNLDESResult):
            raise RuntimeError(
                "%s expected a layered ldes-result and the run returned %s; the engine "
                "was given a Network document." % (caller, type(self._result).__name__))
        return self._result

    def getEnsembleAvg(self):
        """Mean metrics over the LQN element index space.

        Returns:
            (QLN, ULN, RLN, WLN, TLN), each an ``nidx`` vector in the absolute
            LQN index space -- hosts, then tasks, then entries, then activities.
            This is the raw measurement; :meth:`getLNAvgTable` is the same data
            masked and named the way SolverLN and LQNS report it.
        """
        r = self._lnResultOrRun('getEnsembleAvg')
        z = lambda v: np.zeros(r.nidx) if v is None else np.asarray(v, dtype=float).ravel()
        return z(r.QLN), z(r.ULN), z(r.RLN), z(r.WLN), z(r.TLN)

    get_ensemble_avg = getEnsembleAvg

    @staticmethod
    def _hostProcessorMult(lsn, idx: int) -> float:
        """Multiplicity of the processor `idx` ultimately runs on.

        Twin of SolverLDES.hostProcessorMult in the JAR. An infinite server
        reports 1.0 so that its utilization keeps the mean-busy-servers value,
        which is the LINE convention and may exceed 1.
        """
        from ....layered import LayeredNetworkElement
        nidx = int(lsn.nidx)
        cur = int(idx)
        types = np.asarray(lsn.type).ravel()
        parents = np.asarray(lsn.parent).ravel()
        mult = np.asarray(lsn.mult, dtype=float).ravel()
        for _ in range(nidx + 1):
            if cur < 0 or cur >= nidx:
                return 1.0
            if int(types[cur]) == int(LayeredNetworkElement.PROCESSOR):
                m = float(mult[cur])
                return m if (m > 0 and not np.isinf(m)) else 1.0
            p = int(parents[cur])
            if p < 0 or p == cur:
                return 1.0
            cur = p
        return 1.0

    def getLNAvgTable(self) -> pd.DataFrame:
        """Mean metrics of a LayeredNetwork, one row per LQN element.

        The simulated twin of ``SolverLN.get_avg_table``, carrying the same
        columns and the same NaN mask so the two can be compared cell for cell:
        a processor has no queue length, response time or throughput of its own,
        and a task has no response time. ArvR is not measured here.

        Returns:
            A DataFrame with columns Node, NodeType, QLen, Util, RespT, ResidT,
            ArvR, Tput, indexed by the absolute LQN element index.
        """
        from ....layered import LayeredNetworkElement
        r = self._lnResultOrRun('getLNAvgTable')
        lsn = self.model.getStruct()
        types = np.asarray(lsn.type).ravel()
        names = list(np.asarray(lsn.names).ravel())
        isref = np.asarray(lsn.isref).ravel() if lsn.isref is not None else None
        vec = lambda v: (np.full(r.nidx, np.nan) if v is None
                         else np.asarray(v, dtype=float).ravel())
        QLN, ULN, RLN, WLN, TLN = vec(r.QLN), vec(r.ULN), vec(r.RLN), vec(r.WLN), vec(r.TLN)

        def type_name(i):
            t = int(types[i])
            if t == int(LayeredNetworkElement.PROCESSOR):
                return 'Processor'
            if t == int(LayeredNetworkElement.TASK):
                return 'RefTask' if (isref is not None and bool(isref[i])) else 'Task'
            if t == int(LayeredNetworkElement.ENTRY):
                return 'Entry'
            if t == int(LayeredNetworkElement.ACTIVITY):
                return 'Activity'
            return 'Unknown'

        rows = []
        for i in range(int(lsn.nidx)):
            t = int(types[i])
            is_proc = t == int(LayeredNetworkElement.PROCESSOR)
            is_task = t == int(LayeredNetworkElement.TASK)
            # Per-server fraction, as everywhere else in LINE: the engine
            # accumulates mean busy servers over the multiplicity.
            m = self._hostProcessorMult(lsn, i)
            util = ULN[i] / m if m > 1.0 else ULN[i]
            rows.append({
                'Node': names[i] if i < len(names) else str(i),
                'NodeType': type_name(i),
                'QLen': np.nan if is_proc else QLN[i],
                'Util': util,
                'RespT': np.nan if (is_proc or is_task) else RLN[i],
                'ResidT': WLN[i],
                'ArvR': np.nan,
                'Tput': np.nan if is_proc else TLN[i],
            })
        return pd.DataFrame(rows)

    get_ln_avg_table = getLNAvgTable

    def getCdfRespT(self, R=None):
        """
        Empirical response time CDF, measured on the simulated sample path.

        A simulator must report what it OBSERVED. The analytical solvers fall
        back to an exponential law carrying the right mean, which says nothing
        about the tail; the LDES engine records every per-job response time
        (Solver_ssj.responseTimeSamples -> LDESResult.respTimeSamples -> the
        JSON respTimeSamples block, exported under --respt-samples), so the CDF
        here is the ecdf of those samples.

        Ported from MATLAB @SolverLDES/getCdfRespT.m. The return shape is
        SolverJMT's, the other empirical simulator on this side: RD[i][k] is an
        (n x 2) array whose columns are [F(t), t], the convention every
        getCdfRespT follows; swapping the columns yields a CDF that looks like
        a time axis.

        A LayeredNetwork has no (station, class) grid: its response times
        belong to ENTRIES, so it is dispatched to getCdfRespTLN and the
        per-entry laws are returned, exactly as the reference does.

        Args:
            R: optional response time handles, accepted for signature
               compatibility with the analytical solvers and not read.

        Returns:
            RD[station][class], an (n x 2) array of [cdf, time], or None for a
            (station, class) pair the run observed nothing at. A pair with no
            observation is left empty rather than filled with a guess. For a
            LayeredNetwork, the per-entry list of getCdfRespTLN.
        """
        if self._is_layered():
            return self.getCdfRespTLN()
        sn = self._get_network_struct()
        M, K = sn.nstations, sn.nclasses
        RD = [[None for _ in range(K)] for _ in range(M)]

        # The samples are a distributional output, not a trajectory, so they
        # have their own flag; requesting it also forces the jar runner, since
        # the prebuilt AOT binary predates it.
        result = self.runAnalyzer(respt_samples=True)
        samples_all = getattr(result, 'respTimeSamples', None)
        if not samples_all:
            raise ValueError(
                "The LDES run returned no response time samples, so an empirical "
                "CDF cannot be built. Increase options.samples, or use "
                "getPerctRespT(..., 'forktail') for the analytical fork-join tail.")

        for i, station_samples in enumerate(samples_all):
            if i >= M:
                break
            for k, samples in enumerate(station_samples):
                if k >= K or samples is None or len(samples) == 0:
                    continue
                x = np.sort(np.asarray(samples, dtype=float).ravel())
                F = np.arange(1, x.size + 1, dtype=float) / x.size
                # Collapse repeated observations, keeping the LARGEST CDF value
                # at each distinct time: a tie left expanded makes the ecdf
                # multivalued, and an interpolating consumer then reads a
                # quantile off whichever duplicate it happened to hit.
                xu, last_idx = np.unique(x[::-1], return_index=True)
                last_idx = x.size - 1 - last_idx
                RD[i][k] = np.column_stack([F[last_idx], xu])

        return RD

    def getTranCdfPassT(self):
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

    @result.setter
    def result(self, value: Optional[LDESResult]) -> None:
        # NetworkSolver._clearResultStores() (base.py) sets both `_result`
        # and `result` directly, expecting `result` to be assignable (a
        # plain attribute on most backends); without this setter, reset()
        # raises before a fresh solve, which breaks SolverLDES used as a
        # SolverLN per-layer factory (SolverLN.post() calls reset() every
        # outer iteration).
        self._result = value

    def getName(self) -> str:
        """Get solver name."""
        return "LDES"

    def get_name(self) -> str:
        """Get solver name (Python convention)."""
        return self.getName()

    def listValidMethods(self):
        """Valid methods for this solver, SolverLDES.m verbatim.

        'parallel' asks the engine for INDEPENDENT REPLICATIONS and the mean
        over them, which is what its parallel analyzer is; it is not a second
        engine. The CLI argument builder turns the name into --replications,
        taking options.replications when set and 8 otherwise -- the same default
        the SSA parallel analyzer uses for its replica count.
        """
        return ['default', 'parallel']

    list_valid_methods = listValidMethods

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
            'JoinPartial',  # quorum join: fires at the k-th sibling, stragglers discarded on arrival
            'Place', 'Transition',
            'Linkage', 'Enabling', 'Inhibiting', 'Timing', 'Firing', 'Storage',
            'Logger', 'LogTunnel',
            'Cache', 'CacheClassSwitcher', 'CacheRetrieval', 'CacheItemSize',
            'ReplacementStrategy_LRU', 'ReplacementStrategy_FIFO', 'ReplacementStrategy_RR', 'ReplacementStrategy_SFIFO',
            'ReplacementStrategy_HLRU', 'ReplacementStrategy_CLIMB', 'ReplacementStrategy_QLRU',
            'Buffer', 'Region',
            'Exp', 'Erlang', 'HyperExp', 'PH', 'APH', 'Coxian', 'Cox2', 'MAP', 'DMAP', 'MMAP', 'BMAP', 'MMPP2', 'ME', 'RAP', 'Immediate', 'Disabled', 'Replayer', 'Trace',  # Trace is an alias of Replayer
            'Det', 'Uniform', 'Gamma', 'Pareto', 'Weibull', 'Lognormal',
            'Geometric',  # Lattice-valued interarrival/service time on {1,2,...} (Geo/Geo/1 and slotted models)
            'Bernoulli', 'Binomial', 'Poisson',  # Counting distributions; zero atom becomes an immediate interval (continuous mode only)
            'NHPP',
            # Time-inhomogeneous MAP: the piecewise-constant (D0,D1) schedule is
            # simulated exactly by carrying the phase across a breakpoint. PHt service
            # is walked from the SERVICE START epoch, so processor sharing, preemption,
            # load dependence and heterogeneous servers are rejected at runtime.
            'MAPt',
            'PHt',
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
            'RoutingStrategy_SDR',  # Krzesinski (1987) product-form state-dependent routing
            'OpenClass', 'ClosedClass', 'SelfLoopingClass',
            'OpenSignal', 'ClosedSignal',
            'SignalType_NEGATIVE', 'SignalType_REPLY', 'SignalType_CATASTROPHE',
            'SignalBatchRemoval', 'SignalRemovalPolicy',
            'LoadDependence', 'ClassDependence', 'JointDependence', 'Balking', 'Reneging', 'Retrial',
            # Engine simulates the SETUP/DELAYOFF server states.
            'SetupDelayOff',
            # set_breakdown: the server alternates up/down on the breakdownMu/repairMu
            # clocks, a job in service holds its residual work across the outage
            # (preemptive resume), and downServiceRates runs the server at a degraded
            # speed instead of stopping it. Rejected at runtime in slotted mode and with
            # time-inhomogeneous service (MAPt/PHt/NHPP).
            'Breakdown',
            # Queue.add_server_type: the engine keeps one pool per server type,
            # assigns each job a type from the pools compatible with its class and
            # serves it at that pool's own rate, so the pools are an exact
            # sample-path feature rather than a flattened nservers. Only the JVM
            # engine (common/ldes.jar) implements them -- the native C++ binary
            # refuses the model by name (ldes_engine_reject), which is what makes
            # the runner list fall through to the jar.
            'HeteroServers',
            # Source.set_arrival_batch: save_model writes arrivalBatch and the
            # engine reads sn.arrivalbatch, so a batch releases several jobs at
            # one arrival epoch on the sample path.
            'BatchArrival',
            # c-server stations (sn.nservers) and finite buffers with their drop
            # rule (sn.cap/classcap, the 'Buffer' marker above) are simulated
            # directly by both engines.
            'MultiServer', 'FiniteCapacity',
        }

    @staticmethod
    def supports(model) -> bool:
        """Check if model is supported.

        Mirrors MATLAB SolverLDES.supports. No supports() existed anywhere in
        the MRO, so calling it raised AttributeError and the solver had no gate.
        """
        from ...base import supports_via_featureset
        return supports_via_featureset(SolverLDES, model)
