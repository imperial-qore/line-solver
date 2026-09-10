"""
JMT Solver handler - Native Python implementation.

Calls JMT via subprocess.

Port from:



Note: This is a simplified implementation supporting basic queueing networks.
Complex features (caches, transitions, etc.) may require the Java implementation.
"""

import math
import numpy as np
import subprocess
import tempfile
import os
import shutil
import platform
import urllib.request
from dataclasses import dataclass
from typing import Optional, Tuple, List, Dict, Any
from xml.etree import ElementTree as ET
from xml.dom import minidom
import time

from .runner import run_jmt, has_java, result_path_for, JMTBackendError

from ...sn import (
    NetworkStruct,
    NodeType,
    SchedStrategy,
    sn_get_demands_chain,
    sn_deaggregate_chain_results,
    sn_get_arvr_from_tput,
)
from ....constants import ProcessType, PollingType, EventType, JoinStrategy


@dataclass
class SolverJMTOptions:
    """Options for JMT solver."""
    method: str = 'jsim'
    samples: int = 10000
    seed: int = 23000
    max_simulated_time: float = float('inf')
    conf_int: float = 0.99
    max_rel_err: float = 0.03
    verbose: bool = False
    keep: bool = False
    # Backend selection, see api/solvers/jmt/runner.py. rest_url points at a
    # JMT REST server; container overrides the Docker image used when no local
    # JVM exists. Both empty means the local JVM plus common/JMT.jar.
    rest_url: Optional[str] = None
    container: Optional[str] = None
    timeout: float = float('inf')


@dataclass
class SolverJMTReturn:
    """
    Result of JMT solver handler.

    Attributes:
        Q: Mean queue lengths (M x K)
        U: Utilizations (M x K)
        R: Response times (M x K)
        T: Throughputs (M x K)
        A: Arrival rates (M x K)
        W: Waiting times (M x K)
        C: Cycle times (1 x K)
        X: System throughputs (1 x K)
        runtime: Runtime in seconds
        method: Method used
    """
    Q: Optional[np.ndarray] = None
    U: Optional[np.ndarray] = None
    R: Optional[np.ndarray] = None
    T: Optional[np.ndarray] = None
    A: Optional[np.ndarray] = None
    # Finite Capacity Region (FCR) region-level metrics (nregions x nclasses),
    # None when the model has no regions. Ufcr/Afcr are NaN (JMT omits them).
    Qfcr: Optional[np.ndarray] = None
    Ufcr: Optional[np.ndarray] = None
    Rfcr: Optional[np.ndarray] = None
    Wfcr: Optional[np.ndarray] = None
    Tfcr: Optional[np.ndarray] = None
    Afcr: Optional[np.ndarray] = None
    # Region loss metrics for getAvgRegionLossTable (solver-agnostic names): the
    # carried throughput TNfcr = Tfcr, and DropRateNfcr = Afcr - Tfcr (offered
    # minus carried), matching the LDES result field names.
    TNfcr: Optional[np.ndarray] = None
    DropRateNfcr: Optional[np.ndarray] = None
    W: Optional[np.ndarray] = None
    C: Optional[np.ndarray] = None
    X: Optional[np.ndarray] = None
    runtime: float = 0.0
    method: str = "jsim"
    timedOut: bool = False
    # log normalizing constant, reported by the JMVA engine only (its
    # <normconst logValue>); NaN on the simulation path and on the JMVA
    # algorithms that do not compute one.
    logNormConstAggr: float = float('nan')


def _get_jmt_jar_path() -> str:
    """Get path to JMT.jar, downloading if necessary."""
    # Look in common/ directory
    package_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    python_dir = os.path.dirname(package_dir)
    root_dir = os.path.dirname(python_dir)
    common_dir = os.path.join(root_dir, 'common')
    jmt_path = os.path.join(common_dir, 'JMT.jar')

    if os.path.isfile(jmt_path):
        return jmt_path

    # Try to download
    os.makedirs(common_dir, exist_ok=True)
    jmt_url = 'https://line-solver.sourceforge.net/latest/JMT.jar'
    try:
        urllib.request.urlretrieve(jmt_url, jmt_path)
        return jmt_path
    except Exception as e:
        raise RuntimeError(
            f"JMT.jar not found and download failed: {e}\n"
            f"Please manually download from {jmt_url} and place in {common_dir}"
        )


def is_jmt_available() -> bool:
    """Check if JMT is available."""
    # Check for JMT
    from .runner import has_java
    if not has_java():
        return False

    # Check for JMT.jar
    try:
        jmt_path = _get_jmt_jar_path()
        return os.path.isfile(jmt_path)
    except RuntimeError:
        return False


def _get_sched_strategy_class(sched: SchedStrategy) -> str:
    """Map scheduling strategy to JMT QueueGetStrategy class name.

    Note: JMT handles PS/GPS/DPS through PSStrategies, not QueueGetStrategies.
    For QueueGetStrategies, only FCFS and LCFS are available.
    PS and other strategies are handled by the Server section's PSStrategy.
    """
    strategy_map = {
        SchedStrategy.FCFS: "jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy",
        SchedStrategy.LCFS: "jmt.engine.NetStrategies.QueueGetStrategies.LCFSstrategy",
        # PS/SIRO/INF use FCFS queue get strategy - actual scheduling is in PSStrategy
        SchedStrategy.PS: "jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy",
        SchedStrategy.SIRO: "jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy",
        SchedStrategy.INF: "jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy",
    }
    return strategy_map.get(sched, "jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy")


def _get_polling_get_strategy_class(polling_type: PollingType) -> str:
    """Map polling type to JMT QueueGetStrategy class name for polling queues.

    Args:
        polling_type: PollingType enum value

    Returns:
        JMT class path for the polling get strategy
    """
    if polling_type == PollingType.GATED:
        return "jmt.engine.NetStrategies.QueueGetStrategies.GatedPollingGetStrategy"
    elif polling_type == PollingType.EXHAUSTIVE:
        return "jmt.engine.NetStrategies.QueueGetStrategies.ExhaustivePollingGetStrategy"
    elif polling_type == PollingType.KLIMITED:
        return "jmt.engine.NetStrategies.QueueGetStrategies.LimitedPollingGetStrategy"
    elif polling_type == PollingType.DECREMENTING:
        raise ValueError("JMT does not support the decrementing (semiexhaustive) "
                         "polling discipline; use the LDES solver.")
    else:
        return "jmt.engine.NetStrategies.QueueGetStrategies.ExhaustivePollingGetStrategy"


def _write_polling_get_strategy(queue_elem: ET.Element, polling_type: PollingType, polling_k: int = 1):
    """Write polling get strategy to JMT Queue section.

    Args:
        queue_elem: Queue XML section element
        polling_type: PollingType enum value
        polling_k: K value for KLIMITED polling
    """
    strategy_param = ET.SubElement(queue_elem, 'parameter')
    strategy_param.set('classPath', _get_polling_get_strategy_class(polling_type))
    strategy_param.set('name', 'FCFSstrategy')

    # For KLIMITED polling, add the pollingKValue subparameter
    if polling_type == PollingType.KLIMITED:
        polling_k_param = ET.SubElement(strategy_param, 'subParameter')
        polling_k_param.set('classPath', 'java.lang.Integer')
        polling_k_param.set('name', 'pollingKValue')
        value = ET.SubElement(polling_k_param, 'value')
        value.text = str(int(polling_k))


def _write_switchover_service_time_strategy(parent: ET.Element, procid: int, proc, rate: float, scv: float = 1.0):
    """Write service time strategy for switchover distributions.

    Args:
        parent: Parent XML element for the subParameter
        procid: ProcessType ID for the distribution
        proc: Process data (distribution parameters, e.g., list of matrices for PH)
        rate: Service rate (for simple distributions)
        scv: Squared coefficient of variation (for distributions that need it)
    """
    service_time_node = ET.SubElement(parent, 'subParameter')

    if procid == ProcessType.DISABLED:
        service_time_node.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.DisabledServiceTimeStrategy')
        service_time_node.set('name', 'DisabledServiceTimeStrategy')
        return

    if procid == ProcessType.IMMEDIATE:
        service_time_node.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ZeroServiceTimeStrategy')
        service_time_node.set('name', 'ZeroServiceTimeStrategy')
        return

    service_time_node.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy')
    service_time_node.set('name', 'ServiceTimeStrategy')

    # Distribution node
    distr_node = ET.SubElement(service_time_node, 'subParameter')
    distr_par_node = ET.SubElement(service_time_node, 'subParameter')

    if procid == ProcessType.EXP:
        distr_node.set('classPath', 'jmt.engine.random.Exponential')
        distr_node.set('name', 'Exponential')
        distr_par_node.set('classPath', 'jmt.engine.random.ExponentialPar')
        distr_par_node.set('name', 'distrPar')

        lambda_param = ET.SubElement(distr_par_node, 'subParameter')
        lambda_param.set('classPath', 'java.lang.Double')
        lambda_param.set('name', 'lambda')
        value = ET.SubElement(lambda_param, 'value')
        value.text = f'{rate:.12f}'

    elif procid == ProcessType.DET:
        distr_node.set('classPath', 'jmt.engine.random.DeterministicDistr')
        distr_node.set('name', 'Deterministic')
        distr_par_node.set('classPath', 'jmt.engine.random.DeterministicDistrPar')
        distr_par_node.set('name', 'distrPar')

        t_param = ET.SubElement(distr_par_node, 'subParameter')
        t_param.set('classPath', 'java.lang.Double')
        t_param.set('name', 't')
        value = ET.SubElement(t_param, 'value')
        value.text = f'{1.0/rate:.12f}' if rate > 0 else '0.0'

    elif procid == ProcessType.ERLANG:
        phases = len(proc[0]) if proc and isinstance(proc, (list, tuple)) and len(proc) > 0 else 2
        distr_node.set('classPath', 'jmt.engine.random.Erlang')
        distr_node.set('name', 'Erlang')
        distr_par_node.set('classPath', 'jmt.engine.random.ErlangPar')
        distr_par_node.set('name', 'distrPar')

        alpha_param = ET.SubElement(distr_par_node, 'subParameter')
        alpha_param.set('classPath', 'java.lang.Double')
        alpha_param.set('name', 'alpha')
        value = ET.SubElement(alpha_param, 'value')
        value.text = f'{rate * phases:.12f}'

        r_param = ET.SubElement(distr_par_node, 'subParameter')
        r_param.set('classPath', 'java.lang.Long')
        r_param.set('name', 'r')
        value = ET.SubElement(r_param, 'value')
        value.text = str(phases)

    elif procid == ProcessType.HYPEREXP:
        # Recover the full (probs, rates) pair; an n>2 HyperExp is emitted as the
        # equivalent APH rather than being truncated to 2 phases.
        params = _hyperexp_probs_rates(proc)
        if params is None:
            params = (np.array([0.5, 0.5]), np.array([rate, rate]))
        _fill_hyperexp_distribution(distr_node, distr_par_node, params[0], params[1])

    elif procid == ProcessType.PARETO:
        # Pareto distribution - reconstruct shape/scale from scv and rate
        shape = np.sqrt(1.0 + 1.0 / scv) + 1.0
        scale_val = (1.0 / rate) * (shape - 1.0) / shape

        distr_node.set('classPath', 'jmt.engine.random.Pareto')
        distr_node.set('name', 'Pareto')
        distr_par_node.set('classPath', 'jmt.engine.random.ParetoPar')
        distr_par_node.set('name', 'distrPar')

        alpha_param = ET.SubElement(distr_par_node, 'subParameter')
        alpha_param.set('classPath', 'java.lang.Double')
        alpha_param.set('name', 'alpha')
        value = ET.SubElement(alpha_param, 'value')
        value.text = f'{shape:.12f}'

        k_param = ET.SubElement(distr_par_node, 'subParameter')
        k_param.set('classPath', 'java.lang.Double')
        k_param.set('name', 'k')
        value = ET.SubElement(k_param, 'value')
        value.text = f'{scale_val:.12f}'

    elif procid == ProcessType.GAMMA:
        # Gamma distribution
        gamma_alpha = 1.0 / scv
        gamma_beta = scv / rate

        distr_node.set('classPath', 'jmt.engine.random.GammaDistr')
        distr_node.set('name', 'Gamma')
        distr_par_node.set('classPath', 'jmt.engine.random.GammaDistrPar')
        distr_par_node.set('name', 'distrPar')

        alpha_param = ET.SubElement(distr_par_node, 'subParameter')
        alpha_param.set('classPath', 'java.lang.Double')
        alpha_param.set('name', 'alpha')
        value = ET.SubElement(alpha_param, 'value')
        value.text = f'{gamma_alpha:.12f}'

        beta_param = ET.SubElement(distr_par_node, 'subParameter')
        beta_param.set('classPath', 'java.lang.Double')
        beta_param.set('name', 'beta')
        value = ET.SubElement(beta_param, 'value')
        value.text = f'{gamma_beta:.12f}'

    else:
        # Default to exponential for unsupported types
        distr_node.set('classPath', 'jmt.engine.random.Exponential')
        distr_node.set('name', 'Exponential')
        distr_par_node.set('classPath', 'jmt.engine.random.ExponentialPar')
        distr_par_node.set('name', 'distrPar')

        lambda_param = ET.SubElement(distr_par_node, 'subParameter')
        lambda_param.set('classPath', 'java.lang.Double')
        lambda_param.set('name', 'lambda')
        value = ET.SubElement(lambda_param, 'value')
        value.text = f'{rate:.12f}'


def _server_pools(node_idx: int, sn: NetworkStruct, nclasses: int):
    """Server pools and job parallelism of a node, or None when it declares neither.

    JMT's Server section takes classParallelism, serverNames, serversPerServerType,
    serverCompatibilities and schedulingPolicy as one positional block of its
    constructor (jmt.engine.NodeSections.Server), so the five are emitted together
    or not at all, and always after the service strategies. A station declaring
    parallelism alone is therefore given one synthetic pool holding all of its
    servers, since the pool counts, not numberOfServers, size the server pool once
    any pool is declared.

    Returns:
        dict with keys names, counts, compat, policy, parallelism; or None
    """
    from ....lang.base import HeteroSchedPolicy

    param = None
    if sn.nodeparam and node_idx in sn.nodeparam:
        param = sn.nodeparam[node_idx]
    if not isinstance(param, dict):
        return None

    n_types = int(param.get('nservertypes', 0) or 0)
    names = param.get('servertypenames')
    counts = param.get('serverspertype')
    compat = param.get('servercompat')
    has_types = n_types > 0 and names is not None and counts is not None and compat is not None

    declared = param.get('serverparallelism')
    has_parallelism = declared is not None and any(float(x) > 1 for x in np.asarray(declared).ravel())

    if not has_types and not has_parallelism:
        return None

    parallelism = np.ones(nclasses, dtype=int)
    if declared is not None:
        flat = np.asarray(declared).ravel()
        for r in range(min(nclasses, flat.size)):
            parallelism[r] = max(1, int(flat[r]))

    if has_types:
        policy = param.get('heteroschedpolicy') or HeteroSchedPolicy.ORDER
        return {'names': list(names), 'counts': np.asarray(counts).ravel(),
                'compat': np.asarray(compat).reshape(n_types, -1),
                'policy': policy, 'parallelism': parallelism}

    ist = int(sn.nodeToStation[node_idx])
    nservers = sn.nservers[ist] if len(sn.nservers.shape) == 1 else sn.nservers[ist, 0]
    node_name = sn.nodenames[node_idx] if node_idx < len(sn.nodenames) else f'node {node_idx}'
    return {'names': [f'{node_name} - Server Type 1'],
            'counts': np.array([int(nservers)]),
            'compat': np.ones((1, nclasses)),
            'policy': HeteroSchedPolicy.ORDER, 'parallelism': parallelism}


def _write_server_pools(server_elem: ET.Element, node_idx: int, sn: NetworkStruct, classnames: List[str]):
    """Write classParallelism and the heterogeneous pool block of a Server section.

    Mirrors MATLAB saveClassParallelism/saveServerTypeNames/saveServersPerType/
    saveServerCompatibilities/saveHeteroSchedPolicy, in that order.
    """
    K = len(classnames)
    pools = _server_pools(node_idx, sn, K)
    if pools is None:
        return

    par_param = ET.SubElement(server_elem, 'parameter')
    par_param.set('array', 'true')
    par_param.set('classPath', 'java.lang.Integer')
    par_param.set('name', 'classParallelism')
    for r in range(K):
        ref_class = ET.SubElement(par_param, 'refClass')
        ref_class.text = classnames[r]
        sub_param = ET.SubElement(par_param, 'subParameter')
        sub_param.set('classPath', 'java.lang.Integer')
        sub_param.set('name', 'serverParallelism')
        value = ET.SubElement(sub_param, 'value')
        value.text = str(int(pools['parallelism'][r]))

    names_param = ET.SubElement(server_elem, 'parameter')
    names_param.set('array', 'true')
    names_param.set('classPath', 'java.lang.String')
    names_param.set('name', 'serverNames')
    for name in pools['names']:
        sub_param = ET.SubElement(names_param, 'subParameter')
        sub_param.set('classPath', 'java.lang.String')
        sub_param.set('name', 'serverTypesNames')
        value = ET.SubElement(sub_param, 'value')
        value.text = str(name)

    counts_param = ET.SubElement(server_elem, 'parameter')
    counts_param.set('array', 'true')
    counts_param.set('classPath', 'java.lang.Integer')
    counts_param.set('name', 'serversPerServerType')
    for count in pools['counts']:
        sub_param = ET.SubElement(counts_param, 'subParameter')
        sub_param.set('classPath', 'java.lang.Integer')
        sub_param.set('name', 'serverTypesNumOfServers')
        value = ET.SubElement(sub_param, 'value')
        value.text = str(int(count))

    compat_param = ET.SubElement(server_elem, 'parameter')
    compat_param.set('array', 'true')
    compat_param.set('classPath', 'java.lang.Object')
    compat_param.set('name', 'serverCompatibilities')
    for t in range(pools['compat'].shape[0]):
        type_node = ET.SubElement(compat_param, 'subParameter')
        type_node.set('array', 'true')
        type_node.set('classPath', 'java.lang.Boolean')
        type_node.set('name', 'serverTypesCompatibilities')
        for r in range(K):
            class_node = ET.SubElement(type_node, 'subParameter')
            class_node.set('classPath', 'java.lang.Boolean')
            class_node.set('name', 'compatibilities')
            value = ET.SubElement(class_node, 'value')
            value.text = 'true' if pools['compat'][t, r] > 0 else 'false'

    policy_param = ET.SubElement(server_elem, 'parameter')
    policy_param.set('classPath', 'java.lang.String')
    policy_param.set('name', 'schedulingPolicy')
    value = ET.SubElement(policy_param, 'value')
    value.text = pools['policy'].to_jmt_text()

    _warn_hetero_rates(node_idx, sn)


def _warn_hetero_rates(node_idx: int, sn: NetworkStruct):
    """Warn that per-server-type service rates cannot reach the JMT engine.

    JMT keys the ServiceStrategy array of a station by refClass, so its loader
    (jmt.engine.simEngine.SimLoader) keeps one strategy per class however many
    (type, class) entries are written, and every pool of a station ends up serving
    at the class rate. Pool sizes, class compatibilities and the assignment policy
    do cross; set_hetero_service rates do not.
    """
    param = sn.nodeparam[node_idx] if (sn.nodeparam and node_idx in sn.nodeparam) else None
    if not isinstance(param, dict):
        return
    rates = param.get('heterorates')
    if rates is None:
        return
    flat = [float(x) for x in np.asarray(rates).ravel() if float(x) > 0]
    if len(flat) < 2 or max(flat) - min(flat) < 1e-12:
        return
    import warnings
    node_name = sn.nodenames[node_idx] if node_idx < len(sn.nodenames) else f'node {node_idx}'
    warnings.warn(f"JMT keys service strategies by job class, so the per-server-type service "
                  f"rates of station {node_name} cannot be exported; every pool will serve at "
                  f"the class service rate. Use the LDES or CTMC solver for per-type rates.")


def _write_switchover_strategy(server_elem: ET.Element, node_idx: int, sn: NetworkStruct, classnames: List[str]):
    """Write switchover strategy for polling queues.

    Writes the SwitchoverStrategy parameter to the Server section for polling queues.

    Args:
        server_elem: Server XML section element
        node_idx: Node index in the network
        sn: NetworkStruct containing nodeparam with switchover info
        classnames: List of class names
    """
    K = len(classnames)
    is_polling_queue = False
    has_switchover = False

    # Check if this is a polling queue with switchover
    if sn.nodeparam and node_idx in sn.nodeparam:
        nodeparam = sn.nodeparam[node_idx]
        if isinstance(nodeparam, dict):
            # Check if any class has switchover
            for r in range(K):
                if r in nodeparam and isinstance(nodeparam[r], dict):
                    if 'pollingType' in nodeparam[r]:
                        is_polling_queue = True
                    if 'switchoverTime' in nodeparam[r]:
                        has_switchover = True
                    if is_polling_queue:
                        break

    if not is_polling_queue:
        return

    param_node = ET.SubElement(server_elem, 'parameter')
    param_node.set('array', 'true')

    if has_switchover:
        param_node.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategy')
        param_node.set('name', 'SwitchoverStrategy')

        nodeparam = sn.nodeparam[node_idx]
        for r in range(K):
            ref_class = ET.SubElement(param_node, 'refClass')
            ref_class.text = classnames[r]

            # Get switchover distribution for this class
            switchover_time = None
            switchover_proc_id = ProcessType.DISABLED
            rate = 1.0

            if r in nodeparam and isinstance(nodeparam[r], dict):
                switchover_times = nodeparam[r].get('switchoverTime', {})
                switchover_proc_ids = nodeparam[r].get('switchoverProcId', {})

                # For polling, switchover is typically uniform across class transitions
                # Use first available switchover or default
                if switchover_times:
                    first_key = list(switchover_times.keys())[0]
                    switchover_time = switchover_times[first_key]
                    switchover_proc_id = switchover_proc_ids.get(first_key, ProcessType.EXP)

                    # Get rate from distribution
                    if hasattr(switchover_time, 'getMean'):
                        mean = switchover_time.getMean()
                        rate = 1.0 / mean if mean > 0 else 1.0
                    elif hasattr(switchover_time, 'get_mean'):
                        mean = switchover_time.get_mean()
                        rate = 1.0 / mean if mean > 0 else 1.0
                    elif hasattr(switchover_time, '_rate'):
                        rate = switchover_time._rate

            # Get proc data for complex distributions
            proc = None
            if switchover_time and hasattr(switchover_time, 'get_representation'):
                proc = switchover_time.get_representation()
            elif switchover_time and hasattr(switchover_time, '_representation'):
                proc = switchover_time._representation

            _write_switchover_service_time_strategy(param_node, switchover_proc_id, proc, rate)
    else:
        # No switchover - write empty strategy
        param_node.set('classPath', 'java.lang.Object')
        param_node.set('name', 'SwitchoverStrategy')

        for r in range(K):
            ref_class = ET.SubElement(param_node, 'refClass')
            ref_class.text = classnames[r]

            # Empty subparameter for each class with class-to-class structure
            sub_param_row = ET.SubElement(param_node, 'subParameter')
            sub_param_row.set('array', 'true')
            sub_param_row.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategy')
            sub_param_row.set('name', 'SwitchoverStrategy')

            for s in range(K):
                ref_class_s = ET.SubElement(sub_param_row, 'refClass')
                ref_class_s.text = classnames[s]

                # Disabled switchover
                _write_switchover_service_time_strategy(sub_param_row, ProcessType.DISABLED, None, 1.0)


def _write_zero_time_xml(parent: ET.Element) -> None:
    """Append a zero-valued deterministic timing distribution."""
    distr = ET.SubElement(parent, 'subParameter')
    distr.set('classPath', 'jmt.engine.random.DeterministicDistr')
    distr.set('name', 'Deterministic')

    distr_par = ET.SubElement(parent, 'subParameter')
    distr_par.set('classPath', 'jmt.engine.random.DeterministicDistrPar')
    distr_par.set('name', 'distrPar')

    t_node = ET.SubElement(distr_par, 'subParameter')
    t_node.set('classPath', 'java.lang.Double')
    t_node.set('name', 't')
    value = ET.SubElement(t_node, 'value')
    value.text = '0.0'


def _write_delayoff_distribution(parent: ET.Element, dist) -> None:
    """Append a setup/delay-off timing distribution.

    Mirrors appendDistributionXml in MATLAB saveDelayOffStrategy.m: Exp, Erlang
    and Det are exported exactly, Immediate becomes a zero time, and any other
    distribution collapses to a deterministic time at its mean.
    """
    dist_name = dist._name if hasattr(dist, '_name') else type(dist).__name__

    if dist is None or dist_name == 'Immediate':
        _write_zero_time_xml(parent)
        return

    if dist_name == 'Exp':
        distr = ET.SubElement(parent, 'subParameter')
        distr.set('classPath', 'jmt.engine.random.Exponential')
        distr.set('name', 'Exponential')

        distr_par = ET.SubElement(parent, 'subParameter')
        distr_par.set('classPath', 'jmt.engine.random.ExponentialPar')
        distr_par.set('name', 'distrPar')

        lambda_node = ET.SubElement(distr_par, 'subParameter')
        lambda_node.set('classPath', 'java.lang.Double')
        lambda_node.set('name', 'lambda')
        value = ET.SubElement(lambda_node, 'value')
        value.text = '%.12f' % dist.get_rate()
    elif dist_name == 'Erlang':
        distr = ET.SubElement(parent, 'subParameter')
        distr.set('classPath', 'jmt.engine.random.Erlang')
        distr.set('name', 'Erlang')

        distr_par = ET.SubElement(parent, 'subParameter')
        distr_par.set('classPath', 'jmt.engine.random.ErlangPar')
        distr_par.set('name', 'distrPar')

        phases = int(dist.get_number_of_phases())
        alpha_node = ET.SubElement(distr_par, 'subParameter')
        alpha_node.set('classPath', 'java.lang.Double')
        alpha_node.set('name', 'alpha')
        value = ET.SubElement(alpha_node, 'value')
        value.text = '%.12f' % (phases / dist.get_mean())

        r_node = ET.SubElement(distr_par, 'subParameter')
        r_node.set('classPath', 'java.lang.Long')
        r_node.set('name', 'r')
        value = ET.SubElement(r_node, 'value')
        value.text = '%d' % phases
    else:
        # Det, and the deterministic-at-the-mean fallback for every other law.
        distr = ET.SubElement(parent, 'subParameter')
        distr.set('classPath', 'jmt.engine.random.DeterministicDistr')
        distr.set('name', 'Deterministic')

        distr_par = ET.SubElement(parent, 'subParameter')
        distr_par.set('classPath', 'jmt.engine.random.DeterministicDistrPar')
        distr_par.set('name', 'distrPar')

        t_node = ET.SubElement(distr_par, 'subParameter')
        t_node.set('classPath', 'java.lang.Double')
        t_node.set('name', 't')
        value = ET.SubElement(t_node, 'value')
        value.text = '%.12f' % dist.get_mean()


def _write_delayoff_strategy(server_elem: ET.Element, node_idx: int, classnames: List[str],
                             model: Any) -> None:
    """Write delayOffTime and setUpTime strategies for a Queue with delay-off.

    Port of MATLAB saveDelayOffStrategy.m. The queue node object carries the
    setup/delay-off distributions, which do not reach sn.nodeparam in Python.
    """
    if model is None:
        return

    nodes = model.get_nodes() if hasattr(model, 'get_nodes') else getattr(model, 'nodes', None)
    if not nodes or node_idx >= len(nodes):
        return

    node = nodes[node_idx]
    if not hasattr(node, 'is_delay_off_enabled') or not node.is_delay_off_enabled():
        return

    jobclasses = model.get_classes() if hasattr(model, 'get_classes') else getattr(model, 'classes', [])

    for name, getter in (('delayOffTime', 'get_delay_off_time'),
                         ('setUpTime', 'get_setup_time')):
        param_node = ET.SubElement(server_elem, 'parameter')
        param_node.set('array', 'true')
        param_node.set('classPath', 'java.lang.Object')
        param_node.set('name', name)

        for r in range(len(classnames)):
            ref_class = ET.SubElement(param_node, 'refClass')
            ref_class.text = classnames[r]

            # JMT reads Server.java:1280 element [0] of a ServiceStrategy[];
            # the strategy must be wrapped in a single-element array.
            class_row = ET.SubElement(param_node, 'subParameter')
            class_row.set('array', 'true')
            class_row.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategy')
            class_row.set('name', name)

            sub_param = ET.SubElement(class_row, 'subParameter')
            sub_param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy')
            sub_param.set('name', 'ServiceTimeStrategy')

            dist = getattr(node, getter)(jobclasses[r]) if r < len(jobclasses) else None
            if dist is None:
                _write_zero_time_xml(sub_param)
            else:
                _write_delayoff_distribution(sub_param, dist)


def _write_balking_strategy(balking_node: ET.Element, thresholds, nservers) -> None:
    """Populate a JMT Balking element from LINE queue-length balking thresholds.

    Port of saveBalkingStrategy in MATLAB saveImpatience.m. LINE stores a list
    of {minJobs, maxJobs, probability} closed intervals; JMT's Balking reads a
    LoadDependentStrategy whose LDParameter ranges are {from -> probability},
    selecting the last range with from <= queueLength (default 0). The closed
    intervals are translated into from-based breakpoints, with explicit
    0-probability breakpoints at gaps so that queue lengths outside any interval
    do not inherit a neighbour's probability.

    Queue-length convention: LINE evaluates balking against the TOTAL station
    population (in-service + waiting), as do the CTMC/SSA/LDES solvers, whereas
    JMT's Balking evaluates against the number WAITING only. When the servers
    are busy (the regime where balking is meaningful) the two differ by exactly
    the server count S, so every `from` is shifted down by S: JMT waiting w maps
    to LINE total n = w + S.
    """
    entries = []
    for th in (thresholds or []):
        if th is None:
            continue
        entries.append((float(th[0]), float(th[1]), float(th[2])))
    entries.sort(key=lambda e: e[0])

    S = 1
    if nservers is not None and np.isfinite(nservers) and nservers >= 1:
        S = int(nservers)

    froms = []
    probs = []
    for ti, (lo, hi, pr) in enumerate(entries):
        froms.append(max(0, int(lo) - S))
        probs.append(pr)
        if np.isfinite(hi):
            next_lo = entries[ti + 1][0] if ti + 1 < len(entries) else float('inf')
            if hi + 1 < next_lo:
                froms.append(max(0, int(hi) + 1 - S))
                probs.append(0.0)

    # LoadDependentStrategy wrapper (Parameter #1 of Balking)
    ld_strategy = ET.SubElement(balking_node, 'subParameter')
    ld_strategy.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.LoadDependentStrategy')
    ld_strategy.set('name', 'LoadDependentStrategy')

    ld_array = ET.SubElement(ld_strategy, 'subParameter')
    ld_array.set('array', 'true')
    ld_array.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.LDParameter')
    ld_array.set('name', 'LDParameter')

    for from_val, prob in zip(froms, probs):
        range_node = ET.SubElement(ld_array, 'subParameter')
        range_node.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.LDParameter')
        range_node.set('name', 'LDParameter')

        from_node = ET.SubElement(range_node, 'subParameter')
        from_node.set('classPath', 'java.lang.Integer')
        from_node.set('name', 'from')
        ET.SubElement(from_node, 'value').text = '%d' % from_val

        # Dummy distribution; only the function/probability is read for balking.
        distr = ET.SubElement(range_node, 'subParameter')
        distr.set('classPath', 'jmt.engine.random.Exponential')
        distr.set('name', 'Exponential')

        distr_par = ET.SubElement(range_node, 'subParameter')
        distr_par.set('classPath', 'jmt.engine.random.ExponentialPar')
        distr_par.set('name', 'distrPar')
        lambda_node = ET.SubElement(distr_par, 'subParameter')
        lambda_node.set('classPath', 'java.lang.Double')
        lambda_node.set('name', 'lambda')
        ET.SubElement(lambda_node, 'value').text = '1.0'

        func_node = ET.SubElement(range_node, 'subParameter')
        func_node.set('classPath', 'java.lang.String')
        func_node.set('name', 'function')
        ET.SubElement(func_node, 'value').text = '%.12f' % prob

    # priorityActivated flag (Parameter #2 of Balking)
    prio = ET.SubElement(balking_node, 'subParameter')
    prio.set('classPath', 'java.lang.Boolean')
    prio.set('name', 'priorityActivated')
    ET.SubElement(prio, 'value').text = 'false'


def _has_queue_length_balking(sn: NetworkStruct, ist: int, r: int) -> bool:
    """True if station ist balks class r on queue length."""
    strategy = getattr(sn, 'balkingStrategy', None)
    if strategy is None or ist < 0:
        return False
    strategy = np.asarray(strategy)
    if ist >= strategy.shape[0] or r >= strategy.shape[1]:
        return False
    from ....lang.base import BalkingStrategy
    return int(strategy[ist, r]) == int(BalkingStrategy.QUEUE_LENGTH)


def _matlab_num2str(v: float) -> str:
    """MATLAB num2str() of a real scalar.

    The service weights are the one JMT field MATLAB writes with num2str
    (savePreemptiveWeights.m) rather than %.12f, and num2str spells an
    integer-valued double without a decimal point and everything else with
    floor(log10(|v|))+5 significant digits. Writing repr(1.0) = '1.0' where
    MATLAB writes '1' is a document difference for no reason.
    """
    v = float(v)
    if not np.isfinite(v):
        return str(v)
    if v == int(v):
        return '%d' % int(v)
    digits = int(np.floor(np.log10(abs(v)))) + 5
    if digits < 1:
        digits = 1
    return ('%.' + str(digits) + 'g') % v


_DROP_STRATEGY_TEXT = {
    -1: 'waiting queue',   # WAITQ
    1: 'drop',             # DROP
    2: 'BAS blocking',     # BAS
    3: 'BBS blocking',     # BBS
    4: 'RSRD blocking',    # RSRD
    5: 'retrial',          # RETRIAL
    6: 'retrial with limit',
}


def _jmt_is_bas_destination(sn: NetworkStruct, ist: int, r: int) -> bool:
    """True when station ist is the RECEIVING side of a true-BAS relation for class r.

    That is, an arrival of r that finds ist full must block an upstream station
    rather than be lost. LINE accepts the BAS declaration in two places -- on the
    blocking (upstream) station, as cqn_bas_blocking does, or on the full
    destination, as a model read back from JMT does -- and ctmc_ssg resolves both
    into sn.isbasdestination (BUG-83). Reading sn.droprule at the capped station
    sees only the second form, which is what made SolverJMT refuse the first one.
    """
    isbd = getattr(sn, 'isbasdestination', None)
    if isbd is None or ist is None or ist < 0:
        return False
    isbd = np.atleast_2d(np.asarray(isbd))
    if ist >= isbd.shape[0] or r >= isbd.shape[1]:
        return False
    return bool(isbd[ist, r])


#: The only dropStrategy ids JMT's queue section reads: DROP, BAS, WAITQ, RETRIAL.
#: It matches them with a lookupswitch on String.hashCode in
#: jmt/engine/NodeSections/Queue.class (the Storage section of a Place is even
#: narrower and drops 'retrial'); an unrecognized value falls through the default
#: arm with NO flag set, so BBS, RSRD and retrial-with-limit are not approximated,
#: they are IGNORED.
_JMT_READABLE_DROP_IDS = frozenset((1, 2, -1, 5))


def _drop_strategy_text(sn: NetworkStruct, ist: int, r: int) -> str:
    """JMT dropStrategy string for station ist, class r.

    Mirrors MATLAB jmtDropStrategyText.m: an unset slot (no station, NaN, or 0)
    is the bufferless-node field 'drop'; otherwise the DropStrategy id is spelled
    out. The ids are shared across the codebases (see lang/base.py DropStrategy),
    so they are read here as numbers.

    Beyond the id table this resolves the two ways LINE can declare BAS blocking
    onto the one way JMT can read it. JMT's queue section says what happens to an
    arrival that finds THIS buffer full, so it only understands the rule on the
    destination; a WAITQ slot that _jmt_is_bas_destination marks is therefore
    written out as 'BAS blocking'.

    It also keeps the written file VALID: a strategy outside
    _JMT_READABLE_DROP_IDS is spelled 'waiting queue', JMT's own no-limit
    default. That substitution is only ever reached where the rule cannot be
    consulted (infinite size, or a closed capacity equal to the population): a
    buffer that can actually fill under one of those is refused outright by
    _jmt_station_cap_assert.
    """
    droprule = getattr(sn, 'droprule', None)
    if droprule is None or ist is None or ist < 0:
        return 'drop'
    droprule = np.asarray(droprule)
    if droprule.ndim < 2 or ist >= droprule.shape[0] or r >= droprule.shape[1]:
        return 'drop'
    drop_val = droprule[ist, r]
    if np.isnan(drop_val) or int(drop_val) == 0:
        return 'drop'
    if int(drop_val) == -1 and _jmt_is_bas_destination(sn, ist, r):
        return _DROP_STRATEGY_TEXT[2]   # WAITQ slot standing in for upstream-declared BAS
    if int(drop_val) not in _JMT_READABLE_DROP_IDS:
        if int(drop_val) not in _DROP_STRATEGY_TEXT:
            raise ValueError('Unrecognized drop strategy type: %s' % drop_val)
        return _DROP_STRATEGY_TEXT[-1]  # unreachable rule, written as JMT's no-limit default
    return _DROP_STRATEGY_TEXT[int(drop_val)]


def _balking_thresholds(sn: NetworkStruct, ist: int, r: int):
    thresholds = getattr(sn, 'balkingThresholds', None)
    if thresholds is None:
        return None
    try:
        return thresholds[ist][r]
    except (IndexError, TypeError, KeyError):
        return None


def _detect_class_switches(sn: NetworkStruct) -> Dict[Tuple[int, int], np.ndarray]:
    """
    Detect node pairs that require ClassSwitch nodes for JMT.

    For each pair of connected nodes (i, j), check if there's class switching
    in the routing. If the class transition matrix from i to j is not diagonal,
    a ClassSwitch node needs to be inserted.

    Args:
        sn: NetworkStruct object

    Returns:
        Dictionary mapping (source_idx, dest_idx) to class switching matrix (K x K)
        where matrix[r][s] is the probability of switching from class r to class s
        when routing from source to dest.
    """
    if sn.rtnodes is None or sn.connmatrix is None:
        return {}

    K = sn.nclasses
    nnodes = sn.nnodes
    cs_nodes = {}

    # For each connected node pair, check if there's class switching
    for i in range(nnodes):
        # Skip if source node is already a ClassSwitch - it already handles class transitions
        if sn.nodetype is not None and len(sn.nodetype) > i:
            if sn.nodetype[i] == NodeType.CLASSSWITCH:
                continue

        for j in range(nnodes):
            if sn.connmatrix[i, j] <= 0:
                continue

            # Build the class switching matrix for this node pair
            cs_matrix = np.zeros((K, K))
            has_nonzero = False

            for r in range(K):
                row_sum = 0.0
                for s in range(K):
                    src_idx = i * K + r
                    dst_idx = j * K + s
                    if src_idx < sn.rtnodes.shape[0] and dst_idx < sn.rtnodes.shape[1]:
                        prob = sn.rtnodes[src_idx, dst_idx]
                        cs_matrix[r, s] = prob
                        row_sum += prob
                        if prob > 0:
                            has_nonzero = True

                # Normalize row if it has non-zero entries
                if row_sum > 0:
                    cs_matrix[r, :] /= row_sum

            if not has_nonzero:
                continue

            # Check if matrix is not diagonal (i.e., has class switching)
            is_diagonal = True
            for r in range(K):
                for s in range(K):
                    if r != s and cs_matrix[r, s] > 1e-10:
                        is_diagonal = False
                        break
                if not is_diagonal:
                    break

            if not is_diagonal:
                cs_nodes[(i, j)] = cs_matrix

    return cs_nodes


def _write_auto_classswitch_node(
    sim: ET.Element,
    cs_name: str,
    cs_matrix: np.ndarray,
    dest_node: str,
    classnames: List[str]
) -> None:
    """
    Write an auto-generated ClassSwitch node for handling class switching in routing.

    Args:
        sim: Parent XML element (sim)
        cs_name: Name for the ClassSwitch node (e.g., "CS_Queue1_to_Delay")
        cs_matrix: K x K class switching probability matrix
        dest_node: Name of destination node to route to
        classnames: List of class names
    """
    K = len(classnames)

    node_elem = ET.SubElement(sim, 'node')
    node_elem.set('name', cs_name)

    # 1. Queue section (input buffer)
    queue = ET.SubElement(node_elem, 'section')
    queue.set('className', 'Queue')

    size_param = ET.SubElement(queue, 'parameter')
    size_param.set('classPath', 'java.lang.Integer')
    size_param.set('name', 'size')
    value = ET.SubElement(size_param, 'value')
    value.text = '-1'

    # Drop strategies
    drop_strategy = ET.SubElement(queue, 'parameter')
    drop_strategy.set('array', 'true')
    drop_strategy.set('classPath', 'java.lang.String')
    drop_strategy.set('name', 'dropStrategies')

    for r in range(K):
        ref_class = ET.SubElement(drop_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(drop_strategy, 'subParameter')
        sub_param.set('classPath', 'java.lang.String')
        sub_param.set('name', 'dropStrategy')
        value = ET.SubElement(sub_param, 'value')
        value.text = 'drop'

    # Queue get strategy (FCFS)
    strategy_param = ET.SubElement(queue, 'parameter')
    strategy_param.set('classPath', 'jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy')
    strategy_param.set('name', 'FCFSstrategy')

    # Queue put strategy
    put_strategy = ET.SubElement(queue, 'parameter')
    put_strategy.set('array', 'true')
    put_strategy.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategy')
    put_strategy.set('name', 'QueuePutStrategy')

    for r in range(K):
        ref_class = ET.SubElement(put_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(put_strategy, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy')
        sub_param.set('name', 'TailStrategy')

    # 2. ClassSwitch section
    cs_section = ET.SubElement(node_elem, 'section')
    cs_section.set('className', 'ClassSwitch')

    matrix_param = ET.SubElement(cs_section, 'parameter')
    matrix_param.set('array', 'true')
    matrix_param.set('classPath', 'java.lang.Object')
    matrix_param.set('name', 'matrix')

    for r in range(K):
        ref_class = ET.SubElement(matrix_param, 'refClass')
        ref_class.text = classnames[r]

        row_param = ET.SubElement(matrix_param, 'subParameter')
        row_param.set('array', 'true')
        row_param.set('classPath', 'java.lang.Float')
        row_param.set('name', 'row')

        for s in range(K):
            ref_class_col = ET.SubElement(row_param, 'refClass')
            ref_class_col.text = classnames[s]

            cell_param = ET.SubElement(row_param, 'subParameter')
            cell_param.set('classPath', 'java.lang.Float')
            cell_param.set('name', 'cell')
            cell_value = ET.SubElement(cell_param, 'value')
            cell_value.text = f'{cs_matrix[r, s]:.12f}'

    # 3. Router section - always route to destination with Random strategy
    router = ET.SubElement(node_elem, 'section')
    router.set('className', 'Router')

    routing_param = ET.SubElement(router, 'parameter')
    routing_param.set('array', 'true')
    routing_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategy')
    routing_param.set('name', 'RoutingStrategy')

    for r in range(K):
        ref_class = ET.SubElement(routing_param, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(routing_param, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.RandomStrategy')
        sub_param.set('name', 'Random')


def _write_jsim_file(sn: NetworkStruct, model_path: str, options: SolverJMTOptions, model: Any = None) -> None:
    """
    Write the network model to JSIM XML format.

    This is a simplified version supporting basic queueing networks.

    Args:
        sn: NetworkStruct containing network configuration
        model_path: Path to write the JSIM file
        options: Solver options
        model: Optional Network model (for FCR regions)
    """
    M = sn.nstations
    K = sn.nclasses

    # BUG-83: sn.isbasdestination is derived, and unlike MATLAB/JAR/C++ the
    # native struct refresh does not build it, so it is populated here exactly
    # as ssa/serial.py does. _drop_strategy_text and _jmt_station_cap_assert
    # both read it to see BAS declared on the UPSTREAM station.
    from ...state.ctmc_ssg import _populate_isbasblocking
    _populate_isbasblocking(sn)

    from datetime import datetime
    timestamp = datetime.now().strftime('%a %b %d %H:%M:%S %Y')

    # Create sim as root element (matches MATLAB's writeJSIM format)
    # MATLAB uses: simDoc = com.mathworks.xml.XMLUtils.createDocument('sim')
    sim = ET.Element('sim')
    sim.set('xmlns:xsi', 'http://www.w3.org/2001/XMLSchema-instance')
    sim.set('disableStatisticStop', 'true')
    sim.set('logDecimalSeparator', '.')
    sim.set('logDelimiter', ';')

    # Determine logPath: use Logger's filePath if available, otherwise empty string
    # MATLAB/JAR use empty logPath by default, which is important for simulation consistency
    log_path = ''
    if sn.nodeparam is not None:
        for node_idx in sn.nodeparam:
            param = sn.nodeparam[node_idx]
            if hasattr(param, 'filePath') and param.filePath:
                # Use the Logger's filePath as the sim's logPath
                log_path = param.filePath
                break
    sim.set('logPath', log_path)
    sim.set('logReplaceMode', '0')
    sim.set('maxEvents', '-1')
    sim.set('maxSamples', str(options.samples))
    # Only set maxSimulated if it's finite (matches JAR/MATLAB behavior)
    # JMT interprets absence of maxSimulated as unlimited simulation time
    if not np.isinf(options.max_simulated_time):
        sim.set('maxSimulated', f'{options.max_simulated_time:.3f}')
    sim.set('name', os.path.basename(model_path))
    sim.set('polling', '1.0')
    sim.set('seed', str(options.seed))
    sim.set('xsi:noNamespaceSchemaLocation', 'SIMmodeldefinition.xsd')

    # Create class definitions
    njobs = sn.njobs.flatten() if sn.njobs is not None else np.zeros(K)
    classnames = sn.classnames if sn.classnames else [f'Class{i+1}' for i in range(K)]
    nodenames = sn.nodenames if sn.nodenames else [f'Node{i+1}' for i in range(sn.nnodes)]
    refstat = sn.refstat.flatten() if hasattr(sn, 'refstat') and sn.refstat is not None else np.zeros(K, dtype=int)

    # JMT uses higher priority value = higher priority, LINE uses lower value = higher priority
    # We need to invert priorities when exporting to JMT
    max_prio = 0
    if hasattr(sn, 'classprio') and sn.classprio is not None:
        max_prio = int(np.max(sn.classprio))

    for r in range(K):
        class_elem = ET.SubElement(sim, 'userClass')
        class_elem.set('name', classnames[r])

        # Set priority with inversion: LINE uses lower=higher, JMT uses higher=higher
        if hasattr(sn, 'classprio') and sn.classprio is not None:
            classprio_flat = sn.classprio.flatten() if sn.classprio.ndim > 1 else sn.classprio
            line_prio = int(classprio_flat[r])
            jmt_prio = max_prio - line_prio
            class_elem.set('priority', str(jmt_prio))
        else:
            class_elem.set('priority', '0')

        if np.isinf(njobs[r]):
            # Open class - reference source is Source node
            # Find the actual source node name for this class
            source_name = 'Source'  # Default fallback
            if hasattr(sn, 'nodetype') and sn.nodetype is not None:
                for node_idx in range(len(sn.nodetype)):
                    if int(sn.nodetype[node_idx]) == NodeType.SOURCE:
                        source_name = nodenames[node_idx] if node_idx < len(nodenames) else 'Source'
                        break
            class_elem.set('referenceSource', source_name)
            class_elem.set('softDeadline', '0.0')
            class_elem.set('type', 'open')
        else:
            # Closed class - reference source is the station where jobs start
            # refstat contains station indices, convert to node index using stationToNode
            ref_station = int(refstat[r])
            if hasattr(sn, 'stationToNode') and sn.stationToNode is not None and ref_station < len(sn.stationToNode):
                ref_node_idx = int(sn.stationToNode[ref_station])
            else:
                ref_node_idx = ref_station
            ref_name = nodenames[ref_node_idx] if ref_node_idx < len(nodenames) else classnames[r] + '_RefStation'
            class_elem.set('referenceSource', ref_name)
            class_elem.set('type', 'closed')
            class_elem.set('customers', str(int(njobs[r])))
            class_elem.set('softDeadline', '0.0')

    # Detect class switches in routing - JMT requires ClassSwitch nodes for class transitions
    cs_nodes = _detect_class_switches(sn)
    # Map (source, dest) -> ClassSwitch node name
    cs_node_names = {}
    for (src_idx, dst_idx), cs_matrix in cs_nodes.items():
        cs_name = f"CS_{nodenames[src_idx]}_to_{nodenames[dst_idx]}"
        cs_node_names[(src_idx, dst_idx)] = cs_name

    # Create nodes
    nodenames = sn.nodenames if sn.nodenames else [f'Node{i+1}' for i in range(sn.nnodes)]

    for i in range(sn.nnodes):
        node_type = sn.nodetype[i] if sn.nodetype is not None and len(sn.nodetype) > i else NodeType.QUEUE
        node_name = nodenames[i]

        node_elem = ET.SubElement(sim, 'node')
        node_elem.set('name', node_name)

        if node_type == NodeType.SOURCE:
            _write_source_node(node_elem, i, sn, classnames, cs_node_names)
        elif node_type == NodeType.SINK:
            _write_sink_node(node_elem, sn, classnames)
        elif node_type == NodeType.DELAY:
            _write_delay_node(node_elem, i, sn, classnames, cs_node_names)
        elif node_type == NodeType.QUEUE:
            # Check if this is actually a Delay (infinite server) queue by checking SchedStrategy.INF
            ist = int(sn.nodeToStation[i]) if hasattr(sn, 'nodeToStation') and sn.nodeToStation is not None and i < len(sn.nodeToStation) else i
            sched = sn.sched.get(ist, SchedStrategy.FCFS) if sn.sched else SchedStrategy.FCFS
            if sched == SchedStrategy.INF:
                _write_delay_node(node_elem, i, sn, classnames, cs_node_names)
            else:
                _write_queue_node(node_elem, i, sn, classnames, options, cs_node_names, model)
        elif node_type == NodeType.ROUTER:
            _write_router_node(node_elem, i, sn, classnames, cs_node_names)
        elif node_type == NodeType.FORK:
            _write_fork_node(node_elem, i, sn, classnames, cs_node_names)
        elif node_type == NodeType.JOIN:
            _write_join_node(node_elem, i, sn, classnames, cs_node_names)
        elif node_type == NodeType.CLASSSWITCH:
            _write_classswitch_node(node_elem, i, sn, classnames)
        elif node_type == NodeType.PLACE:
            _write_place_node(node_elem, i, sn, classnames)
        elif node_type == NodeType.TRANSITION:
            _write_transition_node(node_elem, i, sn, classnames)
        elif node_type == NodeType.LOGGER:
            _write_logger_node(node_elem, i, sn, classnames, cs_node_names)

    # Create auto-generated ClassSwitch nodes for class switching in routing
    for (src_idx, dst_idx), cs_matrix in cs_nodes.items():
        cs_name = cs_node_names[(src_idx, dst_idx)]
        dest_name = nodenames[dst_idx]
        _write_auto_classswitch_node(sim, cs_name, cs_matrix, dest_name, classnames)

    # Metrics must precede connections per JMT schema; group by metric type
    # (QLen, Util, RespT, Tput, ArvR, Tard), matching JAR/MATLAB order.
    # Q/U/R disabled for Source/Sink; T/A enabled for Source only.
    alpha_str = f'{round(1 - options.conf_int, 10)}'

    # Helper to get station info
    def _get_station_info(ist):
        node_idx = int(sn.stationToNode[ist]) if sn.stationToNode is not None else ist
        node_name = nodenames[node_idx]
        node_type = sn.nodetype[node_idx] if sn.nodetype is not None and len(sn.nodetype) > node_idx else NodeType.QUEUE
        is_source = (node_type == NodeType.SOURCE)
        is_sink = (node_type == NodeType.SINK)
        is_fork = (node_type == NodeType.FORK)
        is_join = (node_type == NodeType.JOIN)
        return node_name, is_source, is_sink, is_fork, is_join

    # 1. Queue Length (Number of Customers) - skip Source and Sink
    for i in range(M):
        node_name, is_source, is_sink, _, _ = _get_station_info(i)
        if is_source or is_sink:
            continue
        for r in range(K):
            metric = ET.SubElement(sim, 'measure')
            metric.set('alpha', alpha_str)
            metric.set('name', f'Performance_{i+1}')
            metric.set('nodeType', 'station')
            metric.set('precision', str(options.max_rel_err))
            metric.set('referenceNode', node_name)
            metric.set('referenceUserClass', classnames[r])
            metric.set('type', 'Number of Customers')
            metric.set('verbose', 'false')

    # 2. Utilization - skip Source, Sink, Fork, Join
    for i in range(M):
        node_name, is_source, is_sink, is_fork, is_join = _get_station_info(i)
        if is_source or is_sink or is_fork or is_join:
            continue
        for r in range(K):
            metric = ET.SubElement(sim, 'measure')
            metric.set('alpha', alpha_str)
            metric.set('name', f'Performance_{i+1}')
            metric.set('nodeType', 'station')
            metric.set('precision', str(options.max_rel_err))
            metric.set('referenceNode', node_name)
            metric.set('referenceUserClass', classnames[r])
            metric.set('type', 'Utilization')
            metric.set('verbose', 'false')

    # 3. Response Time - skip Source and Sink
    for i in range(M):
        node_name, is_source, is_sink, _, _ = _get_station_info(i)
        if is_source or is_sink:
            continue
        for r in range(K):
            metric = ET.SubElement(sim, 'measure')
            metric.set('alpha', alpha_str)
            metric.set('name', f'Performance_{i+1}')
            metric.set('nodeType', 'station')
            metric.set('precision', str(options.max_rel_err))
            metric.set('referenceNode', node_name)
            metric.set('referenceUserClass', classnames[r])
            metric.set('type', 'Response Time')
            metric.set('verbose', 'false')

    # 4. Throughput - all stations including Source, skip Sink
    for i in range(M):
        node_name, _, is_sink, _, _ = _get_station_info(i)
        if is_sink:
            continue
        for r in range(K):
            metric = ET.SubElement(sim, 'measure')
            metric.set('alpha', alpha_str)
            metric.set('name', f'Performance_{i+1}')
            metric.set('nodeType', 'station')
            metric.set('precision', str(options.max_rel_err))
            metric.set('referenceNode', node_name)
            metric.set('referenceUserClass', classnames[r])
            metric.set('type', 'Throughput')
            metric.set('verbose', 'false')

    # 5. Arrival Rate - all stations including Source, skip Sink
    for i in range(M):
        node_name, _, is_sink, _, _ = _get_station_info(i)
        if is_sink:
            continue
        for r in range(K):
            metric = ET.SubElement(sim, 'measure')
            metric.set('alpha', alpha_str)
            metric.set('name', f'Performance_{i+1}')
            metric.set('nodeType', 'station')
            metric.set('precision', str(options.max_rel_err))
            metric.set('referenceNode', node_name)
            metric.set('referenceUserClass', classnames[r])
            metric.set('type', 'Arrival Rate')
            metric.set('verbose', 'false')

    # 6. Tardiness - skip Source and Sink (only if classes have deadlines)
    for i in range(M):
        node_name, is_source, is_sink, _, _ = _get_station_info(i)
        if is_source or is_sink:
            continue
        for r in range(K):
            metric = ET.SubElement(sim, 'measure')
            metric.set('alpha', alpha_str)
            metric.set('name', f'Performance_{i+1}')
            metric.set('nodeType', 'station')
            metric.set('precision', str(options.max_rel_err))
            metric.set('referenceNode', node_name)
            metric.set('referenceUserClass', classnames[r])
            metric.set('type', 'Tardiness')
            metric.set('verbose', 'false')

    # 7. System Tardiness - one per class, no station. MATLAB saveMetrics.m
    # emits the SysTard handles right after the per-station Tard ones, and
    # saveMetric.m blanks nodeType/referenceNode for a system-level index. The
    # measure list is what JMT tests its precision target against, so a missing
    # one is not merely a missing column.
    for r in range(K):
        metric = ET.SubElement(sim, 'measure')
        metric.set('alpha', alpha_str)
        metric.set('name', 'Performance_1')
        metric.set('nodeType', '')
        metric.set('precision', str(options.max_rel_err))
        metric.set('referenceNode', '')
        metric.set('referenceUserClass', classnames[r])
        metric.set('type', 'System Tardiness')
        metric.set('verbose', 'false')

    # FCR metrics: region-aggregate QLen/RespT/ResidT/Tput (no per-class
    # breakdown), named "FCRegion{n}" to match saveRegions; mirrors JAR/MATLAB.
    if model is not None:
        _fcr_regions = []
        if hasattr(model, 'get_regions'):
            _fcr_regions = model.get_regions()
        elif hasattr(model, 'regions'):
            _fcr_regions = model.regions
        for r_idx in range(len(_fcr_regions)):
            fcr_name = f'FCRegion{r_idx + 1}'
            for mtype in ('Number of Customers', 'Response Time', 'Residence Time', 'Throughput'):
                metric = ET.SubElement(sim, 'measure')
                metric.set('alpha', alpha_str)
                metric.set('name', f'FCR_{fcr_name}_{mtype.replace(" ", "")}')
                metric.set('nodeType', 'region')
                metric.set('precision', str(options.max_rel_err))
                metric.set('referenceNode', fcr_name)
                metric.set('referenceUserClass', '')  # FCR metrics are not class-specific
                metric.set('type', mtype)
                metric.set('verbose', 'false')

    # Connections must come after metrics per JMT schema. Class switching
    # routes i -> CS_i_to_j -> j instead of i -> j. Column-major order (j
    # outer, i inner) matches MATLAB's find() behavior.
    if sn.connmatrix is not None:
        for j in range(sn.nnodes):  # columns first (like MATLAB find)
            for i in range(sn.nnodes):  # rows
                if sn.connmatrix[i, j] > 0:
                    if (i, j) in cs_node_names:
                        # Route through ClassSwitch node
                        cs_name = cs_node_names[(i, j)]
                        # Connection: source -> ClassSwitch
                        conn = ET.SubElement(sim, 'connection')
                        conn.set('source', nodenames[i])
                        conn.set('target', cs_name)
                        # Connection: ClassSwitch -> dest
                        conn = ET.SubElement(sim, 'connection')
                        conn.set('source', cs_name)
                        conn.set('target', nodenames[j])
                    else:
                        # Direct connection (no class switching)
                        conn = ET.SubElement(sim, 'connection')
                        conn.set('source', nodenames[i])
                        conn.set('target', nodenames[j])

    # Add blocking regions (FCR - Finite Capacity Regions)
    # Reference: MATLAB saveRegions.m
    if model is not None:
        regions = []
        if hasattr(model, 'get_regions'):
            regions = model.get_regions()
        elif hasattr(model, 'regions'):
            regions = model.regions

        for r_idx, region in enumerate(regions):
            blocking_region = ET.SubElement(sim, 'blockingRegion')
            # JMT-internal region name must match the FCR measure referenceNode
            # above, independent of the LINE display name used in the node table.
            region_name = f'FCRegion{r_idx + 1}'
            blocking_region.set('name', region_name)
            blocking_region.set('type', 'default')

            # 1. regionNode elements - nodes in this region
            region_nodes = region.nodes if hasattr(region, 'nodes') else []
            for node in region_nodes:
                node_name = node.get_name() if hasattr(node, 'get_name') else str(node)
                region_node = ET.SubElement(blocking_region, 'regionNode')
                region_node.set('nodeName', node_name)

            # 2. globalConstraint
            global_constraint = ET.SubElement(blocking_region, 'globalConstraint')
            global_max = region.global_max_jobs if hasattr(region, 'global_max_jobs') else -1
            global_constraint.set('maxJobs', str(global_max))

            # 3. globalMemoryConstraint
            global_mem_constraint = ET.SubElement(blocking_region, 'globalMemoryConstraint')
            global_max_mem = region.global_max_memory if hasattr(region, 'global_max_memory') else -1
            global_mem_constraint.set('maxMemory', str(global_max_mem))

            # Get classes from region or sn
            region_classes = region.classes if hasattr(region, 'classes') else []

            # 4. classConstraint elements
            for job_class in region_classes:
                class_name = job_class.get_name() if hasattr(job_class, 'get_name') else str(job_class)
                class_max_jobs = region.get_class_max_jobs(job_class) if hasattr(region, 'get_class_max_jobs') else -1

                # Only write if not unbounded (-1)
                if class_max_jobs != -1:
                    class_constraint = ET.SubElement(blocking_region, 'classConstraint')
                    class_constraint.set('jobClass', class_name)
                    class_constraint.set('maxJobsPerClass', str(class_max_jobs))

            # 5. classMemoryConstraint elements
            for job_class in region_classes:
                class_name = job_class.get_name() if hasattr(job_class, 'get_name') else str(job_class)
                class_max_mem = region.get_class_max_memory(job_class) if hasattr(region, 'get_class_max_memory') else -1

                # Only write if not unbounded (-1)
                if class_max_mem != -1:
                    class_mem_constraint = ET.SubElement(blocking_region, 'classMemoryConstraint')
                    class_mem_constraint.set('jobClass', class_name)
                    class_mem_constraint.set('maxMemoryPerClass', str(class_max_mem))

            # 6. dropRules elements - always write for each class
            for job_class in region_classes:
                class_name = job_class.get_name() if hasattr(job_class, 'get_name') else str(job_class)
                drop_rule = region.get_drop_rule(job_class) if hasattr(region, 'get_drop_rule') else None

                drop_rules = ET.SubElement(blocking_region, 'dropRules')
                drop_rules.set('jobClass', class_name)

                # Determine if DROP or WAITQ
                if drop_rule is not None:
                    # Check if it's DROP strategy
                    is_drop = False
                    if hasattr(drop_rule, 'name'):
                        is_drop = drop_rule.name == 'DROP'
                    elif hasattr(drop_rule, 'value'):
                        is_drop = drop_rule.value == 1  # DROP = 1
                    elif drop_rule is True:
                        is_drop = True
                    drop_rules.set('dropThisClass', 'true' if is_drop else 'false')
                else:
                    drop_rules.set('dropThisClass', 'false')

            # 7. classSize elements (only if not default value of 1)
            for job_class in region_classes:
                class_name = job_class.get_name() if hasattr(job_class, 'get_name') else str(job_class)
                class_size = region.get_class_size(job_class) if hasattr(region, 'get_class_size') else 1

                if class_size != 1:
                    class_size_elem = ET.SubElement(blocking_region, 'classSize')
                    class_size_elem.set('jobClass', class_name)
                    class_size_elem.set('size', str(class_size))

    # Preload section (MATLAB writeJSIM.m:145-185): closed-network jobs
    # start at their reference stations; SPN Places need initial token populations.
    njobs = sn.njobs.flatten() if sn.njobs is not None else np.zeros(K)

    # Get initial state if available (for SPNs with Places)
    s0 = sn.state if hasattr(sn, 'state') and sn.state is not None else None
    stationToStateful = sn.stationToStateful if hasattr(sn, 'stationToStateful') else None
    # Flatten stationToStateful to ensure proper indexing (it may be 2D)
    if stationToStateful is not None:
        stationToStateful = np.asarray(stationToStateful).flatten()

    # Check if we need preload section
    has_reference_nodes = False
    preload = ET.SubElement(sim, 'preload')

    # For each station (excluding Source and Join nodes)
    for ist in range(M):
        node_idx = int(sn.stationToNode[ist]) if sn.stationToNode is not None else ist
        node_type = int(sn.nodetype[node_idx]) if sn.nodetype is not None else -1

        # Skip Source (0) and Join (5) nodes
        if node_type == 0 or node_type == 5:
            continue

        node_name = nodenames[node_idx] if node_idx < len(nodenames) else f'Node{node_idx}'

        # Get initial population from state (like MATLAB's State.toMarginal)
        # For Places in SPNs, this gives the initial token counts
        nir = np.zeros(K)
        has_explicit_state = False
        if s0 is not None and stationToStateful is not None:
            stateful_idx = int(stationToStateful[ist]) if ist < len(stationToStateful) else -1
            if stateful_idx >= 0 and stateful_idx < len(s0):
                # State.toMarginal as in writeJSIM: the raw state is per-class-per-phase, so its first K entries are not the per-class counts.
                from ...state.marginal import toMarginal
                _, nir_states, _, _ = toMarginal(
                    sn, node_idx, np.asarray(s0[stateful_idx]))
                nir = np.asarray(nir_states)[0, :K].astype(float)
                has_explicit_state = True

        # For closed classes, use reference station logic ONLY if state was not explicitly set
        # If init_from_marginal was called, we trust the state values (even if 0)
        if not has_explicit_state:
            for r in range(K):
                if np.isfinite(njobs[r]) and njobs[r] > 0 and nir[r] == 0:
                    # Check if this station is the reference station for class r
                    ref_idx = int(refstat[r]) if r < len(refstat) else 0
                    if ist == ref_idx:
                        # All jobs of this class start at reference station
                        nir[r] = njobs[r]

        # JAR includes all stations (Queue, Delay) in preload, even at
        # population=0, for proper JMT initialization.
        class_populations = []
        for r in range(K):
            # Open classes: njobs[r]=Inf, nir[r]=0. Closed: njobs[r] finite,
            # nir[r] = jobs at this station.
            class_populations.append((classnames[r], int(round(nir[r]))))

        # Include station in preload - JAR includes all queue stations
        if class_populations:
            has_reference_nodes = True
            station_pop = ET.SubElement(preload, 'stationPopulations')
            station_pop.set('stationName', node_name)

            for class_name, pop in class_populations:
                class_pop = ET.SubElement(station_pop, 'classPopulation')
                class_pop.set('population', str(pop))
                class_pop.set('refClass', class_name)

    # Only keep preload section if we have reference nodes
    if not has_reference_nodes:
        sim.remove(preload)

    # Write to file
    xml_str = ET.tostring(sim, encoding='unicode')
    dom = minidom.parseString(xml_str)
    pretty_xml = dom.toprettyxml(indent='  ')

    with open(model_path, 'w') as f:
        f.write(pretty_xml)


def _write_source_node(node_elem: ET.Element, node_idx: int, sn: NetworkStruct, classnames: List[str],
                       cs_node_names: Optional[Dict[Tuple[int, int], str]] = None):
    """Write source node section."""
    doc = node_elem

    section = ET.SubElement(node_elem, 'section')
    section.set('className', 'RandomSource')

    param = ET.SubElement(section, 'parameter')
    param.set('array', 'true')
    param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategy')
    param.set('name', 'ServiceStrategy')

    ist = int(sn.nodeToStation[node_idx]) if hasattr(sn, 'nodeToStation') and sn.nodeToStation is not None and node_idx < len(sn.nodeToStation) else 0
    K = sn.nclasses

    for r in range(K):
        ref_class = ET.SubElement(param, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(param, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy')
        sub_param.set('name', 'ServiceTimeStrategy')

        njobs = sn.njobs.flatten() if sn.njobs is not None else np.zeros(K)
        if not np.isinf(njobs[r]):
            # Closed class - no arrivals at source
            value = ET.SubElement(sub_param, 'value')
            value.text = 'null'
        else:
            # Open class - check arrival distribution type
            rate = sn.rates[ist, r] if sn.rates is not None and ist < sn.rates.shape[0] and r < sn.rates.shape[1] else 1.0
            if np.isnan(rate) or rate <= 0:
                value = ET.SubElement(sub_param, 'value')
                value.text = 'null'
            else:
                # Get the process type for arrivals
                procid = None
                if hasattr(sn, 'procid') and sn.procid is not None:
                    try:
                        procid = sn.procid[ist][r]
                    except (IndexError, TypeError, KeyError):
                        pass

                # Handle Phase-Type distributions (PH, APH, Coxian, Cox2)
                if procid in (ProcessType.PH, ProcessType.APH, ProcessType.COXIAN, ProcessType.COX2):
                    _write_phase_type_service_distribution(sub_param, sn, ist, r)
                elif procid == ProcessType.ERLANG:
                    # Erlang distribution
                    proc = None
                    if hasattr(sn, 'proc') and sn.proc is not None:
                        try:
                            proc = sn.proc[ist][r]
                        except (IndexError, TypeError, KeyError):
                            pass
                    phases = 2
                    if proc is not None and isinstance(proc, (list, tuple)) and len(proc) > 0:
                        T = np.asarray(proc[0], dtype=np.float64)
                        phases = T.shape[0] if T.ndim >= 1 else 2

                    distr = ET.SubElement(sub_param, 'subParameter')
                    distr.set('classPath', 'jmt.engine.random.Erlang')
                    distr.set('name', 'Erlang')

                    distr_par = ET.SubElement(sub_param, 'subParameter')
                    distr_par.set('classPath', 'jmt.engine.random.ErlangPar')
                    distr_par.set('name', 'distrPar')

                    alpha_param = ET.SubElement(distr_par, 'subParameter')
                    alpha_param.set('classPath', 'java.lang.Double')
                    alpha_param.set('name', 'alpha')
                    value = ET.SubElement(alpha_param, 'value')
                    value.text = f'{rate * phases:.12f}'

                    r_param = ET.SubElement(distr_par, 'subParameter')
                    r_param.set('classPath', 'java.lang.Long')
                    r_param.set('name', 'r')
                    value = ET.SubElement(r_param, 'value')
                    value.text = str(phases)
                elif procid == ProcessType.HYPEREXP:
                    # HyperExponential distribution. An n>2 HyperExp is emitted as
                    # the equivalent APH rather than being truncated to 2 phases.
                    proc = None
                    if hasattr(sn, 'proc') and sn.proc is not None:
                        try:
                            proc = sn.proc[ist][r]
                        except (IndexError, TypeError, KeyError):
                            pass
                    pie = None
                    if hasattr(sn, 'pie') and sn.pie is not None:
                        try:
                            pie = sn.pie[ist][r]
                        except (IndexError, TypeError, KeyError):
                            pass

                    params = _hyperexp_probs_rates(proc, pie)
                    if params is None:
                        params = (np.array([0.5, 0.5]), np.array([rate, rate]))
                    _write_hyperexp_distribution(sub_param, params[0], params[1])
                elif procid == ProcessType.DET:
                    # Deterministic distribution
                    distr = ET.SubElement(sub_param, 'subParameter')
                    distr.set('classPath', 'jmt.engine.random.DeterministicDistr')
                    distr.set('name', 'Deterministic')

                    distr_par = ET.SubElement(sub_param, 'subParameter')
                    distr_par.set('classPath', 'jmt.engine.random.DeterministicDistrPar')
                    distr_par.set('name', 'distrPar')

                    t_param = ET.SubElement(distr_par, 'subParameter')
                    t_param.set('classPath', 'java.lang.Double')
                    t_param.set('name', 't')
                    value = ET.SubElement(t_param, 'value')
                    value.text = f'{1.0 / rate:.12f}' if rate > 0 else '0.000000000000'
                elif procid == ProcessType.PARETO:
                    # Pareto distribution - reconstruct shape/scale from sn.scv and sn.rates
                    scv_val = sn.scv[ist, r] if hasattr(sn, 'scv') and sn.scv is not None else 1.0
                    shape = np.sqrt(1.0 + 1.0 / scv_val) + 1.0
                    scale = (1.0 / rate) * (shape - 1.0) / shape

                    distr = ET.SubElement(sub_param, 'subParameter')
                    distr.set('classPath', 'jmt.engine.random.Pareto')
                    distr.set('name', 'Pareto')

                    distr_par = ET.SubElement(sub_param, 'subParameter')
                    distr_par.set('classPath', 'jmt.engine.random.ParetoPar')
                    distr_par.set('name', 'distrPar')

                    alpha_param = ET.SubElement(distr_par, 'subParameter')
                    alpha_param.set('classPath', 'java.lang.Double')
                    alpha_param.set('name', 'alpha')
                    value = ET.SubElement(alpha_param, 'value')
                    value.text = f'{shape:.12f}'

                    k_param = ET.SubElement(distr_par, 'subParameter')
                    k_param.set('classPath', 'java.lang.Double')
                    k_param.set('name', 'k')
                    value = ET.SubElement(k_param, 'value')
                    value.text = f'{scale:.12f}'
                elif procid == ProcessType.GAMMA:
                    # Gamma distribution
                    scv_val = sn.scv[ist, r] if hasattr(sn, 'scv') and sn.scv is not None else 1.0
                    gamma_alpha = 1.0 / scv_val
                    gamma_beta = scv_val / rate

                    distr = ET.SubElement(sub_param, 'subParameter')
                    distr.set('classPath', 'jmt.engine.random.GammaDistr')
                    distr.set('name', 'Gamma')

                    distr_par = ET.SubElement(sub_param, 'subParameter')
                    distr_par.set('classPath', 'jmt.engine.random.GammaDistrPar')
                    distr_par.set('name', 'distrPar')

                    alpha_param = ET.SubElement(distr_par, 'subParameter')
                    alpha_param.set('classPath', 'java.lang.Double')
                    alpha_param.set('name', 'alpha')
                    value = ET.SubElement(alpha_param, 'value')
                    value.text = f'{gamma_alpha:.12f}'

                    beta_param = ET.SubElement(distr_par, 'subParameter')
                    beta_param.set('classPath', 'java.lang.Double')
                    beta_param.set('name', 'beta')
                    value = ET.SubElement(beta_param, 'value')
                    value.text = f'{gamma_beta:.12f}'
                elif procid == ProcessType.WEIBULL:
                    # Weibull distribution
                    from scipy.special import gamma as gamma_func
                    scv_val = sn.scv[ist, r] if hasattr(sn, 'scv') and sn.scv is not None else 1.0
                    c = np.sqrt(scv_val)
                    rval = c ** (-1.086)
                    weibull_alpha = (1.0 / rate) / gamma_func(1.0 + 1.0 / rval)

                    distr = ET.SubElement(sub_param, 'subParameter')
                    distr.set('classPath', 'jmt.engine.random.Weibull')
                    distr.set('name', 'Weibull')

                    distr_par = ET.SubElement(sub_param, 'subParameter')
                    distr_par.set('classPath', 'jmt.engine.random.WeibullPar')
                    distr_par.set('name', 'distrPar')

                    alpha_param = ET.SubElement(distr_par, 'subParameter')
                    alpha_param.set('classPath', 'java.lang.Double')
                    alpha_param.set('name', 'alpha')
                    value = ET.SubElement(alpha_param, 'value')
                    value.text = f'{weibull_alpha:.12f}'

                    r_param = ET.SubElement(distr_par, 'subParameter')
                    r_param.set('classPath', 'java.lang.Double')
                    r_param.set('name', 'r')
                    value = ET.SubElement(r_param, 'value')
                    value.text = f'{rval:.12f}'
                elif procid == ProcessType.LOGNORMAL:
                    # Lognormal distribution
                    scv_val = sn.scv[ist, r] if hasattr(sn, 'scv') and sn.scv is not None else 1.0
                    c = np.sqrt(scv_val)
                    mu = np.log((1.0 / rate) / np.sqrt(c * c + 1.0))
                    sigma = np.sqrt(np.log(c * c + 1.0))

                    distr = ET.SubElement(sub_param, 'subParameter')
                    distr.set('classPath', 'jmt.engine.random.Lognormal')
                    distr.set('name', 'Lognormal')

                    distr_par = ET.SubElement(sub_param, 'subParameter')
                    distr_par.set('classPath', 'jmt.engine.random.LognormalPar')
                    distr_par.set('name', 'distrPar')

                    mu_param = ET.SubElement(distr_par, 'subParameter')
                    mu_param.set('classPath', 'java.lang.Double')
                    mu_param.set('name', 'mu')
                    value = ET.SubElement(mu_param, 'value')
                    value.text = f'{mu:.12f}'

                    sigma_param = ET.SubElement(distr_par, 'subParameter')
                    sigma_param.set('classPath', 'java.lang.Double')
                    sigma_param.set('name', 'sigma')
                    value = ET.SubElement(sigma_param, 'value')
                    value.text = f'{sigma:.12f}'
                elif procid == ProcessType.UNIFORM:
                    # Uniform distribution
                    scv_val = sn.scv[ist, r] if hasattr(sn, 'scv') and sn.scv is not None else 1.0
                    maxVal = (np.sqrt(12.0 * scv_val / (rate * rate)) + 2.0 / rate) / 2.0
                    minVal = 2.0 / rate - maxVal

                    distr = ET.SubElement(sub_param, 'subParameter')
                    distr.set('classPath', 'jmt.engine.random.Uniform')
                    distr.set('name', 'Uniform')

                    distr_par = ET.SubElement(sub_param, 'subParameter')
                    distr_par.set('classPath', 'jmt.engine.random.UniformPar')
                    distr_par.set('name', 'distrPar')

                    min_param = ET.SubElement(distr_par, 'subParameter')
                    min_param.set('classPath', 'java.lang.Double')
                    min_param.set('name', 'min')
                    value = ET.SubElement(min_param, 'value')
                    value.text = f'{minVal:.12f}'

                    max_param = ET.SubElement(distr_par, 'subParameter')
                    max_param.set('classPath', 'java.lang.Double')
                    max_param.set('name', 'max')
                    value = ET.SubElement(max_param, 'value')
                    value.text = f'{maxVal:.12f}'
                else:
                    # Default: Exponential distribution
                    distr = ET.SubElement(sub_param, 'subParameter')
                    distr.set('classPath', 'jmt.engine.random.Exponential')
                    distr.set('name', 'Exponential')

                    distr_par = ET.SubElement(sub_param, 'subParameter')
                    distr_par.set('classPath', 'jmt.engine.random.ExponentialPar')
                    distr_par.set('name', 'distrPar')

                    lambda_param = ET.SubElement(distr_par, 'subParameter')
                    lambda_param.set('classPath', 'java.lang.Double')
                    lambda_param.set('name', 'lambda')
                    value = ET.SubElement(lambda_param, 'value')
                    value.text = f'{rate:.12f}'

    # ServiceTunnel and Router sections
    tunnel = ET.SubElement(node_elem, 'section')
    tunnel.set('className', 'ServiceTunnel')

    router = ET.SubElement(node_elem, 'section')
    router.set('className', 'Router')
    _write_routing_strategy(router, node_idx, sn, classnames, cs_node_names)


def _write_sink_node(node_elem: ET.Element, sn: NetworkStruct, classnames: List[str]):
    """Write sink node section."""
    section = ET.SubElement(node_elem, 'section')
    section.set('className', 'JobSink')


def _write_delay_node(node_elem: ET.Element, node_idx: int, sn: NetworkStruct, classnames: List[str],
                      cs_node_names: Optional[Dict[Tuple[int, int], str]] = None):
    """Write delay station node.

    Delay nodes (infinite servers) in JMT need:
    1. Queue section (input buffer)
    2. Delay section (service strategy)
    3. Router section (routing)
    """
    K = sn.nclasses
    ist = int(sn.nodeToStation[node_idx]) if hasattr(sn, 'nodeToStation') and sn.nodeToStation is not None and node_idx < len(sn.nodeToStation) else 0

    # Queue section (input buffer)
    queue = ET.SubElement(node_elem, 'section')
    queue.set('className', 'Queue')

    size_param = ET.SubElement(queue, 'parameter')
    size_param.set('classPath', 'java.lang.Integer')
    size_param.set('name', 'size')
    value = ET.SubElement(size_param, 'value')
    value.text = '-1'

    # Drop strategies
    drop_strategy = ET.SubElement(queue, 'parameter')
    drop_strategy.set('array', 'true')
    drop_strategy.set('classPath', 'java.lang.String')
    drop_strategy.set('name', 'dropStrategies')

    for r in range(K):
        ref_class = ET.SubElement(drop_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(drop_strategy, 'subParameter')
        sub_param.set('classPath', 'java.lang.String')
        sub_param.set('name', 'dropStrategy')
        value = ET.SubElement(sub_param, 'value')
        value.text = _drop_strategy_text(sn, ist, r)

    # Queue get strategy (FCFS)
    strategy_param = ET.SubElement(queue, 'parameter')
    strategy_param.set('classPath', 'jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy')
    strategy_param.set('name', 'FCFSstrategy')

    # Queue put strategy
    put_strategy = ET.SubElement(queue, 'parameter')
    put_strategy.set('array', 'true')
    put_strategy.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategy')
    put_strategy.set('name', 'QueuePutStrategy')

    for r in range(K):
        ref_class = ET.SubElement(put_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(put_strategy, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy')
        sub_param.set('name', 'TailStrategy')

    # Impatience section (null for no impatience, required for JMT compatibility)
    impatience_param = ET.SubElement(queue, 'parameter')
    impatience_param.set('array', 'true')
    impatience_param.set('classPath', 'jmt.engine.NetStrategies.ImpatienceStrategies.Impatience')
    impatience_param.set('name', 'Impatience')

    for r in range(K):
        ref_class = ET.SubElement(impatience_param, 'refClass')
        ref_class.text = classnames[r]

        # The ARRAY is named for the abstract Impatience type; each ENTRY names
        # the concrete strategy, which is Reneging even when it carries no
        # distribution (MATLAB saveImpatience.m). A Delay never renege-s, so the
        # value stays null.
        sub_param = ET.SubElement(impatience_param, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.ImpatienceStrategies.Reneging')
        sub_param.set('name', 'Reneging')
        value = ET.SubElement(sub_param, 'value')
        value.text = 'null'

    # Delay section (service)
    server = ET.SubElement(node_elem, 'section')
    server.set('className', 'Delay')

    service_param = ET.SubElement(server, 'parameter')
    service_param.set('array', 'true')
    service_param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategy')
    service_param.set('name', 'ServiceStrategy')

    for r in range(K):
        ref_class = ET.SubElement(service_param, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(service_param, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy')
        sub_param.set('name', 'ServiceTimeStrategy')

        rate = sn.rates[ist, r] if sn.rates is not None and ist < sn.rates.shape[0] and r < sn.rates.shape[1] else 1.0

        # Get process type
        procid = None
        if hasattr(sn, 'procid') and sn.procid is not None:
            if ist < sn.procid.shape[0] and r < sn.procid.shape[1]:
                procid = sn.procid[ist, r]

        if np.isnan(rate) or rate < 0 or procid == ProcessType.DISABLED:
            # Disabled service - use DisabledServiceTimeStrategy like MATLAB
            sub_param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.DisabledServiceTimeStrategy')
            sub_param.set('name', 'DisabledServiceTimeStrategy')
        elif rate == 0 or procid == ProcessType.IMMEDIATE:
            # Immediate service (zero service time)
            sub_param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ZeroServiceTimeStrategy')
            sub_param.set('name', 'ZeroServiceTimeStrategy')
        elif procid in (ProcessType.MAP, ProcessType.MMPP2):
            # MAP/MMPP2 distribution
            _write_map_service_distribution(sub_param, sn, ist, r)
        elif procid == ProcessType.ERLANG:
            # Erlang distribution
            distr = ET.SubElement(sub_param, 'subParameter')
            distr.set('classPath', 'jmt.engine.random.Erlang')
            distr.set('name', 'Erlang')

            distr_par = ET.SubElement(sub_param, 'subParameter')
            distr_par.set('classPath', 'jmt.engine.random.ErlangPar')
            distr_par.set('name', 'distrPar')

            # Get number of phases from sn.proc or sn.phases
            phases = 1
            if hasattr(sn, 'proc') and sn.proc is not None:
                try:
                    proc_entry = sn.proc[ist][r]
                    if isinstance(proc_entry, dict) and 'k' in proc_entry:
                        phases = int(proc_entry['k'])
                except (IndexError, TypeError, KeyError):
                    pass
            if phases == 1 and hasattr(sn, 'phases') and sn.phases is not None:
                if ist < sn.phases.shape[0] and r < sn.phases.shape[1]:
                    phases = int(sn.phases[ist, r])

            # Alpha = rate * phases
            alpha_param = ET.SubElement(distr_par, 'subParameter')
            alpha_param.set('classPath', 'java.lang.Double')
            alpha_param.set('name', 'alpha')
            value = ET.SubElement(alpha_param, 'value')
            value.text = f'{rate * phases:.12f}'

            # r = number of phases
            r_param = ET.SubElement(distr_par, 'subParameter')
            r_param.set('classPath', 'java.lang.Long')
            r_param.set('name', 'r')
            value = ET.SubElement(r_param, 'value')
            value.text = str(phases)
        elif procid == ProcessType.HYPEREXP:
            # see _kb/06-solver-catalog.md (JMT: "phase-type distribution
            # export is shared across service, impatience, retrial")
            params = None
            if hasattr(sn, 'proc') and sn.proc is not None:
                try:
                    pie_entry = None
                    if getattr(sn, 'pie', None) is not None:
                        pie_entry = sn.pie[ist][r]
                    params = _hyperexp_probs_rates(sn.proc[ist][r], pie_entry)
                except (IndexError, TypeError, KeyError):
                    params = None
            if params is None:
                params = (np.array([0.5, 0.5]), np.array([rate, rate]))
            _write_hyperexp_distribution(sub_param, params[0], params[1])
        elif procid in (ProcessType.PH, ProcessType.APH, ProcessType.COXIAN, ProcessType.COX2):
            # Phase-type distributions
            _write_phase_type_service_distribution(sub_param, sn, ist, r)
        elif procid == ProcessType.PARETO:
            # Pareto distribution - reconstruct shape/scale from sn.scv and sn.rates
            # MATLAB formula: shape = sqrt(1+1/sn.scv(i,r))+1; scale = 1/sn.rates(i,r) * (shape-1)/shape
            scv_val = sn.scv[ist, r] if hasattr(sn, 'scv') and sn.scv is not None else 1.0
            shape = np.sqrt(1.0 + 1.0 / scv_val) + 1.0
            scale = (1.0 / rate) * (shape - 1.0) / shape

            distr = ET.SubElement(sub_param, 'subParameter')
            distr.set('classPath', 'jmt.engine.random.Pareto')
            distr.set('name', 'Pareto')

            distr_par = ET.SubElement(sub_param, 'subParameter')
            distr_par.set('classPath', 'jmt.engine.random.ParetoPar')
            distr_par.set('name', 'distrPar')

            alpha_param = ET.SubElement(distr_par, 'subParameter')
            alpha_param.set('classPath', 'java.lang.Double')
            alpha_param.set('name', 'alpha')
            value = ET.SubElement(alpha_param, 'value')
            value.text = f'{shape:.12f}'

            k_param = ET.SubElement(distr_par, 'subParameter')
            k_param.set('classPath', 'java.lang.Double')
            k_param.set('name', 'k')
            value = ET.SubElement(k_param, 'value')
            value.text = f'{scale:.12f}'
        elif procid == ProcessType.GAMMA:
            # Gamma distribution - reconstruct from sn.scv and sn.rates
            scv_val = sn.scv[ist, r] if hasattr(sn, 'scv') and sn.scv is not None else 1.0
            gamma_alpha = 1.0 / scv_val
            gamma_beta = scv_val / rate

            distr = ET.SubElement(sub_param, 'subParameter')
            distr.set('classPath', 'jmt.engine.random.GammaDistr')
            distr.set('name', 'Gamma')

            distr_par = ET.SubElement(sub_param, 'subParameter')
            distr_par.set('classPath', 'jmt.engine.random.GammaDistrPar')
            distr_par.set('name', 'distrPar')

            alpha_param = ET.SubElement(distr_par, 'subParameter')
            alpha_param.set('classPath', 'java.lang.Double')
            alpha_param.set('name', 'alpha')
            value = ET.SubElement(alpha_param, 'value')
            value.text = f'{gamma_alpha:.12f}'

            beta_param = ET.SubElement(distr_par, 'subParameter')
            beta_param.set('classPath', 'java.lang.Double')
            beta_param.set('name', 'beta')
            value = ET.SubElement(beta_param, 'value')
            value.text = f'{gamma_beta:.12f}'
        elif procid == ProcessType.WEIBULL:
            # Weibull distribution - reconstruct from sn.scv and sn.rates
            from scipy.special import gamma as gamma_func
            scv_val = sn.scv[ist, r] if hasattr(sn, 'scv') and sn.scv is not None else 1.0
            c = np.sqrt(scv_val)
            rval = c ** (-1.086)  # Justus approximation (1976)
            weibull_alpha = (1.0 / rate) / gamma_func(1.0 + 1.0 / rval)

            distr = ET.SubElement(sub_param, 'subParameter')
            distr.set('classPath', 'jmt.engine.random.Weibull')
            distr.set('name', 'Weibull')

            distr_par = ET.SubElement(sub_param, 'subParameter')
            distr_par.set('classPath', 'jmt.engine.random.WeibullPar')
            distr_par.set('name', 'distrPar')

            alpha_param = ET.SubElement(distr_par, 'subParameter')
            alpha_param.set('classPath', 'java.lang.Double')
            alpha_param.set('name', 'alpha')
            value = ET.SubElement(alpha_param, 'value')
            value.text = f'{weibull_alpha:.12f}'

            r_param = ET.SubElement(distr_par, 'subParameter')
            r_param.set('classPath', 'java.lang.Double')
            r_param.set('name', 'r')
            value = ET.SubElement(r_param, 'value')
            value.text = f'{rval:.12f}'
        elif procid == ProcessType.LOGNORMAL:
            # Lognormal distribution - reconstruct from sn.scv and sn.rates
            scv_val = sn.scv[ist, r] if hasattr(sn, 'scv') and sn.scv is not None else 1.0
            c = np.sqrt(scv_val)
            mu = np.log((1.0 / rate) / np.sqrt(c * c + 1.0))
            sigma = np.sqrt(np.log(c * c + 1.0))

            distr = ET.SubElement(sub_param, 'subParameter')
            distr.set('classPath', 'jmt.engine.random.Lognormal')
            distr.set('name', 'Lognormal')

            distr_par = ET.SubElement(sub_param, 'subParameter')
            distr_par.set('classPath', 'jmt.engine.random.LognormalPar')
            distr_par.set('name', 'distrPar')

            mu_param = ET.SubElement(distr_par, 'subParameter')
            mu_param.set('classPath', 'java.lang.Double')
            mu_param.set('name', 'mu')
            value = ET.SubElement(mu_param, 'value')
            value.text = f'{mu:.12f}'

            sigma_param = ET.SubElement(distr_par, 'subParameter')
            sigma_param.set('classPath', 'java.lang.Double')
            sigma_param.set('name', 'sigma')
            value = ET.SubElement(sigma_param, 'value')
            value.text = f'{sigma:.12f}'
        elif procid == ProcessType.UNIFORM:
            # Uniform distribution - reconstruct from sn.scv and sn.rates
            scv_val = sn.scv[ist, r] if hasattr(sn, 'scv') and sn.scv is not None else 1.0
            maxVal = (np.sqrt(12.0 * scv_val / (rate * rate)) + 2.0 / rate) / 2.0
            minVal = 2.0 / rate - maxVal

            distr = ET.SubElement(sub_param, 'subParameter')
            distr.set('classPath', 'jmt.engine.random.Uniform')
            distr.set('name', 'Uniform')

            distr_par = ET.SubElement(sub_param, 'subParameter')
            distr_par.set('classPath', 'jmt.engine.random.UniformPar')
            distr_par.set('name', 'distrPar')

            min_param = ET.SubElement(distr_par, 'subParameter')
            min_param.set('classPath', 'java.lang.Double')
            min_param.set('name', 'min')
            value = ET.SubElement(min_param, 'value')
            value.text = f'{minVal:.12f}'

            max_param = ET.SubElement(distr_par, 'subParameter')
            max_param.set('classPath', 'java.lang.Double')
            max_param.set('name', 'max')
            value = ET.SubElement(max_param, 'value')
            value.text = f'{maxVal:.12f}'
        elif procid == ProcessType.DET:
            # Deterministic distribution
            distr = ET.SubElement(sub_param, 'subParameter')
            distr.set('classPath', 'jmt.engine.random.DeterministicDistr')
            distr.set('name', 'Deterministic')

            distr_par = ET.SubElement(sub_param, 'subParameter')
            distr_par.set('classPath', 'jmt.engine.random.DeterministicDistrPar')
            distr_par.set('name', 'distrPar')

            t_param = ET.SubElement(distr_par, 'subParameter')
            t_param.set('classPath', 'java.lang.Double')
            t_param.set('name', 't')
            value = ET.SubElement(t_param, 'value')
            value.text = f'{1.0 / rate:.12f}'
        else:
            # Default: Exponential distribution
            distr = ET.SubElement(sub_param, 'subParameter')
            distr.set('classPath', 'jmt.engine.random.Exponential')
            distr.set('name', 'Exponential')

            distr_par = ET.SubElement(sub_param, 'subParameter')
            distr_par.set('classPath', 'jmt.engine.random.ExponentialPar')
            distr_par.set('name', 'distrPar')

            lambda_param = ET.SubElement(distr_par, 'subParameter')
            lambda_param.set('classPath', 'java.lang.Double')
            lambda_param.set('name', 'lambda')
            value = ET.SubElement(lambda_param, 'value')
            value.text = f'{rate:.12f}'

    # Router section
    router = ET.SubElement(node_elem, 'section')
    router.set('className', 'Router')
    _write_routing_strategy(router, node_idx, sn, classnames, cs_node_names)


def _write_queue_node(node_elem: ET.Element, node_idx: int, sn: NetworkStruct, classnames: List[str],
                      options: SolverJMTOptions, cs_node_names: Optional[Dict[Tuple[int, int], str]] = None,
                      model: Any = None):
    """Write queue station node."""
    K = sn.nclasses
    ist = int(sn.nodeToStation[node_idx]) if hasattr(sn, 'nodeToStation') and sn.nodeToStation is not None and node_idx < len(sn.nodeToStation) else 0

    # Get scheduling strategy
    sched = SchedStrategy.FCFS
    if sn.sched and ist in sn.sched:
        sched = sn.sched[ist]

    # Queue section
    queue = ET.SubElement(node_elem, 'section')
    queue.set('className', 'Queue')

    # Size parameter: JMT "size" = K (total capacity), -1 = infinite.
    # Mirrors MATLAB saveBufferCapacity.m: cap>=sum(njobs) also treated as infinite.
    size_param = ET.SubElement(queue, 'parameter')
    size_param.set('classPath', 'java.lang.Integer')
    size_param.set('name', 'size')
    value = ET.SubElement(size_param, 'value')
    capacity = -1  # Default: infinite capacity
    if hasattr(sn, 'cap') and sn.cap is not None and ist < len(sn.cap):
        cap_val = sn.cap[ist]
        if not np.isinf(cap_val):
            # A capacity the population cannot reach is unbounded, and the test
            # is >= and not ==: refresh_capacity DERIVES sn.cap for a station the
            # user never capped, as sum over the classes served there of the
            # chain population, so a multi-class station gets (#classes) x N --
            # 8 on a two-class model of 4 jobs. Under == only the single-class
            # case matched, and every multi-class one fell through to
            # _jmt_station_cap_assert and was refused as a "binding" buffer
            # nobody declared.
            # An absent njobs establishes no bound, so it must not read as one:
            # inf keeps the capacity on the exported-and-asserted path, where
            # the old `== 0` left it, rather than silently dropping it.
            #
            # THE POPULATION COMPARED AGAINST IS THE ONE THAT CAN REACH ist, not
            # the model total. A class that never visits this station cannot
            # fill it, so counting its jobs makes a capacity that is exactly the
            # reachable population look like a buffer -- which is what a
            # SELF-LOOPING CLASS does. See the MATLAB twin in
            # JMTIO/saveBufferCapacity.m.
            total_jobs = _jmt_reachable_population(sn, ist)
            if cap_val >= total_jobs:
                capacity = -1
            else:
                # Check for infinite servers (delay node) - no buffer needed
                nservers = sn.nservers[ist] if hasattr(sn, 'nservers') and sn.nservers is not None and ist < len(sn.nservers) else 1
                if np.isinf(nservers):
                    capacity = -1
                else:
                    _jmt_station_cap_assert(sn, ist)
                    capacity = int(cap_val)
    value.text = str(capacity)

    # 2. Drop strategies (required)
    drop_strategy = ET.SubElement(queue, 'parameter')
    drop_strategy.set('array', 'true')
    drop_strategy.set('classPath', 'java.lang.String')
    drop_strategy.set('name', 'dropStrategies')

    for r in range(K):
        ref_class = ET.SubElement(drop_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(drop_strategy, 'subParameter')
        sub_param.set('classPath', 'java.lang.String')
        sub_param.set('name', 'dropStrategy')
        value = ET.SubElement(sub_param, 'value')
        # Port of MATLAB saveDropStrategy.m: 0 is not a DropStrategy value, it
        # is the "no rule recorded" slot, which JMT spells 'drop' (the field a
        # bufferless node carries); every other value goes through
        # DropStrategy.toText. WAITQ is -1 in all three codebases, so the
        # earlier `drop_val == 0 -> waiting queue` test never fired and every
        # blocking queue was exported as a DROPPING one instead.
        value.text = _drop_strategy_text(sn, ist, r)

    # 3. Queue get strategy
    # Check if this is a polling queue
    is_polling_queue = False
    polling_type = None
    polling_k = 1

    if sched == SchedStrategy.POLLING and sn.nodeparam and node_idx in sn.nodeparam:
        nodeparam = sn.nodeparam[node_idx]
        if isinstance(nodeparam, dict):
            # Get polling type from first class that has it
            for r in range(K):
                if r in nodeparam and isinstance(nodeparam[r], dict):
                    if 'pollingType' in nodeparam[r]:
                        is_polling_queue = True
                        polling_type = nodeparam[r]['pollingType']
                        polling_par = nodeparam[r].get('pollingPar', [1])
                        polling_k = polling_par[0] if polling_par else 1
                        break

    if is_polling_queue and polling_type is not None:
        _write_polling_get_strategy(queue, polling_type, polling_k)
    else:
        strategy_param = ET.SubElement(queue, 'parameter')
        strategy_param.set('classPath', _get_sched_strategy_class(sched))
        strategy_param.set('name', 'FCFSstrategy')

    # 4. Queue put strategy
    put_strategy = ET.SubElement(queue, 'parameter')
    put_strategy.set('array', 'true')
    put_strategy.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategy')
    put_strategy.set('name', 'QueuePutStrategy')

    for r in range(K):
        ref_class = ET.SubElement(put_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(put_strategy, 'subParameter')
        # Select put strategy based on scheduling discipline
        if sched == SchedStrategy.SIRO:
            sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.RandStrategy')
            sub_param.set('name', 'RandStrategy')
        elif sched == SchedStrategy.LCFS:
            sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.HeadStrategy')
            sub_param.set('name', 'HeadStrategy')
        elif sched == SchedStrategy.LCFSPRIO:
            sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.HeadStrategyPriority')
            sub_param.set('name', 'HeadStrategyPriority')
        elif sched == SchedStrategy.LCFSPR:
            sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.LCFSPRStrategy')
            sub_param.set('name', 'LCFSPRStrategy')
        elif sched == SchedStrategy.LCFSPRPRIO:
            sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.LCFSPRStrategyPriority')
            sub_param.set('name', 'LCFSPRStrategyPriority')
        elif sched == SchedStrategy.LCFSPI:
            sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.LCFSPIStrategy')
            sub_param.set('name', 'LCFSPIStrategy')
        elif sched == SchedStrategy.HOL:
            sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.TailStrategyPriority')
            sub_param.set('name', 'TailStrategyPriority')
        elif sched == SchedStrategy.SJF:
            sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.SJFStrategy')
            sub_param.set('name', 'SJFStrategy')
        elif sched == SchedStrategy.LJF:
            sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.LJFStrategy')
            sub_param.set('name', 'LJFStrategy')
        elif sched == SchedStrategy.SEPT:
            sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.SEPTStrategy')
            sub_param.set('name', 'SEPTStrategy')
        elif sched == SchedStrategy.LEPT:
            sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.LEPTStrategy')
            sub_param.set('name', 'LEPTStrategy')
        else:
            # Default: TailStrategy for FCFS, PS, and other strategies
            sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy')
            sub_param.set('name', 'TailStrategy')

    # see _kb/06-solver-catalog.md (JMT: "impatience (reneging/balking) export")
    impatience_param = ET.SubElement(queue, 'parameter')
    impatience_param.set('array', 'true')
    impatience_param.set('classPath', 'jmt.engine.NetStrategies.ImpatienceStrategies.Impatience')
    impatience_param.set('name', 'Impatience')

    for r in range(K):
        ref_class = ET.SubElement(impatience_param, 'refClass')
        ref_class.text = classnames[r]

        pdist = None
        if getattr(sn, 'impatienceDist', None) is not None                 and 0 <= ist < len(sn.impatienceDist) and r < len(sn.impatienceDist[ist]):
            pdist = sn.impatienceDist[ist][r]
        imp_sub = ET.SubElement(impatience_param, 'subParameter')
        # see _kb/06-solver-catalog.md (JMT: "impatience (reneging/balking) export")
        if _has_queue_length_balking(sn, ist, r):
            imp_sub.set('classPath', 'jmt.engine.NetStrategies.ImpatienceStrategies.Balking')
            imp_sub.set('name', 'Balking')
            nservers_ist = None
            if sn.nservers is not None and ist < np.asarray(sn.nservers).ravel().size:
                nservers_ist = np.asarray(sn.nservers).ravel()[ist]
            _write_balking_strategy(imp_sub, _balking_thresholds(sn, ist, r), nservers_ist)
        elif pdist is not None:
            imp_sub.set('classPath', 'jmt.engine.NetStrategies.ImpatienceStrategies.Reneging')
            imp_sub.set('name', 'Reneging')
            _write_distribution_param(imp_sub, pdist)
        else:
            # No impatience configured: still a Reneging entry, carrying null.
            # Naming the entry after the abstract Impatience type left JMT's
            # SimLoader with a class it cannot instantiate (MATLAB
            # saveImpatience.m emits Reneging in both cases).
            imp_sub.set('classPath', 'jmt.engine.NetStrategies.ImpatienceStrategies.Reneging')
            imp_sub.set('name', 'Reneging')
            value = ET.SubElement(imp_sub, 'value')
            value.text = 'null'

    # Server section
    server = ET.SubElement(node_elem, 'section')
    # Use PSServer for processor sharing variants (PS, DPS, GPS, and priority variants)
    # Use PollingServer variants for POLLING scheduling
    if sched in (SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS,
                 SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO,
                 SchedStrategy.LPS):
        server.set('className', 'PSServer')
    elif sched == SchedStrategy.POLLING and is_polling_queue and polling_type is not None:
        # Use appropriate polling server class based on polling type
        if polling_type == PollingType.GATED:
            server.set('className', 'GatedPollingServer')
        elif polling_type == PollingType.EXHAUSTIVE:
            server.set('className', 'ExhaustivePollingServer')
        elif polling_type == PollingType.KLIMITED:
            server.set('className', 'LimitedPollingServer')
        else:
            server.set('className', 'ExhaustivePollingServer')  # Default
    else:
        server.set('className', 'Server')

    # see _kb/06-solver-catalog.md (JMT: "RNG consumption, timing authority,
    # SPN wiring, servers") -- number of servers via lldscaling
    nservers = 1
    if sn.nservers is not None:
        nservers_val = sn.nservers[ist] if len(sn.nservers.shape) == 1 else sn.nservers[ist, 0]
        if np.isinf(nservers_val):
            nservers = 1000000  # Very large number for infinite servers
        else:
            nservers = int(nservers_val)

    # Check load-dependent scaling for effective number of servers
    if hasattr(sn, 'lldscaling') and sn.lldscaling is not None:
        if ist < sn.lldscaling.shape[0]:
            effective_servers = int(np.max(sn.lldscaling[ist, :]))
            if nservers < 1000000:  # Don't override infinite servers
                nservers = max(nservers, effective_servers)

    servers_param = ET.SubElement(server, 'parameter')
    servers_param.set('classPath', 'java.lang.Integer')
    servers_param.set('name', 'maxJobs')
    value = ET.SubElement(servers_param, 'value')
    value.text = str(nservers)

    # Number of visits
    visits_param = ET.SubElement(server, 'parameter')
    visits_param.set('array', 'true')
    visits_param.set('classPath', 'java.lang.Integer')
    visits_param.set('name', 'numberOfVisits')

    for r in range(K):
        ref_class = ET.SubElement(visits_param, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(visits_param, 'subParameter')
        sub_param.set('classPath', 'java.lang.Integer')
        sub_param.set('name', 'numberOfVisits')
        value = ET.SubElement(sub_param, 'value')
        value.text = '1'

    # Service strategy
    service_param = ET.SubElement(server, 'parameter')
    service_param.set('array', 'true')
    service_param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategy')
    service_param.set('name', 'ServiceStrategy')

    for r in range(K):
        ref_class = ET.SubElement(service_param, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(service_param, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy')
        sub_param.set('name', 'ServiceTimeStrategy')

        rate = sn.rates[ist, r] if sn.rates is not None and ist < sn.rates.shape[0] and r < sn.rates.shape[1] else 1.0

        # Get process type
        procid = None
        if hasattr(sn, 'procid') and sn.procid is not None:
            if ist < sn.procid.shape[0] and r < sn.procid.shape[1]:
                procid = sn.procid[ist, r]

        if np.isnan(rate) or rate < 0 or procid == ProcessType.DISABLED:
            # Disabled service - use DisabledServiceTimeStrategy like MATLAB
            sub_param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.DisabledServiceTimeStrategy')
            sub_param.set('name', 'DisabledServiceTimeStrategy')
        elif rate == 0 or procid == ProcessType.IMMEDIATE:
            # Immediate service (zero service time)
            sub_param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ZeroServiceTimeStrategy')
            sub_param.set('name', 'ZeroServiceTimeStrategy')
        elif procid in (ProcessType.MAP, ProcessType.MMPP2):
            # MAP/MMPP2 distribution
            _write_map_service_distribution(sub_param, sn, ist, r)
        elif procid == ProcessType.ERLANG:
            # Erlang distribution
            distr = ET.SubElement(sub_param, 'subParameter')
            distr.set('classPath', 'jmt.engine.random.Erlang')
            distr.set('name', 'Erlang')

            distr_par = ET.SubElement(sub_param, 'subParameter')
            distr_par.set('classPath', 'jmt.engine.random.ErlangPar')
            distr_par.set('name', 'distrPar')

            # Get number of phases from sn.proc or sn.phases
            phases = 1
            # First try sn.proc[ist][r]['k'] (Python stores Erlang info there)
            if hasattr(sn, 'proc') and sn.proc is not None:
                try:
                    proc_entry = sn.proc[ist][r]
                    if isinstance(proc_entry, dict) and 'k' in proc_entry:
                        phases = int(proc_entry['k'])
                except (IndexError, TypeError, KeyError):
                    pass
            # Fallback to sn.phases if available
            if phases == 1 and hasattr(sn, 'phases') and sn.phases is not None:
                if ist < sn.phases.shape[0] and r < sn.phases.shape[1]:
                    phases = int(sn.phases[ist, r])

            # Alpha = rate * phases (MATLAB: sn.rates(i,r)*sn.phases(i,r))
            alpha_param = ET.SubElement(distr_par, 'subParameter')
            alpha_param.set('classPath', 'java.lang.Double')
            alpha_param.set('name', 'alpha')
            value = ET.SubElement(alpha_param, 'value')
            value.text = f'{rate * phases:.12f}'

            # r = number of phases
            r_param = ET.SubElement(distr_par, 'subParameter')
            r_param.set('classPath', 'java.lang.Long')
            r_param.set('name', 'r')
            value = ET.SubElement(r_param, 'value')
            value.text = str(phases)
        elif procid == ProcessType.HYPEREXP:
            # see _kb/06-solver-catalog.md (JMT: "phase-type distribution
            # export is shared across service, impatience, retrial")
            params = None
            if hasattr(sn, 'proc') and sn.proc is not None:
                try:
                    pie_entry = None
                    if getattr(sn, 'pie', None) is not None:
                        pie_entry = sn.pie[ist][r]
                    params = _hyperexp_probs_rates(sn.proc[ist][r], pie_entry)
                except (IndexError, TypeError, KeyError):
                    params = None
            if params is None:
                params = (np.array([0.5, 0.5]), np.array([rate, rate]))
            _write_hyperexp_distribution(sub_param, params[0], params[1])
        elif procid in (ProcessType.PH, ProcessType.APH, ProcessType.COXIAN, ProcessType.COX2):
            # Phase-type distributions
            _write_phase_type_service_distribution(sub_param, sn, ist, r)
        elif procid == ProcessType.REPLAYER:
            # Replayer distribution - reads service times from a trace file
            distr = ET.SubElement(sub_param, 'subParameter')
            distr.set('classPath', 'jmt.engine.random.Replayer')
            distr.set('name', 'Replayer')

            distr_par = ET.SubElement(sub_param, 'subParameter')
            distr_par.set('classPath', 'jmt.engine.random.ReplayerPar')
            distr_par.set('name', 'distrPar')

            # Get the file path from sn.proc[ist][r]
            file_path = ''
            if hasattr(sn, 'proc') and sn.proc is not None:
                try:
                    proc_entry = sn.proc[ist][r]
                    if isinstance(proc_entry, dict):
                        if 'file_path' in proc_entry:
                            file_path = proc_entry['file_path']
                        elif 'filePath' in proc_entry:
                            file_path = proc_entry['filePath']
                        elif 'fileName' in proc_entry:
                            file_path = proc_entry['fileName']
                except (IndexError, TypeError, KeyError):
                    pass

            file_param = ET.SubElement(distr_par, 'subParameter')
            file_param.set('classPath', 'java.lang.String')
            file_param.set('name', 'fileName')
            value = ET.SubElement(file_param, 'value')
            value.text = str(file_path)
        elif procid == ProcessType.PARETO:
            # Pareto distribution - reconstruct shape/scale from sn.scv and sn.rates
            # MATLAB formula: shape = sqrt(1+1/sn.scv(i,r))+1; scale = 1/sn.rates(i,r) * (shape-1)/shape
            scv_val = sn.scv[ist, r] if hasattr(sn, 'scv') and sn.scv is not None else 1.0
            shape = np.sqrt(1.0 + 1.0 / scv_val) + 1.0
            scale = (1.0 / rate) * (shape - 1.0) / shape

            distr = ET.SubElement(sub_param, 'subParameter')
            distr.set('classPath', 'jmt.engine.random.Pareto')
            distr.set('name', 'Pareto')

            distr_par = ET.SubElement(sub_param, 'subParameter')
            distr_par.set('classPath', 'jmt.engine.random.ParetoPar')
            distr_par.set('name', 'distrPar')

            alpha_param = ET.SubElement(distr_par, 'subParameter')
            alpha_param.set('classPath', 'java.lang.Double')
            alpha_param.set('name', 'alpha')
            value = ET.SubElement(alpha_param, 'value')
            value.text = f'{shape:.12f}'

            k_param = ET.SubElement(distr_par, 'subParameter')
            k_param.set('classPath', 'java.lang.Double')
            k_param.set('name', 'k')
            value = ET.SubElement(k_param, 'value')
            value.text = f'{scale:.12f}'
        elif procid == ProcessType.GAMMA:
            # Gamma distribution - reconstruct from sn.scv and sn.rates
            # MATLAB: alpha = 1/sn.scv(i,r); beta = sn.scv(i,r)/sn.rates(i,r)
            scv_val = sn.scv[ist, r] if hasattr(sn, 'scv') and sn.scv is not None else 1.0
            gamma_alpha = 1.0 / scv_val
            gamma_beta = scv_val / rate

            distr = ET.SubElement(sub_param, 'subParameter')
            distr.set('classPath', 'jmt.engine.random.GammaDistr')
            distr.set('name', 'Gamma')

            distr_par = ET.SubElement(sub_param, 'subParameter')
            distr_par.set('classPath', 'jmt.engine.random.GammaDistrPar')
            distr_par.set('name', 'distrPar')

            alpha_param = ET.SubElement(distr_par, 'subParameter')
            alpha_param.set('classPath', 'java.lang.Double')
            alpha_param.set('name', 'alpha')
            value = ET.SubElement(alpha_param, 'value')
            value.text = f'{gamma_alpha:.12f}'

            beta_param = ET.SubElement(distr_par, 'subParameter')
            beta_param.set('classPath', 'java.lang.Double')
            beta_param.set('name', 'beta')
            value = ET.SubElement(beta_param, 'value')
            value.text = f'{gamma_beta:.12f}'
        elif procid == ProcessType.WEIBULL:
            # Weibull distribution - reconstruct from sn.scv and sn.rates
            # MATLAB: c = sqrt(sn.scv(i,r)); rval = c^(-1.086); alpha = 1/sn.rates(i,r) / gamma(1+1/rval)
            from scipy.special import gamma as gamma_func
            scv_val = sn.scv[ist, r] if hasattr(sn, 'scv') and sn.scv is not None else 1.0
            c = np.sqrt(scv_val)
            rval = c ** (-1.086)  # Justus approximation (1976)
            weibull_alpha = (1.0 / rate) / gamma_func(1.0 + 1.0 / rval)

            distr = ET.SubElement(sub_param, 'subParameter')
            distr.set('classPath', 'jmt.engine.random.Weibull')
            distr.set('name', 'Weibull')

            distr_par = ET.SubElement(sub_param, 'subParameter')
            distr_par.set('classPath', 'jmt.engine.random.WeibullPar')
            distr_par.set('name', 'distrPar')

            alpha_param = ET.SubElement(distr_par, 'subParameter')
            alpha_param.set('classPath', 'java.lang.Double')
            alpha_param.set('name', 'alpha')
            value = ET.SubElement(alpha_param, 'value')
            value.text = f'{weibull_alpha:.12f}'

            r_param = ET.SubElement(distr_par, 'subParameter')
            r_param.set('classPath', 'java.lang.Double')
            r_param.set('name', 'r')
            value = ET.SubElement(r_param, 'value')
            value.text = f'{rval:.12f}'
        elif procid == ProcessType.LOGNORMAL:
            # Lognormal distribution - reconstruct from sn.scv and sn.rates
            # MATLAB: c = sqrt(sn.scv(i,r)); mu = log(1/sn.rates(i,r) / sqrt(c*c+1)); sigma = sqrt(log(c*c+1))
            scv_val = sn.scv[ist, r] if hasattr(sn, 'scv') and sn.scv is not None else 1.0
            c = np.sqrt(scv_val)
            mu = np.log((1.0 / rate) / np.sqrt(c * c + 1.0))
            sigma = np.sqrt(np.log(c * c + 1.0))

            distr = ET.SubElement(sub_param, 'subParameter')
            distr.set('classPath', 'jmt.engine.random.Lognormal')
            distr.set('name', 'Lognormal')

            distr_par = ET.SubElement(sub_param, 'subParameter')
            distr_par.set('classPath', 'jmt.engine.random.LognormalPar')
            distr_par.set('name', 'distrPar')

            mu_param = ET.SubElement(distr_par, 'subParameter')
            mu_param.set('classPath', 'java.lang.Double')
            mu_param.set('name', 'mu')
            value = ET.SubElement(mu_param, 'value')
            value.text = f'{mu:.12f}'

            sigma_param = ET.SubElement(distr_par, 'subParameter')
            sigma_param.set('classPath', 'java.lang.Double')
            sigma_param.set('name', 'sigma')
            value = ET.SubElement(sigma_param, 'value')
            value.text = f'{sigma:.12f}'
        elif procid == ProcessType.UNIFORM:
            # see _kb/06-solver-catalog.md (JMT: "RNG consumption, timing
            # authority, SPN wiring, servers") -- Uniform reconstruction
            scv_val = sn.scv[ist, r] if hasattr(sn, 'scv') and sn.scv is not None else 1.0
            maxVal = (np.sqrt(12.0 * scv_val / (rate * rate)) + 2.0 / rate) / 2.0
            minVal = 2.0 / rate - maxVal

            distr = ET.SubElement(sub_param, 'subParameter')
            distr.set('classPath', 'jmt.engine.random.Uniform')
            distr.set('name', 'Uniform')

            distr_par = ET.SubElement(sub_param, 'subParameter')
            distr_par.set('classPath', 'jmt.engine.random.UniformPar')
            distr_par.set('name', 'distrPar')

            min_param = ET.SubElement(distr_par, 'subParameter')
            min_param.set('classPath', 'java.lang.Double')
            min_param.set('name', 'min')
            value = ET.SubElement(min_param, 'value')
            value.text = f'{minVal:.12f}'

            max_param = ET.SubElement(distr_par, 'subParameter')
            max_param.set('classPath', 'java.lang.Double')
            max_param.set('name', 'max')
            value = ET.SubElement(max_param, 'value')
            value.text = f'{maxVal:.12f}'
        elif procid == ProcessType.DET:
            # Deterministic distribution
            # MATLAB: t = 1/sn.rates(i,r)
            distr = ET.SubElement(sub_param, 'subParameter')
            distr.set('classPath', 'jmt.engine.random.DeterministicDistr')
            distr.set('name', 'Deterministic')

            distr_par = ET.SubElement(sub_param, 'subParameter')
            distr_par.set('classPath', 'jmt.engine.random.DeterministicDistrPar')
            distr_par.set('name', 'distrPar')

            t_param = ET.SubElement(distr_par, 'subParameter')
            t_param.set('classPath', 'java.lang.Double')
            t_param.set('name', 't')
            value = ET.SubElement(t_param, 'value')
            value.text = f'{1.0 / rate:.12f}'
        else:
            # Default: Exponential distribution
            distr = ET.SubElement(sub_param, 'subParameter')
            distr.set('classPath', 'jmt.engine.random.Exponential')
            distr.set('name', 'Exponential')

            distr_par = ET.SubElement(sub_param, 'subParameter')
            distr_par.set('classPath', 'jmt.engine.random.ExponentialPar')
            distr_par.set('name', 'distrPar')

            lambda_param = ET.SubElement(distr_par, 'subParameter')
            lambda_param.set('classPath', 'java.lang.Double')
            lambda_param.set('name', 'lambda')
            value = ET.SubElement(lambda_param, 'value')
            value.text = f'{rate:.12f}'

    # JMT binds the Server constructor by parameter arity: delay-off/setup
    # must follow ServiceStrategy and precede PS weights, as in MATLAB writeJSIM.
    if sched != SchedStrategy.POLLING:
        _write_delayoff_strategy(server, node_idx, classnames, model)

    # Job parallelism and heterogeneous pools, only for the plain Server class:
    # SimLoader picks the constructor by the positional types of the parameters.
    if server.get('className') == 'Server':
        _write_server_pools(server, node_idx, sn, classnames)

    # PSStrategy (for PS/DPS/GPS and priority variants)
    if sched in (SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS,
                 SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO,
                 SchedStrategy.LPS):
        ps_strategy_param = ET.SubElement(server, 'parameter')
        ps_strategy_param.set('array', 'true')
        ps_strategy_param.set('classPath', 'jmt.engine.NetStrategies.PSStrategy')
        ps_strategy_param.set('name', 'PSStrategy')

        for r in range(K):
            ref_class = ET.SubElement(ps_strategy_param, 'refClass')
            ref_class.text = classnames[r]

            sub_param = ET.SubElement(ps_strategy_param, 'subParameter')
            if sched == SchedStrategy.PS:
                sub_param.set('classPath', 'jmt.engine.NetStrategies.PSStrategies.EPSStrategy')
                sub_param.set('name', 'EPSStrategy')
            elif sched == SchedStrategy.DPS:
                sub_param.set('classPath', 'jmt.engine.NetStrategies.PSStrategies.DPSStrategy')
                sub_param.set('name', 'DPSStrategy')
            elif sched == SchedStrategy.GPS:
                sub_param.set('classPath', 'jmt.engine.NetStrategies.PSStrategies.GPSStrategy')
                sub_param.set('name', 'GPSStrategy')
            elif sched == SchedStrategy.PSPRIO:
                sub_param.set('classPath', 'jmt.engine.NetStrategies.PSStrategies.EPSStrategyPriority')
                sub_param.set('name', 'EPSStrategyPriority')
            elif sched == SchedStrategy.DPSPRIO:
                sub_param.set('classPath', 'jmt.engine.NetStrategies.PSStrategies.DPSStrategyPriority')
                sub_param.set('name', 'DPSStrategyPriority')
            elif sched == SchedStrategy.GPSPRIO:
                sub_param.set('classPath', 'jmt.engine.NetStrategies.PSStrategies.GPSStrategyPriority')
                sub_param.set('name', 'GPSStrategyPriority')
            elif sched == SchedStrategy.LPS:
                sub_param.set('classPath', 'jmt.engine.NetStrategies.PSStrategies.EPSStrategy')
                sub_param.set('name', 'EPSStrategy')

    # Service weights (required for PSServer - PS/DPS/GPS and priority variants)
    if sched in (SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS,
                 SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO,
                 SchedStrategy.LPS):
        weights_param = ET.SubElement(server, 'parameter')
        weights_param.set('array', 'true')
        weights_param.set('classPath', 'java.lang.Double')
        weights_param.set('name', 'serviceWeights')

        for r in range(K):
            ref_class = ET.SubElement(weights_param, 'refClass')
            ref_class.text = classnames[r]

            sub_param = ET.SubElement(weights_param, 'subParameter')
            sub_param.set('classPath', 'java.lang.Double')
            sub_param.set('name', 'serviceWeight')
            value = ET.SubElement(sub_param, 'value')
            # Get weight from schedparam (for DPS/GPS and priority variants), default to 1 for PS
            weight = 1.0
            if sched in (SchedStrategy.DPS, SchedStrategy.GPS, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO):
                if hasattr(sn, 'schedparam') and sn.schedparam is not None:
                    if ist < sn.schedparam.shape[0] and r < sn.schedparam.shape[1]:
                        w = sn.schedparam[ist, r]
                        if not np.isnan(w) and w > 0:
                            weight = w
            value.text = _matlab_num2str(weight)

    # Switchover strategy for polling queues
    if sched == SchedStrategy.POLLING:
        _write_switchover_strategy(server, node_idx, sn, classnames)
    else:
        # Check if switchover times are defined for non-polling queues and warn
        has_switchover = False
        if sn.nodeparam and node_idx in sn.nodeparam:
            nodeparam = sn.nodeparam[node_idx]
            if isinstance(nodeparam, dict):
                for r in range(K):
                    if r in nodeparam and isinstance(nodeparam[r], dict):
                        if 'switchoverTime' in nodeparam[r] and nodeparam[r]['switchoverTime']:
                            has_switchover = True
                            break
        if has_switchover:
            import warnings
            node_name = sn.nodenames[node_idx] if hasattr(sn, 'nodenames') and node_idx < len(sn.nodenames) else f"node {node_idx}"
            warnings.warn(f"JMT does not support switchover times for non-polling queues. "
                         f"Switchover times will be ignored for {node_name}.")

    # Router section
    router = ET.SubElement(node_elem, 'section')
    router.set('className', 'Router')
    _write_routing_strategy(router, node_idx, sn, classnames, cs_node_names)


def _write_router_node(node_elem: ET.Element, node_idx: int, sn: NetworkStruct, classnames: List[str],
                       cs_node_names: Optional[Dict[Tuple[int, int], str]] = None):
    """Write router node section.

    Router nodes in JMT need:
    1. Queue (input section - buffer)
    2. ServiceTunnel (middle section - pass-through)
    3. Router (output section with routing strategy)
    """
    K = sn.nclasses

    # 1. Queue section (input buffer)
    queue = ET.SubElement(node_elem, 'section')
    queue.set('className', 'Queue')

    size_param = ET.SubElement(queue, 'parameter')
    size_param.set('classPath', 'java.lang.Integer')
    size_param.set('name', 'size')
    value = ET.SubElement(size_param, 'value')
    value.text = '-1'  # Infinite capacity

    # Drop strategies
    drop_strategy = ET.SubElement(queue, 'parameter')
    drop_strategy.set('array', 'true')
    drop_strategy.set('classPath', 'java.lang.String')
    drop_strategy.set('name', 'dropStrategies')

    for r in range(K):
        ref_class = ET.SubElement(drop_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(drop_strategy, 'subParameter')
        sub_param.set('classPath', 'java.lang.String')
        sub_param.set('name', 'dropStrategy')
        value = ET.SubElement(sub_param, 'value')
        value.text = 'drop'

    # Queue get strategy (FCFS)
    strategy_param = ET.SubElement(queue, 'parameter')
    strategy_param.set('classPath', 'jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy')
    strategy_param.set('name', 'FCFSstrategy')

    # Queue put strategy
    put_strategy = ET.SubElement(queue, 'parameter')
    put_strategy.set('array', 'true')
    put_strategy.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategy')
    put_strategy.set('name', 'QueuePutStrategy')

    for r in range(K):
        ref_class = ET.SubElement(put_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(put_strategy, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy')
        sub_param.set('name', 'TailStrategy')

    # 2. ServiceTunnel section (pass-through)
    tunnel = ET.SubElement(node_elem, 'section')
    tunnel.set('className', 'ServiceTunnel')

    # 3. Router section
    router = ET.SubElement(node_elem, 'section')
    router.set('className', 'Router')
    _write_routing_strategy(router, node_idx, sn, classnames, cs_node_names)


def _write_classswitch_node(node_elem: ET.Element, node_idx: int, sn: NetworkStruct, classnames: List[str]):
    """Write ClassSwitch node section.

    ClassSwitch nodes in JMT need:
    1. Queue (input section - buffer)
    2. ClassSwitch (middle section - class switching matrix)
    3. Router (output section with routing strategy)
    """
    K = sn.nclasses

    # 1. Queue section (input buffer)
    queue = ET.SubElement(node_elem, 'section')
    queue.set('className', 'Queue')

    size_param = ET.SubElement(queue, 'parameter')
    size_param.set('classPath', 'java.lang.Integer')
    size_param.set('name', 'size')
    value = ET.SubElement(size_param, 'value')
    value.text = '-1'  # Infinite capacity

    # Drop strategies
    drop_strategy = ET.SubElement(queue, 'parameter')
    drop_strategy.set('array', 'true')
    drop_strategy.set('classPath', 'java.lang.String')
    drop_strategy.set('name', 'dropStrategies')

    for r in range(K):
        ref_class = ET.SubElement(drop_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(drop_strategy, 'subParameter')
        sub_param.set('classPath', 'java.lang.String')
        sub_param.set('name', 'dropStrategy')
        value = ET.SubElement(sub_param, 'value')
        value.text = 'drop'

    # Queue get strategy (FCFS)
    strategy_param = ET.SubElement(queue, 'parameter')
    strategy_param.set('classPath', 'jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy')
    strategy_param.set('name', 'FCFSstrategy')

    # Queue put strategy
    put_strategy = ET.SubElement(queue, 'parameter')
    put_strategy.set('array', 'true')
    put_strategy.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategy')
    put_strategy.set('name', 'QueuePutStrategy')

    for r in range(K):
        ref_class = ET.SubElement(put_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(put_strategy, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy')
        sub_param.set('name', 'TailStrategy')

    # 2. ClassSwitch section (class switching matrix)
    cs_section = ET.SubElement(node_elem, 'section')
    cs_section.set('className', 'ClassSwitch')

    # Build the class switching matrix
    matrix_param = ET.SubElement(cs_section, 'parameter')
    matrix_param.set('array', 'true')
    matrix_param.set('classPath', 'java.lang.Object')
    matrix_param.set('name', 'matrix')

    # Get connections from this node to determine destination nodes
    if sn.connmatrix is not None:
        conn_i = sn.connmatrix[node_idx, :]
        jset = np.where(conn_i > 0)[0]
    else:
        jset = np.array([])

    for r in range(K):
        ref_class = ET.SubElement(matrix_param, 'refClass')
        ref_class.text = classnames[r]

        row_param = ET.SubElement(matrix_param, 'subParameter')
        row_param.set('array', 'true')
        row_param.set('classPath', 'java.lang.Float')
        row_param.set('name', 'row')

        for s in range(K):
            ref_class_col = ET.SubElement(row_param, 'refClass')
            ref_class_col.text = classnames[s]

            cell_param = ET.SubElement(row_param, 'subParameter')
            cell_param.set('classPath', 'java.lang.Float')
            cell_param.set('name', 'cell')
            cell_value = ET.SubElement(cell_param, 'value')

            # Calculate class switching probability from rtnodes
            val = 0.0
            if sn.rtnodes is not None and len(jset) > 0:
                for j in jset:
                    src_idx = node_idx * K + r
                    dst_idx = int(j) * K + s
                    if src_idx < sn.rtnodes.shape[0] and dst_idx < sn.rtnodes.shape[1]:
                        val += sn.rtnodes[src_idx, dst_idx]
            elif r == s:
                # Default: keep same class
                val = 1.0
            cell_value.text = f'{val:.12f}'

    # ClassSwitch nodes always use RAND strategy, matching MATLAB
    # saveRoutingStrategy.m:18-20.
    router = ET.SubElement(node_elem, 'section')
    router.set('className', 'Router')

    routing_param = ET.SubElement(router, 'parameter')
    routing_param.set('array', 'true')
    routing_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategy')
    routing_param.set('name', 'RoutingStrategy')

    for r in range(K):
        ref_class = ET.SubElement(routing_param, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(routing_param, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.RandomStrategy')
        sub_param.set('name', 'Random')


def _write_place_node(node_elem: ET.Element, node_idx: int, sn: NetworkStruct, classnames: List[str]):
    """Write Place node section for Petri nets.

    Place nodes in JMT use a Storage section.
    """
    K = sn.nclasses

    # Get station index for this node
    ist = -1
    if hasattr(sn, 'nodeToStation') and sn.nodeToStation is not None and node_idx < len(sn.nodeToStation):
        ist = int(sn.nodeToStation[node_idx])

    # Get total capacity from sn.cap (following MATLAB's saveTotalCapacity.m logic)
    total_cap = -1  # Default to infinite
    if ist >= 0 and hasattr(sn, 'cap') and sn.cap is not None and ist < len(sn.cap):
        cap_val = sn.cap[ist]
        if not np.isinf(cap_val):
            total_cap = int(cap_val)

    # Storage section
    storage = ET.SubElement(node_elem, 'section')
    storage.set('className', 'Storage')

    # Total capacity
    capacity_param = ET.SubElement(storage, 'parameter')
    capacity_param.set('classPath', 'java.lang.Integer')
    capacity_param.set('name', 'totalCapacity')
    value = ET.SubElement(capacity_param, 'value')
    value.text = str(total_cap)

    # Place capacities (per-class) - use classcap if available
    place_cap = ET.SubElement(storage, 'parameter')
    place_cap.set('array', 'true')
    place_cap.set('classPath', 'java.lang.Integer')
    place_cap.set('name', 'capacities')

    for r in range(K):
        ref_class = ET.SubElement(place_cap, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(place_cap, 'subParameter')
        sub_param.set('classPath', 'java.lang.Integer')
        sub_param.set('name', 'capacity')
        value = ET.SubElement(sub_param, 'value')
        # Use per-class capacity if available, otherwise use total capacity
        class_cap = -1
        if ist >= 0 and hasattr(sn, 'classcap') and sn.classcap is not None:
            if ist < sn.classcap.shape[0] and r < sn.classcap.shape[1]:
                cc_val = sn.classcap[ist, r]
                if not np.isinf(cc_val):
                    class_cap = int(cc_val)
        value.text = str(class_cap)

    # Drop rules - use 'waiting queue' for SPNs (matches MATLAB DropStrategy.WAITQ)
    drop_rules = ET.SubElement(storage, 'parameter')
    drop_rules.set('array', 'true')
    drop_rules.set('classPath', 'java.lang.String')
    drop_rules.set('name', 'dropRules')

    for r in range(K):
        ref_class = ET.SubElement(drop_rules, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(drop_rules, 'subParameter')
        sub_param.set('classPath', 'java.lang.String')
        sub_param.set('name', 'dropRule')
        value = ET.SubElement(sub_param, 'value')
        value.text = 'waiting queue'

    # Get strategy (FCFS)
    get_strategy = ET.SubElement(storage, 'parameter')
    get_strategy.set('classPath', 'jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy')
    get_strategy.set('name', 'FCFSstrategy')

    # Put strategies - use QueuePutStrategy name to match MATLAB
    put_strategy = ET.SubElement(storage, 'parameter')
    put_strategy.set('array', 'true')
    put_strategy.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategy')
    put_strategy.set('name', 'QueuePutStrategy')

    for r in range(K):
        ref_class = ET.SubElement(put_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(put_strategy, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy')
        sub_param.set('name', 'TailStrategy')

    # 2. ServiceTunnel section (pass-through - tokens don't need service)
    tunnel = ET.SubElement(node_elem, 'section')
    tunnel.set('className', 'ServiceTunnel')

    # 3. Linkage section (for Places, not Router - JMT handles routing via connections)
    linkage = ET.SubElement(node_elem, 'section')
    linkage.set('className', 'Linkage')


def _write_transition_node(node_elem: ET.Element, node_idx: int, sn: NetworkStruct, classnames: List[str]):
    """Write Transition node section for Petri nets.

    Transition nodes in JMT have three sections:
    1. Enabling (enabling and inhibiting conditions)
    2. Timing (mode names, servers, timing strategies)
    3. Firing (firing outcomes)
    """
    K = sn.nclasses
    nodenames = sn.nodenames if sn.nodenames else [f'Node{i+1}' for i in range(sn.nnodes)]

    # Get transition parameters from nodeparam
    trans_param = None
    if hasattr(sn, 'nodeparam') and sn.nodeparam is not None:
        trans_param = sn.nodeparam.get(node_idx)

    # Defaults if no nodeparam
    nmodes = 1
    modenames = ['Mode1']
    firing_prio = [1]
    fire_weight = [1.0]
    nmodeservers = [1]
    enabling = [np.zeros((sn.nnodes, K))]
    inhibiting = [np.full((sn.nnodes, K), np.inf)]
    firing = [np.zeros((sn.nnodes, K))]

    if trans_param is not None:
        nmodes = getattr(trans_param, 'nmodes', 1)
        modenames = getattr(trans_param, 'modenames', ['Mode1'])
        firing_prio = getattr(trans_param, 'firingprio', [1])
        fire_weight = getattr(trans_param, 'fireweight', [1.0])
        nmodeservers = getattr(trans_param, 'nmodeservers', np.array([1]))
        enabling = getattr(trans_param, 'enabling', enabling)
        inhibiting = getattr(trans_param, 'inhibiting', inhibiting)
        firing = getattr(trans_param, 'firing', firing)

    # 1. Enabling section
    enabling_section = ET.SubElement(node_elem, 'section')
    enabling_section.set('className', 'Enabling')

    # Enabling conditions - uses TransitionMatrix structure
    enabling_param = ET.SubElement(enabling_section, 'parameter')
    enabling_param.set('array', 'true')
    enabling_param.set('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix')
    enabling_param.set('name', 'enablingConditions')

    for m in range(nmodes):
        # Each mode has a TransitionMatrix
        mode_matrix = ET.SubElement(enabling_param, 'subParameter')
        mode_matrix.set('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix')
        mode_matrix.set('name', 'enablingCondition')

        # enablingVectors array
        vectors = ET.SubElement(mode_matrix, 'subParameter')
        vectors.set('array', 'true')
        vectors.set('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionVector')
        vectors.set('name', 'enablingVectors')

        for k in range(sn.nnodes):
            if sn.nodetype[k] != NodeType.PLACE:
                continue

            # Check if this place has relevant entries
            has_relevant = False
            for r in range(K):
                en_val = enabling[m][k, r] if m < len(enabling) else 0
                in_val = inhibiting[m][k, r] if m < len(inhibiting) else np.inf
                if (not np.isinf(en_val) and en_val > 0) or (not np.isinf(in_val) and in_val > 0):
                    has_relevant = True
                    break

            if not has_relevant:
                continue

            # Create vector for this place
            vector = ET.SubElement(vectors, 'subParameter')
            vector.set('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionVector')
            vector.set('name', 'enablingVector')

            # Station name
            station_name = ET.SubElement(vector, 'subParameter')
            station_name.set('classPath', 'java.lang.String')
            station_name.set('name', 'stationName')
            value = ET.SubElement(station_name, 'value')
            value.text = nodenames[k]

            # Enabling entries array
            entries = ET.SubElement(vector, 'subParameter')
            entries.set('array', 'true')
            entries.set('classPath', 'java.lang.Integer')
            entries.set('name', 'enablingEntries')

            for r in range(K):
                ref_class = ET.SubElement(entries, 'refClass')
                ref_class.text = classnames[r]

                entry = ET.SubElement(entries, 'subParameter')
                entry.set('classPath', 'java.lang.Integer')
                entry.set('name', 'enablingEntry')
                val = ET.SubElement(entry, 'value')
                en_val = enabling[m][k, r] if m < len(enabling) else 0
                val.text = '-1' if np.isinf(en_val) else str(int(en_val))

    # Inhibiting conditions - MATLAB writes vectors for ALL input places (never skips)
    # and uses '0' for infinite inhibiting values (not '-1' like enabling)
    inhibiting_param = ET.SubElement(enabling_section, 'parameter')
    inhibiting_param.set('array', 'true')
    inhibiting_param.set('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix')
    inhibiting_param.set('name', 'inhibitingConditions')

    # Get input places (nodes connected TO this transition)
    input_places = []
    if sn.connmatrix is not None:
        for k in range(sn.nnodes):
            if sn.connmatrix[k, node_idx] > 0 and sn.nodetype[k] == NodeType.PLACE:
                input_places.append(k)

    for m in range(nmodes):
        mode_matrix = ET.SubElement(inhibiting_param, 'subParameter')
        mode_matrix.set('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix')
        mode_matrix.set('name', 'inhibitingCondition')

        vectors = ET.SubElement(mode_matrix, 'subParameter')
        vectors.set('array', 'true')
        vectors.set('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionVector')
        vectors.set('name', 'inhibitingVectors')

        # Write vectors for ALL input places (no skipping)
        for k in input_places:
            vector = ET.SubElement(vectors, 'subParameter')
            vector.set('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionVector')
            vector.set('name', 'inhibitingVector')

            station_name = ET.SubElement(vector, 'subParameter')
            station_name.set('classPath', 'java.lang.String')
            station_name.set('name', 'stationName')
            value = ET.SubElement(station_name, 'value')
            value.text = nodenames[k]

            entries = ET.SubElement(vector, 'subParameter')
            entries.set('array', 'true')
            entries.set('classPath', 'java.lang.Integer')
            entries.set('name', 'inhibitingEntries')

            for r in range(K):
                ref_class = ET.SubElement(entries, 'refClass')
                ref_class.text = classnames[r]

                entry = ET.SubElement(entries, 'subParameter')
                entry.set('classPath', 'java.lang.Integer')
                entry.set('name', 'inhibitingEntry')
                val = ET.SubElement(entry, 'value')
                in_val = inhibiting[m][k, r] if m < len(inhibiting) else np.inf
                # Use '0' for infinite inhibiting (matches MATLAB), unlike enabling which uses '-1'
                val.text = '0' if np.isinf(in_val) else str(int(in_val))

    # 2. Timing section
    timing_section = ET.SubElement(node_elem, 'section')
    timing_section.set('className', 'Timing')

    # Mode names
    modenames_param = ET.SubElement(timing_section, 'parameter')
    modenames_param.set('array', 'true')
    modenames_param.set('classPath', 'java.lang.String')
    modenames_param.set('name', 'modeNames')

    for m in range(nmodes):
        sub_param = ET.SubElement(modenames_param, 'subParameter')
        sub_param.set('classPath', 'java.lang.String')
        sub_param.set('name', 'modeName')
        value = ET.SubElement(sub_param, 'value')
        value.text = modenames[m]

    # Number of servers
    servers_param = ET.SubElement(timing_section, 'parameter')
    servers_param.set('array', 'true')
    servers_param.set('classPath', 'java.lang.Integer')
    servers_param.set('name', 'numbersOfServers')

    for m in range(nmodes):
        sub_param = ET.SubElement(servers_param, 'subParameter')
        sub_param.set('classPath', 'java.lang.Integer')
        sub_param.set('name', 'numberOfServers')
        value = ET.SubElement(sub_param, 'value')
        nservers_raw = nmodeservers[m] if m < len(nmodeservers) else 1
        if np.isinf(nservers_raw) or nservers_raw >= 1000000:
            value.text = '-1'  # JMT uses -1 for infinite servers
        else:
            value.text = str(int(nservers_raw))

    # Timing strategies
    timing_strat_param = ET.SubElement(timing_section, 'parameter')
    timing_strat_param.set('array', 'true')
    timing_strat_param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategy')
    timing_strat_param.set('name', 'timingStrategies')

    for m in range(nmodes):
        dist = None
        if trans_param is not None and hasattr(trans_param, 'distributions'):
            dists = trans_param.distributions
            if m < len(dists) and dists[m] is not None:
                dist = dists[m]

        # see _kb/06-solver-catalog.md (JMT: "RNG consumption, timing
        # authority, SPN wiring, servers")
        is_immediate = False
        if trans_param is not None and hasattr(trans_param, 'timingstrategies'):
            tss = trans_param.timingstrategies
            if m < len(tss):
                ts = tss[m]
                is_immediate = str(getattr(ts, 'name', ts)).upper() == 'IMMEDIATE'
        if dist is not None and dist.isImmediate():
            is_immediate = True

        if dist is not None and not is_immediate:
            sub_param = ET.SubElement(timing_strat_param, 'subParameter')
            sub_param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy')
            sub_param.set('name', 'timingStrategy')
            _write_distribution_param(sub_param, dist)
        else:
            sub_param = ET.SubElement(timing_strat_param, 'subParameter')
            sub_param.set('classPath', 'jmt.engine.NetStrategies.ServiceStrategies.ZeroServiceTimeStrategy')
            sub_param.set('name', 'ZeroServiceTimeStrategy')

    # Firing priorities
    prio_param = ET.SubElement(timing_section, 'parameter')
    prio_param.set('array', 'true')
    prio_param.set('classPath', 'java.lang.Integer')
    prio_param.set('name', 'firingPriorities')

    for m in range(nmodes):
        sub_param = ET.SubElement(prio_param, 'subParameter')
        sub_param.set('classPath', 'java.lang.Integer')
        sub_param.set('name', 'firingPriority')
        value = ET.SubElement(sub_param, 'value')
        value.text = str(int(firing_prio[m]) if m < len(firing_prio) else 1)

    # Firing weights
    weight_param = ET.SubElement(timing_section, 'parameter')
    weight_param.set('array', 'true')
    weight_param.set('classPath', 'java.lang.Double')
    weight_param.set('name', 'firingWeights')

    for m in range(nmodes):
        sub_param = ET.SubElement(weight_param, 'subParameter')
        sub_param.set('classPath', 'java.lang.Double')
        sub_param.set('name', 'firingWeight')
        value = ET.SubElement(sub_param, 'value')
        value.text = str(float(fire_weight[m]) if m < len(fire_weight) else 1.0)

    # 3. Firing section
    firing_section = ET.SubElement(node_elem, 'section')
    firing_section.set('className', 'Firing')

    # see _kb/06-solver-catalog.md (JMT: "RNG consumption, timing authority,
    # SPN wiring, servers") -- firing outcomes for all output Places
    firing_param = ET.SubElement(firing_section, 'parameter')
    firing_param.set('array', 'true')
    firing_param.set('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix')
    firing_param.set('name', 'firingOutcomes')

    # All output nodes connected from this transition (MATLAB does not filter).
    output_nodes = []
    if sn.connmatrix is not None:
        for k in range(sn.nnodes):
            if sn.connmatrix[node_idx, k] > 0:
                # Include Places and Sink (matching MATLAB behavior)
                if sn.nodetype[k] in (NodeType.PLACE, NodeType.SINK):
                    output_nodes.append(k)

    for m in range(nmodes):
        mode_matrix = ET.SubElement(firing_param, 'subParameter')
        mode_matrix.set('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionMatrix')
        mode_matrix.set('name', 'firingOutcome')

        vectors = ET.SubElement(mode_matrix, 'subParameter')
        vectors.set('array', 'true')
        vectors.set('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionVector')
        vectors.set('name', 'firingVectors')

        # Write vectors for all output nodes
        for k in output_nodes:
            vector = ET.SubElement(vectors, 'subParameter')
            vector.set('classPath', 'jmt.engine.NetStrategies.TransitionUtilities.TransitionVector')
            vector.set('name', 'firingVector')

            station_name = ET.SubElement(vector, 'subParameter')
            station_name.set('classPath', 'java.lang.String')
            station_name.set('name', 'stationName')
            value = ET.SubElement(station_name, 'value')
            value.text = nodenames[k]

            entries = ET.SubElement(vector, 'subParameter')
            entries.set('array', 'true')
            entries.set('classPath', 'java.lang.Integer')
            entries.set('name', 'firingEntries')

            for r in range(K):
                ref_class = ET.SubElement(entries, 'refClass')
                ref_class.text = classnames[r]

                entry = ET.SubElement(entries, 'subParameter')
                entry.set('classPath', 'java.lang.Integer')
                entry.set('name', 'firingEntry')
                val = ET.SubElement(entry, 'value')
                fire_val = firing[m][k, r] if m < len(firing) else 0
                val.text = str(int(fire_val))


def _write_logger_node(node_elem: ET.Element, node_idx: int, sn: NetworkStruct,
                       classnames: List[str], cs_node_names: Optional[Dict[Tuple[int, int], str]] = None):
    """
    Write Logger node for job arrival/departure logging.

    Logger nodes are used by getCdfRespT to collect response time samples
    via transient simulation with logging.

    Logger has three sections:
    1. Queue (input buffer)
    2. LogTunnel (server that logs)
    3. Router (output)
    """
    K = sn.nclasses

    # Get logger parameters from nodeparam
    logger_param = None
    if sn.nodeparam is not None and node_idx in sn.nodeparam:
        logger_param = sn.nodeparam[node_idx]

    # Default values if no params
    file_name = 'log.csv'
    file_path = '/tmp/'
    start_time = 'false'
    logger_name = 'false'
    timestamp = 'true'
    job_id = 'true'
    job_class = 'true'
    time_same_class = 'false'
    time_any_class = 'false'

    if logger_param is not None:
        if hasattr(logger_param, 'fileName'):
            file_name = logger_param.fileName
        if hasattr(logger_param, 'filePath'):
            file_path = logger_param.filePath
        if hasattr(logger_param, 'startTime'):
            start_time = 'true' if logger_param.startTime else 'false'
        if hasattr(logger_param, 'loggerName'):
            logger_name = 'true' if logger_param.loggerName else 'false'
        if hasattr(logger_param, 'timestamp'):
            timestamp = 'true' if logger_param.timestamp else 'false'
        if hasattr(logger_param, 'jobID'):
            job_id = 'true' if logger_param.jobID else 'false'
        if hasattr(logger_param, 'jobClass'):
            job_class = 'true' if logger_param.jobClass else 'false'
        if hasattr(logger_param, 'timeSameClass'):
            time_same_class = 'true' if logger_param.timeSameClass else 'false'
        if hasattr(logger_param, 'timeAnyClass'):
            time_any_class = 'true' if logger_param.timeAnyClass else 'false'

    # Ensure path ends with separator
    import os
    if not file_path.endswith(os.sep):
        file_path = file_path + os.sep

    # 1. Queue section (input buffer) - like a zero-capacity queue
    queue_section = ET.SubElement(node_elem, 'section')
    queue_section.set('className', 'Queue')

    # Size parameter (infinite capacity for logger)
    size_param = ET.SubElement(queue_section, 'parameter')
    size_param.set('classPath', 'java.lang.Integer')
    size_param.set('name', 'size')
    size_val = ET.SubElement(size_param, 'value')
    size_val.text = '-1'  # Infinite capacity

    # Drop strategy
    drop_strat = ET.SubElement(queue_section, 'parameter')
    drop_strat.set('array', 'true')
    drop_strat.set('classPath', 'java.lang.String')
    drop_strat.set('name', 'dropStrategies')

    for r in range(K):
        ref_class = ET.SubElement(drop_strat, 'refClass')
        ref_class.text = classnames[r]
        sub_param = ET.SubElement(drop_strat, 'subParameter')
        sub_param.set('classPath', 'java.lang.String')
        sub_param.set('name', 'dropStrategy')
        val = ET.SubElement(sub_param, 'value')
        val.text = 'drop'

    # Get strategy (FCFS) - no subParameter, classPath directly references FCFSstrategy
    get_strat = ET.SubElement(queue_section, 'parameter')
    get_strat.set('classPath', 'jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy')
    get_strat.set('name', 'FCFSstrategy')

    # Put strategies
    put_strat = ET.SubElement(queue_section, 'parameter')
    put_strat.set('array', 'true')
    put_strat.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategy')
    put_strat.set('name', 'QueuePutStrategy')

    for r in range(K):
        ref_class = ET.SubElement(put_strat, 'refClass')
        ref_class.text = classnames[r]
        sub_param = ET.SubElement(put_strat, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy')
        sub_param.set('name', 'TailStrategy')

    # 2. LogTunnel section
    section = ET.SubElement(node_elem, 'section')
    section.set('className', 'LogTunnel')

    # Logger parameters following MATLAB's saveLogTunnel.m
    params = [
        ('logfileName', 'java.lang.String', file_name),
        ('logfilePath', 'java.lang.String', file_path),
        ('logExecTimestamp', 'java.lang.Boolean', start_time),
        ('logLoggerName', 'java.lang.Boolean', logger_name),
        ('logTimeStamp', 'java.lang.Boolean', timestamp),
        ('logJobID', 'java.lang.Boolean', job_id),
        ('logJobClass', 'java.lang.Boolean', job_class),
        ('logTimeSameClass', 'java.lang.Boolean', time_same_class),
        ('logTimeAnyClass', 'java.lang.Boolean', time_any_class),
        ('numClasses', 'java.lang.Integer', str(K)),
    ]

    for param_name, class_path, param_value in params:
        param = ET.SubElement(section, 'parameter')
        param.set('classPath', class_path)
        param.set('name', param_name)
        value = ET.SubElement(param, 'value')
        value.text = str(param_value)

    # see _kb/06-solver-catalog.md (JMT: "RNG consumption, timing authority,
    # SPN wiring, servers")
    router = ET.SubElement(node_elem, 'section')
    router.set('className', 'Router')

    routing_param = ET.SubElement(router, 'parameter')
    routing_param.set('array', 'true')
    routing_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategy')
    routing_param.set('name', 'RoutingStrategy')

    # Find target nodes from connection matrix
    nodenames = sn.nodenames if sn.nodenames else [f'Node{i+1}' for i in range(sn.nnodes)]
    target_names = []
    if hasattr(sn, 'connmatrix') and sn.connmatrix is not None:
        for j in range(sn.connmatrix.shape[1]):
            if sn.connmatrix[node_idx, j] > 0:
                target_names.append(nodenames[j])

    for r in range(K):
        ref_class = ET.SubElement(routing_param, 'refClass')
        ref_class.text = classnames[r]

        if target_names:
            # Use EmpiricalStrategy with explicit routing probabilities
            sub_param = ET.SubElement(routing_param, 'subParameter')
            sub_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.EmpiricalStrategy')
            sub_param.set('name', 'Probabilities')

            emp_array = ET.SubElement(sub_param, 'subParameter')
            emp_array.set('array', 'true')
            emp_array.set('classPath', 'jmt.engine.random.EmpiricalEntry')
            emp_array.set('name', 'EmpiricalEntryArray')

            # Distribute probability equally among targets
            prob = 1.0 / len(target_names)
            for tgt_name in target_names:
                emp_entry = ET.SubElement(emp_array, 'subParameter')
                emp_entry.set('classPath', 'jmt.engine.random.EmpiricalEntry')
                emp_entry.set('name', 'EmpiricalEntry')

                station_param = ET.SubElement(emp_entry, 'subParameter')
                station_param.set('classPath', 'java.lang.String')
                station_param.set('name', 'stationName')
                station_value = ET.SubElement(station_param, 'value')
                station_value.text = tgt_name

                prob_param = ET.SubElement(emp_entry, 'subParameter')
                prob_param.set('classPath', 'java.lang.Double')
                prob_param.set('name', 'probability')
                prob_value = ET.SubElement(prob_param, 'value')
                prob_value.text = f'{prob:.12f}'
        else:
            # Fallback to RandomStrategy if no connections found
            sub_param = ET.SubElement(routing_param, 'subParameter')
            sub_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.RandomStrategy')
            sub_param.set('name', 'Random')


def _write_distribution_param(parent: ET.Element, dist) -> None:
    """Write distribution parameters for timing strategy."""
    dist_name = dist._name if hasattr(dist, '_name') else type(dist).__name__

    if dist_name == 'Exp':
        distr = ET.SubElement(parent, 'subParameter')
        distr.set('classPath', 'jmt.engine.random.Exponential')
        distr.set('name', 'Exponential')

        distr_par = ET.SubElement(parent, 'subParameter')
        distr_par.set('classPath', 'jmt.engine.random.ExponentialPar')
        distr_par.set('name', 'distrPar')

        lambda_param = ET.SubElement(distr_par, 'subParameter')
        lambda_param.set('classPath', 'java.lang.Double')
        lambda_param.set('name', 'lambda')
        value = ET.SubElement(lambda_param, 'value')
        rate = dist.get_rate() if hasattr(dist, 'get_rate') else 1.0
        value.text = str(rate)
    elif dist_name == 'Erlang':
        distr = ET.SubElement(parent, 'subParameter')
        distr.set('classPath', 'jmt.engine.random.Erlang')
        distr.set('name', 'Erlang')

        distr_par = ET.SubElement(parent, 'subParameter')
        distr_par.set('classPath', 'jmt.engine.random.ErlangPar')
        distr_par.set('name', 'distrPar')

        # Get number of phases - try multiple methods
        phases = 1
        if hasattr(dist, '_phases'):
            phases = dist._phases
        elif hasattr(dist, 'getNumberOfPhases'):
            try:
                phases = dist.getNumberOfPhases()
            except NotImplementedError:
                pass

        # JMT Erlang uses alpha = phase_rate = phases/mean
        # Use _phase_rate directly if available, otherwise compute from mean
        if hasattr(dist, '_phase_rate'):
            alpha = dist._phase_rate
        else:
            mean = dist.get_mean() if hasattr(dist, 'get_mean') else 1.0
            alpha = phases / mean if mean > 0 else 1.0

        alpha_param = ET.SubElement(distr_par, 'subParameter')
        alpha_param.set('classPath', 'java.lang.Double')
        alpha_param.set('name', 'alpha')
        value = ET.SubElement(alpha_param, 'value')
        value.text = f'{alpha:.12f}'

        r_param = ET.SubElement(distr_par, 'subParameter')
        r_param.set('classPath', 'java.lang.Long')
        r_param.set('name', 'r')
        value = ET.SubElement(r_param, 'value')
        value.text = str(phases)
    elif dist_name == 'HyperExp':
        # An n>2 HyperExp is emitted as the equivalent APH, since JMT's
        # HyperExpPar is 2-phase only.
        probs = np.array([0.5, 0.5])
        rates = np.array([1.0, 1.0])
        if hasattr(dist, '_probs') and hasattr(dist, '_rates'):
            probs = np.asarray(dist._probs, dtype=np.float64).flatten()
            rates = np.asarray(dist._rates, dtype=np.float64).flatten()
        _write_hyperexp_distribution(parent, probs, rates)
    elif dist_name in ('Coxian', 'PH', 'APH'):
        # Phase-Type distributions (Coxian, PH, APH) use PhaseTypeDistr format
        # Get alpha (initial probability vector) and T (sub-generator matrix)
        alpha = None
        T = None
        if hasattr(dist, 'alpha'):
            alpha = np.asarray(dist.alpha, dtype=np.float64)
        elif hasattr(dist, '_alpha'):
            alpha = np.asarray(dist._alpha, dtype=np.float64)
        if hasattr(dist, 'T'):
            T = np.asarray(dist.T, dtype=np.float64)
        elif hasattr(dist, '_T'):
            T = np.asarray(dist._T, dtype=np.float64)

        if alpha is not None and T is not None:
            # Ensure alpha is 1D and T is 2D
            alpha = alpha.flatten()
            if T.ndim == 1:
                T = T.reshape(1, -1)
            n_phases = len(alpha)

            # Write distribution element
            distr = ET.SubElement(parent, 'subParameter')
            distr.set('classPath', 'jmt.engine.random.PhaseTypeDistr')
            distr.set('name', 'Phase-Type')

            # Write parameter element
            distr_par = ET.SubElement(parent, 'subParameter')
            distr_par.set('classPath', 'jmt.engine.random.PhaseTypePar')
            distr_par.set('name', 'distrPar')

            # Write alpha (initial probability vector)
            alpha_param = ET.SubElement(distr_par, 'subParameter')
            alpha_param.set('array', 'true')
            alpha_param.set('classPath', 'java.lang.Object')
            alpha_param.set('name', 'alpha')

            alpha_vector = ET.SubElement(alpha_param, 'subParameter')
            alpha_vector.set('array', 'true')
            alpha_vector.set('classPath', 'java.lang.Object')
            alpha_vector.set('name', 'vector')

            for i in range(n_phases):
                entry = ET.SubElement(alpha_vector, 'subParameter')
                entry.set('classPath', 'java.lang.Double')
                entry.set('name', 'entry')
                value = ET.SubElement(entry, 'value')
                value.text = f'{alpha[i]:.12f}'

            # Write T (sub-generator matrix)
            T_param = ET.SubElement(distr_par, 'subParameter')
            T_param.set('array', 'true')
            T_param.set('classPath', 'java.lang.Object')
            T_param.set('name', 'T')

            for i in range(n_phases):
                row_vector = ET.SubElement(T_param, 'subParameter')
                row_vector.set('array', 'true')
                row_vector.set('classPath', 'java.lang.Object')
                row_vector.set('name', 'vector')

                for j in range(n_phases):
                    entry = ET.SubElement(row_vector, 'subParameter')
                    entry.set('classPath', 'java.lang.Double')
                    entry.set('name', 'entry')
                    value = ET.SubElement(entry, 'value')
                    value.text = f'{T[i, j]:.12f}'
        else:
            # Fallback to exponential if PH parameters not available
            distr = ET.SubElement(parent, 'subParameter')
            distr.set('classPath', 'jmt.engine.random.Exponential')
            distr.set('name', 'Exponential')

            distr_par = ET.SubElement(parent, 'subParameter')
            distr_par.set('classPath', 'jmt.engine.random.ExponentialPar')
            distr_par.set('name', 'distrPar')

            lambda_param = ET.SubElement(distr_par, 'subParameter')
            lambda_param.set('classPath', 'java.lang.Double')
            lambda_param.set('name', 'lambda')
            value = ET.SubElement(lambda_param, 'value')
            value.text = '1.0'
    elif dist_name == 'Pareto':
        distr = ET.SubElement(parent, 'subParameter')
        distr.set('classPath', 'jmt.engine.random.Pareto')
        distr.set('name', 'Pareto')

        distr_par = ET.SubElement(parent, 'subParameter')
        distr_par.set('classPath', 'jmt.engine.random.ParetoPar')
        distr_par.set('name', 'distrPar')

        # alpha (shape) / k (scale) read directly from the distribution object.
        alpha = 3.0  # default shape
        k = 1.0  # default scale
        if hasattr(dist, 'alpha'):
            alpha = float(dist.alpha)
        elif hasattr(dist, '_alpha'):
            alpha = float(dist._alpha)
        if hasattr(dist, 'scale'):
            k = float(dist.scale)
        elif hasattr(dist, '_scale'):
            k = float(dist._scale)

        alpha_param = ET.SubElement(distr_par, 'subParameter')
        alpha_param.set('classPath', 'java.lang.Double')
        alpha_param.set('name', 'alpha')
        value = ET.SubElement(alpha_param, 'value')
        value.text = f'{alpha:.12f}'

        k_param = ET.SubElement(distr_par, 'subParameter')
        k_param.set('classPath', 'java.lang.Double')
        k_param.set('name', 'k')
        value = ET.SubElement(k_param, 'value')
        value.text = f'{k:.12f}'
    elif dist_name == 'Replayer':
        # Replayer distribution - reads service times from a trace file
        distr = ET.SubElement(parent, 'subParameter')
        distr.set('classPath', 'jmt.engine.random.Replayer')
        distr.set('name', 'Replayer')

        distr_par = ET.SubElement(parent, 'subParameter')
        distr_par.set('classPath', 'jmt.engine.random.ReplayerPar')
        distr_par.set('name', 'distrPar')

        # Get the file path from the distribution
        file_path = ''
        if hasattr(dist, 'file_path'):
            file_path = dist.file_path
        elif hasattr(dist, '_file_path'):
            file_path = dist._file_path
        elif hasattr(dist, 'filePath'):
            file_path = dist.filePath
        elif hasattr(dist, 'fileName'):
            file_path = dist.fileName

        file_param = ET.SubElement(distr_par, 'subParameter')
        file_param.set('classPath', 'java.lang.String')
        file_param.set('name', 'fileName')
        value = ET.SubElement(file_param, 'value')
        value.text = str(file_path)
    else:
        # Default to exponential with rate 1
        distr = ET.SubElement(parent, 'subParameter')
        distr.set('classPath', 'jmt.engine.random.Exponential')
        distr.set('name', 'Exponential')

        distr_par = ET.SubElement(parent, 'subParameter')
        distr_par.set('classPath', 'jmt.engine.random.ExponentialPar')
        distr_par.set('name', 'distrPar')

        lambda_param = ET.SubElement(distr_par, 'subParameter')
        lambda_param.set('classPath', 'java.lang.Double')
        lambda_param.set('name', 'lambda')
        value = ET.SubElement(lambda_param, 'value')
        value.text = '1.0'


def _hyperexp_probs_rates(proc, pie=None):
    """
    Recover the (probs, rates) pair of an n-phase HyperExp from an sn.proc entry.

    Two storage forms occur across the writers:
      - dict form:   {'probs': [...], 'rates': [...]}
      - MAP form:    (D0, D1) with D0 = -diag(rates) and D1[i,j] = rates[i]*probs[j]

    Args:
        proc: The sn.proc[ist][r] entry.
        pie: Optional sn.pie[ist][r] entry (entry probability vector).

    Returns:
        (probs, rates) as float ndarrays, or None if they cannot be recovered.
    """
    if proc is None:
        return None

    if isinstance(proc, dict):
        if 'probs' in proc and 'rates' in proc:
            probs = np.asarray(proc['probs'], dtype=np.float64).flatten()
            rates = np.asarray(proc['rates'], dtype=np.float64).flatten()
            if len(probs) == len(rates) and len(rates) > 0:
                return probs, rates
        return None

    if isinstance(proc, (list, tuple)) and len(proc) >= 2:
        D0 = np.asarray(proc[0], dtype=np.float64)
        D1 = np.asarray(proc[1], dtype=np.float64)
        if D0.ndim != 2 or D0.shape[0] != D0.shape[1] or D0.shape != D1.shape:
            return None
        rates = -np.diag(D0)
        if np.any(rates <= 0):
            return None
        # Prefer the explicit entry vector when available, otherwise recover the
        # mixing probabilities from D1, whose row i equals rates[i]*probs.
        if pie is not None:
            probs = np.asarray(pie, dtype=np.float64).flatten()
            if len(probs) == len(rates):
                return probs, rates
        probs = D1[0, :] / rates[0]
        return probs, rates

    return None


def _write_hyperexp_distribution(parent: ET.Element, probs, rates) -> None:
    """
    Write a HyperExp distribution to JMT XML.

    JMT's jmt.engine.random.HyperExpPar is 2-phase only. An n>2 phase HyperExp is
    therefore emitted as an equivalent APH (PhaseTypeDistr/PhaseTypePar) built from
    the very same phase-type representation (alpha = probs, T = -diag(rates)), so the
    emitted distribution is identical rather than a moment-matched refit. The n<=2
    case keeps emitting HyperExpPar unchanged.

    Args:
        parent: Parent XML element.
        probs: Mixing probability vector (length n).
        rates: Phase rate vector (length n).
    """
    distr = ET.SubElement(parent, 'subParameter')
    distr_par = ET.SubElement(parent, 'subParameter')
    _fill_hyperexp_distribution(distr, distr_par, probs, rates)


def _fill_hyperexp_distribution(distr: ET.Element, distr_par: ET.Element,
                                probs, rates) -> None:
    """
    Populate pre-created distribution/parameter elements with a HyperExp distribution,
    falling back to the equivalent APH when it has more than 2 phases.

    Args:
        distr: Pre-created distribution element.
        distr_par: Pre-created parameter element.
        probs: Mixing probability vector (length n).
        rates: Phase rate vector (length n).
    """
    probs = np.asarray(probs, dtype=np.float64).flatten()
    rates = np.asarray(rates, dtype=np.float64).flatten()
    n_phases = len(rates)

    if n_phases > 2:
        # Same PH as the HyperExp itself: alpha = probs, T = -diag(rates).
        _fill_phase_type_par(distr, distr_par, probs, -np.diag(rates))
        return

    p = float(probs[0]) if len(probs) > 0 else 0.5
    lambda1 = float(rates[0])
    lambda2 = float(rates[1]) if n_phases > 1 else float(rates[0])

    distr.set('classPath', 'jmt.engine.random.HyperExp')
    distr.set('name', 'Hyperexponential')

    distr_par.set('classPath', 'jmt.engine.random.HyperExpPar')
    distr_par.set('name', 'distrPar')

    # %.12f, not repr: MATLAB saveServiceStrategy.m and JAR SaveHandlers write
    # every JMT distribution parameter at twelve decimals. A full-precision
    # rate here is a DIFFERENT model from the one the other two codebases hand
    # to the same JMT.jar, which at a fixed seed is a different sample path.
    p_param = ET.SubElement(distr_par, 'subParameter')
    p_param.set('classPath', 'java.lang.Double')
    p_param.set('name', 'p')
    value = ET.SubElement(p_param, 'value')
    value.text = '%.12f' % p

    l1_param = ET.SubElement(distr_par, 'subParameter')
    l1_param.set('classPath', 'java.lang.Double')
    l1_param.set('name', 'lambda1')
    value = ET.SubElement(l1_param, 'value')
    value.text = '%.12f' % lambda1

    l2_param = ET.SubElement(distr_par, 'subParameter')
    l2_param.set('classPath', 'java.lang.Double')
    l2_param.set('name', 'lambda2')
    value = ET.SubElement(l2_param, 'value')
    value.text = '%.12f' % lambda2


def _write_phase_type_service_distribution(parent: ET.Element, sn: NetworkStruct, ist: int, r: int) -> None:
    """
    Write Phase-Type service distribution to JMT XML.

    This writes the PhaseTypeDistr format used by JMT for general phase-type
    distributions (PH, APH, Coxian, etc.).

    Args:
        parent: Parent XML element (serviceTimeStrategyNode)
        sn: NetworkStruct containing proc and pie fields
        ist: Station index
        r: Class index
    """
    # Get phase-type representation from proc and pie. sn.proc holds (D0, D1);
    # proc_to_ph returns the PH view of it, so this reader never has to tell
    # (D0, D1) from a legacy [alpha, T] by shape -- both are a pair of arrays
    # and the guess was wrong for every Markovian family since the storage form
    # was unified (see _kb/04-networkstruct.md).
    T = None
    alpha = None

    if hasattr(sn, 'proc') and sn.proc is not None:
        try:
            proc = sn.proc[ist][r]
            if proc is not None:
                from ...sn.proc_form import proc_to_ph
                alpha_candidate, T_candidate = proc_to_ph(proc)
                if T_candidate is not None:
                    T = np.asarray(T_candidate, dtype=np.float64)
                    alpha = np.asarray(alpha_candidate, dtype=np.float64)
                elif isinstance(proc, (list, tuple)) and len(proc) == 1:
                    # Single element - assume it's T
                    T = np.asarray(proc[0], dtype=np.float64)
        except (IndexError, TypeError, KeyError):
            pass

    # If alpha not found in proc, try pie field
    if alpha is None and hasattr(sn, 'pie') and sn.pie is not None:
        try:
            pie = sn.pie[ist][r]
            if pie is not None:
                alpha = np.asarray(pie, dtype=np.float64)
        except (IndexError, TypeError, KeyError):
            pass

    if T is None:
        # Fallback to exponential if phase-type not available
        distr = ET.SubElement(parent, 'subParameter')
        distr.set('classPath', 'jmt.engine.random.Exponential')
        distr.set('name', 'Exponential')

        distr_par = ET.SubElement(parent, 'subParameter')
        distr_par.set('classPath', 'jmt.engine.random.ExponentialPar')
        distr_par.set('name', 'distrPar')

        lambda_param = ET.SubElement(distr_par, 'subParameter')
        lambda_param.set('classPath', 'java.lang.Double')
        lambda_param.set('name', 'lambda')
        value = ET.SubElement(lambda_param, 'value')
        rate = sn.rates[ist, r] if sn.rates is not None and ist < sn.rates.shape[0] and r < sn.rates.shape[1] else 1.0
        value.text = str(rate)
        return

    # Ensure T is 2D matrix
    if T.ndim == 1:
        # Convert 1D array to 2D (single phase)
        T = T.reshape(1, -1) if len(T) > 1 else np.array([[T[0]]])
    elif T.ndim == 0:
        # Scalar - convert to 1x1 matrix
        T = np.array([[float(T)]])

    n_phases = T.shape[0]

    if alpha is None:
        alpha = np.zeros(n_phases)
        alpha[0] = 1.0
    else:
        # Ensure alpha is 1D
        alpha = np.asarray(alpha).flatten()
        # Ensure correct size
        if len(alpha) < n_phases:
            alpha_new = np.zeros(n_phases)
            alpha_new[:len(alpha)] = alpha
            alpha = alpha_new
        elif len(alpha) > n_phases:
            alpha = alpha[:n_phases]

    # Ensure alpha is positive
    alpha = np.abs(alpha)

    _write_phase_type_par(parent, alpha, T)


def _write_phase_type_par(parent: ET.Element, alpha: np.ndarray, T: np.ndarray) -> None:
    """
    Write a PhaseTypeDistr/PhaseTypePar element pair from an explicit (alpha, T) pair.

    Args:
        parent: Parent XML element.
        alpha: Initial probability vector (length n).
        T: Sub-generator matrix (n x n).
    """
    distr = ET.SubElement(parent, 'subParameter')
    distr_par = ET.SubElement(parent, 'subParameter')
    _fill_phase_type_par(distr, distr_par, alpha, T)


def _fill_phase_type_par(distr: ET.Element, distr_par: ET.Element,
                         alpha: np.ndarray, T: np.ndarray) -> None:
    """
    Populate pre-created distribution/parameter elements as a PhaseTypeDistr/PhaseTypePar
    pair from an explicit (alpha, T) pair.

    Args:
        distr: Pre-created distribution element.
        distr_par: Pre-created parameter element.
        alpha: Initial probability vector (length n).
        T: Sub-generator matrix (n x n).
    """
    alpha = np.asarray(alpha, dtype=np.float64).flatten()
    T = np.asarray(T, dtype=np.float64)
    n_phases = T.shape[0]

    # Write distribution element
    distr.set('classPath', 'jmt.engine.random.PhaseTypeDistr')
    distr.set('name', 'Phase-Type')

    # Write parameter element
    distr_par.set('classPath', 'jmt.engine.random.PhaseTypePar')
    distr_par.set('name', 'distrPar')

    # Write alpha (initial probability vector)
    alpha_param = ET.SubElement(distr_par, 'subParameter')
    alpha_param.set('array', 'true')
    alpha_param.set('classPath', 'java.lang.Object')
    alpha_param.set('name', 'alpha')

    alpha_vec = ET.SubElement(alpha_param, 'subParameter')
    alpha_vec.set('array', 'true')
    alpha_vec.set('classPath', 'java.lang.Object')
    alpha_vec.set('name', 'vector')

    for k in range(n_phases):
        entry = ET.SubElement(alpha_vec, 'subParameter')
        entry.set('classPath', 'java.lang.Double')
        entry.set('name', 'entry')
        value = ET.SubElement(entry, 'value')
        value.text = f'{alpha[k]:.12f}'

    # Write T matrix (sub-generator)
    t_param = ET.SubElement(distr_par, 'subParameter')
    t_param.set('array', 'true')
    t_param.set('classPath', 'java.lang.Object')
    t_param.set('name', 'T')

    for k in range(n_phases):
        row_vec = ET.SubElement(t_param, 'subParameter')
        row_vec.set('array', 'true')
        row_vec.set('classPath', 'java.lang.Object')
        row_vec.set('name', 'vector')

        for j in range(n_phases):
            entry = ET.SubElement(row_vec, 'subParameter')
            entry.set('classPath', 'java.lang.Double')
            entry.set('name', 'entry')
            value = ET.SubElement(entry, 'value')
            # MATLAB: if k==j, use -abs(T(k,j)), else use abs(T(k,j))
            if k == j:
                value.text = f'{-abs(T[k, j]):.12f}'
            else:
                value.text = f'{abs(T[k, j]):.12f}'


def _write_map_service_distribution(parent: ET.Element, sn: NetworkStruct, ist: int, r: int) -> None:
    """
    Write MAP/MMPP2 service distribution to JMT XML.

    This writes the MAPDistr format used by JMT for Markov Arrival Processes
    and Markov Modulated Poisson Processes.

    Args:
        parent: Parent XML element (serviceTimeStrategyNode)
        sn: NetworkStruct containing proc field with D0, D1 matrices
        ist: Station index
        r: Class index
    """
    # Get D0 and D1 matrices from proc
    D0 = None
    D1 = None

    if hasattr(sn, 'proc') and sn.proc is not None:
        try:
            proc = sn.proc[ist][r]
            if proc is not None and isinstance(proc, (list, tuple)) and len(proc) >= 2:
                D0 = np.asarray(proc[0], dtype=np.float64)
                D1 = np.asarray(proc[1], dtype=np.float64)
        except (IndexError, TypeError):
            pass

    if D0 is None or D1 is None:
        # Fallback to exponential if matrices not available
        distr = ET.SubElement(parent, 'subParameter')
        distr.set('classPath', 'jmt.engine.random.Exponential')
        distr.set('name', 'Exponential')

        distr_par = ET.SubElement(parent, 'subParameter')
        distr_par.set('classPath', 'jmt.engine.random.ExponentialPar')
        distr_par.set('name', 'distrPar')

        lambda_param = ET.SubElement(distr_par, 'subParameter')
        lambda_param.set('classPath', 'java.lang.Double')
        lambda_param.set('name', 'lambda')
        value = ET.SubElement(lambda_param, 'value')
        rate = sn.rates[ist, r] if sn.rates is not None and ist < sn.rates.shape[0] and r < sn.rates.shape[1] else 1.0
        value.text = str(rate)
        return

    # Get number of phases
    n_phases = D0.shape[0]

    # Write distribution element
    distr = ET.SubElement(parent, 'subParameter')
    distr.set('classPath', 'jmt.engine.random.MAPDistr')
    distr.set('name', 'Burst (MAP)')

    # Write parameter element
    distr_par = ET.SubElement(parent, 'subParameter')
    distr_par.set('classPath', 'jmt.engine.random.MAPPar')
    distr_par.set('name', 'distrPar')

    # Write D0 matrix
    d0_param = ET.SubElement(distr_par, 'subParameter')
    d0_param.set('array', 'true')
    d0_param.set('classPath', 'java.lang.Object')
    d0_param.set('name', 'D0')

    for k in range(n_phases):
        row_param = ET.SubElement(d0_param, 'subParameter')
        row_param.set('array', 'true')
        row_param.set('classPath', 'java.lang.Object')
        row_param.set('name', 'vector')

        for j in range(n_phases):
            entry_param = ET.SubElement(row_param, 'subParameter')
            entry_param.set('classPath', 'java.lang.Double')
            entry_param.set('name', 'entry')
            value = ET.SubElement(entry_param, 'value')
            value.text = f'{D0[k, j]:.12f}'

    # Write D1 matrix
    d1_param = ET.SubElement(distr_par, 'subParameter')
    d1_param.set('array', 'true')
    d1_param.set('classPath', 'java.lang.Object')
    d1_param.set('name', 'D1')

    for k in range(n_phases):
        row_param = ET.SubElement(d1_param, 'subParameter')
        row_param.set('array', 'true')
        row_param.set('classPath', 'java.lang.Object')
        row_param.set('name', 'vector')

        for j in range(n_phases):
            entry_param = ET.SubElement(row_param, 'subParameter')
            entry_param.set('classPath', 'java.lang.Double')
            entry_param.set('name', 'entry')
            value = ET.SubElement(entry_param, 'value')
            value.text = f'{D1[k, j]:.12f}'


def _write_fork_node(node_elem: ET.Element, node_idx: int, sn: NetworkStruct, classnames: List[str],
                     cs_node_names: Optional[Dict[Tuple[int, int], str]] = None):
    """Write Fork node section.

    Fork nodes in JMT have three sections:
    1. Queue (buffer)
    2. ServiceTunnel
    3. Fork (with ForkStrategy)
    """
    K = sn.nclasses

    # Get fanOut (tasks per link), default to 1
    fan_out = 1
    fan_out_link = None
    fan_out_prob = None
    fan_out_dist = None
    if hasattr(sn, 'nodeparam') and sn.nodeparam is not None:
        if node_idx in sn.nodeparam and sn.nodeparam[node_idx] is not None:
            if isinstance(sn.nodeparam[node_idx], dict):
                fan_out = sn.nodeparam[node_idx].get('fanOut', 1)
                fan_out_link = sn.nodeparam[node_idx].get('fanOutLink', None)
                fan_out_prob = sn.nodeparam[node_idx].get('fanOutProb', None)
                fan_out_dist = sn.nodeparam[node_idx].get('fanOutDist', None)

    # 1. Queue section (buffer)
    queue = ET.SubElement(node_elem, 'section')
    queue.set('className', 'Queue')

    size_param = ET.SubElement(queue, 'parameter')
    size_param.set('classPath', 'java.lang.Integer')
    size_param.set('name', 'size')
    value = ET.SubElement(size_param, 'value')
    value.text = '-1'  # Infinite capacity

    # Drop strategies
    drop_strategy = ET.SubElement(queue, 'parameter')
    drop_strategy.set('array', 'true')
    drop_strategy.set('classPath', 'java.lang.String')
    drop_strategy.set('name', 'dropStrategies')

    for r in range(K):
        ref_class = ET.SubElement(drop_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(drop_strategy, 'subParameter')
        sub_param.set('classPath', 'java.lang.String')
        sub_param.set('name', 'dropStrategy')
        value = ET.SubElement(sub_param, 'value')
        value.text = 'drop'

    # Queue get strategy (FCFS)
    strategy_param = ET.SubElement(queue, 'parameter')
    strategy_param.set('classPath', 'jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy')
    strategy_param.set('name', 'FCFSstrategy')

    # Queue put strategy
    put_strategy = ET.SubElement(queue, 'parameter')
    put_strategy.set('array', 'true')
    put_strategy.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategy')
    put_strategy.set('name', 'QueuePutStrategy')

    for r in range(K):
        ref_class = ET.SubElement(put_strategy, 'refClass')
        ref_class.text = classnames[r]

        sub_param = ET.SubElement(put_strategy, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy')
        sub_param.set('name', 'TailStrategy')

    # Impatience section (null for no impatience, required for JMT compatibility)
    impatience_param = ET.SubElement(queue, 'parameter')
    impatience_param.set('array', 'true')
    impatience_param.set('classPath', 'jmt.engine.NetStrategies.ImpatienceStrategies.Impatience')
    impatience_param.set('name', 'Impatience')

    for r in range(K):
        ref_class = ET.SubElement(impatience_param, 'refClass')
        ref_class.text = classnames[r]

        # Concrete strategy name, as at every other buffer: MATLAB's Buffer
        # branch runs the same saveImpatience.m for a Fork.
        sub_param = ET.SubElement(impatience_param, 'subParameter')
        sub_param.set('classPath', 'jmt.engine.NetStrategies.ImpatienceStrategies.Reneging')
        sub_param.set('name', 'Reneging')
        value = ET.SubElement(sub_param, 'value')
        value.text = 'null'

    # 2. ServiceTunnel section
    tunnel = ET.SubElement(node_elem, 'section')
    tunnel.set('className', 'ServiceTunnel')

    # 3. Fork section
    fork = ET.SubElement(node_elem, 'section')
    fork.set('className', 'Fork')

    # jobsPerLink parameter
    jpl_param = ET.SubElement(fork, 'parameter')
    jpl_param.set('classPath', 'java.lang.Integer')
    jpl_param.set('name', 'jobsPerLink')
    value = ET.SubElement(jpl_param, 'value')
    # sn.nodeparam carries fanOut as a float; the parameter is declared
    # java.lang.Integer, so JMT's Integer(String) ctor rejects "1.0" outright.
    value.text = str(int(round(float(fan_out))))

    # block parameter
    block_param = ET.SubElement(fork, 'parameter')
    block_param.set('classPath', 'java.lang.Integer')
    block_param.set('name', 'block')
    value = ET.SubElement(block_param, 'value')
    value.text = '-1'

    # isSimplifiedFork lets JMT ignore the branch list and send one job down
    # every link. That is only the same model when every branch is certain and
    # carries the same number of tasks, so a variable forking level switches it
    # off and makes JMT read the per-branch entries emitted below.
    is_simplified = True
    if fan_out_link is not None:
        taken = fan_out_prob > 0
        is_simplified = (bool(np.all(fan_out_prob[taken] == 1.0))
                         and bool(np.all(fan_out_link[taken] == fan_out))
                         and all(d is None for row in fan_out_dist for d in row))
    simpl_param = ET.SubElement(fork, 'parameter')
    simpl_param.set('classPath', 'java.lang.Boolean')
    simpl_param.set('name', 'isSimplifiedFork')
    value = ET.SubElement(simpl_param, 'value')
    value.text = 'true' if is_simplified else 'false'

    # ForkStrategy parameter
    strategy_param = ET.SubElement(fork, 'parameter')
    strategy_param.set('array', 'true')
    strategy_param.set('classPath', 'jmt.engine.NetStrategies.ForkStrategy')
    strategy_param.set('name', 'ForkStrategy')

    # Find outgoing connections for this fork
    outgoing_nodes = []
    if sn.connmatrix is not None:
        for j in range(sn.nnodes):
            if sn.connmatrix[node_idx, j] > 0:
                outgoing_nodes.append(j)

    nodenames = sn.nodenames if sn.nodenames else [f'Node{i+1}' for i in range(sn.nnodes)]

    # Check which classes actually route through this Fork
    # A class routes through this Fork if there's any incoming routing probability > 0
    rtnodes = sn.rtnodes if hasattr(sn, 'rtnodes') and sn.rtnodes is not None else None
    classes_using_fork = set()
    if rtnodes is not None:
        for r in range(K):
            # Check if any other node routes to this fork for class r
            col_idx = node_idx * K + r
            for i in range(sn.nnodes):
                if i != node_idx:
                    row_idx = i * K + r
                    if row_idx < rtnodes.shape[0] and col_idx < rtnodes.shape[1]:
                        if rtnodes[row_idx, col_idx] > 0:
                            classes_using_fork.add(r)
                            break
    else:
        # Fallback: assume all classes use all forks
        classes_using_fork = set(range(K))

    for r in range(K):
        ref_class = ET.SubElement(strategy_param, 'refClass')
        ref_class.text = classnames[r]

        class_strat = ET.SubElement(strategy_param, 'subParameter')
        class_strat.set('classPath', 'jmt.engine.NetStrategies.ForkStrategies.ProbabilitiesFork')
        class_strat.set('name', 'Branch Probabilities')

        emp_array = ET.SubElement(class_strat, 'subParameter')
        emp_array.set('array', 'true')
        emp_array.set('classPath', 'jmt.engine.NetStrategies.ForkStrategies.OutPath')
        emp_array.set('name', 'EmpiricalEntryArray')

        # One OutPathEntry per outgoing link. This used to emit only the last
        # output, matching a MATLAB/JAR defect that was invisible because
        # isSimplifiedFork makes JMT send one job down every link and ignore the
        # branch list; all four codebases now emit the full list.
        if r in classes_using_fork and outgoing_nodes:
            fork_outputs = list(outgoing_nodes)
        else:
            fork_outputs = []
        for out_node in fork_outputs:
            out_path = ET.SubElement(emp_array, 'subParameter')
            out_path.set('classPath', 'jmt.engine.NetStrategies.ForkStrategies.OutPath')
            out_path.set('name', 'OutPathEntry')

            # outUnitProbability
            emp_entry = ET.SubElement(out_path, 'subParameter')
            emp_entry.set('classPath', 'jmt.engine.random.EmpiricalEntry')
            emp_entry.set('name', 'outUnitProbability')

            # stationName
            station_name = ET.SubElement(emp_entry, 'subParameter')
            station_name.set('classPath', 'java.lang.String')
            station_name.set('name', 'stationName')
            value = ET.SubElement(station_name, 'value')
            value.text = nodenames[out_node]

            # probability: the branch activation probability
            branch_p = 1.0 if fan_out_prob is None else float(fan_out_prob[out_node, r])
            prob_param = ET.SubElement(emp_entry, 'subParameter')
            prob_param.set('classPath', 'java.lang.Double')
            prob_param.set('name', 'probability')
            value = ET.SubElement(prob_param, 'value')
            value.text = repr(branch_p)

            # JobsPerLinkDis is an EmpiricalEntry ARRAY: one entry per point of
            # the jobs-per-link distribution. A deterministic fork emits the
            # single degenerate entry it always did.
            dist = None if fan_out_dist is None else fan_out_dist[out_node][r]
            if dist is not None:
                jpl_points = list(dist.values)
                jpl_probs = list(dist.probs)
            elif fan_out_link is not None:
                jpl_points = [fan_out_link[out_node, r]]
                jpl_probs = [1.0]
            else:
                jpl_points = [fan_out]
                jpl_probs = [1.0]

            jpl_dis = ET.SubElement(out_path, 'subParameter')
            jpl_dis.set('classPath', 'jmt.engine.random.EmpiricalEntry')
            jpl_dis.set('array', 'true')
            jpl_dis.set('name', 'JobsPerLinkDis')

            for point, prob in zip(jpl_points, jpl_probs):
                jpl_entry = ET.SubElement(jpl_dis, 'subParameter')
                jpl_entry.set('classPath', 'jmt.engine.random.EmpiricalEntry')
                jpl_entry.set('name', 'EmpiricalEntry')

                # numbers (jobs per link)
                numbers_param = ET.SubElement(jpl_entry, 'subParameter')
                numbers_param.set('classPath', 'java.lang.String')
                numbers_param.set('name', 'numbers')
                value = ET.SubElement(numbers_param, 'value')
                value.text = str(int(round(float(point))))

                # probability for this distribution point
                prob_param2 = ET.SubElement(jpl_entry, 'subParameter')
                prob_param2.set('classPath', 'java.lang.Double')
                prob_param2.set('name', 'probability')
                value = ET.SubElement(prob_param2, 'value')
                value.text = repr(float(prob))


def _write_join_node(node_elem: ET.Element, node_idx: int, sn: NetworkStruct, classnames: List[str],
                     cs_node_names: Optional[Dict[Tuple[int, int], str]] = None):
    """Write Join node section.

    Join nodes in JMT have three sections:
    1. Join (with JoinStrategy)
    2. ServiceTunnel
    3. Router (dispatcher)
    """
    K = sn.nclasses

    # Get the number of incoming links (fanIn) - this is the number of tasks to wait for
    fan_in = 0
    if sn.connmatrix is not None:
        for i in range(sn.nnodes):
            if sn.connmatrix[i, node_idx] > 0:
                fan_in += 1

    # 1. Join section
    join = ET.SubElement(node_elem, 'section')
    join.set('className', 'Join')

    strategy_param = ET.SubElement(join, 'parameter')
    strategy_param.set('array', 'true')
    strategy_param.set('classPath', 'jmt.engine.NetStrategies.JoinStrategy')
    strategy_param.set('name', 'JoinStrategy')

    # Per-class join rule, as recorded by Network._refresh_node_param
    join_strategy = None
    join_required = None
    if sn.nodeparam is not None and node_idx in sn.nodeparam and isinstance(sn.nodeparam[node_idx], dict):
        join_strategy = sn.nodeparam[node_idx].get('joinStrategy')
        join_required = sn.nodeparam[node_idx].get('joinRequired')

    for r in range(K):
        ref_class = ET.SubElement(strategy_param, 'refClass')
        ref_class.text = classnames[r]

        # A quorum (PARTIAL with 1 <= k < siblings) is a PartialJoin, every
        # other rule a Standard Join. see saveJoinStrategy.m for the reference
        strategy_r = join_strategy[r] if join_strategy is not None and r < len(join_strategy) else JoinStrategy.STD
        required_r = int(join_required[r]) if join_required is not None and r < len(join_required) else -1
        is_quorum = strategy_r != JoinStrategy.STD and 0 < required_r < fan_in

        join_strat = ET.SubElement(strategy_param, 'subParameter')
        if is_quorum:
            join_strat.set('classPath', 'jmt.engine.NetStrategies.JoinStrategies.PartialJoin')
            join_strat.set('name', 'Quorum')
        else:
            join_strat.set('classPath', 'jmt.engine.NetStrategies.JoinStrategies.NormalJoin')
            join_strat.set('name', 'Standard Join')

        req_param = ET.SubElement(join_strat, 'subParameter')
        req_param.set('classPath', 'java.lang.Integer')
        req_param.set('name', 'numRequired')
        value = ET.SubElement(req_param, 'value')
        # -1 on a standard join is JMT's automatic count, taken from the fork
        value.text = str(required_r) if is_quorum else '-1'

    # 2. ServiceTunnel section
    tunnel = ET.SubElement(node_elem, 'section')
    tunnel.set('className', 'ServiceTunnel')

    # 3. Router section
    router = ET.SubElement(node_elem, 'section')
    router.set('className', 'Router')
    _write_routing_strategy(router, node_idx, sn, classnames, cs_node_names)


def _write_routing_strategy(router: ET.Element, node_idx: int, sn: NetworkStruct, classnames: List[str],
                            cs_node_names: Optional[Dict[Tuple[int, int], str]] = None):
    """Write routing strategy for a node.

    Handles different routing strategies (RAND, RROBIN, JSQ, PROB, etc.) by using
    the appropriate JMT strategy class.

    When class switching exists between this node and a destination, routes to
    the ClassSwitch node instead of the direct destination.

    Args:
        router: XML element for the Router section
        node_idx: Index of the current node
        sn: NetworkStruct object
        classnames: List of class names
        cs_node_names: Optional dict mapping (src_idx, dst_idx) to ClassSwitch node names
    """
    from ...sn.network_struct import RoutingStrategy

    if cs_node_names is None:
        cs_node_names = {}

    K = sn.nclasses
    M = sn.nnodes
    nodenames = sn.nodenames if sn.nodenames else [f'Node{i+1}' for i in range(M)]

    # Check if we have routing probability data
    has_rtnodes = hasattr(sn, 'rtnodes') and sn.rtnodes is not None and sn.rtnodes.size > 0
    has_connmatrix = hasattr(sn, 'connmatrix') and sn.connmatrix is not None
    has_routing = hasattr(sn, 'routing') and sn.routing is not None

    param = ET.SubElement(router, 'parameter')
    param.set('array', 'true')
    param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategy')
    param.set('name', 'RoutingStrategy')

    for r in range(K):
        ref_class = ET.SubElement(param, 'refClass')
        ref_class.text = classnames[r]

        # Get routing strategy for this node/class
        strategy = RoutingStrategy.RAND  # Default
        if has_routing and node_idx < sn.routing.shape[0] and r < sn.routing.shape[1]:
            strategy = RoutingStrategy(int(sn.routing[node_idx, r]))

        # Handle different routing strategies
        if strategy == RoutingStrategy.RROBIN:
            # Round Robin strategy
            sub_param = ET.SubElement(param, 'subParameter')
            sub_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.RoundRobinStrategy')
            sub_param.set('name', 'Round Robin')

        elif strategy == RoutingStrategy.WRROBIN:
            # Weighted Round Robin strategy
            sub_param = ET.SubElement(param, 'subParameter')
            sub_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.WeightedRoundRobinStrategy')
            sub_param.set('name', 'Weighted Round Robin')

            weight_array = ET.SubElement(sub_param, 'subParameter')
            weight_array.set('array', 'true')
            weight_array.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.WeightEntry')
            weight_array.set('name', 'WeightEntryArray')

            # Get routing weights for this node/class if available
            weights_key = (node_idx, r)
            has_weights = hasattr(sn, 'routingweights') and sn.routingweights and weights_key in sn.routingweights
            dest_weights = sn.routingweights.get(weights_key, {}) if has_weights else {}

            # Add weight entries for connected nodes
            if has_connmatrix:
                for j in range(M):
                    if sn.connmatrix[node_idx, j] > 0:
                        weight_entry = ET.SubElement(weight_array, 'subParameter')
                        weight_entry.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.WeightEntry')
                        weight_entry.set('name', 'WeightEntry')

                        station_param = ET.SubElement(weight_entry, 'subParameter')
                        station_param.set('classPath', 'java.lang.String')
                        station_param.set('name', 'stationName')
                        station_value = ET.SubElement(station_param, 'value')
                        station_value.text = nodenames[j]

                        weight_param = ET.SubElement(weight_entry, 'subParameter')
                        weight_param.set('classPath', 'java.lang.Integer')
                        weight_param.set('name', 'weight')
                        weight_value = ET.SubElement(weight_param, 'value')
                        # Use actual weight if available, otherwise default to 1
                        weight = int(dest_weights.get(j, 1))
                        weight_value.text = str(weight)

        elif strategy == RoutingStrategy.JSQ:
            # Join Shortest Queue strategy
            sub_param = ET.SubElement(param, 'subParameter')
            sub_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.ShortestQueueLengthRoutingStrategy')
            sub_param.set('name', 'Join the Shortest Queue (JSQ)')

        elif strategy == RoutingStrategy.SQ:
            # Power of K choices strategy
            sub_param = ET.SubElement(param, 'subParameter')
            sub_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.PowerOfKRoutingStrategy')
            sub_param.set('name', 'Power of k')

            k_param = ET.SubElement(sub_param, 'subParameter')
            k_param.set('classPath', 'java.lang.Integer')
            k_param.set('name', 'k')
            k_value = ET.SubElement(k_param, 'value')
            # sn.nodeparam[node][class]['d'] as in MATLAB saveRoutingStrategy; this was hardcoded 2, so SQ(d) was simulated as SQ(2) for every d.
            sq_d = 2
            try:
                sq_d = int(sn.nodeparam[node_idx][r]['d'])
            except (AttributeError, IndexError, TypeError, KeyError):
                pass
            k_value.text = str(sq_d)

            mem_param = ET.SubElement(sub_param, 'subParameter')
            mem_param.set('classPath', 'java.lang.Boolean')
            mem_param.set('name', 'withMemory')
            mem_value = ET.SubElement(mem_param, 'value')
            mem_value.text = 'false'

        elif strategy == RoutingStrategy.RAND:
            # RAND routing uses JMT's RandomStrategy (uniform random among connections)
            # This matches MATLAB's saveRoutingStrategy.m behavior (lines 31-34)
            sub_param = ET.SubElement(param, 'subParameter')
            sub_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.RandomStrategy')
            sub_param.set('name', 'Random')

        elif strategy == RoutingStrategy.PROB:
            # see _kb/06-solver-catalog.md (JMT: "RNG consumption, timing
            # authority, SPN wiring, servers") -- EmpiricalStrategy + class switching
            routing_probs = []
            if has_rtnodes and has_connmatrix:
                row_idx = node_idx * K + r
                if row_idx < sn.rtnodes.shape[0]:
                    for j in range(M):
                        if sn.connmatrix[node_idx, j] > 0:
                            # Check if there's a ClassSwitch node for this edge
                            if (node_idx, j) in cs_node_names:
                                # Sum probabilities across ALL destination classes
                                total_prob = 0.0
                                for s in range(K):
                                    col_idx = j * K + s
                                    if col_idx < sn.rtnodes.shape[1]:
                                        total_prob += sn.rtnodes[row_idx, col_idx]
                                if total_prob > 0:
                                    routing_probs.append((cs_node_names[(node_idx, j)], total_prob))
                            else:
                                # No class switching - use same-class probability
                                col_idx = j * K + r
                                if col_idx < sn.rtnodes.shape[1]:
                                    prob = sn.rtnodes[row_idx, col_idx]
                                    if prob > 0:
                                        routing_probs.append((nodenames[j], prob))

            if routing_probs:
                # Normalize probabilities to sum to 1.0 if needed
                total = sum(p for _, p in routing_probs)
                if total > 0 and abs(total - 1.0) > 1e-10:
                    routing_probs = [(name, prob / total) for name, prob in routing_probs]

                sub_param = ET.SubElement(param, 'subParameter')
                sub_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.EmpiricalStrategy')
                sub_param.set('name', 'Probabilities')

                emp_array = ET.SubElement(sub_param, 'subParameter')
                emp_array.set('array', 'true')
                emp_array.set('classPath', 'jmt.engine.random.EmpiricalEntry')
                emp_array.set('name', 'EmpiricalEntryArray')

                for dest_name, prob in routing_probs:
                    emp_entry = ET.SubElement(emp_array, 'subParameter')
                    emp_entry.set('classPath', 'jmt.engine.random.EmpiricalEntry')
                    emp_entry.set('name', 'EmpiricalEntry')

                    station_param = ET.SubElement(emp_entry, 'subParameter')
                    station_param.set('classPath', 'java.lang.String')
                    station_param.set('name', 'stationName')
                    station_value = ET.SubElement(station_param, 'value')
                    station_value.text = dest_name

                    prob_param = ET.SubElement(emp_entry, 'subParameter')
                    prob_param.set('classPath', 'java.lang.Double')
                    prob_param.set('name', 'probability')
                    prob_value = ET.SubElement(prob_param, 'value')
                    prob_value.text = f'{prob:.12f}'
            else:
                # Fallback to Random if no probabilities defined
                sub_param = ET.SubElement(param, 'subParameter')
                sub_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.RandomStrategy')
                sub_param.set('name', 'Random')

        elif strategy == RoutingStrategy.DISABLED:
            # Disabled routing
            sub_param = ET.SubElement(param, 'subParameter')
            sub_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.DisabledRoutingStrategy')
            sub_param.set('name', 'Disabled')

        else:
            # see _kb/06-solver-catalog.md (JMT: "RNG consumption, timing
            # authority, SPN wiring, servers")
            routing_probs = []
            if has_rtnodes and has_connmatrix:
                row_idx = node_idx * K + r
                if row_idx < sn.rtnodes.shape[0]:
                    # Iterate over all nodes (not just stations) to include Sink
                    n_nodes = sn.nnodes if hasattr(sn, 'nnodes') else sn.connmatrix.shape[1]
                    for j in range(n_nodes):
                        if j < sn.connmatrix.shape[1] and sn.connmatrix[node_idx, j] > 0:
                            # Check if there's a ClassSwitch node for this edge
                            if (node_idx, j) in cs_node_names:
                                # Sum probabilities across ALL destination classes
                                total_prob = 0.0
                                for s in range(K):
                                    col_idx = j * K + s
                                    if col_idx < sn.rtnodes.shape[1]:
                                        total_prob += sn.rtnodes[row_idx, col_idx]
                                if total_prob > 0:
                                    routing_probs.append((cs_node_names[(node_idx, j)], total_prob))
                            else:
                                # No class switching - use same-class probability
                                col_idx = j * K + r
                                if col_idx < sn.rtnodes.shape[1]:
                                    prob = sn.rtnodes[row_idx, col_idx]
                                    if prob > 0:
                                        routing_probs.append((nodenames[j], prob))

            if routing_probs:
                # Normalize to sum 1.0: rtnodes may hold RAND connection weights, not probabilities.
                total = sum(p for _, p in routing_probs)
                if total > 0 and abs(total - 1.0) > 1e-10:
                    routing_probs = [(name, prob / total) for name, prob in routing_probs]

                # Use EmpiricalStrategy with explicit probabilities
                sub_param = ET.SubElement(param, 'subParameter')
                sub_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.EmpiricalStrategy')
                sub_param.set('name', 'Probabilities')

                emp_array = ET.SubElement(sub_param, 'subParameter')
                emp_array.set('array', 'true')
                emp_array.set('classPath', 'jmt.engine.random.EmpiricalEntry')
                emp_array.set('name', 'EmpiricalEntryArray')

                for dest_name, prob in routing_probs:
                    emp_entry = ET.SubElement(emp_array, 'subParameter')
                    emp_entry.set('classPath', 'jmt.engine.random.EmpiricalEntry')
                    emp_entry.set('name', 'EmpiricalEntry')

                    station_param = ET.SubElement(emp_entry, 'subParameter')
                    station_param.set('classPath', 'java.lang.String')
                    station_param.set('name', 'stationName')
                    station_value = ET.SubElement(station_param, 'value')
                    station_value.text = dest_name

                    prob_param = ET.SubElement(emp_entry, 'subParameter')
                    prob_param.set('classPath', 'java.lang.Double')
                    prob_param.set('name', 'probability')
                    prob_value = ET.SubElement(prob_param, 'value')
                    prob_value.text = f'{prob:.12f}'
            else:
                # Fallback to Random strategy
                sub_param = ET.SubElement(param, 'subParameter')
                sub_param.set('classPath', 'jmt.engine.NetStrategies.RoutingStrategies.RandomStrategy')
                sub_param.set('name', 'Random')


def _parse_jsim_results(result_path: str, sn: NetworkStruct) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Parse JMT simulation results from XML output.

    Returns:
        Tuple of (Q, U, R, T, A) matrices
    """
    M = sn.nstations
    K = sn.nclasses

    Q = np.full((M, K), np.nan)
    U = np.full((M, K), np.nan)
    R = np.full((M, K), np.nan)
    T = np.full((M, K), np.nan)
    A = np.full((M, K), np.nan)

    if not os.path.exists(result_path):
        return Q, U, R, T, A

    try:
        tree = ET.parse(result_path)
        root = tree.getroot()

        classnames = sn.classnames if sn.classnames else [f'Class{i+1}' for i in range(K)]
        nodenames = sn.nodenames if sn.nodenames else [f'Node{i+1}' for i in range(sn.nnodes)]

        for measure in root.iter('measure'):
            measure_type = measure.get('measureType', measure.get('type', ''))
            node_name = measure.get('station', measure.get('referenceNode', ''))
            class_name = measure.get('class', measure.get('referenceUserClass', ''))
            mean_value = measure.get('meanValue', '0')
            # successful="false" still carries valid means (precision target not met).

            # Find station index
            station_idx = -1
            for i in range(M):
                node_idx = int(sn.stationToNode[i]) if sn.stationToNode is not None else i
                if node_idx < len(nodenames) and nodenames[node_idx] == node_name:
                    station_idx = i
                    break

            if station_idx < 0:
                continue

            # Find class index
            class_idx = -1
            for r in range(K):
                if r < len(classnames) and classnames[r] == class_name:
                    class_idx = r
                    break

            if class_idx < 0:
                continue

            try:
                value = float(mean_value)
            except ValueError:
                continue

            # For closed classes, filter metrics with insufficient analyzed samples (matches MATLAB getResults.m)
            # A class is considered recurrent only if analyzedSamples > total jobs in its chain
            njobs_flat = sn.njobs.flatten() if sn.njobs is not None else None
            if njobs_flat is not None and class_idx < len(njobs_flat) and not np.isinf(njobs_flat[class_idx]):
                # Find chain containing this class and sum all closed jobs in chain (MATLAB getResults.m line 144)
                chain_jobs = njobs_flat[class_idx]  # default: per-class
                if hasattr(sn, 'chains') and sn.chains is not None:
                    for c in range(sn.nchains):
                        if sn.chains[c, class_idx]:
                            chain_jobs = sum(
                                njobs_flat[r] for r in range(sn.nclasses)
                                if sn.chains[c, r] and not np.isinf(njobs_flat[r])
                            )
                            break
                analyzed_samples = int(measure.get('analyzedSamples', '0'))
                if analyzed_samples <= chain_jobs:
                    continue  # Leave as 0/NaN (starved class)

            if 'Number of Customers' in measure_type or 'QLen' in measure_type:
                Q[station_idx, class_idx] = value
            elif 'Utilization' in measure_type or 'Util' in measure_type:
                U[station_idx, class_idx] = value
            elif 'Response Time' in measure_type or 'RespT' in measure_type:
                R[station_idx, class_idx] = value
            elif 'Throughput' in measure_type or 'Tput' in measure_type:
                T[station_idx, class_idx] = value
            elif 'Arrival Rate' in measure_type or 'ArvR' in measure_type:
                A[station_idx, class_idx] = value

    except Exception as e:
        pass  # Return NaN-filled matrices on parse error

    # Set U to 0 for Fork and Join nodes (utilization is not meaningful for these)
    # This matches MATLAB's behavior where U is initialized to zeros
    for i in range(M):
        node_idx = int(sn.stationToNode[i]) if sn.stationToNode is not None else i
        node_type = sn.nodetype[node_idx] if sn.nodetype is not None and len(sn.nodetype) > node_idx else None
        if node_type == NodeType.FORK or node_type == NodeType.JOIN:
            for r in range(K):
                U[i, r] = 0.0

    return Q, U, R, T, A


def _is_jmva_method(method) -> bool:
    """True for the analytical JMVA methods, false for the simulation ones.

    The accepted spellings are MATLAB's (@SolverJMT/runAnalyzer.m): 'jmva'
    plus the ten algorithm suffixes, each also reachable with a 'jmt.' prefix.
    """
    m = str(method or '').lower()
    return m.startswith('jmva') or m.startswith('jmt.jmva')


def _parse_jmva_results(result_path: str, sn: NetworkStruct):
    """Parse a JMVA result file and disaggregate its chains back into classes.

    JMVA solves the CHAIN-AGGREGATED model that `write_jmva` hands it (one
    'ChainNN' customer class per chain, source stations omitted), so its
    station results are per chain and have to be pushed back onto the classes
    before they mean anything to LINE. This is a port of MATLAB
    `@SolverJMT/getResultsJMVA.m`, whose factors are specific to JMVA's own
    conventions -- in particular its 'Residence time' is per chain visit, not
    per class visit -- and therefore deliberately not the generic
    `sn_deaggregate_chain_results`.

    Stations are matched BY NAME rather than by position. MATLAB indexes the
    result list positionally against the station index, which holds only while
    the model has no Source: `writeJMVA.m` skips Source stations, so on an open
    model every station in the result is read against the wrong row of
    `rates`/`nservers`. Matching by name is right on both.

    Returns:
        (Q, U, R, T, lG) with the four (M x K) matrices and the log normalizing
        constant, NaN when the algorithm reported none.
    """
    M = sn.nstations
    K = sn.nclasses

    Q = np.zeros((M, K))
    U = np.zeros((M, K))
    R = np.zeros((M, K))
    T = np.zeros((M, K))
    lG = float('nan')

    if not os.path.exists(result_path):
        raise RuntimeError("JMT did not output a result file, the analysis has likely failed.")

    root = ET.parse(result_path).getroot()
    alg = root.find('./solutions/algorithm')
    if alg is None:
        raise RuntimeError("JMVA result file carries no solution: %s" % result_path)

    normconst = alg.find('normconst')
    if normconst is not None:
        try:
            lG = float(normconst.get('logValue'))
        except (TypeError, ValueError):
            lG = float('nan')

    demands = sn_get_demands_chain(sn)
    STchain = np.asarray(demands.STchain, dtype=float)
    Vchain = np.asarray(demands.Vchain, dtype=float)
    alpha = np.asarray(demands.alpha, dtype=float)

    rates = np.asarray(sn.rates, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        ST = np.where(rates != 0, 1.0 / rates, 0.0)
    ST = np.where(np.isnan(ST), 0.0, ST)

    nservers = np.asarray(sn.nservers, dtype=float).ravel()
    refstat = np.asarray(sn.refstat).ravel().astype(int)
    station_to_node = np.asarray(sn.stationToNode).ravel().astype(int)
    nodenames = sn.nodenames if sn.nodenames else []

    station_of_name = {}
    for i in range(M):
        node_idx = int(station_to_node[i])
        if node_idx < len(nodenames):
            station_of_name[nodenames[node_idx]] = i

    for statres in alg.findall('stationresults'):
        name = statres.get('station')
        if name not in station_of_name:
            raise RuntimeError("JMVA reported station '%s', which is not in the model" % name)
        i = station_of_name[name]
        for pos, classres in enumerate(statres.findall('classresults')):
            # 'ChainNN' is written by write_jmva in chain order; the name is
            # authoritative and the position is the fallback.
            c = pos
            cname = classres.get('customerclass') or ''
            if cname.lower().startswith('chain'):
                try:
                    c = int(cname[5:]) - 1
                except ValueError:
                    c = pos
            if c < 0 or c >= sn.nchains:
                raise RuntimeError("JMVA reported chain '%s', which is not in the model" % cname)
            inchain = np.asarray(sn.inchain[c]).ravel().astype(int)
            # A multiserver Queue is written as <ldstation servers="1">, and one
            # ldstation switches JMVA to its load-dependent algorithm, whose
            # Utilization is 1-p_i(0) at EVERY station, delay ones included. That
            # is a different random variable from LINE's E[busy servers], so no
            # rescaling recovers it; derive U from the chain throughput, which
            # both JMVA algorithms report alike.
            chain_tput = float('nan')
            for measure in classres.findall('measure'):
                if measure.get('measureType') == 'Throughput':
                    try:
                        chain_tput = float(measure.get('meanValue'))
                    except (TypeError, ValueError):
                        chain_tput = float('nan')
                    break
            for measure in classres.findall('measure'):
                mtype = measure.get('measureType')
                try:
                    value = float(measure.get('meanValue'))
                except (TypeError, ValueError):
                    value = float('nan')
                for k in inchain:
                    # A zero denominator means this class does not load this
                    # station in this chain, so it holds none of the chain's
                    # measure. MATLAB reaches the same rows as 0*Inf = NaN.
                    stc = STchain[i, c]
                    vref = Vchain[int(refstat[k]), c]
                    if mtype == 'Utilization':
                        if vref == 0:
                            continue
                        U[i, k] = ST[i, k] * chain_tput / vref * alpha[i, k]
                        # The divisor is the capacity write_jmva exported, which is
                        # max(nservers, max(lldscaling)): a load-dependent station
                        # carries its c in the scaling and leaves sn.nservers at 1,
                        # so reading nservers alone reported U = c * E[busy]/c.
                        c_eff = nservers[i]
                        lld = getattr(sn, 'lldscaling', None)
                        if lld is not None:
                            lld = np.atleast_2d(np.asarray(lld, dtype=float))
                            if i < lld.shape[0] and lld.shape[1] > 0:
                                c_eff = max(c_eff, float(np.max(lld[i, :])))
                        if np.isfinite(c_eff):
                            U[i, k] /= c_eff
                    elif mtype == 'Throughput':
                        T[i, k] = value * alpha[i, k]
                    elif mtype == 'Number of Customers':
                        if stc == 0 or vref == 0:
                            continue
                        Q[i, k] = value * ST[i, k] / stc / vref * alpha[i, k]
                    elif mtype == 'Residence time':
                        # JMVA reports residence over the chain's visits; LINE
                        # wants response time per visit of THIS class.
                        visits_c = np.asarray(sn.visits[c], dtype=float)
                        vik = visits_c[i, k] if i < visits_c.shape[0] and k < visits_c.shape[1] else 0.0
                        if stc == 0 or vref == 0 or vik == 0:
                            continue
                        R[i, k] = (value / vik) * ST[i, k] / stc / vref * alpha[i, k]

    # A Source carries no queue: JMVA never sees it (writeJMVA.m skips it), and
    # its throughput is the arrival rate the model already declares.
    nodetype = sn.nodetype
    for i in range(M):
        node_idx = int(station_to_node[i])
        if node_idx < len(nodetype) and nodetype[node_idx] == NodeType.SOURCE:
            for r in range(K):
                T[i, r] = 0.0 if np.isnan(rates[i, r]) else rates[i, r]

    return Q, U, R, T, lG


def _solve_jmva(sn: NetworkStruct, options: SolverJMTOptions, jmt_path, start_time) -> SolverJMTReturn:
    """Run one analytical JMVA solve and return its class-level metrics.

    Mirrors the 'jmva*' branch of MATLAB @SolverJMT/runAnalyzer.m: write the
    chain-aggregated .jmva, run the JMT command line in 'mva' mode, parse and
    disaggregate, then derive the arrival rates from the routing matrix, which
    is what MATLAB does here too -- an analytical solve measures nothing.
    """
    from ....solvers.wrappers.solver_qns.jmva_writer import write_jmva

    M = sn.nstations
    K = sn.nclasses
    method = str(getattr(options, 'method', 'jmva') or 'jmva').lower()

    workspace_dir = os.path.join(tempfile.gettempdir(), 'workspace', 'jmva')
    os.makedirs(workspace_dir, exist_ok=True)
    temp_dir = tempfile.mkdtemp(dir=workspace_dir)
    try:
        # 'model.jmva' matches MATLAB getJMVATempPath.m; the result path must
        # stay derived from it (Jmt.java writes args[1] + '-result.jmva').
        model_path = os.path.join(temp_dir, 'model.jmva')
        write_jmva(sn, model_path, {'method': method, 'samples': options.samples})
        print(f"JMT Model: {model_path}")

        returncode, _stdout, stderr = run_jmt(
            'mva', model_path, options.seed, jmt_jar=jmt_path, options=options,
            timeout=getattr(options, 'timeout', None), verbose=bool(options.verbose))
        if returncode != 0:
            raise RuntimeError(f"JMVA analysis failed: {stderr}")

        Q, U, R, T, lG = _parse_jmva_results(result_path_for(model_path, 'mva'), sn)
    finally:
        if not getattr(options, 'keep', False):
            shutil.rmtree(temp_dir, ignore_errors=True)

    W = R.copy()
    A = sn_get_arvr_from_tput(sn, T, T)
    if A is None or np.asarray(A).size == 0:
        A = np.zeros((M, K))

    X = np.zeros((1, K))
    refstat = np.asarray(sn.refstat).ravel().astype(int) if sn.refstat is not None else np.zeros(K, dtype=int)
    for r in range(K):
        if int(refstat[r]) < M:
            X[0, r] = T[int(refstat[r]), r]

    C = np.zeros((1, K))
    njobs = np.asarray(sn.njobs, dtype=float).ravel()
    for r in range(K):
        if not np.isinf(njobs[r]) and X[0, r] > 0:
            C[0, r] = njobs[r] / X[0, r]

    result = SolverJMTReturn(
        Q=Q, U=U, R=R, T=T, A=A, W=W, C=C, X=X,
        runtime=time.time() - start_time,
        method=method,
    )
    result.logNormConstAggr = lG
    return result


def _parse_jsim_fcr_results(result_path, sn):
    """Parse Finite Capacity Region (FCR) measures from JMT output.

    Returns (Qfcr, Ufcr, Rfcr, Wfcr, Tfcr, Afcr), each (nregions x nclasses).
    JMT reports one region-aggregate value per metric (not per class); mirroring
    MATLAB getResults.m, QLen/Tput are split evenly across classes while
    RespT/ResidT are replicated, and Util/ArvR are NaN (JMT does not provide
    them for regions). Regions are matched by the JMT-internal name FCRegion<n>.
    """
    F = int(sn.nregions) if getattr(sn, 'nregions', 0) else 0
    K = sn.nclasses
    Qfcr = np.zeros((F, K)); Rfcr = np.zeros((F, K))
    Wfcr = np.zeros((F, K)); Tfcr = np.zeros((F, K))
    Ufcr = np.full((F, K), np.nan); Afcr = np.full((F, K), np.nan)
    if F == 0 or not os.path.exists(result_path):
        return Qfcr, Ufcr, Rfcr, Wfcr, Tfcr, Afcr
    try:
        root = ET.parse(result_path).getroot()
        for measure in root.iter('measure'):
            node_type = measure.get('nodeType', '')
            ref = measure.get('station', measure.get('referenceNode', ''))
            if node_type != 'region' or not ref.startswith('FCRegion'):
                continue
            try:
                f = int(ref[len('FCRegion'):]) - 1
            except ValueError:
                continue
            if f < 0 or f >= F:
                continue
            try:
                val = float(measure.get('meanValue', '0'))
            except ValueError:
                continue
            mtype = measure.get('measureType', measure.get('type', ''))
            if 'Number of Customers' in mtype or 'QLen' in mtype:
                Qfcr[f, :] = val / K
            elif 'Residence Time' in mtype:
                Wfcr[f, :] = val
            elif 'Response Time' in mtype:
                Rfcr[f, :] = val
            elif 'Throughput' in mtype or 'Tput' in mtype:
                Tfcr[f, :] = val / K
    except Exception:
        pass
    return Qfcr, Ufcr, Rfcr, Wfcr, Tfcr, Afcr


def _compute_region_offered(sn, model, TN):
    """Offered arrival rate into each finite-capacity region by flow balance.

    For region f and class r, the offered rate is the rate at which stations
    OUTSIDE the region route class-r arrivals across the region boundary:
    sum over member stations i, external stations j, classes r' of
    TN[j, r'] * rt[(j,r'), (i,r)]. This is unaffected by the region drop (the
    upstream stations complete regardless), so offered - region throughput is
    the region drop rate. rt is in station-major (station*K + class) order.
    Returns an (nregions x nclasses) array, or None when the routing matrix is
    unavailable.
    """
    F = int(sn.nregions) if getattr(sn, 'nregions', 0) else 0
    if F == 0:
        return None
    rt = getattr(sn, 'rt', None)
    if rt is None:
        return None
    rt = np.asarray(rt, dtype=float)
    M = int(sn.nstations)
    K = int(sn.nclasses)
    if rt.shape != (M * K, M * K):
        return None
    TN = np.asarray(TN, dtype=float)
    members = getattr(sn, 'regionmembers', None)
    offered = np.zeros((F, K))
    for f in range(F):
        if members is not None and len(members) > f and members[f] is not None:
            mask = np.asarray(members[f]).ravel().astype(bool)
            member = set(int(i) for i in range(min(M, mask.size)) if mask[i])
        else:
            member = set()
        if not member:
            continue
        for r in range(K):
            acc = 0.0
            for i in member:
                for j in range(M):
                    if j in member:
                        continue
                    for rp in range(K):
                        w = rt[j * K + rp, i * K + r]
                        if w != 0.0:
                            acc += TN[j, rp] * w
            offered[f, r] = acc
    return offered


def parse_tran_resp_t(arv_file: str, dep_file: str, class_names: List[str] = None) -> Tuple[List[np.ndarray], np.ndarray, np.ndarray]:
    """
    Parse arrival and departure logs to calculate response times.

    Ported from MATLAB's JMTIO.parseTranRespT.

    Args:
        arv_file: Path to arrival log CSV file
        dep_file: Path to departure log CSV file
        class_names: Optional ordered list of class names from the model.
                     If provided, ensures deterministic class-to-index mapping.

    Returns:
        Tuple of (class_resp_t, job_resp_t, job_resp_t_arv_ts) where:
        - class_resp_t: List of response time arrays per class
        - job_resp_t: Response times per job
        - job_resp_t_arv_ts: Arrival timestamps for each response time
    """
    import csv

    # Parse arrival log: timestamp, job_id, class
    job_arv_ts = []
    job_arv_id = []
    job_arv_class = []

    try:
        with open(arv_file, 'r') as f:
            reader = csv.reader(f, delimiter=';')
            next(reader)  # Skip header
            for row in reader:
                if len(row) >= 4:
                    try:
                        job_arv_ts.append(float(row[1]))
                        job_arv_id.append(int(float(row[2])))
                        job_arv_class.append(row[3].strip())
                    except (ValueError, IndexError):
                        continue
    except FileNotFoundError:
        return [], np.array([]), np.array([])

    # Parse departure log: timestamp, job_id, class
    job_dep_ts = []
    job_dep_id = []
    job_dep_class = []

    try:
        with open(dep_file, 'r') as f:
            reader = csv.reader(f, delimiter=';')
            next(reader)  # Skip header
            for row in reader:
                if len(row) >= 4:
                    try:
                        job_dep_ts.append(float(row[1]))
                        job_dep_id.append(int(float(row[2])))
                        job_dep_class.append(row[3].strip())
                    except (ValueError, IndexError):
                        continue
    except FileNotFoundError:
        return [], np.array([]), np.array([])

    if not job_arv_id or not job_dep_id:
        return [], np.array([]), np.array([])

    # Build class name to index mapping
    # Use provided class_names for deterministic ordering; fall back to sorted names
    if class_names is not None:
        all_classes = list(class_names)
        # Add any classes found in logs but not in the model (shouldn't happen normally)
        for cls in set(job_arv_class + job_dep_class):
            if cls not in all_classes:
                all_classes.append(cls)
    else:
        all_classes = sorted(set(job_arv_class + job_dep_class))
    class_to_idx = {name: i for i, name in enumerate(all_classes)}
    num_classes = len(all_classes)

    # Mirrors MATLAB parseTranRespT: combine +1 arrivals/-1 departures per
    # job, sort by timestamp, discard pre-first-arrival data, pair alternately.
    # Group events by job ID
    job_events = {}
    for ts, jid, cls in zip(job_arv_ts, job_arv_id, job_arv_class):
        if jid not in job_events:
            job_events[jid] = []
        # Arrival: +1
        job_events[jid].append((ts, +1, class_to_idx.get(cls, 0)))

    for ts, jid, cls in zip(job_dep_ts, job_dep_id, job_dep_class):
        if jid not in job_events:
            job_events[jid] = []
        # Departure: -1
        job_events[jid].append((ts, -1, class_to_idx.get(cls, 0)))

    # Calculate response times per class
    class_resp_t = [[] for _ in range(num_classes)]
    all_resp_t = []
    all_arv_ts = []

    for jid, events in job_events.items():
        if len(events) < 2:
            continue

        # Sort by timestamp
        events = sorted(events, key=lambda x: x[0])

        # Find first arrival (event_type = +1)
        first_arv_idx = None
        for i, (ts, event_type, cls) in enumerate(events):
            if event_type > 0:
                first_arv_idx = i
                break

        if first_arv_idx is None:
            continue

        # Find last departure (event_type = -1)
        last_dep_idx = None
        for i in range(len(events) - 1, -1, -1):
            if events[i][1] < 0:
                last_dep_idx = i
                break

        if last_dep_idx is None or last_dep_idx <= first_arv_idx:
            continue

        # Keep only events from first arrival to last departure
        events = events[first_arv_idx:last_dep_idx + 1]

        # Pair arrivals with departures alternately
        # Events should now alternate: arv, dep, arv, dep, ...
        i = 0
        while i + 1 < len(events):
            arv_ts, arv_type, arv_cls = events[i]
            dep_ts, dep_type, dep_cls = events[i + 1]

            # Verify this is an arrival followed by departure
            if arv_type > 0 and dep_type < 0:
                resp_t = dep_ts - arv_ts
                if resp_t >= 0:
                    class_resp_t[arv_cls].append(resp_t)
                    all_resp_t.append(resp_t)
                    all_arv_ts.append(arv_ts)
                i += 2
            else:
                # Skip mismatched events (shouldn't happen normally)
                i += 1

    # Convert to numpy arrays
    class_resp_t = [np.array(rt) if rt else np.array([]) for rt in class_resp_t]

    return class_resp_t, np.array(all_resp_t), np.array(all_arv_ts)


def parse_tran_state(arv_file: str, dep_file: str, node_preload,
                     class_names: List[str] = None):
    """
    Reconstruct the per-class queue-length trajectory at a station from its
    arrival/departure logs. Faithful port of MATLAB ``JMTIO.parseTranState``.

    Args:
        arv_file: path to the ``<node>-Arv.csv`` log
        dep_file: path to the ``<node>-Dep.csv`` log
        node_preload: (K,) initial per-class population at the node
        class_names: ordered model class names (for deterministic class indexing)

    Returns:
        Tuple ``(state, evtype, evclass, evjob)`` where ``state`` is an
        ``(E+1, 1+K)`` array (column 0 = timestamps, columns 1..K = cumulative
        per-class queue length after each event, with a leading INIT row
        ``[0, preload]``); ``evtype`` are EventType values, ``evclass`` 0-based
        class indices (NaN for INIT) and ``evjob`` job ids (NaN for INIT).
    """
    import csv

    node_preload = np.asarray(node_preload, dtype=float).ravel()
    K = len(node_preload)

    def _read(path):
        ts, jid, cls = [], [], []
        try:
            with open(path, 'r') as f:
                reader = csv.reader(f, delimiter=';')
                next(reader, None)  # header
                for row in reader:
                    if len(row) >= 4:
                        try:
                            ts.append(float(row[1]))
                            jid.append(int(float(row[2])))
                            cls.append(row[3].strip())
                        except (ValueError, IndexError):
                            continue
        except FileNotFoundError:
            pass
        return ts, jid, cls

    arv_ts, arv_id, arv_cls = _read(arv_file)
    dep_ts, dep_id, dep_cls = _read(dep_file)

    if class_names is not None:
        name_to_idx = {name: i for i, name in enumerate(class_names)}
        K = max(K, len(class_names))
    else:
        allc = sorted(set(arv_cls + dep_cls))
        name_to_idx = {n: i for i, n in enumerate(allc)}
        K = max(K, len(allc))

    if len(node_preload) < K:
        node_preload = np.concatenate([node_preload, np.zeros(K - len(node_preload))])

    nA, nD = len(arv_ts), len(dep_ts)
    E = nA + nD
    state = np.zeros((E, 1 + K))
    evtype = np.empty(E, dtype=object)
    evclass = np.full(E, np.nan)
    evjob = np.full(E, np.nan)

    for i in range(nA):
        c = name_to_idx.get(arv_cls[i], 0)
        state[i, 0] = arv_ts[i]
        state[i, 1 + c] = +1.0
        evtype[i] = EventType.ARV
        evclass[i] = c
        evjob[i] = arv_id[i]
    for i in range(nD):
        c = name_to_idx.get(dep_cls[i], 0)
        state[nA + i, 0] = dep_ts[i]
        state[nA + i, 1 + c] = -1.0
        evtype[nA + i] = EventType.DEP
        evclass[nA + i] = c
        evjob[nA + i] = dep_id[i]

    # sort on timestamps (stable, mirrors MATLAB sortrows on column 1)
    order = np.argsort(state[:, 0], kind='stable')
    state = state[order]
    evtype = evtype[order]
    evclass = evclass[order]
    evjob = evjob[order]

    # prepend the INIT row [0, preload]
    state = np.vstack([np.concatenate([[0.0], node_preload[:K]]), state])
    evtype = np.concatenate([np.array([EventType.INIT], dtype=object), evtype])
    evclass = np.concatenate([[np.nan], evclass])
    evjob = np.concatenate([[np.nan], evjob])

    # Causality correction for zero-gap events, mirroring MATLAB parseTranState.
    ev_inst = np.where(np.diff(state[:, 0]) == 0)[0]
    for ev in ev_inst:
        ej = evjob[ev]
        if np.isnan(ej):
            continue
        prev_idx = [k for k in range(ev) if evjob[k] == ej]
        if not prev_idx:
            continue
        prev_ev = prev_idx[-1]
        nxt_idx = [k for k in range(ev + 1, len(evjob)) if evjob[k] == ej]
        if not nxt_idx:
            continue
        nxt = nxt_idx[0]
        et = evtype[ev]
        if (et == EventType.ARV and evtype[prev_ev] == EventType.ARV) or \
           (et == EventType.DEP and evtype[prev_ev] == EventType.DEP):
            state[[ev, nxt]] = state[[nxt, ev]]
            evjob[ev], evjob[nxt] = evjob[nxt], evjob[ev]
            evtype[ev], evtype[nxt] = evtype[nxt], evtype[ev]

    # cumulative per-class queue length
    for j in range(1, K + 1):
        state[:, j] = np.cumsum(state[:, j])

    return state, evtype, evclass, evjob


def parse_logs(model, is_node_logged: List[bool], metric_type: str = 'RespT') -> Dict:
    """
    Parse JMT log files to extract response time data.

    Ported from MATLAB's JMTIO.parseLogs.

    Args:
        model: Network model with log path
        is_node_logged: Boolean list indicating which nodes were logged
        metric_type: Type of metric to extract ('RespT' or 'QLen')

    Returns:
        Dict mapping (node_idx, class_idx) to response time data
    """
    import os

    log_path = model.get_log_path() if hasattr(model, 'get_log_path') else '/tmp'
    node_names = model.get_node_names() if hasattr(model, 'get_node_names') else []
    class_names = model.get_class_names() if hasattr(model, 'get_class_names') else []
    nnodes = len(node_names) if node_names else len(is_node_logged)
    nclasses = len(class_names) if class_names else 1

    log_data = {}

    for ind in range(nnodes):
        if not is_node_logged[ind]:
            continue

        node_name = node_names[ind] if ind < len(node_names) else f'Node{ind}'
        arv_file = os.path.join(log_path, f"{node_name}-Arv.csv")
        dep_file = os.path.join(log_path, f"{node_name}-Dep.csv")

        if not os.path.exists(arv_file) or not os.path.exists(dep_file):
            continue

        if metric_type == 'RespT':
            class_resp_t, _, _ = parse_tran_resp_t(arv_file, dep_file)

            for r, resp_times in enumerate(class_resp_t):
                if len(resp_times) > 0:
                    log_data[(ind, r)] = {
                        'RespT': resp_times,
                        't': np.arange(len(resp_times))  # Placeholder for timestamps
                    }

    return log_data


def solver_jmt(
    sn: NetworkStruct,
    options: Optional[SolverJMTOptions] = None,
    model: Any = None
) -> SolverJMTReturn:
    """
    JMT solver handler - calls JMT via subprocess.

    Performs discrete-event simulation using JMT by:
    1. Writing the model to JSIM XML format
    2. Calling JMT via subprocess
    3. Parsing the results

    Args:
        sn: Network structure
        options: Solver options
        model: Optional Network model (for FCR regions)

    Returns:
        SolverJMTReturn with all performance metrics

    Raises:
        RuntimeError: If JMT is not available or fails
    """
    start_time = time.time()

    if options is None:
        options = SolverJMTOptions()

    # A local JVM is the default backend but no longer the only one: the runner
    # also speaks to a JMT REST server (options.rest_url) and, when no JVM
    # exists, to the JMT Docker image with the user's consent. Only the local
    # backend needs the jar, so its absence is not fatal by itself.
    # The jar is resolved (and, when absent, downloaded) only for the local
    # backend: a host with no JVM must reach the Docker fallback without first
    # fetching 50MB it will never load.
    rest_url = getattr(options, 'rest_url', None)
    jmt_path = None
    if not rest_url and has_java():
        jmt_path = _get_jmt_jar_path()

    # THE METHOD SELECTS THE ENGINE, and until 2026-08-09 it selected nothing:
    # this handler wrote a .jsim and simulated whatever was asked, so a caller
    # asking for the analytical JMVA silently received a simulation (and the
    # ten jmva* aliases that listValidMethods advertises had no test covering
    # any of them). MATLAB @SolverJMT/runAnalyzer.m branches here.
    if _is_jmva_method(getattr(options, 'method', None)):
        return _solve_jmva(sn, options, jmt_path, start_time)

    M = sn.nstations
    K = sn.nclasses

    # Create temporary directory in /tmp/workspace/jsim/ to match JAR behavior
    workspace_dir = os.path.join(tempfile.gettempdir(), 'workspace', 'jsim')
    os.makedirs(workspace_dir, exist_ok=True)
    temp_dir = tempfile.mkdtemp(dir=workspace_dir)

    try:
        # 'model.jsim' matches MATLAB getJSIMTempPath.m and the JAR; the
        # result path must stay derived from model_path (Jmt.java:128 args[1]+'-result.jsim').
        model_path = os.path.join(temp_dir, 'model.jsim')
        result_path = result_path_for(model_path, 'sim')  # JMT creates <model>-result.jsim

        # Write model to JSIM format
        _write_jsim_file(sn, model_path, options, model)

        # Print model path (matching wrapper behavior)
        print(f"JMT Model: {model_path}")

        # Wall-clock budget: only a FINITE options.timeout is one. An infinite
        # one is the default and means the caller asked for no budget, as in
        # MATLAB (the local JVM arm is untimed) and the JAR
        # (simulationTimeoutSeconds = 0). See run_jmt for why the old 600 s
        # fallback was a wrong answer rather than a safety net.
        _tmo = float(getattr(options, 'timeout', float('inf')))
        _sub_timeout = _tmo if math.isfinite(_tmo) and _tmo > 0 else None

        # Execute through the backend the options select: local JVM by default,
        # a JMT REST server when rest_url is set, the JMT Docker image when no
        # JVM exists and the user consents.
        try:
            returncode, _stdout, _stderr = run_jmt(
                'sim', model_path, options.seed,
                jmt_jar=jmt_path, options=options,
                timeout=_sub_timeout, verbose=bool(options.verbose))
        except subprocess.TimeoutExpired:
            # A KILLED SIMULATION IS NOT A RESULT. Returning an empty table left
            # the analysis reading as COMPLETED with no rows, which every
            # consumer downstream reports as the solver never having run -- the
            # parity harness says "solver JMT missing from the recorded
            # results", indistinguishable from a solver that refused the model.
            # The budget is the caller's, so its expiry is theirs to hear about.
            raise RuntimeError(
                "JMT simulation exceeded the %gs wall-clock budget (options.timeout) "
                "and was terminated after %.1fs, so it produced no result."
                % (_tmo, time.time() - start_time)) from None
        except JMTBackendError as exc:
            raise RuntimeError(str(exc)) from None

        if returncode != 0:
            stderr = _stderr
            raise RuntimeError(f"JMT simulation failed: {stderr}")

        # Parse results (includes arrival rates from JMT)
        Q, U, R, T, A = _parse_jsim_results(result_path, sn)
        Qfcr, Ufcr, Rfcr, Wfcr, Tfcr, Afcr = _parse_jsim_fcr_results(result_path, sn)
        # Region loss. JMT exposes no region drop measure, and its "FCR Arrival
        # Rate" reports the ADMITTED rate (equal to region throughput), so the
        # offered rate must be reconstructed by flow balance: the rate at which
        # external stations route jobs across the region boundary, unaffected by
        # the drop. DropRateNfcr = offered - carried (region Tput), clamped at 0.
        TNfcr = None
        DropRateNfcr = None
        if Tfcr is not None:
            offered = _compute_region_offered(sn, model, T)
            TNfcr = Tfcr.copy()
            if offered is not None:
                DropRateNfcr = np.maximum(0.0, offered - Tfcr)
                # Only DROP-rule (region drop id 1) classes lose jobs; WAITQ and
                # blocking regions back up instead, so zero their loss.
                regionrule = getattr(sn, 'regionrule', None)
                if regionrule is not None:
                    rr = np.atleast_2d(np.asarray(regionrule, dtype=float))
                    if rr.shape == DropRateNfcr.shape:
                        DropRateNfcr[rr != 1] = 0.0

        # Zero mask: when RespT is near zero, zero out QLen and Util
        # (matches JAR NetworkSolver.getAvg() and MATLAB @NetworkSolver/getAvg.m)
        fine_tol = 1e-8  # GlobalConstants.FineTol
        zero_mask = np.where(np.isnan(R), True, R < 10 * fine_tol)
        Q[zero_mask] = 0.0
        U[zero_mask] = 0.0
        R[zero_mask] = 0.0

        # Set source station metrics only if JMT didn't report them
        # For source: Q=0, U=0, R=0, A=0 (T is reported by JMT as simulated throughput)
        rates = np.asarray(sn.rates) if sn.rates is not None else None
        for i in range(M):
            node_idx = int(sn.stationToNode[i]) if sn.stationToNode is not None else i
            if node_idx < len(sn.nodetype) and sn.nodetype[node_idx] == NodeType.SOURCE:
                for r in range(K):
                    # Only set if JMT didn't report (NaN)
                    if np.isnan(Q[i, r]):
                        Q[i, r] = 0.0
                    if np.isnan(U[i, r]):
                        U[i, r] = 0.0
                    if np.isnan(R[i, r]):
                        R[i, r] = 0.0
                    if np.isnan(A[i, r]):
                        A[i, r] = 0.0  # No arrivals to source itself
                    # Only use theoretical rate if JMT didn't report throughput
                    if np.isnan(T[i, r]):
                        if rates is not None and node_idx < rates.shape[0] and r < rates.shape[1]:
                            T[i, r] = rates[node_idx, r]
                        else:
                            T[i, r] = 0.0

        W = R.copy()

        # System throughput (sum at reference stations)
        X = np.zeros((1, K))
        refstat = sn.refstat.flatten() if sn.refstat is not None else np.zeros(K, dtype=int)
        for r in range(K):
            if int(refstat[r]) < M:
                X[0, r] = T[int(refstat[r]), r]

        # Cycle times
        C = np.zeros((1, K))
        njobs = sn.njobs.flatten() if sn.njobs is not None else np.zeros(K)
        for r in range(K):
            if not np.isinf(njobs[r]) and X[0, r] > 0:
                C[0, r] = njobs[r] / X[0, r]

        runtime = time.time() - start_time

        return SolverJMTReturn(
            Q=Q,
            U=U,
            R=R,
            T=T,
            A=A,
            W=W,
            C=C,
            X=X,
            Qfcr=Qfcr,
            Ufcr=Ufcr,
            Rfcr=Rfcr,
            Wfcr=Wfcr,
            Tfcr=Tfcr,
            Afcr=Afcr,
            TNfcr=TNfcr,
            DropRateNfcr=DropRateNfcr,
            runtime=runtime,
            method='jsim'
        )

    finally:
        # Clean up temporary directory (unless keep=True)
        if not getattr(options, 'keep', False):
            shutil.rmtree(temp_dir, ignore_errors=True)


def _jmt_station_cap_assert(sn, ist):
    """Refuse a binding station capacity that a CLOSED class can reach.

    Refused on two counts.

    (1) THE RULE IS ONE JMT CANNOT READ -- BBS, RSRD or retrial-with-limit; see
    ``_JMT_READABLE_DROP_IDS``. Such a value is not approximated, it is IGNORED,
    so the capacity stops being enforced and JMT returns the unconstrained
    answer.

    (2) THE RULE IS WAITQ AND A CLOSED CLASS CAN REACH THE LIMIT, the same reason
    ``jmtClassCapAssert`` (MATLAB ``saveRegions.m``, C++
    ``assert_class_cap_exportable``) refuses the per-class one: JMT cannot hold
    a blocked closed job at its upstream station. Note this is the case where NO
    blocking rule is declared. A model that does declare BAS is exported as JMT
    "BAS blocking", which is the same queueing model, under either declaration
    form -- see ``_jmt_is_bas_destination``. That the limit CAN be reached is the
    caller's to establish and is not retested here: the buffer-capacity writer
    reaches this function only for a capacity strictly below the total
    population, which is the one thing that makes a buffer a buffer.

    That entry advised expressing the limit as the STATION capacity instead;
    measured on 2026-08-19 the advice was wrong, and neither of the two
    strategies a WAITQ station maps onto reproduces the UNDECLARED case:

    ``waiting queue``
        does not enforce ``size`` at all. On a closed 3-queue tandem, N=6,
        Exp(1) FCFS, cap 2 at Q2, JMT returned the UNCONSTRAINED
        [2.03 1.99 1.98], X = 0.750, against the exact [3.6090 0.9711 1.4199],
        X = 0.6522.
    ``BAS blocking``
        enforces it, but completes the service BEFORE blocking, so the blocked
        job moves the instant room frees -- a different queueing model, not a
        rounding: same fixture, [2.871 1.373 1.756], X = 0.7126.

    With no rule declared LINE instead disables the upstream departure while the
    destination is full, which for exponential service is repetitive service (RS)
    and is what SolverCTMC, SolverSSA and SolverLDES all agree on. So THAT model
    is refused rather than exported as either of the two things JMT can say. See
    BUG-81. A declared-BAS model is a different model and is exported, not
    refused: blocking after service is precisely what JMT's "BAS blocking" does.
    """
    cap = np.asarray(sn.cap).flatten()
    if ist >= cap.size or np.isinf(cap[ist]):
        return
    njobs = np.asarray(sn.njobs).flatten()
    rates = np.atleast_2d(np.asarray(sn.rates, dtype=float))
    droprule = np.atleast_2d(np.asarray(sn.droprule))
    for r in range(int(sn.nclasses)):
        if ist < rates.shape[0] and r < rates.shape[1] and np.isnan(rates[ist, r]):
            continue   # class r is not served at this station
        # -1 is DropStrategy.WAITQ; the ids are shared across the codebases and
        # read here as numbers, as _drop_strategy_text does.
        dr = int(droprule[ist, r]) if ist < droprule.shape[0] and r < droprule.shape[1] else -1
        if dr != 0 and dr not in _JMT_READABLE_DROP_IDS:
            # Unmappable for EITHER class type, so this test precedes the open-class skip
            raise ValueError(
                "SolverJMT: station '%s' applies drop strategy \"%s\" to class '%s' and carries "
                "a finite capacity %d it can reach. JMT's queue section reads only \"drop\", "
                "\"waiting queue\", \"BAS blocking\" and \"retrial\"; it does not approximate "
                "anything else, it ignores it, so the capacity would stop being enforced and the "
                "run would return the unconstrained answer. Use SolverCTMC, SolverSSA or "
                "SolverLDES."
                % (sn.nodenames[int(sn.stationToNode[ist])],
                   _DROP_STRATEGY_TEXT.get(dr, str(dr)), sn.classnames[r], int(cap[ist])))
        if r >= njobs.size or np.isinf(njobs[r]):
            continue   # open class: JMT loses its arrivals, as LINE does
        if dr != -1:
            continue   # a mappable declared blocking rule is exported as itself
        if _jmt_is_bas_destination(sn, ist, r):
            # BAS declared on the UPSTREAM station: _drop_strategy_text moves it
            # onto this one, which is where JMT reads it
            continue
        raise ValueError(
            "SolverJMT: station '%s' carries a finite capacity %d that binds for the closed "
            "class '%s'. LINE blocks a closed job that finds no room -- the upstream departure "
            "is disabled and the job stays where it is -- and no JMT drop strategy reproduces "
            "that: \"waiting queue\" does not enforce the size at all, and \"BAS blocking\" "
            "completes the service before blocking, which is a different queueing model. "
            "Use SolverCTMC, SolverSSA or SolverLDES, or declare DropStrategy.BAS if "
            "blocking after service is the model you want, which SolverJMT does export."
            % (sn.nodenames[int(sn.stationToNode[ist])], int(cap[ist]), sn.classnames[r]))


def _jmt_reachable_population(sn, ist):
    """Most jobs that can be present at station ``ist``.

    Read the way ``refresh_capacity`` derives the capacity itself: per CHAIN,
    because a chain's whole population can reach a station that serves any one
    of its classes (class switching moves jobs between them), and a chain none
    of whose classes is served there cannot put a single job on it.

    Deliberately NOT read off ``sn.classcap``, which ``refresh_capacity`` has
    already clamped by the station's own cap: comparing a capacity against a
    quantity derived from it would make every user-declared buffer look
    non-binding. ``inf`` when an open chain is served here, which is what the
    model total gave before and which sends the station to
    ``_jmt_station_cap_assert``, where the open classes are skipped by name.
    """
    if getattr(sn, 'njobs', None) is None:
        return np.inf
    njobs = np.asarray(sn.njobs).ravel()
    inchain = getattr(sn, 'inchain', None)
    rates = getattr(sn, 'rates', None)
    if inchain is None or rates is None:
        return np.sum(njobs)
    total = 0.0
    for c in range(int(getattr(sn, 'nchains', 0))):
        members = np.asarray(inchain[c]).ravel().astype(int)
        if members.size == 0:
            continue
        served = any(not np.isnan(rates[ist, r]) for r in members if r < rates.shape[1])
        if served:
            total += float(np.sum(njobs[members]))
    return total
