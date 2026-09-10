"""
Native Python implementation of LQNS solver.

This module provides a native Python wrapper for the lqns and lqsim command-line
tools that analyze layered queueing networks.
"""

import numpy as np
import subprocess
import tempfile
import shutil
import os
import platform
import xml.etree.ElementTree as ET
import json
import urllib.request
import urllib.error
from dataclasses import dataclass, field
from ....constants import default_verbose
from typing import Optional, Dict, Any, List, Tuple
import pandas as pd

from ....api.io.logging import line_debug, line_ack
from ...base import Solver
from ....layered import LayeredNetworkElement

# The kinds lqn.type stores, under the names parseXMLResults.m uses for them.
LQN_HOST = int(LayeredNetworkElement.PROCESSOR)
LQN_TASK = int(LayeredNetworkElement.TASK)
LQN_ENTRY = int(LayeredNetworkElement.ENTRY)
LQN_ACTIVITY = int(LayeredNetworkElement.ACTIVITY)
LQN_CALL = int(LayeredNetworkElement.CALL)


@dataclass
class LQNSOptions:
    """Options for the LQNS solver."""
    method: str = 'default'  # default, lqns, srvn, exactmva, srvn.exactmva, sim, lqsim
    multiserver: str = 'rolia'  # rolia, conway, etc.
    samples: int = 10000  # For simulation methods
    verbose: bool = field(default_factory=default_verbose)
    keep: bool = True  # Keep temporary files
    seed: int = 23000
    remote: bool = False  # Enable remote execution via REST API
    remote_url: str = 'http://localhost:8080'  # URL of lqns-rest server
    # see _kb/06-solver-catalog.md (Wrappers) for the lang toggle semantics
    lang: str = field(default_factory=lambda: os.environ.get('LINE_SOLVER_LANG', 'python'))


@dataclass
class LQNSResult:
    """Result from LQNS solver."""
    PN: np.ndarray = field(default_factory=lambda: np.array([]))  # Processor utilization
    SN: np.ndarray = field(default_factory=lambda: np.array([]))  # Phase 1 service times
    TN: np.ndarray = field(default_factory=lambda: np.array([]))  # Throughputs
    UN: np.ndarray = field(default_factory=lambda: np.array([]))  # Utilizations
    RN: np.ndarray = field(default_factory=lambda: np.array([]))  # Processor waiting
    QN: np.ndarray = field(default_factory=lambda: np.array([]))  # Queue lengths
    AN: np.ndarray = field(default_factory=lambda: np.array([]))  # Arrival rates
    WN: np.ndarray = field(default_factory=lambda: np.array([]))  # Residence times
    runtime: float = 0.0
    method: str = 'default'
    iterations: int = 0


def _snap_to_tenth(values):
    """Snap values that are within CoarseTol (relative) of an exact tenth.

    lqns reports quantities that are analytically an exact multiple of 0.1
    (a deterministic 0.4 service time, say) with a few digits of iteration
    noise, as 0.39990. The MATLAB LQNS wrapper removes that noise in
    getAvgTable before tabulating, and the Python row must apply the same
    step or the two rows disagree on values neither solver actually computed
    differently.

    References:
        MATLAB: matlab/src/solvers/wrappers/LQNS/@SolverLQNS/SolverLQNS.m
    """
    from ....constants import GlobalConstants
    arr = np.asarray(values, dtype=float).copy()
    if arr.size == 0:
        return arr
    scaled = arr * 10.0
    snapped = np.round(scaled)
    with np.errstate(invalid='ignore'):
        # The tolerance is relative to the scaled value, exactly as in MATLAB;
        # NaN and non-positive entries fail the comparison and pass through.
        mask = np.abs(scaled - snapped) < GlobalConstants.CoarseTol * scaled
    arr[mask] = snapped[mask] / 10.0
    return arr


class SolverLQNS(Solver):
    """
    Native Python LQNS solver using external lqns/lqsim tools.

    This solver wraps the LQNS (Layered Queueing Network Solver) command-line tools
    to analyze layered queueing networks. The CLI is called directly without going
    through SolverLQNS.

    Supported methods:
    - default/lqns: Standard LQNS analytical solver
    - srvn: SRVN layering
    - exactmva: Exact MVA algorithm
    - srvn.exactmva: SRVN with exact MVA
    - sim/lqsim: Simulation-based solver

    Requirements:
        The 'lqns' and 'lqsim' commands must be available in the system PATH.
        Install from: http://www.sce.carleton.ca/rads/lqns/

        LINE ships no LQNS binary and runs none from a container image: LQNS is
        distributed under an evaluation agreement that forbids redistribution.

        Alternatively, point LINE at a host that already runs LQNS:
          options.config.remote = True; options.config.remote_url = 'http://localhost:8080'

    Note:
        Model serialization to LQNX format
        method. The native aspect is in the CLI execution and result parsing.

    Example:
        >>> lqn = LayeredNetwork('model')
        >>> # ... build model ...
        >>> solver = SolverLQNS(lqn, LQNSOptions(method='lqns'))
        >>> result = solver.runAnalyzer()
    """

    def __init__(self, model, options: Optional[LQNSOptions] = None, **kwargs):
        """
        Initialize the LQNS solver.

        Args:
            model: LayeredNetwork model
            options: Optional LQNSOptions configuration
            **kwargs: Additional parameters (keep, etc.) for compatibility
        """
        self.model = model
        if options is not None:
            self.options = options
        else:
            # Route recognized kwargs into options; unknown ones are ignored.
            import dataclasses
            valid = {f.name for f in dataclasses.fields(LQNSOptions)}
            self.options = LQNSOptions(**{k: v for k, v in kwargs.items() if k in valid})
        self._result: Optional[LQNSResult] = None
        self._lqn = None  # Will be set during analysis
        self._keep = kwargs.get('keep', True)  # For compatibility

    @staticmethod
    def isAvailable() -> bool:
        """
        Check if lqns can be run: a native binary is on the PATH.

        LINE never runs LQNS from a container image, because its licence forbids
        redistribution. To exercise a containerised build, put a shim on the PATH
        with run-tests.sh --lqns-docker.

        Returns:
            True if a native lqns binary is available, False otherwise
        """
        return SolverLQNS._has_native_lqns()

    @staticmethod
    def _has_native_lqns() -> bool:
        """Check if a native lqns binary is available in the system PATH."""
        try:
            if platform.system() == 'Windows':
                process = subprocess.Popen(
                    ['lqns', '--help'],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE
                )
            else:
                process = subprocess.Popen(
                    ['lqns', '--help'],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE
                )

            try:
                stdout, stderr = process.communicate(timeout=5)
                output = (stdout + stderr).decode().lower()
                return 'lqns' in output or 'usage' in output
            except subprocess.TimeoutExpired:
                process.kill()
                return False

        except (FileNotFoundError, OSError):
            return False

    @staticmethod
    def listValidMethods() -> List[str]:
        """List valid methods for the LQNS solver."""
        return ['default', 'lqns', 'srvn', 'exactmva', 'srvn.exactmva', 'sim', 'lqsim', 'lqnsdefault']

    def isStochasticMethod(self, method):
        """The lqsim simulator is stochastic; the analytical lqns/srvn methods
        are deterministic.
        """
        return str(method).lower() in ('sim', 'lqsim') if method else False

    is_stochastic_method = isStochasticMethod

    def getStruct(self):
        """Get the LayeredNetworkStruct."""
        if hasattr(self.model, 'getStruct'):
            return self.model.getStruct()
        return None

    def runAnalyzer(self) -> LQNSResult:
        """
        Run the LQNS analysis.

        Returns:
            LQNSResult containing performance metrics

        Raises:
            RuntimeError: If lqns is not available or fails
        """
        line_ack('LQNS', self.options.verbose)
        multiserver = self.options.multiserver or 'default'

        # Solver console: SolverLQNS extends Solver, not NetworkSolver, so it
        # never reaches the shared run hook in solvers/base.py and opens its
        # own run here, as the MATLAB wrapper does.
        from line_solver.api.io import console as _console
        with _console.run_scope(self, self.options):
            return self._runAnalyzerBody(multiserver)

    def _runAnalyzerBody(self, multiserver: str) -> LQNSResult:
        """Write the model, run the binary, parse what came back."""
        from line_solver.api.io import console as _console
        line_debug("LQNS: starting (method=%s, multiserver=%s)", self.options.method, multiserver, options=self.options)

        import time
        start_time = time.time()

        # LINE_WORKSPACE_ROOT relocates the staging dir; run-tests.sh sets it when
        # a solver is wrapped in a container, so the path is bind-mountable.
        workspace_root = os.environ.get('LINE_WORKSPACE_ROOT', '').strip()
        if workspace_root:
            base = os.path.join(workspace_root, 'line_workspace', 'lqns')
            os.makedirs(base, exist_ok=True)
            temp_dir = tempfile.mkdtemp(prefix='lqns_', dir=base)
        else:
            temp_dir = tempfile.mkdtemp(prefix='lqns_')

        try:
            # Write model to LQNX format
            lqnx_file = os.path.join(temp_dir, 'model.lqnx')

            # Check if model is native LayeredNetwork
            model_type = type(self.model).__name__
            _console.step('writing the LQN model to %s', lqnx_file)
            if model_type == 'LayeredNetwork' or hasattr(self.model, 'processors'):
                # Use native writeXML
                self.model.writeXML(lqnx_file, False)
            elif hasattr(self.model, 'writeXML'):
                # Use LayeredNetwork writeXML
                self.model.writeXML(lqnx_file, False)
            else:
                raise RuntimeError(f"Model type {model_type} does not support writeXML")

            # Check for remote execution
            if self.options.remote:
                line_debug("LQNS: using remote execution at %s", self.options.remote_url, options=self.options)
                if self.options.verbose and not _console.is_active():
                    print(f"Using remote LQNS at: {self.options.remote_url}")
                _console.step('running the lqns service at %s', self.options.remote_url)
                self._run_remote_lqns(lqnx_file)
            else:
                # Build and execute command locally
                cmd = self._build_command(lqnx_file)
                line_debug("LQNS: using local execution, command: %s", cmd, options=self.options)

                # the console already reports the command through the routed
                # line_debug above, so printing it again would duplicate it
                if self.options.verbose and not _console.is_active():
                    print(f"LQNS command: {cmd}")

                _console.step('running the lqns binary as a subprocess')

                # Execute
                if platform.system() == 'Windows':
                    process = subprocess.Popen(
                        cmd.split(),
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT
                    )
                else:
                    process = subprocess.Popen(
                        cmd.split(),
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT
                    )

                stdout, _ = process.communicate()
                output = stdout.decode()

                if process.returncode != 0:
                    raise RuntimeError(
                        f"LQNS/LQSIM did not terminate correctly.\n"
                        f"Exit code: {process.returncode}\n"
                        f"Output: {output}"
                    )

            # Parse results from .lqxo file
            lqxo_file = lqnx_file.replace('.lqnx', '.lqxo')
            _console.step('parsing the lqns XML results')
            self._parse_xml_results(lqxo_file)

        finally:
            # Clean up
            if not self.options.keep:
                shutil.rmtree(temp_dir, ignore_errors=True)

        runtime = time.time() - start_time
        if self._result:
            self._result.runtime = runtime

        return self._result


    def _build_command(self, filename: str) -> str:
        """Build the LQNS/LQSIM command line."""
        method = self.options.method.lower()

        # Verbose flags
        verbose_flag = '' if self.options.verbose else '-a -w'

        # Multiserver policy
        praqma_flag = ''
        if method not in ('lqsim', 'sim'):
            pol = self.options.multiserver or 'rolia'
            praqma_flag = f'-Pmultiserver={pol}'

        common = f'{verbose_flag} {praqma_flag} -Pstop-on-message-loss=false -x {filename}'

        if method == 'srvn':
            cmd = f'lqns {common} -Playering=srvn'
        elif method == 'exactmva':
            cmd = f'lqns {common} -Pmva=exact'
        elif method == 'srvn.exactmva':
            cmd = f'lqns {common} -Playering=srvn -Pmva=exact'
        elif method in ('sim', 'lqsim'):
            cmd = f'lqsim {common} -A {self.options.samples}'
        elif method == 'lqnsdefault':
            cmd = f'lqns {verbose_flag} {praqma_flag} -x {filename}'
        else:  # default or lqns
            cmd = f'lqns {common}'

        # Clean up multiple spaces
        cmd = ' '.join(cmd.split())
        return cmd

    def _run_remote_lqns(self, lqnx_file: str) -> None:
        """
        Execute LQNS via remote REST API.

        Args:
            lqnx_file: Path to the LQNX model file
        """
        # Read LQNX model file
        with open(lqnx_file, 'r', encoding='utf-8') as f:
            model_content = f.read()

        # Determine endpoint based on method
        base_url = self.options.remote_url.rstrip('/')
        method = self.options.method.lower()
        if method in ('sim', 'lqsim'):
            endpoint = f"{base_url}/api/v1/solve/lqsim"
        else:
            endpoint = f"{base_url}/api/v1/solve/lqns"

        # Build request
        request_data = {
            'model': {
                'content': model_content,
                'base64': False
            },
            'options': {
                'include_raw_output': True
            }
        }

        # Add pragmas based on method
        if method not in ('sim', 'lqsim'):
            multiserver = self.options.multiserver or 'rolia'
            if multiserver == 'default':
                multiserver = 'rolia'

            request_data['options']['pragmas'] = {
                'multiserver': multiserver,
                'stop_on_message_loss': False
            }

            # Method-specific pragmas
            if method == 'srvn':
                request_data['options']['pragmas']['layering'] = 'srvn'
            elif method == 'exactmva':
                request_data['options']['pragmas']['mva'] = 'exact'
            elif method == 'srvn.exactmva':
                request_data['options']['pragmas']['layering'] = 'srvn'
                request_data['options']['pragmas']['mva'] = 'exact'
        else:
            # LQSIM options
            request_data['options']['blocks'] = 30
            if self.options.samples > 0:
                request_data['options']['run_time'] = self.options.samples

        # Make HTTP request
        json_data = json.dumps(request_data).encode('utf-8')
        req = urllib.request.Request(
            endpoint,
            data=json_data,
            headers={
                'Content-Type': 'application/json',
                'Accept': 'application/json'
            },
            method='POST'
        )

        try:
            with urllib.request.urlopen(req, timeout=300) as response:
                response_data = json.loads(response.read().decode('utf-8'))
        except urllib.error.HTTPError as e:
            error_body = e.read().decode('utf-8') if e.fp else ''
            raise RuntimeError(f"Remote LQNS returned error: {e.code} - {error_body}")
        except urllib.error.URLError as e:
            raise RuntimeError(f"Remote LQNS connection failed: {e.reason}")

        # Check response status
        status = response_data.get('status', '')
        if status in ('error', 'failed'):
            error_msg = response_data.get('error', 'Unknown error')
            raise RuntimeError(f"Remote solver returned error: {error_msg}")

        # Extract and write LQXO content
        raw_output = response_data.get('raw_output', {})
        lqxo_content = raw_output.get('lqxo', '') if isinstance(raw_output, dict) else ''

        if lqxo_content:
            lqxo_file = lqnx_file.replace('.lqnx', '.lqxo')
            with open(lqxo_file, 'w', encoding='utf-8') as f:
                f.write(lqxo_content)
        else:
            raise RuntimeError("Remote solver did not return LQXO output")

    def _parse_xml_results(self, filename: str) -> None:
        """
        Parse LQNS XML output file (.lqxo).

        Extracts utilization, service times, throughputs, etc. from the output.
        """
        # Get LayeredNetworkStruct for node info
        self._lqn = self.getStruct()

        if self._lqn is None:
            raise RuntimeError("Cannot get LayeredNetworkStruct from model")

        # Wait for file to exist
        import time
        max_wait = 10  # seconds
        waited = 0
        while not os.path.exists(filename) and waited < max_wait:
            time.sleep(0.01)
            waited += 0.01

        if not os.path.exists(filename):
            raise RuntimeError(f"LQNS output file not found: {filename}")

        # Parse XML
        tree = ET.parse(filename)
        root = tree.getroot()

        # Try to get names from model first
        names = self._get_node_names()
        types = self._get_node_types()

        # If names is empty, build names mapping from the XML itself
        if not names:
            names, types = self._build_names_from_xml(root)

        # Get structure info
        try:
            # Native Python struct has nidx and ncalls as attributes
            num_nodes = int(self._lqn.nidx)
            num_calls = int(self._lqn.ncalls) if hasattr(self._lqn, 'ncalls') else 0
        except:
            # Use the number of names we found
            num_nodes = max(names.keys()) if names else 100
            num_calls = 100

        # Initialize result matrices
        utilization = np.full(num_nodes, np.nan)
        phase1_util = np.full(num_nodes, np.nan)
        phase2_util = np.full(num_nodes, np.nan)
        phase1_st = np.full(num_nodes, np.nan)
        phase2_st = np.full(num_nodes, np.nan)
        throughput = np.full(num_nodes, np.nan)
        proc_waiting = np.full(num_nodes, np.nan)
        proc_util = np.full(num_nodes, np.nan)
        edges_waiting = np.full(num_calls, np.nan)

        iterations = 0

        # Store names for later use in getAvgTable
        self._parsed_names = names

        # Parse solver-params for iterations
        for solver_params in root.findall('.//solver-params'):
            for result_general in solver_params.findall('result-general'):
                iter_str = result_general.get('iterations', '0')
                try:
                    iterations = int(iter_str)
                except ValueError:
                    pass

        # Parse processors
        for proc_elem in root.findall('.//processor'):
            proc_name = proc_elem.get('name')
            proc_pos = self._find_lqn_elem(names, types, proc_name, LQN_HOST)

            # Get processor utilization
            for proc_result in proc_elem.findall('result-processor'):
                util_str = proc_result.get('utilization', '')
                if util_str and proc_pos >= 0 and proc_pos < len(proc_util):
                    proc_util[proc_pos] = float(util_str)

            # Parse tasks
            for task_elem in proc_elem.findall('task'):
                task_name = task_elem.get('name')
                task_pos = self._find_lqn_elem(names, types, task_name, LQN_TASK)

                for task_result in task_elem.findall('result-task'):
                    if task_pos >= 0 and task_pos < len(utilization):
                        utilization[task_pos] = float(task_result.get('utilization', 0))
                        p1u = task_result.get('phase1-utilization', '')
                        if p1u:
                            phase1_util[task_pos] = float(p1u)
                        p2u = task_result.get('phase2-utilization', '')
                        if p2u:
                            phase2_util[task_pos] = float(p2u)
                        throughput[task_pos] = float(task_result.get('throughput', 0))
                        proc_util[task_pos] = float(task_result.get('proc-utilization', 0))

                # Parse entries
                for entry_elem in task_elem.findall('.//entry'):
                    entry_name = entry_elem.get('name')
                    entry_pos = self._find_lqn_elem(names, types, entry_name, LQN_ENTRY)

                    for entry_result in entry_elem.findall('result-entry'):
                        if entry_pos >= 0 and entry_pos < len(utilization):
                            utilization[entry_pos] = float(entry_result.get('utilization', 0))
                            p1u = entry_result.get('phase1-utilization', '')
                            if p1u:
                                phase1_util[entry_pos] = float(p1u)
                            p2u = entry_result.get('phase2-utilization', '')
                            if p2u:
                                phase2_util[entry_pos] = float(p2u)
                            p1st = entry_result.get('phase1-service-time', '')
                            if p1st:
                                phase1_st[entry_pos] = float(p1st)
                            p2st = entry_result.get('phase2-service-time', '')
                            if p2st:
                                phase2_st[entry_pos] = float(p2st)
                            throughput[entry_pos] = float(entry_result.get('throughput', 0))
                            proc_util[entry_pos] = float(entry_result.get('proc-utilization', 0))

                    # Parse entry-phase-activities (PH1PH2 format): fill the
                    # phase activity rows, otherwise left NaN (JLINE parity)
                    entry_tput = throughput[entry_pos] if 0 <= entry_pos < len(throughput) else np.nan
                    for epa in entry_elem.findall('entry-phase-activities'):
                        for activity in epa.findall('activity'):
                            act_name = activity.get('name')
                            act_pos = self._find_lqn_elem(names, types, act_name, LQN_ACTIVITY)
                            if act_pos < 0 or act_pos >= len(utilization):
                                continue
                            for act_result in activity.findall('result-activity'):
                                u = act_result.get('utilization', '')
                                if u:
                                    utilization[act_pos] = float(u)
                                st = act_result.get('service-time', '')
                                if st:
                                    phase1_st[act_pos] = float(st)
                                pw = act_result.get('proc-waiting', '')
                                if pw:
                                    proc_waiting[act_pos] = float(pw)
                                t = act_result.get('throughput', '')
                                if t:
                                    throughput[act_pos] = float(t)
                                else:
                                    # LQNS omits throughput here; each phase executes once per entry invocation
                                    throughput[act_pos] = entry_tput
                                pu = act_result.get('proc-utilization', '')
                                hd = activity.get('host-demand-mean', '')
                                if pu:
                                    proc_util[act_pos] = float(pu)
                                elif hd and not np.isnan(entry_tput):
                                    # LQNS omits proc-util here: per-phase value = entry throughput * phase host demand
                                    proc_util[act_pos] = entry_tput * float(hd)

        # Parse task-activities
        for task_acts in root.findall('.//task-activities'):
            for activity in task_acts.findall('activity'):
                # Activities under task-activities (already filtered by findall)
                act_name = activity.get('name')
                act_pos = self._find_lqn_elem(names, types, act_name, LQN_ACTIVITY)

                for act_result in activity.findall('result-activity'):
                    if act_pos >= 0 and act_pos < len(utilization):
                        utilization[act_pos] = float(act_result.get('utilization', 0))
                        st_raw = act_result.get('service-time', '')
                        # absent means UNREPORTED, not zero -- MATLAB's str2double('')
                        # is NaN and the unanimity rule below must read the same thing
                        phase1_st[act_pos] = float(st_raw) if st_raw else np.nan
                        throughput[act_pos] = float(act_result.get('throughput', 0))
                        pw = act_result.get('proc-waiting', '')
                        if pw:
                            proc_waiting[act_pos] = float(pw)
                        proc_util[act_pos] = float(act_result.get('proc-utilization', 0))

        # Processor utilization of an entry, aggregated from its activity graph.
        # lqns credits host work to whichever level carries the host demand: in
        # the activity-graph form an entry declares none, so lqns reports
        # result-entry proc-utilization as a literal 0 and the work sits on the
        # result-activity rows. The entry value is the sum over the activities
        # reachable from the entry within its own task, which is what actsof
        # holds. In PH1PH2 form the same sum runs over the phase activities and
        # reproduces the value lqns reports there, so no form test is needed. An
        # entry with no activities, or any activity lqns left unreported, keeps
        # the raw attribute rather than a partial sum.
        lqn = self._lqn
        if lqn is not None and getattr(lqn, 'actsof', None):
            for eoff in range(int(lqn.nentries)):
                eidx = int(lqn.eshift) + eoff
                acts = [a for a in lqn.actsof.get(eidx, []) if a < num_nodes]
                if not acts or len(acts) != len(lqn.actsof.get(eidx, [])):
                    continue
                pu_acts = proc_util[acts]
                if not np.any(np.isnan(pu_acts)):
                    proc_util[eidx] = float(np.sum(pu_acts))

        # Phase-1 service time of an entry lqns never invoked.
        # lqns omits phase1-service-time from result-entry exactly when the entry's
        # throughput is zero: nothing was served, so there is no per-invocation mean
        # to report. LINE then carried a NaN where the table says an entry HAS a
        # response time and every other solver reports one, breaking the NaN mask --
        # see _kb/06-solver-catalog.md. The value is taken from the activity rows,
        # and ONLY where they are unanimous: if every activity reachable from the
        # entry reports a zero service time then every aggregation law agrees on
        # zero -- the serial sum, the branch-weighted mean of an OrFork, the order
        # statistic of an AndFork -- so the derivation does not depend on which one
        # applies.
        # It is deliberately NOT generalised the way proc_util is above.
        # Utilizations add over an activity graph; response times do not. Measured
        # over the example corpus, sum(actsof) reproduces phase1-service-time on
        # serial chains only and misses it wherever the graph branches
        # (lqn_workflows `Entry`: 12.5667 reported against 8.5667 summed,
        # lqn_fork_open_arrival `SE`: 0.841667 against 1.0), so a summed fallback
        # would answer with a number lqns contradicts. An entry whose activities are
        # unreported, absent, or not all zero keeps NaN.
        if lqn is not None and getattr(lqn, 'actsof', None):
            for eoff in range(int(lqn.nentries)):
                eidx = int(lqn.eshift) + eoff
                if not np.isnan(phase1_st[eidx]):
                    continue
                acts = [a for a in lqn.actsof.get(eidx, []) if a < num_nodes]
                if not acts or len(acts) != len(lqn.actsof.get(eidx, [])):
                    continue
                st_acts = phase1_st[acts]
                if not np.any(np.isnan(st_acts)) and np.all(st_acts == 0):
                    phase1_st[eidx] = 0.0

        # Build result
        self._result = LQNSResult(
            PN=proc_util,
            SN=phase1_st,
            TN=throughput,
            UN=utilization,
            RN=proc_waiting,
            QN=np.full_like(proc_waiting, np.nan),
            AN=np.full(num_nodes, np.nan),
            WN=np.full(num_nodes, np.nan),
            runtime=0.0,
            method=self.options.method,
            iterations=iterations
        )

    def _build_names_from_xml(self, root) -> Tuple[Dict[int, str], Dict[int, int]]:
        """
        Build the node name and node kind mappings by walking the .lqxo itself.

        This is the fallback used when the struct carries no names. The kind is
        the tag being walked, so it is known exactly here and is returned
        alongside the name; see _find_lqn_elem for why the pair is the key.
        """
        names = {}
        types = {}
        idx = 0  # 0-based, in step with the struct index space

        # Parse processors
        for proc_elem in root.findall('.//processor'):
            proc_name = proc_elem.get('name')
            if proc_name:
                names[idx] = proc_name
                types[idx] = LQN_HOST
                idx += 1

            # Parse tasks
            for task_elem in proc_elem.findall('task'):
                task_name = task_elem.get('name')
                if task_name:
                    names[idx] = task_name
                    types[idx] = LQN_TASK
                    idx += 1

                # Parse entries
                for entry_elem in task_elem.findall('.//entry'):
                    entry_name = entry_elem.get('name')
                    if entry_name:
                        names[idx] = entry_name
                        types[idx] = LQN_ENTRY
                        idx += 1

                # Parse activities from task-activities
                for task_acts in task_elem.findall('task-activities'):
                    for activity in task_acts.findall('activity'):
                        act_name = activity.get('name')
                        if act_name:
                            names[idx] = act_name
                            types[idx] = LQN_ACTIVITY
                            idx += 1

        return names, types

    def _get_node_names(self) -> Dict[int, str]:
        """Get mapping of node index to node name."""
        names = {}
        try:
            lqn = self._lqn
            if lqn is None:
                return names
            # Native Python struct has names as numpy array
            if hasattr(lqn, 'names') and lqn.names is not None:
                # THE STRUCT INDEX SPACE IS 0-BASED (778978b66), so element 0 is
                # a real node -- the first processor -- and skipping it drops
                # that node from every result and shifts each remaining one down
                # a slot. On lqn_twotasks that lost the P1 row outright and left
                # a table that still READ correctly, because the same map labels
                # the rows it mis-indexes.
                for idx in range(len(lqn.names)):
                    name = lqn.names[idx]
                    if name is not None and name != '':
                        names[idx] = str(name)
        except Exception:
            pass
        return names

    def _get_node_types(self) -> Dict[int, int]:
        """Get mapping of node index to LayeredNetworkElement kind."""
        types = {}
        lqn = self._lqn
        if lqn is None or not hasattr(lqn, 'type') or lqn.type is None:
            return types
        # lqn.type is 0-based, in step with lqn.names.
        for idx in range(len(lqn.type)):
            types[idx] = int(lqn.type[idx])
        return types

    def _find_lqn_elem(self, names: Dict[int, str], types: Dict[int, int],
                       target: str, elem_type: int) -> int:
        """
        0-based position of the element called TARGET whose kind is ELEM_TYPE,
        or -1 when no such element exists.

        THE KIND IS PART OF THE KEY, and has to be. A LINE-generated layered
        model routinely gives a processor, its task and that task's entry the
        SAME name, and lqn.names holds all three, so a name-only lookup returns
        whichever one it meets first and the .lqxo rows for the other two are
        written into it: one result file then yields three different wrong
        answers. The document states which kind each row describes -- it is the
        tag being read -- so the ambiguity does not have to exist. On a model
        whose names are unique this agrees element for element with the
        name-only lookup it replaces.
        """
        for pos, name in names.items():
            if name == target and types.get(pos) == elem_type:
                return pos
        return -1

    def getAvg(self) -> LQNSResult:
        """Get average performance metrics."""
        if self._result is None:
            # runAnalyzer, NOT the NetworkSolver _ensureAvgResults funnel:
            # SolverLQNS extends Solver, not NetworkSolver, so that helper is
            # not inherited, and its MAP/MMPP random-environment gate has no
            # meaning for a LayeredNetwork anyway.
            self.runAnalyzer()

        # Copy result and swap QN/UN/RN
        result = LQNSResult(
            QN=self._result.UN.copy(),
            UN=self._result.PN.copy(),
            RN=self._result.SN.copy(),
            TN=self._result.TN.copy(),
            PN=self._result.PN.copy(),
            SN=self._result.SN.copy(),
            AN=self._result.AN.copy(),
            WN=self._result.WN.copy(),
            runtime=self._result.runtime,
            method=self._result.method,
            iterations=self._result.iterations
        )

        # UN is lqns' proc-utilization, verbatim for hosts, tasks and
        # activities and aggregated over the activity graph for entries, which
        # lqns itself reports as 0 in the activity-graph form. Both lqns and LN
        # report the processor utilization summed over the host's servers, so no
        # rescaling by the host multiplicity applies.
        return result

    def getAvgTable(self) -> pd.DataFrame:
        """
        Get average performance metrics table.

        Returns:
            pandas.DataFrame with layered network performance metrics
        """
        # see _kb/06-solver-catalog.md (Wrappers: "lang='java' opt-in JAR delegation")
        if getattr(self.options, 'lang', 'python') == 'java':
            from ...jar_dispatch import ln_avg_table_via_jar
            return ln_avg_table_via_jar(self)

        result = self.getAvg()

        # see _kb/06-solver-catalog.md (Wrappers: "LQNS snap-to-tenth")
        QN = _snap_to_tenth(result.QN)
        UN = _snap_to_tenth(result.UN)
        RN = _snap_to_tenth(result.RN)
        TN = _snap_to_tenth(result.TN)

        # Get node info - prefer parsed names from XML
        if hasattr(self, '_parsed_names') and self._parsed_names:
            names = self._parsed_names
        else:
            names = self._get_node_names()

        # Build rows for all nodes (preserves NaN for values not computed)
        rows = []
        for idx, name in names.items():
            i = idx
            if i < len(QN):
                # Determine node type from the model structure
                node_type = self._get_node_type(idx)

                rows.append({
                    'Node': name,
                    'NodeType': node_type,
                    'QLen': QN[i],
                    'Util': UN[i],
                    'RespT': RN[i],
                    'ResidT': np.nan,  # LQNS doesn't compute ResidT
                    'ArvR': np.nan,    # LQNS doesn't compute ArvR
                    'Tput': TN[i],
                })

        df = pd.DataFrame(rows)

        if not self._table_silent and len(df) > 0:
            print(df.to_string(index=False))

        # IndexedTable gives MATLAB-style 5-sig-fig formatting; a bare
        # DataFrame's fixed 5-decimal repr drops a sig fig below 0.1.
        from line_solver.indexed_table import IndexedTable
        return IndexedTable(df)

    def avg_table(self) -> pd.DataFrame:
        """Alias for getAvgTable() for API consistency."""
        return self.getAvgTable()

    # Alias for snake_case naming convention
    get_avg_table = avg_table

    def _get_node_type(self, idx: int) -> str:
        """Get node type name for a given index."""
        try:
            lqn = self.getStruct()
            if lqn is None:
                return 'Unknown'

            # Check if model has type information
            # Note: lqn.type is 0-based, matching lqn.names
            if hasattr(lqn, 'type') and idx < len(lqn.type):
                elem_type = int(lqn.type[idx])
                # Native Python uses integer values:
                # 0=PROCESSOR, 1=TASK, 2=ENTRY, 3=ACTIVITY
                if elem_type == 0:
                    return 'Processor'
                elif elem_type == 1:
                    if hasattr(lqn, 'isref') and lqn.isref[idx, 0]:
                        return 'RefTask'
                    return 'Task'
                elif elem_type == 2:
                    return 'Entry'
                elif elem_type == 3:
                    return 'Activity'
                elif elem_type == 4:
                    return 'Call'
        except Exception:
            pass
        return 'Unknown'

    def getRawAvgTables(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Get raw average tables including call metrics.

        Returns:
            Tuple of (avg_table, call_avg_table) DataFrames
        """
        result = self.getAvg()

        # Get node info
        if hasattr(self, '_parsed_names') and self._parsed_names:
            names = self._parsed_names
        else:
            names = self._get_node_names()

        # Build average table rows (preserves NaN for values not computed)
        rows = []
        for idx, name in names.items():
            i = idx
            if i < len(result.QN):
                node_type = self._get_node_type(idx)
                rows.append({
                    'Node': name,
                    'NodeType': node_type,
                    'QLen': result.QN[i],
                    'Util': result.UN[i],
                    'RespT': result.RN[i],
                    'ResidT': np.nan,  # LQNS doesn't compute ResidT
                    'ArvR': np.nan,    # LQNS doesn't compute ArvR
                    'Tput': result.TN[i],
                })

        avg_table = pd.DataFrame(rows)

        # Build call average table (inter-entry calls)
        # This table shows call statistics between entries
        call_rows = []
        # For now, return an empty call table - full implementation would
        # parse call statistics from LQNS output
        call_avg_table = pd.DataFrame(call_rows, columns=['From', 'To', 'CallRate', 'Wait'])

        return avg_table, call_avg_table

    # Alias
    get_raw_avg_tables = getRawAvgTables

    @staticmethod
    def getFeatureSet() -> set:
        """Get set of features supported per layer by the LQNS solver.

        Returns the canonical feature names (mirrors MATLAB
        SolverLQNS.supports and the JAR SolverLQNS.getFeatureSet).
        """
        return {
            'Sink', 'Source', 'Queue',
            'Coxian', 'Erlang', 'Exp', 'HyperExp',
            'Buffer', 'Server', 'JobSink', 'RandomSource', 'ServiceTunnel',
            'SchedStrategy_PS', 'SchedStrategy_FCFS',
            'ClosedClass',
        }

    @staticmethod
    def supports(model) -> bool:
        """Check if model is supported.

        Mirrors MATLAB SolverLQNS.supports. No supports() existed anywhere in
        the MRO, so calling it raised AttributeError and the solver had no gate.
        """
        from ...base import supports_via_featureset
        return supports_via_featureset(SolverLQNS, model)

    @staticmethod
    def defaultOptions() -> 'LQNSOptions':
        """Get default solver options.

        Returns:
            LQNSOptions with default configuration
        """
        return LQNSOptions()

    # Aliases
    run_analyzer = runAnalyzer
    get_avg = getAvg
    is_available = isAvailable
    list_valid_methods = listValidMethods
    default_options = defaultOptions
