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
from typing import Optional, Dict, Any, List, Tuple
import pandas as pd


@dataclass
class LQNSOptions:
    """Options for the LQNS solver."""
    method: str = 'default'  # default, lqns, srvn, exactmva, srvn.exactmva, sim, lqsim
    multiserver: str = 'rolia'  # rolia, conway, etc.
    samples: int = 10000  # For simulation methods
    verbose: bool = False
    keep: bool = True  # Keep temporary files
    seed: int = 23000
    remote: bool = False  # Enable remote execution via REST API
    remote_url: str = 'http://localhost:8080'  # URL of lqns-rest server


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


class SolverLQNS:
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

        Alternatively, use remote execution via Docker:
          1. Pull and run: docker run -d -p 8080:8080 imperialqore/lqns-rest:latest
          2. Configure: options.config.remote = True; options.config.remote_url = 'http://localhost:8080'

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
        self.options = options or LQNSOptions()
        self._result: Optional[LQNSResult] = None
        self._table_silent = False
        self._lqn = None  # Will be set during analysis
        self._keep = kwargs.get('keep', True)  # For compatibility

    @staticmethod
    def isAvailable() -> bool:
        """
        Check if lqns is available in the system PATH.

        Returns:
            True if lqns is available, False otherwise
        """
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
        import time
        start_time = time.time()

        # Create temporary directory
        temp_dir = tempfile.mkdtemp(prefix='lqns_')

        try:
            # Write model to LQNX format
            lqnx_file = os.path.join(temp_dir, 'model.lqnx')

            # Check if model is native LayeredNetwork
            model_type = type(self.model).__name__
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
                if self.options.verbose:
                    print(f"Using remote LQNS at: {self.options.remote_url}")
                self._run_remote_lqns(lqnx_file)
            else:
                # Build and execute command locally
                cmd = self._build_command(lqnx_file)

                if self.options.verbose:
                    print(f"LQNS command: {cmd}")

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

        # If names is empty, build names mapping from the XML itself
        if not names:
            names = self._build_names_from_xml(root)

        # Get structure info
        try:
            num_nodes = int(self._lqn.obj.getNidx())
            num_calls = int(self._lqn.obj.getNcalls())
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
            proc_pos = self._find_name_position(names, proc_name)

            # Get processor utilization
            for proc_result in proc_elem.findall('result-processor'):
                util_str = proc_result.get('utilization', '')
                if util_str and proc_pos >= 0 and proc_pos < len(proc_util):
                    proc_util[proc_pos] = float(util_str)

            # Parse tasks
            for task_elem in proc_elem.findall('task'):
                task_name = task_elem.get('name')
                task_pos = self._find_name_position(names, task_name)

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
                    entry_pos = self._find_name_position(names, entry_name)

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

        # Parse task-activities
        for task_acts in root.findall('.//task-activities'):
            for activity in task_acts.findall('activity'):
                # Activities under task-activities (already filtered by findall)
                act_name = activity.get('name')
                act_pos = self._find_name_position(names, act_name)

                for act_result in activity.findall('result-activity'):
                    if act_pos >= 0 and act_pos < len(utilization):
                        utilization[act_pos] = float(act_result.get('utilization', 0))
                        phase1_st[act_pos] = float(act_result.get('service-time', 0))
                        throughput[act_pos] = float(act_result.get('throughput', 0))
                        pw = act_result.get('proc-waiting', '')
                        if pw:
                            proc_waiting[act_pos] = float(pw)
                        proc_util[act_pos] = float(act_result.get('proc-utilization', 0))

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

    def _build_names_from_xml(self, root) -> Dict[int, str]:
        """Build node names mapping by parsing the LQXO XML structure."""
        names = {}
        idx = 1  # 1-based indexing

        # Parse processors
        for proc_elem in root.findall('.//processor'):
            proc_name = proc_elem.get('name')
            if proc_name:
                names[idx] = proc_name
                idx += 1

            # Parse tasks
            for task_elem in proc_elem.findall('task'):
                task_name = task_elem.get('name')
                if task_name:
                    names[idx] = task_name
                    idx += 1

                # Parse entries
                for entry_elem in task_elem.findall('.//entry'):
                    entry_name = entry_elem.get('name')
                    if entry_name:
                        names[idx] = entry_name
                        idx += 1

                # Parse activities from task-activities
                for task_acts in task_elem.findall('task-activities'):
                    for activity in task_acts.findall('activity'):
                        act_name = activity.get('name')
                        if act_name:
                            names[idx] = act_name
                            idx += 1

        return names

    def _get_node_names(self) -> Dict[int, str]:
        """Get mapping of node index to node name."""
        names = {}
        try:
            names_dict = self._lqn.obj.getNames()
            # Convert Map to Python dict
            for key in java_names.keySet():
                names[int(key)] = str(java_names.get(key))
        except Exception:
            pass
        return names

    def _find_name_position(self, names: Dict[int, str], target: str) -> int:
        """Find position of a name in the names dictionary."""
        for pos, name in names.items():
            if name == target:
                return pos - 1  # 0-based index
        return -1

    def getAvg(self) -> LQNSResult:
        """Get average performance metrics."""
        if self._result is None:
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
        return result

    def getAvgTable(self) -> pd.DataFrame:
        """
        Get average performance metrics table.

        Returns:
            pandas.DataFrame with layered network performance metrics
        """
        result = self.getAvg()

        # Get node info - prefer parsed names from XML
        if hasattr(self, '_parsed_names') and self._parsed_names:
            names = self._parsed_names
        else:
            names = self._get_node_names()

        # Build rows only for nodes that have data (non-NaN values)
        rows = []
        for idx, name in names.items():
            i = idx - 1  # Convert to 0-based index
            if i < len(result.QN):
                # Only include nodes with at least some non-zero data
                qlen = result.QN[i] if not np.isnan(result.QN[i]) else 0
                util = result.UN[i] if not np.isnan(result.UN[i]) else 0
                respt = result.RN[i] if not np.isnan(result.RN[i]) else 0
                tput = result.TN[i] if not np.isnan(result.TN[i]) else 0

                rows.append({
                    'Node': name,
                    'QLen': qlen,
                    'Util': util,
                    'RespT': respt,
                    'Tput': tput,
                })

        df = pd.DataFrame(rows)

        if not self._table_silent and len(df) > 0:
            print(df.to_string(index=False))

        return df

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

        # Build average table rows
        rows = []
        for idx, name in names.items():
            i = idx - 1
            if i < len(result.QN):
                rows.append({
                    'Node': name,
                    'QLen': result.QN[i] if not np.isnan(result.QN[i]) else 0,
                    'Util': result.UN[i] if not np.isnan(result.UN[i]) else 0,
                    'RespT': result.RN[i] if not np.isnan(result.RN[i]) else 0,
                    'Tput': result.TN[i] if not np.isnan(result.TN[i]) else 0,
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
    def defaultOptions() -> 'LQNSOptions':
        """Get default solver options.

        Returns:
            LQNSOptions with default configuration
        """
        return LQNSOptions()

    # Aliases
    run_analyzer = runAnalyzer
    get_avg = getAvg
    get_avg_table = getAvgTable
    avg_table = getAvgTable
    is_available = isAvailable
    list_valid_methods = listValidMethods
    default_options = defaultOptions
