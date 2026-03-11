"""
Native Python implementation of QNS solver.

This module provides a native Python wrapper for the qnsolver command-line tool
that analyzes queueing networks using various multiserver approximation methods.
"""

import numpy as np
import subprocess
import tempfile
import shutil
import os
import platform
from dataclasses import dataclass, field
from typing import Optional, Dict, Any, List, Tuple
import pandas as pd

from ...api.sn import NetworkStruct, NodeType
from ...api.io.logging import line_debug
from .jmva_writer import write_jmva
from ..base import NetworkSolver


@dataclass
class QNSOptions:
    """Options for the QNS solver."""
    method: str = 'default'  # conway, rolia, zhou, suri, reiser, schmidt
    multiserver: str = 'default'
    samples: int = 10000
    verbose: bool = False
    keep: bool = False  # Keep temporary files


@dataclass
class QNSResult:
    """Result from QNS solver."""
    QN: np.ndarray = field(default_factory=lambda: np.array([]))  # Queue lengths [M x K]
    UN: np.ndarray = field(default_factory=lambda: np.array([]))  # Utilizations [M x K]
    RN: np.ndarray = field(default_factory=lambda: np.array([]))  # Response times [M x K]
    TN: np.ndarray = field(default_factory=lambda: np.array([]))  # Throughputs [M x K]
    AN: np.ndarray = field(default_factory=lambda: np.array([]))  # Arrival rates [M x K]
    WN: np.ndarray = field(default_factory=lambda: np.array([]))  # Residence times [M x K]
    CN: np.ndarray = field(default_factory=lambda: np.array([]))  # System response times [1 x C]
    XN: np.ndarray = field(default_factory=lambda: np.array([]))  # System throughputs [1 x C]
    runtime: float = 0.0
    method: str = 'default'
    iter: int = 0


class SolverQNS(NetworkSolver):
    """
    Native Python QNS solver using external qnsolver tool.

    This solver wraps the qnsolver command-line tool from the LQNS toolkit
    to analyze queueing networks using various multiserver approximation methods.

    Supported methods:
    - default: Uses Conway approximation
    - conway: Conway's approximation
    - rolia: Rolia's method
    - zhou: Zhou's approximation
    - suri: Suri's approximation
    - reiser: Reiser's method
    - schmidt: Schmidt's method

    Requirements:
        The 'qnsolver' command must be available in the system PATH.
        Install from: http://www.sce.carleton.ca/rads/lqns/

    Example:
        >>> solver = SolverQNS(sn, QNSOptions(method='conway'))
        >>> result = solver.runAnalyzer()
        >>> print(result.QN)  # Queue lengths
    """

    def __init__(self, model_or_sn, options: Optional[QNSOptions] = None, **kwargs):
        """
        Initialize the QNS solver.

        Args:
            model_or_sn: Network model or NetworkStruct containing the queueing network
            options: Optional QNSOptions configuration
            **kwargs: Additional options (method, multiserver, samples, verbose, keep)
        """
        # Handle model vs sn
        if hasattr(model_or_sn, 'getStruct'):
            # It's a Network model
            self.model = model_or_sn
            self.sn = model_or_sn.getStruct()
        else:
            # It's already a NetworkStruct
            self.model = None
            self.sn = model_or_sn

        # Handle options
        if options is not None:
            self.options = options
        else:
            # Build options from kwargs
            method = kwargs.pop('method', 'default')
            qns_kwargs = {k: v for k, v in kwargs.items() if k in ['multiserver', 'samples', 'verbose', 'keep']}
            self.options = QNSOptions(method=method, **qns_kwargs)

        self._result: Optional[QNSResult] = None

        # Map method to multiserver config
        if self.options.method in ('conway', 'rolia', 'zhou', 'suri', 'reiser', 'schmidt'):
            self.options.multiserver = self.options.method

    @staticmethod
    def isAvailable() -> bool:
        """
        Check if qnsolver is available in the system PATH.

        Returns:
            True if qnsolver is available, False otherwise
        """
        try:
            if platform.system() == 'Windows':
                result = subprocess.run(
                    ['cmd', '/c', 'qnsolver -h'],
                    capture_output=True,
                    timeout=5
                )
            else:
                result = subprocess.run(
                    ['sh', '-c', 'qnsolver -h'],
                    capture_output=True,
                    timeout=5
                )

            # Check if command was found
            output = result.stdout.decode() + result.stderr.decode()
            output_lower = output.lower()

            if 'command not found' in output_lower or 'not recognized' in output_lower:
                return False

            # If we got any output, the command exists
            return len(output) > 0 or result.returncode == 0

        except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
            return False

    @staticmethod
    def listValidMethods() -> List[str]:
        """List valid methods for the QNS solver."""
        return ['default', 'conway', 'rolia', 'zhou', 'suri', 'reiser', 'schmidt']

    def getName(self) -> str:
        """Get solver name."""
        return 'QNS'

    def get_name(self) -> str:
        """Get solver name (snake_case alias)."""
        return self.getName()

    def runAnalyzer(self) -> QNSResult:
        """
        Run the QNS analysis.

        Returns:
            QNSResult containing performance metrics

        Raises:
            RuntimeError: If qnsolver is not available or fails
        """
        import time
        start_time = time.time()

        M = self.sn.nstations
        K = self.sn.nclasses
        C = self.sn.nchains

        # Initialize result matrices
        QN = np.zeros((M, K))
        UN = np.zeros((M, K))
        RN = np.zeros((M, K))
        TN = np.zeros((M, K))
        AN = np.zeros((M, K))
        WN = np.zeros((M, K))
        CN = np.zeros((1, C))
        XN = np.zeros((1, C))

        actual_method = self.options.method
        multiserver = self.options.config.get('multiserver', 'default') if hasattr(self.options, 'config') and self.options.config else 'default'
        line_debug("QNS: starting (method=%s, multiserver=%s)", self.options.method, multiserver, options=self.options)

        # Determine which path to use (matches MATLAB runAnalyzer lines 46-97)
        from ...api.sn.predicates import sn_has_product_form
        has_pf = sn_has_product_form(self.sn)
        has_open = self.model.has_open_classes() if self.model is not None and hasattr(self.model, 'has_open_classes') else False

        if has_pf or has_open:
            # Product-form or open: use qnsolver directly
            line_debug("QNS: product-form or open model, routing to qns_analyzer", options=self.options)
            QN, UN, RN, TN, AN, WN, CN, XN, actual_method = self._run_qnsolver_path(M, K, C)
        else:
            # Non-product-form closed: convert to LQN and solve via SolverLQNS
            line_debug("QNS: non-product-form closed model, converting to LQN and using SolverLQNS", options=self.options)
            QN, UN, RN, TN, AN, WN, CN, XN, actual_method = self._run_lqns_path(M, K, C)

        runtime = time.time() - start_time

        # Handle default method naming (matches MATLAB lines 108-112)
        method = self.options.method
        if method == 'default':
            method = f'default/{actual_method}'
        else:
            method = actual_method

        self._result = QNSResult(
            QN=QN, UN=UN, RN=RN, TN=TN, AN=AN, WN=WN,
            CN=CN, XN=XN, runtime=runtime, method=method, iter=0
        )

        return self._result

    def _run_qnsolver_path(self, M, K, C):
        """Run the qnsolver (external tool) path for product-form/open networks."""
        QN = np.zeros((M, K))
        UN = np.zeros((M, K))
        RN = np.zeros((M, K))
        TN = np.zeros((M, K))
        AN = np.zeros((M, K))
        WN = np.zeros((M, K))
        CN = np.zeros((1, C))
        XN = np.zeros((1, C))

        # Create temporary directory
        temp_dir = tempfile.mkdtemp(prefix='qns_')

        try:
            # Write model to JMVA format
            model_file = os.path.join(temp_dir, 'model.jmva')
            write_jmva(self.sn, model_file, {
                'method': self.options.method,
                'samples': self.options.samples,
            })

            # Prepare paths
            result_file = os.path.join(temp_dir, 'result.jmva')
            log_file = os.path.join(temp_dir, 'console.out')

            # Build command
            cmd = self._build_command(model_file, result_file, log_file)

            if self.options.verbose:
                print(f"SolverQNS command: {cmd}")

            # Execute command
            if platform.system() == 'Windows':
                process = subprocess.run(
                    ['cmd', '/c', cmd],
                    cwd=temp_dir,
                    capture_output=True
                )
            else:
                process = subprocess.run(
                    ['sh', '-c', cmd],
                    cwd=temp_dir,
                    capture_output=True
                )

            if process.returncode != 0:
                log_content = ''
                if os.path.exists(log_file):
                    with open(log_file, 'r') as f:
                        log_content = f.read()
                raise RuntimeError(
                    f"QNS solver failed with exit code: {process.returncode}\n"
                    f"Log: {log_content}\n"
                    f"Stderr: {process.stderr.decode()}"
                )

            # Parse results
            Uchain, Qchain, Wchain, Tchain = self._parse_results(result_file, C)

            # Get demands per chain
            Lchain, STchain, Vchain, Nchain, alpha = self._get_demands_chain()

            # Calculate system throughput for each chain
            Xchain = np.zeros((1, C))
            refstat = self.sn.refstat.flatten() if self.sn.refstat is not None else np.zeros(K, dtype=int)
            for c in range(C):
                ref_idx = int(refstat[c]) if c < len(refstat) else 0
                if ref_idx < Tchain.shape[0]:
                    Xchain[0, c] = Tchain[ref_idx, c]

            # Response times per chain
            Rchain = np.nan_to_num(Wchain, nan=0.0)

            # Adjust utilizations for multi-server stations
            for i in range(M):
                if i < len(self.sn.nservers):
                    servers = self.sn.nservers[i]
                    if not np.isinf(servers) and servers > 1:
                        for c in range(C):
                            Uchain[i, c] = Uchain[i, c] / servers

            # Deaggregate chain results to station-class results
            from ...api.sn.deaggregate import sn_deaggregate_chain_results
            deagg = sn_deaggregate_chain_results(
                self.sn, Lchain, None, STchain, Vchain, alpha,
                None, None, Rchain, Tchain, None, Xchain
            )

            QN[:, :] = deagg.Q[:M, :K]
            UN[:, :] = deagg.U[:M, :K]
            RN[:, :] = deagg.R[:M, :K]
            TN[:, :] = deagg.T[:M, :K]
            WN[:, :] = deagg.R[:M, :K]
            CN = deagg.C
            XN = deagg.X

            # Calculate arrival rates from throughputs
            AN = self._get_arrival_rates(TN)

            actual_method = self._actual_method

        finally:
            # Clean up temporary directory
            if not self.options.keep:
                shutil.rmtree(temp_dir, ignore_errors=True)

        return QN, UN, RN, TN, AN, WN, CN, XN, actual_method

    def _run_lqns_path(self, M, K, C):
        """Run the QN2LQN + SolverLQNS path for non-product-form closed networks."""
        from ...api.io.converters import qn2lqn
        from ..solver_lqns.solver_lqns import SolverLQNS, LQNSOptions

        # Convert QN to LQN
        lqnmodel = qn2lqn(self.model)
        lqn = lqnmodel.getStruct()

        # Configure LQNS options (matches MATLAB runAnalyzer lines 57-78)
        actual_method = self.options.method
        multiserver = self.options.method  # method name IS the multiserver method
        if multiserver == 'default':
            multiserver = 'rolia'
            actual_method = 'rolia'

        lqns_options = LQNSOptions(method='default', multiserver=multiserver)

        # Create and run SolverLQNS
        solver = SolverLQNS(lqnmodel, lqns_options)
        avg_table = solver.getAvgTable()

        # Extract results using lqn.ashift indexing (matches MATLAB lines 81-93)
        QN = np.zeros((M, K))
        UN = np.zeros((M, K))
        RN = np.zeros((M, K))
        TN = np.zeros((M, K))
        WN = np.zeros((M, K))

        ashift = lqn.ashift if hasattr(lqn, 'ashift') else 0

        # avg_table is a DataFrame with columns: Node, NodeType, QLen, Util, RespT, ResidT, ArvR, Tput
        # Index by activity position: t = ashift + r + i * nclasses (0-based)
        nrows = len(avg_table) if avg_table is not None else 0

        for r in range(K):
            for i in range(M):
                t = ashift + r + i * K
                if t < nrows:
                    QN[i, r] = avg_table.iloc[t]['QLen']
                    if not np.isinf(self.sn.nservers[i]):
                        UN[i, r] = avg_table.iloc[t]['Util'] / self.sn.nservers[i]
                    else:
                        UN[i, r] = avg_table.iloc[t]['Util']
                    RN[i, r] = avg_table.iloc[t]['RespT']
                    WN[i, r] = avg_table.iloc[t]['ResidT'] if not np.isnan(avg_table.iloc[t]['ResidT']) else avg_table.iloc[t]['RespT']
                    TN[i, r] = avg_table.iloc[t]['Tput']

        AN = self._get_arrival_rates(TN)
        CN = np.zeros((1, C))
        XN = np.zeros((1, C))

        return QN, UN, RN, TN, AN, WN, CN, XN, actual_method

    def _build_command(self, model_file: str, result_file: str, log_file: str) -> str:
        """Build the qnsolver command line."""
        cmd_parts = ['qnsolver']
        cmd_parts.append(f'-l "{model_file}"')

        self._actual_method = self.options.method

        # Add multiserver method if needed
        if self._has_multi_server() and self.options.multiserver:
            ms = self.options.multiserver.lower()
            if ms in ('default', 'conway'):
                cmd_parts.append('-mconway')
                self._actual_method = 'conway'
            elif ms == 'reiser':
                cmd_parts.append('-mreiser')
                self._actual_method = 'reiser'
            elif ms == 'rolia':
                cmd_parts.append('-mrolia')
                self._actual_method = 'rolia'
            elif ms == 'zhou':
                cmd_parts.append('-mzhou')
                self._actual_method = 'zhou'
            # Note: 'suri' and 'schmidt' don't have CLI flags in original MATLAB

        cmd_parts.append(f'-o "{result_file}"')
        cmd_parts.append(f'> "{log_file}" 2>&1')

        return ' '.join(cmd_parts)

    def _has_multi_server(self) -> bool:
        """Check if the model has multi-server stations."""
        if self.sn.nservers is None:
            return False
        for ns in self.sn.nservers:
            if ns > 1 and not np.isinf(ns):
                return True
        return False

    def _parse_results(self, result_file: str, nchains: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Parse qnsolver output file.

        Returns:
            Tuple of (Uchain, Qchain, Wchain, Tchain) matrices
        """
        M = self.sn.nstations
        Uchain = np.zeros((M, nchains))
        Qchain = np.zeros((M, nchains))
        Wchain = np.zeros((M, nchains))
        Tchain = np.zeros((M, nchains))

        if not os.path.exists(result_file):
            raise RuntimeError(f"QNS result file not found: {result_file}")

        with open(result_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or '$' in line:
                    continue
                if ',' not in line:
                    continue

                # Parse based on number of chains
                if nchains == 1:
                    parsed = self._parse_dollar_output_single_class(line, nchains)
                else:
                    parsed = self._parse_dollar_output(line, nchains)

                if parsed is None:
                    continue

                stat_name, Q, W, U, T = parsed

                # Find station index
                station_idx = -1
                for i in range(M):
                    node_idx = int(self.sn.stationToNode[i])
                    if self.sn.nodenames[node_idx] == stat_name:
                        station_idx = i
                        break

                if station_idx != -1:
                    for c in range(nchains):
                        Qchain[station_idx, c] = Q[c]
                        Wchain[station_idx, c] = W[c]
                        Uchain[station_idx, c] = U[c]
                        Tchain[station_idx, c] = T[c]

        return Uchain, Qchain, Wchain, Tchain

    def _parse_dollar_output(self, line: str, nchains: int) -> Optional[Tuple[str, List[float], List[float], List[float], List[float]]]:
        """
        Parse multi-class output line.

        Format: Station, $Q(Chain01), ..., $Q, $R(Chain01), ..., $R, $U(Chain01), ..., $U, $X(Chain01), ..., $X
        """
        parts = line.replace(' ', '').split(',')
        expected_cols = 1 + 4 * (nchains + 1)
        if len(parts) < expected_cols:
            return None

        stat_name = parts[0]
        Q = [0.0] * nchains
        W = [0.0] * nchains
        U = [0.0] * nchains
        T = [0.0] * nchains

        ptr = 1

        # Q values
        for r in range(nchains):
            try:
                Q[r] = float(parts[ptr + r])
            except (ValueError, IndexError):
                Q[r] = 0.0
        ptr += nchains + 1  # Skip aggregate

        # R values (map to W)
        for r in range(nchains):
            try:
                W[r] = float(parts[ptr + r])
            except (ValueError, IndexError):
                W[r] = 0.0
        ptr += nchains + 1

        # U values
        for r in range(nchains):
            try:
                U[r] = float(parts[ptr + r])
            except (ValueError, IndexError):
                U[r] = 0.0
        ptr += nchains + 1

        # X values (map to T)
        for r in range(nchains):
            try:
                T[r] = float(parts[ptr + r])
            except (ValueError, IndexError):
                T[r] = 0.0

        return stat_name, Q, W, U, T

    def _parse_dollar_output_single_class(self, line: str, nchains: int) -> Optional[Tuple[str, List[float], List[float], List[float], List[float]]]:
        """
        Parse single-class output line.

        Format: Station, $Q, $R, $U, $X
        """
        parts = line.replace(' ', '').split(',')
        if len(parts) < 5:
            return None

        stat_name = parts[0]
        Q = [0.0] * nchains
        W = [0.0] * nchains
        U = [0.0] * nchains
        T = [0.0] * nchains

        ptr = 1
        for r in range(nchains):
            try:
                Q[r] = float(parts[ptr + r])
            except (ValueError, IndexError):
                Q[r] = 0.0
        ptr += 1

        for r in range(nchains):
            try:
                W[r] = float(parts[ptr + r])
            except (ValueError, IndexError):
                W[r] = 0.0
        ptr += 1

        for r in range(nchains):
            try:
                U[r] = float(parts[ptr + r])
            except (ValueError, IndexError):
                U[r] = 0.0
        ptr += 1

        for r in range(nchains):
            try:
                T[r] = float(parts[ptr + r])
            except (ValueError, IndexError):
                T[r] = 0.0

        return stat_name, Q, W, U, T

    def _get_demands_chain(self):
        """Calculate demands per chain using proper chain aggregation."""
        from ...api.sn.demands import sn_get_demands_chain
        result = sn_get_demands_chain(self.sn)
        return result.Lchain, result.STchain, result.Vchain, result.Nchain, result.alpha

    def _deaggregate_chain_results(self, Dchain, STchain, Vchain, alpha, Qchain, Uchain, Rchain, Tchain, Xchain):
        """
        Deaggregate chain-level results to station-class level.

        This converts results from [M x C] to [M x K] format.
        """
        M = self.sn.nstations
        K = self.sn.nclasses
        C = self.sn.nchains

        Q = np.zeros((M, K))
        U = np.zeros((M, K))
        R = np.zeros((M, K))
        T = np.zeros((M, K))

        # For each class, find its chain and copy results
        for k in range(K):
            for c in range(C):
                if alpha[k, c] > 0:
                    for i in range(M):
                        Q[i, k] = Qchain[i, c]
                        U[i, k] = Uchain[i, c]
                        R[i, k] = Rchain[i, c]
                        T[i, k] = Tchain[i, c]
                    break

        C_out = np.sum(R, axis=0, keepdims=True)
        X_out = Xchain

        return Q, U, R, T, C_out, X_out

    def _get_arrival_rates(self, TN: np.ndarray) -> np.ndarray:
        """Calculate arrival rates from throughputs."""
        M, K = TN.shape
        AN = np.zeros((M, K))

        # For open classes, arrival rate equals throughput at source
        for k in range(K):
            if k < len(self.sn.njobs) and np.isinf(self.sn.njobs[k]):
                # Open class - use throughput as arrival rate
                for i in range(M):
                    AN[i, k] = TN[i, k]
            else:
                # Closed class - arrival rate equals throughput
                for i in range(M):
                    AN[i, k] = TN[i, k]

        return AN

    def getAvgTable(self) -> pd.DataFrame:
        """
        Get comprehensive average performance metrics table.

        Returns:
            pandas.DataFrame with columns: Station, JobClass, QLen, Util, RespT, ResidT, ArvR, Tput
        """
        if self._result is None:
            self.runAnalyzer()

        result = self._result
        M = result.QN.shape[0]
        K = result.QN.shape[1] if len(result.QN.shape) > 1 else 1

        rows = []
        for i in range(M):
            node_idx = int(self.sn.stationToNode[i])
            station_name = self.sn.nodenames[node_idx]

            for k in range(K):
                class_name = self.sn.classnames[k] if k < len(self.sn.classnames) else f'Class{k+1}'

                qlen = result.QN[i, k] if K > 1 else result.QN[i]
                util = result.UN[i, k] if K > 1 else result.UN[i]
                respt = result.RN[i, k] if K > 1 else result.RN[i]
                residt = result.WN[i, k] if K > 1 else result.WN[i]
                arvr = result.AN[i, k] if K > 1 else result.AN[i]
                tput = result.TN[i, k] if K > 1 else result.TN[i]

                # Filter out rows where all metrics are zero (matching MATLAB behavior)
                metrics = [qlen, util, respt, residt, arvr, tput]
                has_significant_value = any(
                    (not np.isnan(v) and v > 0) for v in metrics
                )
                if not has_significant_value:
                    continue

                rows.append({
                    'Station': station_name,
                    'JobClass': class_name,
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

    # =========================================================================
    # Core Metric Getters
    # =========================================================================

    def getAvg(self):
        """Get all average metrics at once.

        Returns:
            Tuple of (Q, U, R, T, A, W)
        """
        if self._result is None:
            self.runAnalyzer()
        r = self._result
        return r.QN.copy(), r.UN.copy(), r.RN.copy(), r.TN.copy(), r.AN.copy(), r.WN.copy()

    def getAvgQLen(self) -> np.ndarray:
        """Get average queue lengths (M x K)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.QN.copy()

    def getAvgUtil(self) -> np.ndarray:
        """Get average utilizations (M x K)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.UN.copy()

    def getAvgRespT(self) -> np.ndarray:
        """Get average response times (M x K)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.RN.copy()

    def getAvgResidT(self) -> np.ndarray:
        """Get average residence times (M x K)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.WN.copy()

    def getAvgWaitT(self) -> np.ndarray:
        """Get average waiting times (M x K)."""
        if self._result is None:
            self.runAnalyzer()
        R = self._result.RN.copy()
        if hasattr(self.sn, 'rates') and self.sn.rates is not None:
            rates = np.asarray(self.sn.rates)
            S = np.zeros_like(rates)
            nonzero = rates > 0
            S[nonzero] = 1.0 / rates[nonzero]
            W = R - S
            W = np.maximum(W, 0.0)
            return W
        return R

    def getAvgTput(self) -> np.ndarray:
        """Get average throughputs (M x K)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.TN.copy()

    def getAvgArvR(self) -> np.ndarray:
        """Get average arrival rates (M x K)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.AN.copy()

    # =========================================================================
    # System-Level Methods
    # =========================================================================

    def getAvgSysRespT(self) -> np.ndarray:
        """Get system response times (1 x C)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.CN.flatten().copy()

    def getAvgSysTput(self) -> np.ndarray:
        """Get system throughputs (1 x C)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.XN.flatten().copy()

    def getAvgSys(self):
        """Get system-level average metrics.

        Returns:
            Tuple of (CN, XN) - system response times and throughputs
        """
        return self.getAvgSysRespT(), self.getAvgSysTput()

    def getAvgSysTable(self) -> pd.DataFrame:
        """Get system-level metrics as DataFrame."""
        CN, XN = self.getAvgSys()
        nchains = len(CN)
        rows = []
        for c in range(nchains):
            rows.append({
                'Chain': f'Chain{c + 1}',
                'SysRespT': CN[c],
                'SysTput': XN[c],
            })
        return pd.DataFrame(rows)

    # =========================================================================
    # Chain-Level Methods
    # =========================================================================

    def _get_chains(self) -> List[List[int]]:
        """Get chain-to-class mapping from network structure."""
        if hasattr(self.sn, 'chains') and self.sn.chains is not None:
            chains_arr = np.asarray(self.sn.chains)
            nchains = self.sn.nchains if hasattr(self.sn, 'nchains') else 1
            if chains_arr.ndim == 1:
                if len(chains_arr) == 0:
                    return [[k] for k in range(self.sn.nclasses)]
                nchains = max(nchains, int(np.max(chains_arr)) + 1)
                chains = [[] for _ in range(nchains)]
                for k in range(self.sn.nclasses):
                    if k < len(chains_arr):
                        c = int(chains_arr[k])
                        if 0 <= c < nchains:
                            chains[c].append(k)
                chains = [c for c in chains if c]
                return chains if chains else [[k for k in range(self.sn.nclasses)]]
            else:
                chains = []
                for c in range(nchains):
                    chain_classes = []
                    for k in range(self.sn.nclasses):
                        if c < chains_arr.shape[0] and k < chains_arr.shape[1]:
                            if chains_arr[c, k] > 0:
                                chain_classes.append(k)
                    chains.append(chain_classes)
                chains = [c for c in chains if c]
                return chains if chains else [[k for k in range(self.sn.nclasses)]]
        return [[k] for k in range(self.sn.nclasses)]

    def getAvgQLenChain(self) -> np.ndarray:
        """Get average queue lengths aggregated by chain."""
        Q = self.getAvgQLen()
        chains = self._get_chains()
        nstations = Q.shape[0]
        nchains = len(chains)
        QN_chain = np.zeros((nstations, nchains))
        for c, cc in enumerate(chains):
            if cc:
                QN_chain[:, c] = np.sum(Q[:, cc], axis=1)
        return QN_chain

    def getAvgUtilChain(self) -> np.ndarray:
        """Get average utilizations aggregated by chain."""
        U = self.getAvgUtil()
        chains = self._get_chains()
        nstations = U.shape[0]
        nchains = len(chains)
        UN_chain = np.zeros((nstations, nchains))
        for c, cc in enumerate(chains):
            if cc:
                UN_chain[:, c] = np.sum(U[:, cc], axis=1)
        return UN_chain

    def getAvgRespTChain(self) -> np.ndarray:
        """Get average response times aggregated by chain."""
        R = self.getAvgRespT()
        chains = self._get_chains()
        nstations = R.shape[0]
        nchains = len(chains)
        RN_chain = np.zeros((nstations, nchains))
        for c, cc in enumerate(chains):
            if cc:
                RN_chain[:, c] = np.mean(R[:, cc], axis=1)
        return RN_chain

    def getAvgResidTChain(self) -> np.ndarray:
        """Get average residence times aggregated by chain."""
        return self.getAvgRespTChain()

    def getAvgTputChain(self) -> np.ndarray:
        """Get average throughputs aggregated by chain."""
        T = self.getAvgTput()
        chains = self._get_chains()
        nstations = T.shape[0]
        nchains = len(chains)
        TN_chain = np.zeros((nstations, nchains))
        for c, cc in enumerate(chains):
            if cc:
                TN_chain[:, c] = np.sum(T[:, cc], axis=1)
        return TN_chain

    def getAvgArvRChain(self) -> np.ndarray:
        """Get average arrival rates aggregated by chain."""
        return self.getAvgTputChain()

    def getAvgChain(self):
        """Get all average metrics aggregated by chain.

        Returns:
            Tuple of (QN, UN, RN, WN, AN, TN) aggregated by chain
        """
        return (self.getAvgQLenChain(), self.getAvgUtilChain(),
                self.getAvgRespTChain(), self.getAvgResidTChain(),
                self.getAvgArvRChain(), self.getAvgTputChain())

    def getAvgChainTable(self) -> pd.DataFrame:
        """Get average metrics by chain as DataFrame."""
        QN, UN, RN, WN, AN, TN = self.getAvgChain()
        nstations, nchains = QN.shape
        rows = []
        for i in range(nstations):
            node_idx = int(self.sn.stationToNode[i])
            station_name = self.sn.nodenames[node_idx]
            for c in range(nchains):
                rows.append({
                    'Station': station_name,
                    'Chain': f'Chain{c + 1}',
                    'QLen': QN[i, c],
                    'Util': UN[i, c],
                    'RespT': RN[i, c],
                    'ResidT': WN[i, c],
                    'ArvR': AN[i, c],
                    'Tput': TN[i, c],
                })
        return pd.DataFrame(rows)

    # =========================================================================
    # Node-Level Methods
    # =========================================================================

    def getAvgNode(self):
        """Get average metrics per node.

        Returns:
            Tuple of (QNn, UNn, RNn, WNn, ANn, TNn) - node-level metrics
        """
        if self._result is None:
            self.runAnalyzer()

        sn = self.sn
        I = sn.nnodes
        M = sn.nstations
        K = sn.nclasses

        QN = self._result.QN
        UN = self._result.UN
        RN = self._result.RN
        TN = self._result.TN
        WN = self._result.WN

        QNn = np.zeros((I, K))
        UNn = np.zeros((I, K))
        RNn = np.zeros((I, K))
        WNn = np.zeros((I, K))
        TNn = np.zeros((I, K))
        ANn = np.zeros((I, K))

        for ist in range(M):
            ind = int(sn.stationToNode[ist])
            if 0 <= ind < I:
                QNn[ind, :] = QN[ist, :]
                UNn[ind, :] = UN[ist, :]
                RNn[ind, :] = RN[ist, :]
                WNn[ind, :] = WN[ist, :]
                TNn[ind, :] = TN[ist, :]
                ANn[ind, :] = self._result.AN[ist, :]

        return QNn, UNn, RNn, WNn, ANn, TNn

    def getAvgNodeTable(self) -> pd.DataFrame:
        """Get average metrics by node as DataFrame."""
        QNn, UNn, RNn, WNn, ANn, TNn = self.getAvgNode()
        sn = self.sn
        nodenames = list(sn.nodenames) if hasattr(sn, 'nodenames') and sn.nodenames else []
        rows = []
        for node_idx in range(sn.nnodes):
            node_name = nodenames[node_idx] if node_idx < len(nodenames) else f'Node{node_idx}'
            for r in range(sn.nclasses):
                class_name = sn.classnames[r] if r < len(sn.classnames) else f'Class{r + 1}'
                if (abs(QNn[node_idx, r]) < 1e-10 and abs(UNn[node_idx, r]) < 1e-10 and
                        abs(RNn[node_idx, r]) < 1e-10 and abs(ANn[node_idx, r]) < 1e-10 and
                        abs(TNn[node_idx, r]) < 1e-10):
                    continue
                rows.append({
                    'Node': node_name,
                    'JobClass': class_name,
                    'QLen': QNn[node_idx, r],
                    'Util': UNn[node_idx, r],
                    'RespT': RNn[node_idx, r],
                    'ResidT': WNn[node_idx, r],
                    'ArvR': ANn[node_idx, r],
                    'Tput': TNn[node_idx, r],
                })
        df = pd.DataFrame(rows)
        if not self._table_silent:
            print(df.to_string(index=False))
        return df

    def getAvgNodeChain(self):
        """Get average metrics by node and chain."""
        return self.getAvgChain()

    def getAvgNodeChainTable(self) -> pd.DataFrame:
        """Get average metrics by node and chain as DataFrame."""
        return self.getAvgChainTable()

    def getAvgNodeQLenChain(self) -> np.ndarray:
        """Get average queue lengths by node aggregated by chain."""
        return self.getAvgQLenChain()

    def getAvgNodeUtilChain(self) -> np.ndarray:
        """Get average utilizations by node aggregated by chain."""
        return self.getAvgUtilChain()

    def getAvgNodeRespTChain(self) -> np.ndarray:
        """Get average response times by node aggregated by chain."""
        return self.getAvgRespTChain()

    def getAvgNodeResidTChain(self) -> np.ndarray:
        """Get average residence times by node aggregated by chain."""
        return self.getAvgResidTChain()

    def getAvgNodeTputChain(self) -> np.ndarray:
        """Get average throughputs by node aggregated by chain."""
        return self.getAvgTputChain()

    def getAvgNodeArvRChain(self) -> np.ndarray:
        """Get average arrival rates by node aggregated by chain."""
        return self.getAvgArvRChain()

    # =========================================================================
    # Utility / Static Methods
    # =========================================================================

    @staticmethod
    def getFeatureSet() -> set:
        """Get supported features."""
        return {
            'Source', 'Sink', 'Queue', 'Delay', 'DelayStation',
            'Exp', 'Erlang', 'HyperExp', 'PH', 'APH', 'Cox2',
            'InfiniteServer', 'SharedServer', 'Buffer', 'Dispatcher',
            'Server', 'ServiceTunnel', 'RandomSource', 'JobSink',
            'SchedStrategy_INF', 'SchedStrategy_PS', 'SchedStrategy_FCFS',
            'RoutingStrategy_PROB', 'RoutingStrategy_RAND',
            'ClosedClass', 'OpenClass',
        }

    @staticmethod
    def supports(model) -> bool:
        """Check if model is supported."""
        try:
            if hasattr(model, 'nstations'):
                nstations = model.nstations
            elif hasattr(model, 'getNumberOfStations'):
                nstations = model.getNumberOfStations()
            else:
                return False
            if hasattr(model, 'nclasses'):
                nclasses = model.nclasses
            elif hasattr(model, 'getNumberOfClasses'):
                nclasses = model.getNumberOfClasses()
            else:
                return False
            return nstations > 0 and nclasses > 0
        except Exception:
            return False

    @staticmethod
    def defaultOptions():
        """Get default solver options."""
        return QNSOptions()

    # =========================================================================
    # Aliases
    # =========================================================================

    run_analyzer = runAnalyzer
    is_available = isAvailable
    list_valid_methods = listValidMethods
    get_avg = getAvg
    get_avg_qlen = getAvgQLen
    get_avg_util = getAvgUtil
    get_avg_respt = getAvgRespT
    get_avg_residt = getAvgResidT
    get_avg_waitt = getAvgWaitT
    get_avg_tput = getAvgTput
    get_avg_arvr = getAvgArvR
    get_avg_sys_respt = getAvgSysRespT
    get_avg_sys_tput = getAvgSysTput
    get_avg_sys = getAvgSys
    get_avg_sys_table = getAvgSysTable
    get_avg_chain = getAvgChain
    get_avg_chain_table = getAvgChainTable
    get_avg_node = getAvgNode
    get_avg_node_table = getAvgNodeTable
