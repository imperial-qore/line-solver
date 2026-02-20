"""
JMTIO - Unified I/O handler for JMT model files.

Provides both read and write functionality for JMT XML format files
(JSIM and JMVA). This class provides a consistent interface matching
the Java and MATLAB implementations.

Port from:

    /matlab/src/solvers/JMT/@JMTIO/

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import numpy as np
import os
from dataclasses import dataclass
from typing import Optional, Tuple, List, Dict, Any
from xml.etree import ElementTree as ET
from xml.dom import minidom

from ...sn import NetworkStruct


@dataclass
class JMTIOOptions:
    """Options for JMTIO operations."""
    samples: int = 10000
    seed: int = 23000
    max_simulated_time: float = float('inf')
    conf_int: float = 0.99
    max_rel_err: float = 0.03


class JMTIO:
    """
    Unified I/O handler for JMT files.

    Provides both read and write functionality for JMT XML format files
    (JSIM and JMVA) from LINE network models.

    Example:
        jmtio = JMTIO(options)
        # Write
        jmtio.write_jsim(sn, "model.jsim")
        # Read
        Q, U, R, T = jmtio.parse_jsim_results("model-result.jsim", sn)
    """

    def __init__(self, options: Optional[JMTIOOptions] = None):
        """
        Create a JMTIO instance.

        Args:
            options: I/O options (samples, seed, etc.)
        """
        self.options = options or JMTIOOptions()

    # ==================== WRITE METHODS ====================

    def write_jsim(self, sn: NetworkStruct, output_path: str) -> str:
        """
        Write a network model to JSIM XML format.

        Args:
            sn: Network structure
            output_path: Output file path

        Returns:
            Path to the written file
        """
        from .handler import _write_jsim_file, SolverJMTOptions

        opts = SolverJMTOptions(
            samples=self.options.samples,
            seed=self.options.seed,
            max_simulated_time=self.options.max_simulated_time,
            conf_int=self.options.conf_int,
            max_rel_err=self.options.max_rel_err
        )
        _write_jsim_file(sn, output_path, opts)
        return output_path

    # ==================== READ METHODS ====================

    def parse_jsim_results(self, result_path: str, sn: NetworkStruct) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Parse JMT simulation results from XML output.

        Args:
            result_path: Path to result file
            sn: Network structure

        Returns:
            Tuple of (Q, U, R, T) matrices:
            - Q: Queue lengths (M x K)
            - U: Utilizations (M x K)
            - R: Response times (M x K)
            - T: Throughputs (M x K)
        """
        from .handler import _parse_jsim_results
        return _parse_jsim_results(result_path, sn)

    def parse_csv_log(self, file_path: str) -> List[List[str]]:
        """
        Parse a CSV log file from JMT.

        Args:
            file_path: Path to the CSV log file

        Returns:
            List of parsed rows, each as a list of strings
        """
        data = []
        with open(file_path, 'r') as f:
            first_line = True
            for line in f:
                if first_line:
                    first_line = False
                    continue
                parts = line.strip().split(';')
                data.append(parts)
        return data

    def parse_tran_respt(self, arv_data: List[List[str]], dep_data: List[List[str]],
                          classnames: List[str]) -> Dict[int, Tuple[np.ndarray, np.ndarray]]:
        """
        Parse transient response time data from arrival and departure logs.

        Args:
            arv_data: Arrival log data
            dep_data: Departure log data
            classnames: List of class names

        Returns:
            Dict mapping class index to (timestamps, response_times) tuples
        """
        nclasses = len(classnames)

        # Build maps for job tracking
        job_arv_times: Dict[int, float] = {}
        job_classes: Dict[int, int] = {}

        # Process arrivals
        for row in arv_data:
            if len(row) >= 4:
                try:
                    timestamp = float(row[1])
                    job_id = int(float(row[2]))
                    class_name = row[3].strip()
                    class_idx = classnames.index(class_name) if class_name in classnames else -1

                    job_arv_times[job_id] = timestamp
                    job_classes[job_id] = class_idx
                except (ValueError, IndexError):
                    continue

        # Build response time lists per class
        resp_times: List[List[float]] = [[] for _ in range(nclasses)]
        arv_timestamps: List[List[float]] = [[] for _ in range(nclasses)]

        # Process departures and compute response times
        for row in dep_data:
            if len(row) >= 4:
                try:
                    dep_timestamp = float(row[1])
                    job_id = int(float(row[2]))

                    if job_id in job_arv_times:
                        arv_time = job_arv_times[job_id]
                        class_idx = job_classes[job_id]
                        resp_time = dep_timestamp - arv_time

                        if 0 <= class_idx < nclasses:
                            resp_times[class_idx].append(resp_time)
                            arv_timestamps[class_idx].append(arv_time)
                except (ValueError, IndexError):
                    continue

        # Convert to numpy arrays
        result = {}
        for r in range(nclasses):
            if resp_times[r]:
                result[r] = (
                    np.array(arv_timestamps[r]),
                    np.array(resp_times[r])
                )

        return result

    def parse_tran_state(self, arv_data: List[List[str]], dep_data: List[List[str]],
                         node_preload: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Parse transient state data from log files.

        Args:
            arv_data: Arrival log data
            dep_data: Departure log data
            node_preload: Initial state at the node

        Returns:
            Tuple of (state, evtype, evclass, evjob) arrays
        """
        timestamps = []
        event_types = []
        event_classes = []
        event_jobs = []

        # Process arrivals (event type = 1)
        for row in arv_data:
            if len(row) >= 4:
                try:
                    timestamp = float(row[1])
                    job_id = int(float(row[2]))
                    timestamps.append(timestamp)
                    event_types.append(1)
                    event_classes.append(0)
                    event_jobs.append(job_id)
                except (ValueError, IndexError):
                    continue

        # Process departures (event type = 2)
        for row in dep_data:
            if len(row) >= 4:
                try:
                    timestamp = float(row[1])
                    job_id = int(float(row[2]))
                    timestamps.append(timestamp)
                    event_types.append(2)
                    event_classes.append(0)
                    event_jobs.append(job_id)
                except (ValueError, IndexError):
                    continue

        n = len(timestamps)
        ncols = (node_preload.shape[1] + 1) if node_preload is not None else 2

        state = np.zeros((n, ncols))
        evtype = np.zeros((n, 1))
        evclass = np.zeros((n, 1))
        evjob = np.zeros((n, 1))

        for i in range(n):
            state[i, 0] = timestamps[i]
            evtype[i, 0] = event_types[i]
            evclass[i, 0] = event_classes[i]
            evjob[i, 0] = event_jobs[i]

        return state, evtype, evclass, evjob

    def parse_logs(self, sn: NetworkStruct, is_node_logged: List[bool],
                   log_path: str) -> Dict[int, Dict[int, Tuple[np.ndarray, np.ndarray]]]:
        """
        Parse JMT log files for all logged nodes.

        Args:
            sn: Network structure
            is_node_logged: Boolean array indicating which nodes have logging
            log_path: Path to log directory

        Returns:
            Nested dict: node_idx -> class_idx -> (timestamps, values)
        """
        result: Dict[int, Dict[int, Tuple[np.ndarray, np.ndarray]]] = {}

        if not log_path:
            return result

        classnames = sn.classnames if sn.classnames else [f'Class{i+1}' for i in range(sn.nclasses)]
        nodenames = sn.nodenames if sn.nodenames else [f'Node{i+1}' for i in range(sn.nnodes)]

        for ind in range(sn.nnodes):
            if is_node_logged[ind]:
                node_name = nodenames[ind]
                log_file_arv = os.path.join(log_path, f'{node_name}-Arv.csv')
                log_file_dep = os.path.join(log_path, f'{node_name}-Dep.csv')

                if os.path.exists(log_file_arv) and os.path.exists(log_file_dep):
                    try:
                        arv_data = self.parse_csv_log(log_file_arv)
                        dep_data = self.parse_csv_log(log_file_dep)
                        node_resp_t = self.parse_tran_respt(arv_data, dep_data, classnames)
                        result[ind] = node_resp_t
                    except Exception:
                        continue

        return result
