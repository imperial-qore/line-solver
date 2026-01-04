"""
DES solver options and result dataclasses.

This module provides configuration options and result containers for the
SimPy-based Discrete Event Simulation (DES) solver.
"""

from dataclasses import dataclass, field
from typing import Optional, Dict, Any
import numpy as np


@dataclass
class DESOptions:
    """
    Configuration options for Discrete Event Simulation (DES) solver.

    DES uses the SimPy library for discrete-event simulation that generates
    sample paths of the queueing network state evolution over time.

    Key characteristics:
    - Discrete event simulation approach using SimPy library
    - Handles multiclass Jackson queueing networks
    - Supports steady-state analysis
    - Provides statistical estimates with confidence intervals
    - Uses event-count based stopping (samples = max service completions)
    - Default: 200,000 service completion events

    Default Configuration (based on benchmark analysis):
        samples     = 200,000   (good accuracy/speed balance, ~5% mean error)
        tranfilter  = "mser5"   (adaptive warmup detection)
        warmupfrac  = 0.20      (20% warmup for fixed filter)
        cimethod    = "obm"     (overlapping batch means)
        obmoverlap  = 0.50      (50% overlap)
        cnvgon      = false     (fixed event count)

    For heavy-tailed workloads (high SCV), consider increasing samples to 500,000+.
    """

    # Event-count based stopping
    samples: int = 200_000
    """Maximum number of service completion events. Default: 200,000"""

    seed: int = 23000
    """Random seed for reproducibility. Use -1 for random seed."""

    # Convergence options
    cnvgon: bool = False
    """Enable convergence-based stopping. Default: False"""

    cnvgtol: float = 0.05
    """Convergence tolerance - stop when (CI half-width / mean) < this. Default: 0.05 (5%)"""

    cnvgbatch: int = 20
    """Minimum number of batches before checking convergence. Default: 20"""

    cnvgchk: int = 0
    """Events between convergence checks. 0 = auto (samples/50). Default: 0"""

    # Transient detection options
    tranfilter: str = "mser5"
    """Transient filter method: 'mser5' (auto), 'fixed' (fraction), 'none'. Default: 'mser5'"""

    mserbatch: int = 5
    """Batch size for MSER transient detection. Default: 5 (MSER-5)"""

    warmupfrac: float = 0.2
    """Warmup fraction for fixed filter (0.0 to 1.0). Default: 0.2 (20%)"""

    # Confidence interval options
    cimethod: str = "obm"
    """CI computation method: 'obm' (overlapping batch means), 'bm', 'none'. Default: 'obm'"""

    obmoverlap: float = 0.5
    """Overlap fraction for OBM (0.0 to 1.0). Default: 0.5 (50%)"""

    ciminbatch: int = 10
    """Minimum batch size for CI computation. Default: 10"""

    ciminobs: int = 100
    """Minimum post-warmup observations for CI. Default: 100"""

    confint: float = 0.95
    """Confidence level for intervals. Default: 0.95 (95%)"""

    # Output control
    verbose: str = "silent"
    """Verbosity level: 'silent', 'std', 'debug'. Default: 'silent'"""

    # Method selection
    method: str = "default"
    """Simulation method: 'default' or 'simpy'. Default: 'default'"""

    # Transient analysis options
    transient: bool = False
    """Enable transient analysis (collect time-series data). Default: False"""

    transient_interval: float = 1.0
    """Time interval between transient samples. Default: 1.0"""

    transient_max_samples: int = 10000
    """Maximum transient samples to keep. Default: 10000"""

    timespan: Optional[tuple] = None
    """
    Time horizon for transient analysis as (start_time, end_time).
    When set with finite end_time, simulation runs in transient mode
    (no warmup period, collects time-series data).
    Example: timespan=(0.0, 100.0) runs simulation for 100 time units.
    Default: None (steady-state mode)
    """

    def validate(self) -> None:
        """Validate option values and raise ValueError if invalid."""
        if self.tranfilter not in ("mser5", "fixed", "none"):
            raise ValueError("tranfilter must be 'mser5', 'fixed', or 'none'")
        if self.cimethod not in ("obm", "bm", "none"):
            raise ValueError("cimethod must be 'obm', 'bm', or 'none'")
        if self.mserbatch < 1:
            raise ValueError("mserbatch must be at least 1")
        if not (0.0 <= self.warmupfrac < 1.0):
            raise ValueError("warmupfrac must be in [0, 1)")
        if not (0.0 <= self.obmoverlap <= 1.0):
            raise ValueError("obmoverlap must be in [0, 1]")
        if self.ciminbatch < 2:
            raise ValueError("ciminbatch must be at least 2")
        if self.ciminobs < 10:
            raise ValueError("ciminobs must be at least 10")

    def copy(self) -> 'DESOptions':
        """Create a deep copy of this options object."""
        return DESOptions(
            samples=self.samples,
            seed=self.seed,
            cnvgon=self.cnvgon,
            cnvgtol=self.cnvgtol,
            cnvgbatch=self.cnvgbatch,
            cnvgchk=self.cnvgchk,
            tranfilter=self.tranfilter,
            mserbatch=self.mserbatch,
            warmupfrac=self.warmupfrac,
            cimethod=self.cimethod,
            obmoverlap=self.obmoverlap,
            ciminbatch=self.ciminbatch,
            ciminobs=self.ciminobs,
            confint=self.confint,
            verbose=self.verbose,
            method=self.method,
            transient=self.transient,
            transient_interval=self.transient_interval,
            transient_max_samples=self.transient_max_samples,
            timespan=self.timespan,
        )

    def is_transient_mode(self) -> bool:
        """
        Check if simulation should run in transient mode.

        Transient mode is enabled when:
        1. transient=True, or
        2. timespan is specified with finite end time

        Returns:
            True if transient mode is enabled
        """
        if self.transient:
            return True
        if self.timespan is not None and len(self.timespan) >= 2:
            return np.isfinite(self.timespan[1])
        return False


@dataclass
class DESResult:
    """
    Result container for Discrete Event Simulation (DES) solver computations.

    All metrics are stored as numpy arrays with dimensions [stations x classes].
    Confidence interval half-widths have the same dimensions as corresponding metrics.

    Performance Metrics:
        QN: Average queue lengths (number of jobs)
        UN: Server utilization (fraction of busy time)
        RN: Response times (time from arrival to departure)
        TN: Throughput (jobs completed per time unit)
        AN: Arrival rates at each station
        WN: Residence times (visit count * response time)
        CN: Visit counts (for residence time calculations)
        XN: System throughput per class

    Transient Analysis Results (if applicable):
        QNt: Queue length evolution over time
        UNt: Utilization evolution over time
        TNt: Throughput evolution over time
        t: Time points for transient measurements
    """

    # Mean performance metrics [stations x classes]
    QN: Optional[np.ndarray] = None
    """Average queue lengths [stations x classes]"""

    UN: Optional[np.ndarray] = None
    """Server utilizations [stations x classes]"""

    RN: Optional[np.ndarray] = None
    """Response times [stations x classes]"""

    TN: Optional[np.ndarray] = None
    """Throughputs [stations x classes]"""

    AN: Optional[np.ndarray] = None
    """Arrival rates [stations x classes]"""

    WN: Optional[np.ndarray] = None
    """Residence times [stations x classes]"""

    CN: Optional[np.ndarray] = None
    """Visit counts [stations x classes]"""

    XN: Optional[np.ndarray] = None
    """System throughput per class [1 x classes]"""

    # Confidence interval half-widths [stations x classes]
    QNCI: Optional[np.ndarray] = None
    """CI half-widths for queue lengths"""

    UNCI: Optional[np.ndarray] = None
    """CI half-widths for utilizations"""

    RNCI: Optional[np.ndarray] = None
    """CI half-widths for response times"""

    TNCI: Optional[np.ndarray] = None
    """CI half-widths for throughputs"""

    ANCI: Optional[np.ndarray] = None
    """CI half-widths for arrival rates"""

    WNCI: Optional[np.ndarray] = None
    """CI half-widths for residence times"""

    # Relative precision (CI half-width / mean) [stations x classes]
    QNRelPrec: Optional[np.ndarray] = None
    """Relative precision for queue lengths"""

    UNRelPrec: Optional[np.ndarray] = None
    """Relative precision for utilizations"""

    RNRelPrec: Optional[np.ndarray] = None
    """Relative precision for response times"""

    TNRelPrec: Optional[np.ndarray] = None
    """Relative precision for throughputs"""

    # Transient analysis results
    QNt: Optional[np.ndarray] = None
    """Queue lengths over time [time_points x stations x classes]"""

    UNt: Optional[np.ndarray] = None
    """Utilizations over time [time_points x stations x classes]"""

    TNt: Optional[np.ndarray] = None
    """Throughputs over time [time_points x stations x classes]"""

    t: Optional[np.ndarray] = None
    """Time points for transient measurements"""

    # Metadata
    method: str = "simpy"
    """Simulation method used"""

    runtime: float = 0.0
    """Execution time in seconds"""

    converged: bool = False
    """Whether simulation converged before reaching max events/time"""

    stopping_reason: str = ""
    """Stopping reason: 'convergence', 'max_events', or 'max_time'"""

    convergence_batches: int = 0
    """Number of batches used for final CI estimation"""

    mser_truncation_batch: int = 0
    """MSER-5 truncation batch index (0 = no truncation)"""

    warmup_end_time: float = 0.0
    """Simulation time when warmup ended"""

    total_events: int = 0
    """Total number of service completion events processed"""

    simulation_time: float = 0.0
    """Total simulation time (virtual time, not wall clock)"""

    # Response time samples for CDF computation
    resp_time_samples: Optional[Dict[int, Dict[int, list]]] = None
    """Response time samples [station_idx][class_idx] -> list of samples"""

    pass_time_samples: Optional[Dict[int, Dict[int, list]]] = None
    """Passage time samples [station_idx][class_idx] -> list of samples (same as resp_time for now)"""
