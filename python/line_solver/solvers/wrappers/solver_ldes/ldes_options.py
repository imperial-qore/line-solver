"""
LDES solver options and result dataclasses.

This module provides configuration options and result containers for the
LDES (LINE Discrete Event Simulator) solver, which runs via subprocess
calling ldes.jar.
"""

from dataclasses import dataclass, field
from typing import Optional, Dict, Any, List
import numpy as np

from ....constants import default_verbose


def _default_verbose_level() -> str:
    """Inherit the banner verbosity from GlobalConstants ('std' or 'silent')."""
    return 'std' if default_verbose() else 'silent'


@dataclass
class LDESOptions:
    """
    Configuration options for LDES solver.

    Options map to CLI flags of `java -jar ldes.jar solve`.

    Key characteristics:
    - Discrete event simulation via ldes.jar subprocess
    - Handles multiclass Jackson queueing networks
    - Supports steady-state analysis
    - Provides statistical estimates with confidence intervals
    - Uses event-count based stopping (samples = max service completions)
    - Default: 200,000 service completion events
    """

    samples: int = 200_000
    """Maximum number of service completion events (-s). Default: 200,000.
    Deprecated alias of events; when events is set it overrides this value."""

    events: Optional[int] = None
    """DES event budget (service completion events). None (default) = unset:
    the solver falls back to samples. When set, overrides samples."""

    seed: int = 23000
    """Random seed for reproducibility (--seed). Use -1 for random seed."""

    method: str = "default"
    """Simulation method (--method). Default: 'default'"""

    cnvgon: bool = False
    """Enable convergence-based stopping (--cnvgon). Default: False"""

    cnvgtol: float = 0.05
    """Convergence tolerance (--cnvgtol). Default: 0.05 (5%)"""

    tranfilter: str = "mser5"
    """Transient filter method (--tranfilter): 'mser5', 'fixed', 'none'. Default: 'mser5'"""

    mserbatch: int = 5
    """MSER batch size (--mserbatch). Only used when tranfilter='mser5'. Default: 5"""

    warmupfrac: float = 0.2
    """Warmup fraction for fixed filter (--warmupfrac). Default: 0.2 (20%)"""

    cimethod: str = "obm"
    """CI computation method (--cimethod): 'obm', 'bm', 'spectral', 'none'. Default: 'obm'.
    'spectral' uses Heidelberger-Welch spectral analysis (log-periodogram regression) to
    account for autocorrelation between batches, giving wider but more honest CIs."""

    obmoverlap: float = 0.5
    """Overlap fraction for OBM confidence intervals (--obmoverlap). Only used when
    cimethod='obm'; 0.0 reduces it to non-overlapping batch means. Default: 0.5"""

    ciminbatch: int = 10
    """Minimum batch size for CI computation (--ciminbatch); the batch size used is
    max(ciminbatch, sqrt(n)) on n observations. Default: 10"""

    ciminobs: int = 100
    """Minimum post-warmup observations below which no CI is reported (--ciminobs).
    Default: 100"""

    spectral_low_freq_frac: float = 0.25
    """Fraction of lowest frequencies for spectral log-periodogram regression
    (--spectrallowfreqfrac). Only used when cimethod='spectral'. Default: 0.25"""

    cnvgbatch: int = 20
    """Minimum batches before the first convergence check (--cnvgbatch). Only used
    when cnvgon=True. Default: 20"""

    cnvgchk: int = 0
    """Events between convergence checks (--cnvgchk); 0 selects samples/50. Only
    used when cnvgon=True. Default: 0 (auto)"""

    slotted: bool = False
    """Run on a discrete time scale (--slotted). Every sampled interarrival and
    service time must fall on the slot lattice; a non-lattice sample is an error,
    not something the engine rounds. Intra-slot events are ordered by a fixed
    phase and NHPP / processor-sharing models are rejected. Default: False"""

    slot_length: float = 1.0
    """Slot length in model time units (--slotlength). Only meaningful when
    slotted is True. Default: 1.0"""

    replications: Optional[int] = None
    """Number of replications (--replications). Default: None (use ldes.jar default)"""

    numthreads: Optional[int] = None
    """Number of threads (--numthreads). Default: None (use ldes.jar default)"""

    verbose: str = field(default_factory=_default_verbose_level)
    """Verbosity level: 'silent', 'std', 'debug'. Default inherits GlobalConstants."""

    timespan: Optional[List[float]] = None
    """Time horizon [T0, T1] for transient analysis. None = steady-state."""

    timeout: float = float('inf')
    """Wall-clock time budget in seconds (--maxtime). inf = no budget."""

    init_sol: Optional[np.ndarray] = None
    """Warm-start initial placement (--initsol): station-major vector
    [st0_cl0, st0_cl1, ..., stM-1_clK-1]. None = default initialization
    (closed jobs at their reference stations)."""

    rest_url: Optional[str] = None
    """Base URL of an LDES REST server (the imperialqore/ldes container), e.g.
    'http://localhost:8080'. When set, the model is POSTed to <url>/api/v1/solve
    instead of being solved by a local native binary or ldes.jar. None = local."""

    def copy(self) -> 'LDESOptions':
        """Create a deep copy of this options object."""
        return LDESOptions(
            samples=self.samples,
            events=self.events,
            seed=self.seed,
            method=self.method,
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
            spectral_low_freq_frac=self.spectral_low_freq_frac,
            slotted=self.slotted,
            slot_length=self.slot_length,
            replications=self.replications,
            numthreads=self.numthreads,
            verbose=self.verbose,
            timespan=list(self.timespan) if self.timespan else None,
            timeout=self.timeout,
            init_sol=np.array(self.init_sol, copy=True) if self.init_sol is not None else None,
            rest_url=self.rest_url,
        )


@dataclass
class LDESResult:
    """
    Result container for LDES solver computations.

    All metrics are stored as numpy arrays with dimensions [stations x classes].
    Populated from the JSON output of ldes.jar.
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

    # Transient time-series data [stations][classes] -> numpy array (numTimePoints x 2)
    QNt: Optional[List[List[Optional[np.ndarray]]]] = None
    """Queue length time series [stations][classes], each (numTimePoints x 2) with [value, time]"""

    UNt: Optional[List[List[Optional[np.ndarray]]]] = None
    """Utilization time series [stations][classes]"""

    TNt: Optional[List[List[Optional[np.ndarray]]]] = None
    """Throughput time series [stations][classes]"""

    t: Optional[np.ndarray] = None
    """Time vector for transient analysis (numTimePoints x 1)"""

    respTimeSamples: Optional[List[List[Optional[List[float]]]]] = None
    """Response time samples [stations][classes] -> list of individual observations"""

    # Metadata
    method: str = "default"
    """Simulation method used"""

    runtime: float = 0.0
    """Execution time in seconds"""

    converged: bool = False
    """Whether simulation converged before reaching max events/time"""

    stopping_reason: str = ""
    """Stopping reason: 'convergence', 'max_events' (service-completion budget,
    options.samples/events, reached), 'max_sim_events' (all-event cap,
    maxSimEvents, reached), or 'max_time' (wall-clock budget)."""

    timedOut: bool = False
    """True if the wall-clock time budget (options.timeout) was exceeded."""

    convergence_batches: int = 0
    """Number of batches used for final CI estimation"""

    cache_metrics: Optional[dict] = None
    """Per-cache hit/miss ratios and expected latency, keyed by cache node name."""

    state_histogram_space: Optional[np.ndarray] = None
    """Exact joint-state residence-time histogram states [nstates x nstations*nclasses]
    in the CTMC stateSpaceAggr layout (present when run with --export-histogram)."""

    state_histogram_time: Optional[np.ndarray] = None
    """Residence time per state in state_histogram_space [nstates]."""

    state_trajectory_space: Optional[np.ndarray] = None
    """Integer joint-state trajectory along the sampled path [npoints x nstations*nclasses],
    downsampled; used by get_tran_reward."""

    state_trajectory_time: Optional[np.ndarray] = None
    """Time points for state_trajectory_space [npoints]."""

    # Finite capacity region (FCR) metrics [regions x classes]
    QNfcr: Optional[np.ndarray] = None
    """FCR mean queue lengths (raw in-region job count) [regions x classes]"""

    UNfcr: Optional[np.ndarray] = None
    """FCR mean utilizations (NaN; not applicable) [regions x classes]"""

    RNfcr: Optional[np.ndarray] = None
    """FCR mean response times [regions x classes]"""

    TNfcr: Optional[np.ndarray] = None
    """FCR mean throughputs [regions x classes]"""

    ANfcr: Optional[np.ndarray] = None
    """FCR mean arrival rates [regions x classes]"""

    WNfcr: Optional[np.ndarray] = None
    """FCR mean residence times [regions x classes]"""

    WeightNfcr: Optional[np.ndarray] = None
    """FCR mean total weight (weighted occupation); row sum = region aggregate [regions x classes]"""

    MemOccNfcr: Optional[np.ndarray] = None
    """FCR mean memory occupation; row sum = region aggregate [regions x classes]"""

    DropRateJoin: Optional[np.ndarray] = None
    """Fork-Join quorum sibling-drop rate [stations x classes]: rate of forked
    siblings discarded because they arrive after a quorum/PARTIAL Join already
    fired. Non-zero only at Join station rows; folded into getAvgLossTable."""
