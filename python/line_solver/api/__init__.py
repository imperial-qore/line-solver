"""
Native Python implementations of LINE API algorithms.

This module provides pure Python/NumPy/SciPy implementations of the
queueing network analysis algorithms, replacing the JPype-based Java wrapper.

Implemented Modules:
    sn: Service network utilities (predicates, transforms, 69+ functions)
    qsys: Single-queue analysis (M/M/1, M/M/k, M/G/1, G/G/1 approx, 17 functions)
    mc: Markov chain analysis (CTMC, DTMC solvers, 21 functions)
    pfqn: Product-form queueing network algorithms (MVA family, 64 functions)
    mam: Matrix-analytic methods (QBD, MAP, PH, APH fitting)
    cache: Cache system analysis (TTL, LRU, FIFO, 25 functions)
    aoi: Age of Information analysis (FCFS, LCFS, LST, 23 functions)
    polling: Polling system analysis (gated, exhaustive, k-limited)
    lossn: Loss network analysis (Erlang B, Erlang C)
    measures: Statistical distance measures (KL, JS, Wasserstein, etc.)
    map: MAP-driven queue analysis (MAP/M/1-PS)
    trace: Trace analysis functions (statistics, correlation, IDI/IDC)
    mmdp: Markov-Modulated Deterministic Process (fluid queue modeling)
    npfqn: Non-product-form approximations (traffic merge/split, nonexp approx)
    lsn: Layered stochastic network utilities (max multiplicity)
    mapqn: MAP queueing network bounds (LP-based)
    me: Maximum Entropy methods (open queueing networks)
    io: I/O utilities (logging, XML I/O, code generation, model adapters)

Usage:
    from line_solver.api import qsys
    result = qsys.qsys_mm1(0.5, 1.0)

    from line_solver.api import mc
    pi = mc.ctmc_solve(Q)
"""

# Fully implemented modules
__all__ = [
    'sn',        # Service network utilities (69 functions)
    'qsys',      # Single-queue analysis (17 functions)
    'mc',        # Markov chain analysis (21 functions)
    'pfqn',      # Product-form queueing networks (64 functions)
    'mam',       # Matrix-analytic methods (14 functions)
    'cache',     # Cache system analysis (25 functions)
    'aoi',       # Age of Information analysis (23 functions)
    'polling',   # Polling system analysis (3 functions)
    'lossn',     # Loss network analysis (3 functions)
    'measures',  # Statistical distance measures (57 functions)
    'map',       # MAP-driven queue analysis (1 function)
    'trace',     # Trace analysis (16 functions)
    'kpctoolbox',  # KPC-Toolbox (Markov chains, APH, MMPP)
    'm3a',       # M3A compression algorithms for MMAPs
    'perm',      # Matrix permanent computation
    'smc',       # Structured Markov Chain solvers (QBD, M/G/1, G/I/M/1)
    'lti',       # Laplace Transform Inversion (Euler, Talbot, Gaver-Stehfest)
    'qmam',      # Queue analysis with Matrix-Analytic Methods
    'mom',       # Method of Moments solver
    'butools',   # BuTools library (cloned from GitHub)
    'mmdp',      # Markov-Modulated Deterministic Process (fluid queues)
    'fes',       # Flow-Equivalent Server analysis (3 functions)
    'fj',        # Fork-Join topology analysis (3 functions)
    'npfqn',     # Non-product-form approximations (5 functions)
    'lsn',       # Layered stochastic network utilities (1 function)
    'mapqn',     # MAP queueing network bounds (2 algorithms)
    'me',        # Maximum Entropy methods (1 function)
    'io',        # I/O and model transformation utilities
    'wf',        # Workflow analysis (pattern detection, AUTO integration)
]

# Import implemented modules
from . import sn
from . import qsys
from . import mc
from . import pfqn
from . import mam
from . import cache
from . import aoi
from . import polling
from . import lossn
from . import measures
from . import map
from . import trace
from . import kpctoolbox
from . import m3a
from . import perm
from . import smc
from . import lti
from . import qmam
from . import mom
from . import butools
from . import mmdp
from . import fes
from . import fj
from . import npfqn
from . import lsn
from . import mapqn
from . import me
from . import io
from . import wf
