"""
Native Python implementations of LINE API algorithms.

This module provides pure Python/NumPy/SciPy implementations of the
queueing network analysis algorithms.

Implemented Modules:
    sn: Service network utilities (predicates, transforms, 69+ functions)
    qsys: Single-queue analysis (M/M/1, M/M/k, M/G/1, G/G/1 approx, 17 functions)
    mc: Markov chain analysis (CTMC, DTMC solvers, 21 functions)
    pfqn: Product-form queueing network algorithms (MVA family, 64 functions)
    mam: Matrix-analytic methods (QBD, MAP, PH, APH fitting)
    cache: Cache system analysis (TTL, LRU, FIFO, 25 functions)
    aoi: Age of Information analysis (FCFS, LCFS, LST, 23 functions)
    polling: Polling system analysis (gated, exhaustive, k-limited)
    moment: Moment conversions for discrete distributions (raw, central, factorial,
        binomial, upward-factorial, negative-binomial; 20 functions)
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
    'da',        # Decomposition-aggregation toolkit (fixed-point driver)
    'polling',   # Polling system analysis (3 functions)
    'lossn',     # Loss network analysis (3 functions)
    'measures',  # Statistical distance measures (57 functions)
    'map',       # MAP-driven queue analysis (1 function)
    'trace',     # Trace analysis (16 functions)
    'perm',      # Matrix permanent computation
    'lti',       # Laplace Transform Inversion (Euler, Talbot, Gaver-Stehfest)
    'mom',       # Method of Moments solver
    'moment',    # Moment conversions for discrete distributions (20 functions)
    'mmdp',      # Markov-Modulated Deterministic Process (fluid queues)
    'fes',       # Flow-Equivalent Server analysis (3 functions)
    'npfqn',     # Non-product-form approximations (5 functions)
    'sum',       # Summation method SUM/ESUM and closing method (2 functions)
    'lsn',       # Layered stochastic network utilities (1 function)
    'mapqn',     # MAP queueing network bounds (2 algorithms)
    'me',        # Maximum Entropy methods (1 function)
    'io',        # I/O and model transformation utilities
    'wf',        # Workflow analysis (pattern detection, AUTO integration)
]

# Import implemented modules (LINE-native API)
from . import sn
from . import qsys
from . import mc
from . import pfqn
from . import mam
from . import cache
from . import da
from . import polling
from . import lossn
from . import measures
from . import map
from . import trace
from . import perm
from . import lti
from . import mom
from . import mmdp
from . import fes
from . import npfqn
from . import sum
from . import lsn
from . import mapqn
from . import me
from . import io
from . import wf

# Re-export modules moved to lib/ for backward compatibility.
# sys.modules aliasing ensures deep imports (e.g. line_solver.api.butools.ph)
# resolve to the new lib/ locations even though api/{pkg}/ dirs were removed.
import sys as _sys
import line_solver.lib.thirdparty.aoi as aoi
import line_solver.lib.thirdparty.butools as butools
import line_solver.lib.thirdparty.fj as fj
import line_solver.lib.thirdparty.qmam as qmam
import line_solver.lib.thirdparty.smc as smc
import line_solver.lib.m3a as m3a
import line_solver.lib.kpctoolbox as kpctoolbox

for _alias, _target in [
    ('line_solver.api.aoi', 'line_solver.lib.thirdparty.aoi'),
    ('line_solver.api.butools', 'line_solver.lib.thirdparty.butools'),
    ('line_solver.api.fj', 'line_solver.lib.thirdparty.fj'),
    ('line_solver.api.qmam', 'line_solver.lib.thirdparty.qmam'),
    ('line_solver.api.smc', 'line_solver.lib.thirdparty.smc'),
    ('line_solver.api.m3a', 'line_solver.lib.m3a'),
    ('line_solver.api.kpctoolbox', 'line_solver.lib.kpctoolbox'),
]:
    if _target in _sys.modules:
        _sys.modules[_alias] = _sys.modules[_target]
