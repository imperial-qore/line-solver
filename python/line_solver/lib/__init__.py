"""
LINE Library API - Comprehensive low-level utility functions.

This package exposes the complete set of lib/ functions,
providing direct access to foundational implementations for:

- Trace processing and analysis
- Phase-type distributions (PH, ME, APH, DPH)
- Markov chain analysis (CTMC, DTMC)
- Markovian arrival processes (MAP, MMAP)
- Moment fitting and computation
- Laplace transform inversion
- EM-based parameter estimation
- Queueing-specific analytical methods
- Advanced moment approximation techniques

Library Modules
===============

**Distribution Fitting & Analysis:**
- butools      - Comprehensive BUTools library (PH, ME, DPH, moments, Markov chains)
- kpctoolbox   - KPC toolkit (trace, MMPP, KPC fitting, MVPH, CTMC/DTMC)
- phasetype    - Phase-type specific operations (legacy, overlaps with butools)
- markov       - Markov chain utilities (legacy, overlaps with butools)
- empht        - EM-based phase-type fitting (moved to line-apps.git)

**Process Analysis:**
- trace        - Single and multi-trace analysis (legacy, overlaps with kpctoolbox)
- mvph         - Multivariate phase-type analysis (legacy, overlaps with kpctoolbox)
- smc          - Stationary Markov chain solvers (QBD, GI/M/1, M/G/1-type)

**Approximation & Fitting:**
- m3a          - 3rd moment approximation compression
- mom          - Moment-based solver
- lti          - Laplace transform inversion (Talbot, Gaverstehfest, Euler, etc.)

**Specialized Queueing:**
- qmam         - Queueing Markov analytical methods (MAP/MAP/1, PH/PH/1, etc.)

**Primary Modules**
-------------------

The main packages (butools, kpctoolbox, smc, lti, etc.) provide the most
comprehensive coverage. Legacy modules (phasetype, markov, trace, mvph)
are maintained for backward compatibility but their functions are also
available in the primary modules with consistent naming.

Usage
-----

Import specific functions:
    from line_solver.lib import butools, kpctoolbox
    alpha, A = butools.lib_butools_ph_from_moments([1.0, 2.0, 6.0])
    trace_mean = kpctoolbox.lib_kpc_trace_mean([0.5, 0.6, 0.4, 0.7])

Or import entire modules:
    import line_solver.lib.empht as empht
    result = empht.lib_empht_fit_aph(data, n_phases=2)

Available Modules
-----------------

All modules are submodules of line_solver.lib:

- line_solver.lib.butools     (128+ functions)
- line_solver.lib.kpctoolbox  (40+ functions)
- line_solver.lib.qmam        (7 functions)
- line_solver.lib.smc         (20+ functions)
- line_solver.lib.m3a         (17+ functions)
- line_solver.lib.lti         (35+ functions)
- line_solver.lib.mom         (8 functions)
- line_solver.lib.trace       (20+ functions)
- line_solver.lib.phasetype   (13 functions)
- line_solver.lib.markov      (8+ functions)
- line_solver.lib.mvph        (10+ functions)

Function Naming Convention
--------------------------

All wrapped functions follow the pattern:
    lib_<package>_<operation>[_<variant>]

Examples:
- lib_butools_ph_moments(alpha, A)
- lib_kpc_trace_mean(trace)
- lib_empht_fit_aph(data, n_phases)
- lib_qmam_ct_map_map_1_steady_state(D0, D1)

Package Structure
-----------------

All library modules are implemented in pure Python using NumPy.

Dependencies
------------

All lib functions require:
- NumPy for array handling
- line_solver package initialization
"""

# Lazy imports - modules are imported on demand
__all__ = [
    'butools',
    'kpctoolbox',
    'qmam',
    'smc',
    'lti',
    'm3a',
    'mom',
    'trace',
    'phasetype',
    'markov',
    'mvph',
]

# This allows direct access via: from line_solver.lib import butools
# Individual functions still require full module path:
#   from line_solver.lib.butools import lib_butools_ph_moments

def __getattr__(name):
    """Lazy loading of submodules."""
    if name in __all__:
        import importlib
        return importlib.import_module(f'.{name}', __package__)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
