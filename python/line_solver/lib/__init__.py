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

**Process Analysis:**
- smc          - Stationary Markov chain solvers (QBD, GI/M/1, M/G/1-type)
- mapdist      - Distance measures between MAPs and D-MAPs

**Approximation & Fitting:**
- m3a          - 3rd moment approximation compression
- mom          - Moment-based solver
- lti          - Laplace transform inversion (Talbot, Gaverstehfest, Euler, etc.)

**Specialized Queueing:**
- qmam         - Queueing Markov analytical methods (MAP/MAP/1, PH/PH/1, etc.)

**Mean Field Approximation:**
- rmftool      - Refined Mean Field approximation (vendored from ngast/rmf_tool, MIT License)

Five thin wrappers over the retired JPype backend -- phasetype, markov, trace,
mvph and a shadowed second `smc` -- were listed here until 2026-08-17. Every one
of their bodies had lost the line that called the JAR, so none of the five had
parsed since `4becd99c4`; they are removed rather than reimplemented, because
butools and kpctoolbox already cover what they claimed. `empht` went to
line-apps.git earlier and is likewise gone.

Usage
-----

Import specific functions:
    from line_solver.lib import butools, kpctoolbox
    alpha, A = butools.lib_butools_ph_from_moments([1.0, 2.0, 6.0])
    trace_mean = kpctoolbox.lib_kpc_trace_mean([0.5, 0.6, 0.4, 0.7])

Or import entire modules:
    import line_solver.lib.lti as lti
    f = lti.lib_lti_talbot(laplace_func, t=1.0)

Available Modules
-----------------

All modules are submodules of line_solver.lib:

- line_solver.lib.butools     (128+ functions)
- line_solver.lib.kpctoolbox  (40+ functions)
- line_solver.lib.qmam        (7 functions)
- line_solver.lib.smc         (20+ functions)
- line_solver.lib.m3a         (17+ functions)
- line_solver.lib.lti         (35+ functions)
- line_solver.lib.mapdist     (MAP/D-MAP distance measures)
- line_solver.lib.rmftool     (refined mean field, vendored)

Function Naming Convention
--------------------------

All wrapped functions follow the pattern:
    lib_<package>_<operation>[_<variant>]

Examples:
- lib_butools_ph_moments(alpha, A)
- lib_kpc_trace_mean(trace)
- lib_lti_talbot(laplace_func, t)
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

# Direct submodules (LINE-authored libraries)
# These are packages under lib/ directly
__all__ = [
    'kpctoolbox',  # KPC toolkit (from matlab/lib/kpctoolbox/)
    'm3a',         # M3A compression (from matlab/lib/m3a/)
    'lti',         # Laplace transform inversion
    'mapdist',     # MAP distribution utilities
    'rmftool',     # Refined Mean Field (vendored, MIT License)
    'thirdparty',  # Third-party library ports
]

# Convenience aliases for thirdparty modules
# so that 'from line_solver.lib import butools' still works
_thirdparty_aliases = {
    'aoi': 'thirdparty.aoi',
    'butools': 'thirdparty.butools',
    'fj': 'thirdparty.fj',
    'hurst_estimators': 'thirdparty.hurst_estimators',
    'iltcme': 'thirdparty.iltcme',
    'mapmsg': 'thirdparty.mapmsg',
    'qmam': 'thirdparty.qmam',
    'smc': 'thirdparty.smc',
}

def __getattr__(name):
    """Lazy loading of submodules, with thirdparty aliases."""
    import importlib
    if name in __all__:
        return importlib.import_module(f'.{name}', __package__)
    if name in _thirdparty_aliases:
        return importlib.import_module(f'.{_thirdparty_aliases[name]}', __package__)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
