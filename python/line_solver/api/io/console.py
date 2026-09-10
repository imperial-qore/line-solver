"""Running progress log of a LINE solver run (the "solver console").

The console narrates what a solver is doing while it does it: reading the
model, compiling the network structure, computing chains, visits and demands,
resolving the method, iterating, and closing with the figures of merit. Each
line carries the elapsed time since the run started::

    [   0.014s] compiling the network structure of model 'cqn'
    [   0.031s]   solving the routing DTMC for the visit ratios
    [   0.052s] recognized a closed queueing network: 3 stations, 1 class, 1 chain
    [   0.061s]   AMVA sweep 10: queue-length residual 1.11e-01, X = 1.2225

It prints no tables: the result table stays the caller's own ``getAvgTable``.
THE CONSOLE IS ``VerboseLevel.DEBUG``: it narrates exactly when the run is at
DEBUG and is silent at every lower level, and it never alters a numerical
result. There is no separate console switch -- the console was one, until it
became clear that a running progress log IS what a debug verbosity is for, and
two switches for one channel only let a session ask for DEBUG and get nothing.

Typical use: ``LineLogger().set_verbose(VerboseLevel.DEBUG)`` for the session,
or ``verbose=VerboseLevel.DEBUG`` on one solver.

Nested runs (an inner solver driven by SolverLN, SolverENV or UQ) do not
narrate: only the outermost run writes, so several inner solves cannot
interleave their lines. Use :func:`is_active` to suppress a legacy print,
:func:`owns_log` to emit.

References:
    MATLAB: matlab/src/io/LineConsole.m
"""

import re
import sys
import time

from line_solver.constants import VerboseLevel

__all__ = [
    'begin_run', 'run_scope', 'close_run', 'reset',
    'is_active', 'owns_log', 'writes', 'wanted',
    'step', 'substep', 'compile_detail', 'compiling', 'detail', 'loop', 'iter_line',
    'push_quiet', 'quiet_scope',
]

# ---------------------------------------------------------------- state

_DEPTH = 0                # nesting level of open solver runs
_ACTIVE = False           # True while an outermost run is narrating
_TAG = ''                 # solver name of the open run
_MODEL_NAME = ''          # name of the model under study
_T0 = None                # perf_counter at the start of the open run
_TSETUP = None            # seconds spent before the analysis began
_QUIET = 0                # >0 while structure-compile detail is suppressed
_MUTED = 0                # >0 while a silenced run is executing
_FORCE_DETAIL = False     # True while the run compiles its own model
_LAST_LOOP = ''           # text of the iteration loop last announced
_SHOWN = 0                # iteration lines written for the current loop
_TRUNC = False            # True once that loop's lines were cut short
_CLOCK = None             # console clock, restarted by each run's opening line
_PENDING = []             # legacy completion lines held until the run closes

# routed line_debug bookkeeping
_D_LAST = ''
_D_SHAPES = []
_D_SHAPE_COUNT = []
_D_TOTAL = 0

_NUMBER = re.compile(r'[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?')

_MAX_ITER_LINES = 30      # per loop
_MAX_DETAIL_LINES = 200   # per run
_MAX_PER_SHAPE = 3        # routed messages that differ only in their numbers


def reset():
    """Forget any open run (used after an interrupted solve)."""
    global _DEPTH, _ACTIVE, _TAG, _MODEL_NAME, _T0, _TSETUP, _QUIET, _MUTED
    global _FORCE_DETAIL, _LAST_LOOP, _SHOWN, _TRUNC
    del _PENDING[:]
    _DEPTH = 0
    _ACTIVE = False
    _TAG = ''
    _MODEL_NAME = ''
    _T0 = None
    _TSETUP = None
    _QUIET = 0
    _MUTED = 0
    _FORCE_DETAIL = False
    _LAST_LOOP = ''
    _SHOWN = 0
    _TRUNC = False
    _reset_detail()


def _reset_detail():
    global _D_LAST, _D_SHAPES, _D_SHAPE_COUNT, _D_TOTAL
    _D_LAST = ''
    _D_SHAPES = []
    _D_SHAPE_COUNT = []
    _D_TOTAL = 0


def _global_verbose():
    """The session verbosity, read from the logging singleton."""
    try:
        from line_solver.api.io.logging import LineLogger
        return LineLogger().verbose
    except Exception:
        return VerboseLevel.STD


def _is_silent(verbose):
    """True when a verbosity value asks for silence."""
    if verbose is None:
        return False
    if isinstance(verbose, VerboseLevel):
        return verbose == VerboseLevel.SILENT
    if isinstance(verbose, bool):
        return not verbose
    try:
        return int(verbose) <= VerboseLevel.SILENT.value
    except (TypeError, ValueError):
        return False


def _is_debug(level):
    """True when a verbosity LEVEL asks for DEBUG.

    A bool is not a level (see :func:`wanted`) and never reaches here from that
    path; treated as False so a stray one cannot switch the console on by
    itself.
    """
    if level is None or isinstance(level, bool):
        return False
    if isinstance(level, VerboseLevel):
        return level == VerboseLevel.DEBUG
    try:
        return int(level) >= VerboseLevel.DEBUG.value
    except (TypeError, ValueError):
        return False


def wanted(options=None):
    """Resolve whether a run with these options should narrate.

    THE CONSOLE IS DEBUG, and this is the whole rule: a run narrates when it is
    at ``VerboseLevel.DEBUG`` and at no lower level. The run's own
    ``options.verbose`` decides when it carries a LEVEL, otherwise the session
    level does; there is no third switch that could put the two out of step.

    A NATIVE SOLVER OPTION SPELLS ITS VERBOSITY AS A BOOL
    (``constants.default_verbose``), and a bool carries no third state, so it
    may VETO the console but never switch it on: False silences the run, True
    defers to the session level. Reading it numerically instead would make
    ``True == 1 == STD`` pass a ``>= DEBUG`` test only by accident of the
    encoding -- and, worse, would leave the console unreachable for every native
    solver run, since that is the only type this field ever carries.
    """
    session = _is_debug(_global_verbose())
    if options is not None:
        verbose = options.get('verbose') if isinstance(options, dict) else getattr(options, 'verbose', None)
        if verbose is not None:
            if isinstance(verbose, bool):
                return session and verbose
            return _is_debug(verbose)
    return session


def is_active():
    """True while a run is narrating (gates SUPPRESSION of legacy prints)."""
    return _ACTIVE


def defer_print(text):
    """Queue a legacy line for after the run closes.

    The console owns the log while a run narrates, so a solver's standard
    completion message is held here and written once the closing DONE line has
    gone out, reading as it would with the console off.
    """
    if not _ACTIVE:
        print(text)
        return
    if _DEPTH > 1:
        return  # a nested run does not narrate, and does not report
    _PENDING.append(text)


def owns_log():
    """True only inside the OUTERMOST open run (gates EMISSION)."""
    return _ACTIVE and _DEPTH <= 1 and _MUTED == 0


def writes():
    """True when a progress line should be printed.

    Either the outermost run is narrating, or no run is open and the console is
    switched on for the session -- the second case is what lets model
    construction narrate before any solver exists.
    """
    if _MUTED > 0:
        return False
    if _DEPTH > 0:
        return owns_log()
    return wanted()


# ---------------------------------------------------------------- emission

def _emit(indent, text):
    global _CLOCK
    if _CLOCK is None:
        _CLOCK = time.perf_counter()
    # a top-level row opens with a capital, an indented substep stays lowercase
    if not indent and text:
        text = text[0].upper() + text[1:]
    sys.stdout.write('[%8.3fs] %s%s\n' % (time.perf_counter() - _CLOCK, indent, text))


def step(fmt, *args):
    """Write one progress line."""
    if not writes():
        return
    _emit('', fmt % args if args else fmt)


def substep(fmt, *args):
    """Write one indented progress line."""
    if not writes():
        return
    _emit('  ', fmt % args if args else fmt)


def compile_detail(fmt, *args):
    """One stage line of a structure compile.

    Silenced inside a :func:`push_quiet` scope, and inside an open run, where
    the structures being compiled are those of auxiliary models (SolverLN
    layers) rather than of the model under study. :func:`begin_run` lifts the
    second rule while it compiles its own model.
    """
    if _QUIET > 0 or (_DEPTH > 0 and not _FORCE_DETAIL):
        return
    substep(fmt, *args)


def compiling(name):
    """Announce the compilation of a model structure."""
    global _CLOCK
    if _DEPTH > 0 and name != _MODEL_NAME:
        # an ensemble rebuilds the same submodel once per stage or per
        # iteration, so these go through detail() and collapse to one line
        detail("refreshing the auxiliary submodel '%s'" % name)
    else:
        if _DEPTH == 0:  # a compile outside any run opens its own timeline
            _CLOCK = time.perf_counter()
        step("compiling the network structure of model '%s'", name)


def detail(text):
    """Report a solver's own debug message as a substep.

    ``line_debug`` routes here while a run narrates. Consecutive repeats are
    dropped, at most three messages of the same SHAPE (the text with its
    numbers masked) are reported, and the channel is capped per run, since a
    message inside a loop would otherwise bury the narration.
    """
    global _D_LAST, _D_TOTAL
    if not owns_log():
        return
    text = (text or '').strip()
    if not text or text == _D_LAST:
        return
    shape = _NUMBER.sub('#', text)
    if shape in _D_SHAPES:
        idx = _D_SHAPES.index(shape)
        _D_SHAPE_COUNT[idx] += 1
        if _D_SHAPE_COUNT[idx] > _MAX_PER_SHAPE:
            return
    else:
        _D_SHAPES.append(shape)
        _D_SHAPE_COUNT.append(1)
    _D_LAST = text
    _D_TOTAL += 1
    if _D_TOTAL > _MAX_DETAIL_LINES:
        if _D_TOTAL == _MAX_DETAIL_LINES + 1:
            _emit('  ', 'further solver detail not reported')
        return
    # a message that already reads as a sentence keeps its own wording
    substep('%s', _lower_first(text))


def _lower_first(text):
    """Lowercase the first letter unless the word opens with an acronym."""
    if len(text) >= 2 and text[:2] != text[:2].upper():
        return text[0].lower() + text[1:]
    return text


def loop(fmt, *args):
    """Announce an iteration loop and reset its reporting budget.

    Re-announcing the SAME text (a solver that restarts its loop, as the fluid
    integration passes do) neither reprints the header nor refills the budget.
    """
    global _LAST_LOOP, _SHOWN, _TRUNC
    if not owns_log():
        return
    text = fmt % args if args else fmt
    if text == _LAST_LOOP:
        return
    _LAST_LOOP = text
    _SHOWN = 0
    _TRUNC = False
    _emit('', text)


def iter_line(k, fmt, *args):
    """Report iteration ``k`` of the current loop.

    Lines are decimated: the first 20 iterations report in full, then every
    10th, and the loop stops reporting after 30 lines, so that a long run
    cannot bury the rest of the narration.
    """
    global _SHOWN, _TRUNC
    if not owns_log():
        return
    if k > 20 and k % 10 != 0:
        return
    if _SHOWN >= _MAX_ITER_LINES:
        if not _TRUNC:
            _TRUNC = True
            _emit('  ', 'further iterations of this loop not reported')
        return
    _SHOWN += 1
    _emit('  ', fmt % args if args else fmt)


# ---------------------------------------------------------------- quiet scope

def push_quiet():
    """Suppress structure-compile detail until :func:`pop_quiet`."""
    global _QUIET
    _QUIET += 1


def pop_quiet():
    """End one :func:`push_quiet` scope."""
    global _QUIET
    _QUIET = max(0, _QUIET - 1)


class quiet_scope(object):
    """Context manager form of :func:`push_quiet`.

    Used where many auxiliary models are compiled in a row (the SolverLN layer
    builders): each keeps its one headline, and drops its stage breakdown.
    """

    def __enter__(self):
        push_quiet()
        return self

    def __exit__(self, *exc):
        pop_quiet()
        return False


# ---------------------------------------------------------------- lifecycle

def begin_run(solver, options):
    """Open a console run. Returns True when this call opened the outermost one.

    The caller MUST pair this with :func:`close_run` in a finally block, or use
    :func:`run_scope`.
    """
    global _DEPTH, _ACTIVE, _TAG, _MODEL_NAME, _T0, _TSETUP, _MUTED
    global _LAST_LOOP, _SHOWN, _TRUNC
    if _DEPTH == 0 and not wanted(options):
        # A run that must stay silent MUTES the console for its whole
        # duration, so that nothing it triggers -- a structure compile, a step
        # line of its own -- leaks out at depth 0.
        _MUTED += 1
        return False
    _DEPTH += 1
    if _DEPTH > 1:  # nested run: the outer analyzer owns the log
        return False
    _ACTIVE = True
    _TAG = _solver_tag(solver)
    _MODEL_NAME = _model_name(solver)
    _T0 = time.perf_counter()
    _TSETUP = None
    _LAST_LOOP = ''
    _SHOWN = 0
    _TRUNC = False
    _reset_detail()
    _opening_lines(solver, options)
    _TSETUP = time.perf_counter() - _T0
    return True


def close_run(solver=None):
    """Close the innermost open run, writing the closing lines."""
    global _DEPTH, _MUTED
    if _DEPTH <= 0:
        if _MUTED > 0:
            _MUTED = max(0, _MUTED - 1)
        return
    _DEPTH -= 1
    if _DEPTH > 0 or not _ACTIVE:
        return
    if solver is not None:
        _closing_lines(solver)
    # the standard completion message follows the closing DONE line
    for line in _PENDING:
        print(line)
    # the mute of an enclosing silenced run outlives this run's own state
    muted = _MUTED
    reset()
    _MUTED = muted


class run_scope(object):
    """Context manager that opens a console run and closes it on exit.

    Closing happens on an exception too, so a failed analysis still reports
    what it had reached.
    """

    def __init__(self, solver, options):
        self.solver = solver
        self.options = options

    def __enter__(self):
        begin_run(self.solver, self.options)
        return self

    def __exit__(self, *exc):
        close_run(self.solver)
        return False


# ---------------------------------------------------------------- narration

def _opening_lines(solver, options):
    global _CLOCK
    _CLOCK = time.perf_counter()  # each run's timeline starts at zero
    if writes():  # the opening row is set off from whatever preceded it
        sys.stdout.write('\n')
    step("LINE %s: Solver%s starting on model '%s' (lang %s)",
         _version_string(), _TAG, _MODEL_NAME, _lang(options))
    _read_model(solver)
    _compile_struct(solver)
    _recognize_model(solver)
    _report_method(solver, options)
    _presolve(solver)


def _closing_lines(solver):
    res = _result_of(solver)
    if res is None:
        _ensemble_closing(solver)
        return
    method = res.get('method', 'default')
    iters = res.get('iter', None)
    mtype = _method_type(_TAG, method)
    if iters is None or (isinstance(iters, float) and iters != iters) or iters <= 1:
        step('solved by %s (%s)', method, mtype)
    else:
        step('solved by %s (%s) in %d iterations', method, mtype, int(round(iters)))
    _result_lines(solver, res)
    total = time.perf_counter() - _T0
    if _TSETUP is None:
        step('DONE in %.4f s', total)
    else:
        step('DONE in %.4f s (setup %.4f s, analysis %.4f s)',
             total, _TSETUP, max(0.0, total - _TSETUP))


def _ensemble_closing(solver):
    errs = getattr(solver, 'maxitererr', None)
    if errs is not None and len(errs) > 0:
        nonzero = [e for e in errs if e > 0]
        last = nonzero[-1] if nonzero else 0.0
        tol = getattr(getattr(solver, 'options', None), 'iter_tol', float('nan'))
        step('fixed point reached after %d iterations, final error %.3e against tolerance %.3e',
             len(errs), last, tol)
    else:
        step('solved')
    step('DONE in %.4f s', time.perf_counter() - _T0)


def _read_model(solver):
    model = _model_of(solver)
    if model is None or not hasattr(model, 'getNumberOfNodes'):
        return
    try:
        step('reading the model: %s, %s',
             _plural(model.getNumberOfNodes(), 'node'),
             _plural(model.getNumberOfClasses(), 'job class', 'job classes'))
    except Exception:
        pass


def _compile_struct(solver):
    global _FORCE_DETAIL
    model = _model_of(solver)
    if model is None or not hasattr(model, 'getStruct'):
        return
    if _has_compiled_struct(model):
        step('network structure already compiled, reusing it')
        return
    _FORCE_DETAIL = True  # this compile is the run's own model
    try:
        model.getStruct()
    finally:
        _FORCE_DETAIL = False


def _has_compiled_struct(model):
    sn = getattr(model, 'sn', None)
    if sn is None:
        return False
    has = getattr(model, 'hasStruct', None)
    if has is None:
        return True
    return bool(has)


def _recognize_model(solver):
    sn = _struct_of(solver)
    if sn is None:
        return
    if not hasattr(sn, 'nstations'):  # layered
        if hasattr(sn, 'nhosts'):
            step('layered queueing network: %d hosts, %d tasks, %d entries, %d activities',
                 sn.nhosts, sn.ntasks, sn.nentries, sn.nacts)
        return
    kind = _model_kind(sn)
    step('recognized %s %s: %s, %s, %s', _article(kind), kind,
         _plural(sn.nstations, 'station'),
         _plural(sn.nclasses, 'class', 'classes'),
         _plural(sn.nchains, 'chain'))
    sched = _sched_mix(sn)
    if sched:
        substep('scheduling: %s', sched)
    pop = _population_line(sn)
    if pop:
        substep('populations: %s', pop)


def _report_method(solver, options):
    import numpy as np
    requested = _opt(options, 'method', 'default')
    resolved = requested
    try:
        resolved = solver.resolveMethod(options)
    except Exception:
        pass
    tol = _opt(options, 'tol', float('nan'))
    itermax = _opt(options, 'iter_max', _opt(options, 'max_iter', 0))
    cap = '' if not itermax else ', iteration cap %d' % int(itermax)
    # a wrapper whose binary owns the stopping rule carries no tolerance of
    # its own, and reporting it as nan states a parameter that does not exist
    tolstr = '' if not np.isfinite(tol) else ', tolerance %g' % tol
    if resolved == requested:
        step("method '%s'%s%s", resolved, tolstr, cap)
    else:
        step("method '%s' resolves to '%s'%s%s",
             requested, resolved, tolstr, cap)
    if _TAG in ('SSA', 'LDES', 'JMT'):
        samples = _opt(options, 'samples', None)
        if samples is not None:
            substep('sample budget %g, seed %s', samples, _opt(options, 'seed', 'unset'))


def _presolve(solver):
    sn = _struct_of(solver)
    if sn is None or not hasattr(sn, 'nstations'):
        return
    try:
        import numpy as np
        from line_solver.api.sn import sn_get_demands_chain
        from line_solver.constants import SchedStrategy
    except Exception:
        return
    if getattr(sn, 'nchains', 0) == 0 or getattr(sn, 'nstations', 0) == 0:
        return
    try:
        demands = sn_get_demands_chain(sn)
        Lchain, Nchain = demands.Lchain, demands.Nchain
    except Exception:
        return
    step('computing service demands per chain')
    sched = _sched_list(sn)
    is_delay = [_is_sched(s, 'INF') for s in sched]
    is_source = [_is_sched(s, 'EXT') for s in sched]
    Lchain = np.asarray(Lchain)
    for c in range(int(sn.nchains)):
        L = np.array(Lchain[:, c], dtype=float).ravel()
        L[~np.isfinite(L)] = 0.0
        Z = float(sum(L[i] for i in range(len(L)) if i < len(is_delay) and is_delay[i]))
        Lq = L.copy()
        for i in range(len(Lq)):
            if i < len(is_delay) and (is_delay[i] or is_source[i]):
                Lq[i] = 0.0
        if Lq.size == 0 or Lq.max() <= 0:
            substep('chain %d has no queueing demand (delay only)', c + 1)
            continue
        bidx = int(Lq.argmax())
        dmax = float(Lq[bidx])
        bname = _station_name(sn, bidx)
        nc = float(np.asarray(Nchain).ravel()[c])
        if np.isfinite(nc) and nc > 0:
            substep('chain %d closed, N = %g, Z = %g: bottleneck %s at D = %g, '
                    'so X <= %g and the knee is at N* = %.3g',
                    c + 1, nc, Z, bname, dmax, 1.0 / dmax, (float(Lq.sum()) + Z) / dmax)
        else:
            lam = _open_rate(sn, c)
            substep('chain %d open, lambda = %g: bottleneck %s at D = %g, utilization %.4g',
                    c + 1, lam, bname, dmax, lam * dmax)


def _result_lines(solver, res):
    import numpy as np
    sn = _struct_of(solver)
    X = res.get('X', None)
    C = res.get('C', None)
    if X is not None and np.size(X) > 0:
        xs = np.asarray(X, dtype=float).sum(axis=0)
        if C is not None and np.size(C) > 0:
            cs = np.asarray(C, dtype=float).sum(axis=0)
            substep('system throughput %s, system response time %s', _vec(xs), _vec(cs))
        else:
            substep('system throughput %s', _vec(xs))
    U = res.get('U', None)
    if U is not None and np.size(U) > 0 and sn is not None and hasattr(sn, 'nstations'):
        from line_solver.constants import SchedStrategy
        Ust = np.asarray(U, dtype=float).sum(axis=1)
        sched = _sched_list(sn)
        best, bidx = -np.inf, None
        for i in range(min(len(Ust), len(sched))):
            # a delay station reports jobs in service, not a busy fraction,
            # and a Source has no server at all
            if _is_sched(sched[i], 'INF') or _is_sched(sched[i], 'EXT'):
                continue
            if Ust[i] > best:
                best, bidx = Ust[i], i
        if bidx is not None:
            substep('busiest queueing station %s at utilization %.4f',
                    _station_name(sn, bidx), best)
    Q = res.get('Q', None)
    if Q is not None and np.size(Q) > 0:
        # an LQN store leaves rows at NaN (the LQNS wrapper leaves QN entirely
        # unset), so sum the rows that carry a queue length and say nothing
        # when none of them does -- a nansum over all-NaN is 0, not a total
        Qa = np.asarray(Q, dtype=float)
        if np.any(np.isfinite(Qa)):
            substep('mean jobs in the network %.4f', float(np.nansum(Qa)))
    lg = res.get('lG', None)
    if lg is not None and np.isfinite(lg):
        substep('normalizing constant log G = %.6g', float(lg))


# ---------------------------------------------------------------- helpers

def _solver_tag(solver):
    name = type(solver).__name__
    return name[6:] if name.startswith('Solver') else name


def _model_of(solver):
    return getattr(solver, 'model', None)


def _model_name(solver):
    model = _model_of(solver)
    if model is None:
        return '(unnamed)'
    for attr in ('getName', 'get_name'):
        fn = getattr(model, attr, None)
        if callable(fn):
            try:
                name = fn()
                if name:
                    return str(name)
            except Exception:
                pass
    return str(getattr(model, 'name', '(unnamed)'))


def _struct_of(solver):
    model = _model_of(solver)
    if model is None or not hasattr(model, 'getStruct'):
        return None
    try:
        return model.getStruct()
    except Exception:
        return None


def _result_of(solver):
    """The analyzer result as a dict keyed Q/U/R/T/C/X/method/iter, or None.

    The python store is a flat dict with QN/UN/RN/TN/XN/CN keys, unlike the
    MATLAB result.Avg block, so the names are mapped here rather than at every
    reader.
    """
    res = getattr(solver, 'result', None)
    if res is None:
        res = getattr(solver, '_result', None)
    if res is None:
        return None
    if not isinstance(res, dict):
        avg = getattr(res, 'Avg', None)
        if avg is not None and isinstance(avg, dict):
            res = avg
        else:
            # a dataclass store (SolverCTMCReturn and friends) carries the
            # metrics as attributes under either the short or the long name
            src = avg if avg is not None else res
            res = {}
            for k in ('Q', 'U', 'R', 'T', 'C', 'X', 'QN', 'UN', 'RN', 'TN', 'CN', 'XN',
                      'method', 'iter', 'runtime', 'lG', 'lognormconst'):
                if hasattr(src, k):
                    res[k] = getattr(src, k)
    out = {}
    for short, long in (('Q', 'QN'), ('U', 'UN'), ('R', 'RN'),
                        ('T', 'TN'), ('C', 'CN'), ('X', 'XN')):
        out[short] = res.get(short, res.get(long))
    out['method'] = res.get('method', 'default')
    out['iter'] = res.get('iter')
    out['runtime'] = res.get('runtime')
    # the normalizing constant, when the solver recorded one (NC, and MVA's
    # asymptotic estimate): reported from the result rather than from each of
    # the analyzer's several exits
    out['lG'] = res.get('lG', res.get('lognormconst'))
    return out


def _opt(options, name, default):
    if options is None:
        return default
    if isinstance(options, dict):
        val = options.get(name, default)
    else:
        val = getattr(options, name, default)
    return default if val is None else val


def _lang(options):
    return str(_opt(options, 'lang', 'python'))


def _version_string():
    try:
        from line_solver.constants import GlobalConstants
        return str(GlobalConstants.Version)
    except Exception:
        return ''


def _method_type(tag, method):
    try:
        from line_solver.solvers.base import method_type
        return method_type(tag, method)
    except Exception:
        return 'unknown'


def _model_kind(sn):
    from line_solver.constants import NodeType
    import numpy as np
    nodetype = list(getattr(sn, 'nodetype', []))
    if NodeType.Transition in nodetype or NodeType.Place in nodetype:
        return 'stochastic Petri net'
    if NodeType.Cache in nodetype:
        base = 'caching network'
    elif NodeType.Fork in nodetype:
        base = 'fork-join network'
    else:
        base = 'queueing network'
    njobs = np.asarray(getattr(sn, 'njobs', []), dtype=float).ravel()
    nopen = int(np.sum(np.isinf(njobs)))
    if nopen == 0:
        return 'closed ' + base
    if nopen == njobs.size:
        return 'open ' + base
    return 'mixed ' + base


def _sched_mix(sn):
    sched = _sched_list(sn)
    if not sched:
        return ''
    seen, parts = [], []
    for s in sched:
        if s in seen:
            continue
        seen.append(s)
        parts.append('%s x%d' % (_sched_text(s).upper(), sched.count(s)))
    return ', '.join(parts)


def _sched_list(sn):
    """sn.sched as a plain list: it is a dict keyed by station index here."""
    sched = getattr(sn, 'sched', None)
    if sched is None:
        return []
    if isinstance(sched, dict):
        return [sched[k] for k in sorted(sched.keys())]
    return list(sched)


def _is_sched(value, name):
    """True when a sched entry is the named strategy, as an enum or a value."""
    from line_solver.constants import SchedStrategy
    target = getattr(SchedStrategy, name)
    if value is target:
        return True
    try:
        return int(getattr(value, 'value', value)) == int(target.value)
    except (TypeError, ValueError):
        return False


def _sched_text(s):
    name = getattr(s, 'name', None)
    return str(name) if name is not None else str(s)


def _population_line(sn):
    import numpy as np
    njobs = np.asarray(getattr(sn, 'njobs', []), dtype=float).ravel()
    parts = []
    for r in range(njobs.size):
        if np.isinf(njobs[r]):
            parts.append('%s open' % _class_name(sn, r))
        else:
            parts.append('%s N=%g' % (_class_name(sn, r), njobs[r]))
    return ', '.join(parts)


def _open_rate(sn, c):
    import numpy as np
    from line_solver.constants import SchedStrategy
    rates = np.asarray(getattr(sn, 'rates', []), dtype=float)
    sched = _sched_list(sn)
    inchain = getattr(sn, 'inchain', None)
    if inchain is None or rates.size == 0:
        return 0.0
    lam = 0.0
    for r in np.asarray(inchain[c]).ravel().astype(int):
        for i in range(min(rates.shape[0], len(sched))):
            if _is_sched(sched[i], 'EXT') and np.isfinite(rates[i, r]):
                lam += float(rates[i, r])
    return lam


def _station_name(sn, i):
    names = getattr(sn, 'nodenames', None)
    to_node = getattr(sn, 'stationToNode', None)
    try:
        return str(names[int(to_node[i])])
    except Exception:
        return 'station %d' % (i + 1)


def _class_name(sn, r):
    names = getattr(sn, 'classnames', None)
    try:
        return str(names[r])
    except Exception:
        return 'class %d' % (r + 1)


def _article(word):
    return 'an' if word[:1].lower() in 'aeiou' else 'a'


def _plural(n, singular, plural=None):
    n = int(n)
    if plural is None:
        plural = singular + 's'
    return '%d %s' % (n, singular if n == 1 else plural)


def _vec(v):
    import numpy as np
    v = np.asarray(v, dtype=float).ravel()
    if v.size == 1:
        return '%.4f' % float(v[0])
    return '[' + ' '.join('%.4f' % float(x) for x in v) + ']'
