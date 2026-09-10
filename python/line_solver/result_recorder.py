"""Capture the result tables a script produces, with the solver that produced them.

WHY THIS EXISTS. Cross-codebase parity is asserted against one shared golden per
example (`goldens/baselines/*.json`), keyed by SOLVER NAME. Until 2026-08-19 the
only way to recover that key was to scrape the banner an example printed above
each table -- tens of thousands of lines of regex, one dialect per codebase, and
a value truncated to whatever the printer showed. The recorder supplies the same
attribution BY CONSTRUCTION: a getter knows which solver it belongs to, which
method that solver resolved, and the full-precision values it is returning.

It is off unless asked for, and asking is one of:

    LINE_RECORD_RESULTS=1 python3 my_example.py      # env, at import time
    from line_solver.result_recorder import recorder
    recorder.enable()                                 # or in-process

With it off, `install()` is never called and no getter is wrapped, so an
ordinary run carries no cost and no behaviour change at all.

WHAT IS RECORDED. One entry per OUTERMOST getter call -- the facade layer calls
through to a native solver's getter of the same name, and recording both would
double every table. Each entry carries the solver's golden label, the method it
resolved, which view was asked for ('avg', 'node', 'chain', ...), and the table
as plain rows of Python floats. Nothing is rounded here: quantizing to the
golden's printed precision is the comparator's job, and doing it at the source
would throw away the digits a future full-precision golden needs.

COVERAGE IS SWEPT, NOT LISTED. `install()` walks every Solver subclass that has
been imported and wraps every attribute whose name is a result-table getter, so a
getter added to a solver later is recorded the day it lands. A hand-written hook
table was tried first and was wrong within the hour: `SolverMVA.getAvgSysTable`
is a class attribute REBOUND to `NetworkSolver.getAvgSysTable`, so wrapping the
base class left five solvers unrecorded while the table claimed they were
covered. An unrecorded table reads downstream as a solver that produced nothing
-- a parity failure with no defect behind it -- which is exactly the kind of gap
a list maintained by hand produces.
"""
import functools
import os
import re

__all__ = ['recorder', 'ResultRecorder', 'Record', 'install', 'solver_label',
           'solver_classes', 'METHOD_VIEW', 'GETTER_PATTERN',
           'SCALAR_GETTERS', 'OBJECT_SCALAR_GETTERS',
           'METHOD_FAMILY_LABELS', 'qualified_label', 'auto_selection_method']


def _camel_to_snake(name):
    return re.sub(r'(?<!^)(?=[A-Z])', '_', name).lower()


# Solver class name -> the label the shared goldens key that solver by. These
# are the names the printed banners carried, which is what the goldens were
# generated from, so the map is a transcription and not a new naming scheme.
SOLVER_LABELS = {
    'SolverMVA': 'MVA',
    'SolverNC': 'NC',
    'SolverCTMC': 'CTMC',
    'SolverSSA': 'SSA',
    'SolverFLD': 'FLD',
    'SolverFluid': 'FLD',
    'SolverMAM': 'MAM',
    'SolverJMT': 'JMT',
    'SolverLDES': 'LDES',
    'SolverLQNS': 'LQNS',
    'SolverQNS': 'QNS',
    'SolverBA': 'BA',
    'SolverAG': 'AG',
    'SolverAUTO': 'AUTO',
    'SolverLN': 'LN',
    'SolverENV': 'ENV',
    'SolverUQ': 'UQ',
}

# Solvers that RUN ANOTHER SOLVER, and are therefore qualified by it. The member
# is not decoration: MVA layers and NC layers are different fixed points, and an
# environment answer scored against one stage's standalone table would report the
# coupling itself as a defect.
ENSEMBLE_SOLVERS = ('LN', 'ENV', 'UQ')

# The method-family method name `findSolver()` reports -> the label the goldens key
# that family by. `SolverAUTO(model, 'mva.schmidt')` is the documented way to ask
# for one runnable (family, method) pair -- it is literally the method name in the
# Method column of `model.findSolver()` -- and a solve asked for that way is
# recorded under the QUALIFIED label `MVA:schmidt`, which is the golden key for a
# non-default method.
#
# ONLY THE AUTO SPELLING IS QUALIFIED, and that is what keeps the 209 existing
# baselines untouched: `SolverMVA(model, 'exact')` still records as plain `MVA`,
# because a bare golden key means "this family's default path AS THE EXAMPLE
# DROVE IT" and cqn_oneline pins `exact` while keying `MVA`. Re-labelling direct
# constructions would move every such golden onto a key no example produces.
METHOD_FAMILY_LABELS = {
    'mva': 'MVA', 'nc': 'NC', 'ctmc': 'CTMC', 'fluid': 'FLD', 'mam': 'MAM',
    'ba': 'BA', 'ssa': 'SSA', 'ldes': 'LDES', 'jmt': 'JMT', 'ln': 'LN',
    'env': 'ENV', 'uq': 'UQ', 'lqns': 'LQNS', 'qns': 'QNS', 'ag': 'AG',
}


def qualified_label(method):
    """`MVA:schmidt` for the method name `mva.schmidt`, or None.

    None when the method name names no family this project keys a golden by, which is
    the honest answer for a bare method name: `AUTO(model, 'exact')` asks AUTO
    to choose, and the family it chooses is not knowable from the method name.
    """
    if not method or '.' not in method:
        return None
    family, rest = method.split('.', 1)
    label = METHOD_FAMILY_LABELS.get(family.lower())
    return '%s:%s' % (label, rest) if label and rest else None


def auto_selection_method(solver_obj):
    """The `family.method` method name AUTO was pinned to, or ''.

    SolverAUTO keeps it on `options.selection_method` rather than on
    `options.method`: `method` is what it hands the family solver it builds, so
    it is reset per delegation and reads 'default' by the time a getter runs.
    Reading the wrong one is not a cosmetic mistake -- it labels every swept
    solve 'AUTO' and the sweep produces no qualified table at all.
    """
    options = getattr(solver_obj, 'options', None)
    if options is None:
        return ''
    for attr in ('selection_method', 'selectionMethod', 'method'):
        value = getattr(options, attr, None)
        if isinstance(value, str) and value:
            return value
    return ''


LABEL_VALUES = frozenset(SOLVER_LABELS.values())


# The two label columns of each view, in the order the goldens carry them
# (Station/JobClass). A view whose table names its columns differently is read
# off the DataFrame itself; this map only says which columns are LABELS rather
# than metrics, so an unlisted view still records correctly.
VIEW_LABELS = {
    'avg': ('Station', 'JobClass'),
    'node': ('Node', 'JobClass'),
    'chain': ('Station', 'Chain'),
    'nodechain': ('Node', 'Chain'),
    'sys': ('Chain', 'JobClass'),
    'cache': ('Cache', 'Item'),
    'item': ('Item', 'JobClass'),
    'loss': ('Station', 'JobClass'),
    'regionloss': ('Region', 'JobClass'),
    'region': ('Region', 'JobClass'),
    'orbit': ('Station', 'JobClass'),
}

# Columns that are never metrics. Anything else numeric in a recorded table is
# kept, so a solver that reports a metric the goldens do not carry is recorded
# rather than dropped -- the comparator iterates the GOLDEN, so an extra column
# costs nothing and a missing one cannot be recovered later.
LABEL_COLUMNS = ('Station', 'Node', 'JobClass', 'Chain', 'Item', 'Cache',
                 'Region', 'Class', 'Name', 'NodeType', 'Type')


class Record(object):
    """One table, as the getter returned it, with what produced it."""

    __slots__ = ('solver', 'method', 'view', 'labels', 'rows', 'seq', 'derived')

    def __init__(self, solver, method, view, labels, rows, seq, derived=False):
        self.solver = solver
        self.method = method
        self.view = view
        self.labels = labels      # the two label column names, as found
        self.rows = rows          # [{col: value}], values are float or str
        self.seq = seq            # call order, 0-based
        # A DERIVED value is compared at the golden's own written precision
        # rather than at the shared five printed digits: these goldens were
        # written by the example's own format string, not by a result table.
        self.derived = derived

    def as_dict(self):
        return {'solver': self.solver, 'method': self.method, 'view': self.view,
                'labels': list(self.labels), 'rows': self.rows, 'seq': self.seq,
                'derived': self.derived}

    def __repr__(self):
        return '<Record %s(%s) %s %d rows>' % (self.solver, self.method,
                                               self.view, len(self.rows))


class ResultRecorder(object):
    """The buffer, and the on/off flag. One per process; see `recorder`."""

    def __init__(self):
        self.enabled = False
        self.records = []
        self.notes = []
        self._depth = 0
        self._seq = 0

    def enable(self):
        install()
        self.enabled = True

    def disable(self):
        self.enabled = False

    def reset(self):
        self.records = []
        self.notes = []
        self._depth = 0
        self._seq = 0

    def note(self, text):
        """Record something the run reported that is not a table.

        Used for a solver's own refusal ('this engine does not implement X'),
        which a consumer must be able to tell apart from a table that simply
        never arrived: the first is a fact about the port and is a named skip,
        the second is a failure.
        """
        self.notes.append(text)

    def capture(self, solver_obj, table, view):
        label = solver_label(solver_obj)
        if label is None:
            return
        rows, labels = table_rows(table, view)
        if not rows:
            return
        self.records.append(Record(label, solver_method(solver_obj), view,
                                   labels, rows, self._seq))
        self._seq += 1

    def capture_scalar(self, key, quantity, value, method='default'):
        """Record a derived quantity as a one-row table.

        `key` is the solver label ('CTMC') or the object's own key ('WF'), and
        `quantity` names what was computed ('probSysAggr', 'mean'). The value is
        flattened to a list of floats, because several of these getters return a
        pair or an array and the caller downstream decides which element the
        golden means.
        """
        values = _as_floats(value)
        if not values:
            return
        rows = [{'Quantity': quantity, 'Index': str(i), 'QLen': v}
                for i, v in enumerate(values)]
        self.records.append(Record(key, method, 'scalar',
                                   ('Quantity', 'Index'), rows, self._seq))
        self.records[-1].derived = True
        self._seq += 1

    def as_dict(self):
        return {'records': [r.as_dict() for r in self.records],
                'notes': list(self.notes)}


recorder = ResultRecorder()


def solver_label(solver_obj):
    """The golden's key for this solver, or None when it is not one we key by.

    A LAYERED or ENVIRONMENT solve is qualified by the member solver it drove --
    'LN(NC)', 'ENV(FLD)' -- because that is how several goldens spell it and
    because the two are genuinely different computations. The comparator still
    reconciles the qualified and bare spellings against the golden's own key;
    recording the member is what gives it the evidence to do so safely.
    """
    name = type(solver_obj).__name__
    base = SOLVER_LABELS.get(name)
    if base is None:
        # A subclass of a known solver (the `LINE` alias of SolverAUTO, the
        # `ENV` alias of SolverENV) keys as its parent.
        for klass in type(solver_obj).__mro__[1:]:
            base = SOLVER_LABELS.get(klass.__name__)
            if base is not None:
                break
    if base is None:
        return None
    # AUTO ASKED FOR ONE PAIR IS THAT PAIR. `SolverAUTO(model, 'mva.schmidt')`
    # delegates to exactly the family and method the method name names, so recording
    # it under 'AUTO' would key it by the meta-solver rather than by what ran --
    # and 'AUTO' is never goldened (it reports under whichever family it PICKED,
    # so its table would be scored against another family's numbers). The
    # qualified label is what a swept golden is keyed by; see qualified_label.
    if base == 'AUTO':
        return qualified_label(auto_selection_method(solver_obj)) or base
    if base not in ENSEMBLE_SOLVERS:
        return base
    # AN ENSEMBLE'S OWN METHOD NAMES IT WHEN IT HAS ONE, because that is the
    # distinction the goldens draw. `lqn_moment3` solves the SAME model with the
    # same NC layers twice, default and `moment3`, and prints the two under
    # `LN Results:` and `LN(moment3) Results:` -- its golden holds the first, and
    # they differ by 360x on T1. Labelling both `LN(NC)` would let the second
    # overwrite the first and report the default solve as a 99.7% error.
    method = solver_method(solver_obj)
    if method and method != 'default' and method.upper() not in LABEL_VALUES:
        return '%s(%s)' % (base, method)
    member = _member_solver(solver_obj)
    return '%s(%s)' % (base, member) if member else base


def _member_solver(solver_obj):
    """The layer/stage solver an ensemble solver actually ran, as a golden label.

    THE FACTORY IS THE DECLARATION, and it hides the class: `LN(model, lambda m:
    NC(m, opts))` names NC nowhere an attribute can be read. SolverLN already
    solves that for its own dispatcher with `probe_layer_solver_name`, which
    applies the factory to a layer and names the product, so that is what is
    asked here. Getting it wrong is not a spelling difference -- MVA layers and
    NC layers are different fixed points (lcq_threehosts: cache hit 0.5 against
    0.48331).
    """
    probe = getattr(solver_obj, 'probe_layer_solver_name', None)
    if callable(probe):
        try:
            name = probe()
        except Exception:
            name = None
        if name:
            # The probe answers with the class name on one path ('SolverNC') and
            # with the bare label on another ('NC'), so both spellings resolve.
            name = str(name)
            label = SOLVER_LABELS.get(name) or SOLVER_LABELS.get('Solver' + name)
            if label:
                return label
            if name.upper() in set(SOLVER_LABELS.values()):
                return name.upper()
    # SolverENV builds one stage solver per environment stage in its
    # constructor, so the stage engine is readable directly. It matters for the
    # same reason the layer solver does: the goldens key an environment solve by
    # the stage solver the example ran ('FLD' on renv_node_breakdown), because
    # that is the only solver name anywhere near the table.
    # SolverUQ is handed the solver CLASS it runs over each alternative model.
    klass = getattr(solver_obj, 'solver_class', None)
    if klass is not None:
        label = SOLVER_LABELS.get(getattr(klass, '__name__', ''))
        if label:
            return label
    for stage in (getattr(solver_obj, '_solvers', None) or ()):
        if stage is None:
            continue
        label = SOLVER_LABELS.get(type(stage).__name__)
        if label is None:
            for klass in type(stage).__mro__[1:]:
                label = SOLVER_LABELS.get(klass.__name__)
                if label is not None:
                    break
        if label:
            return label
    options = getattr(solver_obj, 'options', None)
    if options is None:
        return ''
    for attr in ('method', 'solver', 'layer_solver', 'stage_solver'):
        val = getattr(options, attr, None)
        if isinstance(val, str) and val:
            method_name = val.split('.')[-1].strip().upper()
            if method_name in ('MVA', 'NC', 'COMOM', 'CTMC', 'FLUID', 'FLD', 'SSA',
                         'LQNS', 'MAM'):
                return 'FLD' if method_name == 'FLUID' else method_name
    return ''


def solver_method(solver_obj):
    """The method this solver resolved, as a plain string.

    `MAM(dec.source)` and `MAM(inap)` differ by 236% on the same model, so the
    method is not decoration: it is what makes a recorded table attributable to
    the run the golden was generated from.
    """
    options = getattr(solver_obj, 'options', None)
    method = getattr(options, 'method', None) if options is not None else None
    if isinstance(method, str) and method:
        return method
    # A pinned AUTO keeps its `family.method` method name elsewhere (see
    # auto_selection_method), and that method name IS the method that ran. Consulted
    # only when `method` is unset, so nothing that already had one changes.
    method_name = auto_selection_method(solver_obj)
    return method_name if method_name else 'default'


def _as_floats(value):
    """Every float in a returned scalar, in order; [] when there is none.

    These getters return a float, a (value, error) pair, or an array, depending
    on the solver. Flattening rather than insisting on one shape is what lets a
    single record carry all three; the consumer names the element it wants.
    """
    if value is None or isinstance(value, bool):
        return []
    if isinstance(value, float) or isinstance(value, int):
        return [float(value)]
    try:
        import numpy as np
        if isinstance(value, np.ndarray):
            return [float(v) for v in value.ravel()]
    except ImportError:
        pass
    if isinstance(value, (list, tuple)):
        out = []
        for item in value:
            out.extend(_as_floats(item))
        return out
    try:
        return [float(value)]
    except (TypeError, ValueError):
        return []


def _dataframe(table):
    """The pandas frame behind a returned table, or None if it is not one.

    Handles the IndexedTable wrapper (which holds it on `.data`), a bare
    DataFrame, and anything else by declining.
    """
    data = getattr(table, 'data', None)
    if data is not None and hasattr(data, 'columns'):
        return data
    if hasattr(table, 'columns') and hasattr(table, 'itertuples'):
        return table
    return None


def table_rows(table, view):
    """(rows, label_columns) for a returned table; ([], ()) when it is not one.

    Values come out as Python floats and strings, at full precision. A cell that
    will not convert stays a string, which is right for a label column and
    harmless for anything else: the comparator only reads the metrics its golden
    names.
    """
    df = _dataframe(table)
    if df is None or len(df) == 0:
        return [], ()
    cols = [str(c) for c in df.columns]
    labels = tuple(c for c in cols if c in LABEL_COLUMNS)
    if not labels:
        # No named label column: fall back to the view's declared pair, keeping
        # only the ones this table actually has.
        labels = tuple(c for c in VIEW_LABELS.get(view, ()) if c in cols)
    rows = []
    for _, row in df.iterrows():
        out = {}
        for col in cols:
            val = row[col]
            if col in labels:
                out[col] = str(val)
                continue
            try:
                out[col] = float(val)
            except (TypeError, ValueError):
                out[col] = str(val)
        rows.append(out)
    return rows, labels


# ---------------------------------------------------------------------------
# Derived scalars.
#
# EIGHTEEN GOLDENS HOLD A QUANTITY NO RESULT TABLE CARRIES: a state probability,
# a workflow's phase-type moments, a cache hit ratio. They are computed by a
# library call the example makes and then printed, so the same construction that
# records a table records them -- the value, the call that produced it, and the
# solver it belongs to.
#
# They are recorded as one-row tables so that everything downstream handles a
# single shape. `parity/derived.py` is what maps them onto the key a particular
# golden uses, because THAT part is example-specific: `statepr_aggr`'s golden
# holds one probability where the example computes three, an artefact of the
# scraper having taken the first bare number it saw.
# ---------------------------------------------------------------------------

# Scalar getters worth recording, and the quantity each returns. A getter absent
# from here is not recorded: unlike the table getters, most methods on a solver
# return something that is not a result at all.
SCALAR_GETTERS = {
    'getProb': 'prob',
    'getProbAggr': 'probAggr',
    'getProbSys': 'probSys',
    'getProbSysAggr': 'probSysAggr',
    'getProbMarg': 'probMarg',
    'getProbSysMarg': 'probSysMarg',
    'getProbNormConstAggr': 'probNormConstAggr',
    'getProbNormConst': 'probNormConst',
    'getHitRatio': 'hitRatio',
    'getMissRatio': 'missRatio',
}
# Both snake_case spellings the package carries: `get_prob_sys_aggr` and the
# bare `prob_sys_aggr`. They are separate function objects on several classes,
# not aliases, so each must be named.
for _name, _q in list(SCALAR_GETTERS.items()):
    _snake = _camel_to_snake(_name)
    SCALAR_GETTERS[_snake] = _q
    if _snake.startswith('get_'):
        SCALAR_GETTERS[_snake[4:]] = _q

# Getters on a non-solver object, recorded under a key of their own because the
# quantity is the OBJECT's and no solver produced it. A Workflow's phase-type
# moments are the whole content of the seven `wf_*` goldens.
OBJECT_SCALAR_GETTERS = {
    ('Workflow', 'getMean'): ('WF', 'mean'),
    ('Workflow', 'getSCV'): ('WF', 'SCV'),
    ('Workflow', 'toPH'): ('WF', 'phases'),
}


# ---------------------------------------------------------------------------
# Installation.
# ---------------------------------------------------------------------------

# Method name -> the view it returns. Both spellings of every getter are here
# because the package carries both, often as separate function objects rather
# than aliases of one.
_VIEW_BY_METHOD = {
    'avg': ('getAvgTable', 'avgTable', 'avg_table', 'get_avg_table'),
    'node': ('getAvgNodeTable', 'avg_node_table', 'get_avg_node_table'),
    'chain': ('getAvgChainTable', 'avg_chain_table', 'get_avg_chain_table'),
    'nodechain': ('getAvgNodeChainTable', 'avg_node_chain_table',
                  'get_avg_node_chain_table'),
    'sys': ('getAvgSysTable', 'avg_sys_table', 'get_avg_sys_table'),
    'cache': ('getAvgCacheTable', 'avg_cache_table', 'get_avg_cache_table'),
    'item': ('getAvgItemTable', 'avg_item_table', 'get_avg_item_table'),
    'loss': ('getAvgLossTable', 'avg_loss_table', 'get_avg_loss_table'),
    'region': ('getAvgRegionTable', 'avg_region_table', 'get_avg_region_table'),
    'regionloss': ('getAvgRegionLossTable', 'avg_region_loss_table',
                   'get_avg_region_loss_table'),
    'orbit': ('getAvgOrbitTable', 'avg_orbit_table', 'get_avg_orbit_table'),
    # Single-metric views. They are recorded like any other: they merge into the
    # same solver's table under the same keys, so an example that asks only for
    # throughputs still contributes the column it asked for instead of nothing.
    'qlen': ('getAvgQLenTable', 'avg_qlen_table', 'get_avg_qlen_table'),
    'util': ('getAvgUtilTable', 'avg_util_table', 'get_avg_util_table'),
    'respt': ('getAvgRespTTable', 'avg_respt_table', 'get_avg_respt_table'),
    'tput': ('getAvgTputTable', 'avg_tput_table', 'get_avg_tput_table'),
}

METHOD_VIEW = dict((name, view)
                   for view, names in _VIEW_BY_METHOD.items() for name in names)

# Anything matching this and NOT in METHOD_VIEW is a getter nobody has classified
# yet. `install()` records it under a view named after the method, so it is
# captured rather than dropped, and the coverage test names it so the view can
# be declared deliberately.
GETTER_PATTERN = re.compile(
    r'^(getAvg\w*Table|get_avg\w*table|avg\w*_table|avgTable)$')

_installed = False


def solver_classes():
    """Every imported Solver class, base first.

    Walks `__subclasses__` from the roots rather than the package tree: a class
    that has not been imported cannot have produced a table, and walking the
    tree would import optional wrappers (and their third-party dependencies)
    that this process never asked for.
    """
    roots = []
    try:
        from line_solver.solvers.base import Solver, NetworkSolver
        roots += [Solver, NetworkSolver]
    except ImportError:
        pass
    for modname, clsname in (('line_solver.solvers.base', 'EnsembleSolver'),
                             ('line_solver.environment', 'SolverENV')):
        try:
            mod = __import__(modname, fromlist=[clsname])
        except ImportError:
            continue
        klass = getattr(mod, clsname, None)
        if klass is not None:
            roots.append(klass)

    seen = []
    stack = list(roots)
    while stack:
        klass = stack.pop()
        if any(klass is k for k in seen):
            continue
        seen.append(klass)
        try:
            stack.extend(klass.__subclasses__())
        except TypeError:
            pass
    return seen


def _wrap_scalar(fn, quantity):
    """Record the scalar this getter returns, attributed to its solver.

    TOP-LEVEL CALLS ONLY, on the same rule as the tables: a solve calls these
    getters internally (a layered fixed point evaluates a per-layer probability
    on every pass), and recording those would bury the one value the example
    asked for under hundreds it never saw.
    """
    @functools.wraps(fn)
    def wrapper(self, *args, **kwargs):
        if not recorder.enabled:
            return fn(self, *args, **kwargs)
        top = recorder._depth == 0
        recorder._depth += 1
        try:
            out = fn(self, *args, **kwargs)
        finally:
            recorder._depth -= 1
        if top:
            label = solver_label(self)
            if label is not None:
                recorder.capture_scalar(label, quantity, out,
                                        solver_method(self))
        return out
    wrapper.__line_recorded__ = quantity
    wrapper.__line_original__ = fn
    return wrapper


def _wrap_object_scalar(fn, key, quantity):
    """Record a quantity an ordinary model object computes (Workflow's PH).

    `toPH` returns (alpha, T) and the golden asks for the PHASE COUNT, so the
    number of phases is what is recorded -- taking the first float of a 3x3
    generator would record a rate.

    TOP-LEVEL CALLS ONLY. A layered solve composes an activity graph into a
    phase-type distribution for every entry of every layer, so `lqn_twotasks`
    alone calls `Workflow.toPH` 159 times inside its solve. Those are the
    solver's working, not a result the example asked for.
    """
    @functools.wraps(fn)
    def wrapper(self, *args, **kwargs):
        if not recorder.enabled:
            return fn(self, *args, **kwargs)
        top = recorder._depth == 0
        recorder._depth += 1
        try:
            out = fn(self, *args, **kwargs)
        finally:
            recorder._depth -= 1
        if top:
            value = out
            if quantity == 'phases':
                value = _phase_count(out)
            recorder.capture_scalar(key, quantity, value)
        return out
    wrapper.__line_recorded__ = quantity
    wrapper.__line_original__ = fn
    return wrapper


def _phase_count(ph):
    """The order of a phase-type representation returned as (alpha, T)."""
    try:
        _alpha, sub = ph
        return float(sub.shape[0])
    except (TypeError, ValueError, AttributeError, IndexError):
        return None


def _wrap(fn, view):
    @functools.wraps(fn)
    def wrapper(self, *args, **kwargs):
        if not recorder.enabled:
            return fn(self, *args, **kwargs)
        # OUTERMOST CALL ONLY. A facade getter calls through to a getter of the
        # same name on the solver it wraps, and several snake_case spellings are
        # thin forwards to the camelCase one, so recording at every level would
        # enter the same table two or three times under the same key.
        recorder._depth += 1
        try:
            out = fn(self, *args, **kwargs)
        finally:
            recorder._depth -= 1
        if recorder._depth == 0 and out is not None:
            recorder.capture(self, out, view)
        return out
    wrapper.__line_recorded__ = view
    wrapper.__line_original__ = fn
    return wrapper


def install():
    """Wrap every result-table getter on every imported Solver class.

    Idempotent, and safe to call again after new solver modules are imported --
    a second call picks up whatever the first could not see.

    ONE WRAPPER PER UNDERLYING FUNCTION, rebound everywhere it appears. Solver
    classes rebind each other's getters freely (`getAvgSysTable =
    NetworkSolver.getAvgSysTable`, `avg_sys_table = getAvgSysTable`), so wrapping
    per (class, name) would wrap one copy and leave the others calling the
    original -- silently unrecorded.
    """
    global _installed
    _installed = True
    classes = solver_classes()
    wrapped = {}
    for klass in classes:
        for name, value in list(vars(klass).items()):
            if not GETTER_PATTERN.match(name) or not callable(value):
                continue
            if getattr(value, '__line_recorded__', None) is not None:
                continue
            key = id(value)
            if key not in wrapped:
                wrapped[key] = (value, _wrap(value, METHOD_VIEW.get(name, name)))
            setattr(klass, name, wrapped[key][1])
    # SCALAR getters: a derived quantity no result table carries. Same
    # one-wrapper-per-function rule, for the same reason.
    for klass in classes:
        for name, value in list(vars(klass).items()):
            quantity = SCALAR_GETTERS.get(name)
            if quantity is None or not callable(value):
                continue
            if getattr(value, '__line_recorded__', None) is not None:
                continue
            key = id(value)
            if key not in wrapped:
                wrapped[key] = (value, _wrap_scalar(value, quantity))
            setattr(klass, name, wrapped[key][1])

    # Quantities an ordinary model object computes. There is one today (a
    # Workflow's phase-type moments, which are the whole content of the seven
    # `wf_*` goldens) and the table is keyed by class name so a second costs one
    # line.
    for (clsname, method), (key, quantity) in OBJECT_SCALAR_GETTERS.items():
        klass = _model_class(clsname)
        if klass is None:
            continue
        fn = klass.__dict__.get(method)
        if fn is None or getattr(fn, '__line_recorded__', None) is not None:
            continue
        setattr(klass, method, _wrap_object_scalar(fn, key, quantity))

    # A getter rebound onto ANOTHER class under a name this sweep did not visit
    # (an alias with a different spelling) still points at the original; rebind
    # those too, so the two spellings cannot disagree about whether they record.
    originals = dict((id(orig), wrap) for orig, wrap in wrapped.values())
    for klass in classes:
        for name, value in list(vars(klass).items()):
            wrap = originals.get(id(value))
            if wrap is not None and value is not wrap:
                setattr(klass, name, wrap)


def _model_class(clsname):
    """A model class by name, or None when this build has not imported it."""
    import importlib
    for modname in ('line_solver.lang.workflow', 'line_solver.lang'):
        try:
            mod = importlib.import_module(modname)
        except ImportError:
            continue
        klass = getattr(mod, clsname, None)
        if klass is not None:
            return klass
    return None


if os.environ.get('LINE_RECORD_RESULTS'):
    recorder.enable()
