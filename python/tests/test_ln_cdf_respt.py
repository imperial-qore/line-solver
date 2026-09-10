"""
SolverLN.getCdfRespT: the per-entry response time distribution.

The law exists only under the moment3 update pass, which fits an APH to the
first three moments of every activity and call term and convolves them. The
getter therefore re-runs the ensemble under moment3 when the solver was built
for another method, exactly as MATLAB @SolverLN/getCdfRespT.m does, and returns
one (n, 2) [F(t), t] table per entry in the ENTRY-LOCAL index space.

The model is the one of examples/basic/layeredModel/lqn_moment3: a reference
task T1 whose single activity A1 (Exp(1)) calls E2 and E3 on an infinite task
T2, with E2 served by the serial chain A20 -> A21 -> A22 and E3 by A3.

Entry 0 (E1) is the case that used to be lost: both ports carried MATLAB's
1-based `0 < e <= nentries` guard into a 0-based index space, so E1's law was
never stored and the getter had nothing to return.
"""

import numpy as np
import pytest

from line_solver import (LayeredNetwork, Processor, Task, Entry, Activity,
                         ActivityPrecedence, Erlang, Exp, SchedStrategy,
                         SolverLN, SolverMVA)
from line_solver.solvers.solver_ln.solver_ln import SolverLNOptions


def _build():
    m = LayeredNetwork('myLayeredModel')
    P1 = Processor(m, 'P1', 1, SchedStrategy.PS)
    T1 = Task(m, 'T1', 100, SchedStrategy.REF).on(P1)
    E1 = Entry(m, 'E1').on(T1)
    P2 = Processor(m, 'P2', 1, SchedStrategy.PS)
    T2 = Task(m, 'T2', 1, SchedStrategy.INF).on(P2)
    E2 = Entry(m, 'E2').on(T2)
    E3 = Entry(m, 'E3').on(T2)
    T1.set_think_time(Erlang.fit_mean_and_order(10, 1))
    Activity(m, 'A1', Exp(1)).on(T1).bound_to(E1).synch_call(E2).synch_call(E3, 1)
    A20 = Activity(m, 'A20', Exp(1)).on(T2).bound_to(E2)
    A21 = Activity(m, 'A21', Exp(1)).on(T2)
    A22 = Activity(m, 'A22', Exp(1)).on(T2).replies_to(E2)
    T2.add_precedence(ActivityPrecedence.serial([A20, A21, A22]))
    Activity(m, 'A3', Exp(1)).on(T2).bound_to(E3).replies_to(E3)
    return m


def _solver(method=None):
    """MVA layers, as the MATLAB reference example uses.

    NOT SolverNC: its getCdfRespT is the tagged-job STDF law, which is defined
    for FCFS stations only and returns [] with a warning on the PS/INF layers an
    LQN is mostly made of. The repo is then legitimately empty and every entry
    law collapses -- an upstream limit of that solver, not of the getter.
    """
    opts = SolverLNOptions()
    opts.verbose = False
    if method is not None:
        opts.method = method
    return SolverLN(_build(), lambda m: SolverMVA(m, verbose=False), opts)


def _assert_is_cdf_table(cdf):
    """[F(t), t], NOT [t, F(t)]: column 0 is the probability, column 1 the time."""
    assert isinstance(cdf, np.ndarray)
    assert cdf.ndim == 2 and cdf.shape[1] == 2
    F, t = cdf[:, 0], cdf[:, 1]
    assert np.all(np.diff(t) >= 0), "the time column must be nondecreasing"
    assert t[0] == pytest.approx(0.0)
    assert np.all(F >= -1e-12) and np.all(F <= 1 + 1e-12)
    assert np.all(np.diff(F) >= -1e-9), "the CDF column must be nondecreasing"
    assert F[0] == pytest.approx(0.0, abs=1e-9)
    assert F[-1] > 0.99, "the grid must reach into the tail"


def test_moment3_returns_one_law_per_entry():
    s = _solver('moment3')
    s.get_avg_table()
    RD = s.getCdfRespT()
    assert len(RD) == s.lqn.nentries == 3
    # EVERY entry, not just the first. E1 is the one the off-by-one used to
    # drop; E2 and E3 are the ones the layer-CDF contract mismatch starved.
    for e in range(len(RD)):
        _assert_is_cdf_table(RD[e])


def test_layer_cdf_repo_accepts_both_getter_contracts():
    """SolverFLD returns RD[station][class] arrays; SolverMVA/SolverNC return a
    flat list of {station, class, t, p} dicts. Both must reach the entry
    assembly -- the guard this replaced dropped every dict in silence."""
    s = _solver('moment3')
    s.get_avg_table()

    class _Layer(object):
        class _Sn(object):
            nstations, nclasses = 2, 2

        def getStruct(self):
            return self._Sn()

    t = np.array([0.0, 1.0, 2.0])
    p = np.array([0.0, 0.5, 1.0])
    flat = [{'station': 2, 'class': 1, 't': t, 'p': p}]
    RD = s._nested_respt_cdf(flat, _Layer())
    assert RD[0][0] is None and RD[0][1] is None and RD[1][1] is None
    cell = RD[1][0]                       # 1-based on the wire -> 0-based here
    assert cell.shape == (3, 2)
    assert np.allclose(cell[:, 0], p)     # column 0 is the CDF
    assert np.allclose(cell[:, 1], t)     # column 1 is the time

    nested = [[np.column_stack([p, t]), None], [None, None]]
    assert s._nested_respt_cdf(nested, _Layer()) is nested

    # A shape nobody publishes is an error, not a silent drop.
    with pytest.raises(TypeError, match="unrecognised getCdfRespT return shape"):
        s._nested_respt_cdf(["not a contract"], _Layer())
    with pytest.raises(ValueError, match="outside the layer"):
        s._nested_respt_cdf([{'station': 9, 'class': 1, 't': t, 'p': p}], _Layer())
    with pytest.raises(ValueError, match="times and"):
        s._nested_respt_cdf([{'station': 1, 'class': 1, 't': t, 'p': p[:2]}], _Layer())


def test_caller_entry_convolves_its_call_terms():
    """E1's law must carry A1's own service AND both synchronous calls.

    When the call terms are missing the entry reports its activity's bare host
    demand (~Exp(1), mean 1.07) while the activity row beside it carries the
    full 390 -- which is the tell that the layer CDF repo is starved.
    """
    s = _solver('moment3')
    table = s.get_avg_table()
    names = [str(n) for n in np.asarray(table['Node'])]
    respt = [float(v) for v in np.asarray(table['RespT'])]
    e1, a1 = respt[names.index('E1')], respt[names.index('A1')]
    assert e1 > 100.0, "E1 collapsed to its bare host demand: call terms lost"
    assert e1 == pytest.approx(a1, rel=0.05)


def test_law_agrees_with_the_reported_entry_response_time():
    """The tabulated law and the entry's mean are the SAME quantity.

    They part company when the interlock rescale of update_populations is let
    at the entry servt after the moment pass has formed it -- BUG-97, whose
    tell is a service time BELOW that of the activity the entry holds.
    """
    s = _solver('moment3')
    table = s.get_avg_table()
    RD = s.getCdfRespT()
    F, t = RD[0][:, 0], RD[0][:, 1]
    # E[X] = int (1-F) dt over the tabulated grid.
    mean_from_cdf = float(np.trapezoid(1.0 - F, t)) if hasattr(np, 'trapezoid') \
        else float(np.trapz(1.0 - F, t))
    assert mean_from_cdf == pytest.approx(s.entryproc[0].getMean(), rel=2e-2)

    names = [str(n) for n in np.asarray(table['Node'])]
    respt = [float(v) for v in np.asarray(table['RespT'])]
    # The AvgTable CARRIES 4 significant figures, not just displays them (388.8
    # against the law's 388.804925), so the tolerance is set by that rounding --
    # measured 1.3e-5, 8.4e-6 and 3.4e-4 on the three entries. It stays orders of
    # magnitude below the FACTOR a visit-ratio error would introduce, which is
    # what this assertion is here to catch.
    assert respt[names.index('E1')] == pytest.approx(s.entryproc[0].getMean(), rel=1e-3)
    # A1 holds Exp(1) plus what it waits on; the entry cannot be faster.
    assert s.entryproc[0].getMean() > 1.0


def test_a_second_solve_still_reports_metrics():
    """The moment3 pass is terminal WITHIN ONE SOLVE only.

    The flag that makes it terminal is scoped to one iterate(); left standing it
    short-circuits converged() at it=0 on the next solve, the loop body never
    runs, and every metric comes back zero.
    """
    s = _solver('moment3')
    first = s.get_avg_table()
    s.getCdfRespT()
    second = s.get_avg_table()
    for table in (first, second):
        names = [str(n) for n in np.asarray(table['Node'])]
        respt = [float(v) for v in np.asarray(table['RespT'])]
        assert respt[names.index('E1')] > 1.0


def test_getter_reruns_under_moment3_and_restores_the_method():
    """The mean-based update forms no law, so the getter flips the method itself."""
    s = _solver('srvn.cs')
    s.get_avg_table()
    assert s.lnmethod == 'srvn.cs'
    RD = s.getCdfRespT()
    _assert_is_cdf_table(RD[0])
    assert s.lnmethod == 'srvn.cs', "the caller's method must be put back"
    assert getattr(s.options, 'method') == 'srvn.cs'


def test_ph_encoding_is_refused_by_name():
    """srvn.ph replaces the activity graph by a composed server law, so there is
    no routing to re-run the distribution pass over. Refuse rather than rebuild
    the wrong topology."""
    s = _solver('srvn.ph')
    s.get_avg_table()
    if not s._is_ph_encoding():
        pytest.skip("the alias fell back to a routing encoding on this model")
    with pytest.raises(ValueError, match="routing encoding"):
        s.getCdfRespT()


def test_aph_evalcdf_call_forms():
    """APH.evalCDF mirrors MATLAB's three forms; the no-arg one is what the
    moment pass calls, and it used to raise TypeError into a bare except."""
    s = _solver('moment3')
    s.get_avg_table()
    aph = s.entryproc[0]
    grid = aph.evalCDF()
    _assert_is_cdf_table(grid)
    assert grid.shape[0] == 500
    scalar = aph.evalCDF(float(grid[100, 1]))
    assert isinstance(scalar, float)
    assert scalar == pytest.approx(float(grid[100, 0]), abs=1e-9)
    vector = aph.evalCDF(grid[:10, 1])
    assert np.asarray(vector).shape == (10,)
    assert np.allclose(vector, grid[:10, 0])
