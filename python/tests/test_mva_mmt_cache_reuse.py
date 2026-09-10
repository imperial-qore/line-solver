"""Regression tests for the reuse of SolverMVA's cached MMT transformation.

SolverMVA caches the fork-join (MMT) transformation of a model so that an outer
fixed point -- SolverLN re-solving a layer once per iteration -- does not pay for
a full deep copy of the model every time. The transformed model freezes the
service rates it was built with, so reusing it is only correct if those rates are
re-fed from the base model first, and if the slots the transformation owns (the
joins it turns into zero-service delays, which the fork loop then overwrites with
synchronisation delays) are reset to their initial value.

Getting this wrong is silent: it produces plausible numbers that are simply
wrong. These tests pin the invariant that a reused solver agrees exactly with a
cold one.

The reuse tests below discard the fork warm start (options.config.fj_warmstart,
SolverMVA._fj_fork_lambda) before each re-solve. Like options.init_sol, that
iterate deliberately survives reset() in all three codebases, so leaving it in
place makes the MMT fixed point start from the previous converged point and land
on a different point of the same tolerance ball. That is the warm start working
as designed, not stale cache state, and it would mask the invariant under test.
"""
import os

import numpy as np
import pytest

from line_solver import (Network, Delay, Queue, Fork, Join, ClosedClass, Exp,
                         Erlang, SchedStrategy, MVA, GlobalConstants,
                         VerboseLevel)

GlobalConstants.set_verbose(VerboseLevel.SILENT)

METRICS = ['QLen', 'Util', 'Tput']


def _build_fj(rate2, erlang=False):
    """Closed fork-join network whose Q2 service is the parameter under test."""
    model = Network('fj')
    d = Delay(model, 'Clients')
    q1 = Queue(model, 'Q1', SchedStrategy.PS)
    q2 = Queue(model, 'Q2', SchedStrategy.PS)
    f = Fork(model, 'F')
    j = Join(model, 'J', f)
    c = ClosedClass(model, 'C', 4, d)
    d.set_service(c, Exp(1.0))
    q1.set_service(c, Exp(2.0))
    q2.set_service(c, Erlang.fit_mean_and_order(1.0 / rate2, 3) if erlang else Exp(rate2))
    P = model.init_routing_matrix()
    P.set(c, c, d, f, 1.0)
    P.set(c, c, f, q1, 1.0)
    P.set(c, c, f, q2, 1.0)
    P.set(c, c, q1, j, 1.0)
    P.set(c, c, q2, j, 1.0)
    P.set(c, c, j, d, 1.0)
    model.link(P)
    return model, q2, c


def _avg(solver):
    return np.asarray(solver.get_avg_table()[METRICS].values, dtype=float)


def test_reused_solver_matches_cold_solver_across_rate_changes():
    """The SolverLN pattern: reset() then re-solve after only the rates changed."""
    rates = [2.0, 3.0, 1.5, 5.0, 2.5]

    cold = []
    for r in rates:
        model, _, _ = _build_fj(r)
        cold.append(_avg(MVA(model)))

    model, q2, c = _build_fj(rates[0])
    solver = MVA(model)
    for i, r in enumerate(rates):
        q2.set_service(c, Exp(r))
        solver.reset()
        solver.reset_fork_warm_start()
        reused = _avg(solver)
        np.testing.assert_array_equal(
            reused, cold[i],
            err_msg=("reused MMT transformation disagrees with a cold solve at "
                     "rate2=%s; the cached nonfjmodel is holding stale state" % r))


def test_reused_solver_invalidates_on_structural_change():
    """A phase-count change is structural, so the cached transformation of the old
    topology must be discarded rather than re-fed."""
    model, q2, c = _build_fj(2.0)
    solver = MVA(model)
    solver.get_avg_table()

    q2.set_service(c, Erlang.fit_mean_and_order(1.0 / 3.0, 3))
    solver.reset()
    solver.reset_fork_warm_start()
    reused = _avg(solver)

    cold_model, _, _ = _build_fj(3.0, erlang=True)
    np.testing.assert_array_equal(
        reused, _avg(MVA(cold_model)),
        err_msg="cached MMT transformation survived a structural change")


@pytest.mark.skipif(
    os.environ.get('LINE_SOLVER_LANG') == 'java',
    reason="inspects SolverMVA._mmt_cache, a native-Python internal that the "
           "java-dispatch path (solve delegated to jline.jar) does not populate")
def test_join_service_is_reset_not_inherited():
    """The join delay is owned by the fork loop, which overwrites it with the
    current sync delay each pass. On reuse it must be restored to Immediate, not
    left holding the previous solve's converged value."""
    from line_solver.io.model_adapter import ModelAdapter

    model, q2, c = _build_fj(2.0)
    solver = MVA(model)
    solver.get_avg_table()
    cache = getattr(solver, '_mmt_cache', None)
    assert cache is not None, "expected the MMT transformation to be cached"

    mmt_result = cache['mmt_result']
    assert mmt_result.immediate_slots, "join slots were not recorded as mmt-owned"
    # No mmt-owned slot may claim a base-model provenance, or a refresh would
    # overwrite the transformation's own Immediate with a base service.
    assert not (set(mmt_result.immediate_slots) & set(mmt_result.service_src)), \
        "a slot is both mmt-owned and base-derived"

    assert ModelAdapter.refresh_services_from_base(mmt_result)
    nonfj = mmt_result.nonfjmodel
    for (i, k) in mmt_result.immediate_slots:
        svc = nonfj._nodes[i].get_service(nonfj._classes[k])
        assert svc is not None and svc.isImmediate(), (
            "join slot (%d,%d) was not restored to Immediate on refresh" % (i, k))
