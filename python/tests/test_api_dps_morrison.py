"""
Morrison's heavy-usage expansion for a closed think+DPS network, and its NC solver route.

The pinned values are the MATLAB reference (matlab/src/api/npfqn/npfqn_dps_morrison.m) and agree
with the JAR and C++ ports to the digits written here; the same models appear in the test of each
codebase, so a divergence in any of them fails somewhere.

Reference: J.A. Morrison, Queueing Systems 9 (1991) 191-214.
"""

import numpy as np
import pytest

from line_solver import (ClosedClass, Delay, Erlang, Exp, GlobalConstants, NC, Network, Queue,
                         SchedStrategy, VerboseLevel)
from line_solver.api.npfqn import npfqn_dps_morrison

GlobalConstants.set_verbose(VerboseLevel.SILENT)

TOL = 1e-9

# M2: two classes, K = 100 each, rho = 0.9 -- the regime the expansion is derived for.
M2 = dict(N=[100.0, 100.0], Z=[1.0, 0.5], S=[0.006, 0.0015], w=[1.0, 4.0])
# M3: three classes, unequal populations and weights, rho = 1.06 (appendix A's saturated side).
M3 = dict(N=[20.0, 12.0, 28.0], Z=[1.0, 0.5, 2.0], S=[0.02, 0.01, 0.03], w=[1.0, 4.0, 2.0])


def test_two_classes_matches_reference():
    res = npfqn_dps_morrison(**M2)
    assert res.rho == pytest.approx(0.9, abs=1e-12)
    assert res.a == pytest.approx(0.1, abs=1e-12)
    assert res.Q == pytest.approx([4.053964005705, 0.930378969124], abs=TOL)
    assert res.R == pytest.approx([0.042452089190, 0.004666835854], abs=TOL)
    assert res.X == pytest.approx([95.946035994295, 198.139242061752], abs=1e-8)
    assert res.Qlead == pytest.approx([4.373155763276, 0.546644470409], abs=TOL)
    assert res.sigma == pytest.approx([0.084297520661, -0.337190082645], abs=TOL)


def test_three_classes_matches_reference():
    res = npfqn_dps_morrison(**M3)
    assert res.rho == pytest.approx(1.06, abs=1e-12)
    assert res.Q == pytest.approx([3.213021210443, 0.918773941158, 2.645208279899], abs=TOL)
    assert res.R == pytest.approx([0.208463525546, 0.039776387080, 0.202390704352], abs=TOL)
    assert res.sigma == pytest.approx([0.497718870654, -0.218459289315, -0.258992817331], abs=TOL)


def test_littles_law_and_the_two_response_time_forms():
    """Throughput closes the think station EXACTLY, which is why the analyzers report R = Q/T.

    Morrison's RESULT 2 (4.17) is the EXPANDED ratio (4.11)/(4.15), so it differs from Q/X by
    higher-order terms -- 22% at K=6, 0.5% at K=100, 0.09% at K=400. The two forms agree only
    asymptotically, and that convergence is what is asserted here."""
    res = npfqn_dps_morrison(**M3)
    N = np.asarray(M3['N']); Z = np.asarray(M3['Z'])
    assert res.X == pytest.approx((N - res.Q) / Z, abs=1e-10)

    def gap(kernel):
        return np.max(np.abs(kernel.R - kernel.Q / kernel.X) / kernel.R)

    g6 = gap(npfqn_dps_morrison([6.0, 6.0], M2['Z'], [0.075, 0.0375], M2['w']))
    g100 = gap(npfqn_dps_morrison(**M2))
    g400 = gap(npfqn_dps_morrison([400.0, 400.0], M2['Z'], [0.0015, 0.000375], M2['w']))
    assert g400 < g100 < g6
    assert g100 < 0.01


def test_equal_weights_collapse_the_correction():
    """Equal weights make the network product-form, so Morrison's (4.18)-(4.21) collapse:
    D = B, Q = C, M = H = I, L = J, hence R = S = V = U = 0, A = -K, delta = 0 and sigma = 0."""
    res = npfqn_dps_morrison(M2['N'], M2['Z'], M2['S'], [1.0, 1.0])
    assert res.cD == pytest.approx(res.cB, rel=1e-9)
    assert res.cQ == pytest.approx(res.cC, rel=1e-9)
    assert res.cI == pytest.approx(res.cH, rel=1e-9)
    assert res.cM == pytest.approx(res.cH, rel=1e-9)
    assert res.cL == pytest.approx(res.cJ, rel=1e-9)
    assert res.cR == pytest.approx(0.0, abs=1e-9)
    assert res.cS == pytest.approx(0.0, abs=1e-9)
    assert res.cV == pytest.approx(0.0, abs=1e-9)
    assert res.cU == pytest.approx(0.0, abs=1e-9)
    assert res.cA == pytest.approx(-res.cK, rel=1e-9)
    assert res.delta == pytest.approx(0.0, abs=1e-9)
    assert res.sigma == pytest.approx([0.0, 0.0], abs=1e-9)
    assert res.Q == pytest.approx([3.737598812660, 1.937595955666], abs=TOL)


def test_accuracy_improves_with_population():
    """An asymptotic expansion has to get relatively better as the populations grow at fixed
    usage. Both models run at rho = 0.9; the exact values are from the CTMC of eq. (2.1)."""
    exact_small = np.array([1.178868777433, 0.746116560046])
    exact_big = np.array([4.082236085180, 0.862829749692])
    small = npfqn_dps_morrison([6.0, 6.0], M2['Z'], [0.075, 0.0375], M2['w'])
    big = npfqn_dps_morrison(**M2)
    e_small = np.abs(small.Q - exact_small) / exact_small
    e_big = np.abs(big.Q - exact_big) / exact_big
    assert np.all(e_big < e_small)
    assert np.all(e_big < 0.1)


def test_invalid_inputs_rejected():
    with pytest.raises(ValueError):
        npfqn_dps_morrison(M2['N'], M2['Z'], M2['S'], [1.0, -4.0])
    with pytest.raises(ValueError):
        npfqn_dps_morrison([0.0, 100.0], M2['Z'], M2['S'], M2['w'])
    with pytest.raises(ValueError):
        npfqn_dps_morrison(M2['N'], M2['Z'], M2['S'], [1.0])


def _think_dps(n1=100, n2=100, s1=0.006, s2=0.0015, w2=4.0, dist1=None):
    model = Network('morrison')
    delay = Delay(model, 'Think')
    queue = Queue(model, 'DPS', SchedStrategy.DPS)
    c1 = ClosedClass(model, 'C1', n1, delay, 0)
    c2 = ClosedClass(model, 'C2', n2, delay, 0)
    delay.set_service(c1, Exp(1.0))
    delay.set_service(c2, Exp(2.0))
    queue.set_service(c1, dist1 if dist1 is not None else Exp(1.0 / s1), 1.0)
    queue.set_service(c2, Exp(1.0 / s2), w2)
    P = model.init_routing_matrix()
    for c in (c1, c2):
        P.set(c, c, delay, queue, 1.0)
        P.set(c, c, queue, delay, 1.0)
    model.link(P)
    return model


def _qlen(table, station):
    df = table.to_dataframe() if hasattr(table, 'to_dataframe') else table
    sel = df[df['Station'].astype(str).str.strip() == station]
    return sel['QLen'].to_numpy(dtype=float)


def test_solver_nc_routes_to_morrison():
    """SolverNC answers the think+DPS shape by Morrison, conserving the population exactly."""
    solver = NC(_think_dps(), verbose=False)
    table = solver.avg_table()
    q_dps = _qlen(table, 'DPS')
    q_think = _qlen(table, 'Think')
    assert q_dps == pytest.approx([4.053964005705, 0.930378969124], abs=1e-6)
    assert (q_dps + q_think) == pytest.approx([100.0, 100.0], abs=1e-9)


def test_morrison_is_advertised_and_other_methods_refused():
    model = _think_dps()
    assert 'morrison' in NC(model, verbose=False).listValidMethods()
    with pytest.raises(Exception):
        NC(_think_dps(), method='exact', verbose=False).avg_table()


def test_dps_outside_shape_is_refused():
    """A DPS station outside Morrison's shape is refused rather than approximated."""
    model = _think_dps(dist1=Erlang.fit_mean_and_order(0.006, 3))
    with pytest.raises(Exception):
        NC(model, verbose=False).avg_table()


def _plain_ps():
    """A product-form PS model: no DPS station anywhere."""
    model = Network('ps')
    delay = Delay(model, 'Think')
    queue = Queue(model, 'PS', SchedStrategy.PS)
    c1 = ClosedClass(model, 'C1', 6, delay, 0)
    c2 = ClosedClass(model, 'C2', 6, delay, 0)
    delay.set_service(c1, Exp(1.0))
    delay.set_service(c2, Exp(2.0))
    queue.set_service(c1, Exp(3.0))
    queue.set_service(c2, Exp(5.0))
    P = model.init_routing_matrix()
    for c in (c1, c2):
        P.set(c, c, delay, queue, 1.0)
        P.set(c, c, queue, delay, 1.0)
    model.link(P)
    return model


def test_morrison_named_on_a_non_dps_model_is_refused():
    """'morrison' is in listValidMethods, so without an explicit arm it would pass the method
    gate, match no route, and fall through to the ordinary normalizing-constant path -- answering
    a product-form model UNDER THE CALLER'S LABEL. It must be refused by name instead."""
    with pytest.raises(ValueError, match='morrison'):
        NC(_plain_ps(), method='morrison', verbose=False).avg_table()


def test_analyzer_refuses_a_model_off_the_shape_when_called_directly():
    """The analyzer re-checks the shape rather than trusting its caller."""
    from line_solver.solvers.solver_nc.solver_nc_dps_analyzer import solver_nc_dps_analyzer
    model = _plain_ps()
    with pytest.raises(ValueError):
        solver_nc_dps_analyzer(model.getStruct(), NC(model, verbose=False).options)


def test_default_still_solves_the_ps_model_exactly():
    """The gate must not disturb the ordinary product-form path."""
    q = _qlen(NC(_plain_ps(), verbose=False).avg_table(), 'PS')
    assert np.all(np.isfinite(q)) and np.all(q > 0)
