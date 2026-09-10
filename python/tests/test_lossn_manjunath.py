"""
Validation of the exact Manjunath-Sikdar transform for loss networks
(lossn_manjunath) and of the exact path through solver_nc_lossn_analyzer.

The primary oracle is an independently written enumeration of the admissible
set: the transform's whole point is that it never enumerates, so a
coefficient-domain bug cannot hide behind a shared traversal. Single-link cases
are additionally pinned to the Erlang-B recursion, a closed form with no code in
common with either.

The end-to-end expectations are the MATLAB reference implementation's output for
the same model, reproduced to 16 digits because the transform is exact: there is
no tolerance to hide behind, and any drift is a defect rather than noise.
"""

import itertools
import math

import numpy as np
import pytest

from line_solver import (Delay, Exp, NC, Network, OpenClass, Sink, Source)
from line_solver.api.lossn import lossn_manjunath

# MATLAB SolverNC(fcr_lossn) with method 'default' -> 'default/lossn.exact'
LAMBDA1, LAMBDA2, MU1, MU2 = 0.3, 0.2, 1.0, 0.8
GMAX, C1MAX, C2MAX = 5, 3, 3
MATLAB_Q = (0.2989813602605746, 0.2494742965193094)
MATLAB_X = (0.2989813602605746, 0.1995794372154475)
MATLAB_LOSS = (0.003395465798084696, 0.002102813922762459)
MATLAB_LG = 0.5495940110776684


def _dispatched(result):
    """The method the analyzer actually ran.

    The reported name carries the REQUESTED method as a prefix when the two
    differ, exactly as MATLAB reports 'default/lossn.exact' for this model. The
    native python result records the resolved half alone, so the prefix is
    stripped here rather than asserted on: what these tests discriminate is
    which of the three loss-network paths answered, not how the request was
    spelled.
    """
    return str(getattr(result, 'method', None) or '').rsplit('/', 1)[-1]


def erlang_b(a, c):
    """Erlang B by the numerically stable recursion."""
    b = 1.0
    for k in range(1, int(c) + 1):
        b = a * b / (k + a * b)
    return b


def enumerate_exact(nu, A, C):
    """
    Independent oracle: g(C), E[n_r] and blocking by direct enumeration.

    A route in no constraint row has an untruncated Poisson marginal, so it never
    blocks, carries its full load, and factors exp(nu_r) out of g(C); summing a
    truncated box over it would be a different network.
    """
    nu = np.asarray(nu, float)
    A = np.atleast_2d(np.asarray(A, float))
    C = np.asarray(C, float)
    R = len(nu)
    free = ~np.any(A != 0, axis=0)
    N = []
    for r in range(R):
        pos = A[:, r] > 0
        N.append(int(np.floor(np.min(C[pos] / A[pos, r]))) if np.any(pos) else 0)
    G = 0.0
    En = np.zeros(R)
    acc = np.zeros(R)
    for n in itertools.product(*[range(v + 1) for v in N]):
        n = np.array(n)
        if np.any(A @ n > C + 1e-12):
            continue
        w = np.prod([nu[r] ** n[r] / math.factorial(n[r]) for r in range(R)])
        G += w
        En += w * n
        for r in range(R):
            if np.all(A @ (n + np.eye(R, dtype=int)[r]) <= C + 1e-12):
                acc[r] += w
    qlen = En / G
    loss = 1.0 - acc / G
    qlen[free] = nu[free]
    loss[free] = 0.0
    return qlen, loss, math.log(G) + float(np.sum(nu[free]))


CASES = [
    # single link, two classes with different circuit requirements
    ([1.5, 0.7], [[1.0, 2.0]], [5.0]),
    # the shape a FiniteCapacityRegion produces: a global cap plus per-class caps
    ([2.0, 1.0], [[1.0, 1.0], [1.0, 0.0], [0.0, 1.0]], [6.0, 4.0, 3.0]),
    # overlapping rows, so more than one link is live at once and the interleaved
    # elimination order actually does something
    ([1.0, 2.0, 0.5], [[1.0, 1.0, 0.0], [0.0, 1.0, 1.0], [1.0, 0.0, 1.0]], [4.0, 5.0, 3.0]),
    # (4, 8) <= 20 reduces to (1, 2) <= 5: the dimension shrinks from 21
    # coefficients to 6 and the answer must not move
    ([3.0, 1.0], [[4.0, 8.0]], [20.0]),
    # class 2 appears in no row: it never blocks and contributes exp(nu_2) to g(C)
    ([1.0, 2.5], [[1.0, 0.0]], [4.0]),
    # nu = 0 has no logarithm to shift by, hence its own branch
    ([0.0, 2.0], [[1.0, 1.0]], [3.0]),
    # a single class 2 call needs 5 units of a link with 3: blocked in every state
    ([1.0, 1.0], [[1.0, 5.0]], [3.0]),
]


@pytest.mark.parametrize("nu,A,C", CASES)
def test_matches_direct_enumeration(nu, A, C):
    q, loss, lG, niter = lossn_manjunath(nu, A, C)
    eq, el, elG = enumerate_exact(nu, A, C)
    np.testing.assert_allclose(q, eq, rtol=0, atol=1e-10)
    np.testing.assert_allclose(loss, el, rtol=0, atol=1e-10)
    assert abs(lG - elG) < 1e-10
    assert niter == 1, "the transform is direct"


def test_class_that_cannot_fit_is_fully_blocked():
    q, loss, _lG, _n = lossn_manjunath([1.0, 1.0], [[1.0, 5.0]], [3.0])
    assert loss[1] == pytest.approx(1.0, abs=1e-12)
    assert q[1] == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("a,c", [(1.0, 1), (5.0, 10), (30.0, 25), (200.0, 210)])
def test_single_link_is_erlang_b(a, c):
    q, loss, _lG, _n = lossn_manjunath([a], [[1.0]], [c])
    b = erlang_b(a, c)
    assert loss[0] == pytest.approx(b, abs=1e-11)
    assert q[0] == pytest.approx(a * (1 - b), abs=1e-9)


def test_heavy_load_does_not_overflow():
    """
    nu = C = 900 is where the pre-fix reference returned nan: nu**n/n! peaks at
    exp(nu)/sqrt(2 pi nu), so forming the terms before rescaling overflows to inf,
    and a rescaling step guarded on finiteness then declines to run, sending the
    inf through to g.
    """
    q, loss, lG, _n = lossn_manjunath([900.0], [[1.0]], [900.0])
    assert np.isfinite(loss[0]) and np.isfinite(q[0]) and np.isfinite(lG)
    b = erlang_b(900.0, 900)
    assert loss[0] == pytest.approx(b, abs=1e-12)
    assert q[0] == pytest.approx(900.0 * (1 - b), abs=1e-9)


@pytest.mark.parametrize("A,C", [([[1.5]], [3.0]), ([[1.0]], [2.5])])
def test_fractional_input_refused_by_name(A, C):
    # The residue argument counts whole units of capacity, so a fractional entry
    # is refused rather than rounded into a different network.
    with pytest.raises(ValueError, match="lossn_mci"):
        lossn_manjunath([1.0], A, C)


def test_live_state_cap_refused():
    # Peak memory is the product of (C_j+1) over the simultaneously live links, so
    # an oversized region is refused rather than allowed to exhaust memory.
    with pytest.raises(RuntimeError, match="lossn_mci"):
        lossn_manjunath([1.0, 1.0], [[1.0, 1.0]], [50.0], max_live_states=4)


def _fcr_lossn(drop1=True, drop2=True):
    model = Network('FCR Loss Network')
    source = Source(model, 'Source')
    delay = Delay(model, 'Delay')
    sink = Sink(model, 'Sink')
    c1 = OpenClass(model, 'Class1', 0)
    c2 = OpenClass(model, 'Class2', 1)
    source.set_arrival(c1, Exp.fit_rate(LAMBDA1))
    source.set_arrival(c2, Exp.fit_rate(LAMBDA2))
    delay.set_service(c1, Exp.fit_rate(MU1))
    delay.set_service(c2, Exp.fit_rate(MU2))
    P = model.init_routing_matrix()
    P.set(c1, c1, source, delay, 1.0)
    P.set(c1, c1, delay, sink, 1.0)
    P.set(c2, c2, source, delay, 1.0)
    P.set(c2, c2, delay, sink, 1.0)
    model.link(P)
    fcr = model.add_region(delay)
    fcr.set_global_max_jobs(GMAX)
    fcr.set_class_max_jobs(c1, C1MAX)
    fcr.set_class_max_jobs(c2, C2MAX)
    fcr.set_drop_rule(c1, drop1)
    fcr.set_drop_rule(c2, drop2)
    return model


@pytest.mark.parametrize("drop1,drop2", [(True, False), (False, True), (False, False)])
def test_mixed_drop_rules_are_refused(drop1, drop2):
    """
    The loss-network path requires EVERY class to be dropped, as the reference
    tests with all(regionrule(1,:) == DROP).

    A mixed region is not a loss network at all: its blocked class occupies the
    region while it waits, so the per-class loss probabilities Kelly's truncation
    implies are not the ones the model implies, and there is no truncated product
    form to evaluate. The failure mode being guarded against is a wrong NUMBER,
    not an error -- a dispatch that checked only class 0 would solve
    (DROP, WAITQ) as though the second class were dropped too.

    Verified against MATLAB SolverNC on all four configurations: only [1 1]
    solves. Note the DIAGNOSTIC differs by codebase: MATLAB and the JAR name the
    WAITQ policy, whereas here the feature gate fires first and reports the
    generic unsupported-feature message. The refusal itself is what matters, and
    it agrees -- so the pattern below admits both wordings, which is also what
    lets this run unchanged under LINE_SOLVER_LANG=java, where the JAR raises.
    """
    with pytest.raises(RuntimeError, match="not supported|does not support"):
        NC(_fcr_lossn(drop1=drop1, drop2=drop2)).get_avg_table()


def test_all_drop_still_solves_after_rewriting_the_rules():
    # The companion to the refusal above: setting both rules explicitly to DROP
    # must still take the exact path, so the all-classes test is not simply
    # rejecting everything.
    solver = NC(_fcr_lossn(drop1=True, drop2=True))
    solver.get_avg_table()
    assert _dispatched(solver._result) == 'lossn.exact'
    np.testing.assert_allclose(np.asarray(solver._result.Q)[1, :], MATLAB_Q,
                               rtol=0, atol=1e-9)


def test_lossn_manjunath_on_the_region_rule_matches_matlab():
    # The rule the analyzer assembles for _fcr_lossn: a global cap plus one
    # per-class cap each, on the SCALED offered load lambda_r / mu_r.
    nu = [LAMBDA1 / MU1, LAMBDA2 / MU2]
    q, loss, lG, _n = lossn_manjunath(nu, [[1, 1], [1, 0], [0, 1]], [GMAX, C1MAX, C2MAX])
    np.testing.assert_allclose(q, MATLAB_Q, rtol=0, atol=1e-12)
    np.testing.assert_allclose(loss, MATLAB_LOSS, rtol=0, atol=1e-12)
    assert abs(lG - MATLAB_LG) < 1e-12


@pytest.mark.parametrize("method", ['default', 'ms'])
def test_solver_nc_exact_path_matches_matlab(method):
    solver = NC(_fcr_lossn(), method=method)
    solver.get_avg_table()
    r = solver._result
    assert _dispatched(r) == 'lossn.exact'
    # QLen is the CARRIED load E[n_r], so Q is QLen and the throughput is
    # lambda_r (1 - Loss_r); dividing QLen by mu again would double-apply
    # Little's law. Station order is Source then Delay.
    np.testing.assert_allclose(np.asarray(r.Q)[1, :], MATLAB_Q, rtol=0, atol=1e-9)
    np.testing.assert_allclose(np.asarray(r.U)[1, :], MATLAB_Q, rtol=0, atol=1e-9)
    np.testing.assert_allclose(np.asarray(r.R)[1, :], [1 / MU1, 1 / MU2], rtol=0, atol=1e-9)
    np.testing.assert_allclose(np.asarray(r.T)[1, :], MATLAB_X, rtol=0, atol=1e-9)
    np.testing.assert_allclose(np.asarray(r.X).ravel(), MATLAB_X, rtol=0, atol=1e-9)
    assert abs(r.lG - MATLAB_LG) < 1e-9
    assert r.it == 1, "the transform is direct"


def test_solver_nc_offered_load_is_scaled_not_the_arrival_rate():
    """
    nu_r is the arrival rate times the mean holding time in the region,
    lambda_r V_r / mu_r, not the bare arrival rate. Class 2 has mu = 0.8, so the
    bare rate would give nu_2 = 0.2 instead of 0.25 and report the blocking of a
    different network -- which is what the pre-port analyzer did.
    """
    bare_q, _l, _g, _n = lossn_manjunath([LAMBDA1, LAMBDA2], [[1, 1], [1, 0], [0, 1]],
                                  [GMAX, C1MAX, C2MAX])
    solver = NC(_fcr_lossn())
    solver.get_avg_table()
    got = np.asarray(solver._result.Q)[1, :]
    np.testing.assert_allclose(got, MATLAB_Q, rtol=0, atol=1e-9)
    assert abs(got[1] - bare_q[1]) > 1e-3, \
        "the unscaled offered load must not reproduce the reference answer"


def test_solver_nc_erlangfp_is_close_but_not_exact():
    # The reduced-load approximation treats the links as independent, so it must
    # land near the truth without reproducing it. Asserting BOTH bounds is what
    # keeps the test honest: an 'erlangfp' that silently dispatched to the
    # transform would pass a one-sided closeness check.
    solver = NC(_fcr_lossn(), method='erlangfp')
    solver.get_avg_table()
    r = solver._result
    assert _dispatched(r) == 'lossn.erlangfp'
    got = np.asarray(r.Q)[1, :]
    np.testing.assert_allclose(got, MATLAB_Q, rtol=0, atol=5e-3)
    assert np.max(np.abs(got - np.asarray(MATLAB_Q))) > 1e-12


def test_solver_nc_mci_brackets_the_exact_answer():
    solver = NC(_fcr_lossn(), method='mci', samples=200000, seed=42)
    solver.get_avg_table()
    r = solver._result
    assert _dispatched(r) == 'lossn.mci'
    np.testing.assert_allclose(np.asarray(r.Q)[1, :], MATLAB_Q, rtol=0, atol=5e-3)
    assert abs(r.lG - MATLAB_LG) < 5e-3
