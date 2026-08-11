"""
Tests for the Maximum Entropy Method (Kouvatsos 1994) in SolverNC.

Validates the me_oqn fixed-point algorithm and the SolverNC method='mem'
wrapper against closed-form results: M/M/1, M/M/1 with Bernoulli feedback
(Jackson exact), GE/M/1 (paper closed form), tandem M/M/1, M/M/c (Erlang-C
exact) and the GE/GE/inf building block.
"""
import math
import warnings

import numpy as np
import pytest

from line_solver import (ClosedClass, Delay, Exp, HyperExp, Network,
                         OpenClass, Queue, Sink, SolverNC, Source)
from line_solver.api.me.me_cqn import me_cqn
from line_solver.api.me.me_oqn import me_oqn

TOL = 1e-8


def _no_routing(M, R):
    return np.zeros((M, M, R))


class TestMeOqn:
    def test_mm1(self):
        L, W, Ca, Cd, lam, rho, _ = me_oqn(
            1, 1, [[0.6]], [[1.0]], [[1.0]], [[1.0]], _no_routing(1, 1))
        assert L[0, 0] == pytest.approx(1.5, abs=TOL)
        assert Cd[0, 0] == pytest.approx(1.0, abs=TOL)
        assert W[0, 0] == pytest.approx(2.5, abs=TOL)

    def test_mm1_feedback(self):
        # M/M/1 with 50% Bernoulli feedback: composite service exponential,
        # hence Jackson-exact L = rho/(1-rho) = 1 at rho = 0.5
        P = _no_routing(1, 1)
        P[0, 0, 0] = 0.5
        L, W, Ca, Cd, lam, rho, _ = me_oqn(
            1, 1, [[1.0]], [[1.0]], [[4.0]], [[1.0]], P)
        assert L[0, 0] == pytest.approx(1.0, abs=TOL)
        assert lam[0, 0] == pytest.approx(2.0, abs=TOL)  # visit-inclusive
        assert rho[0, 0] == pytest.approx(0.5, abs=TOL)

    def test_gem1(self):
        # GE/M/1 with Ca=5, rho=0.5: paper closed form
        # L = rho/2*(Ca+1) + rho^2*(Ca+Cs)/(2*(1-rho)) = 3.0
        L, _, _, _, _, _, _ = me_oqn(
            1, 1, [[1.0]], [[5.0]], [[2.0]], [[1.0]], _no_routing(1, 1))
        assert L[0, 0] == pytest.approx(3.0, abs=TOL)

    def test_tandem_mm1(self):
        P = _no_routing(2, 1)
        P[0, 1, 0] = 1.0
        L, _, _, _, _, _, _ = me_oqn(
            2, 1, [[0.5], [0.0]], [[1.0], [1.0]], [[1.0], [0.8]],
            [[1.0], [1.0]], P)
        assert L[0, 0] == pytest.approx(1.0, abs=TOL)
        assert L[1, 0] == pytest.approx(5.0 / 3.0, abs=TOL)

    def test_mmc_erlang_c(self):
        # M/M/3 with lambda=2, mu=1: Erlang-C exact
        c, a = 3, 2.0
        rc = a / c
        p0 = 1.0 / (sum(a ** n / math.factorial(n) for n in range(c))
                    + a ** c / math.factorial(c) / (1 - rc))
        L_exact = a + a ** c / math.factorial(c) / (1 - rc) * p0 * rc / (1 - rc)
        L, _, _, Cd, _, rho, _ = me_oqn(
            1, 1, [[2.0]], [[1.0]], [[1.0]], [[1.0]], _no_routing(1, 1), c=[3])
        assert L[0, 0] == pytest.approx(L_exact, abs=TOL)
        assert Cd[0, 0] == pytest.approx(1.0, abs=TOL)
        assert rho[0, 0] == pytest.approx(rc, abs=TOL)

    def test_infinite_server(self):
        # GE/GE/inf: L = lambda/mu, departures inherit the arrival scv
        L, _, _, Cd, _, _, _ = me_oqn(
            1, 1, [[3.0]], [[1.0]], [[2.0]], [[1.0]], _no_routing(1, 1),
            c=[np.inf])
        assert L[0, 0] == pytest.approx(1.5, abs=TOL)
        assert Cd[0, 0] == pytest.approx(1.0, abs=TOL)

    def test_two_class_tandem_parity(self):
        # Two-class tandem; reference values from the MATLAB implementation
        P = _no_routing(2, 2)
        P[0, 1, 0] = 1.0
        P[0, 1, 1] = 1.0
        L, _, _, _, _, _, _ = me_oqn(
            2, 2, [[0.3, 0.2], [0.0, 0.0]], np.ones((2, 2)),
            [[1.0, 1.0], [0.8, 0.8]], np.ones((2, 2)), P)
        assert L[0, 0] == pytest.approx(0.6, abs=TOL)
        assert L[0, 1] == pytest.approx(0.4, abs=TOL)
        assert L[1, 0] == pytest.approx(1.24, abs=TOL)
        assert L[1, 1] == pytest.approx(0.8266666666666667, abs=TOL)


class TestSolverNCMem:
    def _qlen(self, solver, row):
        q = np.asarray(solver.getAvgQLen())
        return q[row, 0] if q.ndim > 1 else q[row]

    def test_mem_mm1(self):
        model = Network('mm1')
        source = Source(model, 'Source')
        queue = Queue(model, 'Queue')
        sink = Sink(model, 'Sink')
        oclass = OpenClass(model, 'Class1')
        source.setArrival(oclass, Exp(0.6))
        queue.setService(oclass, Exp(1.0))
        model.link(Network.serialRouting(source, queue, sink))
        solver = SolverNC(model, method='mem')
        assert self._qlen(solver, 1) == pytest.approx(1.5, abs=1e-6)

    def test_mem_mmc_and_delay(self):
        # M/M/3 (Erlang-C exact L=26/9) followed by an IS delay (exact L=1)
        model = Network('mmc_is')
        source = Source(model, 'Source')
        queue = Queue(model, 'Queue')
        queue.setNumberOfServers(3)
        delay = Delay(model, 'Delay')
        sink = Sink(model, 'Sink')
        oclass = OpenClass(model, 'Class1')
        source.setArrival(oclass, Exp(2.0))
        queue.setService(oclass, Exp(1.0))
        delay.setService(oclass, Exp(2.0))
        model.link(Network.serialRouting(source, queue, delay, sink))
        solver = SolverNC(model, method='mem')
        assert self._qlen(solver, 1) == pytest.approx(26.0 / 9.0, abs=1e-6)
        assert self._qlen(solver, 2) == pytest.approx(1.0, abs=1e-6)

    def test_mem_feedback(self):
        # M/M/1 with 50% Bernoulli feedback: Jackson-exact L=1 at rho=0.5
        model = Network('fb')
        source = Source(model, 'Source')
        queue = Queue(model, 'Queue')
        sink = Sink(model, 'Sink')
        oclass = OpenClass(model, 'Class1')
        source.setArrival(oclass, Exp(1.0))
        queue.setService(oclass, Exp(4.0))
        P = model.initRoutingMatrix()
        P.set(oclass, oclass, source, queue, 1.0)
        P.set(oclass, oclass, queue, queue, 0.5)
        P.set(oclass, oclass, queue, sink, 0.5)
        model.link(P)
        solver = SolverNC(model, method='mem')
        assert self._qlen(solver, 1) == pytest.approx(1.0, abs=1e-6)
        t = np.asarray(solver.getAvgTput())
        tput = t[1, 0] if t.ndim > 1 else t[1]
        assert tput == pytest.approx(2.0, abs=1e-6)  # visit-inclusive

    def test_default_keeps_product_form_path(self):
        # All-exponential open model: default must keep the exact
        # product-form normalizing-constant path
        model = Network('mm1pf')
        source = Source(model, 'Source')
        queue = Queue(model, 'Queue')
        sink = Sink(model, 'Sink')
        oclass = OpenClass(model, 'Class1')
        source.setArrival(oclass, Exp(0.8))
        queue.setService(oclass, Exp(1.0))
        model.link(Network.serialRouting(source, queue, sink))
        solver = SolverNC(model)  # method='default'
        assert self._qlen(solver, 1) == pytest.approx(4.0, abs=1e-6)
        assert getattr(solver._result, 'method', None) != 'mem'

    def test_mem_rejects_priority(self):
        from line_solver import SchedStrategy
        model = Network('prio')
        source = Source(model, 'Source')
        queue = Queue(model, 'Queue', SchedStrategy.HOL)
        sink = Sink(model, 'Sink')
        oclass = OpenClass(model, 'Class1')
        source.setArrival(oclass, Exp(0.5))
        queue.setService(oclass, Exp(1.0))
        model.link(Network.serialRouting(source, queue, sink))
        solver = SolverNC(model, method='mem')
        with pytest.raises(Exception):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                solver.getAvgQLen()


if __name__ == '__main__':
    import sys
    sys.exit(pytest.main([__file__, '-v']))


class TestMeCqn:
    def test_cyclic_exact(self):
        # Cyclic 2-station M/M/1, N=3: closed ME product form reduces to
        # the exact BCMP solution on Markovian single-class networks
        P = np.zeros((2, 2, 1))
        P[0, 1, 0] = 1.0
        P[1, 0, 0] = 1.0
        L, W, Ca, Cd, lam, rho, X, _ = me_cqn(
            2, 1, [3], [[1.0], [2.0]], [[1.0], [1.0]], P)
        assert L[0, 0] == pytest.approx(2.2667, abs=1e-3)
        assert L[1, 0] == pytest.approx(0.7333, abs=1e-3)
        assert X[0] == pytest.approx(0.9333, abs=1e-3)
        assert np.sum(L) == pytest.approx(3.0, abs=1e-8)  # population exact

    def test_repairmen_exact(self):
        P = np.zeros((2, 2, 1))
        P[0, 1, 0] = 1.0
        P[1, 0, 0] = 1.0
        L, _, _, _, _, _, X, _ = me_cqn(
            2, 1, [4], [[1.0], [1.25]], [[1.0], [1.0]], P, c=[np.inf, 1])
        assert L[0, 0] == pytest.approx(1.2132, abs=1e-3)
        assert L[1, 0] == pytest.approx(2.7868, abs=1e-3)

    def test_h2_reference(self):
        # GE-type reference values (MATLAB/JAR parity)
        P = np.zeros((2, 2, 1))
        P[0, 1, 0] = 1.0
        P[1, 0, 0] = 1.0
        L, _, _, _, _, _, X, _ = me_cqn(
            2, 1, [3], [[1.0], [2.0]], [[4.0], [1.0]], P)
        assert L[0, 0] == pytest.approx(2.018982, abs=1e-4)
        assert X[0] == pytest.approx(0.820026, abs=1e-4)

    def test_two_class_reference(self):
        # Class-dependent rates; 0.1% off exact CTMC (MATLAB/JAR parity)
        P = np.zeros((2, 2, 2))
        P[0, 1, :] = 1.0
        P[1, 0, :] = 1.0
        L, _, _, _, _, _, X, _ = me_cqn(
            2, 2, [2, 2], [[1.0, 0.8], [2.0, 1.5]], np.ones((2, 2)), P)
        assert L[0, 0] == pytest.approx(1.5690, abs=1e-3)
        assert L[0, 1] == pytest.approx(1.5434, abs=1e-3)
        assert X[0] == pytest.approx(0.4446, abs=1e-3)
        assert X[1] == pytest.approx(0.4140, abs=1e-3)


class TestSolverNCMemClosed:
    def _qlen(self, solver, row):
        q = np.asarray(solver.getAvgQLen())
        return q[row, 0] if q.ndim > 1 else q[row]

    def test_mem_closed_cyclic(self):
        model = Network('cyc')
        q1 = Queue(model, 'Q1')
        q2 = Queue(model, 'Q2')
        cclass = ClosedClass(model, 'C', 3, q1)
        q1.setService(cclass, Exp(1.0))
        q2.setService(cclass, Exp(2.0))
        model.link(Network.serialRouting(q1, q2))
        solver = SolverNC(model, method='mem')
        assert self._qlen(solver, 0) == pytest.approx(2.2667, abs=1e-3)
        assert self._qlen(solver, 1) == pytest.approx(0.7333, abs=1e-3)

    def test_mem_closed_repairmen(self):
        model = Network('rep')
        d1 = Delay(model, 'D1')
        q2 = Queue(model, 'Q2')
        cclass = ClosedClass(model, 'C', 4, d1)
        d1.setService(cclass, Exp(1.0))
        q2.setService(cclass, Exp(1.25))
        model.link(Network.serialRouting(d1, q2))
        solver = SolverNC(model, method='mem')
        assert self._qlen(solver, 0) == pytest.approx(1.2132, abs=1e-3)
        assert self._qlen(solver, 1) == pytest.approx(2.7868, abs=1e-3)

    def test_mem_closed_default_not_routed(self):
        # Closed models keep the exact normalizing-constant default path
        model = Network('cdef')
        q1 = Queue(model, 'Q1')
        q2 = Queue(model, 'Q2')
        cclass = ClosedClass(model, 'C', 3, q1)
        q1.setService(cclass, HyperExp.fitMeanAndSCV(1.0, 4.0))
        q2.setService(cclass, Exp(2.0))
        model.link(Network.serialRouting(q1, q2))
        solver = SolverNC(model)  # method='default'
        solver.getAvgQLen()
        assert getattr(solver._result, 'method', None) != 'mem'

    def test_mem_closed_rejects_multiserver(self):
        model = Network('ms')
        q1 = Queue(model, 'Q1')
        q2 = Queue(model, 'Q2')
        q2.setNumberOfServers(2)
        cclass = ClosedClass(model, 'C', 3, q1)
        q1.setService(cclass, Exp(1.0))
        q2.setService(cclass, Exp(2.0))
        model.link(Network.serialRouting(q1, q2))
        solver = SolverNC(model, method='mem')
        with pytest.raises(Exception):
            solver.getAvgQLen()


class TestSolverNCMemMixed:
    def test_mem_mixed_product_form(self):
        # Open class (Poisson 0.3) traverses Q1->Q2; closed class (N=2)
        # cycles Q1<->Q2; exponential services. Exact mixed MVA reference:
        # L=[1.0822 1.5252; 0.2603 0.4748], X=[0.3, 0.6249]
        from line_solver import SchedStrategy
        model = Network('mix')
        source = Source(model, 'Source')
        q1 = Queue(model, 'Q1', SchedStrategy.PS)
        q2 = Queue(model, 'Q2', SchedStrategy.PS)
        sink = Sink(model, 'Sink')
        oclass = OpenClass(model, 'O')
        cclass = ClosedClass(model, 'C', 2, q1)
        source.setArrival(oclass, Exp(0.3))
        q1.setService(oclass, Exp(1.0))
        q1.setService(cclass, Exp(1.0))
        q2.setService(oclass, Exp(2.0))
        q2.setService(cclass, Exp(2.0))
        P = model.initRoutingMatrix()
        P.set(oclass, oclass, source, q1, 1.0)
        P.set(oclass, oclass, q1, q2, 1.0)
        P.set(oclass, oclass, q2, sink, 1.0)
        P.set(cclass, cclass, q1, q2, 1.0)
        P.set(cclass, cclass, q2, q1, 1.0)
        model.link(P)
        solver = SolverNC(model, method='mem')
        QN = np.asarray(solver.getAvgQLen())
        ref = np.array([[0, 0], [1.0822, 1.5252], [0.2603, 0.4748]])
        assert np.max(np.abs(QN - ref)) < 1e-3


class TestNCMemDefaultSelector:
    """Benchmark-calibrated default selector: open non-Markovian models route
    to MEM except bursty arrivals at low load; closed models route to MEM for
    hypoexponential services or class-dependent FCFS rates only."""

    def _method(self, model):
        solver = SolverNC(model)
        solver.getAvgQLen()
        return getattr(solver._result, 'method', None)

    def test_open_bursty_low_load_not_routed(self):
        from line_solver import APH
        model = Network('bursty')
        source = Source(model, 'Source')
        queue = Queue(model, 'Queue')
        sink = Sink(model, 'Sink')
        oclass = OpenClass(model, 'C')
        source.setArrival(oclass, APH.fitMeanAndSCV(3.0, 64.0))
        queue.setService(oclass, Exp(1.0))
        model.link(Network.serialRouting(source, queue, sink))
        assert self._method(model) != 'mem'

    def test_closed_hyperexponential_not_routed(self):
        model = Network('h2c')
        q1 = Queue(model, 'Q1')
        q2 = Queue(model, 'Q2')
        cclass = ClosedClass(model, 'C', 3, q1)
        q1.setService(cclass, HyperExp.fitMeanAndSCV(1.0, 4.0))
        q2.setService(cclass, Exp(2.0))
        model.link(Network.serialRouting(q1, q2))
        assert self._method(model) != 'mem'
