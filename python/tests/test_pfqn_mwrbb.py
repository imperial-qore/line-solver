"""
Validates the Majumdar-Woodside robust box bounds (pfqn_mwrbb and the
mwba.upper/mwba.lower SolverBA methods) against the numerical results of
S. Majumdar and C.M. Woodside, "Robust bounds and throughput guarantees for
closed multiclass queueing networks", Performance Evaluation 32 (1998) 101-136.

Reference model QNM1 (Section 3.1): two FIFO devices, two closed classes,
N1=2, N2=3, zero think time, V=[[1,1],[1,10]], S=[[0.9,0.9],[0.1,0.1]]
(rows = devices, cols = classes). Paper robust box bounds:
X1-=0.4, X1+=0.74075, X2-=0.37037, X2+=0.71111.
"""

import numpy as np
import pytest

from line_solver.api.pfqn import pfqn_mwrbb
from line_solver import Network, Queue, ClosedClass, Exp, SchedStrategy, SolverBA

TOL = 1e-4


def _qnm1_params():
    V = np.array([[1.0, 1.0], [1.0, 10.0]])
    S = np.array([[0.9, 0.9], [0.1, 0.1]])
    N = np.array([2.0, 3.0])
    Z = np.array([0.0, 0.0])
    return V, S, N, Z


def test_pfqn_mwrbb_qnm1():
    V, S, N, Z = _qnm1_params()
    Xlo, Xup, _ = pfqn_mwrbb(V, S, N, Z)
    assert Xlo[0] == pytest.approx(0.40000, abs=TOL)
    assert Xup[0] == pytest.approx(0.74074, abs=TOL)
    assert Xlo[1] == pytest.approx(0.37037, abs=TOL)
    assert Xup[1] == pytest.approx(0.71111, abs=TOL)
    assert np.all(Xlo <= Xup + TOL)


def test_pfqn_mwrbb_mqnm2_mixed_disciplines():
    # MQNM2 Table 2: 4 classes, 4 devices, N=1 each. Devices: 2=FIFO, 3=NPP,
    # 4=PS; device 1 varies. Priorities (lower=higher): [0,1,1,2].
    prio = np.array([0, 1, 1, 2])
    N = np.array([1, 1, 1, 1])
    Z = np.array([2, 5, 5, 20])
    V = np.array([[1, 1, 1, 1 / 3], [0, 1, 1, 1 / 3],
                  [0, 1, 1, 1 / 3], [0, 4, 4, 2]], float)
    S = np.array([[1, 2, 3, 10], [0, 2, 3, 10],
                  [0, 2, 3, 10], [0, 2, 3, 10]], float)
    codes = [0, 2, 3, 1]  # FIFO, NPP, PP, PS at device 1
    exp = np.array([[0.2376, 0.01778, 0.01404, 0.00832],
                    [0.2376, 0.01061, 0.008333, 0.002459],
                    [0.3333, 0.01212, 0.009524, 0.002459],
                    [0.2376, 0.01844, 0.01376, 0.007666]])
    for di in range(4):
        sched = np.array([codes[di], 0, 2, 1])
        Xlo, Xup, _ = pfqn_mwrbb(V, S, N, Z, sched, prio)
        assert np.max(np.abs(Xlo - exp[di])) < 1e-3
        # upper bounds discipline-independent
        assert Xup == pytest.approx([0.3333, 0.05263, 0.03846, 0.02], abs=1e-4)


def test_pfqn_mwrbb_mqnm1_population_sweep():
    # MQNM1 Table 1: device 1 = preemptive priority, device 2 = FIFO,
    # device 3 = NPP, device 4 = PS; N1=N3=N4=1 and N2 swept. Paper values.
    prio = np.array([0, 1, 1, 2])
    Z = np.array([2, 5, 5, 20])
    V = np.array([[1, 1, 1, 1 / 3], [0, 1, 1, 1 / 3],
                  [0, 1, 1, 1 / 3], [0, 4, 4, 2]], float)
    S = np.array([[1, 2, 3, 10], [0, 2, 3, 10],
                  [0, 2, 3, 10], [0, 2, 3, 10]], float)
    sched = np.array([3, 0, 2, 1])  # PP, FIFO, NPP, PS
    exp = {
        1: ([0.3333, 0.01212, 0.009524, 0.002459], [0.3333, 0.05263, 0.03846, 0.02]),
        2: ([0.3333, 0.01839, 0.007207, 0.0001321], [0.3333, 0.1053, 0.03846, 0.02]),
        3: ([0.3333, 0.02222, 0.005952, 0.0], [0.3333, 0.1161, 0.03846, 0.02]),
    }
    for n2, (elo, eup) in exp.items():
        N = np.array([1, n2, 1, 1])
        Xlo, Xup, _ = pfqn_mwrbb(V, S, N, Z, sched, prio)
        assert Xlo == pytest.approx(elo, abs=1e-4)
        assert Xup == pytest.approx(eup, abs=1e-4)


def test_pfqn_mwrbb_aba_discipline_looser_than_fifo():
    # QNM1 with all stations set to the ABA full-contention code (4) must give
    # a lower (looser) throughput guarantee than FIFO (0), since ABA forces the
    # arrival-contention probability P_cm = 1 (no rate-ratio discount).
    V, S, N, Z = _qnm1_params()
    xlo_fifo, _, _ = pfqn_mwrbb(V, S, N, Z, np.array([0, 0]))
    xlo_aba, _, _ = pfqn_mwrbb(V, S, N, Z, np.array([4, 4]))
    assert xlo_aba[1] == pytest.approx(0.31579, abs=TOL)   # class 2 ABA lower
    assert np.all(xlo_aba <= xlo_fifo + 1e-12)


def test_pfqn_mwrbb_single_class_reduces_to_asymptotic():
    # single class must reduce to Muntz-Wong asymptotic bounds
    D = np.array([0.9, 0.1])
    Np, Zc = 4.0, 1.0
    V = np.array([[1.0], [1.0]])
    S = D.reshape(-1, 1)
    Xlo, Xup, _ = pfqn_mwrbb(V, S, np.array([Np]), np.array([Zc]))
    x_up_ref = min(1.0 / D.max(), Np / (Zc + D.sum()))
    x_lo_ref = Np / (Zc + Np * D.sum())
    assert Xup[0] == pytest.approx(x_up_ref, abs=TOL)
    assert Xlo[0] == pytest.approx(x_lo_ref, abs=TOL)


def _build_qnm1_model():
    model = Network('QNM1')
    A = Queue(model, 'Q1', SchedStrategy.FCFS)
    B = Queue(model, 'Q2', SchedStrategy.FCFS)
    c1 = ClosedClass(model, 'C1', 2, A, 0)
    c2 = ClosedClass(model, 'C2', 3, A, 0)
    A.setService(c1, Exp(1 / 0.9))
    B.setService(c1, Exp(1 / 0.1))
    A.setService(c2, Exp(1 / 0.9))
    B.setService(c2, Exp(1 / 0.1))
    P = model.initRoutingMatrix()
    P.set(c1, c1, A, B, 1.0)
    P.set(c1, c1, B, A, 1.0)
    P.set(c2, c2, A, B, 1.0)
    P.set(c2, c2, B, B, 0.9)
    P.set(c2, c2, B, A, 0.1)
    model.link(P)
    return model


def test_solverba_mwba_qnm1_end_to_end():
    model = _build_qnm1_model()
    # Bound methods live in SolverBA; SolverMVA rejects them by design in
    # all three codebases (MATLAB @SolverMVA/runAnalyzer line_error).
    Tup = SolverBA(model, method='mwba.upper').getAvgTable()
    Tlo = SolverBA(model, method='mwba.lower').getAvgTable()

    # class throughput = Tput at the reference station Q1 (rows 0=Q1/C1, 1=Q1/C2)
    def tput(tbl, irow):
        return float(np.asarray(tbl['Tput'])[irow])

    assert tput(Tup, 0) == pytest.approx(0.74074, abs=TOL)
    assert tput(Tup, 1) == pytest.approx(0.71111, abs=TOL)
    assert tput(Tlo, 0) == pytest.approx(0.40000, abs=TOL)
    assert tput(Tlo, 1) == pytest.approx(0.37037, abs=TOL)
