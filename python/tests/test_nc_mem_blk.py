"""
Tests for the finite-buffer Maximum Entropy Method in SolverNC.

Validates the censored GE/GE/c/K;N building block of Kouvatsos (1994),
Section 4.1, against the exact Markovian queue-length distributions, and the
transfer-blocking network algorithm of Tahilramani, Manjunath and Bose (1999)
against the results published in their Table 1.
"""
import math

import numpy as np
import pytest

from line_solver import (DropStrategy, Exp, Network, OpenClass, Queue,
                         SchedStrategy, Sink, SolverNC, Source)
from line_solver.api.me.me_gegecn import me_gegecn
from line_solver.api.me.me_oqn import me_oqn
from line_solver.api.me.me_oqn_blk import RULE_BAS, RULE_LOSS, me_oqn_blk

TOL = 1e-10


@pytest.mark.parametrize('rho', [0.3, 0.8, 1.0, 1.5],
                         ids=['rho0.3', 'rho0.8', 'rho1.0', 'rho1.5'])
@pytest.mark.parametrize('nbuf', [1, 2, 5, 10],
                         ids=['nbuf1', 'nbuf2', 'nbuf5', 'nbuf10'])
def test_censored_single_server_is_exact(rho, nbuf):
    """The censored GE/GE/1/0;N reduces to M/M/1/N on Markovian streams."""
    p, L, U, PB, _ = me_gegecn(rho, 1.0, 1.0, 1.0, 1, 0, nbuf)
    n = np.arange(nbuf + 1)
    un = rho ** n
    pex = un / np.sum(un)
    assert np.max(np.abs(p - pex)) < TOL
    assert abs(L - float(np.sum(n * pex))) < TOL
    # A Poisson stream is blocked exactly when the buffer is full (PASTA)
    assert abs(PB - pex[-1]) < TOL
    assert abs(U - (1 - pex[0])) < TOL


@pytest.mark.parametrize('c', [2, 3, 5], ids=['c2', 'c3', 'c5'])
def test_censored_multiserver_is_exact(c):
    """The censored GE/GE/c/0;N reduces to M/M/c/N on Markovian streams."""
    for nbuf in (c, c + 3, c + 8):
        lam = 0.7 * c
        p, L, U, PB, _ = me_gegecn(lam, 1.0, 1.0, 1.0, c, 0, nbuf)
        n = np.arange(nbuf + 1)
        un = np.zeros(nbuf + 1)
        for k in range(nbuf + 1):
            if k <= c:
                un[k] = lam ** k / math.factorial(k)
            else:
                un[k] = lam ** c / math.factorial(c) * (lam / c) ** (k - c)
        pex = un / np.sum(un)
        assert np.max(np.abs(p - pex)) < TOL
        assert abs(L - float(np.sum(n * pex))) < TOL
        assert abs(U - float(np.sum(np.minimum(n, c) * pex)) / c) < TOL
        assert abs(PB - pex[-1]) < TOL


@pytest.mark.parametrize('rho', [0.5, 0.9, 1.2],
                         ids=['rho0.5', 'rho0.9', 'rho1.2'])
def test_single_station_with_loss(rho):
    """The network algorithm on a single M/M/1/N with loss is exact."""
    nbuf = 4
    Q, W, T, U, _, _, PBa, _, _ = me_oqn_blk(
        1, [rho], [1.0], [1.0], [1.0], [[0.0]], [1.0], [nbuf], [RULE_LOSS])
    n = np.arange(nbuf + 1)
    pex = rho ** n / np.sum(rho ** n)
    assert abs(Q[0] - float(np.sum(n * pex))) < 1e-6
    assert abs(T[0] - rho * (1 - pex[-1])) < 1e-5
    assert abs(PBa[0] - pex[-1]) < 1e-5


def test_unbounded_reduces_to_me_oqn():
    """With every buffer unbounded the algorithm reproduces me_oqn."""
    M = 3
    lam0 = np.array([1.0, 0.0, 0.0])
    ca0 = np.ones(3)
    mu = np.array([3.0, 2.5, 2.0])
    cs = np.array([1.0, 2.0, 0.5])
    P = np.array([[0.0, 0.6, 0.4], [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]])
    c = np.ones(3)
    Q = me_oqn_blk(M, lam0, ca0, mu, cs, P, c, [np.inf] * 3, [RULE_LOSS] * 3)[0]
    L = me_oqn(M, 1, lam0.reshape(3, 1), ca0.reshape(3, 1), mu.reshape(3, 1),
               cs.reshape(3, 1), P.reshape(3, 3, 1), c, np.zeros(3, dtype=bool))[0]
    assert np.max(np.abs(Q - L.ravel())) < 1e-10


def test_transfer_blocking_against_published_table():
    """Feed-forward network of Tahilramani, Manjunath and Bose (1999), Table 1.

    The published values come from an independent implementation of the same
    algorithm, so agreement is asserted at the level of the algorithm rather
    than bitwise: 5 percent on the queue lengths, 2 percent on the throughputs.
    """
    Q, W, T, U, _, _, _, _, _ = me_oqn_blk(
        3, [1.5, 0, 0], [2.0, 1, 1], [2, 2, 2], [1, 1, 1],
        [[0, 0.4, 0.4], [0, 0, 0.5], [0, 0, 0]], [1, 3, 1], [5, 4, 3],
        [RULE_BAS] * 3)
    k_paper = np.array([1.702, 0.261, 0.578])
    t_paper = np.array([1.284, 0.514, 0.771])
    assert np.max(np.abs(Q - k_paper) / k_paper) < 0.05
    assert np.max(np.abs(T - t_paper) / t_paper) < 0.02


@pytest.mark.parametrize('nbuf,qlen,tput', [(2, 0.85245902, 0.59016393),
                                            (4, 1.56306521, 0.70252261),
                                            (6, 2.14243372, 0.74692668)],
                         ids=['nbuf2', 'nbuf4', 'nbuf6'])
def test_solver_nc_mem_accepts_finite_capacity(nbuf, qlen, tput):
    """SolverNC method='mem' solves M/M/1/N exactly, the CTMC values being
    the reference."""
    model = Network('mm1n')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue1', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass = OpenClass(model, 'Class1')
    source.setArrival(oclass, Exp(0.8))
    queue.setService(oclass, Exp(1.0))
    queue.setCapacity(nbuf)
    model.link(Network.serialRouting(source, queue, sink))
    table = SolverNC(model, method='mem').getAvgTable()
    assert abs(float(table.QLen[1]) - qlen) < 1e-6
    assert abs(float(table.Tput[1]) - tput) < 1e-5


def test_hypoexponential_is_rejected():
    """The GE domain is enforced rather than approximated."""
    with pytest.raises(ValueError):
        me_gegecn(0.5, 1.0, 1.0, 0.5, 1, 0, 4)
    with pytest.raises(ValueError):
        me_oqn_blk(1, [0.5], [1.0], [1.0], [0.5], [[0.0]], [1.0], [4], [RULE_LOSS])
