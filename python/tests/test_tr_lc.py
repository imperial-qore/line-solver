"""Load concealment as a model TRANSFORMATION (options.config['transform']='lc').

Birman-Kogan Algorithm 2 exists twice in the tree and the two are not the same
computation. ``pfqn_bklc`` is the KERNEL: it sweeps chains on a demand matrix
with pfqn_mva as the inner single-chain solve. ``transform='lc'`` is the
TRANSFORMATION: the same sweep, but each single-chain subproblem is a real
single-class Network solved by whichever solver was called.

THE ORACLE IS THE KERNEL'S OWN FIXED POINT. On a model whose per-chain
subproblem is single-class, PS and exponential, the chain aggregation is exact
and the concealment is the identity operation on demands, so the two must reach
the SAME fixed point in the SAME number of sweeps. Agreement on the fixed point
alone is not enough: a Jacobi sweep converges to the same point in a different
number of sweeps, and sweep count is what the four codebases are pinned on.

The MATLAB twin is line-test.git test_tr_lc.m, the C++ twin is
cpp/tests/test_tr_lc.cpp and the JAR twin is LcStrategyTest.java.
"""

import numpy as np
import pytest

from line_solver import (ClosedClass, Delay, Exp, Network, Queue, SchedStrategy,
                         SolverCTMC, SolverMVA)
from line_solver.api.pfqn import pfqn_bklc

#: The demand matrix, populations and think times the model below reduces to.
L = np.array([[0.0, 0.0], [0.5, 1.0 / 3], [1.0 / 3, 2.0 / 3], [5.0 / 9, 0.4]])
N = np.array([4.0, 3.0])
Z = np.array([1.0, 0.5])


def build():
    """Delay -> Q1 -> Q2 -> Q3 -> Delay, two closed classes, one chain each."""
    m = Network('lc')
    d = Delay(m, 'Think')
    q1 = Queue(m, 'Q1', SchedStrategy.PS)
    q2 = Queue(m, 'Q2', SchedStrategy.PS)
    q3 = Queue(m, 'Q3', SchedStrategy.PS)
    c1 = ClosedClass(m, 'C1', 4, d, 0)
    c2 = ClosedClass(m, 'C2', 3, d, 0)
    d.setService(c1, Exp(1.0))
    d.setService(c2, Exp(2.0))
    q1.setService(c1, Exp(2.0))
    q1.setService(c2, Exp(3.0))
    q2.setService(c1, Exp(3.0))
    q2.setService(c2, Exp(1.5))
    q3.setService(c1, Exp(1.8))
    q3.setService(c2, Exp(2.5))
    P = m.initRoutingMatrix()
    P.set(c1, c1, Network.serialRouting(d, q1, q2, q3))
    P.set(c2, c2, Network.serialRouting(d, q1, q2, q3))
    m.link(P)
    return m


def kernel():
    out = pfqn_bklc(L, N, Z, 'mva', 1e-10, 1000)
    return np.asarray(out[0]).ravel(), int(out[3])


def test_the_model_reduces_to_the_demands_the_kernel_is_given():
    sn = build().getStruct()
    assert int(sn.nclasses) == 2
    assert int(sn.nchains) == 2


def test_transformation_reaches_the_kernel_fixed_point():
    Xk, _ = kernel()
    Xt = np.asarray(SolverCTMC(build(), config={'transform': 'lc'}).getAvgSysTput()).ravel()
    assert Xt.shape == Xk.shape
    assert np.allclose(Xt, Xk, rtol=1e-8), (Xt, Xk)


def test_the_sweep_is_gauss_seidel_not_jacobi():
    # Same fixed point in the same number of sweeps. Jacobi coupling reaches the
    # same point at a DIFFERENT sweep count, which is what this catches.
    _, itk = kernel()
    s = SolverCTMC(build(), config={'transform': 'lc'})
    s.getAvgSysTput()
    assert int(s._result.iter) == itk


def test_the_transformed_solve_names_itself():
    s = SolverCTMC(build(), config={'transform': 'lc'})
    s.getAvgSysTput()
    assert 'lc' in str(s._result.method)


def test_a_second_solver_reaches_the_same_fixed_point():
    """The transformation is hosted by the solver, not owned by one.

    SolverMVA and SolverCTMC drive the SAME strategy; both are exact on a
    single-class PS network, so both must land on the kernel's fixed point in
    the kernel's number of sweeps. SolverNC also hosts it and reaches the same
    decomposition, but its single-chain solve is itself approximate (it lands
    within ~4e-6) and costs ~20s per sweep in python, so it is deliberately not
    asserted here.
    """
    Xk, itk = kernel()
    s = SolverMVA(build(), config={'transform': 'lc'})
    X = np.asarray(s.getAvgSysTput()).ravel()
    assert np.allclose(X, Xk, rtol=1e-8), (X, Xk)
    assert int(s._result['iter']) == itk


def test_an_unknown_token_is_refused_by_name():
    with pytest.raises(ValueError, match='not a known model transformation'):
        SolverCTMC(build(), config={'transform': 'bogus'}).getAvgSysTput()
