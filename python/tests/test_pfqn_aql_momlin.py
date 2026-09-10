"""AQL and the moment-linearizer branch of getMomentTable.

pfqn_aql was a Schweitzer refinement of pfqn_bs under the AQL name until
2026-08-01, so the reference values below are the guard against that: they come
from MATLAB pfqn_aql([0.5 0.3;0.4 0.6],[4 3],[1 0.5]) printed at %.12g. The
solver-level tests check the dispatch that was wired in the same change.
"""

import numpy as np
import pytest

from line_solver import (ClosedClass, Delay, Exp, Network, Queue,
                         SchedStrategy, SolverMVA)
from line_solver.api.pfqn import pfqn_aql

# MATLAB pfqn_aql on L = [0.5 0.3; 0.4 0.6], N = [4 3], Z = [1 0.5]
AQL_L = np.array([[0.5, 0.3], [0.4, 0.6]])
AQL_N = np.array([4.0, 3.0])
AQL_Z = np.array([1.0, 0.5])
AQL_X = np.array([1.01427076478, 0.852510135342])
AQL_Q = np.array([[1.37097188005, 0.770649024465],
                  [1.61475736928, 1.80309587945]])
AQL_U = np.array([[0.50713538239, 0.255753040603],
                  [0.405708305912, 0.511506081205]])


def test_pfqn_aql_matches_matlab_reference():
    XN, CN, QN, UN, RN, TN, AN = pfqn_aql(AQL_L, AQL_N, AQL_Z)
    assert np.allclose(XN.flatten(), AQL_X, rtol=1e-10)
    assert np.allclose(QN, AQL_Q, rtol=1e-10)
    assert np.allclose(UN, AQL_U, rtol=1e-10)


def test_pfqn_aql_conserves_the_closed_population():
    XN, CN, QN, UN, RN, TN, AN = pfqn_aql(AQL_L, AQL_N, AQL_Z)
    # tol defaults to 1e-7 on the relative queue-length change, so the
    # population identity holds to about that much, not to machine precision
    assert np.allclose(QN.sum(axis=0) + XN.flatten() * AQL_Z, AQL_N, atol=1e-6)


def _two_class_model():
    model = Network('aqlModel')
    d = Delay(model, 'Think')
    q1 = Queue(model, 'Q1', SchedStrategy.PS)
    q2 = Queue(model, 'Q2', SchedStrategy.PS)
    c1 = ClosedClass(model, 'C1', 4, d)
    c2 = ClosedClass(model, 'C2', 3, d)
    d.setService(c1, Exp(1))
    d.setService(c2, Exp(2))
    q1.setService(c1, Exp(2))
    q1.setService(c2, Exp(3))
    q2.setService(c1, Exp(2.5))
    q2.setService(c2, Exp(1.7))
    model.link(Network.serialRouting(d, q1, q2))
    return model


def test_solver_mva_dispatches_aql():
    model = _two_class_model()
    solver = SolverMVA(model, 'aql')
    table = solver.getAvgTable()
    assert 'aql' in solver.listValidMethods()
    qlen = np.asarray(table['QLen'], dtype=float)
    # MATLAB SolverMVA(model,'aql') on the same model
    assert np.allclose(sorted(np.round(qlen, 4)),
                       sorted([1.0097, 0.42181, 1.4391, 0.87775, 1.5513, 1.7004]),
                       atol=1e-3)


def test_aql_is_not_advertised_on_multiserver():
    model = Network('msModel')
    d = Delay(model, 'Think')
    q = Queue(model, 'Q1', SchedStrategy.PS)
    q.setNumberOfServers(3)
    c = ClosedClass(model, 'C1', 5, d)
    d.setService(c, Exp(1))
    q.setService(c, Exp(2))
    model.link(Network.serialRouting(d, q))
    methods = SolverMVA(model).listValidMethods()
    assert 'aql' not in methods and 'amva.aql' not in methods
    assert 'tay' not in methods and 'amva.tay' not in methods


def test_moment_table_momlin_branch():
    model = _two_class_model()
    solver = SolverMVA(model)
    exact_table, exact_mom = solver.getMomentTable(2, 'exact')
    ml_table, ml_mom = solver.getMomentTable(2, 'momlin')
    assert getattr(ml_mom['qlen'], 'method', None) == 'momlin'
    assert getattr(exact_mom['qlen'], 'method', None) is None
    # the cross-station tensor is the momlin branch's own output
    assert np.shape(ml_mom['qlen'].QCovFull) == (2, 2, 2, 2)
    # same shape, and an approximation of the same quantity: the means agree to
    # a few percent on this small model, which is the AMVA error, not a bug
    qe = np.asarray(exact_table['QLen'], dtype=float)
    qm = np.asarray(ml_table['QLen'], dtype=float)
    assert qe.shape == qm.shape
    assert np.allclose(qe, qm, rtol=0.12)


def test_moment_table_rejects_unknown_method():
    with pytest.raises(ValueError):
        SolverMVA(_two_class_model()).getMomentTable(2, 'bogus')
