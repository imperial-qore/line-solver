"""The load-dependent QRF arms on delay, multiserver and load-dependent models.

WHAT CHANGED AND WHY IT IS NOT AN APPROXIMATION
-----------------------------------------------
``qrf.mmi.ld`` and ``qrf.mmi.linear`` carry a scaling ``alpha(i,n)`` that
multiplies every rate out of station i at population n, completions and
background phase changes alike. That is exactly the rate law of an infinite
server (``alpha = n``), of a c-server station (``alpha = min(n,c)``) and of
limited load dependence, so deriving alpha from the model makes the relaxed
chain the model's OWN chain rather than an approximation of it. The two arms
therefore serve models the rest of the QRF family still refuses.

THE ORACLE IS EXACTNESS, NOT CLOSENESS
--------------------------------------
At M = 2 the pairwise joint of a closed chain is fully determined by the
marginal, so the QRF polytope is tight and the answer must EQUAL the exact CTMC
one (see the same argument in test_solver_ba_qrf_path.py). Every two-station
fixture here is asserted at 1e-9 against SolverCTMC, which pins the alpha
derivation, the BN readout and the utilization normalizer at once. A looser
tolerance would let a wrong alpha through.

The N = 1 fixture is a second exact oracle that needs no CTMC: ``min(1,c) = 1``,
so a c-server model at population 1 is the c = 1 model, and the two runs must
agree to the digit.
"""
import numpy as np
import pytest

from line_solver import (Network, Queue, Delay, ClosedClass, SchedStrategy, Exp,
                         Erlang, SolverBA, SolverCTMC)
from line_solver.api.sn.sn_to_qrf_alpha import sn_to_qrf_alpha

LD_ARMS = ['qrf.mmi.ld', 'qrf.mmi.linear']
ALPHA_FREE_ARMS = ['qr', 'qrf.mmi', 'qrf.mem', 'qrf.bethe']


def _cqn(N, servers=(1, 1), rates=(1.0, 1.5), delay0=False, lld=None, phases=(1, 1)):
    """Two stations in a cycle; station 0 is a Delay when delay0."""
    m = Network('ld')
    nodes = []
    for i in range(len(rates)):
        if delay0 and i == 0:
            nodes.append(Delay(m, 'D%d' % i))
        else:
            q = Queue(m, 'Q%d' % i, SchedStrategy.FCFS)
            q.setNumberOfServers(int(servers[i]))
            nodes.append(q)
    cls = ClosedClass(m, 'C', N, nodes[0], 0)
    for i, r in enumerate(rates):
        if phases[i] > 1:
            nodes[i].setService(cls, Erlang.fitMeanAndOrder(1.0 / r, phases[i]))
        else:
            nodes[i].setService(cls, Exp(r))
    if lld is not None:
        nodes[lld[0]].setLoadDependence(np.asarray(lld[1], dtype=float))
    m.link(Network.serialRouting(*nodes))
    return m


def _avg(solver):
    return (np.asarray(solver.getAvgQLen(), dtype=float).ravel(),
            np.asarray(solver.getAvgUtil(), dtype=float).ravel(),
            np.asarray(solver.getAvgTput(), dtype=float).ravel())


def test_alpha_is_min_n_c_at_a_multiserver_station():
    alpha, msg, ld, peak = sn_to_qrf_alpha(_cqn(4, servers=(3, 1)).getStruct())
    assert msg == '' and ld
    assert np.allclose(alpha[0], [1, 2, 3, 3])
    assert np.allclose(alpha[1], [1, 1, 1, 1])
    # The normalizer is the DECLARED peak, so three servers even though the
    # reachable alpha would peak at 3 here too; the N < c case is below.
    assert np.allclose(peak, [3, 1])


def test_alpha_peak_is_the_declared_servers_not_the_reachable_maximum():
    # c = 3 with N = 2: alpha reaches 2, but the station still has 3 servers and
    # LINE reports U = T*S/c. Normalizing by max(alpha) would overstate U by 3/2.
    alpha, msg, _ld, peak = sn_to_qrf_alpha(_cqn(2, servers=(3, 1)).getStruct())
    assert msg == ''
    assert np.allclose(alpha[0], [1, 2])
    assert np.allclose(peak, [3, 1])


def test_alpha_is_n_at_a_delay():
    alpha, msg, ld, peak = sn_to_qrf_alpha(_cqn(3, delay0=True).getStruct())
    assert msg == '' and ld
    assert np.allclose(alpha[0], [1, 2, 3])
    assert np.isinf(peak[0]) and peak[1] == 1


def test_alpha_composes_lld_with_the_server_count():
    m = _cqn(3, servers=(2, 1), lld=(0, [1.0, 1.5, 2.0]))
    alpha, msg, ld, peak = sn_to_qrf_alpha(m.getStruct())
    assert msg == '' and ld
    # min(n,2) times the lld column, and the peak is the PRODUCT c * max(lld).
    assert np.allclose(alpha[0], [1 * 1.0, 2 * 1.5, 2 * 2.0])
    assert np.allclose(peak, [4.0, 1.0])


def test_alpha_is_all_ones_on_a_single_server_model():
    """The load-independent model must reach the ld arms unchanged, which is
    what keeps this change a no-op there."""
    alpha, msg, ld, peak = sn_to_qrf_alpha(_cqn(3).getStruct())
    assert msg == '' and not ld
    assert np.allclose(alpha, 1.0)
    assert np.allclose(peak, 1.0)


@pytest.mark.parametrize('kw', [dict(servers=(2, 1)), dict(delay0=True)])
def test_phasetype_where_several_jobs_are_served_at_once_is_refused(kw):
    """The QRF local state carries ONE phase per station, which describes one
    job in service and no more, so a PH multiserver would answer a different
    chain and the number would bound nothing."""
    m = _cqn(2, phases=(2, 1), **kw)
    _alpha, msg, ld, _peak = sn_to_qrf_alpha(m.getStruct())
    assert 'one phase per station' in msg
    assert ld, 'ld must stay set through the refusal'
    with pytest.raises(ValueError, match='one phase per station'):
        SolverBA(m, method='qrf.mmi.ld').getAvgQLen()


@pytest.mark.parametrize('meth', ALPHA_FREE_ARMS)
@pytest.mark.parametrize('kw', [dict(servers=(2, 1)), dict(delay0=True),
                                dict(lld=(0, [1.0, 2.0]))])
def test_alpha_free_arms_refuse_and_name_the_arms_that_serve(meth, kw):
    m = _cqn(2, **kw)
    with pytest.raises(ValueError, match="Use 'qrf.mmi.ld' or 'qrf.mmi.linear'"):
        SolverBA(m, method=meth).getAvgQLen()


@pytest.mark.parametrize('meth', LD_ARMS)
@pytest.mark.parametrize('label,kw', [
    ('delay', dict(N=3, delay0=True)),
    ('c=2', dict(N=3, servers=(2, 1))),
    ('c=3', dict(N=4, servers=(3, 1))),
    ('lld', dict(N=3, lld=(0, [1.0, 1.5, 2.0]))),
    ('single-server', dict(N=3)),
])
def test_two_station_answers_are_exact(meth, label, kw):
    """M = 2 makes the polytope tight, so these are equalities, not bounds."""
    m = _cqn(**kw)
    q, u, t = _avg(SolverBA(m, method=meth))
    qc, uc, tc = _avg(SolverCTMC(m))
    assert np.allclose(q, qc, atol=1e-9), '%s QLen' % label
    assert np.allclose(u, uc, atol=1e-9), '%s Util' % label
    assert np.allclose(t, tc, atol=1e-9), '%s Tput' % label


@pytest.mark.parametrize('meth', LD_ARMS)
def test_multiserver_at_population_one_is_the_single_server_model(meth):
    """min(1,c) = 1, so the two chains are identical and only the utilization
    normalizer differs. No CTMC needed: the oracle is the other run."""
    q1, u1, t1 = _avg(SolverBA(_cqn(1, servers=(1, 1)), method=meth))
    q3, u3, t3 = _avg(SolverBA(_cqn(1, servers=(3, 1)), method=meth))
    assert np.allclose(q1, q3, atol=1e-9)
    assert np.allclose(t1, t3, atol=1e-9)
    assert np.allclose(u3[0] * 3.0, u1[0], atol=1e-9)


def test_list_valid_methods_keeps_the_ld_arms_on_a_delay_model():
    """A caller enumerating the list must still see the only two bound methods
    the model has; dropping them with the rest of the family would hide them."""
    methods = list(SolverBA(_cqn(2, delay0=True)).listValidMethods())
    for m in LD_ARMS:
        assert m in methods
    for m in ALPHA_FREE_ARMS + ['qrf.bas', 'qrf.rsrd']:
        assert m not in methods
