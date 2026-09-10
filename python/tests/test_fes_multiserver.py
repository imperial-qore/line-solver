"""The FES of a subnetwork holding a multiserver station must satisfy Norton.

`fes_build_isolated` reads the real server count into `mi` and the throughput
and metric kernels handed it straight to `pfqn_mva`. But `pfqn_mva`'s `mi` is
NOT a server count: it enters only as the additive term of the residence-time
recursion, `C(i,s) = L(i,s)*(mi(i) + Qarv)`. For `mi=1` that is the ordinary
arrival theorem; for a c-server station it INFLATES the residence time by c
instead of adding c servers, i.e. it makes the station slower, not faster. The
FES table was therefore silently wrong wherever the aggregated subset contained
a multiserver, and wrong in a direction that still looks like a plausible
queueing model. Both kernels now go through `pfqn_mvams`, which forwards to
`pfqn_mva` when every station is a single server and to the load-dependent
recursion with `mu(i,n)=min(n,S(i))` when one is not.

Norton's theorem is the oracle: replacing a subset of stations by its
flow-equivalent server must leave every throughput unchanged in a closed
product-form network. It is exact under the convolution (`pfqn_conv` reads
`mu_{r,i}(n)` off the same class-dependence handle), so `SolverNC` is the
reference here; AMVA-QD reads the handle approximately and is not.

The single-server case is kept alongside as the control: it passed before the
fix and must keep passing, which is what pins the change to the multiserver
branch alone.
"""
import numpy as np
import pytest

from line_solver import (Network, Queue, Delay, ClosedClass, Exp, SchedStrategy,
                         SolverNC)
from line_solver.io.model_adapter import aggregate_fes


def _build(servers_at_q2):
    """Closed 2-class network: Think -> Q1 (PS, 1 server) -> Q2 (FCFS, c servers).

    Q2 has equal rates across classes, so it is BCMP type 1 and the whole model
    stays product form whatever c is.
    """
    net = Network('fes_multiserver')
    d = Delay(net, 'Think')
    q1 = Queue(net, 'Q1', SchedStrategy.PS)
    q1.set_number_of_servers(1)
    q2 = Queue(net, 'Q2', SchedStrategy.FCFS)
    q2.set_number_of_servers(servers_at_q2)
    c1 = ClosedClass(net, 'C1', 3, d, 0)
    c2 = ClosedClass(net, 'C2', 2, d, 0)
    d.set_service(c1, Exp(1.0))
    d.set_service(c2, Exp(0.5))
    q1.set_service(c1, Exp(2.0))
    q1.set_service(c2, Exp(1.25))
    q2.set_service(c1, Exp(1.0))
    q2.set_service(c2, Exp(1.0))
    P = net.init_routing_matrix()
    for r in range(2):
        P[r][r] = net.serial_routing([d, q1, q2])
    net.link(P)
    return net, q1, q2


def _tput(model):
    return np.asarray(SolverNC(model, method='exact', verbose=False).getAvgTput())


@pytest.mark.parametrize('servers', [1, 3])
def test_fes_aggregation_is_norton_exact(servers):
    net, q1, q2 = _build(servers)
    Xorig = _tput(net)[1, :]

    res = aggregate_fes(net, [q1, q2])
    fes_model = res['fes_model'] if isinstance(res, dict) else res[0]
    Xfes = _tput(fes_model)[1, :]

    assert np.allclose(Xfes, Xorig, atol=1e-10), (
        'Norton violated with %d server(s) at Q2: FES %s vs original %s'
        % (servers, Xfes, Xorig))


def test_multiserver_fes_differs_from_single_server():
    """The multiserver branch must actually be exercised.

    A kernel that ignores the server count returns the same aggregate for c=1
    and c=3; the defect's signature was exactly that kind of insensitivity, so
    the assertion is that the two aggregates DIFFER.
    """
    net1, a1, b1 = _build(1)
    net3, a3, b3 = _build(3)
    r1 = aggregate_fes(net1, [a1, b1])
    r3 = aggregate_fes(net3, [a3, b3])
    m1 = r1['fes_model'] if isinstance(r1, dict) else r1[0]
    m3 = r3['fes_model'] if isinstance(r3, dict) else r3[0]
    X1 = _tput(m1)[1, :]
    X3 = _tput(m3)[1, :]
    assert not np.allclose(X1, X3, atol=1e-6), (
        'the FES is insensitive to the server count: c=1 and c=3 both give %s' % X1)


def test_pfqn_mva_mi_is_not_a_server_count():
    """Pin the semantics the two kernels now avoid relying on.

    `pfqn_mva(L, N, Z, mi)` with mi=c is NOT an M/M/c: it is a single server
    whose residence time carries the factor c. `pfqn_mvams(..., S=c)` is the
    M/M/c. On one 3-server station of unit demand holding one job, the exact
    throughput is 1.0; `pfqn_mva` with mi=3 answers 1/3.
    """
    from line_solver.api.pfqn.mva import pfqn_mva
    from line_solver.api.pfqn.mvald import pfqn_mvams

    L = np.array([[1.0, 1.0]])
    N = np.array([1.0, 0.0])
    Z = np.zeros(2)

    X_ms = np.asarray(pfqn_mvams(np.zeros(2), L, N, Z,
                                 np.ones(1), np.array([3.0]))[0]).ravel()
    assert np.isclose(X_ms[0], 1.0, atol=1e-10), X_ms

    X_mi = np.asarray(pfqn_mva(L, N, Z, np.array([3.0]))[0]).ravel()
    assert np.isclose(X_mi[0], 1.0 / 3.0, atol=1e-10), X_mi
