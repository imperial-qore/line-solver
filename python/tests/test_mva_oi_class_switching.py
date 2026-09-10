"""The exact OI mean-value analyzer must REFUSE class switching, not answer it.

`nc_is_oi_model` (and its python twin `find_oi_stations`) tests STATION SHAPES:
an OI/PAS station with an all-zero swap graph plus product-form neighbours. It
says nothing about chains, so a class-switching model passes it and reaches the
analyzer. The CMVA recursion there is driven by the per-class population vector
`sn.njobs`, which class switching makes meaningless: a class that only ever
appears mid-chain carries `njobs = 0`, so the OI station is analyzed as if empty
and the whole population is parked in the delay.

`solver_nc_oi_analyzer` always carried the guard; `solver_mva_oi_analyzer` did
not, and returned that degenerate answer SILENTLY -- zero throughput at the OI
station where the exact answer is nonzero. The model below is the smallest shape
that exhibits it, and is also the shape every cache model takes (read -> hit |
miss is class switching by construction).

The restriction is an implementation limit, not a theorem: the OI product form
survives class switching, but the convolution would have to run on the CHAIN
population over a per-class count lattice rather than on a fixed per-class N.
"""
import numpy as np
import pytest

from line_solver import (Network, Queue, Delay, ClosedClass, Exp, SchedStrategy,
                         SolverCTMC, SolverMVA, SolverNC)

N = 4        # closed population
Z = 1.0      # think time
MU_HIT = 4.0
MU_MISS = 1.0
P_HIT = 0.6


def _model():
    """Delay -> OI queue, with Read switching to Hit or Miss on the way in."""
    net = Network('oi_class_switching')
    think = Delay(net, 'Think')
    srv = Queue(net, 'Srv', SchedStrategy.OI)
    read = ClosedClass(net, 'Read', N, think, 0)
    hit = ClosedClass(net, 'Hit', 0, think, 0)
    miss = ClosedClass(net, 'Miss', 0, think, 0)
    think.set_service(read, Exp(1.0 / Z))

    # python-native hands svcRateFun 0-BASED class indices (MATLAB 1-based);
    # see _kb/07-cross-language-parity.md.
    sigma = [0.0, MU_HIT, MU_MISS]

    def ordrate(c):
        c = np.atleast_1d(np.asarray(c, dtype=int)).ravel()
        return float(sum(sigma[int(j)] for j in c if 0 <= int(j) < len(sigma)))

    srv.set_service(ordrate)
    srv.set_number_of_servers(N)

    P = net.init_routing_matrix()
    P.set(read, hit, think, srv, P_HIT)
    P.set(read, miss, think, srv, 1.0 - P_HIT)
    P.set(hit, read, srv, think, 1.0)
    P.set(miss, read, srv, think, 1.0)
    net.link(P)
    return net


def _refused(solver_factory):
    with pytest.raises(Exception) as exc:
        solver_factory().getAvgTput()
    return str(exc.value)


def test_mva_refuses_class_switching_rather_than_returning_an_empty_station():
    msg = _refused(lambda: SolverMVA(_model(), method='exact', verbose=False))
    assert 'one class per chain' in msg, msg


def test_nc_refuses_the_same_model_the_same_way():
    msg = _refused(lambda: SolverNC(_model(), method='exact', verbose=False))
    assert 'one class per chain' in msg, msg


def test_ctmc_still_solves_it_and_the_oi_station_is_not_empty():
    """The model itself is fine -- it is the mean-value recursion that cannot
    take it. CTMC walks the ordered microstate and calls mu(c) directly, so it
    is the oracle showing the refused answer would not have been zero."""
    X = np.asarray(SolverCTMC(_model(), verbose=False).getAvgTput())
    assert X[1, 1] > 1e-6 and X[1, 2] > 1e-6, X
