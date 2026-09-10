"""SolverNC must reach its exact OI analyzer, and must refuse the methods that cannot.

An order-independent station carries a rank rate mu(n) of the whole per-class
occupancy vector, so `hasProductFormSolution()` reads false. Two consequences
were live:

- `method='exact'` was rejected by the product-form guard BEFORE `runAnalyzer`
  could route the model to `solver_nc_oi_analyzer` (exact, `pfqn_ncoi`). The
  guard now carries the same OI/PAS exemption `mvaDispatch` already had via its
  `hasOIorPAS` test.
- every other NC method reads `sn.rates`, which holds only the single-job rate
  `mu([r])`: the rank rate was SILENTLY DROPPED and the answer was that of an
  ordinary queue -- a wrong number, not a coarse one. `comom` returned
  [0.535, 0.465] where the answer is [0.985, 1.264]. Those methods are refused
  by name now. `is`/`sampling` stay admissible because
  `solver_nc_pas_is_analyzer` does read mu(n).

The oracle is SolverCTMC, which walks the ordered microstate and calls the rate
function directly.
"""
import numpy as np
import pytest

from line_solver import (Network, Queue, Delay, ClosedClass, Exp, SchedStrategy,
                         SolverCTMC, SolverNC, SolverMVA)

CSRV = 3
CAP1 = 1     # at most one class-1 job may occupy a server
CAP2 = 2


def _mu(n):
    """MSCCC-style rank rate: per-class concurrency caps inside a CSRV pool."""
    return float(min(CSRV, min(n[0], CAP1) + min(n[1], CAP2)))


def _model():
    net = Network('oi_dispatch')
    d = Delay(net, 'Think')
    q = Queue(net, 'AGG', SchedStrategy.OI)
    c1 = ClosedClass(net, 'C1', 3, d, 0)
    c2 = ClosedClass(net, 'C2', 3, d, 0)
    d.set_service(c1, Exp(1.0))
    d.set_service(c2, Exp(0.5))

    def ordrate(c):
        # NOTE: python-native hands svcRateFun 0-BASED class indices
        # (MATLAB hands 1-based). See _kb/07-cross-language-parity.md.
        c = np.atleast_1d(np.asarray(c, dtype=int)).ravel()
        return _mu([int(np.sum(c == 0)), int(np.sum(c == 1))])

    q.set_service(ordrate)
    q.set_number_of_servers(1)
    q.set_cap(6)
    P = net.init_routing_matrix()
    for r in range(2):
        P[r][r] = net.serial_routing([d, q])
    net.link(P)
    return net


@pytest.fixture(scope='module')
def reference():
    return np.asarray(SolverCTMC(_model(), verbose=False).getAvgTput())[1, :]


@pytest.mark.parametrize('method', ['default', 'exact'])
def test_nc_reaches_the_exact_oi_analyzer(reference, method):
    X = np.asarray(SolverNC(_model(), method=method, verbose=False).getAvgTput())[1, :]
    assert np.allclose(X, reference, atol=1e-10), (method, X, reference)


def test_mva_agrees(reference):
    X = np.asarray(SolverMVA(_model(), method='exact', verbose=False).getAvgTput())[1, :]
    assert np.allclose(X, reference, atol=1e-10), (X, reference)


def test_comom_is_refused_rather_than_silently_wrong():
    with pytest.raises(Exception) as exc:
        SolverNC(_model(), method='comom', verbose=False).getAvgTput()
    assert 'order-independent' in str(exc.value), str(exc.value)
