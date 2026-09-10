"""Routing reducibility and its repair, against the MATLAB reference.

MATLAB `@MNetwork/isRoutingErgodic`, `getReducibilityInfo`,
`getAbsorbingStations` and `makeErgodic` on the fixture D -> Q1 -> Q2 -> Q2:
Q2 absorbs, D and Q1 are transient, and there are three strongly connected
components.
"""

import warnings

import numpy as np

from line_solver import (ClosedClass, Delay, Exp, Network, Queue, SchedStrategy)


def _reducible():
    m = Network('red')
    d = Delay(m, 'D')
    q1 = Queue(m, 'Q1', SchedStrategy.PS)
    q2 = Queue(m, 'Q2', SchedStrategy.PS)
    c = ClosedClass(m, 'C', 2, d, 0)
    d.setService(c, Exp(1.0))
    q1.setService(c, Exp(2.0))
    q2.setService(c, Exp(3.0))
    P = m.initRoutingMatrix()
    M = np.zeros((3, 3))
    M[0, 1] = 1.0
    M[1, 2] = 1.0
    M[2, 2] = 1.0
    P.set(c, c, M)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        m.link(P)
    return m


def _ergodic():
    m = Network('erg')
    d = Delay(m, 'D')
    q = Queue(m, 'Q', SchedStrategy.PS)
    c = ClosedClass(m, 'C', 2, d, 0)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    m.link(Network.serialRouting(d, q))
    return m


def test_reducible_structure_matches_matlab():
    is_erg, info = _reducible().is_routing_ergodic()
    assert is_erg is False
    assert info['isReducible'] is True
    assert info['numSCCs'] == 3
    assert info['absorbingStations'] == ['Q2']
    assert sorted(info['transientStations']) == ['D', 'Q1']


def test_reducibility_info_suggests_one_fix_per_absorbing_station():
    info = _reducible().get_reducibility_info()
    assert info['isRoutingErgodic'] is False
    assert info['suggestedFixes'] == [
        'Route jobs from Q2 back to D (e.g., P{class}(Q2, D) = 1.0)']


def test_absorbing_stations_are_returned_as_objects():
    stations, idxs = _reducible().get_absorbing_stations()
    assert [s.name for s in stations] == ['Q2']
    # get_node_index is 1-based, as in MATLAB
    assert idxs == [3]


def test_make_ergodic_redirects_the_absorbing_row():
    P = _reducible().make_ergodic()
    # MATLAB returns [0 1 0; 0 0 1; 1 0 0]
    want = np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]])
    assert np.allclose(P[0][0], want)


def test_an_ergodic_routing_is_left_alone():
    m = _ergodic()
    assert m.is_routing_ergodic()[0] is True
    info = m.get_reducibility_info()
    assert info['isRoutingErgodic'] is True
    assert info['suggestedFixes'] == []
    assert m.get_absorbing_stations()[0] == []
