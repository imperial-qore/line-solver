"""
MDD-rec and the SPN measures built on it: mdd_rec, mdd_rec_marginal,
spn_rec_enabled, spn_metrics, spn_sinvariants, spn_conv.

S. Balsamo, A. Marin, I. Stojic, "Computation of the normalising constant for
product-form models of distributed systems with synchronisation", FGCS 111
(2020) 475-490.

THREE INDEPENDENT ORACLES, because a normalising constant is a single number
that a wrong recursion can still produce plausibly:

 1. EXPLICIT SUM. G is also the sum over the enumerated reachable set of the
    product of the g_l. The diagram walk and the enumeration share no code path,
    so agreement pins the recursion itself.
 2. THE CONVOLUTION. spn_conv decomposes {m : S m = V} instead of walking the
    diagram; on this net the two must agree exactly (FGCS Sec. 5.2).
 3. THE OTHER CODEBASES. The marginals must reproduce MATLAB's and the C++
    port's queue lengths on the same net, and the mode throughputs must match
    the level aggregation's X, which comes from a different algorithm on a
    different descriptor (mdd_mcd on the QN form of the same net).

THE MODEL is the 3-place cyclic net at N = 4 with firing rates {1, 1.5, 2}: a
Gordon-Newell network, so its product form is known in closed form,
g_l(n) = (1/mu_l)^n with unit visit ratios, and no product-form TEST is needed to
obtain the g_l -- which is the part the paper itself declares out of scope.
"""

import numpy as np
import pytest

from line_solver import ClosedClass, Exp, Network, Place, Transition
from line_solver.api.mdd import (mdd_descriptor, mdd_mcd, mdd_reachset, mdd_rec,
                                 mdd_rec_marginal, mdd_rec_masked)
from line_solver.api.spn import (spn_conv, spn_mdd, spn_metrics, spn_rec_enabled,
                                 spn_sinvariants)

RATES = [1.0, 1.5, 2.0]
NJOBS = 4
# MATLAB and the C++ port on this net, at %.12f.
QLEN_REF = [2.249874392899, 1.069837548149, 0.680288058952]


def cyclic_spn(ntokens=NJOBS):
    """P0 -> T0 -> P1 -> T1 -> P2 -> T2 -> P0, one token class."""
    model = Network('spn')
    places = [Place(model, 'P%d' % i) for i in range(3)]
    trans = [Transition(model, 'T%d' % i) for i in range(3)]
    cls = ClosedClass(model, 'Class1', ntokens, places[0])
    for i in range(3):
        mode = trans[i].add_mode('fire')
        trans[i].set_distribution(mode, Exp(RATES[i]))
        trans[i].set_enabling_conditions(mode, cls, places[i], 1)
        trans[i].set_firing_outcome(mode, cls, places[(i + 1) % 3], 1)
    R = model.init_routing_matrix()
    for i in range(3):
        R.set(cls, cls, places[i], trans[i], 1.0)
        R.set(cls, cls, trans[i], places[(i + 1) % 3], 1.0)
    model.link(R)
    places[0].set_state([ntokens])
    places[1].set_state([0])
    places[2].set_state([0])
    return model


def gordon_newell_g(domain):
    """g_l(n) = (v_l/mu_l)^n, unit visit ratios on a cycle."""
    return [[RATES[l] ** (-n) for n in range(int(domain[l]))] for l in range(len(domain))]


@pytest.fixture(scope='module')
def net():
    mdds, desc, info = spn_mdd(cyclic_spn())
    return mdds, desc, info, gordon_newell_g(mdds.domain)


def test_mdd_rec_matches_the_explicit_sum(net):
    mdds, _, info, g = net
    states = info['mdd'].enumerate()
    assert len(states) == 15                       # C(4+2,2), the closed lattice
    expected = sum(float(np.prod([g[l][s[l]] for l in range(int(mdds.K))])) for s in states)
    assert mdd_rec(mdds, g) == pytest.approx(expected, rel=1e-12)


def test_an_empty_mask_is_the_unmasked_recursion(net):
    mdds, _, _, g = net
    mask = [np.ones(int(mdds.domain[l]), dtype=bool) for l in range(int(mdds.K))]
    assert mdd_rec_masked(mdds, g, mask) == pytest.approx(mdd_rec(mdds, g), rel=1e-12)


def test_the_marginals_partition_the_normalising_constant(net):
    mdds, _, info, g = net
    G = mdd_rec(mdds, g)
    for l in range(int(info['nplacelevels'])):
        mass = mdd_rec_marginal(mdds, g, l)
        assert mass.sum() == pytest.approx(G, rel=1e-12)


def test_sinvariants_find_the_token_conservation_law(net):
    _, _, _, _ = net
    inv = spn_sinvariants(cyclic_spn().getStruct())
    # the only minimal-support invariant of a cycle is "the tokens are conserved"
    assert inv['S'] == [[1, 1, 1]]
    assert inv['V'] == [NJOBS]
    assert inv['m0'] == [NJOBS, 0, 0]


def test_the_convolution_agrees_with_the_diagram_walk(net):
    mdds, _, _, g = net
    inv = spn_sinvariants(cyclic_spn().getStruct())
    assert spn_conv(inv, g) == pytest.approx(mdd_rec(mdds, g), rel=1e-12)


def test_metrics_reproduce_the_matlab_queue_lengths(net):
    mdds, _, info, g = net
    met = spn_metrics(mdds, g, info)
    assert met['tokens'] == pytest.approx(QLEN_REF, abs=1e-9)
    assert sum(met['tokens']) == pytest.approx(NJOBS, abs=1e-12)
    for l in range(int(info['nplacelevels'])):
        assert met['placeUtil'][l] == pytest.approx(1.0 - met['marginal'][l][0], rel=1e-12)


def test_mode_throughput_matches_the_level_aggregation(net):
    mdds, _, info, g = net
    met = spn_metrics(mdds, g, info)
    P = np.zeros((3, 3))
    for i in range(3):
        P[i, (i + 1) % 3] = 1.0
    qdesc = mdd_descriptor(RATES, P, [1.0, 1.0, 1.0], NJOBS)
    qdiag = mdd_reachset(qdesc['domain'], qdesc['init'], qdesc['nextfun'])
    X = np.ravel(np.asarray(mdd_mcd(qdiag.to_struct(), qdesc)['X'], dtype=float))
    # a cycle carries one flow, so every transition sees the same throughput
    assert met['modeTput'] == pytest.approx(X[0], rel=1e-6)
    assert met['placeTput'] == pytest.approx(X[0], rel=1e-6)


def test_enabling_degrees_are_a_distribution(net):
    mdds, _, info, g = net
    G = mdd_rec(mdds, g)
    met = spn_metrics(mdds, g, info)
    for e, mde in enumerate(info['modes']):
        en = spn_rec_enabled(mdds, g, mde, int(info['nplacelevels']))
        assert en['maxDegree'] == NJOBS               # one token per firing set
        assert en['ge'][0] == pytest.approx(G, rel=1e-12)
        assert sum(en['eq']) == pytest.approx(G, rel=1e-12)
        assert en['ge'][1] / G == pytest.approx(met['modeUtil'][e], rel=1e-12)
        # P(e >= k) is non-increasing by construction; a mask that grew the set
        # would be a silently wrong measure rather than an error
        assert np.all(np.diff(en['ge']) <= 1e-15)


def test_a_mode_with_no_input_place_is_refused(net):
    mdds, _, info, g = net
    mde = dict(info['modes'][0])
    mde['enab'] = np.zeros_like(np.asarray(mde['enab'], dtype=float))
    with pytest.raises(Exception):
        spn_rec_enabled(mdds, g, mde, int(info['nplacelevels']))


def test_a_shape_mismatch_in_g_is_refused(net):
    mdds, _, _, g = net
    with pytest.raises(Exception):
        mdd_rec(mdds, g[:-1])
