"""
spn_pf: derive the product form of a stochastic Petri net, and the SolverNC
route that consumes it.

The MDD-rec paper (Balsamo-Marin-Stojic, FGCS 111 (2020) 475-490) takes the g_l
as GIVEN -- deriving them is declared out of scope in its Sec. 3.2 -- so
everything the api/spn functions do downstream was, until spn_pf, unreachable
from a solver. These tests pin the derivation itself.

FOUR INDEPENDENT ORACLES, because a product form that is merely plausible is
worse than none:

 1. GLOBAL BALANCE. pi built from the derived g_l must satisfy pi Q = 0 on the
    generator assembled independently from the net's own rate law. This is the
    definition of stationarity and shares no code path with the derivation.
 2. SolverCTMC. The measures must reproduce the explicit-generator solve.
 3. THE CERTIFICATE. Deficiency, linkage classes, stoichiometric rank and weak
    reversibility are structural facts about the net that can be read off by
    hand on these small examples.
 4. THE OTHER CODEBASES. y, G and the token counts are pinned at 12 decimals.

THE MODELS are the 3-place cyclic net at N = 4 with rates {1, 1.5, 2}, whose
Gordon-Newell factors are known in closed form and whose token counts are the
value the four codebases already agree on for spn_rec; the same net under
INFINITE-SERVER firing, which is mass action rather than a constant rate and so
selects the other psi; and a FORK-JOIN net, whose marking is not a conserved job
population at all -- a mode consumes one token and produces two -- which is the
case the MDD-rec paper exists to serve and which SolverCTMC cannot solve.
"""

from math import factorial

import numpy as np
import pytest

from line_solver import ClosedClass, Exp, Network, Place, SolverCTMC, SolverNC, Transition
from line_solver.api.spn import spn_metrics, spn_pf

RATES = [1.0, 1.5, 2.0]
NJOBS = 4

# MATLAB, the JAR and the C++ port on these nets, at %.12f.
CYCLIC_Y = [1.442249570307, 0.961499713538, 0.721124785154]
CYCLIC_G = 19.934426352559
CYCLIC_TOKENS = [2.249874392899, 1.069837548149, 0.680288058952]
CYCLIC_MODEX = [0.869201138838, 0.869201138838, 0.869201138838]

FJ_Y = [1.028383350947, 1.381974961646, 1.381974961646, 0.703630713806]
FJ_G = 20.320471787308
FJ_TOKENS = [0.710262429604, 1.864198441415, 1.864198441415, 0.425539128981]
FJ_MODEX = [0.607360882508, 0.607360882508, 0.607360882508]


def cyclic_spn(ntokens=NJOBS, servers=1):
    """P0 -> T0 -> P1 -> T1 -> P2 -> T2 -> P0, one token class.

    servers=1 makes every mode fire at its rate constant, which is psi = 1;
    servers=inf makes it fire at rate lambda*m_p, which is mass action.
    """
    model = Network('spn')
    places = [Place(model, 'P%d' % i) for i in range(3)]
    trans = [Transition(model, 'T%d' % i) for i in range(3)]
    cls = ClosedClass(model, 'Class1', ntokens, places[0])
    for i in range(3):
        mode = trans[i].add_mode('fire')
        trans[i].set_distribution(mode, Exp(RATES[i]))
        trans[i].set_enabling_conditions(mode, cls, places[i], 1)
        trans[i].set_firing_outcome(mode, cls, places[(i + 1) % 3], 1)
        if servers != 1:
            trans[i].set_number_of_servers(mode, servers)
    R = model.init_routing_matrix()
    for i in range(3):
        R.set(cls, cls, places[i], trans[i], 1.0)
        R.set(cls, cls, trans[i], places[(i + 1) % 3], 1.0)
    model.link(R)
    for i, v in enumerate([ntokens, 0, 0]):
        places[i].set_state([v])
    return model


def forkjoin_spn(ntokens=3, lf=1.3, lj=0.7, lb=1.9):
    """P0 -(Tf)-> P1 + P2 -(Tj)-> P3 -(Tb)-> P0.

    Tf consumes ONE token and produces TWO, Tj the reverse, so the marking is
    not a conserved population and the net has no queueing-network counterpart.
    Its place invariant is (2, 1, 1, 2).
    """
    model = Network('fj')
    P = [Place(model, 'P%d' % i) for i in range(4)]
    Tf, Tj, Tb = Transition(model, 'Tf'), Transition(model, 'Tj'), Transition(model, 'Tb')
    cls = ClosedClass(model, 'C', ntokens, P[0])
    m = Tf.add_mode('f')
    Tf.set_distribution(m, Exp(lf))
    Tf.set_enabling_conditions(m, cls, P[0], 1)
    Tf.set_firing_outcome(m, cls, P[1], 1)
    Tf.set_firing_outcome(m, cls, P[2], 1)
    m = Tj.add_mode('j')
    Tj.set_distribution(m, Exp(lj))
    Tj.set_enabling_conditions(m, cls, P[1], 1)
    Tj.set_enabling_conditions(m, cls, P[2], 1)
    Tj.set_firing_outcome(m, cls, P[3], 1)
    m = Tb.add_mode('b')
    Tb.set_distribution(m, Exp(lb))
    Tb.set_enabling_conditions(m, cls, P[3], 1)
    Tb.set_firing_outcome(m, cls, P[0], 1)
    R = model.init_routing_matrix()
    R.set(cls, cls, P[0], Tf, 1.0)
    R.set(cls, cls, Tf, P[1], 1.0)
    R.set(cls, cls, Tf, P[2], 1.0)
    R.set(cls, cls, P[1], Tj, 1.0)
    R.set(cls, cls, P[2], Tj, 1.0)
    R.set(cls, cls, Tj, P[3], 1.0)
    R.set(cls, cls, P[3], Tb, 1.0)
    R.set(cls, cls, Tb, P[0], 1.0)
    model.link(R)
    for i, v in enumerate([ntokens, 0, 0, 0]):
        P[i].set_state([v])
    return model


def _balance_residual(pf):
    """max |pi Q| over the enumerated reachable set.

    The generator is assembled here from the net's own rate law, independently
    of the derivation, so agreement is evidence and not a tautology.
    """
    info = pf['info']
    L = int(info['nplacelevels'])
    md = info['modes']
    st = info['mdd'].enumerate()[:, :L].astype(float)
    idx = {tuple(r): i for i, r in enumerate(st)}
    n = len(st)
    Q = np.zeros((n, n))
    for i, m in enumerate(st):
        for e in range(len(md)):
            enab = np.asarray(md[e]['enab'])
            fire = np.asarray(md[e]['fire'])
            if np.any(m < enab):
                continue
            rate = float(np.asarray(md[e]['D1']).ravel()[0])
            if pf['kind'] == 'massaction':
                for l in np.nonzero(enab > 0)[0]:
                    for j in range(int(enab[l])):
                        rate *= (m[l] - j)
            else:
                deg = min([np.floor(m[l] / enab[l]) for l in np.nonzero(enab > 0)[0]] or [1.0])
                rate *= min(deg, md[e]['srv'])
            t = tuple(m - enab + fire)
            assert t in idx, 'firing leaves the reachable set: %s' % (t,)
            Q[i, idx[t]] += rate
    for i in range(n):
        Q[i, i] -= Q[i].sum()
    pi = np.array([np.prod([pf['g'][l][int(m[l])] for l in range(L)]) for m in st])
    pi /= pi.sum()
    return float(np.max(np.abs(pi @ Q)))


def test_cyclic_product_form_is_the_gordon_newell_one():
    pf = spn_pf(cyclic_spn())
    assert pf['kind'] == 'geometric'
    # y is proportional to 1/mu: the closed-form Gordon-Newell factors, up to the
    # gauge the minimum-norm solution fixes.
    ratios = np.asarray(pf['y']) / pf['y'][0]
    np.testing.assert_allclose(ratios, [1.0, RATES[0] / RATES[1], RATES[0] / RATES[2]],
                               rtol=1e-12)
    np.testing.assert_allclose(pf['y'], CYCLIC_Y, rtol=1e-11)


def test_cyclic_certificate_is_the_structural_one():
    pf = spn_pf(cyclic_spn())
    # three complexes {e0, e1, e2}, one linkage class (the cycle is strongly
    # connected), stoichiometric rank 2, so deficiency 3 - 1 - 2 = 0
    assert pf['complexes'].shape[0] == 3
    assert pf['linkage'] == 1
    assert pf['srank'] == 2
    assert pf['deficiency'] == 0
    assert pf['weaklyreversible']
    assert pf['residual'] < 1e-12


def test_the_derived_law_satisfies_global_balance():
    for mk in (cyclic_spn, forkjoin_spn, lambda: cyclic_spn(servers=np.inf)):
        assert _balance_residual(spn_pf(mk())) < 1e-12


def test_infinite_server_firing_selects_mass_action():
    pf = spn_pf(cyclic_spn(servers=np.inf))
    assert pf['kind'] == 'massaction'
    # the factors carry the 1/k! that psi = prod 1/m! puts there
    for l in range(3):
        k = np.arange(len(pf['g'][l]))
        np.testing.assert_allclose(
            pf['g'][l], pf['y'][l] ** k / np.array([float(factorial(int(j))) for j in k]),
            rtol=1e-12)


def test_measures_reproduce_the_explicit_generator():
    pf = spn_pf(cyclic_spn())
    met = spn_metrics(pf['mdds'], pf['g'], pf['info'])
    ctmc = SolverCTMC(cyclic_spn())
    np.testing.assert_allclose(met['tokens'], np.ravel(ctmc.get_avg_qlen()), atol=1e-12)
    np.testing.assert_allclose(met['placeTput'], np.ravel(ctmc.get_avg_tput()), atol=1e-12)


def test_cyclic_values_are_pinned_across_the_codebases():
    pf = spn_pf(cyclic_spn())
    met = spn_metrics(pf['mdds'], pf['g'], pf['info'])
    assert met['G'] == pytest.approx(CYCLIC_G, rel=1e-11)
    np.testing.assert_allclose(met['tokens'], CYCLIC_TOKENS, rtol=1e-11)
    np.testing.assert_allclose(met['modeTput'], CYCLIC_MODEX, rtol=1e-11)


def test_forkjoin_values_are_pinned_across_the_codebases():
    pf = spn_pf(forkjoin_spn())
    met = spn_metrics(pf['mdds'], pf['g'], pf['info'])
    assert pf['deficiency'] == 0 and pf['weaklyreversible']
    np.testing.assert_allclose(pf['y'], FJ_Y, rtol=1e-11)
    assert met['G'] == pytest.approx(FJ_G, rel=1e-11)
    np.testing.assert_allclose(met['tokens'], FJ_TOKENS, rtol=1e-11)
    np.testing.assert_allclose(met['modeTput'], FJ_MODEX, rtol=1e-11)


def test_the_fork_and_the_join_carry_the_same_flow():
    # every token that forks must later join and return, so on this cycle the
    # three modes share one throughput -- a conservation law the derivation was
    # never told about
    met = spn_metrics(*(lambda pf: (pf['mdds'], pf['g'], pf['info']))(spn_pf(forkjoin_spn())))
    np.testing.assert_allclose(met['modeTput'], met['modeTput'][0] * np.ones(3), rtol=1e-11)


def test_solver_nc_solves_a_petri_net():
    nc = SolverNC(cyclic_spn())
    ctmc = SolverCTMC(cyclic_spn())
    np.testing.assert_allclose(np.ravel(nc.get_avg_qlen()), np.ravel(ctmc.get_avg_qlen()),
                               atol=1e-12)
    np.testing.assert_allclose(np.ravel(nc.get_avg_tput()), np.ravel(ctmc.get_avg_tput()),
                               atol=1e-12)
    # UN follows LINE's INF-station convention, U = Q, not the paper's
    # 1 - P(m = 0); the latter rides on the certificate
    np.testing.assert_allclose(np.ravel(nc.get_avg_util()), np.ravel(ctmc.get_avg_util()),
                               atol=1e-12)


def test_solver_ctmc_mdd_agrees_on_a_product_form_net():
    # the level aggregation is EXACT on a product-form model (paper Sec. 5), so
    # here it must reproduce the explicit generator rather than approximate it
    ex = SolverCTMC(cyclic_spn())
    md = SolverCTMC(cyclic_spn(), method='mdd')
    np.testing.assert_allclose(np.ravel(md.get_avg_qlen()), np.ravel(ex.get_avg_qlen()),
                               atol=1e-9)
    np.testing.assert_allclose(np.ravel(md.get_avg_tput()), np.ravel(ex.get_avg_tput()),
                               atol=1e-9)


def test_an_inhibitor_arc_is_refused_by_name():
    model = cyclic_spn()
    trans = [n for n in model._nodes if isinstance(n, Transition)]
    cls = model.get_classes()[0]
    places = [n for n in model._nodes if isinstance(n, Place)]
    trans[0].set_inhibiting_conditions(trans[0].get_modes()[0], cls, places[2], 2)
    with pytest.raises(Exception, match='inhibitor'):
        spn_pf(model)
