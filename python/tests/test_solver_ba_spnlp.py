"""SolverBA 'spnlp.*': moment-relaxation LP bounds for stochastic Petri nets.

THE ACCEPTANCE TEST IS THE REFERENCE'S OWN TABLE 2. Liu (1998) publishes four
bound columns for the four-server production line of its Fig. 2b, on five rate
vectors, and this suite asserts all four to three decimals. Those numbers are
what a transcription error in any constraint family would move: the polytope has
twelve of them and the published optimum is a function of the whole set, so a
bracket check alone would not catch a family that was dropped or mis-signed.

The CTMC checks are the other half. A bound that agreed with the paper and still
failed to contain the exact answer would be wrong whatever it agreed with.

Reference: Z. Liu, "Performance Analysis of Stochastic Timed Petri Nets Using
Linear Programming Approach", IEEE Trans. Software Engineering 24(11), 1998,
1014-1030.
"""

import numpy as np
import pytest

from line_solver import (ClosedClass, Erlang, Exp, Network, Pareto, Place,
                         SolverBA, SolverCTMC, SolverNC, Transition)
from line_solver.api.spn import spn_lpbnd


# --------------------------------------------------------------------------
def liu98_line(mu):
    """Fig. 2b: four servers, blocking before service, buffers of 3, 2 and 4.

    (p5, p2), (p4, p1) and (p3, p0) are the three buffer pairs, each conserved
    at its capacity, so the net is a strongly connected marked graph and every
    transition carries the same throughput.
    """
    model = Network('liu98')
    p5, p4, p3 = Place(model, 'p5'), Place(model, 'p4'), Place(model, 'p3')
    p2, p1, p0 = Place(model, 'p2'), Place(model, 'p1'), Place(model, 'p0')
    t1, t2 = Transition(model, 't1'), Transition(model, 't2')
    t3, t4 = Transition(model, 't3'), Transition(model, 't4')
    jc = ClosedClass(model, 'Class1', 9, p2, 0)
    m = t1.addMode('m1')
    t1.setDistribution(m, Exp(mu[0]))
    t1.setEnablingConditions(m, jc, p2, 1)
    t1.setFiringOutcome(m, jc, p5, 1)
    m = t2.addMode('m2')
    t2.setDistribution(m, Exp(mu[1]))
    t2.setEnablingConditions(m, jc, p5, 1)
    t2.setEnablingConditions(m, jc, p1, 1)
    t2.setFiringOutcome(m, jc, p4, 1)
    t2.setFiringOutcome(m, jc, p2, 1)
    m = t3.addMode('m3')
    t3.setDistribution(m, Exp(mu[2]))
    t3.setEnablingConditions(m, jc, p4, 1)
    t3.setEnablingConditions(m, jc, p0, 1)
    t3.setFiringOutcome(m, jc, p3, 1)
    t3.setFiringOutcome(m, jc, p1, 1)
    m = t4.addMode('m4')
    t4.setDistribution(m, Exp(mu[3]))
    t4.setEnablingConditions(m, jc, p3, 1)
    t4.setFiringOutcome(m, jc, p0, 1)
    R = model.initRoutingMatrix()
    for a, b in [(p2, t1), (p5, t2), (p1, t2), (p4, t3), (p0, t3), (p3, t4)]:
        R.set(jc, jc, a, b, 1.0)
    for a, b in [(t1, p5), (t2, p4), (t2, p2), (t3, p3), (t3, p1), (t4, p0)]:
        R.set(jc, jc, a, b, 1.0)
    model.link(R)
    for pl, v in [(p5, 0), (p4, 0), (p3, 0), (p2, 3), (p1, 2), (p0, 4)]:
        pl.setState(v)
    return model


def inhibiting_spn(n=4):
    """Three places, four modes, one inhibitor arc; the token count is conserved."""
    model = Network('spn')
    p1, p2, p3 = Place(model, 'P1'), Place(model, 'P2'), Place(model, 'P3')
    t1, t2, t3 = Transition(model, 'T1'), Transition(model, 'T2'), Transition(model, 'T3')
    jc = ClosedClass(model, 'Class1', n, p1, 0)
    m = t1.addMode('Mode1')
    t1.setDistribution(m, Exp(2))
    t1.setEnablingConditions(m, jc, p1, 2)
    t1.setFiringOutcome(m, jc, p2, 2)
    m = t1.addMode('Mode2')
    t1.setDistribution(m, Exp(1))
    t1.setEnablingConditions(m, jc, p1, 1)
    t1.setFiringOutcome(m, jc, p3, 1)
    m = t2.addMode('Mode3')
    t2.setDistribution(m, Exp(4))
    t2.setEnablingConditions(m, jc, p2, 1)
    t2.setFiringOutcome(m, jc, p1, 1)
    m = t3.addMode('Mode4')
    t3.setDistribution(m, Exp(1))
    t3.setEnablingConditions(m, jc, p3, 3)
    t3.setInhibitingConditions(m, jc, p2, 1)
    t3.setFiringOutcome(m, jc, p1, 3)
    R = model.initRoutingMatrix()
    for a, b in [(p1, t1), (p2, t2), (p2, t3), (p3, t3),
                 (t1, p2), (t1, p3), (t2, p1), (t3, p1)]:
        R.set(jc, jc, a, b, 1.0)
    model.link(R)
    p1.setState(n)
    p2.setState(0)
    p3.setState(0)
    return model


def cyclic_spn(n=4, rates=(1.3, 0.7, 1.9)):
    """P0 -> T0 -> P1 -> T1 -> P2 -> T2 -> P0, one class."""
    model = Network('cyc')
    P = [Place(model, 'P%d' % i) for i in range(3)]
    T = [Transition(model, 'T%d' % i) for i in range(3)]
    jc = ClosedClass(model, 'Class1', n, P[0])
    for i in range(3):
        m = T[i].addMode('fire')
        T[i].setDistribution(m, Exp(rates[i]))
        T[i].setEnablingConditions(m, jc, P[i], 1)
        T[i].setFiringOutcome(m, jc, P[(i + 1) % 3], 1)
    R = model.initRoutingMatrix()
    for i in range(3):
        R.set(jc, jc, P[i], T[i], 1.0)
        R.set(jc, jc, T[i], P[(i + 1) % 3], 1.0)
    model.link(R)
    for i, v in enumerate([n, 0, 0]):
        P[i].setState(v)
    return model


def forkjoin_spn(n=3, lf=1.3, lj=0.7, lb=1.9):
    """P0 -(Tf)-> P1 + P2 -(Tj)-> P3 -(Tb)-> P0.

    Tf consumes one token and produces two, so the marking is not a conserved
    population; the place invariant is (2, 1, 1, 2). This is what exercises the
    WEIGHTED form of the invariant family -- the reference writes that family
    for unweighted cycles, and an implementation that only summed token counts
    would produce an unbounded polytope here.
    """
    model = Network('fj')
    P = [Place(model, 'P%d' % i) for i in range(4)]
    Tf, Tj, Tb = Transition(model, 'Tf'), Transition(model, 'Tj'), Transition(model, 'Tb')
    jc = ClosedClass(model, 'C', n, P[0])
    m = Tf.addMode('f')
    Tf.setDistribution(m, Exp(lf))
    Tf.setEnablingConditions(m, jc, P[0], 1)
    Tf.setFiringOutcome(m, jc, P[1], 1)
    Tf.setFiringOutcome(m, jc, P[2], 1)
    m = Tj.addMode('j')
    Tj.setDistribution(m, Exp(lj))
    Tj.setEnablingConditions(m, jc, P[1], 1)
    Tj.setEnablingConditions(m, jc, P[2], 1)
    Tj.setFiringOutcome(m, jc, P[3], 1)
    m = Tb.addMode('b')
    Tb.setDistribution(m, Exp(lb))
    Tb.setEnablingConditions(m, jc, P[3], 1)
    Tb.setFiringOutcome(m, jc, P[0], 1)
    R = model.initRoutingMatrix()
    for a, b in [(P[0], Tf), (Tf, P[1]), (Tf, P[2]), (P[1], Tj), (P[2], Tj),
                 (Tj, P[3]), (P[3], Tb), (Tb, P[0])]:
        R.set(jc, jc, a, b, 1.0)
    model.link(R)
    for i, v in enumerate([n, 0, 0, 0]):
        P[i].setState(v)
    return model


# --------------------------------------------------------------------------
# Liu (1998) Table 2, p. 1023: the total throughput of the production line.
# Columns, in order: mu, l.b., u.b.2, o.l.b., o.u.b. The published u.b.1 is the
# upper side TIGHTENED by the subnet-throughput theorems (Thms 1 and 2), which
# are deliberately not implemented, so u.b.2 is the column to match.
LIU_TABLE2 = [
    ((1.000, 1.25, 2.00, 0.50), 1.165, 2.000, 0.930, 2.000),
    ((1.000, 1.25, 2.00, 2.50), 1.829, 3.529, 1.481, 4.000),
    ((1.000, 1.25, 1.25, 2.50), 1.581, 3.333, 1.333, 4.000),
    ((1.000, 1.25, 1.25, 1.00), 1.359, 3.333, 1.111, 4.000),
    ((1.111, 1.111, 1.111, 1.111), 1.350, 2.963, 1.111, 4.444),
]


@pytest.mark.parametrize('mu,lb,ub,olb,oub', LIU_TABLE2)
def test_liu98_table2_upper(mu, lb, ub, olb, oub):
    """The Markovian upper side reproduces the published u.b.2 column."""
    bnd = spn_lpbnd(liu98_line(mu).getStruct(), {'markovian': True})
    assert float(bnd['modeTput'][1].sum()) == pytest.approx(ub, abs=5e-4)


@pytest.mark.parametrize('mu,lb,ub,olb,oub', LIU_TABLE2)
def test_liu98_table2_lower_needs_liveness(mu, lb, ub, olb, oub):
    """The published l.b. column is the lower side WITH the liveness rows.

    They are opt-in because they hold only on a live net, and this one is: a
    strongly connected marked graph with a token on every cycle. Without them
    the lower side falls back to the operational value, which is the second
    assertion here and is why the default is not a silent loss.
    """
    sn = liu98_line(mu).getStruct()
    live = spn_lpbnd(sn, {'markovian': True, 'assumelive': True})
    assert float(live['modeTput'][0].sum()) == pytest.approx(lb, abs=5e-4)
    plain = spn_lpbnd(sn, {'markovian': True})
    assert float(plain['modeTput'][0].sum()) <= float(live['modeTput'][0].sum()) + 1e-9


@pytest.mark.parametrize('mu,lb,ub,olb,oub', LIU_TABLE2)
def test_liu98_table2_operational(mu, lb, ub, olb, oub):
    """Both operational columns, which need no Markovian family at all."""
    bnd = spn_lpbnd(liu98_line(mu).getStruct(), {'markovian': False, 'assumelive': True})
    assert float(bnd['modeTput'][0].sum()) == pytest.approx(olb, abs=5e-4)
    assert float(bnd['modeTput'][1].sum()) == pytest.approx(oub, abs=5e-4)


# --------------------------------------------------------------------------
def _exact_by_place(table):
    """QLen and Tput keyed by station name.

    getAvgTable drops a (station, class) row whose metrics are all zero, so the
    table is not a full grid and must be joined by name rather than by position.
    """
    q, t = {}, {}
    for name, ql, tp in zip(table.Station, table.QLen, table.Tput):
        q[str(name)] = float(ql)
        t[str(name)] = float(tp)
    return q, t


# THE ORACLE IS PER NET, and the fork-join case is why. SolverCTMC returns a
# degenerate answer on forkjoin_spn -- all three tokens parked at P0 and zero
# throughput -- while SolverNC 'rec' solves it, which is also the pairing
# matlab/examples/basic/stochPetriNet/spn_productform_nc.m uses: it checks the
# cyclic net against CTMC and the fork-join one against NC alone. That is a
# pre-existing CTMC question and not one this bound can answer, so each net is
# checked against the solver that solves it.
@pytest.mark.parametrize('build,oracle', [(inhibiting_spn, SolverCTMC),
                                          (cyclic_spn, SolverCTMC),
                                          (forkjoin_spn, SolverNC)])
@pytest.mark.parametrize('markovian', [True, False])
def test_bracket_contains_exact(build, oracle, markovian):
    """lower <= exact <= upper, per place."""
    model = build()
    q, t = _exact_by_place(oracle(model).getAvgTable())
    bnd = spn_lpbnd(build().getStruct(), {'markovian': markovian})
    for l, name in enumerate(bnd['levelname']):
        place = name.split('.')[0]
        if place not in q:
            continue
        assert bnd['tokens'][0, l] <= q[place] + 1e-6
        assert q[place] <= bnd['tokens'][1, l] + 1e-6
        assert bnd['placeTput'][0, l] <= t[place] + 1e-6
        assert t[place] <= bnd['placeTput'][1, l] + 1e-6


@pytest.mark.parametrize('build', [inhibiting_spn, cyclic_spn, forkjoin_spn])
def test_markovian_is_at_least_as_tight_as_operational(build):
    """The Markovian polytope is the operational one plus rows, so it nests."""
    sn = build().getStruct()
    mk = spn_lpbnd(sn, {'markovian': True})
    op = spn_lpbnd(sn, {'markovian': False})
    for l in range(mk['nplacelevels']):
        assert op['tokens'][0, l] <= mk['tokens'][0, l] + 1e-9
        assert mk['tokens'][1, l] <= op['tokens'][1, l] + 1e-9


@pytest.mark.parametrize('build', [inhibiting_spn, cyclic_spn, forkjoin_spn])
def test_bracket_is_not_vacuous(build):
    """A [0, B] bracket is a missing constraint family, not a loose bound.

    The QRF lesson recorded in _kb/03-api-layer.md: never accept the whole box
    as "loose but valid", because it means the polytope does not constrain the
    objective at all.
    """
    bnd = spn_lpbnd(build().getStruct(), {'markovian': True})
    for l in range(bnd['nplacelevels']):
        assert np.isfinite(bnd['bound'][l])
        width = bnd['tokens'][1, l] - bnd['tokens'][0, l]
        assert width < bnd['bound'][l] - 1e-9


def test_degenerate_single_place_is_exact():
    """One place, one mode, a self-loop: nothing can move, so the bracket collapses."""
    model = Network('one')
    p = Place(model, 'P')
    t = Transition(model, 'T')
    jc = ClosedClass(model, 'C', 3, p)
    m = t.addMode('fire')
    t.setDistribution(m, Exp(2.0))
    t.setEnablingConditions(m, jc, p, 1)
    t.setFiringOutcome(m, jc, p, 1)
    R = model.initRoutingMatrix()
    R.set(jc, jc, p, t, 1.0)
    R.set(jc, jc, t, p, 1.0)
    model.link(R)
    p.setState(3)
    bnd = spn_lpbnd(model.getStruct(), {'markovian': True})
    assert bnd['tokens'][0, 0] == pytest.approx(3.0, abs=1e-9)
    assert bnd['tokens'][1, 0] == pytest.approx(3.0, abs=1e-9)


# --------------------------------------------------------------------------
def test_solver_ba_offers_only_spnlp_on_a_petri_net():
    """Every other family is parameterized by demands, which a marking is not."""
    valid = SolverBA(inhibiting_spn()).list_valid_methods()
    assert sorted(valid) == ['spnlp.lower', 'spnlp.op.lower', 'spnlp.op.upper', 'spnlp.upper']


def test_solver_ba_hides_spnlp_off_a_petri_net():
    from line_solver import Delay, Queue, SchedStrategy
    m = Network('cqn')
    d = Delay(m, 'Z')
    q = Queue(m, 'Q', SchedStrategy.PS)
    c = ClosedClass(m, 'C', 4, d, 0)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    m.link(Network.serialRouting(d, q))
    assert not [x for x in SolverBA(m).list_valid_methods() if x.startswith('spnlp')]


@pytest.mark.parametrize('method', ['spnlp.upper', 'spnlp.lower',
                                    'spnlp.op.upper', 'spnlp.op.lower'])
def test_solver_ba_end_to_end(method):
    """The four tokens run and land the reported side on the right place rows."""
    model = inhibiting_spn()
    s = SolverBA(model, method)
    s.runAnalyzer()
    r = s._result
    ex = np.asarray(SolverCTMC(inhibiting_spn()).getAvgTable().QLen, dtype=float)
    qn = np.asarray(r['QN'], dtype=float).ravel()
    assert qn.size == ex.size
    # U = Q at an INF station, which is what a Place is; see the analyzer note
    assert np.allclose(np.asarray(r['UN'], dtype=float).ravel(), qn)
    if method.endswith('upper'):
        assert np.all(qn >= ex - 1e-6)
    else:
        assert np.all(qn <= ex + 1e-6)


def test_solver_ba_citations_name_the_paper():
    s = SolverBA(inhibiting_spn(), 'spnlp.upper')
    s.runAnalyzer()
    text = ' '.join(str(c) for c in s.citations())
    assert 'Liu' in text and 'Petri' in text


# --------------------------------------------------------------------------
def test_phase_type_refused_by_markovian_accepted_by_operational():
    """An Erlang mode has no marking-only state, but it does have a mean."""
    model = cyclic_spn()
    # rebuild with one Erlang mode
    model = Network('cycph')
    P = [Place(model, 'P%d' % i) for i in range(3)]
    T = [Transition(model, 'T%d' % i) for i in range(3)]
    jc = ClosedClass(model, 'Class1', 4, P[0])
    for i in range(3):
        m = T[i].addMode('fire')
        T[i].setDistribution(m, Erlang.fitMeanAndOrder(1.0, 2) if i == 0 else Exp(1.0))
        T[i].setEnablingConditions(m, jc, P[i], 1)
        T[i].setFiringOutcome(m, jc, P[(i + 1) % 3], 1)
    R = model.initRoutingMatrix()
    for i in range(3):
        R.set(jc, jc, P[i], T[i], 1.0)
        R.set(jc, jc, T[i], P[(i + 1) % 3], 1.0)
    model.link(R)
    for i, v in enumerate([4, 0, 0]):
        P[i].setState(v)
    sn = model.getStruct()
    with pytest.raises(Exception, match='phase-type'):
        spn_lpbnd(sn, {'markovian': True})
    bnd = spn_lpbnd(sn, {'markovian': False})
    assert np.isfinite(bnd['tokens']).all()


def test_non_phase_type_refused_by_both():
    """A Pareto mode reaches sn as (shape, scale) with no mean to read."""
    model = Network('par')
    p = Place(model, 'P')
    t = Transition(model, 'T')
    jc = ClosedClass(model, 'C', 2, p)
    m = t.addMode('fire')
    t.setDistribution(m, Pareto(2.5, 1.0))
    t.setEnablingConditions(m, jc, p, 1)
    t.setFiringOutcome(m, jc, p, 1)
    R = model.initRoutingMatrix()
    R.set(jc, jc, p, t, 1.0)
    R.set(jc, jc, t, p, 1.0)
    model.link(R)
    p.setState(2)
    sn = model.getStruct()
    for markovian in (True, False):
        with pytest.raises(Exception, match='not phase-type'):
            spn_lpbnd(sn, {'markovian': markovian})


def test_infinite_server_mode_refused_by_name():
    model = Network('is')
    p = Place(model, 'P')
    t = Transition(model, 'T')
    jc = ClosedClass(model, 'C', 2, p)
    m = t.addMode('fire')
    t.setDistribution(m, Exp(1.0))
    t.setNumberOfServers(m, float('inf'))
    t.setEnablingConditions(m, jc, p, 1)
    t.setFiringOutcome(m, jc, p, 1)
    R = model.initRoutingMatrix()
    R.set(jc, jc, p, t, 1.0)
    R.set(jc, jc, t, p, 1.0)
    model.link(R)
    p.setState(2)
    with pytest.raises(Exception, match='servers'):
        spn_lpbnd(model.getStruct(), {'markovian': True})
