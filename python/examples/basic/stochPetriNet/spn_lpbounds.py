"""Bound a stochastic Petri net by linear programming.

SolverBA's 'spnlp' family is the first BOUNDING route LINE offers for a Petri
net. SolverCTMC builds the explicit generator, SolverSSA and SolverLDES
simulate, SolverFLD fluidises, and SolverNC 'rec' needs a product form; this one
needs none of that. It relaxes the stationary chain to a MOMENT POLYTOPE -- the
uniformized evolution equation written for E[X_p], E[X_p^2] and E[X_p1 X_p2],
plus behavioural and probabilistic inequalities -- and then minimises and
maximises each reported measure over it. Every stationary point of the true
chain satisfies every row, so the two optima bracket the exact value.

    spn_lpbnd        assembles the polytope and solves the LPs
    spn_sinvariants  supplies the weighted place invariants, which are both the
                     conservation equalities and the a priori per-place caps

Reference: Z. Liu, "Performance Analysis of Stochastic Timed Petri Nets Using
Linear Programming Approach", IEEE Trans. Software Engineering 24(11), 1998,
1014-1030.
"""

from line_solver import (ClosedClass, Exp, Network, Place, SolverBA, SolverCTMC,
                         Transition)
from line_solver.api.spn import spn_lpbnd


def prodline(mu):
    """The reference's Fig. 2b: four servers, blocking before service.

    Server i cannot start until the downstream buffer has a free slot. The
    buffers hold 3, 2 and 4, and each is a conserved pair of places --
    (p5,p2), (p4,p1), (p3,p0) -- so the net is a strongly connected marked graph
    and all four transitions carry the same throughput.
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


def inhibiting(n=4):
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


# ---------------------------------------------------------------------------
# The reference's own Fig. 2b, and its Table 2.
#
# Table 2 reports four bound columns on five rate vectors. All four are
# reproduced below to the three decimals it prints. The one column not
# reproduced is its u.b.1, which is the upper side further tightened by the
# subnet-throughput theorems (its Thms 1 and 2); those are not implemented, so
# u.b.2 is the column to compare against.
MUS = [(1, 1.25, 2, 0.5), (1, 1.25, 2, 2.5), (1, 1.25, 1.25, 2.5),
       (1, 1.25, 1.25, 1), (1.111, 1.111, 1.111, 1.111)]
PUB = [(1.165, 1.951, 2.000, 0.930, 2.000),
       (1.829, 2.978, 3.529, 1.481, 4.000),
       (1.581, 2.873, 3.333, 1.333, 4.000),
       (1.359, 2.757, 3.333, 1.111, 4.000),
       (1.350, 2.667, 2.963, 1.111, 4.444)]

print('\nLiu (1998) Table 2: total throughput of the production line')
print('%-5s %-31s %-19s %-19s' % ('case', 'Markovian LP', 'published', 'operational LP'))
print('%-5s %9s %9s %9s   %9s %9s   %9s %9s'
      % ('', 'lower', 'simul', 'upper', 'l.b.', 'u.b.2', 'o.l.b.', 'o.u.b.'))
for c, mu in enumerate(MUS):
    sn = prodline(mu).getStruct()
    # The liveness rows of the reference's Table 1 are OPT-IN, because they hold
    # only on a live net and spn_lpbnd cannot certify liveness. This one is
    # live: a strongly connected marked graph with a token on every cycle. They
    # are the whole of the lower side, so the published l.b. needs them.
    lo = spn_lpbnd(sn, {'markovian': True, 'assumelive': True})
    up = spn_lpbnd(sn, {'markovian': True})
    op = spn_lpbnd(sn, {'markovian': False, 'assumelive': True})
    print('%-5d %9.4f %9.4f %9.4f   %9.3f %9.3f   %9.4f %9.4f'
          % (c + 1, lo['modeTput'][0].sum(), PUB[c][1], up['modeTput'][1].sum(),
             PUB[c][0], PUB[c][2], op['modeTput'][0].sum(), op['modeTput'][1].sum()))

# ---------------------------------------------------------------------------
# Through SolverBA, on a net with an inhibitor arc.
#
# The four method names are spnlp.upper, spnlp.lower and their spnlp.op.*
# counterparts, which drop the second-moment, covariance and Little's-law
# families and so need only a mean firing time rather than an exponential one.
# They are the only family SolverBA offers on a Petri net, and the only one it
# withholds off a Petri net: every other family is parameterized by demands and
# a population, which a marking is not.
print('\nmethods offered on this net: %s'
      % ', '.join(SolverBA(inhibiting()).list_valid_methods()))

ex = SolverCTMC(inhibiting()).getAvgTable()
lo = SolverBA(inhibiting(), 'spnlp.lower').getAvgTable()
up = SolverBA(inhibiting(), 'spnlp.upper').getAvgTable()
print('\nMean tokens per place, exact between the two sides:')
print('%-8s %10s %10s %10s' % ('place', 'lower', 'exact', 'upper'))
for i in range(len(ex)):
    print('%-8s %10.5f %10.5f %10.5f'
          % (ex.Station.iloc[i], lo.QLen.iloc[i], ex.QLen.iloc[i], up.QLen.iloc[i]))

# U = Q at a Place, which LINE models as an INF station -- the same convention
# SolverCTMC and SolverNC report. The paper's place utilization 1 - P(m = 0) is
# a different quantity and is not this column.
print('\nsolver.citations():')
SolverBA(inhibiting(), 'spnlp.upper').citations(display=True)
