"""
Solve a stochastic Petri net analytically with SolverNC.

SolverNC's 'rec' method is the first ANALYTICAL route LINE offers for a Petri
net: SolverCTMC builds the explicit generator, SolverSSA and SolverLDES
simulate, SolverFLD fluidises. It works in three steps, each with its own
reference:

    spn_pf      decides whether the net has a product form and derives the
                per-place factors g_l, by complex balance
                (Coleman-Henderson-Taylor, Perform. Eval. 26(3), 1996)
    mdd_rec     evaluates G = sum over the reachable set of prod_l g_l(m_l) by
                ONE memoised walk of the decision diagram holding that set
                (Balsamo-Marin-Stojic, FGCS 111 (2020) 475-490)
    spn_metrics reads the mean tokens, the utilisations and the throughputs off
                masked walks of the same diagram

Two nets are solved here. The first is a closed cycle, which a queueing network
could also express. The second FORKS: its transition Tf consumes one token and
produces two, and Tj consumes two and produces one, so the marking is not a
conserved job population and there is no queueing-network counterpart -- which
is exactly the limitation the MDD-rec paper opens with.
"""

import numpy as np

from line_solver import ClosedClass, Exp, Network, Place, SolverCTMC, SolverNC, Transition
from line_solver.api.spn import spn_metrics, spn_pf


def cyclic(ntokens=4):
    rates = [1.0, 1.5, 2.0]
    model = Network('spn')
    places = [Place(model, 'P%d' % i) for i in range(3)]
    trans = [Transition(model, 'T%d' % i) for i in range(3)]
    cls = ClosedClass(model, 'Class1', ntokens, places[0])
    for i in range(3):
        mode = trans[i].add_mode('fire')
        trans[i].set_distribution(mode, Exp(rates[i]))
        trans[i].set_number_of_servers(mode, 1)
        trans[i].set_enabling_conditions(mode, cls, places[i], 1)
        trans[i].set_firing_outcome(mode, cls, places[(i + 1) % 3], 1)
    R = model.init_routing_matrix()
    for i in range(3):
        R.set(cls, cls, places[i], trans[i], 1.0)
        R.set(cls, cls, trans[i], places[(i + 1) % 3], 1.0)
    model.link(R)
    for i, v in enumerate([ntokens, 0, 0]):
        places[i].set_state([v])
    return model


def forkjoin(ntokens=3):
    model = Network('fj')
    P = [Place(model, 'P%d' % i) for i in range(4)]
    Tf, Tj, Tb = Transition(model, 'Tf'), Transition(model, 'Tj'), Transition(model, 'Tb')
    cls = ClosedClass(model, 'C', ntokens, P[0])
    m = Tf.add_mode('f')
    Tf.set_distribution(m, Exp(1.3))
    Tf.set_number_of_servers(m, 1)
    Tf.set_enabling_conditions(m, cls, P[0], 1)
    Tf.set_firing_outcome(m, cls, P[1], 1)
    Tf.set_firing_outcome(m, cls, P[2], 1)
    m = Tj.add_mode('j')
    Tj.set_distribution(m, Exp(0.7))
    Tj.set_number_of_servers(m, 1)
    Tj.set_enabling_conditions(m, cls, P[1], 1)
    Tj.set_enabling_conditions(m, cls, P[2], 1)
    Tj.set_firing_outcome(m, cls, P[3], 1)
    m = Tb.add_mode('b')
    Tb.set_distribution(m, Exp(1.9))
    Tb.set_number_of_servers(m, 1)
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


# --- A closed cycle, where the exact CTMC gives the reference
print('\n--- 3-place cyclic net, N = 4 ---')
print(SolverNC(cyclic()).get_avg_table())
print('The same net through the explicit generator, for comparison:')
print(SolverCTMC(cyclic()).get_avg_table())

# The certificate the derivation produced. Deficiency zero plus weak
# reversibility is what Feinberg's theorem needs for a positive complex-balanced
# point to exist at ANY choice of rates.
pf = spn_pf(cyclic(), {'verbose': True})
print('product form: %s, deficiency %d, %d linkage classes, rank %d'
      % (pf['kind'], pf['deficiency'], pf['linkage'], pf['srank']))

# --- A fork-join net, which has no queueing-network form at all
print('\n--- fork-join net, P0 -> P1+P2 -> P3 -> P0, 3 tokens at P0 ---')
print(SolverNC(forkjoin()).get_avg_table())

pfj = spn_pf(forkjoin())
met = spn_metrics(pfj['mdds'], pfj['g'], pfj['info'])
# Every token that forks must later join and return, so the three modes share
# one throughput -- a flow-conservation law nothing in the derivation was told.
print('mode throughputs: %s (they must all agree)' % np.round(met['modeTput'], 8))
# The place invariant of this net is 2*m0 + m1 + m2 + 2*m3 = 6, not the token
# count, which is what "not a conserved population" means concretely.
print('place invariant 2*m0 + m1 + m2 + 2*m3 = %.6f'
      % float(np.array([2, 1, 1, 2]) @ np.asarray(met['tokens'])))
