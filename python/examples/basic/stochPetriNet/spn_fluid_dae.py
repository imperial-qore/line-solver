"""
Fluid (Mean-Field) Analysis of a Stochastic Petri Net, with SolverFLD

A GSPN is a density-dependent Markov population process: the marking is the
population, a transition mode is a reaction, and the firing rate
lambda*min(enabling degree, servers) is the same min() non-linearity the
min-normal closure of SolverFLD exists to smooth. The 'dae' method is the one
that can carry it, because a Petri net needs three things stated as EQUATIONS
rather than integrated:

  * the P-invariants, which hold to solver tolerance instead of integrator
    tolerance -- and supply the rank the drift Jacobian is missing;
  * the firing FLOW of an immediate transition, an algebraic unknown pinned by
    the constraint that its input place holds no mass;
  * a bounded place, a linear inequality on the marking.

SolverFLD(model) resolves to 'dae' on any model holding a Transition node, so
no method has to be named. Unlike every other solver of a Petri net in LINE it
also returns a SECOND MOMENT: the marking covariance of the linear noise
approximation, through getMoments().

Twin of `matlab/examples/basic/stochPetriNet/spn_fluid_dae.m`.
"""

from line_solver import *


def spn_fluid_exact():
    """A closed net whose fluid answer is EXACT.

    Every mode is infinite-server with one input arc, so min(m/w, Inf) = m and
    the drift is LINEAR: the fluid mean is then the exact mean, and the
    covariance the exact covariance (a binomial marking).
    """
    exact = Network('spn_fluid_exact')
    P1 = Place(exact, 'P1')
    P2 = Place(exact, 'P2')
    T1 = Transition(exact, 'T1')
    T2 = Transition(exact, 'T2')
    jc = ClosedClass(exact, 'Class1', 4, P1, 0)

    m1 = T1.add_mode('Mode1')
    T1.set_number_of_servers(m1, float('inf'))
    T1.set_distribution(m1, Exp(2))
    T1.set_enabling_conditions(m1, jc, P1, 1)
    T1.set_firing_outcome(m1, jc, P2, 1)

    m2 = T2.add_mode('Mode2')
    T2.set_number_of_servers(m2, float('inf'))
    T2.set_distribution(m2, Exp(3))
    T2.set_enabling_conditions(m2, jc, P2, 1)
    T2.set_firing_outcome(m2, jc, P1, 1)

    R = exact.init_routing_matrix()
    R.set(jc, jc, P1, T1, 1.0)
    R.set(jc, jc, T1, P2, 1.0)
    R.set(jc, jc, P2, T2, 1.0)
    R.set(jc, jc, T2, P1, 1.0)
    exact.link(R)

    P1.set_state(jc.get_population())
    P2.set_state(0)
    return exact, P1, P2


def spn_fluid_immediate():
    """An immediate transition, as an algebraic flow.

    P1 -T1-> P2 -(immediate)-> P3 -T2-> P1. The vanishing place P2 holds exactly
    zero mass and the net answers as the reduced two-place net does, which is
    what the algebraic flow buys: an approximation of the immediate transition
    by a large finite rate would only approach it.
    """
    imm = Network('spn_fluid_immediate')
    Q1 = Place(imm, 'P1')
    Q2 = Place(imm, 'P2')
    Q3 = Place(imm, 'P3')
    U1 = Transition(imm, 'T1')
    Ui = Transition(imm, 'Ti')
    U3 = Transition(imm, 'T3')
    jq = ClosedClass(imm, 'Class1', 4, Q1, 0)

    a1 = U1.add_mode('M1')
    U1.set_distribution(a1, Exp(2))
    U1.set_enabling_conditions(a1, jq, Q1, 1)
    U1.set_firing_outcome(a1, jq, Q2, 1)

    ai = Ui.add_mode('Mi')
    Ui.set_distribution(ai, Immediate())
    Ui.set_timing_strategy(ai, TimingStrategy.IMMEDIATE)
    Ui.set_enabling_conditions(ai, jq, Q2, 1)
    Ui.set_firing_outcome(ai, jq, Q3, 1)

    a3 = U3.add_mode('M3')
    U3.set_distribution(a3, Exp(3))
    U3.set_enabling_conditions(a3, jq, Q3, 1)
    U3.set_firing_outcome(a3, jq, Q1, 1)

    Ri = imm.init_routing_matrix()
    Ri.set(jq, jq, Q1, U1, 1.0)
    Ri.set(jq, jq, U1, Q2, 1.0)
    Ri.set(jq, jq, Q2, Ui, 1.0)
    Ri.set(jq, jq, Ui, Q3, 1.0)
    Ri.set(jq, jq, Q3, U3, 1.0)
    Ri.set(jq, jq, U3, Q1, 1.0)
    imm.link(Ri)

    Q1.set_state(jq.get_population())
    Q2.set_state(0)
    Q3.set_state(0)
    return imm


def spn_fluid_reduced():
    """The same net with the immediate transition eliminated by hand."""
    red = Network('spn_fluid_reduced')
    W1 = Place(red, 'P1')
    W3 = Place(red, 'P3')
    V1 = Transition(red, 'T1')
    V3 = Transition(red, 'T3')
    jr = ClosedClass(red, 'Class1', 4, W1, 0)

    b1 = V1.add_mode('M1')
    V1.set_distribution(b1, Exp(2))
    V1.set_enabling_conditions(b1, jr, W1, 1)
    V1.set_firing_outcome(b1, jr, W3, 1)

    b3 = V3.add_mode('M3')
    V3.set_distribution(b3, Exp(3))
    V3.set_enabling_conditions(b3, jr, W3, 1)
    V3.set_firing_outcome(b3, jr, W1, 1)

    Rr = red.init_routing_matrix()
    Rr.set(jr, jr, W1, V1, 1.0)
    Rr.set(jr, jr, V1, W3, 1.0)
    Rr.set(jr, jr, W3, V3, 1.0)
    Rr.set(jr, jr, V3, W1, 1.0)
    red.link(Rr)

    W1.set_state(jr.get_population())
    W3.set_state(0)
    return red


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    exact, P1, P2 = spn_fluid_exact()
    fld = FLD(exact)
    print(fld.getAvgTable())
    print(CTMC(exact, cutoff=6).getAvgTable())

    # The exact marking is Binomial(4, 3/5), so the variance is 4*0.6*0.4 = 0.96.
    mom = fld.getMoments()
    petri = mom['petri']
    print('marking variance: %.6f %.6f  (exact 0.96)' % (
        petri['markingVar'][P1.get_index0(), 0],
        petri['markingVar'][P2.get_index0(), 0]))
    print('invariant "%s" = %g, error %.2e' % (
        petri['invariantLabel'][0],
        petri['invariantValue'][0],
        petri['invariantError'][0]))

    print(FLD(spn_fluid_immediate()).getAvgTable())
    print(FLD(spn_fluid_reduced()).getAvgTable())
