"""
Colored Generalized Stochastic Petri Net (CGSPN)

This example demonstrates:
- Token colors: a color is a job class, so a place holds one marking per color
- One transition mode per color, with color-dependent rates
- An immediate transition (weights arbitrate between colors) and a timed one
- Inhibiting arcs making the server a mutual-exclusion resource

Two colors, Gold and Silver, circulate between a Buffer place and a Server
place. Admission is immediate, so the server is never idle: the firing weights
split the completions 2:1 in favour of Gold, hence X_Gold = 2*X_Silver,
U_c = X_c/mu_c and U_Gold + U_Silver = 1.
"""

from line_solver import *


def build_model() -> Network:
    """Build and return the colored GSPN model."""
    model = Network('ColoredGSPN')

    # Declare every node before parameterizing any of it: the enabling,
    # inhibiting and firing matrices are indexed by (node, class) over the
    # nodes that exist when they are first set.
    buf = Place(model, 'Buffer')
    srv = Place(model, 'Server')
    admit = Transition(model, 'admit')
    serve = Transition(model, 'serve')

    # Two token colors, both starting in the buffer
    gold = ClosedClass(model, 'Gold', 2, buf, 0)
    silver = ClosedClass(model, 'Silver', 2, buf, 0)

    # Immediate transition: one mode per color, admitting a token only when the
    # server holds no token of either color. The weights decide the color mix.
    mg = admit.add_mode('gold')
    admit.set_timing_strategy(mg, TimingStrategy.IMMEDIATE)
    admit.set_enabling_conditions(mg, gold, buf, 1)
    admit.set_inhibiting_conditions(mg, gold, srv, 1)
    admit.set_inhibiting_conditions(mg, silver, srv, 1)
    admit.set_firing_outcome(mg, gold, srv, 1)
    admit.set_firing_weights(mg, 2.0)

    ms = admit.add_mode('silver')
    admit.set_timing_strategy(ms, TimingStrategy.IMMEDIATE)
    admit.set_enabling_conditions(ms, silver, buf, 1)
    admit.set_inhibiting_conditions(ms, gold, srv, 1)
    admit.set_inhibiting_conditions(ms, silver, srv, 1)
    admit.set_firing_outcome(ms, silver, srv, 1)
    admit.set_firing_weights(ms, 1.0)

    # Timed transition: one mode per color, with color-dependent service rates
    sg = serve.add_mode('gold')
    serve.set_distribution(sg, Exp(3.0))
    serve.set_enabling_conditions(sg, gold, srv, 1)
    serve.set_firing_outcome(sg, gold, buf, 1)

    ss = serve.add_mode('silver')
    serve.set_distribution(ss, Exp(1.5))
    serve.set_enabling_conditions(ss, silver, srv, 1)
    serve.set_firing_outcome(ss, silver, buf, 1)

    # Topology: both colors follow Buffer -> admit -> Server -> serve -> Buffer
    R = model.init_routing_matrix()
    for c in (gold, silver):
        R.set(c, c, buf, admit, 1.0)
        R.set(c, c, admit, srv, 1.0)
        R.set(c, c, srv, serve, 1.0)
        R.set(c, c, serve, buf, 1.0)
    model.link(R)

    # Initial marking: two tokens of each color in the buffer, server empty
    buf.set_state([2, 2])
    srv.set_state([0, 0])

    return model


spn_colored_gspn = build_model


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = build_model()

    # Exact solution of the underlying CTMC, and a simulation cross-check
    print(SolverCTMC(model, cutoff=4).avg_table())
    print(SolverSSA(model, seed=23000, samples=200000).avg_table())
