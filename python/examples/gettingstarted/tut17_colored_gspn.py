"""Tutorial 17: colored generalized stochastic Petri net."""

from line_solver import *


model = Network("ColoredGSPN")
buffer = Place(model, "Buffer")
server = Place(model, "Server")
admit = Transition(model, "admit")
serve = Transition(model, "serve")
gold = ClosedClass(model, "Gold", 2, buffer, 0)
silver = ClosedClass(model, "Silver", 2, buffer, 0)

for token_class, name, weight in ((gold, "gold", 2.0), (silver, "silver", 1.0)):
    mode = admit.add_mode(name)
    admit.set_timing_strategy(mode, TimingStrategy.IMMEDIATE)
    admit.set_enabling_conditions(mode, token_class, buffer, 1)
    admit.set_inhibiting_conditions(mode, gold, server, 1)
    admit.set_inhibiting_conditions(mode, silver, server, 1)
    admit.set_firing_outcome(mode, token_class, server, 1)
    admit.set_firing_weights(mode, weight)

for token_class, name, rate in ((gold, "gold", 3.0), (silver, "silver", 1.5)):
    mode = serve.add_mode(name)
    serve.set_distribution(mode, Exp(rate))
    serve.set_enabling_conditions(mode, token_class, server, 1)
    serve.set_firing_outcome(mode, token_class, buffer, 1)

routing = model.init_routing_matrix()
for token_class in (gold, silver):
    routing.set(token_class, token_class, buffer, admit, 1.0)
    routing.set(token_class, token_class, admit, server, 1.0)
    routing.set(token_class, token_class, server, serve, 1.0)
    routing.set(token_class, token_class, serve, buffer, 1.0)
model.link(routing)
buffer.set_state([2, 2])
server.set_state([0, 0])

print(SolverCTMC(model, cutoff=4).avg_table())
print(SolverSSA(model, seed=23000, samples=200000).avg_table())
