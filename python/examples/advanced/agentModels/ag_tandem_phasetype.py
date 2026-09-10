"""
Open Tandem Queue with Phase-Type Service, solved by the agent-based solver

This example demonstrates that the agent-based methods represent a phase-type service
law exactly rather than collapsing it to its mean rate. Each component of the
agent decomposition is a QBD whose level is the queue length and whose phase is
the pair (arrival phase, service phase), so the first station -- an isolated
M/PH/1, since it sees the Poisson source directly -- comes out at the
Pollaczek-Khinchine mean whatever the reversed-rate iteration does. Both
stations here have the same mean service time and differ only in their
variability, which is exactly what an M/M/1 reading cannot see.

References: Casale and Harrison, "AutoCAT: Automated Product-Form Solution of
            Stochastic Models", Stochastic Models 27, 2013
            Neuts, "Matrix-Geometric Solutions in Stochastic Models", 1981

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

from line_solver import *

# Parameters
arrival_rate = 0.5   # Arrival rate
mean_service = 1.0   # Mean service time at both queues
scv1 = 0.5           # Erlang-2 service at queue 1 (less variable than exponential)
scv2 = 4.0           # HyperExp service at queue 2 (more variable)

# Create model: Source -> Queue1 -> Queue2 -> Sink
model = Network('Tandem-MPH1')

source = Source(model, 'Source')
queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
sink = Sink(model, 'Sink')

oclass = OpenClass(model, 'Class1')
source.setArrival(oclass, Exp(arrival_rate))
queue1.setService(oclass, Erlang.fitMeanAndSCV(mean_service, scv1))
queue2.setService(oclass, HyperExp.fitMeanAndSCV(mean_service, scv2))

model.link(Network.serial_routing([source, queue1, queue2, sink]))

# The exact M/G/1 mean at Queue1, which sees the Poisson source directly
rho = arrival_rate * mean_service
pk1 = rho + rho ** 2 * (1 + scv1) / (2 * (1 - rho))

print('=== Open Tandem with Phase-Type Service ===\n')
print('Queue1 is an isolated M/Er2/1, so its exact mean queue length is the')
print(f'Pollaczek-Khinchine value {pk1:.6f}. The M/M/1 reading would be '
      f'{rho / (1 - rho):.6f}.\n')

print('AG (method=inap):')
print(AG(model, 'inap').get_avg_table())

# 'inapinf' additionally drops the maxStates truncation, solving each open
# component on its infinite state space through Neuts' rate matrix R. That
# matters most at Queue2, whose service law has the heavier tail.
print('\nAG (method=inapinf):')
print(AG(model, 'inapinf').get_avg_table())
