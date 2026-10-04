"""
Stochastic network calculus on a feed-forward network: envelope propagation,
blind multiplexing, and the model classes the family refuses.

The bound of snc_delay_quantile.py is a single station. This example is what the
network calculus adds on top of it: a departure envelope carries a flow to the
next hop, cross traffic is subtracted from a shared server, and both operations
cost burstiness, which is visible as the bound loosening downstream and under
sharing.

The reference throughout is the exact Jackson-network solution, since every
station here is an M/M/1 with Poisson input by Burke's theorem. Twin of
matlab/examples/advanced/networkCalculus/snc_tandem_multiclass.m.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

from line_solver import *
from line_solver.api.snc import (snc_conv, snc_env_poisson, snc_perc_delay,
                                 snc_srv_exp)

print('=== Stochastic network calculus: tandem and shared server ===\n')

lam = 0.6
rates = [1.5, 1.2, 1.0]          # service rates, the last is the bottleneck

model = Network('SncTandem')
source = Source(model, 'Source')
q1 = Queue(model, 'Q1', SchedStrategy.FCFS)
q2 = Queue(model, 'Q2', SchedStrategy.FCFS)
q3 = Queue(model, 'Q3', SchedStrategy.FCFS)
sink = Sink(model, 'Sink')
jobclass = OpenClass(model, 'Class1')
source.setArrival(jobclass, Exp(lam))
q1.setService(jobclass, Exp(rates[0]))
q2.setService(jobclass, Exp(rates[1]))
q3.setService(jobclass, Exp(rates[2]))
model.link(Network.serialRouting(source, q1, q2, q3, sink))

tandem = SolverBA(model, 'snc.upper').getAvgTable()
print(tandem)

# Each hop replaces the arrival envelope by the DEPARTURE envelope of the station
# upstream, which carries the burst the server has added. The exact answer does not
# degrade this way -- by Burke's theorem the departure process of an M/M/1 is again
# Poisson -- so the ratio to the exact response time grows hop by hop. The bound
# stays valid; it is the price of assuming nothing about the departure process
# beyond its envelope.
print('\n%-8s %10s %10s %8s' % ('station', 'R bound', 'R exact', 'ratio'))
for i in range(3):
    exact = 1.0 / (rates[i] - lam)
    bound = float(tandem['RespT'].iloc[i + 1])
    print('%-8s %10.4f %10.4f %8.2f'
          % (tandem['Station'].iloc[i + 1], bound, exact, bound / exact))

# Summing the per-station bounds pays the burst term at every hop. Concatenating
# the three service envelopes with snc_conv first and bounding the composed element
# once pays it only once, which is the classical result of the network calculus and
# is worth several tens of percent here.
arv = lambda theta: snc_env_poisson(lam, theta)


def tandem_service(theta):
    """Min-plus concatenation of the three service envelopes into one element."""
    sigma, rho = snc_srv_exp(rates[0], theta)
    for rate in rates[1:]:
        s2, r2 = snc_srv_exp(rate, theta)
        sigma, rho = snc_conv(sigma, rho, s2, r2, theta)
    return sigma, rho


d_end, theta_end = snc_perc_delay(arv, tandem_service, 1e-3)
d_hop = sum(snc_perc_delay(arv, (lambda r: lambda t: snc_srv_exp(r, t))(rate), 1e-3 / 3)[0]
            for rate in rates)
print('\nend-to-end delay quantile at eps=1e-3')
print('  concatenated (snc_conv) : %8.4f  at theta = %.4f' % (d_end, theta_end))
print('  summed per hop          : %8.4f' % d_hop)
print('  pay bursts once saves   : %7.1f%%' % (100.0 * (1.0 - d_end / d_hop)))

# A class sharing a station sees the server minus whatever the other classes take
# from it: snc_leftover subtracts the cross-flow arrival envelope from the service
# envelope. The result holds for ANY work-conserving discipline at that station,
# which is why it is well above the FCFS answer -- it also covers the policy that
# serves the other class first whenever it can.
shared = Network('SncShared')
src2 = Source(shared, 'Source')
qs = Queue(shared, 'Shared', SchedStrategy.FCFS)
snk2 = Sink(shared, 'Sink')
class_a = OpenClass(shared, 'ClassA')
class_b = OpenClass(shared, 'ClassB')
src2.setArrival(class_a, Exp(0.3))
src2.setArrival(class_b, Exp(0.3))
qs.setService(class_a, Exp(1.0))
qs.setService(class_b, Exp(1.0))
P = shared.initRoutingMatrix()
P.set(class_a, Network.serialRouting(src2, qs, snk2))
P.set(class_b, Network.serialRouting(src2, qs, snk2))
shared.link(P)

shared_solver = SolverBA(shared, 'snc.upper')
print()
print(shared_solver.getAvgTable())
print()
print(shared_solver.getPercTable(1e-3))
print('exact per-class response time (aggregate M/M/1, lam=0.6, mu=1): %.4f'
      % (1.0 / (1.0 - 0.6)))


# The elementary envelope algebra has real limits, and the analyzer states them
# rather than returning a number that looks plausible. Each refusal below is a
# modelling assumption of the calculus, not an implementation gap.
def show_refusal(what, build):
    try:
        build()
        print('  %-45s NO REFUSAL (unexpected)' % what)
    except Exception as err:                                  # noqa: BLE001
        print('  %-45s %s' % (what + ':', err))


def closed_model():
    m = Network('Closed')
    delay = Delay(m, 'Think')
    queue = Queue(m, 'Q', SchedStrategy.PS)
    jobs = ClosedClass(m, 'C', 3, delay)
    delay.setService(jobs, Exp(1.0))
    queue.setService(jobs, Exp(2.0))
    m.link(Network.serialRouting(delay, queue))
    SolverBA(m, 'snc.upper').getAvgTable()


def split_model():
    m = Network('Split')
    src = Source(m, 'Source')
    qa = Queue(m, 'QA', SchedStrategy.FCFS)
    qb = Queue(m, 'QB', SchedStrategy.FCFS)
    snk = Sink(m, 'Sink')
    c = OpenClass(m, 'C')
    src.setArrival(c, Exp(0.4))
    qa.setService(c, Exp(1.0))
    qb.setService(c, Exp(1.0))
    P = m.initRoutingMatrix()
    P.set(c, c, src, qa, 1.0)
    P.set(c, c, qa, qb, 0.5)
    P.set(c, c, qa, snk, 0.5)
    P.set(c, c, qb, snk, 1.0)
    m.link(P)
    SolverBA(m, 'snc.upper').getAvgTable()


def unequal_model():
    m = Network('Unequal')
    src = Source(m, 'Source')
    queue = Queue(m, 'Q', SchedStrategy.FCFS)
    snk = Sink(m, 'Sink')
    ca = OpenClass(m, 'A')
    cb = OpenClass(m, 'B')
    src.setArrival(ca, Exp(0.2))
    src.setArrival(cb, Exp(0.2))
    queue.setService(ca, Exp(1.0))
    queue.setService(cb, Exp(2.0))
    P = m.initRoutingMatrix()
    P.set(ca, Network.serialRouting(src, queue, snk))
    P.set(cb, Network.serialRouting(src, queue, snk))
    m.link(P)
    SolverBA(m, 'snc.upper').getAvgTable()


print('\nrefusals:')
show_refusal('closed network', closed_model)
show_refusal('probabilistic split downstream of the Source', split_model)
show_refusal('unequal service rates at a shared station', unequal_model)
