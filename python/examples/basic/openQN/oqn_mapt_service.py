"""
Open queue with a MAP_t SERVICE process, simulated by LDES.

A MAP_t service is the only process LDES simulates that is both non-renewal and
non-stationary. Segment k covers [breakpoints[k], breakpoints[k+1]) and carries
the pair (D0[k], D1[k]); D1 fires a service COMPLETION and D0 only moves the
modulating phase. The engine walks that phase process forward from the instant
the job ENTERS SERVICE -- not from the instant it arrived, which may fall under
a different segment -- and resumes from the phase the previous completion left
behind, so successive service times are correlated exactly as for an ordinary
MAP.

The schedule below is an MMPP whose two environment states hold the same sojourn
rates throughout (Q_env) while their completion intensities change from segment
to segment, so the server alternates between a fast and a slow mode AND the pair
of modes itself changes with the clock.

Because the walk is a wall-clock one, it is exact only where service runs
continuously at unit rate once started. LDES therefore refuses a MAP_t (or PH_t)
service under processor sharing, under a preemptive discipline, with load
dependence and with heterogeneous server pools, rather than returning a number
from a sample path that does not model any of those. Use INF or a non-preemptive
FCFS/LCFS-family discipline, as here.
"""

import numpy as np

from line_solver import (LDES, MAP, MAPt, Exp, Network, OpenClass, Queue,
                         SchedStrategy, Sink, Source)


def oqn_mapt_service():
    model = Network('model')

    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')

    jobclass = OpenClass(model, 'OpenClass', 0)

    # Three segments held for 1, 1, 2 time units, repeating with period 4.
    breakpoints = [0.0, 1.0, 2.0, 4.0]
    Q_env = np.array([[-1.0, 1.0], [2.0, -2.0]])
    gamma = -np.diag(Q_env)
    mu_values = np.array([[2.0, 8.0, 4.0],
                          [4.0, 2.0, 8.0]])
    D0 = []
    D1 = []
    for segment in range(len(breakpoints) - 1):
        a = Q_env.copy()
        a[0, 0] = -gamma[0] - mu_values[0, segment]
        a[1, 1] = -gamma[1] - mu_values[1, segment]
        D0.append(a)
        D1.append(np.diag(mu_values[:, segment]))

    source.setArrival(jobclass, Exp(10.0))
    queue.setService(jobclass, MAPt(breakpoints, D0, D1, True))
    queue.setCapacity(5)
    queue.setNumberOfServers(1)

    model.link(Network.serialRouting(source, queue, sink))
    return model


def _flat_and_homogeneous():
    """The degenerate pair: a constant MAP_t schedule and the MAP it must equal."""
    Q_env = np.array([[-1.0, 1.0], [2.0, -2.0]])
    D0bar = Q_env - np.diag([5.0, 10.0])
    D1bar = np.diag([5.0, 10.0])

    flat = Network('flat')
    src = Source(flat, 'Source')
    q = Queue(flat, 'Queue', SchedStrategy.FCFS)
    snk = Sink(flat, 'Sink')
    cls = OpenClass(flat, 'OpenClass', 0)
    src.setArrival(cls, Exp(1.0))
    q.setService(cls, MAPt([0.0, 0.25, 0.5], [D0bar, D0bar], [D1bar, D1bar], True))
    flat.link(Network.serialRouting(src, q, snk))

    homog = Network('homog')
    src2 = Source(homog, 'Source')
    q2 = Queue(homog, 'Queue', SchedStrategy.FCFS)
    snk2 = Sink(homog, 'Sink')
    cls2 = OpenClass(homog, 'OpenClass', 0)
    src2.setArrival(cls2, Exp(1.0))
    q2.setService(cls2, MAP(D0bar, D1bar))
    homog.link(Network.serialRouting(src2, q2, snk2))
    return flat, homog


if __name__ == "__main__":
    # Built by the function above rather than inline: a model that exists only
    # under __main__ exposes nothing on import, so the JAVA and C++ parity rows
    # cannot export it and SKIP every solver -- a row that reads as coverage
    # while asserting nothing.
    model = oqn_mapt_service()

    # The offered load exceeds what the server can clear, so the finite buffer
    # drops: Tput falls short of the arrival rate and QLen sits just under the
    # capacity. Both are signatures that the service process is being simulated
    # -- a service that sampled as zero would report QLen = Util = 0, Tput = 10.
    print(LDES(model, seed=1234, samples=1000000).getAvgTable())

    # A MAP_t whose segments are identical carries no time dependence, so it
    # must reproduce the ordinary MAP with those matrices. This is the
    # degeneracy that pins the schedule machinery: it exercises every boundary
    # crossing and still has to land on the time-homogeneous answer.
    flat, homog = _flat_and_homogeneous()
    t_flat = LDES(flat, seed=1234, samples=400000).getAvgTable()
    t_homog = LDES(homog, seed=1234, samples=400000).getAvgTable()
    print('\nconstant schedule vs ordinary MAP (must agree):')
    print(t_flat)
    print(t_homog)
