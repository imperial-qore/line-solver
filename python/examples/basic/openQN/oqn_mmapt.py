"""
Open queueing network whose arrival stream is a MARKED, time-inhomogeneous MAP.

MMAPt crosses the two axes of MAPt and MarkedMAP: arrivals are labelled with
one of K marks, AND the matrices that generate them are functions of the wall
clock. Segment j carries D0[j] together with the K blocks D1k[c][j], and
D0[j] + sum_c D1k[c][j] is a generator in every segment. The aggregate
sum_c D1k[c][j] is the D1 of the underlying MAPt, so hiding the marks recovers
exactly that process.

At a Source the mark SELECTS THE CLASS of the arriving job: one modulating
chain drives every marked class, and only the first one carries the stream.
That is what makes this different from declaring K independent arrival
processes -- the classes are correlated through the shared phase, and the class
MIX shifts with the schedule even when the total rate does not.

The schedule below is built so that only the mix moves. Both segments carry an
aggregate rate of 4, so the total arrival stream is statistically identical
throughout; what changes is the split, 9:1 towards the first class in the
morning segment and 1:9 towards the second in the evening one. A model that
read the marks off the time-averaged matrices would report a flat 2:2 split.

Reductions worth knowing: with K = 1 an MMAPt IS the MAPt with the same
matrices, and with identical segments it IS the stationary MMAP. MPHt is the
phase-type twin, stored lowered to this same form.

Reference for the marked structure: Q.-M. He, "The versatility of MMAP[K] and
the MMAP[K]/G[K]/1 queue", Queueing Systems 38(4), 2001; for the
time-inhomogeneous one: Y. M. Ko, J. Pender, "Diffusion limits for the
(MAP_t/Ph_t/inf)^N queueing network", Oper. Res. Lett. 45(3), 2017.
"""

import numpy as np

from line_solver import (Exp, MMAPt, Network, OpenClass, Queue, SchedStrategy,
                         Sink, SolverLDES, Source)


def oqn_mmapt():
    model = Network('model')

    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')

    morningClass = OpenClass(model, 'Morning', 0)
    eveningClass = OpenClass(model, 'Evening', 0)

    # Two twelve-hour segments repeating on a daily cycle. One phase, so the
    # aggregate stream is Poisson at rate 4 throughout and ONLY the mix moves.
    breakpoints = [0.0, 12.0, 24.0]
    D0 = [np.array([[-4.0]]), np.array([[-4.0]])]
    D1k = [[np.array([[3.6]]), np.array([[0.4]])],   # mark 1 -> Morning
           [np.array([[0.4]]), np.array([[3.6]])]]   # mark 2 -> Evening
    arrival = MMAPt(breakpoints, D0, D1k, True)

    source.setMarkedArrival(arrival, [morningClass, eveningClass])
    queue.setService(morningClass, Exp(8.0))
    queue.setService(eveningClass, Exp(8.0))

    model.link(Network.serialRouting(source, queue, sink))

    # The arrival process rides along: the marked schedule is what the rate
    # block below probes, and find_model searches a returned tuple element-wise,
    # so handing both back keeps the model exportable.
    return model, arrival


if __name__ == "__main__":
    model, arrival = oqn_mmapt()

    print('Aggregate arrival rate over the cycle: %.4f'
          % arrival.getTimeAverageRate())
    lam = arrival.getTimeAverageMarkRates()
    print('Per-mark rates over the cycle:         %.4f %.4f  '
          '(they sum to the aggregate)' % (lam[0], lam[1]))
    print('The cycle average is a 2:2 split, so a run that shows anything else '
          'is reading\nthe marks from the segment in force rather than from the '
          'average.\n')

    # LDES simulates the marked schedule directly: the mark is decided by WHICH
    # block's transition fired, inside the same competing-transitions draw that
    # ends the interval, so no extra randomness is introduced by the labelling.
    print(SolverLDES(model, samples=4e5, seed=23000).getAvgTable())
