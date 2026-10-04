"""
Open queueing network whose arrival stream is a BATCH, MARKED,
time-inhomogeneous MAP.

BMMAPt crosses the three axes LINE models for an arrival stream. A block is
indexed by segment, mark and batch size: D1kb[c][b][j] holds the rates that, in
segment j, release a batch of b+1 jobs ALL carrying mark c+1, and
D0[j] + sum_c sum_b D1kb[c][b][j] is a generator in every segment. At a Source
the mark selects the CLASS of the arriving jobs and the batch size says HOW MANY
of them arrive at that instant.

Two derived levels come for free and are what keep the family legible: summing
over b gives the MMAPt of the marks alone, and summing over c as well gives the
MAPt of the epochs alone. So a consumer that ignores batches sees exactly the
marked schedule, and one that ignores marks too sees exactly the unmarked one.

The schedule below moves BOTH the class mix and the batch size, which is the
point: it is the combination no existing family expresses. Every segment fires
epochs at rate 4, so the epoch stream is statistically identical throughout;
what changes is which class the batch carries and how big it is. In the morning
segment the Premium class arrives in PAIRS and Economy singly; in the evening
segment that reverses. The JOB rate therefore differs from the EPOCH rate, and
differs per class within a segment even though the epoch rate does not -- a
model that read either label off the time-averaged matrices would miss all of it.

Reductions worth knowing: with every batch size 1 a BMMAPt IS the MMAPt with the
same blocks, sample path for sample path on one seed; with K = 1 it is the
unmarked batch schedule; with both it is the MAPt; and with identical segments it
is the stationary BMAP.

References: Q.-M. He, "The versatility of MMAP[K] and the MMAP[K]/G[K]/1 queue",
Queueing Systems 38(4), 2001, for the marked structure; D. M. Lucantoni, "New
results on the single server queue with a batch Markovian arrival process",
Stochastic Models 7(1), 1991, for the batch structure; Y. M. Ko and J. Pender,
"Diffusion limits for the (MAP_t/Ph_t/inf)^N queueing network", Oper. Res. Lett.
45(3), 2017, for the time-inhomogeneous one.
"""

from line_solver import (BMMAPt, Exp, Network, OpenClass, Queue, SchedStrategy,
                         Sink, Source, SolverLDES)


def oqn_bmmapt():
    model = Network('model')

    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')

    premium_class = OpenClass(model, 'Premium', 0)
    economy_class = OpenClass(model, 'Economy', 0)

    # Two twelve-hour segments repeating on a daily cycle. One phase, so the
    # epoch stream is Poisson at rate 4 throughout and only the labels move.
    breakpoints = [0.0, 12.0, 24.0]
    D0 = [[[-4.0]], [[-4.0]]]

    # D1kb[mark][batch][segment]. The batch axis is DENSE: a mark that never
    # releases a batch of that size still declares a zero block for it, exactly
    # as BMAP's [D0, D1, ..., Dk] does.
    #
    #   morning (segment 1): Premium in PAIRS at 3.6, Economy singly at 0.4
    #   evening (segment 2): Premium singly at 0.4,   Economy in PAIRS at 3.6
    D1kb = [[[[[0.0]], [[0.4]]],    # mark 1 (Premium), batch 1
             [[[3.6]], [[0.0]]]],   # mark 1 (Premium), batch 2
            [[[[0.4]], [[0.0]]],    # mark 2 (Economy), batch 1
             [[[0.0]], [[3.6]]]]]   # mark 2 (Economy), batch 2
    arrival = BMMAPt(breakpoints, D0, D1kb, True)

    source.setMarkedArrival(arrival, [premium_class, economy_class])
    queue.setService(premium_class, Exp(16.0))
    queue.setService(economy_class, Exp(16.0))

    model.link(Network.serialRouting(source, queue, sink))
    return model, arrival


if __name__ == "__main__":
    # Built by the function above rather than inline: a model that exists only
    # under __main__ exposes nothing on import, so the JAVA and C++ parity rows
    # cannot export it and SKIP every solver.
    model, arrival = oqn_bmmapt()

    print('Epoch (batch) rate over the cycle: %.4f' % arrival.getTimeAverageRate())
    print('JOB rate over the cycle:           %.4f  (epochs times the mean batch)'
          % arrival.getTimeAverageJobRate())
    print('Mean batch size:                   %.4f' % arrival.getMeanBatchSize())
    lam = arrival.getTimeAverageMarkRates()
    print('Per-mark EPOCH rates:              %.4f %.4f  (they sum to the epoch rate)'
          % (lam[0], lam[1]))
    jobs = arrival.getTimeAverageMarkJobRates()
    print('Per-mark JOB rates:                %.4f %.4f  (they sum to the job rate)'
          % (jobs[0], jobs[1]))
    print('The two classes are symmetric over the cycle, so the averages match; the '
          'run below\nshows the SEGMENTS, where they do not.\n')

    # LDES simulates the batch marked schedule directly: the mark and the batch
    # size are both decided by WHICH block's transition fired, inside the same
    # competing-transitions draw that ends the interval, so neither label costs
    # an extra random draw and the reduction to MMAPt at batch size 1 is exact.
    solver = SolverLDES(model, samples=4e5, seed=23000)
    print(solver.getAvgTable())
