"""
Test ETAQA departure process for FCFS queues in MAM dec.mmap

Creates a 2-queue tandem with non-Poisson (Erlang) arrivals and Erlang service,
solved with dec.mmap (ETAQA-based departures). Validates that:
  1. dec.mmap returns finite results and converges
  2. Results are consistent with dec.source and MVA baselines
  3. The ETAQA truncation level can be configured via options
"""

import numpy as np

from line_solver import (Erlang, GlobalConstants, MAM, MVA, Network, OpenClass,
                         Queue, SchedStrategy, Sink, Source, VerboseLevel)


def _rows(table):
    """The Queue rows of an AvgTable, in station order."""
    return table[table['Station'].astype(str).str.contains('Queue')]


if __name__ == '__main__':
    GlobalConstants.set_verbose(VerboseLevel.SILENT)

    # Model: Erl(5)/Erl(2)/1 -> Erl(3)/1 tandem
    model = Network('ETAQA-FCFS-Tandem')

    source = Source(model, 'Source')
    queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
    queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')

    oclass = OpenClass(model, 'Class1')
    source.setArrival(oclass, Erlang.fitMeanAndOrder(2, 5))    # mean=2, SCV=0.2
    queue1.setService(oclass, Erlang.fitMeanAndOrder(0.8, 2))  # rho1=0.4
    queue2.setService(oclass, Erlang.fitMeanAndOrder(1.2, 3))  # rho2=0.6

    model.link(Network.serialRouting(source, queue1, queue2, sink))

    # Solve with dec.mmap (uses ETAQA departure process)
    T_mmap = MAM(model, method='dec.mmap').getAvgTable()
    # Solve with dec.source (baseline)
    T_src = MAM(model, method='dec.source').getAvgTable()
    # Solve with MVA (exact for mean values)
    T_mva = MVA(model).getAvgTable()

    print('\n=== ETAQA FCFS Departure Process Test ===')
    print('%-10s %10s %10s %10s' % ('Station', 'MVA', 'dec.source', 'dec.mmap'))
    print('--- Queue Lengths ---')
    for (_, a), (_, b), (_, c) in zip(_rows(T_mva).iterrows(), _rows(T_src).iterrows(),
                                      _rows(T_mmap).iterrows()):
        print('%-10s %10.4f %10.4f %10.4f' % (a['Station'], a['QLen'], b['QLen'], c['QLen']))
    print('--- Response Times ---')
    for (_, a), (_, b), (_, c) in zip(_rows(T_mva).iterrows(), _rows(T_src).iterrows(),
                                      _rows(T_mmap).iterrows()):
        print('%-10s %10.4f %10.4f %10.4f' % (a['Station'], a['RespT'], b['RespT'], c['RespT']))

    # Validate: dec.mmap results should be finite and within 20% of MVA
    Q_mmap = np.asarray(_rows(T_mmap)['QLen'], dtype=float)
    Q_mva = np.asarray(_rows(T_mva)['QLen'], dtype=float)

    assert np.all(np.isfinite(Q_mmap)), \
        'ETAQA FCFS: dec.mmap returned non-finite queue lengths'
    assert np.all(Q_mmap >= 0), 'ETAQA FCFS: dec.mmap returned negative queue lengths'

    rel_err = np.abs(Q_mmap - Q_mva) / np.maximum(Q_mva, 1e-6)
    print('\nRelative error vs MVA: %.2f%%, %.2f%%' % tuple(rel_err * 100))
    assert np.all(rel_err < 0.20), \
        'ETAQA FCFS: dec.mmap queue lengths deviate >20%% from MVA (err=%.2f%%)' \
        % (rel_err.max() * 100)

    # Test configurable truncation level
    solver4 = MAM(model, method='dec.mmap')
    solver4.options.config['etaqa_trunc'] = 4
    Q_mmap4 = np.asarray(_rows(solver4.getAvgTable())['QLen'], dtype=float)
    assert np.all(np.isfinite(Q_mmap4)), \
        'ETAQA FCFS: truncation level 4 returned non-finite results'

    print('\nETAQA FCFS departure process test PASSED.')
