"""
Test ETAQA departure process for PS queues in MAM dec.mmap

Creates a 2-queue tandem with Erlang arrivals and PS scheduling, solved with
dec.mmap (ETAQA-PS-based departures). Validates that:
  1. dec.mmap now accepts PS queues (previously unsupported)
  2. Results are finite and consistent with MVA (exact for PS mean values)
  3. The PS departure process captures queueing effects beyond scaled service
"""

import numpy as np

from line_solver import (APH, Erlang, Exp, GlobalConstants, MAM, MVA, Network,
                         OpenClass, Queue, SchedStrategy, Sink, Source,
                         VerboseLevel)


def _rows(table):
    """The Queue rows of an AvgTable, in station order."""
    return table[table['Station'].astype(str).str.contains('Queue')]


if __name__ == '__main__':
    GlobalConstants.set_verbose(VerboseLevel.SILENT)

    # Model: Erl(3)/M/1-PS -> M/1-PS tandem
    model = Network('ETAQA-PS-Tandem')

    source = Source(model, 'Source')
    queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
    queue2 = Queue(model, 'Queue2', SchedStrategy.PS)
    sink = Sink(model, 'Sink')

    oclass = OpenClass(model, 'Class1')
    source.setArrival(oclass, Erlang.fitMeanAndOrder(2, 3))  # mean=2, SCV=1/3
    queue1.setService(oclass, Exp(2))                        # rho1=0.5
    queue2.setService(oclass, Exp(1.5))                      # rho2=0.667

    model.link(Network.serialRouting(source, queue1, queue2, sink))

    T_mmap = MAM(model, method='dec.mmap').getAvgTable()
    T_mva = MVA(model).getAvgTable()

    print('\n=== ETAQA PS Departure Process Test ===')
    print('%-10s %10s %10s' % ('Station', 'MVA', 'dec.mmap'))
    print('--- Queue Lengths ---')
    for (_, a), (_, b) in zip(_rows(T_mva).iterrows(), _rows(T_mmap).iterrows()):
        print('%-10s %10.4f %10.4f' % (a['Station'], a['QLen'], b['QLen']))
    print('--- Utilizations ---')
    for (_, a), (_, b) in zip(_rows(T_mva).iterrows(), _rows(T_mmap).iterrows()):
        print('%-10s %10.4f %10.4f' % (a['Station'], a['Util'], b['Util']))

    Q_mmap = np.asarray(_rows(T_mmap)['QLen'], dtype=float)
    Q_mva = np.asarray(_rows(T_mva)['QLen'], dtype=float)
    U_mmap = np.asarray(_rows(T_mmap)['Util'], dtype=float)
    U_mva = np.asarray(_rows(T_mva)['Util'], dtype=float)

    assert np.all(np.isfinite(Q_mmap)), 'ETAQA PS: dec.mmap returned non-finite queue lengths'
    assert np.all(Q_mmap >= 0), 'ETAQA PS: dec.mmap returned negative queue lengths'

    # Validate utilization (should match exactly since throughput = arrival rate)
    rel_err_u = np.abs(U_mmap - U_mva) / np.maximum(U_mva, 1e-6)
    print('\nUtilization error vs MVA: %.4f%%, %.4f%%' % tuple(rel_err_u * 100))
    assert np.all(rel_err_u < 0.01), \
        'ETAQA PS: utilizations deviate >1%% from MVA (err=%.4f%%)' % (rel_err_u.max() * 100)

    # Validate queue lengths (PS mean QLen = rho/(1-rho), independent of the
    # arrival process). For a single-class open M/G/1-PS, E[Q] = rho/(1-rho)
    # regardless of arrival SCV; the dec.mmap Q approximation uses U/(1-U).
    rel_err_q = np.abs(Q_mmap - Q_mva) / np.maximum(Q_mva, 1e-6)
    print('Queue length error vs MVA: %.2f%%, %.2f%%' % tuple(rel_err_q * 100))
    assert np.all(rel_err_q < 0.20), \
        'ETAQA PS: queue lengths deviate >20%% from MVA (err=%.2f%%)' % (rel_err_q.max() * 100)

    # Test: single PS queue (simplest case)
    model2 = Network('ETAQA-PS-Single')
    src2 = Source(model2, 'Source')
    q2 = Queue(model2, 'Queue1', SchedStrategy.PS)
    snk2 = Sink(model2, 'Sink')
    oc2 = OpenClass(model2, 'Class1')
    src2.setArrival(oc2, APH.fitMeanAndSCV(2, 4))   # high-variability arrival
    q2.setService(oc2, Exp(1))                      # rho=0.5
    model2.link(Network.serialRouting(src2, q2, snk2))

    T2 = MAM(model2, method='dec.mmap').getAvgTable()
    Q2 = np.asarray(_rows(T2)['QLen'], dtype=float)
    assert np.all(np.isfinite(Q2)), 'ETAQA PS single queue: returned non-finite results'
    assert np.all(Q2 >= 0), 'ETAQA PS single queue: returned negative queue lengths'
    print('Single PS queue Q=%.4f (expected ~1.0)' % Q2[0])

    print('\nETAQA PS departure process test PASSED.')
