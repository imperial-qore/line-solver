"""
Test ETAQA FCFS departure process: MAM dec.mmap vs LDES

Erl(5)/Erl(2)/1 -> Erl(3)/1 FCFS tandem
Compares MAM (dec.mmap with ETAQA departures) against LDES simulation.
"""

from line_solver import (Erlang, GlobalConstants, LDES, MAM, MVA, Network,
                         OpenClass, Queue, SchedStrategy, Sink, Source,
                         VerboseLevel)


def _rows(table):
    return table[table['Station'].astype(str).str.contains('Queue')]


if __name__ == '__main__':
    GlobalConstants.set_verbose(VerboseLevel.SILENT)

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

    T_mmap = MAM(model, method='dec.mmap').getAvgTable()
    T_src = MAM(model, method='dec.source').getAvgTable()
    T_mva = MVA(model).getAvgTable()
    T_ldes = LDES(model, seed=23000, samples=int(5e5)).getAvgTable()

    print('\n=== ETAQA FCFS: MAM vs LDES ===')
    print('%-10s %10s %10s %10s %10s'
          % ('Station', 'LDES', 'MVA', 'dec.source', 'dec.mmap'))
    for label in ('QLen', 'RespT', 'Util'):
        print('--- %s ---' % {'QLen': 'Queue Lengths', 'RespT': 'Response Times',
                              'Util': 'Utilizations'}[label])
        for (_, l), (_, v), (_, s), (_, m) in zip(_rows(T_ldes).iterrows(),
                                                  _rows(T_mva).iterrows(),
                                                  _rows(T_src).iterrows(),
                                                  _rows(T_mmap).iterrows()):
            print('%-10s %10.4f %10.4f %10.4f %10.4f'
                  % (l['Station'], l[label], v[label], s[label], m[label]))
