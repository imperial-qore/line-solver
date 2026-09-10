"""
Test ETAQA PS departure process: MAM dec.mmap vs LDES

Erl(3)/M/1-PS -> M/1-PS tandem
Compares MAM (dec.mmap with ETAQA PS departures) against LDES simulation.
"""

from line_solver import (Erlang, Exp, GlobalConstants, LDES, MAM, MVA, Network,
                         OpenClass, Queue, SchedStrategy, Sink, Source,
                         VerboseLevel)


def _rows(table):
    return table[table['Station'].astype(str).str.contains('Queue')]


if __name__ == '__main__':
    GlobalConstants.set_verbose(VerboseLevel.SILENT)

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
    T_ldes = LDES(model, seed=23000, samples=int(5e5)).getAvgTable()

    print('\n=== ETAQA PS: MAM vs LDES ===')
    print('%-10s %10s %10s %10s' % ('Station', 'LDES', 'MVA', 'dec.mmap'))
    for label in ('QLen', 'RespT', 'Util'):
        print('--- %s ---' % {'QLen': 'Queue Lengths', 'RespT': 'Response Times',
                              'Util': 'Utilizations'}[label])
        for (_, l), (_, v), (_, m) in zip(_rows(T_ldes).iterrows(),
                                          _rows(T_mva).iterrows(),
                                          _rows(T_mmap).iterrows()):
            print('%-10s %10.4f %10.4f %10.4f'
                  % (l['Station'], l[label], v[label], m[label]))
