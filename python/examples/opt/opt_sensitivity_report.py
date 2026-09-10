"""
opt_sensitivity_report  Read back the analytic d(metric)/d(service rate) table
that LineEvaluator attaches to every flat-model evaluation.

compute_model_sensitivities fills a three-level nested mapping (metric kind ->
metric key -> parameter key -> derivative). The optimizer consumes it
internally; this example walks it directly, which is also how one inspects which
station a rate change moves.

Native Python returns that mapping as a plain nested dict rather than MATLAB's
SensitivityData object, so the walk indexes it directly instead of calling
forKind; the three levels and their meaning are the same. The keys are tuples
here -- (station, class) for a metric, ('rate', station, class) for a parameter
-- where MATLAB joins them into strings, so they are formatted on the way out.
"""

from line_solver import (Exp, GlobalConstants, LineEvaluator, Network,
                         OpenClass, Queue, SchedStrategy, Sink, Source,
                         VerboseLevel)


def _key(k):
    """MATLAB joins a composite key with '||'; do the same for the report."""
    if isinstance(k, (tuple, list)):
        return '||'.join(str(p) for p in k)
    return str(k)


if __name__ == '__main__':
    GlobalConstants.set_verbose(VerboseLevel.SILENT)

    model = Network('Tandem')
    source = Source(model, 'Arrivals')
    q1 = Queue(model, 'Frontend', SchedStrategy.FCFS)
    q2 = Queue(model, 'Backend', SchedStrategy.FCFS)
    sink = Sink(model, 'Departures')
    jobs = OpenClass(model, 'Jobs')
    source.setArrival(jobs, Exp(1.0))
    q1.setService(jobs, Exp(2.0))
    q2.setService(jobs, Exp(1.5))
    model.link(Network.serialRouting(source, q1, q2, sink))

    evaluator = LineEvaluator(model, [], [])
    result = evaluator.evaluateValues({})

    print('Solver used : %s' % result.solverUsed)
    print('Feasible    : %d' % int(result.feasible))
    print('Frontend Q  : %.4f' % result.getQueueLength('Frontend'))
    print('Backend  Q  : %.4f' % result.getQueueLength('Backend'))

    sens = result.sensitivities
    if not sens:
        print('\nNo analytic sensitivities available for this model.')
        raise SystemExit(0)

    for kind in ('QLen', 'RespT', 'Util', 'Tput'):
        by_metric = sens.get(kind)
        if not by_metric:
            continue
        print('\n=== d(%s)/d(rate) ===' % kind)
        # sort: a dict preserves insertion order, so sort for a stable report
        for mk in sorted(by_metric, key=_key):
            by_param = by_metric[mk]
            for pk in sorted(by_param, key=_key):
                print('  %-24s %-32s % .6f' % (_key(mk), _key(pk), by_param[pk]))
