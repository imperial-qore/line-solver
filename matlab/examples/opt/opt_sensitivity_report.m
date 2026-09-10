% opt_sensitivity_report  Read back the analytic d(metric)/d(service rate)
% table that opt.LineEvaluator attaches to every flat-model evaluation.
%
% opt.sens.computeModelSensitivities fills an opt.SensitivityData, a three-level
% nested dictionary (metric kind -> metric key -> parameter key -> derivative).
% The optimizer consumes it internally; this example walks it directly via
% forKind, which is also how one inspects which station a rate change moves.

model = Network('Tandem');
source = Source(model, 'Arrivals');
q1 = Queue(model, 'Frontend', SchedStrategy.FCFS);
q2 = Queue(model, 'Backend', SchedStrategy.FCFS);
sink = Sink(model, 'Departures');
jobs = OpenClass(model, 'Jobs');
source.setArrival(jobs, Exp(1.0));
q1.setService(jobs, Exp(2.0));
q2.setService(jobs, Exp(1.5));
model.link(Network.serialRouting(source, q1, q2, sink));

evaluator = opt.LineEvaluator(model, {}, {});
result = evaluator.evaluateValues(configureDictionary('string', 'cell'));

fprintf('Solver used : %s\n', result.solverUsed);
fprintf('Feasible    : %d\n', result.feasible);
fprintf('Frontend Q  : %.4f\n', result.getQueueLength('Frontend'));
fprintf('Backend  Q  : %.4f\n', result.getQueueLength('Backend'));

sens = result.sensitivities;
if isempty(sens)
    fprintf('\nNo analytic sensitivities available for this model.\n');
    return
end

for kind = {'QLen', 'RespT', 'Util', 'Tput'}
    byMetric = sens.forKind(kind{1});
    if isempty(byMetric)
        continue
    end
    fprintf('\n=== d(%s)/d(rate) ===\n', kind{1});
    % sort: a dictionary preserves insertion order, so sort for a stable report
    mk = sort(keys(byMetric));
    for i = 1:numel(mk)
        byParam = byMetric{mk(i)};
        pk = sort(keys(byParam));
        for j = 1:numel(pk)
            fprintf('  %-24s %-32s % .6f\n', mk(i), pk(j), byParam(pk(j)));
        end
    end
end
