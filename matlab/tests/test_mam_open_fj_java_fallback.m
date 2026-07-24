% Regression test for Java-backed SolverMAM dec.source.mmap on an open
% fork-join model.

model = Network('fj_basic_open_regression');

source = Source(model, 'Source');
queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS);
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS);
fork = Fork(model, 'Fork');
join = Join(model, 'Join', fork);
sink = Sink(model, 'Sink');

jobclass = OpenClass(model, 'class1');

source.setArrival(jobclass, Exp(0.05));
queue1.setService(jobclass, Exp(1.0));
queue2.setService(jobclass, Exp(2.0));

P = zeros(6);
P(source, fork) = 1;
P(fork, queue1) = 1;
P(fork, queue2) = 1;
P(queue1, join) = 1;
P(queue2, join) = 1;
P(join, sink) = 1;
model.link(P);

options = SolverMAM.defaultOptions;
options.lang = 'java';
options.method = 'dec.source.mmap';
options.verbose = VerboseLevel.SILENT;

solver = SolverMAM(model, options);
avgTable = solver.getAvgTable;
joinRow = strcmp(string(avgTable.Station), "Join") & strcmp(string(avgTable.JobClass), "class1");

assert(~isempty(avgTable), ...
    'SolverMAM did not return results for dec.source.mmap open fork-join model.');
assert(isfield(solver.result, 'Avg') && isfield(solver.result.Avg, 'T') && ...
    ~isempty(solver.result.Avg.T) && any(solver.result.Avg.T(:) > 0), ...
    'SolverMAM returned empty or zero throughput for dec.source.mmap open fork-join model.');
assert(any(joinRow), ...
    'SolverMAM did not report a Join row for dec.source.mmap open fork-join model.');
assert(avgTable.Tput(joinRow) > 0, ...
    'SolverMAM returned zero throughput at the Join for dec.source.mmap open fork-join model.');
assert(avgTable.RespT(joinRow) > 0, ...
    'SolverMAM returned zero synchronization delay at the Join for dec.source.mmap open fork-join model.');
