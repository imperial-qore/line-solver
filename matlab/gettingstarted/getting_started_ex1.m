% https://github.com/imperial-qore/line-solver/wiki/Getting-started#example-1-a-mm1-queue
model = Network('M/M/1');
%% Block 1: nodes
source = Source(model, 'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');
%% Block 2: classes
jobclass = OpenClass(model, 'Class1');
source.setArrival(jobclass, Exp(1));
queue.setService(jobclass, Exp(2));
%% Block 3: topology
model.link(Network.serialRouting(source,queue,sink));
%% Block 4: solution
AvgTable = SolverJMT(model,'seed',23000,'samples',10000).getAvgTable
%% select a particular table row
ARow = tget(AvgTable, queue, jobclass) % this is also valid
%% select a particular table row by node and class label
ARow = tget(AvgTable, 'Queue', 'Class1')
%% export to JMT
%model.jsimgView
