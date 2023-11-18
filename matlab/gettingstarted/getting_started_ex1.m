% https://github.com/imperial-qore/line-solver/wiki/Getting-started#example-1-a-mm1-queue
model = Network('M/M/1');
%% Block 1: nodes
source = Source(model, 'mySource');
queue = Queue(model, 'myQueue', SchedStrategy.FCFS);
sink = Sink(model, 'mySink');
%% Block 2: classes
oclass = OpenClass(model, 'myClass');
source.setArrival(oclass, Exp(1));
queue.setService(oclass, Exp(2));
%% Block 3: topology
model.link(Network.serialRouting(source,queue,sink));
%% Block 4: solution
AvgTable = SolverJMT(model,'seed',23000,'samples',10000).getAvgTable
%% select a particular table row
ARow = tget(AvgTable, queue, oclass) % this is also valid
%% select a particular table row by node and class label
ARow = tget(AvgTable, 'myQueue', 'myClass')
%% export to JMT
%model.jsimgView
