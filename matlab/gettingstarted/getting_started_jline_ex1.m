jmodel = JNetwork('M/M/1');
%% Block 1: nodes
source = Source(jmodel, 'mySource');
queue = Queue(jmodel, 'myQueue', SchedStrategy.FCFS);
sink = Sink(jmodel, 'mySink');
%% Block 2: classes
oclass = OpenClass(jmodel, 'myClass');
source.setArrival(oclass, Exp(1));
queue.setService(oclass, Exp(2));
%% Block 3: topology
jmodel.link(Network.serialRouting(source,queue,sink));
%% Block 4: solution
AvgTable = SolverMVA(jmodel).getAvgTable
%% select a particular table row
%ARow = tget(AvgTable, queue, oclass) % this is also valid
%% select a particular table row by node and class label
%ARow = tget(AvgTable, 'myQueue', 'myClass')
%% export to JMT
%model.jsimgView
