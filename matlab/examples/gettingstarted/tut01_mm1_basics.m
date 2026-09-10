% Example 1: A M/M/1 queue
GlobalConstants.setVerbose(VerboseLevel.STD);
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
AvgTable = JMT(model,'seed',23000,'samples',10000).avgTable();
%% select a particular table row using direct indexing
ARow = AvgTable(queue, jobclass)
%% select a particular table row using string names
ARow = AvgTable.filterBy(queue)  % returns all rows for 'Queue'
%% export to JMT
%model.jsimgView
