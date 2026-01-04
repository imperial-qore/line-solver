function model = gallery_mdk(k)
if nargin<1
    k = 2;
end
model = Network('M/D/k');
%% Block 1: nodes
source = Source(model, 'mySource');
queue = Queue(model, 'myQueue', SchedStrategy.FCFS);
sink = Sink(model, 'mySink');
%% Block 2: classes
oclass = OpenClass(model, 'myClass');
source.setArrival(oclass, Exp.fitMean(1));
queue.setService(oclass, Det(2/k));
queue.setNumberOfServers(k);
%% Block 3: topology
model.link(Network.serialRouting(source,queue,sink));
end