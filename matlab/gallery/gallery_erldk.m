function model = gallery_erldk(k)
if nargin<1
    k = 2;
end
model = Network('Erl/D/k');
%% Block 1: nodes
source = Source(model, 'mySource');
queue = Queue(model, 'myQueue', SchedStrategy.FCFS);
sink = Sink(model, 'mySink');
%% Block 2: classes
oclass = OpenClass(model, 'myClass');
source.setArrival(oclass, Erlang.fitMeanAndOrder(1,5));
queue.setService(oclass, Det(2/k));
queue.setNumberOfServers(k);
%% Block 3: topology
model.link(Network.serialRouting(source,queue,sink));
end