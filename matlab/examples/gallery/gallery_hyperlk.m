function model = gallery_hyperlk(k)
if nargin<1
    k = 2;
end
model = Network('Hyper/Erl/k');
%% Block 1: nodes
source = Source(model, 'mySource');
queue = Queue(model, 'myQueue', SchedStrategy.FCFS);
sink = Sink(model, 'mySink');
%% Block 2: classes
oclass = OpenClass(model, 'myClass');
source.setArrival(oclass, HyperExp.fitMeanAndSCVBalanced(1/1.8,4));
queue.setService(oclass, Erlang.fitMeanAndSCV(1,0.25));
queue.setNumberOfServers(k);
%% Block 3: topology
model.link(Network.serialRouting(source,queue,sink));
end