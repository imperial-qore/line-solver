function model = gallery_mhypk(k)
if ~exist('k','var')
    k = 2;
end
model = Network(['M/Hyper/',num2str(k)]);
%% Block 1: nodes
source = Source(model, 'mySource');
queue = Queue(model, 'myQueue', SchedStrategy.FCFS);
sink = Sink(model, 'mySink');
%% Block 2: classes
oclass = OpenClass(model, 'myClass');
source.setArrival(oclass, Exp(1));
queue.setService(oclass, Coxian.fitMeanAndSCV(0.5,4));
queue.setNumberOfServers(k);
%% Block 3: topology
model.link(Network.serialRouting(source,queue,sink));
end