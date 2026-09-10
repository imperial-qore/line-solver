function model = gallery_replayerm1(filename)
if ~exist('filename','var')
    filename = which('example_trace.txt');
end
model = Network('Trace/M/1');
%% Block 1: nodes
source = Source(model, 'mySource');
queue = Queue(model, 'myQueue', SchedStrategy.FCFS);
sink = Sink(model, 'mySink');
%% Block 2: classes
oclass = OpenClass(model, 'myClass');
replayer = Replayer(filename);
source.setArrival(oclass, replayer);
queue.setService(oclass, Exp(3/replayer.getMean));
%% Block 3: topology
model.link(Network.serialRouting(source,queue,sink));
end