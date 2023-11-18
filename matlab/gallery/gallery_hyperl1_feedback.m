function model = gallery_hyperl1_feedback
model = Network('Hyper/Erl/1-Feedback');
%% Block 1: nodes
source = Source(model, 'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');
%% Block 2: classes
oclass1 = OpenClass(model, 'Class1');
source.setArrival(oclass1, HyperExp.fitMeanAndSCV(1,64));
queue.setService(oclass1, Erlang.fitMeanAndOrder(0.5,5));
%% Block 3: topology
P = model.initRoutingMatrix;
P{oclass1,oclass1}(source,queue)=1;
P{oclass1,oclass1}(queue,queue)=0.9;
P{oclass1,oclass1}(queue,sink)=0.1;
model.link(P);
end