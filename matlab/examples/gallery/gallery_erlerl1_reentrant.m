function model = gallery_erlerl1_reentrant
model = Network('Erl/Erl/1-Reentrant');
n=5;
%% Block 1: nodes
source = Source(model, 'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');
%% Block 2: classes
oclass1 = OpenClass(model, 'Class1');
oclass2 = OpenClass(model, 'Class2');
source.setArrival(oclass1, Erlang.fitMeanAndOrder(1,n));
source.setArrival(oclass2, Disabled());
queue.setService(oclass1, Erlang.fitMeanAndOrder(0.5,n));
queue.setService(oclass2, Exp(3));
%% Block 3: topology
P = model.initRoutingMatrix;
P{oclass1,oclass1}(source,queue)=1;
P{oclass1,oclass2}(queue,queue)=0.5;
P{oclass2,oclass2}(queue,sink)=0.5;
model.link(P);
end