function model = gallery_parm1
model = Network('Par/M/1');
%% Block 1: nodes
source = Source(model, 'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');
%% Block 2: classes
oclass1 = OpenClass(model, 'Class1');
source.setArrival(oclass1, Pareto.fitMeanAndSCV(1,64));
queue.setService(oclass1, Exp(2));
%% Block 3: topology
P = model.initRoutingMatrix;
P{oclass1,oclass1}(source,queue)=1;
P{oclass1,oclass1}(queue,sink)=1;
model.link(P);
end