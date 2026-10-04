% This model examplifies how to specify models with multiple sinks (virtual
% sink)
model = Network('model');

source = Source(model, 'Source');
queue = QueueingStation(model, 'Queue1', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');
vsink1 = Router(model, 'VSink1');
vsink2 = Router(model, 'VSink2');

ocl1 = OpenClass(model, 'Class1');
ocl2 = OpenClass(model, 'Class2');

source.setArrival(ocl1,Exp(1.0));
queue.setService(ocl1,Exp(100.0));

source.setArrival(ocl2,Exp(1.0));
queue.setService(ocl2,Exp(100.0));

P = model.initRoutingMatrix;

P{ocl1}(source,queue) = 1.0;
P{ocl1}(queue,vsink1) = 0.6;
P{ocl1}(queue,vsink2) = 0.4;
P{ocl1}(vsink1,sink) = 1.0;
P{ocl1}(vsink2,sink) = 1.0;

P{ocl2}(source,queue) = 1.0;
P{ocl2}(queue,vsink1) = 0.1;
P{ocl2}(queue,vsink2) = 0.9;
P{ocl2}(vsink1,sink) = 1.0;
P{ocl2}(vsink2,sink) = 1.0;

model.link(P);

% We use getAvgNodeTable to see the throuhgputs of sink1 and sink2
AvgTable = MVA(model).getAvgTable
AvgNodeTable = MVA(model).getAvgNodeTable

AvgTableMAM = MAM(model).getAvgTable
AvgNodeTableMAM = MAM(model).getAvgNodeTable

AvgTableNC = NC(model).getAvgTable
AvgNodeTableNC = NC(model).getAvgNodeTable
