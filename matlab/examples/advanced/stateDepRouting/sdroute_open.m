clear node jobclass solver AvgTable

model = Network('myModel');

% Block 1: nodes
source = Source(model, 'Source');
router = Router(model, 'Router');
queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS);
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');

% Block 2: classes
oclass = OpenClass(model, 'Class1');
source.setArrival(oclass, Exp(1));
queue1.setService(oclass, Exp(2));
queue2.setService(oclass, Exp(2));

% Block 3: topology
model.addLink(source, router);
model.addLink(router, queue1);
model.addLink(router, queue2);
model.addLink(queue1, sink);
model.addLink(queue2, sink);

router.setRouting(oclass, RoutingStrategy.RROBIN);

solver = {};
solver{1} = JMT(model,'seed',23000);
solver{2} = CTMC(model,'cutoff',5);

AvgTable = {};
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());
    AvgTable{s} = solver{s}.getAvgNodeTable();
    AvgTable{s}
end
