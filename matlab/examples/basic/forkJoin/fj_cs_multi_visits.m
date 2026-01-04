% fork-join with multiple visits within the same chain
clear solver AvgTable
model = Network('model');

source = Source(model,'Source');
queue1 = Queue(model,'Queue1',SchedStrategy.PS);
queue2 = Queue(model,'Queue2',SchedStrategy.PS);
fork = Fork(model,'Fork');
join = Join(model,'Join', fork);
sink = Sink(model,'Sink');

jobclass1 = OpenClass(model, 'class1');
jobclass2 = OpenClass(model, 'class2');

source.setArrival(jobclass1, Exp(0.1));
queue1.setService(jobclass1, Exp(1.0));
queue2.setService(jobclass1, Exp(1.0));
queue1.setService(jobclass2, Exp(1.0));
queue2.setService(jobclass2, Exp(1.0));

P = model.initRoutingMatrix;
P{jobclass1,jobclass1}(source,fork) = 1;
P{jobclass1,jobclass1}(fork,queue1) = 1.0;
P{jobclass1,jobclass1}(fork,queue2) = 1.0;
P{jobclass1,jobclass1}(queue1,join) = 1.0;
P{jobclass1,jobclass1}(queue2,join) = 1.0;

% now loop back in class 2
P{jobclass1,jobclass2}(join,fork) = 1.0;
P{jobclass2,jobclass2}(fork,queue1) = 1.0;
P{jobclass2,jobclass2}(fork,queue2) = 1.0;
P{jobclass2,jobclass2}(queue1,join) = 1.0;
P{jobclass2,jobclass2}(queue2,join) = 1.0;
P{jobclass2,jobclass2}(join,sink) = 1.0;

model.link(P);

solver{1} = MVA(model);
solver{end+1} = JMT(model,'seed',23000);
%solver{end+1} = DES(model,'seed',23000);

for s=1:length(solver)
    AvgTable{s} = solver{s}.getAvgTable;
    AvgTable{s}
end
