clear solver AvgTable;

%clearvars -except exampleName;
model = Network('model');

% source = Source(model,'Source');
delay = Delay(model, 'Delay1');
queue1 = Queue(model,'Queue1',SchedStrategy.FCFS);
queue2 = Queue(model,'Queue2',SchedStrategy.FCFS);
fork = Fork(model,'Fork');
join = Join(model,'Join', fork);
% sink = Sink(model,'Sink');

jobclass1 = ClosedClass(model, 'class1', 10, delay, 0);
jobclass2 = ClosedClass(model, 'class2', 10, delay, 0);

% source.setArrival(jobclass1, Exp(0.5));
queue1.setService(jobclass1, Exp(1.0));
queue2.setService(jobclass1, Exp(2.0));
delay.setService(jobclass1, Exp(0.5));

queue1.setService(jobclass2, Exp(1.0));
queue2.setService(jobclass2, Exp(2.0));
delay.setService(jobclass2, Exp(0.2));

P = model.initRoutingMatrix;
% P(source,fork) = 1;
P{jobclass1, jobclass1}(delay, fork) = 1.0;
P{jobclass1, jobclass1}(fork,queue1) = 1.0;
P{jobclass1, jobclass1}(fork,queue2) = 1.0;
P{jobclass1, jobclass1}(queue1,join) = 1.0;
P{jobclass1, jobclass1}(queue2,join) = 1.0;
P{jobclass1, jobclass1}(join,delay) = 1.0;

P{jobclass2, jobclass2}(delay, fork) = 1.0;
P{jobclass2, jobclass2}(fork,queue1) = 1.0;
P{jobclass2, jobclass2}(fork,queue2) = 1.0;
P{jobclass2, jobclass2}(queue1,join) = 1.0;
P{jobclass2, jobclass2}(queue2,join) = 1.0;
P{jobclass2, jobclass2}(join,delay) = 1.0;

model.link(P);

solver = {};
solver{end+1} = JMT(model,'seed',23000);
solver{end+1} = MVA(model);
%solver{end+1} = DES(model,'seed',23000);

AvgTable = {};
for s=1:length(solver)
    AvgTable{end+1} = solver{s}.getAvgTable;
    AvgTable{s}
end
