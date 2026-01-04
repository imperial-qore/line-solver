clear solver AvgTable;

%clearvars -except exampleName;
model = Network('model');

% source = Source(model,'Source');
delay = Delay(model, 'Delay1');
queue1 = Queue(model,'Queue1',SchedStrategy.FCFS);
queue2 = Queue(model,'Queue2',SchedStrategy.FCFS);
queue3 = Queue(model,'Queue3',SchedStrategy.FCFS);
queue4 = Queue(model,'Queue4',SchedStrategy.FCFS);
queue5 = Queue(model,'Queue5',SchedStrategy.FCFS);
fork = Fork(model,'Fork');
join = Join(model,'Join', fork);
% sink = Sink(model,'Sink');

jobclass1 = ClosedClass(model, 'class1', 10, delay, 0);

% source.setArrival(jobclass1, Exp(0.5));
queue1.setService(jobclass1, Exp(1.0));
queue2.setService(jobclass1, Exp(2.0));
queue3.setService(jobclass1, Exp(1.0));
queue4.setService(jobclass1, Exp(3.0));
queue5.setService(jobclass1, Exp(0.8));
delay.setService(jobclass1, Exp(0.5));

P = zeros(5);
% P(source,fork) = 1;
P(delay, fork) = 1.0;
P(fork,queue1) = 1.0;
P(fork,queue2) = 1.0;
P(queue1,queue4) = 1.0;
P(queue4,queue5) = 1.0;
P(queue5,join) = 1.0;
P(queue2,queue3) = 1.0;
P(queue3,join) = 1.0;
P(join,delay) = 1.0;

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
