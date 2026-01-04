clear solver AvgTable;

%clearvars -except exampleName;
model = Network('model');

% source = Source(model,'Source');
delay = Delay(model, 'Delay1');
queue1 = Queue(model,'Queue1',SchedStrategy.FCFS);
queue2 = Queue(model,'Queue2',SchedStrategy.FCFS);
fork = Fork(model,'Fork');
join = Join(model,'Join', fork);
queue3 = Queue(model,'Queue3',SchedStrategy.FCFS);
queue4 = Queue(model,'Queue4',SchedStrategy.FCFS);
fork2 = Fork(model,'Fork2');
join2 = Join(model,'Join2', fork2);
% sink = Sink(model,'Sink');

jobclass1 = ClosedClass(model, 'class1', 1, delay, 0);

% source.setArrival(jobclass1, Exp(0.5));
queue1.setService(jobclass1, Exp(1.0));
queue2.setService(jobclass1, Exp(1.0));
delay.setService(jobclass1, Exp(0.5));
queue3.setService(jobclass1, Exp(2.0));
queue4.setService(jobclass1, Exp(2.0));

P = zeros(9);
% P(source,fork) = 1;
P(delay, fork) = 1.0;
P(fork,queue1) = 1.0;
P(fork,queue2) = 1.0;
P(queue1,fork2) = 1.0;
P(fork2, queue3) = 1.0;
P(fork2, queue4) = 1.0;
P(queue3, join2) = 1.0;
P(queue4, join2) = 1.0;
P(join2, join) = 1.0;
P(queue2,join) = 1.0;
P(join,delay) = 1.0;

model.link(P);

solver = {};
solver{end+1} = JMT(model,'seed',23000);
mva_options = MVA.defaultOptions;
mva_options.config.fork_join = 'fjt';
solver{end+1} = MVA(model, mva_options);
%solver{end+1} = DES(model,'seed',23000);

AvgTable = {};
for s=1:length(solver)
    AvgTable{end+1} = solver{s}.getAvgTable;
    AvgTable{s}
end
