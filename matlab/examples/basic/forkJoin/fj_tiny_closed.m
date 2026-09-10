% Minimal closed fork-join model solved exactly by the native CTMC
% fork-join implementation (tag-augmented state space) and validated
% against the closed form for N=1:
%   X = 1/(Z + E[max(S1,S2)]),  E[max] = 1/mu1 + 1/mu2 - 1/(mu1+mu2)
clear solver AvgTable;
model = Network('model');

delay = Delay(model,'Delay1');
queue1 = Queue(model,'Queue1',SchedStrategy.FCFS);
queue2 = Queue(model,'Queue2',SchedStrategy.FCFS);
fork = Fork(model,'Fork1');
join = Join(model,'Join1',fork);

jobclass1 = ClosedClass(model,'Class1',1,delay);

delay.setService(jobclass1, Exp(1.0));
queue1.setService(jobclass1, Exp(2.0));
queue2.setService(jobclass1, Exp(3.0));

P = zeros(5);
P(delay,fork) = 1;
P(fork,queue1) = 1.0;
P(fork,queue2) = 1.0;
P(queue1,join) = 1.0;
P(queue2,join) = 1.0;
P(join,delay) = 1.0;

model.link(P);

solver = {};
solver{end+1} = CTMC(model,'exact');
solver{end+1} = SSA(model,'seed',23000,'samples',5e4);
solver{end+1} = JMT(model,'seed',23000);

AvgTable = {};
for s=1:length(solver)
    AvgTable{end+1} = solver{s}.getAvgTable;
    AvgTable{s}
end

% closed-form check (N=1): X = 1/(1 + 1/2 + 1/3 - 1/5) = 0.612244...
