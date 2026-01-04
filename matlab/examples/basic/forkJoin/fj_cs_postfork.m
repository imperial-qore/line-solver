clear solver AvgTable
model = Network('model');

delay = Delay(model,'Delay');
fork1 = Fork(model,'Fork1');
join1 = Join(model,'Join1',fork1);
queue1 = Queue(model,'Queue1',SchedStrategy.PS);
queue2 = Queue(model,'Queue2',SchedStrategy.PS);

jobclass1 = ClosedClass(model, 'class1', 1, delay);
jobclass2 = ClosedClass(model, 'class2', 1, delay);

delay.setService(jobclass1, Exp(0.25));
queue1.setService(jobclass1, Exp(2.0));
queue2.setService(jobclass1, Exp(2.0));

delay.setService(jobclass2, Exp(0.25));
queue1.setService(jobclass2, Exp(2.0));
queue2.setService(jobclass2, Exp(2.0));

M = model.getNumberOfNodes;
R = model.getNumberOfClasses;

P = model.initRoutingMatrix;
P{jobclass1,jobclass1}(delay,fork1) = 1.0;
P{jobclass1,jobclass1}(fork1,queue1) = 1.0;
P{jobclass1,jobclass1}(fork1,queue2) = 1.0;
P{jobclass1,jobclass1}(queue1,join1) = 1.0;
P{jobclass1,jobclass1}(queue2,join1) = 1.0;
P{jobclass1,jobclass1}(join1,delay) = 1.0;

P{jobclass2,jobclass2}(delay,fork1) = 1.0;
P{jobclass2,jobclass1}(fork1,queue1) = 1.0;
P{jobclass2,jobclass1}(fork1,queue2) = 1.0;

model.link(P);

solver{1} = MVA(model);
solver{end+1} = JMT(model,'seed',23000); % JMT has a bug on this one
%solver{end+1} = DES(model,'seed',23000);

for s=1:length(solver)
    AvgTable{s} = solver{s}.getAvgTable;
    AvgTable{s}
end
