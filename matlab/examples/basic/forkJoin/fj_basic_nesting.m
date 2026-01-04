clear solver AvgTable;

model = Network('model');

delay = Delay(model,'Delay');
fork1 = Fork(model,'Fork1');
fork1.setTasksPerLink(1);
fork11 = Fork(model,'Fork1_1');
fork11.setTasksPerLink(2);
join1 = Join(model,'Join1',fork1);
join11 = Join(model,'Join1_1',fork11);
queue1 = Queue(model,'Queue1',SchedStrategy.PS);
queue2 = Queue(model,'Queue2',SchedStrategy.PS);

jobclass1 = ClosedClass(model, 'class1', 5, delay);
jobclass2 = ClosedClass(model, 'class2', 2, delay);

delay.setService(jobclass1, Exp(0.25));
queue1.setService(jobclass1, Exp(1.0));
queue2.setService(jobclass1, Exp(0.75));

delay.setService(jobclass2, Exp(0.25));
queue1.setService(jobclass2, Exp(2.0));
queue2.setService(jobclass2, Exp(2.0));

M = model.getNumberOfNodes;
R = model.getNumberOfClasses;

P = cellzeros(R,M);
P{jobclass1,jobclass1}(delay,fork1) = 1.0;
P{jobclass1,jobclass1}(fork1,queue1) = 1.0;
P{jobclass1,jobclass1}(fork1,queue2) = 1.0;
P{jobclass1,jobclass1}(queue1,join1) = 1.0;
P{jobclass1,jobclass1}(queue2,join1) = 1.0;
P{jobclass1,jobclass1}(join1,delay) = 1.0;

P{jobclass2,jobclass2}(delay,fork11) = 1.0;
P{jobclass2,jobclass2}(fork11,fork1) = 1.0;
P{jobclass2,jobclass2}(fork1,queue1) = 1.0;
P{jobclass2,jobclass2}(fork1,queue2) = 1.0;
P{jobclass2,jobclass2}(queue1,join1) = 1.0;
P{jobclass2,jobclass2}(queue2,join1) = 1.0;
P{jobclass2,jobclass2}(join1,join11) = 1.0;
P{jobclass2,jobclass2}(join11,delay) = 1.0;

model.link(P);


solver = {};
solver{end+1} = JMT(model,'seed',23000);
solver{end+1} = MVA(model);
%solver{end+1} = DES(model,'seed',23000);

AvgTable = {};
for s=1:length(solver)
    AvgTable{end+1} = solver{s}.getAvgTable;
    AvgTable{end}
end

%% Corresponding open model
% model = Network('model');
%
% source = Source(model,'Source');
% queue1 = Queue(model,'Queue1',SchedStrategy.PS);
% queue2 = Queue(model,'Queue2',SchedStrategy.PS);
% fork1 = Fork(model,'Fork1');
% fork1.setTasksPerLink(1);
% fork11 = Fork(model,'Fork1_1');
% fork11.setTasksPerLink(2);
% join1 = Join(model,'Join1');
% join11 = Join(model,'Join1_1');
% sink = Sink(model,'Sink');
%
% jobclass1 = OpenClass(model, 'class1');
% jobclass2 = OpenClass(model, 'class2');
%
% source.setArrival(jobclass1, Exp(0.25));
% queue1.setService(jobclass1, Exp(1.0));
% queue2.setService(jobclass1, Exp(0.75));
%
% source.setArrival(jobclass2, Exp(0.25));
% queue1.setService(jobclass2, Exp(2.0));
% queue2.setService(jobclass2, Exp(2.0));
%
% M = model.getNumberOfNodes;
% R = model.getNumberOfClasses;
% P = cellzeros(R,M);
% P{jobclass1,jobclass1}(source,fork1) = 1;
% P{jobclass1,jobclass1}(fork1,queue1) = 1.0;
% P{jobclass1,jobclass1}(fork1,queue2) = 1.0;
% P{jobclass1,jobclass1}(queue1,join1) = 1.0;
% P{jobclass1,jobclass1}(queue2,join1) = 1.0;
% P{jobclass1,jobclass1}(join1,sink) = 1.0;
%
% P{jobclass2,jobclass2}(source,fork11) = 1;
% P{jobclass2,jobclass2}(fork11,fork1) = 1.0;
% P{jobclass2,jobclass2}(fork1,queue1) = 1.0;
% P{jobclass2,jobclass2}(fork1,queue2) = 1.0;
% P{jobclass2,jobclass2}(queue1,join1) = 1.0;
% P{jobclass2,jobclass2}(queue2,join1) = 1.0;
% P{jobclass2,jobclass2}(join1,join11) = 1.0;
% P{jobclass2,jobclass2}(join11,sink) = 1.0;
%
% model.link(P);
% JMT(model,'keep',true).getAvgNodeTable