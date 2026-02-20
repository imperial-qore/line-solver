clear node jobclass solver AvgTable;

model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.GPSPRIO);

jobclass{1} = ClosedClass(model, 'Class1', 6, node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', 4, node{1}, 1);
jobclass{3} = ClosedClass(model, 'Class3', 4, node{1}, 1);
jobclass{4} = ClosedClass(model, 'Class4', 1, node{1}, 2);

node{1}.setService(jobclass{1}, Erlang(3.0,2));
node{1}.setService(jobclass{2}, Exp(1.0));
node{1}.setService(jobclass{3}, Exp(1.0));
node{1}.setService(jobclass{4}, Exp(2.0));

w1=12; node{2}.setService(jobclass{1}, Exp(30), w1);
w2=3; node{2}.setService(jobclass{2}, Exp(2), w2);
w3=5; node{2}.setService(jobclass{3}, Exp(12), w3);
w4=1; node{2}.setService(jobclass{4}, Exp(1), w4);

M = model.getNumberOfStations();
K = model.getNumberOfClasses();

P = model.initRoutingMatrix;
P{1} = Network.serialRouting(node);
P{2} = Network.serialRouting(node);
P{3} = Network.serialRouting(node);
P{4} = Network.serialRouting(node);
model.link(P);
%%
% This part illustrates the execution of different solvers
solver = {};
solver{end+1} = CTMC(model);
solver{end+1} = JMT(model,'seed',23000,'samples',3e4,'keep',true);
%solver{end+1} = SSA(model,'seed',23000,'samples',3e4);
% solver{end+1} = FLD(model);
% solver{end+1} = MVA(model,'exact');
% solver{end+1} = NC(model,'exact');
% solver{end+1} = MAM(model);
% solver{end+1} = LINE(model);
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',strrep(solver{s}.getName(),'Solver',''));
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
