clear node jobclass solver AvgTable AvgChainTable AvgSysChainTable;

model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
node{2}.setNumberOfServers(2);

jobclass{1} = ClosedClass(model, 'Class1', 2, node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', 0, node{1}, 0);
jobclass{3} = ClosedClass(model, 'Class3', 1, node{1}, 0);

node{1}.setService(jobclass{1}, Erlang(3,2));
node{1}.setService(jobclass{2}, HyperExp(0.5,3.0,10.0));
node{1}.setService(jobclass{3}, Exp(1));

node{2}.setService(jobclass{1}, HyperExp(0.1,1.0,10.0));
node{2}.setService(jobclass{2}, Exp(2));
node{2}.setService(jobclass{3}, Exp(3));

M = model.getNumberOfStations();
K = model.getNumberOfClasses();

P = model.initRoutingMatrix;
P{1,1} = [0.3,0.1; 0.2,0];
P{1,2} = [0.6,0; 0.8,0];
P{1,3} = zeros(M);

P{2,1} = [0,0; 1,0];
P{2,2} = [0,1; 0,0];
P{2,3} = zeros(M);

P{3,1} = zeros(M);
P{3,2} = zeros(M);
P{3,3} = circul(M);

model.link(P);
%%
% This part illustrates the execution of different solvers
solver{1} = CTMC(model);
solver{end+1} = JMT(model,'seed',23000,'samples',5e3,'verbose',true);
solver{end+1} = SSA(model,'seed',23000,'samples',5e3,'verbose',true);
solver{end+1} = FLD(model);
solver{end+1} = MVA(model);
solver{end+1} = NC(model,'exact');
solver{end+1} = MAM(model);
solver{end+1} = LINE(model,'seed',23000);
solver{end+1} = DES(model,'seed',23000,'samples',5e3,'verbose',true);

for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());
    AvgTable{s} = solver{s}.getAvgTable();
    AvgChainTable{s} = solver{s}.getAvgChainTable();
    AvgSysChainTable{s} = solver{s}.getAvgSysTable();
    AvgTable{s}
    AvgChainTable{s}
    AvgSysChainTable{s}
end