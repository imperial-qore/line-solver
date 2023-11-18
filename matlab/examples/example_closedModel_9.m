clear node jobclass;

model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.FCFS);
node{3} = Queue(model, 'Queue2', SchedStrategy.FCFS);
jobclass{1} = ClosedClass(model, 'Class1', 10, node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', 10, node{1}, 0);

node{1}.setService(jobclass{1}, Exp.fitMean(1.0)); % mean = 1
node{2}.setService(jobclass{1}, Exp.fitMean(1.5)); % mean = 1.5
node{3}.setService(jobclass{1}, Exp.fitMean(3.0)); % mean = 2.0

node{1}.setService(jobclass{2}, Exp.fitMean(1.0)); % mean = 1
node{2}.setService(jobclass{2}, Exp.fitMean(1.5)); % mean = 1.5
node{3}.setService(jobclass{2}, Exp.fitMean(3.0)); % mean = 2.0

P = model.initRoutingMatrix;
P{1} = Network.serialRouting(node);
P{2} = Network.serialRouting(node);
%P{jobclass{1},jobclass{1}} = [0.7,0.3,0;1.0,0,0;1.0,0,0];
%P{jobclass{2},jobclass{2}} = [0.7,0,0.3;1.0,0,0;1.0,0,0];

% This may be alternatively specified as:
%P{1}(node{1},[node{1},node{2}]) = [0.7,0.3];
%P{1}(node{2},[node{1},node{2}]) = [1.0,0];

model.link(P);

% solver{1} = SolverCTMC(model,'keep',true);
% solver{end+1} = SolverJMT(model,'seed',23000,'verbose',true,'keep',true);
% solver{end+1} = SolverSSA(model,'seed',23000,'verbose',true,'samples',5e3);
% solver{end+1} = SolverFluid(model);
options = SolverMVA.defaultOptions;
%options.method = 'lin';
solver{1} = SolverMVA(model, options);
% solver{end+1} = SolverNC(model,'exact');
% solver{end+1} = SolverMAM(model);
% solver{end+1} = LINE(model);

AvgTable = cell(1,length(solver));
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
