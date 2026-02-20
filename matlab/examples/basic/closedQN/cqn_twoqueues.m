clear node solver jobclass;
% reducible routing matrix
model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.FCFS);
node{3} = Queue(model, 'Queue2', SchedStrategy.FCFS);

jobclass{1} = ClosedClass(model, 'Class1', 1, node{1}, 0);

node{1}.setService(jobclass{1}, Exp.fitMean(1.0)); % mean = 1
node{2}.setService(jobclass{1}, Exp.fitMean(1.5)); % mean = 1.5
node{3}.setService(jobclass{1}, Exp.fitMean(3.0)); % mean = 2.0

P = model.initRoutingMatrix;
P{1}(1,1) = 0.2;
P{1}(1,2) = 0.3;
P{1}(1,3) = 0.5;
P{1}(2,2) = 1.0;
P{1}(3,3) = 1.0;
model.link(P);

options = MVA.defaultOptions;
solver{1} = MVA(model, options);

AvgTable = cell(1,length(solver));
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',strrep(solver{s}.getName(),'Solver',''));
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
