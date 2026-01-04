clear solver node jobclass;

model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.FCFS);
jobclass{1} = ClosedClass(model, 'Class1', 10, node{1}, 0);

node{1}.setService(jobclass{1}, Exp.fitMean(1.0)); % mean = 1
node{2}.setService(jobclass{1}, Exp.fitMean(1.5)); % mean = 1.5

P = model.initRoutingMatrix;
P{1} = [0.7,0.3;1.0,0];

% This may be alternatively specified as:
%P{1}(node{1},[node{1},node{2}]) = [0.7,0.3];
%P{1}(node{2},[node{1},node{2}]) = [1.0,0];

model.link(P);

solver{1} = CTMC(model,'keep',true);
solver{end+1} = JMT(model,'seed',23000,'verbose',true,'keep',true);
solver{end+1} = SSA(model,'seed',23000,'verbose',true,'samples',5e3);
solver{end+1} = FLD(model);
solver{end+1} = MVA(model);
solver{end+1} = NC(model,'exact');
solver{end+1} = MAM(model);
solver{end+1} = LINE(model);
solver{end+1} = DES(model,'seed',23000,'verbose',true);

AvgTable = cell(1,length(solver));
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
