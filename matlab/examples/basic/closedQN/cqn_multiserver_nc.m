clear solver node jobclass;

% Closed Delay + multiserver FCFS queue.
% Exercises the exact multiserver path of NC: with method 'exact'
% the c-server station is converted to a load-dependent station with rate
% min(n,c) (runAnalyzer.m), solved via comomld/pfqn_comomrm_ld. CTMC is the
% exact ground truth; MVA provides an additional cross-check.

model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.FCFS);
node{2}.setNumServers(3);

jobclass{1} = ClosedClass(model, 'Class1', 5, node{1}, 0);

node{1}.setService(jobclass{1}, Exp.fitMean(1.0)); % mean = 1
node{2}.setService(jobclass{1}, Exp.fitMean(0.8)); % mean = 0.8

P = model.initRoutingMatrix;
P{1} = [0,1.0;1.0,0];
model.link(P);

solver{1} = CTMC(model,'exact');
solver{end+1} = MVA(model);
solver{end+1} = NC(model,'exact');

AvgTable = cell(1,length(solver));
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',strrep(solver{s}.getName(),'Solver',''));
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
