clear node jobclass 

model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
node{3} = Queue(model, 'Queue2', SchedStrategy.PS);
node{3}.setNumServers(2);

N=[2,0];
jobclass{1} = ClosedClass(model, 'Class1', N(1), node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', N(2), node{1}, 0);

node{1}.setService(jobclass{1}, Exp(1));
node{1}.setService(jobclass{2}, Exp(1));

node{2}.setService(jobclass{1}, Exp(3));
node{2}.setService(jobclass{2}, Exp(4));

node{3}.setService(jobclass{1}, Exp(1));
node{3}.setService(jobclass{2}, Exp(3));

K = length(jobclass);
P = model.initRoutingMatrix;

P{1,1} = [0,1,0; 0,0,0; 1,0,0];
P{1,2} = [0,0,0; 0,0,1; 0,0,0];
P{2,1} = [0,1,0; 0,0,0; 1,0,0];
P{2,2} = [0,0,0; 0,0,1; 0,0,0];

model.link(P);

M = model.getNumberOfStations();
K = model.getNumberOfClasses();

%This part illustrates the execution of different solvers
fprintf(1,'This example illustrates the calculation of probabilities via normalizing constants.\n')

options = Solver.defaultOptions;
options.verbose=1;
options.samples=2e4;
options.seed = 23000;
% set a custom initial state
n=[ 0,0;
    1,0;
    0,1];

for i=1:M
    node{i}.setState(n(i,:));
end
state = model.getState;

%% getProbAggr
% marginal probabilities for the aggregated state space where station i
% is specified by a tuple (nir), r=1,...R, R being the number of classes.
solver = CTMC(model,options);
Pr = solver.getProbAggr(node{M});
fprintf(1,'Station %d is in state %s with probability %d\n',i,mat2str(state{i}),Pr);
Pmarga_ctmc = Pr

solver = NC(model,options);
Pr = solver.getProbAggr(node{M});
fprintf(1,'Station %d is in state %s with probability %d\n',i,mat2str(state{i}),Pr);
Pmarga_nc = Pr

solver = SSA(model,options);
Pr = solver.getProbAggr(node{M});
fprintf(1,'Station %d is in state %s with probability %d\n',i,mat2str(state{i}),Pr);
Pmarga_ssa = Pr

solver = JMT(model,options);
Pr = solver.getProbAggr(node{M});
fprintf(1,'Station %d is in state %s with probability %d\n',i,mat2str(state{i}),Pr);
Pmarga_jmt = Pr

%% getProb
% marginal probabilities for the detailed state space, which tracks also
% the phases of service
solver = CTMC(model,options);
Pr = solver.getProb(node{M});
fprintf(1,'Station %d is in state %s with probability %d\n',i,mat2str(state{i}),Pr);
Pmarg_ctmc = Pr

% solver = NC(model,options);
% Pr = solver.getProb(node{M});
% Pmarg_nc = Pr

solver = SSA(model,options);
Pr = solver.getProb(node{M});
fprintf(1,'Station %d is in state %s with probability %d\n',i,mat2str(state{i}),Pr);
Pmarg_ssa = Pr

%% getProbSysAggr
% joint state probabilities for the aggregated state space where station i
% is specified by a tuple (nir), r=1,...R, R being the number of classes.
solver = CTMC(model,options);
Pr = solver.getProbSysAggr();
Pjointa_ctmc = Pr

solver = NC(model,options);
Pr = solver.getProbSysAggr();
Pjointa_nc = Pr

solver = SSA(model,options);
Pr = solver.getProbSysAggr();
Pjointa_ssa = Pr

solver = JMT(model,options);
Pr = solver.getProbSysAggr();
Pjointa_jmt = Pr

%% getProbSys
% joint state probabilities for the detailed state space, which tracks also
% the phases of service
solver = CTMC(model,options);
Pr = solver.getProbSys();
Pjoint_ctmc = Pr

% solver = NC(model,options);
% Pr = solver.getProbSys();
% Pjoint_nc = Pr

solver = SSA(model,options);
Pr = solver.getProbSys();
Pjoint_ssa = Pr

% options.samples = 3e5;
% solver = JMT(model,options);
% Pr = solver.getProbSys();
% Pjoint_jmt = Pr

