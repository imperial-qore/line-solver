clear node jobclass solver AvgTable

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
P = cell(K,K);

P{1} = [0,1,0; 0,0,1; 1,0,0];
P{2} = [0,1,0; 0,0,1; 1,0,0];

model.link(P);

M = model.getNumberOfStations();
K = model.getNumberOfClasses();

% This part illustrates the execution of different solvers
fprintf(1,'This example illustrates the calculation of probabilities via normalizing constants.\n')

options = Solver.defaultOptions;
options.verbose=1;
options.seed=23000;
%options.samples=1e5;

n=[-1,-1;
   -1,-1;
    0, 0];  % rows set to -1 are ignored in the calculation of the marginal probabilities

for i=1:M
    node{i}.setState(n(i,:));
end
state = model.getState;

%%
i=M;
solver = CTMC(model,options);
Pr = solver.getProbAggr(node{i});
fprintf(1,'Station %d is in state %s with probability %d\n',i,mat2str(state{i}),Pr);
Pmarg_ctmc = Pr

solver = NC(model,options);
Pr = solver.getProbAggr(node{i});
fprintf(1,'Station %d is in state %s with probability %d\n',i,mat2str(state{i}),Pr);
Pmarg_nc = Pr
