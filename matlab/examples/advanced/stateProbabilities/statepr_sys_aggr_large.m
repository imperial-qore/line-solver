clear node jobclass 

model = Network('model');

node{1} = Queue(model, 'Queue1', SchedStrategy.PS);
node{2} = Queue(model, 'Queue2', SchedStrategy.PS);
node{3} = Queue(model, 'Queue3', SchedStrategy.PS);

node{3}.setNumServers(3);

N=[1,1,1,1];
jobclass{1} = ClosedClass(model, 'Class1', N(1), node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', N(2), node{1}, 0);
jobclass{3} = ClosedClass(model, 'Class3', N(3), node{1}, 0);
jobclass{4} = ClosedClass(model, 'Class4', N(4), node{1}, 0);

node{1}.setService(jobclass{1}, Exp(1));
node{1}.setService(jobclass{2}, Exp(2));
node{1}.setService(jobclass{3}, Exp(1));
node{1}.setService(jobclass{4}, Exp(1));

node{2}.setService(jobclass{1}, Exp(3));
node{2}.setService(jobclass{2}, Exp(4));
node{2}.setService(jobclass{3}, Exp(5));
node{2}.setService(jobclass{4}, Exp(1));

node{3}.setService(jobclass{1}, Exp(1));
node{3}.setService(jobclass{2}, Exp(3));
node{3}.setService(jobclass{3}, Exp(5));
node{3}.setService(jobclass{4}, Exp(2));

K = length(jobclass);
P = cell(K,K);

P{1,1} = [0,1,0; 0,0,1; 0,0,0];
P{1,2} = [0,0,0; 0,0,0; 1,0,0];
P{1,3} = [0,0,0; 0,0,0; 0,0,0];
P{1,4} = [0,0,0; 0,0,0; 0,0,0];

P{2,1} = [0,0,0; 0,0,0; 1,0,0];
P{2,2} = [0,1,0; 0,0,1; 0,0,0];
P{2,3} = [0,0,0; 0,0,0; 0,0,0];
P{2,4} = [0,0,0; 0,0,0; 0,0,0];

P{3,1} = [0,0,0; 0,0,0; 0,0,0];
P{3,2} = [0,0,0; 0,0,0; 0,0,0];
P{3,3} = [0,1,0; 0,0,1; 0,0,0];
P{3,4} = [0,0,0; 0,0,0; 1,0,0];

P{4,1} = [0,0,0; 0,0,0; 0,0,0];
P{4,2} = [0,0,0; 0,0,0; 0,0,0];
P{4,3} = [0,0,0; 0,0,0; 1,0,0];
P{4,4} = [0,0,1; 0,0,0; 0,0,0];

%%pause
model.link(P);

M = model.getNumberOfStations();
K = model.getNumberOfClasses();

solver={};
options = Solver.defaultOptions;
options.verbose=1;
%options.samples=1e5;

n=[0,0,0,0;
   0,0,0,0;
   N(1),N(2),N(3),N(4)];
for i=1:M
    node{i}.setState(n(i,:));
end
state = model.getState;

tic;
solver = LINE.load('ctmc',model);
Pr_ctmc = solver.getProbSysAggr()

toc
tic;
solver = JMT(model,'seed',23000,'samples',1e5);
Pr_jmt = solver.getProbSysAggr()

toc
%%
tic;
solver = NC(model,'seed',23000,'method','ls','samples',10000);
Pr_nc = solver.getProbSysAggr()

toc
