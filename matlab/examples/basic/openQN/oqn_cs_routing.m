clear node jobclass solver AvgTable

%model = JMT2LINE('example_openModel_3.jsimg');
model = Network('myModel');

%% Block 1: nodes
node{1} = Source(model, 'Source 1');
node{2} = Queue(model, 'Queue 1', SchedStrategy.PS);
node{3} = ClassSwitch(model, 'ClassSwitch 1'); % Class switching is embedded in the routing matrix P 
node{4} = Sink(model, 'Sink 1');
node{5} = Queue(model, 'Queue 2', SchedStrategy.PS);

%% Block 2: classes
jobclass{1} = OpenClass(model, 'Class A', 0);
jobclass{2} = OpenClass(model, 'Class B', 0);
jobclass{3} = OpenClass(model, 'Class C', 0);

node{1}.setArrival(jobclass{1}, Exp.fitMean(0.500000)); % (Source 1,Class A)
node{1}.setArrival(jobclass{2}, Exp.fitMean(1.000000)); % (Source 1,Class B)
node{1}.setArrival(jobclass{3}, Disabled.getInstance()); % (Source 1,Class C)
node{2}.setService(jobclass{1}, Exp.fitMean(0.200000)); % (Queue 1,Class A)
node{2}.setService(jobclass{2}, Exp.fitMean(0.300000)); % (Queue 1,Class B)
node{2}.setService(jobclass{3}, Exp.fitMean(0.333333)); % (Queue 1,Class C)
node{5}.setService(jobclass{1}, Exp.fitMean(1.000000)); % (Queue 2,Class A)
node{5}.setService(jobclass{2}, Exp.fitMean(1.000000)); % (Queue 2,Class B)
node{5}.setService(jobclass{3}, Exp.fitMean(0.150000)); % (Queue 2,Class C)

%% Block 3: topology
C = node{3}.initClassSwitchMatrix();
C = eye(length(C));
node{3}.setClassSwitchingMatrix(C);

P = model.initRoutingMatrix(); % initialize routing matrix 
P{1,1}(1,2) = 1; % (Source 1,Class A) -> (Queue 1,Class A)
P{1,1}(2,3) = 1; % (Queue 1,Class A) -> (ClassSwitch 1,Class A)
P{1,1}(5,4) = 1; % (Queue 2,Class A) -> (Sink 1,Class A)
P{1,3}(3,5) = 1; % (ClassSwitch 1,Class A) -> (Queue 2,Class C)
P{2,2}(1,2) = 1; % (Source 1,Class B) -> (Queue 1,Class B)
P{2,2}(2,3) = 1; % (Queue 1,Class B) -> (ClassSwitch 1,Class B)
P{2,2}(5,4) = 1; % (Queue 2,Class B) -> (Sink 1,Class B)
P{2,3}(3,5) = 1; % (ClassSwitch 1,Class B) -> (Queue 2,Class C)
P{3,3}(1,2) = 1; % (Source 1,Class C) -> (Queue 1,Class C)
P{3,3}(2,3) = 1; % (Queue 1,Class C) -> (ClassSwitch 1,Class C)
P{3,3}(3,5) = 1; % (ClassSwitch 1,Class C) -> (Queue 2,Class C)
P{3,3}(5,4) = 1; % (Queue 2,Class C) -> (Sink 1,Class C)
model.link(P);

options = Solver.defaultOptions;
options.keep = true;
options.verbose = 1;
options.cutoff = [1,1,0;3,3,0;0,0,3]; % works well with 7
options.seed = 23000;
options.samples = 1e5;

%disp('This example shows the execution of the solver on a 1-class 2-node open model.')
% This part illustrates the execution of different solvers
solver={};
solver{end+1} = CTMC(model,options); % CTMC is infinite on this model
solver{end+1} = FLD(model,options);
solver{end+1} = MVA(model,options);
solver{end+1} = MAM(model,options);
solver{end+1} = NC(model,options);
solver{end+1} = JMT(model,options);
solver{end+1} = SSA(model,options);
solver{end+1} = DES(model,options);
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',strrep(solver{s}.getName(),'Solver',''));
    AvgTable{s} = solver{s}.getAvgTable()
    AvgTable{s}
end
