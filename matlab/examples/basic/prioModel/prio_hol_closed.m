clear node jobclass solver AvgTable

%% multiclass example with PS, SIRO, FCFS, HOL priority
model = Network('MyNetwork');

%% Block 1: nodes
node{1} = Delay(model, 'SlowDelay');
node{2} = Queue(model, 'FCFSQueue', SchedStrategy.FCFS);
node{3} = Queue(model, 'SIROQueue', SchedStrategy.SIRO);
node{4} = Queue(model, 'PSQueue', SchedStrategy.PS);
node{5} = Queue(model, 'HOLQueue', SchedStrategy.HOL);
node{6} = Delay(model, 'FastDelay');

%% Block 2: classes
jobclass{1} = ClosedClass(model, 'Class1', 18, node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', 18, node{1}, 1);
jobclass{3} = ClosedClass(model, 'Class3', 18, node{1}, 0);

node{1}.setService(jobclass{1}, Exp.fitMean(10.000000)); 
node{1}.setService(jobclass{2}, Exp.fitMean(10.000000)); 
node{1}.setService(jobclass{3}, Exp.fitMean(10.000000)); 

node{2}.setService(jobclass{1}, Exp.fitMean(0.300000)); 
node{2}.setService(jobclass{2}, Exp.fitMean(0.500000)); 
node{2}.setService(jobclass{3}, Exp.fitMean(0.600000)); 

node{3}.setService(jobclass{1}, Exp.fitMean(1.100000)); 
node{3}.setService(jobclass{2}, Exp.fitMean(1.300000)); 
node{3}.setService(jobclass{3}, Exp.fitMean(1.500000)); 

node{4}.setService(jobclass{1}, Exp.fitMean(1.000000)); 
node{4}.setService(jobclass{2}, Exp.fitMean(1.100000)); 
node{4}.setService(jobclass{3}, Exp.fitMean(1.900000)); 

node{5}.setService(jobclass{1}, Exp.fitMean(2.500000)); 
node{5}.setService(jobclass{2}, Exp.fitMean(1.900000)); 
node{5}.setService(jobclass{3}, Exp.fitMean(4.300000)); 

node{6}.setService(jobclass{1}, Exp.fitMean(1.000000)); 
node{6}.setService(jobclass{2}, Exp.fitMean(1.000000)); 
node{6}.setService(jobclass{3}, Exp.fitMean(1.000000)); 

%% Block 3: topology
P = model.initRoutingMatrix(); % initialize routing matrix 
P{1,1}(1,2) = 1; % (Source of Customers,Class0) -> (WebServer,Class0)
P{1,1}(2,3) = 2.500000e-01; % (WebServer,Class0) -> (Storage1,Class0)
P{1,1}(2,4) = 2.500000e-01; % (WebServer,Class0) -> (Storage2,Class0)
P{1,1}(2,5) = 2.500000e-01; % (WebServer,Class0) -> (Storage3,Class0)
P{1,1}(2,6) = 2.500000e-01; % (WebServer,Class0) -> (Out,Class0)
P{1,1}(3,2) = 1; % (Storage1,Class0) -> (WebServer,Class0)
P{1,1}(4,2) = 1; % (Storage2,Class0) -> (WebServer,Class0)
P{1,1}(5,2) = 1; % (Storage3,Class0) -> (WebServer,Class0)
P{2,2}(1,2) = 1; % (Source of Customers,Class1) -> (WebServer,Class1)
P{2,2}(2,3) = 2.500000e-01; % (WebServer,Class1) -> (Storage1,Class1)
P{2,2}(2,4) = 2.500000e-01; % (WebServer,Class1) -> (Storage2,Class1)
P{2,2}(2,5) = 2.500000e-01; % (WebServer,Class1) -> (Storage3,Class1)
P{2,2}(2,6) = 2.500000e-01; % (WebServer,Class1) -> (Out,Class1)
P{2,2}(3,2) = 1; % (Storage1,Class1) -> (WebServer,Class1)
P{2,2}(4,2) = 1; % (Storage2,Class1) -> (WebServer,Class1)
P{2,2}(5,2) = 1; % (Storage3,Class1) -> (WebServer,Class1)
P{3,3}(1,2) = 1; % (Source of Customers,Class2) -> (WebServer,Class2)
P{3,3}(2,3) = 2.500000e-01; % (WebServer,Class2) -> (Storage1,Class2)
P{3,3}(2,4) = 2.500000e-01; % (WebServer,Class2) -> (Storage2,Class2)
P{3,3}(2,5) = 2.500000e-01; % (WebServer,Class2) -> (Storage3,Class2)
P{3,3}(2,6) = 2.500000e-01; % (WebServer,Class2) -> (Out,Class2)
P{3,3}(3,2) = 1; % (Storage1,Class2) -> (WebServer,Class2)
P{3,3}(4,2) = 1; % (Storage2,Class2) -> (WebServer,Class2)
P{3,3}(5,2) = 1; % (Storage3,Class2) -> (WebServer,Class2)
P{1,1}(6,1) = 1;
P{2,2}(6,1) = 1;
P{3,3}(6,1) = 1;
model.link(P);

%%
options = lineDefaults;
options.seed = 23000;
options.cutoff=1;
options.samples=1e4;
solver={};
%solver{end+1} = CTMC(model,options); % CTMC is infinite on this model
%solver{end+1} = FLD(model,options);
solver{end+1} = MVA(model,options);
%solver{end+1} = MAM(model,options);
solver{end+1} = JMT(model,options);
solver{end+1} = SSA(model,options);
%solver{end+1} = NC(model,options);
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',strrep(solver{s}.getName(),'Solver',''));
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
