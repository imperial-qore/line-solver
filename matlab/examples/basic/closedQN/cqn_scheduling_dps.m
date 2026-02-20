clear node jobclass1

model = Network('model');

node1 = Delay(model, 'Delay');
node2 = Queue(model, 'Queue1', SchedStrategy.PS);
node3 = Queue(model, 'Queue2', SchedStrategy.DPS);

jobclass1 = ClosedClass(model, 'Class1', 2, node1, 0);
jobclass2 = ClosedClass(model, 'Class2', 1, node1, 0);

node1.setService(jobclass1, Exp(3));
node1.setService(jobclass2, Exp(0.5));

% this weight assigment is ignored since the node is PS
w1 = 5; node2.setService(jobclass1, Exp(0.1), w1);
w2 = 1; node2.setService(jobclass2, Exp(1), w2);

% this is not ignored since the node is DPS
w1 = 1; node3.setService(jobclass1, Exp(0.1), w1);
w2 = 5; node3.setService(jobclass2, Exp(1), w2);

P = model.initRoutingMatrix();

P.set(jobclass1, jobclass1, node1, node2, 0.300000); % (Delay,Class1) -> (Queue1,Class1)
P.set(jobclass1, jobclass1, node1, node3, 0.700000); % (Delay,Class1) -> (Queue2,Class1)
P.set(jobclass1, jobclass1, node2, node1, 1.000000); % (Queue1,Class1) -> (Delay,Class1)
P.set(jobclass1, jobclass1, node3, node1, 1.000000); % (Queue2,Class1) -> (Delay,Class1)

P.set(jobclass2, jobclass2, node1, node2, 0.700000); % (Delay,Class2) -> (Queue1,Class2)
P.set(jobclass2, jobclass2, node1, node3, 0.300000); % (Delay,Class2) -> (Queue2,Class2)
P.set(jobclass2, jobclass2, node2, node1, 1.000000); % (Queue1,Class2) -> (Delay,Class2)
P.set(jobclass2, jobclass2, node3, node1, 1.000000); % (Queue2,Class2) -> (Delay,Class2)

model.link(P);

% This part illustrates the execution of different solvers

solver={};
options = Solver.defaultOptions;
options.verbose=1;
options.samples=1e4;
options.seed = 23000;
solver{end+1} = CTMC(model,options);
solver{end+1} = JMT(model,options);
solver{end+1} = FLD(model,options);
solver{end+1} = MVA(model,options);
solver{end+1} = DES(model,options);
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',strrep(solver{s}.getName(),'Solver',''));
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
