clear node jobclass AvgTable

%%clearvars -except exampleName; 
model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
node{3} = Queue(model, 'Queue2', SchedStrategy.PS);

jobclass{1} = ClosedClass(model, 'Class1', 1, node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', 1, node{1}, 0);

node{1}.setService(jobclass{1}, Erlang(3,2));
node{1}.setService(jobclass{2}, HyperExp(0.5,3.0,10.0));

node{2}.setService(jobclass{1}, HyperExp(0.1,1.0,10.0));
node{2}.setService(jobclass{2}, MMPP2(1,2,3,4));

node{3}.setService(jobclass{1}, HyperExp(0.1,1.0,10.0));
node{3}.setService(jobclass{2}, Erlang(1,2));

model.addLink(node{1}, node{1});
model.addLink(node{1}, node{2});
model.addLink(node{1}, node{3});
model.addLink(node{2}, node{1});
model.addLink(node{3}, node{1});

node{1}.setProbRouting(jobclass{1}, node{1}, 0.0)
node{1}.setProbRouting(jobclass{1}, node{2}, 0.3)
node{1}.setProbRouting(jobclass{1}, node{3}, 0.7)
node{2}.setProbRouting(jobclass{1}, node{1}, 1.0)
node{3}.setProbRouting(jobclass{1}, node{1}, 1.0)

node{1}.setRouting(jobclass{2},RoutingStrategy.RAND);
node{2}.setRouting(jobclass{2},RoutingStrategy.RAND);
node{3}.setRouting(jobclass{2},RoutingStrategy.RAND);

solver={};
options = Solver.defaultOptions;
options.seed = 23000;
%solver{end+1} = CTMC(model,options);
solver{end+1} = JMT(model,options);
solver{end+1} = DES(model,options);
%solver{end+1} = SSA(model,options);
%solver{end+1} = FLD(model,options);
%solver{end+1} = MVA(model,options);
%solver{end+1} = NC(model,options);
%solver{end+1} = MAM(model,options);

for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
