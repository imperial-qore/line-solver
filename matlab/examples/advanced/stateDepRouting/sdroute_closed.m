clear node jobclass solver AvgTable

model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
node{3} = Queue(model, 'Queue2', SchedStrategy.PS);

jobclass{1} = ClosedClass(model, 'Class1', 1, node{1}, 0);

node{1}.setService(jobclass{1}, HyperExp.fitMeanAndSCV(1,25));
node{2}.setService(jobclass{1}, Exp(1));
node{3}.setService(jobclass{1}, Exp(2));

model.addLink(node{1}, node{1});
model.addLink(node{1}, node{2});
model.addLink(node{1}, node{3});
model.addLink(node{2}, node{1});
model.addLink(node{3}, node{1});

node{1}.setRouting(jobclass{1},RoutingStrategy.RROBIN);
node{2}.setProbRouting(jobclass{1}, node{1}, 1.0);
node{3}.setProbRouting(jobclass{1}, node{1}, 1.0);

solver={};
solver{end+1} = CTMC(model,'keep',true);
solver{end+1} = JMT(model,'samples',1e5,'seed',23000);
solver{end+1} = SSA(model,'verbose',true,'samples',1e4,'seed',23000);

AvgTable = {};
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',strrep(solver{s}.getName(),'Solver',''));
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
