clear node jobclass solver; 
model = Network('model');

node{1} = ClassSwitch(model,'CS',[0,1;1,0]);
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
node{3} = Queue(model, 'Queue2', SchedStrategy.PS);
node{4} = Queue(model, 'Delay',  SchedStrategy.INF);

jobclass{1} = ClosedClass(model, 'Class1', 15, node{4}, 0);
jobclass{2} = ClosedClass(model, 'Class2',  5, node{4}, 0);

node{2}.setService(jobclass{1}, Exp.fitMean(1.5)); % mean = 1.5
node{2}.setService(jobclass{2}, Erlang.fitMeanAndOrder(1.5,2)); % mean = 1.5

node{3}.setService(jobclass{1}, Erlang.fitMeanAndOrder(1.5,2)); % mean = 1.5
node{3}.setService(jobclass{2}, Exp.fitMean(1.5)); % mean = 1.5

node{4}.setService(jobclass{1}, Exp.fitMean(1.0)); % mean = 1
node{4}.setService(jobclass{2}, Exp.fitMean(1.0)); % mean = 1

model.addLink(node{2}, node{1});
model.addLink(node{3}, node{1});
model.addLink(node{1}, node{4});
model.addLink(node{4}, node{2});
model.addLink(node{4}, node{3});

node{1}.setRouting(jobclass{1},RoutingStrategy.RAND);
node{2}.setRouting(jobclass{1},RoutingStrategy.RAND);
node{3}.setRouting(jobclass{1},RoutingStrategy.RAND);
node{4}.setRouting(jobclass{1},RoutingStrategy.RROBIN);

node{1}.setRouting(jobclass{2},RoutingStrategy.RAND);
node{2}.setRouting(jobclass{2},RoutingStrategy.RAND);
node{3}.setRouting(jobclass{2},RoutingStrategy.RAND);
node{4}.setRouting(jobclass{2},RoutingStrategy.WRROBIN,node{2},1);
node{4}.setRouting(jobclass{2},RoutingStrategy.WRROBIN,node{3},2);

simoptions = Solver.defaultOptions; 
simoptions.verbose = true;
simoptions.seed = 23000;
solver = {};
solver{end+1} = JMT(model, simoptions);
solver{end+1} = DES(model, simoptions);
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());    
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
