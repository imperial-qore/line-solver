clear node jobclass solver AvgTable

model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
node{3} = Queue(model, 'Queue2', SchedStrategy.FCFS);

jobclass{1} = ClosedClass(model, 'Class1', 1, node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', 2, node{1}, 0);

% renewal
map21 = APH([1,0],[-2,2; 0,-0.5]);
map22 = PH(1,-1); 
% non-renewal
map31 = MAP([-20,0; 0,-1],[0,20;0.8,0.2]);
map32 = MAP([-4, 3; 4, -6],[1,0;0,2]);
%map32 = MMPP2(1,2,3,4).toMAP();

node{1}.setService(jobclass{1}, HyperExp.fitMeanAndSCV(1,25));
node{2}.setService(jobclass{1}, map21);
node{3}.setService(jobclass{1}, map31);

node{1}.setService(jobclass{2}, HyperExp.fitMeanAndSCV(1,25));
node{2}.setService(jobclass{2}, map22);
node{3}.setService(jobclass{2}, map32);

model.addLink(node{1}, node{1});
model.addLink(node{1}, node{2});
model.addLink(node{1}, node{3});
model.addLink(node{2}, node{1});
model.addLink(node{3}, node{1});

node{1}.setRouting(jobclass{1},RoutingStrategy.RROBIN);
node{2}.setProbRouting(jobclass{1}, node{1}, 1.0)
node{3}.setProbRouting(jobclass{1}, node{1}, 1.0)

node{1}.setRouting(jobclass{2},RoutingStrategy.RROBIN);
node{2}.setProbRouting(jobclass{2}, node{1}, 1.0)
node{3}.setProbRouting(jobclass{2}, node{1}, 1.0)

solver={};
solver{end+1} = CTMC(model,'keep',false);
solver{end+1} = JMT(model,'samples',1e5,'seed',23000);
solver{end+1} = SSA(model,'verbose',true,'samples',5e3,'seed',23000);

AvgTable = {};
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
