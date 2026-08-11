clear node jobclass solver;
model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.FCFS);
node{2}.setNumServers(3);

% Default: scheduling is set as FCFS everywhere, routing as Random
jobclass{1} = ClosedClass(model, 'Class1', 4, node{1}, 0);
jobclass{2} = ClosedClass(model, 'Class2', 2, node{1}, 0);

node{1}.setService(jobclass{1}, Exp(1));
node{1}.setService(jobclass{2}, Exp(1));

node{2}.setService(jobclass{1}, Exp(1));
node{2}.setService(jobclass{2}, Exp(10));

myP = model.initRoutingMatrix;
myP{1} = Network.serialRouting(node);
myP{2} = Network.serialRouting(node);
model.link(myP);
%
options = QNS.defaultOptions;
%options.verbose = false;

solver{1} = CTMC(model);
solver{end+1} =QNS(model,'conway');
solver{end+1} =QNS(model,'reiser');
solver{end+1} =QNS(model,'rolia');
solver{end+1} =QNS(model,'zhou');
options = MVA.defaultOptions;
options.config.multiserver = 'softmin';
solver{end+1} = MVA(model, options);
options.config.multiserver = 'seidmann';
solver{end+1} = MVA(model, options);
solver{end+1} = NC(model)

AvgTable = cell(1,length(solver));
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',strrep(solver{s}.getName(),'Solver',''));    
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end

