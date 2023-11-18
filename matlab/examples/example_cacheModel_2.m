clear solver AvgTable;

model = Network('model');

n = 5; % number of items
m = 2; % cache capacity

delay = Delay(model, 'Delay');
cacheNode = Cache(model, 'Cache', n, m, ReplacementStrategy.FIFO);

jobClass = ClosedClass(model, 'JobClass', 1, delay, 0);
hitClass = ClosedClass(model, 'HitClass', 0, delay, 0);
missClass = ClosedClass(model, 'MissClass', 0, delay, 0);

delay.setService(jobClass, Exp(1));

pAccess = DiscreteSampler((1/n)*ones(1,n));  % uniform item references
cacheNode.setRead(jobClass, pAccess);

cacheNode.setHitClass(jobClass, hitClass);
cacheNode.setMissClass(jobClass, missClass);

P = model.initRoutingMatrix;

P{jobClass, jobClass}(delay, cacheNode) =  1.0;
P{hitClass, jobClass}(cacheNode, delay) =  1.0;
P{missClass, jobClass}(cacheNode, delay) =  1.0;

model.link(P);

solver{1} = SolverCTMC(model,'keep',false,'seed',1);
AvgTable{1} = solver{1}.getAvgNodeTable; AvgTable{1}

model.reset;
solver{2} = SolverSSA(model,'samples',1e4,'verbose',true,'method','serial','seed',1);
AvgTable{2} = solver{2}.getAvgNodeTable; AvgTable{2}

solver{3} = SolverMVA(model,'seed',1);
AvgTable{3} = solver{3}.getAvgNodeTable; AvgTable{3}
