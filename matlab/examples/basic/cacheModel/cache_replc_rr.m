clear solver AvgTable;

model = Network('model');

n = 5; % number of items
m = 2; % cache capacity 

source = Source(model, 'Source');
cacheNode = Cache(model, 'Cache', n, m, ReplacementStrategy.RR);
sink = Sink(model, 'Sink');

jobClass = OpenClass(model, 'InitClass', 0);
hitClass = OpenClass(model, 'HitClass', 0);
missClass = OpenClass(model, 'MissClass', 0);

source.setArrival(jobClass, Exp(2));

%pAccess = DiscreteSampler((1/n)*ones(1,n));  % uniform item references
pAccess = Zipf(1.4, n);  % Zipf-like item references
cacheNode.setRead(jobClass, pAccess);


cacheNode.setHitClass(jobClass, hitClass);
cacheNode.setMissClass(jobClass, missClass);

P = model.initRoutingMatrix;

P{jobClass, jobClass}(source, cacheNode) =  1.0;
P{hitClass, hitClass}(cacheNode, sink) =  1.0;
P{missClass, missClass}(cacheNode, sink) =  1.0;

model.link(P);
%%
solver{1} = CTMC(model,'keep',false,'cutoff',1);
AvgTable{1} = solver{1}.getAvgNodeTable; AvgTable{1}

model.reset;
solver{2} = SSA(model,'samples',1e4,'verbose',true,'method','serial','seed',23000);
AvgTable{2} = solver{2}.getAvgNodeTable; AvgTable{2}

model.reset;
solver{3} = MVA(model);
AvgTable{3} = solver{3}.getAvgNodeTable; AvgTable{3}

model.reset;
solver{4} = NC(model);
AvgTable{4} = solver{4}.getAvgNodeTable; AvgTable{4}

hitRatio=cacheNode.getHitRatio
missRatio=cacheNode.getMissRatio