clear solver AvgTable;

% CLIMB (transposition) replacement: on a hit an item moves up one position;
% on a miss it enters at the tail. Exact in CTMC; simulated in SSA/LDES.
% Not product-form, so MVA/NC/FLD reject it (use CTMC/SSA/LDES).

model = Network('model');

n = 5; % number of items
m = 2; % cache capacity

delay = Delay(model, 'Delay');
cacheNode = Cache(model, 'Cache', n, m, ReplacementStrategy.CLIMB);

jobClass = ClosedClass(model, 'JobClass', 1, delay, 0);
hitClass = ClosedClass(model, 'HitClass', 0, delay, 0);
missClass = ClosedClass(model, 'MissClass', 0, delay, 0);

delay.setService(jobClass, Exp(1));

pAccess = Zipf(1.2, n);  % Zipf-like item references
cacheNode.setRead(jobClass, pAccess);

cacheNode.setHitClass(jobClass, hitClass);
cacheNode.setMissClass(jobClass, missClass);

P = model.initRoutingMatrix;
P{jobClass, jobClass}(delay, cacheNode) =  1.0;
P{hitClass, jobClass}(cacheNode, delay) =  1.0;
P{missClass, jobClass}(cacheNode, delay) =  1.0;
model.link(P);

solver{1} = CTMC(model,'keep',false);
AvgTable{1} = solver{1}.getAvgNodeTable; AvgTable{1}

model.reset;
solver{2} = SSA(model,'samples',1e4,'method','serial','seed',23000);
AvgTable{2} = solver{2}.getAvgNodeTable; AvgTable{2}
