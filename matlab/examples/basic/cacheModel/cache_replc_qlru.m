clear solver AvgTable;

% q-LRU replacement: on a miss the item is admitted (LRU head insert) with
% probability q, otherwise it passes through uncached. Admission filtering can
% raise the hit ratio over plain LRU under skewed popularity. Exact in CTMC;
% simulated in SSA/LDES. Not product-form (MVA/NC/FLD reject it).

model = Network('model');

n = 5; % number of items
m = 2; % cache capacity

delay = Delay(model, 'Delay');
cacheNode = Cache(model, 'Cache', n, m, ReplacementStrategy.QLRU);
cacheNode.setAdmissionProb(0.5); % admit a missed item with probability q=0.5

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

solver{1} = CTMC(model, 'exact','keep',false);
AvgTable{1} = solver{1}.getAvgNodeTable; AvgTable{1}

model.reset;
solver{2} = SSA(model,'samples',1e4,'method','serial','seed',23000);
AvgTable{2} = solver{2}.getAvgNodeTable; AvgTable{2}
