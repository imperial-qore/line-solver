% Cache with per-item storage costs (sizes) and per-list cost caps.
%
% Each item i carries a storage cost sigma_i and list j may hold items of
% total cost at most k_j. A promotion that would breach a cap serves the
% request without changing the cache state. SolverNC evaluates the
% constrained normalizing constant E(m,k) of Casale-Gast (IEEE/ACM ToN
% 29(2), 2021), Sec. IX.

clear solver AvgTable;

model = Network('model');

n = 6;          % number of items
m = [1 1];      % two lists, one item each
sizes = [1 1 1 2 2 2];  % small and large items
caps  = [2 1];  % list 2 admits small items only

source = Source(model, 'Source');
cacheNode = Cache(model, 'Cache', n, m, ReplacementStrategy.RR);
sink = Sink(model, 'Sink');

jobClass = OpenClass(model, 'InitClass', 0);
hitClass = OpenClass(model, 'HitClass', 0);
missClass = OpenClass(model, 'MissClass', 0);

source.setArrival(jobClass, Exp(2));
cacheNode.setRead(jobClass, DiscreteSampler((1/n)*ones(1,n)));
cacheNode.setItemSizes(sizes);
cacheNode.setCostCaps(caps);

cacheNode.setHitClass(jobClass, hitClass);
cacheNode.setMissClass(jobClass, missClass);

P = model.initRoutingMatrix;
P{jobClass, jobClass}(source, cacheNode) = 1.0;
P{hitClass, hitClass}(cacheNode, sink) = 1.0;
P{missClass, missClass}(cacheNode, sink) = 1.0;
model.link(P);

%%
solver = NC(model,'method','exact');
AvgTable = solver.getAvgNodeTable; AvgTable
CacheTable = solver.getAvgCacheTable; CacheTable
ItemTable = solver.getAvgItemTable; ItemTable

hitRatio = cacheNode.getHitRatio
listCost = cacheNode.getListCost
