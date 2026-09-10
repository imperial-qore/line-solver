function model = gallery_cache_lru()
% GALLERY_CACHE_LRU Closed cache model with LRU replacement (n=5 items, m=2 slots)
model = Network('Cache-LRU');
%% Block 1: nodes
n = 5; % number of items
m = 2; % cache capacity
delay = Delay(model, 'Delay');
cacheNode = Cache(model, 'Cache', n, m, ReplacementStrategy.LRU);
%% Block 2: classes
jobClass = ClosedClass(model, 'JobClass', 1, delay, 0);
hitClass = ClosedClass(model, 'HitClass', 0, delay, 0);
missClass = ClosedClass(model, 'MissClass', 0, delay, 0);
delay.setService(jobClass, Exp(1));
pAccess = DiscreteSampler((1/n)*ones(1,n));  % uniform item references
cacheNode.setRead(jobClass, pAccess);
cacheNode.setHitClass(jobClass, hitClass);
cacheNode.setMissClass(jobClass, missClass);
%% Block 3: topology
P = model.initRoutingMatrix();
P{jobClass, jobClass}(delay, cacheNode) = 1.0;
P{hitClass, jobClass}(cacheNode, delay) = 1.0;
P{missClass, jobClass}(cacheNode, delay) = 1.0;
model.link(P);
end
