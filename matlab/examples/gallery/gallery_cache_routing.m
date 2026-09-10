function model = gallery_cache_routing()
% GALLERY_CACHE_ROUTING Open cache with hit/miss routed to distinct queues
model = Network('Cache-Routing');
%% Block 1: nodes
n = 4; % number of items
m = 2; % cache capacity
source = Source(model, 'Source');
cacheNode = Cache(model, 'Cache', n, m, ReplacementStrategy.LRU);
hitQueue = Queue(model, 'HitQueue', SchedStrategy.FCFS);
missQueue = Queue(model, 'MissQueue', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');
%% Block 2: classes
jobClass = OpenClass(model, 'InitClass', 0);
hitClass = OpenClass(model, 'HitClass', 0);
missClass = OpenClass(model, 'MissClass', 0);
source.setArrival(jobClass, Exp(1));
hitQueue.setService(hitClass, Exp(2.0));
missQueue.setService(missClass, Exp(1.0));
pAccess = DiscreteSampler((1/n)*ones(1,n));  % uniform item references
cacheNode.setRead(jobClass, pAccess);
cacheNode.setHitClass(jobClass, hitClass);
cacheNode.setMissClass(jobClass, missClass);
%% Block 3: topology
P = model.initRoutingMatrix();
P{jobClass, jobClass}(source, cacheNode) = 1.0;
P{hitClass, hitClass}(cacheNode, hitQueue) = 1.0;
P{hitClass, hitClass}(hitQueue, sink) = 1.0;
P{missClass, missClass}(cacheNode, missQueue) = 1.0;
P{missClass, missClass}(missQueue, sink) = 1.0;
model.link(P);
end
